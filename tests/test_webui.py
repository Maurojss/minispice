"""Tests for the browser GUI (webui module)."""

import json
import re
import threading
import urllib.request

import pytest

from minispice.webui import (
    _csv_bytes,
    _dat_payload,
    _fmt,
    _histogram_png,
    _plot,
    _result_payload,
    create_server,
)


def _get_page(base: str) -> str:
    with urllib.request.urlopen(f"{base}/", timeout=10) as response:
        return response.read().decode()


def test_page_script_is_valid_js():
    """Guards against the whole script dying from a syntax error."""
    esprima = pytest.importorskip("esprima")

    server = create_server()
    threading.Thread(target=server.serve_forever, daemon=True).start()
    host, port = server.server_address[:2]
    try:
        page = _get_page(f"http://{host}:{port}/")
    finally:
        server.shutdown()
        server.server_close()
    assert "__DEFAULT__" not in page and "__EXAMPLES__" not in page
    script = re.search(r"<script>(.*)</script>", page, re.S).group(1)
    esprima.parseScript(script)  # raises on syntax error


def test_fmt():
    assert _fmt(0.993228) == "0.993228"
    assert _fmt(0) == "0"
    assert _fmt(-0.0) == "0"


def test_parse_gui_args():
    from minispice.webui import DEFAULT_PORT, _parse_gui_args

    assert _parse_gui_args(None) == (DEFAULT_PORT, [])
    assert _parse_gui_args(["-gui"]) == (DEFAULT_PORT, ["-gui"])
    assert _parse_gui_args(["--port", "9999", "x.cir"]) == (9999, ["x.cir"])
    assert _parse_gui_args(["--port=5123", "x.cir"]) == (5123, ["x.cir"])
    assert _parse_gui_args(["x.cir"]) == (DEFAULT_PORT, ["x.cir"])


def _post(url: str, body: dict) -> dict:
    request = urllib.request.Request(
        url,
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=10) as response:
        return json.loads(response.read())


def test_schematic_and_zoom(tmp_path):
    server = create_server()
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address[:2]
    base = f"http://{host}:{port}"
    try:
        netlist = "t\nV1 in 0 1\nR1 in out 1k\nC1 out 0 1u\n.tran 10u 1m\n.end\n"
        data = _post(f"{base}/run", {"netlist": netlist})
        assert data["ok"] and data["schematic"]
        run = data["run"]

        # schematic endpoint serves a PNG
        with urllib.request.urlopen(f"{base}/schematic?run={run}", timeout=10) as r:
            assert r.headers["Content-Type"] == "image/png"
            assert r.read(4) == b"\x89PNG"

        # zoomed plot only contains the first half of the sweep
        with urllib.request.urlopen(
            f"{base}/plot?run={run}&i=0&sigs=V:out&x0=0&x1=5e-4", timeout=10
        ) as r:
            assert r.read(4) == b"\x89PNG"

        # save=1 marks the image as an attachment
        with urllib.request.urlopen(
            f"{base}/plot?run={run}&i=0&sigs=V:out&save=1", timeout=10
        ) as r:
            assert "attachment" in r.headers["Content-Disposition"]

        # /dat serves full arrays for the browser cursor
        with urllib.request.urlopen(f"{base}/dat?run={run}&i=0", timeout=10) as r:
            dat = json.loads(r.read())
        assert dat["axis"] == "time"
        assert len(dat["x"]) == 101
        assert "V:out" in dat["series"]
        assert abs(dat["series"]["V:out"][-1] - 1.0) < 1e-3   # DC remains constant
    finally:
        server.shutdown()
        server.server_close()


def test_server_smoke(tmp_path):
    server = create_server()
    host, port = server.server_address[:2]
    base = f"http://{host}:{port}"
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        # main page renders
        with urllib.request.urlopen(f"{base}/", timeout=10) as r:
            assert r.status == 200
            page = r.read().decode()
            assert "MiniSPICE" in page and "__DEFAULT__" not in page
            assert "__EXAMPLES__" not in page

        # example download
        with urllib.request.urlopen(f"{base}/example?name=rc_lowpass.cir", timeout=10) as r:
            assert b"PULSE" in r.read()

        # run the default RC circuit
        data = _post(f"{base}/run", {"netlist": "RC\nV1 in 0 1\nR1 in out 1k\nC1 out 0 1u\n.tran 10u 1m\n.end\n"})
        assert data["ok"]
        run = data["run"]
        assert [r["type"] for r in data["results"]] == ["series"]

        # plot endpoint
        with urllib.request.urlopen(f"{base}/plot?run={run}&i=0&sigs=V:in,V:out", timeout=10) as r:
            assert r.headers["Content-Type"] == "image/png"
            assert r.read(4) == b"\x89PNG"

        # csv endpoint
        with urllib.request.urlopen(f"{base}/csv?run={run}&i=0", timeout=10) as r:
            csv_text = r.read().decode()
        assert csv_text.startswith("time,V(in),V(out)")
        assert "1.0" in csv_text

        # parse error becomes a clean error payload
        failure = _post(f"{base}/run", {"netlist": "BROKEN\nX1 a b\n.end\n"})
        assert not failure["ok"] and failure["error"]

        # empty netlist rejected
        empty = _post(f"{base}/run", {"netlist": "  \n"})
        assert not empty["ok"]

        # unknown route / unknown run
        with pytest.raises(urllib.error.HTTPError):
            urllib.request.urlopen(f"{base}/plot?run=9999", timeout=10)
    finally:
        server.shutdown()
        server.server_close()


def test_result_payload_and_csv_shapes():
    from minispice import NetlistParser, Simulation

    netlist = "R\nVIN in 0 5\nR1 in out 100\nR2 out 0 100\n.op\n.dc VIN 0 5 1\n.end\n"
    results = Simulation(NetlistParser(text=netlist).parse()).run()
    assert len(results) == 2

    op_payload = _result_payload(0, results[0])
    assert op_payload["type"] == "op"
    labels = [row[0] for row in op_payload["rows"]]
    assert "V(in)" in labels and "I(VIN)" in labels

    dc_payload = _result_payload(1, results[1])
    assert dc_payload["type"] == "series"
    assert dc_payload["head"][0] == "VIN"
    assert len(dc_payload["rows"]) == 6

    png = __import__("io").BytesIO()
    _plot(results[1], ["V:out"], png)
    assert png.getvalue()[:4] == b"\x89PNG"

    csv_text = _csv_bytes(results[1]).decode()
    assert csv_text.startswith("VIN,V(in),V(out)")


def test_noise_payload_and_csv_shapes():
    import io

    from minispice import NetlistParser, Simulation

    netlist = (
        "N\nV1 in 0 dc 0 ac 1\nR1 in out 1k\nC1 out 0 1u\n"
        ".noise V(out) V1 dec 5 1 1e6\n.end\n"
    )
    result = Simulation(NetlistParser(text=netlist).parse()).run()[0]

    payload = _result_payload(0, result)
    assert payload["head"] == ["freq (Hz)", "N:out", "N:in"]
    assert len(payload["rows"]) == 31
    assert [signal["label"] for signal in payload["signals"]] == ["N:out", "N:in"]
    assert "integrated output noise" in payload["info"]

    dat = _dat_payload(result)
    assert dat["axis"] == "freq (Hz)"
    assert set(dat["series"]) == {"N:out", "N:in"}

    csv_text = _csv_bytes(result).decode()
    assert csv_text.startswith(
        "freq,output_psd,output_density,input_psd,input_density"
    )

    png = io.BytesIO()
    _plot(result, ["N:out", "N:in"], png)
    assert png.getvalue()[:4] == b"\x89PNG"


def test_step_result_titles_include_case():
    from minispice import NetlistParser, Simulation

    netlist = (
        "S\nV1 in 0 1\n.param R=1k\nR1 in out {R}\nR2 out 0 1k\n"
        ".step param R LIST 1k 2k\n.op\n.end\n"
    )
    results = Simulation(NetlistParser(text=netlist).parse()).run()
    assert _result_payload(0, results[0])["title"] == "Operating point [r=1000]"
    assert _result_payload(1, results[1])["title"] == "Operating point [r=2000]"


def test_noise_plot_and_dat_endpoints():
    server = create_server()
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address[:2]
    base = f"http://{host}:{port}"
    try:
        netlist = (
            "N\nV1 in 0 dc 0 ac 1\nR1 in out 1k\nC1 out 0 1u\n"
            ".noise V(out) V1 dec 5 1 1e6\n.end\n"
        )
        data = _post(f"{base}/run", {"netlist": netlist})
        assert data["ok"]
        run = data["run"]
        with urllib.request.urlopen(
            f"{base}/plot?run={run}&i=0&sigs=N:out,N:in", timeout=10
        ) as response:
            assert response.read(4) == b"\x89PNG"
        with urllib.request.urlopen(f"{base}/dat?run={run}&i=0", timeout=10) as response:
            dat = json.loads(response.read())
        assert set(dat["series"]) == {"N:out", "N:in"}
    finally:
        server.shutdown()
        server.server_close()


def test_histogram_png_renders_limits():
    png = _histogram_png(
        [0.48, 0.5, 0.52, 0.49],
        "vout",
        {"mean": 0.4975, "limit_min": 0.45, "limit_max": 0.55},
    )
    assert png[:4] == b"\x89PNG"


def test_monte_carlo_payload_has_cases_yield_and_histogram():
    server = create_server()
    threading.Thread(target=server.serve_forever, daemon=True).start()
    host, port = server.server_address[:2]
    base = f"http://{host}:{port}"
    try:
        netlist = (
            "MC GUI\n"
            "V1 in 0 1\n"
            "R1 in out 1k TOL=0.1\n"
            "R2 out 0 1k TOL=0.1\n"
            ".dc V1 0.5 1.5 0.5\n"
            ".meas dc vout MAX V(out) MIN=0.4 MAX=0.9\n"
            ".mc 6 SEED=2\n.end\n"
        )
        data = _post(f"{base}/run", {"netlist": netlist})
        assert data["ok"]
        assert data["cases"] == [{"label": "nominal", "indices": list(range(6))}]
        stats = data["mc"]["measures"]["vout"]
        assert stats["total"] == 6
        assert stats["yield"] == pytest.approx(1.0)
        assert data["mc"]["by_case"]["nominal"]["vout"]["pass"] == 6

        run = data["run"]
        with urllib.request.urlopen(
            f"{base}/hist?run={run}&measure=vout&case=nominal", timeout=10
        ) as response:
            assert response.read(4) == b"\x89PNG"
    finally:
        server.shutdown()
        server.server_close()


def test_step_payload_groups_cases():
    server = create_server()
    threading.Thread(target=server.serve_forever, daemon=True).start()
    host, port = server.server_address[:2]
    base = f"http://{host}:{port}"
    try:
        netlist = (
            "Step GUI\n"
            ".param R=1k\n"
            "V1 in 0 1\n"
            "R1 in out {R}\n"
            "R2 out 0 1k\n"
            ".step param R LIST 1k 2k\n"
            ".op\n.end\n"
        )
        data = _post(f"{base}/run", {"netlist": netlist})
        assert data["ok"]
        assert [case["label"] for case in data["cases"]] == ["r=1000", "r=2000"]
        assert len(data["results"]) == 2
        assert data["mc"] is None
    finally:
        server.shutdown()
        server.server_close()
