"""Tests for advanced analyses: .step, Monte Carlo, noise, FFT, B sources."""

import math

import numpy as np
import pytest

from minispice import NetlistParser, NoiseResult, Simulation


def _run(text: str):
    return Simulation(NetlistParser(text).parse()).run()


def test_step_param_reparses_expressions():
    results = _run(
        "Step param\n"
        ".param R=1k\n"
        "V1 in 0 1\n"
        "R1 in out {R}\n"
        "R2 out 0 1k\n"
        ".step param R LIST 1k 2k 3k\n"
        ".op\n.end\n"
    )
    assert [r.case["r"] for r in results] == [1000.0, 2000.0, 3000.0]
    assert [r.nodes["out"] for r in results] == pytest.approx([0.5, 1 / 3, 0.25])


def test_step_temp_runs_every_temperature():
    results = _run(
        "Step temp\n"
        "V1 in 0 1\n"
        "R1 in out 1k\n"
        "D1 out 0 dmod\n"
        ".model dmod D(IS=1e-14 N=1)\n"
        ".step temp 27 127 50\n"
        ".op\n.end\n"
    )
    assert [r.case["temp"] for r in results] == [27.0, 77.0, 127.0]
    voltages = [r.nodes["out"] for r in results]
    assert voltages[0] > voltages[1] > voltages[2]


def test_step_source_values():
    results = _run(
        "Step source\n"
        "V1 in 0 1\n"
        "R1 in 0 1k\n"
        ".step V1 LIST 0 1 2\n"
        ".op\n.end\n"
    )
    assert [r.case["V1"] for r in results] == [0.0, 1.0, 2.0]
    assert [r.nodes["in"] for r in results] == pytest.approx([0.0, 1.0, 2.0])


def test_monte_carlo_summary_and_measurements():
    simulation = Simulation(
        NetlistParser(
            "Monte Carlo\n"
            "V1 in 0 1\n"
            "R1 in out 1k TOL=0.1\n"
            "R2 out 0 1k TOL=0.1\n"
            ".dc V1 0 1 0.5\n"
            ".meas dc voutmax MAX V(out)\n"
            ".mc 30 SEED=7\n"
            ".end\n"
        ).parse()
    )
    results = simulation.run()
    assert len(results) == 30
    assert [r.case["mc"] for r in results] == list(range(1, 31))
    stats = simulation.mc_summary
    assert stats["V(out)"]["mean"] == pytest.approx(0.5, abs=0.02)
    assert stats["V(out)"]["std"] > 0.0
    assert stats["V(out)"]["min"] <= stats["V(out)"]["mean"] <= stats["V(out)"]["max"]
    assert "meas voutmax" in stats
    assert stats["V(in)"]["std"] == pytest.approx(0.0)


def test_behavioral_sources_dc_tran_and_ac():
    op = _run(
        "Behavioral\n"
        "V1 in 0 1\n"
        "B1 out 0 V={2*V(in)}\n"
        "R1 out 0 1k\n"
        "B2 gout 0 I={V(in)/1k}\n"
        "R2 gout 0 1k\n"
        ".op\n.end\n"
    )[0]
    assert op.nodes["out"] == pytest.approx(2.0)
    assert op.nodes["gout"] == pytest.approx(-1.0)

    ac = _run(
        "Behavioral AC\n"
        "V1 in 0 dc 0 ac 1\n"
        "B1 out 0 V={2*V(in)}\n"
        "R1 out 0 1k\n"
        ".ac lin 2 10 100\n.end\n"
    )[0]
    assert np.allclose(ac.nodes["out"], 2.0 + 0j)

    tran = _run(
        "Behavioral tran\n"
        "V1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\n"
        "B1 out 0 V={V(in,0)}\n"
        "R1 out 0 1k\n"
        ".tran 100u 1m\n.end\n"
    )[0]
    assert tran.nodes["out"][0] == pytest.approx(0.0, abs=1e-9)
    assert tran.nodes["out"][-1] == pytest.approx(1.0, abs=1e-3)

    ramp = _run(
        "Behavioral time\n"
        "V1 in 0 1\n"
        "B1 out 0 V={V(in)*time}\n"
        "R1 out 0 1k\n"
        ".tran 100u 1m\n.end\n"
    )[0]
    assert ramp.nodes["out"][0] == pytest.approx(0.0, abs=1e-12)
    assert ramp.nodes["out"][-1] == pytest.approx(1e-3, rel=1e-3)


def test_resistor_noise_matches_kT_over_C():
    result = _run(
        "RC noise\n"
        "V1 in 0 dc 0 ac 1\n"
        "R1 in out 1k\n"
        "C1 out 0 1u\n"
        ".noise V(out) V1 dec 20 1 1e9\n.end\n"
    )[0]
    assert isinstance(result, NoiseResult)
    kT = 1.380649e-23 * 300.15
    assert result.output_psd[0] == pytest.approx(4 * kT * 1e3, rel=2e-3)
    assert result.total_rms == pytest.approx(math.sqrt(kT / 1e-6), rel=0.03)
    assert result.input_psd is not None
    assert result.input_unit == "V"
    assert result.input_psd[0] == pytest.approx(4 * kT * 1e3, rel=2e-3)
    assert result.input_rms > result.total_rms


def test_noise_current_source_input_referred_units():
    result = _run(
        "RC noise, current input\n"
        "I1 0 out dc 0 ac 1\n"
        "R1 out 0 1k\n"
        "C1 out 0 1u\n"
        ".noise V(out) I1 dec 20 1 1e9\n.end\n"
    )[0]
    kT = 1.380649e-23 * 300.15
    assert result.input_unit == "A"
    assert result.input_psd is not None
    assert result.output_psd[0] == pytest.approx(4 * kT * 1e3, rel=2e-3)
    assert result.input_psd[0] == pytest.approx(4 * kT / 1e3, rel=2e-3)


def test_monte_carlo_noise_summary():
    simulation = Simulation(
        NetlistParser(
            "MC noise\n"
            "V1 in 0 dc 0 ac 1\n"
            "R1 in out 1k TOL=0.2\n"
            "C1 out 0 1u\n"
            ".noise V(out) V1 dec 5 1 1e6\n"
            ".mc 4 SEED=1\n.end\n"
        ).parse()
    )
    results = simulation.run()
    assert len(results) == 4
    assert [r.case["mc"] for r in results] == [1, 2, 3, 4]
    assert "noise_out_rms" in simulation.mc_summary
    assert "noise_in_rms" in simulation.mc_summary
    assert simulation.mc_summary["noise_out_rms"]["std"] > 0.0


def test_transient_fft_and_four_measurements():
    result = _run(
        "FFT\n"
        "V1 out 0 SIN(0 1 1k)\n"
        "R1 out 0 1k\n"
        ".tran 1u 10m\n"
        ".four 1k V(out)\n.end\n"
    )[0]
    freqs, spectrum = result.fft("v(out)")
    peak = freqs[int(np.argmax(spectrum))]
    assert peak == pytest.approx(1e3, rel=0.02)
    assert spectrum.max() == pytest.approx(1.0, rel=0.02)
    assert result.measurements["V_out_f1"] == pytest.approx(1.0, rel=0.05)
    assert result.measurements["thd_V_out"] < 0.02


def test_bjt_high_injection_knee_lowers_beta():
    def beta(model: str) -> float:
        result = _run(
            "BJT knee\n"
            "VCC vcc 0 5\n"
            "I1 0 b 10u\n"
            f"Q1 c b 0 qmod\n"
            "RC vcc c 1k\n"
            f".model qmod {model}\n"
            ".op\n.end\n"
        )[0]
        return float(result.details["Q1"]["beta"])

    ideal = beta("NPN(IS=1e-15 BF=100)")
    knee = beta("NPN(IS=1e-15 BF=100 IKF=0.5m)")
    assert ideal == pytest.approx(100.0, rel=2e-2)
    assert knee < 0.85 * ideal


def test_bjt_nonideal_base_current_lowers_beta():
    def beta(model: str) -> float:
        result = _run(
            "BJT ise\n"
            "VCC vcc 0 5\n"
            "VIN in 0 0.65\n"
            "RB in b 100k\n"
            "RC vcc c 1k\n"
            "Q1 c b 0 qmod\n"
            f".model qmod {model}\n"
            ".op\n.end\n"
        )[0]
        return float(result.details["Q1"]["beta"])

    assert beta("NPN(IS=1e-15 BF=100)") > 100.0
    assert beta("NPN(IS=1e-15 BF=100 ISE=1e-14 NE=1.5)") < 95.0


def test_bjt_area_scales_currents_and_knee():
    def collector(area: float, model: str = "NPN(IS=1e-15 BF=100 IKF=0.5m)") -> float:
        result = _run(
            "BJT area\n"
            "VCC vcc 0 5\n"
            "VIN b 0 0.65\n"
            "RC vcc c 1k\n"
            f"Q1 c b 0 qmod {area}\n"
            f".model qmod {model}\n"
            ".op\n.end\n"
        )[0]
        return float(result.currents["Q1"])

    assert collector(2.0) == pytest.approx(2.0 * collector(1.0), rel=1e-3)


def _bjt_op(model: str, area: float = 1.0):
    text = (
        "BJT parasitics\n"
        "VCC vcc 0 5\n"
        "VIN in 0 0.7\n"
        "RB in b 100k\n"
        "RC vcc c 1k\n"
        f"Q1 c b 0 qmod {area}\n"
        f".model qmod {model}\n"
        ".op\n.end\n"
    )
    simulation = Simulation(NetlistParser(text).parse())
    return simulation, simulation.run()[0]


def test_bjt_series_resistances_drop_voltage_at_internal_nodes():
    plain_sim, plain = _bjt_op("NPN(IS=1e-15 BF=100)")
    sim, result = _bjt_op("NPN(IS=1e-15 BF=100 RB=100k RC=100 RE=10)")

    assert plain_sim.internal_nodes == []
    assert sim.internal_nodes == ["q1#c", "q1#b", "q1#e"]
    assert all("#" not in name for name in result.nodes)

    # The parasitic base resistance lowers the intrinsic Vbe and Ic.
    assert result.details["Q1"]["vbe"] < plain.details["Q1"]["vbe"] - 0.005
    assert result.currents["Q1"] < plain.currents["Q1"]

    # Series collector/emitter resistors show up as voltage drops that are
    # hidden from the public node list but kept in the solution vector.
    ic = float(result.currents["Q1"])
    ie = -float(result.details["Q1"]["ie"])
    vc = result.vector[sim.index["q1#c"]]
    ve = result.vector[sim.index["q1#e"]]
    assert result.nodes["c"] - vc == pytest.approx(ic * 100.0, rel=1e-4)
    assert ve == pytest.approx(ie * 10.0, rel=1e-4)


def test_bjt_parasitic_resistances_scale_with_area():
    sim, _ = _bjt_op("NPN(IS=1e-15 BF=100 RB=100k)", area=2.0)
    base_resistor = next(el for el in sim._bjt_parasitics if el.name == "Q1#rb")
    assert base_resistor.value == pytest.approx(50e3)


def test_bjt_junction_capacitances_roll_off_ac_gain():
    def response(model: str) -> tuple[np.ndarray, np.ndarray]:
        result = _run(
            "BJT ac\n"
            "VCC vcc 0 5\n"
            "VIN in 0 0.7 ac 1\n"
            "RB in b 10k\n"
            "RC vcc c 1k\n"
            "Q1 c b 0 qmod\n"
            f".model qmod {model}\n"
            ".ac dec 10 1 1e10\n.end\n"
        )[0]
        return result.freqs, np.abs(result.nodes["c"])

    freqs, plain = response("NPN(IS=1e-15 BF=100)")
    _, capped = response("NPN(IS=1e-15 BF=100 CJC=10p)")
    index = int(np.argmin(np.abs(freqs - 1e8)))
    assert capped[0] == pytest.approx(plain[0], rel=1e-6)
    assert capped[index] < 0.99 * plain[index]


def test_bjt_charge_storage_slows_transient_edges():
    text = (
        "BJT tran\n"
        "VCC vcc 0 5\n"
        "VIN in 0 PULSE(0 0.8 10n 1n 1n 1u 2u)\n"
        "RB in b 10k\n"
        "RC vcc c 10k\n"
        "Q1 c b 0 qmod\n"
        ".model qmod NPN(IS=1e-15 BF=100 {extra})\n"
        ".tran 1n 3u\n.end\n"
    )

    def fall_time(extra: str) -> tuple[float, np.ndarray, np.ndarray]:
        result = _run(text.format(extra=extra))[0]
        y = result.nodes["c"]
        hits = np.where(y <= 2.5)[0]
        return (result.time[hits[0]] if len(hits) else math.nan), result.time, y

    fast, _, _ = fall_time("")
    slow, _, _ = fall_time("CJE=50p")
    assert slow > 20.0 * fast

    # CJC is multiplied by the Miller gain, keeping the collector high.
    miller, _, _ = fall_time("CJC=50p")
    assert math.isnan(miller)

    _, _, stored = fall_time("TF=50n")
    _, time, plain = fall_time("")
    after = time > 1.02e-6
    recovered = time[np.where(after & (stored >= 4.0))[0][0]]
    baseline = time[np.where(after & (plain >= 4.0))[0][0]]
    assert recovered > baseline + 50e-9


def test_bjt_parasitic_members_are_parsed():
    circuit = NetlistParser(
        "BJT model\n"
        ".model qmod NPN(RB=10 RC=2 RE=1 CJE=1p CJC=2p TF=3n TR=4n)\n"
        ".end\n"
    ).parse()
    model = circuit.models["qmod"]
    assert (model.rb, model.rc, model.re) == (10.0, 2.0, 1.0)
    assert (model.cje, model.cjc, model.tf, model.tr) == pytest.approx(
        (1e-12, 2e-12, 3e-9, 4e-9)
    )



def test_sparse_linear_solver_path_is_correct():
    simulation = Simulation(NetlistParser("Sparse\nV1 in 0 1\n.end\n").parse())
    size = 80
    matrix = np.eye(size)
    rhs = np.arange(size, dtype=float)
    assert np.allclose(simulation._solve_linear(matrix, rhs), rhs)

    lines = ["Sparse ladder", "V1 in 0 1", "R0 in n1 1k"]
    lines += [f"R{i} n{i} n{i + 1} 1k" for i in range(1, 70)]
    lines += ["R70 n70 0 1k", ".op", ".end"]
    ladder = NetlistParser("\n".join(lines)).parse()
    assert Simulation(ladder).size >= 64
    result = Simulation(ladder).run()[0]
    assert result.nodes["n70"] == pytest.approx(1.0 / 71.0, rel=1e-5)


def test_step_and_monte_carlo_combine():
    simulation = Simulation(
        NetlistParser(
            "Step + Monte Carlo\n"
            "V1 in 0 1\n"
            ".param R=1k\n"
            "R1 in out {R}\n"
            "R2 out 0 1k\n"
            ".step param R LIST 1k 2k\n"
            ".mc 10 SEED=3\n"
            ".op\n.end\n"
        ).parse()
    )
    results = simulation.run()
    assert len(results) == 20
    assert {r.case["r"] for r in results} == {1000.0, 2000.0}
    assert [r.case["mc"] for r in results if r.case["r"] == 1000.0] == list(
        range(1, 11)
    )
    expected = {1000.0: 0.5, 2000.0: 1 / 3}
    for stepped, nominal in expected.items():
        subset = [r.nodes["out"] for r in results if r.case["r"] == stepped]
        assert np.mean(subset) == pytest.approx(nominal, abs=0.03)
        assert np.std(subset) > 0.0
    assert "V(out)" in simulation.mc_summary


def test_measure_limits_and_yield_summary():
    circuit = NetlistParser(
        "MC yield\n"
        "V1 in 0 1\n"
        "R1 in out 1k TOL=0.1\n"
        "R2 out 0 1k TOL=0.1\n"
        ".dc V1 0.5 1.5 0.5\n"
        ".meas dc vout MAX V(out) MIN=0.4 MAX=0.9\n"
        ".mc 8 SEED=4\n.end\n"
    ).parse()
    spec = circuit.measures[0]
    assert (spec.limit_min, spec.limit_max) == (0.4, 0.9)

    simulation = Simulation(circuit)
    results = simulation.run()
    assert len(results) == 8
    assert simulation.yield_summary["vout"]["total"] == 8
    assert simulation.yield_summary["vout"]["yield"] == pytest.approx(1.0)
    assert simulation.yield_by_case["nominal"]["vout"]["pass"] == 8
    assert all(np.isfinite(result.measurements["vout"]) for result in results)


def test_measure_rejects_inverted_limits():
    from minispice import NetlistError

    with pytest.raises(NetlistError):
        NetlistParser(
            "bad\nV1 in 0 1\nR1 in out 1k\n.dc V1 0 1 1\n"
            ".meas dc vout MAX V(out) MIN=0.9 MAX=0.1\n.end\n"
        ).parse()


def test_transient_edge_width_and_cross_measurements():
    result = _run(
        "RC measures\n"
        "V1 in 0 PULSE(0 1 0 1n 1n 5m 10m)\n"
        "R1 in out 1k\n"
        "C1 out 0 1u\n"
        ".tran 10u 20m\n"
        ".meas tran tr RISE V(out)\n"
        ".meas tran tf FALL V(out) FROM=4m TO=20m\n"
        ".meas tran tw WIDTH V(out) VAL=0.5\n"
        ".meas tran tc CROSS V(out) VAL=0.5 RISE=1\n"
        ".end\n"
    )[0]
    tau = 1e-3
    assert result.measurements["tr"] == pytest.approx(tau * math.log(9), rel=0.03)
    assert result.measurements["tf"] == pytest.approx(tau * math.log(9), rel=0.03)
    assert result.measurements["tw"] == pytest.approx(5e-3, rel=0.01)
    assert result.measurements["tc"] == pytest.approx(tau * math.log(2), rel=0.02)


def test_ac_gain_phase_bandwidth_and_loopgain():
    result = _run(
        "RC ac measures\n"
        "V1 in 0 dc 0 ac 1\n"
        "R1 in out 1k\n"
        "C1 out 0 1u\n"
        ".ac dec 50 1 1e5\n"
        ".meas ac gain GAIN V(out) AT=1\n"
        ".meas ac ph PHASE V(out) AT=159.15\n"
        ".meas ac bw BW V(out)\n"
        ".end\n"
    )[0]
    fc = 1.0 / (2 * math.pi * 1e3 * 1e-6)
    assert result.measurements["gain"] == pytest.approx(0.0, abs=1e-3)
    assert result.measurements["ph"] == pytest.approx(-45.0, abs=0.2)
    assert result.measurements["bw"] == pytest.approx(fc, rel=0.01)

    loop = _run(
        "Loop gain\n"
        "V1 in 0 dc 0 ac 1\n"
        "E1 out 0 in 0 2\n"
        ".ac lin 3 1 1e3\n"
        ".meas ac lg LOOPGAIN V(out) V(in) AT=1k\n"
        ".end\n"
    )[0]
    assert loop.measurements["lg"] == pytest.approx(20 * math.log10(2), rel=1e-6)


def test_advanced_measure_validation():
    from minispice import NetlistError

    with pytest.raises(NetlistError):
        NetlistParser(
            "bad\nV1 in 0 ac 1\nR1 in out 1k\nC1 out 0 1u\n"
            ".ac dec 5 1 1e3\n.meas ac tr RISE V(out)\n.end\n"
        ).parse()
    with pytest.raises(NetlistError):
        NetlistParser(
            "bad\nV1 in 0 1\nR1 in out 1k\n.tran 1u 1m\n"
            ".meas tran g GAIN V(out) AT=1k\n.end\n"
        ).parse()


def test_expression_if_func_and_comparisons():
    circuit = NetlistParser(
        "Expressions\n"
        ".func sq(x) = {x*x}\n"
        ".func gain(x) = sq(x) / 2\n"
        ".param A = 3\n"
        ".param B = {gain(A)}\n"
        ".param C = {if(B > 4, 10, 20)}\n"
        ".param D = {if(B > 100, 1/0, 99)}\n"
        ".param E = {if(A > 1 and A < 5, 7, 8)}\n"
        ".end\n"
    ).parse()
    assert circuit.params["b"] == pytest.approx(4.5)
    assert circuit.params["c"] == 10.0
    # if() is lazy, so the impossible branch never raises.
    assert circuit.params["d"] == 99.0
    assert circuit.params["e"] == 7.0
    assert circuit.functions["sq"] == (["x"], "x*x")


def test_func_works_in_devices_and_behavioral_sources():
    result = _run(
        "Func usage\n"
        ".func twice(x) = {2*x}\n"
        "V1 in 0 1\n"
        "R1 in out {twice(1k)}\n"
        "R2 out 0 1k\n"
        "B1 buf 0 V={twice(V(out))}\n"
        "R3 buf 0 1k\n"
        ".op\n.end\n"
    )[0]
    assert result.nodes["out"] == pytest.approx(1.0 / 3.0)
    assert result.nodes["buf"] == pytest.approx(2.0 / 3.0)


def test_gauss_expression_is_seeded_and_reproducible():
    text = (
        "Gauss param\n"
        ".param R = {gauss(1k, 10)}\n"
        "V1 in 0 1\n"
        "R1 in out {R}\n"
        "R2 out 0 1k\n"
        ".op\n.end\n"
    )
    first = NetlistParser(text).parse().params["r"]
    second = NetlistParser(text).parse().params["r"]
    seeded_a = NetlistParser(text, expr_seed=5).parse().params["r"]
    seeded_b = NetlistParser(text, expr_seed=5).parse().params["r"]
    assert first == pytest.approx(second)
    assert seeded_a == pytest.approx(seeded_b)
    assert seeded_a != pytest.approx(first)


def test_gauss_is_unavailable_at_runtime():
    from minispice import NetlistError

    circuit = NetlistParser(
        "B source gauss\n"
        "V1 in 0 1\n"
        "B1 out 0 V={gauss(0, 1)}\n"
        ".op\n.end\n"
    ).parse()
    with pytest.raises(NetlistError):
        Simulation(circuit).run()


def _gauss_divider(rho: float, runs: int = 150) -> tuple[float, float]:
    circuit = NetlistParser(
        "Gaussian MC\n"
        "V1 in 0 1\n"
        "R1 in out 1k TOL=0.1\n"
        "R2 out 0 1k TOL=0.1\n"
        ".dc V1 0.9 1.0 0.1\n"
        ".meas dc vout MAX V(out)\n"
        f".mc {runs} SEED=7 DIST=GAUSS SIGMA=0.05 RHO={rho}\n.end\n"
    ).parse()
    results = Simulation(circuit).run()
    values = np.array([result.measurements["vout"] for result in results])
    return float(np.mean(values)), float(np.std(values))


def test_gaussian_monte_carlo_statistics_and_correlation():
    mean, spread = _gauss_divider(0.0)
    assert mean == pytest.approx(0.5, abs=0.005)
    assert 0.005 < spread < 0.03

    # RHO=1 moves both resistors together, so the divider ratio is frozen.
    correlated_mean, correlated_spread = _gauss_divider(1.0, runs=20)
    assert correlated_mean == pytest.approx(0.5, abs=1e-6)
    assert correlated_spread == pytest.approx(0.0, abs=1e-6)


def test_gaussian_option_validation():
    from minispice import NetlistError

    with pytest.raises(NetlistError):
        NetlistParser(
            "bad\nV1 in 0 1\nR1 in out 1k\n.op\n"
            ".mc 5 RHO=0.5\n.end\n"
        ).parse()
    with pytest.raises(NetlistError):
        NetlistParser(
            "bad\nV1 in 0 1\nR1 in out 1k\n.op\n"
            ".mc 5 DIST=NOPE\n.end\n"
        ).parse()


def test_step_and_gaussian_monte_carlo_combine():
    circuit = NetlistParser(
        "Step gaussian MC\n"
        ".param R=1k\n"
        "V1 in 0 1\n"
        "R1 in out {R}\n"
        "R2 out 0 1k\n"
        ".dc V1 0.9 1.0 0.1\n"
        ".meas dc vout MAX V(out) MIN=0.2 MAX=0.8\n"
        ".step param R LIST 1k 2k\n"
        ".mc 8 SEED=3 DIST=GAUSS\n.end\n"
    ).parse()
    simulation = Simulation(circuit)
    results = simulation.run()
    assert len(results) == 16
    assert [r.case["mc"] for r in results if r.case["r"] == 1000.0] == list(
        range(1, 9)
    )
    assert set(simulation.yield_by_case) == {"r=1000", "r=2000"}
    assert simulation.yield_summary["vout"]["total"] == 16


def test_sensitivity_ranking_and_worst_case_corners():
    circuit = NetlistParser(
        "Divider sensitivity\n"
        "V1 in 0 1\n"
        "R1 in out 1k TOL=0.1\n"
        "R2 out 0 1k TOL=0.1\n"
        ".sens V(out) TOP=5 DEVICES=2\n"
        ".op\n.end\n"
    ).parse()
    spec = circuit.sensitivity[0]
    assert spec.output == "V(out)"
    assert spec.top == 5 and spec.devices == 2

    simulation = Simulation(circuit)
    simulation.run()
    report = simulation.sensitivity_report[0]
    assert report["nominal"] == pytest.approx(0.5)
    entries = {entry["name"]: entry for entry in report["entries"]}
    assert entries["R1"]["sensitivity"] == pytest.approx(-0.5, abs=1e-4)
    assert entries["R2"]["sensitivity"] == pytest.approx(0.5, abs=1e-4)

    worst = report["worst_case"]
    assert worst["min"]["value"] == pytest.approx(0.45, rel=1e-5)
    assert worst["max"]["value"] == pytest.approx(0.55, rel=1e-5)
    assert worst["min"]["case"] == {"R1": pytest.approx(1100.0), "R2": pytest.approx(900.0)}


def test_sensitivity_validation():
    from minispice import NetlistError

    with pytest.raises(NetlistError):
        NetlistParser(
            "bad\nV1 in 0 1\nR1 in out 1k\n.sens\n.op\n.end\n"
        ).parse()
    with pytest.raises(NetlistError):
        NetlistParser(
            "bad\nV1 in 0 1\nR1 in out 1k\n.sens V(out) BOGUS=1\n.op\n.end\n"
        ).parse()


def test_source_stepping_fallback_converges(monkeypatch):
    from minispice.solver import ConvergenceError

    original = Simulation._newton
    state = {"started": False}

    def flaky(self, ctx, x0=None):
        if ctx.source_scale < 1.0:
            state["started"] = True
            return original(self, ctx, x0)
        if state["started"]:
            return original(self, ctx, x0)
        raise ConvergenceError("forced direct failure")

    monkeypatch.setattr(Simulation, "_newton", flaky)
    simulation = Simulation(
        NetlistParser(
            "Fallback\nV1 in 0 1\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n"
        ).parse()
    )
    result = simulation.operating_point()
    assert state["started"] is True
    assert result.nodes["out"] == pytest.approx(0.5)


def test_convergence_failure_reports_residuals(monkeypatch):
    from minispice.solver import ConvergenceError

    def always_fail(self, ctx, x0=None):
        raise ConvergenceError("forced")

    monkeypatch.setattr(Simulation, "_newton", always_fail)
    simulation = Simulation(
        NetlistParser(
            "No convergence\nV1 in 0 1\nR1 in out 1k\nR2 out 0 1k\n.op\n.end\n"
        ).parse()
    )
    with pytest.raises(ConvergenceError) as excinfo:
        simulation.operating_point()
    message = str(excinfo.value)
    assert "source/gmin stepping failed" in message
    assert "largest residuals:" in message
