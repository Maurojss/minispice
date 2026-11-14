"""Browser-based graphical interface for MiniSPICE.

Runs a tiny local web server (standard library only) that serves an HTML
page with a netlist editor, runs the declared analyses and renders the
results as tables, an inline matplotlib plot of the selected signals and
CSV download links.

Usage:
    minispice --gui                # open the browser UI
    python -m minispice.webui      # same, without installing
"""

from __future__ import annotations

import csv
import html
import io
import json
import shutil
import subprocess
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np

from .parser import NetlistError, NetlistParser
from .results import AcResult, DcSweep, NoiseResult, OperatingPoint, Transient
from .solver import ConvergenceError, Simulation

_ROW_PREVIEW = 60


def _examples_dir() -> Path:
    """Locate the examples/ folder (CWD first, then the repo layout)."""
    candidates = [Path.cwd() / "examples", Path(__file__).resolve().parents[2] / "examples"]
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    return candidates[0]

_PAGE = """<!DOCTYPE html>
<html lang="pt-BR">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MiniSPICE</title>
<style>
  :root { --bg:#12141a; --panel:#1c1f27; --line:#2c3040; --text:#d8dce6;
          --dim:#8a91a3; --accent:#4f8cff; --ok:#3fb96f; --err:#e05561; }
  * { box-sizing: border-box; }
  body { margin:0; background:var(--bg); color:var(--text);
         font:14px/1.5 system-ui, sans-serif; }
  header { display:flex; align-items:center; gap:10px; padding:12px 18px;
           background:var(--panel); border-bottom:1px solid var(--line); }
  header h1 { margin:0; font-size:17px; }
  header .sub { color:var(--dim); font-size:12px; }
  button { background:var(--accent); color:#fff; border:0; border-radius:6px;
           padding:8px 16px; font-size:14px; cursor:pointer; }
  button.secondary { background:#33384a; }
  button:disabled { opacity:.5; cursor:default; }
  main { display:grid; grid-template-columns: minmax(300px, 34%) 1fr;
         gap:14px; padding:14px 18px; height:calc(100vh - 53px); }
  .col { display:flex; flex-direction:column; gap:10px; min-height:0; }
  .editwrap { position:relative; flex:1; min-height:0; }
  #hl, #netlist {
    position:absolute; inset:0; margin:0; padding:10px;
    font-family:ui-monospace, monospace; font-size:13px; line-height:1.45;
    white-space:pre; border:1px solid var(--line); border-radius:6px;
    overflow:auto; background:var(--panel); tab-size:2;
  }
  #hl { color:var(--dim); pointer-events:none; user-select:none;
        overflow:hidden; }
  #netlist { color:transparent; caret-color:#fff; resize:none;
             outline-color:var(--accent); }
  #netlist::selection { background:#2b4a80; color:transparent; }
  .tk-c { color:#5f8f5f; }  .tk-k { color:#b48cff; font-weight:600; }
  .tk-n { color:#e8a06a; }  .tk-dev { color:#6fb7ff; }
  .zoomrow { display:flex; gap:6px; align-items:center; flex-wrap:wrap;
             margin-top:6px; font-size:12.5px; color:var(--dim); }
  .zoomrow input { width:76px; background:#0d0f14; color:var(--text);
                   border:1px solid var(--line); border-radius:4px; padding:3px 6px; }
  .plotwrap { position:relative; margin-top:8px; font-size:12.5px; color:var(--text); }
  .cursorcv { position:absolute; inset:0; width:100%; height:100%;
              cursor:crosshair; }
  .cursorread { position:absolute; left:6px; top:6px; pointer-events:none;
                background:rgba(13,15,20,.88); border:1px solid var(--line);
                border-radius:5px; padding:4px 9px;
                font:11.5px ui-monospace, monospace; white-space:pre; }
  #results { overflow:auto; display:grid; align-items:start;
             grid-template-columns:repeat(auto-fit, minmax(440px, 1fr)); gap:14px; }
  .resultbody { display:grid; grid-column:1/-1; align-items:start;
                grid-template-columns:repeat(auto-fit, minmax(440px, 1fr)); gap:14px; }
  .casebar { grid-column:1/-1; display:flex; gap:12px; align-items:center;
             flex-wrap:wrap; }
  .yield { color:var(--ok); font-weight:600; font-size:13px; margin:2px 0; }
  .histcard h3 { color:var(--text); }
  .card { background:var(--panel); border:1px solid var(--line);
          border-radius:8px; padding:12px 14px; overflow:auto; }
  .card h3 { margin:0 0 6px; font-size:14px; font-weight:600; }
  img.plot { width:100%; height:auto; background:#0d0f14;
             border-radius:6px; margin-top:8px; }
  .toolbar { display:flex; gap:8px; align-items:center; flex-wrap:wrap; }
  #status { font-size:12px; color:var(--dim); min-height:18px; }
  #status.err { color:var(--err); }  #status.ok { color:var(--ok); }
  table { border-collapse:collapse; font-family:ui-monospace, monospace;
          font-size:12.5px; margin-top:8px; }
  th, td { padding:2px 10px; text-align:right; border-bottom:1px solid var(--line); }
  th { color:var(--dim); font-weight:500; }
  td:first-child, th:first-child { text-align:left; }
  select { background:var(--panel); color:var(--text);
           border:1px solid var(--line); border-radius:6px; padding:8px;
           font-family:ui-monospace, monospace; font-size:13px; }
  .sigbox { display:flex; flex-wrap:wrap; gap:4px 12px; margin-top:6px;
            font-size:12.5px; }
  .sigbox label { cursor:pointer; }
  .note { color:var(--dim); font-size:12px; margin-top:6px; }
  a { color:var(--accent); text-decoration:none; }
  .warn { color:#e0a855; font-size:12.5px; }
</style>
</head>
<body>
<header>
  <h1>MiniSPICE</h1><span class="sub">simulador de circuitos</span>
  <span style="flex:1"></span>
  <button id="run">▶ Simular</button>
</header>
<main>
  <div class="col">
    <b style="font-size:13px">Netlist</b>
    <div class="toolbar">
      <select id="examples" style="flex:1"></select>
      <button class="secondary" id="load">Carregar</button>
      <button class="secondary" id="openBtn">Abrir…
        <input type="file" id="openFile" accept=".cir,.sp,.net,.txt" style="display:none">
      </button>
      <button class="secondary" id="save">Salvar…</button>
    </div>
    <div class="editwrap">
      <pre id="hl" aria-hidden="true"></pre>
      <textarea id="netlist" spellcheck="false">__DEFAULT__</textarea>
    </div>
    <div id="status"></div>
  </div>
  <div id="results" class="col">
    <div class="card" style="color:var(--dim)">
      Escolha ou edite uma netlist à esquerda e clique em <b>Simular</b>.
    </div>
  </div>
</main>
<script>
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g,
  c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const fmt = v => { const s = (+v).toPrecision(6);
  return s.includes('.') ? s.replace(/\\.?0+$/, '') : s; };

/* ---------- syntax highlighting (pre behind the transparent textarea) --- */
function hlLine(l){
  const s = l.trim();
  if (s.startsWith('*')) return `<span class="tk-c">${esc(l)}</span>`;
  if (s.startsWith('.')) return `<span class="tk-k">${esc(l)}</span>`;
  let out = '', last = 0, first = true;
  for (const m of l.matchAll(/(\\S+)/g)) {
    out += esc(l.slice(last, m.index));
    const w = m[0];
    let cls = '';
    if (/^[-+]?[.\\d]/.test(w)) cls = 'tk-n';
    else if (first) cls = 'tk-dev';
    out += cls ? `<span class="${cls}">${esc(w)}</span>` : esc(w);
    last = m.index + w.length; first = false;
  }
  return out + esc(l.slice(last));
}
function hl(){
  const ta = $('netlist'), pre = $('hl');
  pre.innerHTML = ta.value.replace(/\\t/g, '  ').split('\\n').map(hlLine).join('\\n');
  pre.scrollTop = ta.scrollTop; pre.scrollLeft = ta.scrollLeft;
}
$('netlist').addEventListener('input', hl);
$('netlist').addEventListener('scroll', () => {
  $('hl').scrollTop = $('netlist').scrollTop;
  $('hl').scrollLeft = $('netlist').scrollLeft;
});
hl();

const EXAMPLES = __EXAMPLES__;
const exSel = $('examples');
EXAMPLES.forEach(name => {
  const o = document.createElement('option');
  o.value = name; o.textContent = name;
  exSel.appendChild(o);
});

$('load').onclick = () => {
  fetch('/example?name=' + encodeURIComponent(exSel.value))
    .then(r => { if (!r.ok) throw new Error('erro ' + r.status); return r.text(); })
    .then(t => { $('netlist').value = t; hl(); setStatus('Exemplo carregado: ' + exSel.value, 'ok'); })
    .catch(e => setStatus('Falha ao carregar: ' + e.message, 'err'));
};
$('openBtn').onclick = () => $('openFile').click();
$('openFile').onchange = e => {
  const f = e.target.files[0]; if (!f) return;
  f.text().then(t => { $('netlist').value = t; hl();
    setStatus('Arquivo carregado: ' + f.name, 'ok'); });
};
$('save').onclick = () => {
  const blob = new Blob([$('netlist').value], {type: 'text/plain'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = (exSel.value || 'netlist.cir');
  a.click(); URL.revokeObjectURL(a.href);
};

function setStatus(msg, cls){ $('status').textContent = msg; $('status').className = cls || ''; }

$('run').onclick = async () => {
  $('run').disabled = true; setStatus('Simulando…');
  const t0 = performance.now();
  try {
    const r = await fetch('/run', {method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({netlist: $('netlist').value})});
    const data = await r.json();
    if (!data.ok) { $('results').innerHTML = ''; setStatus('erro: ' + data.error, 'err'); }
    else {
      render(data);
      setStatus('Simulação concluída em ' + ((performance.now() - t0) / 1000).toFixed(2) + ' s.', 'ok');
    }
  } catch (e) { setStatus('erro: ' + e.message, 'err'); }
  $('run').disabled = false;
};

function render(data){
  const box = $('results'); box.innerHTML = '';
  window.__msState = {data: data, group: (data.cases && data.cases.length)
      ? data.cases[0].label : null, run: 0};
  for (const w of data.warnings)
    box.insertAdjacentHTML('beforeend', `<div class="card warn">⚠ ${esc(w)}</div>`);
  if (data.schematic) {
    const sch = document.createElement('div'); sch.className = 'card';
    sch.innerHTML = `<h3>Circuito</h3>` +
      `<img class="plot" src="/schematic?run=${data.run}">` +
      `<div class="note"><a href="/schematic?run=${data.run}" download="minispice_esquema.png">Baixar esquema PNG</a></div>`;
    box.appendChild(sch);
  }
  if (!data.results.length)
    box.insertAdjacentHTML('beforeend',
      '<div class="card warn">Netlist sem análise (.op / .dc / .tran).</div>');
  if ((data.cases && data.cases.length > 1) || data.mc) {
    const bar = document.createElement('div'); bar.className = 'card casebar';
    const opts = (data.cases || []).map((c, i) =>
      `<option value="${esc(c.label)}">${esc(c.label)}</option>`).join('');
    bar.innerHTML = `<b>Caso</b> <select id="caseSel">${opts}</select>` +
      (data.mc ? ` <span class="note">Monte Carlo · ${data.mc.runs} execuções por caso</span>` : '');
    box.appendChild(bar);
    bar.querySelector('#caseSel').onchange = e => {
      window.__msState.group = e.target.value;
      window.__msState.run = 0;
      renderBody();
    };
  }
  const body = document.createElement('div');
  body.id = 'resultbody'; body.className = 'resultbody';
  box.appendChild(body);
  renderBody();
}

function renderBody(){
  const st = window.__msState, data = st.data;
  const box = $('resultbody'); box.innerHTML = '';
  const group = (data.cases || []).find(c => c.label === st.group);
  const idxs = group ? group.indices : data.results.map((_, i) => i);
  if (data.mc) {
    const mc = data.mc;
    const stats = (group && mc.by_case[group.label]) || mc.measures || {};
    for (const name of Object.keys(stats)) {
      const s = stats[name];
      if (s.mean === null || s.mean === undefined) continue;
      const card = document.createElement('div'); card.className = 'card histcard';
      let yieldTxt = '';
      if (s.yield !== null && s.yield !== undefined) {
        yieldTxt = `<div class="yield">yield ${(s.yield * 100).toFixed(1)}% ` +
          `(${s.pass}/${s.total})</div>`;
      }
      const lim = [];
      if (s.limit_min !== null) lim.push('min ' + fmt(s.limit_min));
      if (s.limit_max !== null) lim.push('max ' + fmt(s.limit_max));
      card.innerHTML = `<h3>${esc(name)}</h3>${yieldTxt}` +
        `<img class="plot" src="/hist?run=${data.run}&measure=${encodeURIComponent(name)}` +
        `&case=${encodeURIComponent(group ? group.label : '')}">` +
        `<div class="note">média ${fmt(s.mean)} · σ ${fmt(s.std)} · ` +
        `${fmt(s.min)} … ${fmt(s.max)}` + (lim.length ? ' · limites ' + lim.join(', ') : '') + `</div>`;
      box.appendChild(card);
    }
    if (idxs.length > 1) {
      const pick = document.createElement('div'); pick.className = 'card casebar';
      const opts = idxs.map((i, k) =>
        `<option value="${k}">execução ${k + 1}</option>`).join('');
      pick.innerHTML = `<b>Execução</b> <select id="runSel">${opts}</select>`;
      box.appendChild(pick);
      pick.querySelector('#runSel').onchange = e => {
        st.run = +e.target.value; renderBody();
      };
      const chosen = idxs[Math.min(st.run, idxs.length - 1)];
      box.appendChild(resultCard(data.results[chosen], data.run));
    } else {
      idxs.forEach(i => box.appendChild(resultCard(data.results[i], data.run)));
    }
  } else {
    idxs.forEach(i => box.appendChild(resultCard(data.results[i], data.run)));
  }
}

function resultCard(r, run){
  const card = document.createElement('div'); card.className = 'card';
  if (r.type === 'op') {
    card.innerHTML = `<h3>${esc(r.title)}</h3>`;
    const table = document.createElement('table');
    r.rows.forEach(([k, v]) => {
      const tr = document.createElement('tr');
      tr.innerHTML = `<td>${esc(k)}</td><td>${esc(v)}</td>`;
      table.appendChild(tr);
    });
    card.appendChild(table);
  } else {
    const chosen = new Set();
    card.innerHTML = `<h3>${esc(r.title)}</h3><span class="sigbox"></span>`;
    const picker = card.querySelector('.sigbox');
    r.signals.forEach(s => {
      const cb = document.createElement('input');
      cb.type = 'checkbox'; cb.checked = s.sel;
      cb.onchange = () => {
        if (cb.checked) chosen.add(s.label); else chosen.delete(s.label);
        draw(); drawTable();
      };
      if (s.sel) chosen.add(s.label);
      const lb = document.createElement('label');
      lb.appendChild(cb); lb.appendChild(document.createTextNode(s.label));
      picker.appendChild(lb);
    });
    if (!chosen.size) { r.signals.slice(0, 2).forEach(s => chosen.add(s.label)); }

    // zoom controls (analysis x-axis: time or sweep source)
    const zoom = document.createElement('div'); zoom.className = 'zoomrow';
    zoom.innerHTML = `zoom: [ <input id="x0${r.i}" placeholder="de" value="">` +
      ` … <input id="x1${r.i}" placeholder="até" value=""> ] ${esc(r.head[0])}` +
      ` <button class="secondary" style="padding:3px 10px;font-size:12px">Ajustar</button>`;
    const x0 = zoom.querySelector(`#x0${r.i}`), x1 = zoom.querySelector(`#x1${r.i}`);
    const zoomBtn = zoom.querySelector('button');
    card.appendChild(zoom);

    const img = document.createElement('img'); img.className = 'plot';
    const aSave = document.createElement('a'); aSave.textContent = 'Baixar PNG';
    aSave.download = 'minispice_grafico.png';

    const q = () => {
      let u = '/plot?run=' + run + '&i=' + r.i + '&sigs=' + encodeURIComponent([...chosen].join(','));
      if (x0.value && !isNaN(+x0.value)) u += '&x0=' + (+x0.value);
      if (x1.value && !isNaN(+x1.value)) u += '&x1=' + (+x1.value);
      return u;
    };

    // ---- interactive cursor: canvas overlay on top of the plot image ----
    // axes fractions must match webui._AX_* constants exactly
    const CYCLE = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                   '#9467bd', '#8c564b', '#e377c2', '#7f7f7f'];
    const AXL = 0.075, AXR = 0.985, AXB = 0.115, AXT = 0.955;
    const fmtSci = v => (+v).toPrecision(4).replace(/\\.?0+e/, 'e').replace(/\\.?0+$/, '');

    const wrap = document.createElement('div'); wrap.className = 'plotwrap';
    const cv = document.createElement('canvas'); cv.className = 'cursorcv';
    const readout = document.createElement('div'); readout.className = 'cursorread';
    wrap.appendChild(img); wrap.appendChild(cv); wrap.appendChild(readout);
    card.appendChild(wrap);

    let cursorData = null, cursorAxis = null, cursorDisp = null, cursorX = -1;
    function cursorFetch(){
      fetch('/dat?run=' + run + '&i=' + r.i)
        .then(resp => resp.json())
        .then(d => { cursorData = d.series; cursorAxis = d.x; cursorDisp = d.disp || null; overlay(); })
        .catch(() => {});
    }
    function overlay(){
      if (!cv.width || cv.width !== img.clientWidth) {
        cv.width = img.clientWidth; cv.height = img.clientHeight;
      }
      const ctx = cv.getContext('2d'); ctx.clearRect(0, 0, cv.width, cv.height);
      if (!cursorData || !cursorAxis || cursorX < 0 || !img.complete || !img.naturalWidth) return;
      const axL = AXL * cv.width, axR = AXR * cv.width;
      const axT = (1 - AXT) * cv.height, axB = (1 - AXB) * cv.height;
      const uFrac = Math.min(1, Math.max(0, (cursorX - axL) / (axR - axL)));
      const xs = cursorAxis;
      const xTarget = xs[0] + uFrac * (xs[xs.length - 1] - xs[0]);
      let lo = 0, hi = xs.length - 1;            // nearest sample (binary search)
      while (hi - lo > 1) { const mid = (lo + hi) >> 1; xs[mid] < xTarget ? lo = mid : hi = mid; }
      const frac = (xTarget - xs[lo]) / ((xs[hi] - xs[lo]) || 1);
      const px = axL + uFrac * (axR - axL);
      ctx.strokeStyle = '#e05561'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(px, axT); ctx.lineTo(px, axB); ctx.stroke();

      const sel = [...chosen];
      let ymin = Infinity, ymax = -Infinity;
      sel.forEach(l => { const a = cursorData[l]; if (!a) return;
        for (const v of a) { if (v < ymin) ymin = v; if (v > ymax) ymax = v; }
      });
      const nearest = frac > 0.5 ? hi : lo;
      const axisTxt = cursorDisp ? cursorDisp[nearest] : fmtSci(xTarget);
      const lines = [`${r.head[0]} = ${axisTxt}`];
      sel.forEach((label, i) => {
        const a = cursorData[label]; if (!a) return;
        const v = a[lo] + frac * (a[hi] - a[lo]);
        const yFrac = (v - ymin) / ((ymax - ymin) || 1);
        const py = axB - yFrac * (axB - axT);
        ctx.fillStyle = CYCLE[sel.indexOf(label) % CYCLE.length];
        ctx.beginPath(); ctx.arc(px, py, 3.5, 0, 2 * Math.PI); ctx.fill();
        lines.push(`${label} = ${fmtSci(v)}`);
      });
      readout.innerHTML = lines.map(l => `<div>${esc(l)}</div>`).join('');
    }
    cv.addEventListener('mousemove', ev => {
      cursorX = ev.clientX - cv.getBoundingClientRect().left; overlay();
    });
    cv.addEventListener('mouseleave', () => { cursorX = -1; overlay(); readout.innerHTML = ''; });
    window.addEventListener('resize', () => overlay());

    zoomBtn.onclick = () => draw();
    const draw = () => {
      img.src = q() + '&t=' + Date.now();
      aSave.href = q() + '&save=1';
    };
    img.addEventListener('load', () => { cv.width = img.clientWidth; overlay(); });
    img.addEventListener('load', cursorFetch);
    draw();

    const tbl = document.createElement('table');
    function drawTable(){
      tbl.innerHTML = '<tr>' + r.head.map(h => '<th>' + esc(h) + '</th>').join('') + '</tr>' +
        r.rows.map(row => '<tr>' + row.map((v, j) => {
          if (j > 0 && !chosen.has(r.head[j])) return '';
          return '<td>' + esc(v) + '</td>';
        }).join('') + '</tr>').join('');
    }
    drawTable();
    card.appendChild(tbl);
    const aSaveLink = document.createElement('div'); aSaveLink.className = 'note';
    aSaveLink.appendChild(aSave);
    card.appendChild(aSaveLink);
  }
  const info = document.createElement('div'); info.className = 'note';
  info.textContent = r.info + ' · ';
  const a = document.createElement('a');
  a.href = '/csv?run=' + run + '&i=' + r.i; a.download = r.i + '.csv';
  a.textContent = 'Baixar CSV';
  info.appendChild(a);
  card.appendChild(info);
  return card;
}
</script>
</body>
</html>"""


# ------------------------------------------------------------------ rendering


def _fmt(value: float) -> str:
    """Compact human-readable number for display/CSV heads."""
    text = f"{value:.6g}"
    return "0" if text in {"-0", "0"} else text


def _fmt_freq(value: float) -> str:
    """Frequency with an engineering suffix (1.2k, 3MEG...)."""
    text = f"{value:.4g}"
    if abs(value) >= 1e6:
        text = f"{value / 1e6:.4g}MEG"
    elif abs(value) >= 1e3:
        text = f"{value / 1e3:.4g}k"
    return f"{text}Hz"


def _signals(result: DcSweep | Transient | AcResult | NoiseResult) -> list[tuple[str, str]]:
    """(label, element-name) pairs for plots: V:<node> / I:<branch>."""
    if isinstance(result, NoiseResult):
        labels = [("N:out", "out")]
        if result.input_psd is not None:
            labels.append(("N:in", "in"))
        return labels
    return [
        *(("V:" + n, n) for n in result.nodes),
        *(("I:" + n, n) for n in result.currents),
    ]


def _axis_name(result: DcSweep | Transient | AcResult | NoiseResult) -> str:
    if isinstance(result, (AcResult, NoiseResult)):
        return "freq (Hz)"
    return result.source if isinstance(result, DcSweep) else "time"


def _measure_info(result: object) -> str:
    measurements = getattr(result, "measurements", {})
    if not measurements:
        return ""
    parts = [
        f"{name} = {'failed' if not np.isfinite(value) else _fmt(value)}"
        for name, value in measurements.items()
    ]
    return " · ".join(parts)


def _case_suffix(result: object) -> str:
    case = getattr(result, "case", {})
    if not case:
        return ""
    return " [" + ", ".join(f"{name}={value:g}" for name, value in case.items()) + "]"


def _subtitle(result: DcSweep | Transient | AcResult | NoiseResult) -> str:
    if isinstance(result, NoiseResult):
        text = f"integrated output noise: {result.total_rms:.4g} V rms"
        if result.input_psd is not None:
            text += (
                f" · input: {result.input_rms:.4g} {result.input_unit} rms"
            )
        return text
    if isinstance(result, AcResult):
        return (
            f"{len(result.freqs)} points · f = "
            f"{_fmt_freq(result.freqs[0])} .. {_fmt_freq(result.freqs[-1])}"
        )
    if isinstance(result, DcSweep):
        return f"{len(result.sweep)} points · sweep of {result.source}"
    return f"{len(result.time)} points · t = {_fmt(result.time[0])} .. {_fmt(result.time[-1])} s"


# Fixed axes placement so the browser can map pixels <-> data exactly for
# the interactive cursor overlay (fractions of the rendered figure).
_AX_LEFT, _AX_RIGHT = 0.075, 0.985
_AX_BOTTOM, _AX_TOP = 0.115, 0.955


def _plot(
    result: DcSweep | Transient | AcResult | NoiseResult,
    labels: list[str],
    path: str | io.BytesIO,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if isinstance(result, NoiseResult):
        fig, ax = plt.subplots(figsize=(8, 4.2), dpi=115)
        fig.patch.set_facecolor("#0d0f14")
        ax.set_facecolor("#0d0f14")
        density = np.sqrt(np.maximum(result.output_psd, 0.0))
        ax.loglog(result.freqs, np.where(density > 0, density, np.nan), lw=1.6, label="N:out")
        if result.input_psd is not None:
            input_density = np.sqrt(np.maximum(result.input_psd, 0.0))
            ax.loglog(
                result.freqs,
                np.where(input_density > 0, input_density, np.nan),
                lw=1.2,
                label=f"N:in ({result.input_unit})",
            )
        ax.set_xlabel("freq (Hz)", fontsize=10)
        ax.set_ylabel("V/sqrt(Hz)", fontsize=10)
        ax.grid(True, which="both", color="#262a38", lw=0.7)
        for spine in ax.spines.values():
            spine.set_color("#3a3f52")
        ax.tick_params(colors="#8a91a3", labelsize=9)
        fig.subplots_adjust(left=_AX_LEFT, right=_AX_RIGHT, top=_AX_TOP, bottom=_AX_BOTTOM)
        if isinstance(path, str):
            fig.savefig(path, facecolor="#0d0f14")
        else:
            fig.savefig(path, format="png", facecolor="#0d0f14")
            path.seek(0)
        plt.close(fig)
        return

    pairs = _signals(result)
    names = {label: element for label, element in pairs}
    unknown = [label for label in labels if label not in names]
    if unknown:
        raise KeyError(f"unknown signal(s): {', '.join(unknown)}")

    if isinstance(result, AcResult):
        x = result.freqs
        ylabel = "Magnitude (dB)"
    elif isinstance(result, DcSweep):
        x = result.sweep
        ylabel = "V / A"
    else:
        x = result.time
        ylabel = "V / A"
    fig, ax = plt.subplots(figsize=(8, 4.2), dpi=115)
    fig.patch.set_facecolor("#0d0f14")
    ax.set_facecolor("#0d0f14")
    for label in labels:
        if isinstance(result, AcResult):
            series = 20 * np.log10(np.maximum(np.abs(result.nodes[names[label]]), 1e-30))
        elif label.startswith("V:"):
            series = result.nodes[names[label]]
        else:
            series = result.currents[names[label]]
        ax.plot(x, series, lw=1.6, label=label)
    for spine in ax.spines.values():
        spine.set_color("#3a3f52")
    ax.tick_params(colors="#8a91a3", labelsize=9)
    ax.grid(True, which="both", color="#262a38", lw=0.7)
    ax.set_xlabel(_axis_name(result), fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.margins(x=0, y=0)  # xlim/ylim == data range: exact pixel <-> data cursor
    if isinstance(result, AcResult):
        ax.set_xscale("log")
        ax.set_xlim(x[0], x[-1])
    if labels:
        legend = ax.legend(facecolor="#1c1f27", edgecolor="#3a3f52", fontsize=9)
        for text in legend.get_texts():
            text.set_color("#d8dce6")
    fig.subplots_adjust(left=_AX_LEFT, right=_AX_RIGHT, top=_AX_TOP, bottom=_AX_BOTTOM)
    if isinstance(path, str):
        fig.savefig(path, facecolor="#0d0f14")
    else:
        fig.savefig(path, format="png", facecolor="#0d0f14")
        path.seek(0)
    plt.close(fig)


def _db_array(result: AcResult, name: str) -> np.ndarray:
    """dB magnitude series, clamped so log(0) never poisons the payload."""
    return 20 * np.log10(np.maximum(np.abs(result.nodes[name]), 1e-30))


def _dat_payload(result: DcSweep | Transient | AcResult | NoiseResult) -> dict:
    """Full arrays for the browser-side cursor (JSON-encodable).

    For AC the x axis is log10(f) (matches the logarithmic plot axes so the
    linear pixel mapping stays exact) and every series is already in dB;
    `disp` carries formatted frequency strings for the readout.
    """
    if isinstance(result, NoiseResult):
        series = {"N:out": np.sqrt(np.maximum(result.output_psd, 0.0)).tolist()}
        if result.input_psd is not None:
            series["N:in"] = np.sqrt(np.maximum(result.input_psd, 0.0)).tolist()
        return {
            "axis": _axis_name(result),
            "x": np.log10(result.freqs).tolist(),
            "disp": [_fmt_freq(f) for f in result.freqs],
            "series": series,
        }
    if isinstance(result, AcResult):
        return {
            "axis": _axis_name(result),
            "x": np.log10(result.freqs).tolist(),
            "disp": [_fmt_freq(f) for f in result.freqs],
            "series": {
                label: (
                    _db_array(result, name)
                    if label.startswith("V:")
                    else 20 * np.log10(np.maximum(np.abs(result.currents[name]), 1e-30))
                ).tolist()
                for label, name in _signals(result)
            },
        }
    x = result.sweep if isinstance(result, DcSweep) else result.time
    return {
        "axis": _axis_name(result),
        "x": x.tolist(),
        "series": {
            label: (
                result.nodes[name] if label.startswith("V:") else result.currents[name]
            ).tolist()
            for label, name in _signals(result)
        },
    }


def _histogram_png(
    values: list[float], name: str, stats: dict
) -> bytes:
    """Histogram of one measured value across Monte Carlo runs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    fig, ax = plt.subplots(figsize=(4.8, 2.8), dpi=115)
    fig.patch.set_facecolor("#0d0f14")
    ax.set_facecolor("#0d0f14")
    bins = min(30, max(5, int(np.sqrt(max(data.size, 1))) + 2))
    ax.hist(data, bins=bins, color="#4f8cff", edgecolor="#0d0f14", alpha=0.9)
    low = stats.get("limit_min")
    high = stats.get("limit_max")
    if low is not None:
        ax.axvline(float(low), color="#e0a855", lw=1.2, ls="--", label="limit")
    if high is not None:
        ax.axvline(float(high), color="#e0a855", lw=1.2, ls="--")
    mean = stats.get("mean")
    if mean is not None:
        ax.axvline(float(mean), color="#3fb96f", lw=1.2, label="mean")
    ax.set_title(name, color="#d8dce6", fontsize=10)
    ax.set_xlabel("value", color="#8a91a3", fontsize=9)
    ax.set_ylabel("runs", color="#8a91a3", fontsize=9)
    ax.grid(True, color="#262a38", lw=0.6, axis="y")
    for spine in ax.spines.values():
        spine.set_color("#3a3f52")
    ax.tick_params(colors="#8a91a3", labelsize=8)
    if low is not None or high is not None or mean is not None:
        legend = ax.legend(fontsize=7, frameon=False, labelcolor="#8a91a3")
        if legend is not None:
            legend.set_bbox_to_anchor((1.0, 1.0))
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor="#0d0f14")
    plt.close(fig)
    return buf.getvalue()


def _csv_bytes(
    result: OperatingPoint | DcSweep | Transient | AcResult | NoiseResult,
) -> bytes:
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    if isinstance(result, OperatingPoint):
        writer.writerow(["signal", "value"])
        for name, value in result.nodes.items():
            writer.writerow([f"V({name})", repr(value)])
        for name, value in result.currents.items():
            writer.writerow([f"I({name})", repr(value)])
        for device, info in result.details.items():
            writer.writerow([f"{device}.region", str(info.get("region", ""))])
            for key, value in info.items():
                if isinstance(value, (int, float)):
                    writer.writerow([f"{device}.{key}", repr(float(value))])
        for name, value in result.measurements.items():
            writer.writerow([f"meas {name}", repr(value)])
    elif isinstance(result, NoiseResult):
        has_input = result.input_psd is not None
        header = ["freq", "output_psd", "output_density"]
        if has_input:
            header += ["input_psd", "input_density"]
        writer.writerow(header)
        for index, (freq, psd) in enumerate(zip(result.freqs, result.output_psd)):
            row = [
                repr(float(freq)),
                repr(float(psd)),
                repr(float(np.sqrt(max(psd, 0.0)))),
            ]
            if has_input and result.input_psd is not None:
                in_psd = float(result.input_psd[index])
                row += [repr(in_psd), repr(float(np.sqrt(max(in_psd, 0.0))))]
            writer.writerow(row)
    elif isinstance(result, AcResult):
        pairs = _signals(result)
        writer.writerow(
            ["freq", *(f"{lbl} dB" for lbl, _ in pairs), *(f"{lbl} deg" for lbl, _ in pairs)]
        )
        phasors: dict[str, list[complex]] = {}
        phasors.update({f"V:{n}": result.nodes[n] for n in result.nodes})
        phasors.update({f"I:{n}": result.currents[n] for n in result.currents})
        for i, freq in enumerate(result.freqs):
            mags = [20 * np.log10(max(abs(phasors[lbl][i]), 1e-30)) for lbl, _ in pairs]
            phases = [np.rad2deg(np.angle(phasors[lbl][i])) for lbl, _ in pairs]
            writer.writerow([repr(float(freq)), *(repr(float(v)) for v in mags),
                             *(repr(float(v)) for v in phases)])
    else:
        x = result.sweep if isinstance(result, DcSweep) else result.time
        writer.writerow(
            [_axis_name(result), *(f"V({n})" for n in result.nodes), *(f"I({n})" for n in result.currents)]
        )
        for i, xv in enumerate(x):
            writer.writerow([
                repr(float(xv)),
                *(repr(float(series[i])) for series in result.nodes.values()),
                *(repr(float(series[i])) for series in result.currents.values()),
            ])
    return buf.getvalue().encode()


def _preview_indices(count: int) -> list[int | None]:
    if count > _ROW_PREVIEW:
        return [*range(_ROW_PREVIEW // 2), None, *range(count - _ROW_PREVIEW // 2, count)]
    return list(range(count))


def _ac_rows(
    result: AcResult, pairs: list[tuple[str, str]], j: int | None
) -> list[str]:
    if j is None:
        return ["⋮"] * (1 + len(pairs))
    phasors: dict[str, np.ndarray] = {}
    phasors.update({f"V:{n}": result.nodes[n] for n in result.nodes})
    phasors.update({f"I:{n}": result.currents[n] for n in result.currents})
    values = [
        _fmt(20 * np.log10(max(abs(phasors[label][j]), 1e-30)))
        for label, _ in pairs
    ]
    return [_fmt_freq(result.freqs[j]), *values]


def _result_payload(
    index: int, result: OperatingPoint | DcSweep | Transient | AcResult
) -> dict:
    """JSON description of one result for the web UI."""
    if isinstance(result, OperatingPoint):
        rows = [
            *(["V(" + n + ")", _fmt(v)] for n, v in result.nodes.items()),
            *(["I(" + n + ")", _fmt(v)] for n, v in result.currents.items()),
        ]
        for device, info in result.details.items():
            rows.append([f"{device}.region", str(info.get("region", ""))])
            for key in ("id", "vd", "vgs", "vds", "vbs", "vth", "vov", "gm", "gds", "gmb", "gd"):
                if key in info:
                    rows.append([f"{device}.{key}", _fmt(float(info[key]))])
        info_text = f"Newton iterations: {result.iterations}"
        measured = _measure_info(result)
        if measured:
            info_text += " · " + measured
        return {
            "type": "op",
            "i": index,
            "title": "Operating point" + _case_suffix(result),
            "head": [],
            "rows": rows,
            "info": info_text,
        }
    if isinstance(result, NoiseResult):
        density = np.sqrt(np.maximum(result.output_psd, 0.0))
        if result.input_psd is not None:
            input_density = np.sqrt(np.maximum(result.input_psd, 0.0))
            rows = [
                [_fmt_freq(f), _fmt(float(v)), _fmt(float(i))]
                for f, v, i in zip(result.freqs, density, input_density)
            ]
            head = [_axis_name(result), "N:out", "N:in"]
            signals = [
                {"label": "N:out", "sel": True},
                {"label": "N:in", "sel": True},
            ]
        else:
            rows = [
                [_fmt_freq(f), _fmt(float(v))]
                for f, v in zip(result.freqs, density)
            ]
            head = [_axis_name(result), "N:out"]
            signals = [{"label": "N:out", "sel": True}]
        return {
            "type": "series",
            "i": index,
            "title": "Noise analysis" + _case_suffix(result),
            "head": head,
            "rows": rows,
            "signals": signals,
            "info": _subtitle(result),
        }
    if isinstance(result, AcResult):
        x = result.freqs
        pairs = _signals(result)
        head = [_axis_name(result), *(label for label, _ in pairs)]
        signals = [{"label": label, "sel": i < 3} for i, (label, _) in enumerate(pairs)]
        rows = [_ac_rows(result, pairs, j) for j in _preview_indices(len(x))]
        return {
            "type": "series",
            "i": index,
            "title": "AC analysis" + _case_suffix(result),
            "head": head,
            "rows": rows,
            "signals": signals,
            "info": _subtitle(result)
            + ((" · " + _measure_info(result)) if result.measurements else ""),
        }
    x = result.sweep if isinstance(result, DcSweep) else result.time
    pairs = _signals(result)
    head = [_axis_name(result), *(label for label, _ in pairs)]
    signals = [
        {"label": label, "sel": i < 3}
        for i, (label, _) in enumerate(pairs)
    ]

    rows: list[list[str]] = []
    for j in _preview_indices(len(x)):
        if j is None:
            rows.append(["⋮"] * len(head))
            continue
        values = [
            _fmt(result.nodes[name][j]) if label.startswith("V:") else _fmt(result.currents[name][j])
            for label, name in pairs
        ]
        rows.append([_fmt(x[j]), *values])
    kind = "DC sweep of " + result.source if isinstance(result, DcSweep) else "Transient"
    return {
        "type": "series",
        "i": index,
        "title": kind + _case_suffix(result),
        "head": head,
        "rows": rows,
        "signals": signals,
        "info": _subtitle(result)
        + ((" · " + _measure_info(result)) if result.measurements else ""),
    }


def _slice(
    result: DcSweep | Transient | AcResult | NoiseResult, qs: dict[str, str]
) -> DcSweep | Transient | AcResult | NoiseResult:
    """Return a copy of the series limited to x0 <= x <= x1 (zoom)."""
    try:
        x0, x1 = float(qs.get("x0", "")), float(qs.get("x1", ""))
    except ValueError:
        return result
    if isinstance(result, NoiseResult):
        x = result.freqs
        mask = (x >= min(x0, x1)) & (x <= max(x0, x1))
        if not mask.any():
            return result
        return NoiseResult(
            freqs=x[mask],
            output_psd=result.output_psd[mask],
            output_name=result.output_name,
            input_source=result.input_source,
            input_unit=result.input_unit,
            total_rms=result.total_rms,
            input_psd=None if result.input_psd is None else result.input_psd[mask],
            input_rms=result.input_rms,
            measurements=dict(result.measurements),
            case=dict(result.case),
        )
    if isinstance(result, AcResult):
        x = result.freqs
        mask = (x >= min(x0, x1)) & (x <= max(x0, x1))
        if not mask.any():
            return result
        nodes = {name: series[mask] for name, series in result.nodes.items()}
        currents = {name: series[mask] for name, series in result.currents.items()}
        return AcResult(
            x[mask], nodes, currents, measurements=dict(result.measurements)
        )
    x = result.sweep if isinstance(result, DcSweep) else result.time
    mask = (x >= min(x0, x1)) & (x <= max(x0, x1))
    if not mask.any():
        return result
    nodes = {name: series[mask] for name, series in result.nodes.items()}
    currents = {name: series[mask] for name, series in result.currents.items()}
    finite = np.asarray(x)[mask]
    if isinstance(result, DcSweep):
        return DcSweep(
            result.source,
            finite,
            nodes,
            currents,
            measurements=dict(result.measurements),
        )
    return Transient(
        finite, nodes, currents, measurements=dict(result.measurements)
    )


# -------------------------------------------------------------------- server


class _Store:
    """Keeps the results of the last runs in memory, keyed by run id."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._runs: dict[int, dict] = {}

    def insert(
        self,
        results: list,
        circuit,
        schematic: bytes | None,
        cases: list[dict] | None = None,
        mc: dict | None = None,
    ) -> int:
        with self._lock:
            run_id = (max(self._runs) + 1) if self._runs else 1
            self._runs[run_id] = {
                "results": results,
                "circuit": circuit,
                "schematic": schematic,
                "cases": cases or [],
                "mc": mc,
            }
            return run_id

    def get(self, run_id: int) -> dict:
        with self._lock:
            if run_id not in self._runs:
                raise KeyError(run_id)
            return self._runs[run_id]


_STORE = _Store()


class _Handler(BaseHTTPRequestHandler):
    server_version = "MiniSPICE/0.1"
    protocol_version = "HTTP/1.1"

    def log_message(self, *args) -> None:  # silence the request log
        pass

    def _send(self, code: int, body: bytes, ctype: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj: dict, code: int = 200) -> None:
        self._send(code, json.dumps(obj).encode(), "application/json")

    def _qs(self) -> dict[str, str]:
        return {k: v[0] for k, v in parse_qs(urlparse(self.path).query).items()}

    def _route(self) -> None:
        url = urlparse(self.path)
        try:
            if self.command == "GET":
                routes = {
                    "/": lambda: self._page(),
                    "/example": lambda: self._example(self._qs().get("name", "")),
                    "/plot": lambda: self._plot(self._qs()),
                    "/csv": lambda: self._csv(self._qs()),
                    "/schematic": lambda: self._schematic(self._qs()),
                    "/dat": lambda: self._dat(self._qs()),
                    "/hist": lambda: self._hist(self._qs()),
                }
            else:
                routes = {"/run": self._run}
            handler = routes.get(url.path)
            if handler is None:
                self._json({"ok": False, "error": "not found"}, 404)
            else:
                handler()
        except (NetlistError, ConvergenceError, ValueError) as exc:
            self._json({"ok": False, "error": str(exc)})
        except KeyError as exc:
            self._json({"ok": False, "error": f"invalid request: {exc}"}, 400)

    do_GET = _route
    do_POST = _route

    # ---- GET /
    def _page(self) -> None:
        default = _examples_dir() / "rc_lowpass.cir"
        text = default.read_text(encoding="utf-8") if default.exists() else "* netlist\n.end\n"
        examples = sorted(p.name for p in _examples_dir().glob("*.cir"))
        page = _PAGE.replace("__DEFAULT__", html.escape(text, quote=False))
        page = page.replace("__EXAMPLES__", json.dumps(examples))
        self._send(200, page.encode(), "text/html; charset=utf-8")

    # ---- GET /example?name=...
    def _example(self, name: str) -> None:
        path = _examples_dir() / name
        if ".." in Path(name).parts or not path.is_file():
            self._json({"ok": False, "error": f"unknown example: {name}"}, 404)
            return
        self._send(200, path.read_bytes(), "text/plain; charset=utf-8")

    # ---- POST /run  {"netlist": "..."}
    def _run(self) -> None:
        length = int(self.headers.get("Content-Length", "0"))
        body = json.loads(self.rfile.read(length) or b"{}")
        text = body.get("netlist")
        if not isinstance(text, str) or not text.strip():
            self._json({"ok": False, "error": "netlist vazia"})
            return

        circuit = NetlistParser(text=text).parse()
        simulation = Simulation(circuit)
        results = simulation.run()

        schematic: bytes | None = None
        try:
            from .schematic import schematic_png

            buf = io.BytesIO()
            if schematic_png(circuit, buf):
                schematic = buf.getvalue()
        except Exception:
            schematic = None  # the schematic is a bonus, never breaks the run

        mc = None
        if circuit.mc_runs > 0:
            mc = {
                "runs": circuit.mc_runs,
                "measures": simulation.yield_summary,
                "by_case": simulation.yield_by_case,
            }
        run_id = _STORE.insert(
            results, circuit, schematic, simulation.case_groups, mc
        )
        payload = {
            "ok": True,
            "run": run_id,
            "title": circuit.title,
            "warnings": circuit.warnings,
            "schematic": schematic is not None,
            "cases": simulation.case_groups,
            "mc": mc,
            "results": [_result_payload(i, r) for i, r in enumerate(results)],
        }
        self._json(payload)

    # ---- GET /schematic?run=
    def _schematic(self, qs: dict[str, str]) -> None:
        run = _STORE.get(int(qs["run"]))
        if not run["schematic"]:
            self._json({"ok": False, "error": "no schematic available"}, 404)
            return
        self._send(200, run["schematic"], "image/png")

    # ---- GET /plot?run=&i=&sigs=V:out,I:dd[&x0=&x1=&save=1]
    def _plot(self, qs: dict[str, str]) -> None:
        run = _STORE.get(int(qs["run"]))
        result = run["results"][int(qs.get("i", "0"))]
        if not isinstance(result, (DcSweep, Transient, AcResult, NoiseResult)):
            self._json({"ok": False, "error": "no series to plot"}, 400)
            return
        result = _slice(result, qs)
        labels = [s for s in qs.get("sigs", "").split(",") if s]
        buf = io.BytesIO()
        _plot(result, labels, buf)
        self.send_response(200)
        self.send_header("Content-Type", "image/png")
        self.send_header("Content-Length", str(len(buf.getvalue())))
        if qs.get("save"):
            self.send_header(
                "Content-Disposition", f'attachment; filename="minispice_plot_{qs["run"]}_{qs.get("i", "0")}.png"'
            )
        self.end_headers()
        self.wfile.write(buf.getvalue())

    # ---- GET /csv?run=&i=
    def _csv(self, qs: dict[str, str]) -> None:
        run = _STORE.get(int(qs["run"]))
        result = run["results"][int(qs.get("i", "0"))]
        self._send(200, _csv_bytes(result), "text/csv; charset=utf-8")

    # ---- GET /dat?run=&i=   (full arrays for the browser cursor)
    def _dat(self, qs: dict[str, str]) -> None:
        run = _STORE.get(int(qs["run"]))
        result = run["results"][int(qs.get("i", "0"))]
        if not isinstance(result, (DcSweep, Transient, AcResult, NoiseResult)):
            self._json({"ok": False, "error": "no series data"}, 400)
            return
        self._json(_dat_payload(result))

    # ---- GET /hist?run=&measure=&case=
    def _hist(self, qs: dict[str, str]) -> None:
        run = _STORE.get(int(qs["run"]))
        measure = qs.get("measure", "")
        label = qs.get("case") or None
        indices = list(range(len(run["results"])))
        for group in run.get("cases", []):
            if group["label"] == label:
                indices = list(group["indices"])
                break
        stats_by_case = run.get("mc") or {}
        stats = (
            stats_by_case.get("by_case", {}).get(label)
            if label
            else stats_by_case.get("measures")
        ) or {}
        values = [
            run["results"][i].measurements.get(measure, np.nan)
            for i in indices
            if measure in getattr(run["results"][i], "measurements", {})
        ]
        if not values:
            self._json({"ok": False, "error": f"no .meas {measure!r}"}, 404)
            return
        png = _histogram_png(values, measure, stats.get(measure, {}))
        self._send(200, png, "image/png")


def create_server(port: int = 0) -> ThreadingHTTPServer:
    """Create (not start) the GUI server on the given or a free localhost port."""
    return ThreadingHTTPServer(("127.0.0.1", port), _Handler)


def _open_quiet(url: str) -> tuple[bool, str]:
    """Open the URL in the system browser without leaking its console noise
    (Chromium logs zygote errors into the terminal otherwise).

    Returns (opened, hint) where hint explains a failure.
    """
    for opener in ("xdg-open", "gio"):
        binary = shutil.which(opener)
        if not binary:
            continue
        try:
            subprocess.Popen(
                [binary, url],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                stdin=subprocess.DEVNULL,
                start_new_session=True,
            )
            return True, ""
        except OSError:
            continue
    try:
        import webbrowser

        return bool(webbrowser.open(url)), ""
    except Exception as exc:
        return False, str(exc)


DEFAULT_PORT = 8080


def _parse_gui_args(argv: list[str] | None) -> tuple[int, list[str]]:
    """Split GUI argv into (port, netlist args). Accepts --port N / --port=N."""
    port = DEFAULT_PORT
    remaining: list[str] = []
    tokens = list(argv or [])
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token in ("--no-browser", "-n"):
            remaining.append(token)
        elif token == "--port" and i + 1 < len(tokens):
            port = int(tokens[i + 1])
            i += 1
        elif token.startswith("--port="):
            port = int(token.split("=", 1)[1])
        else:
            remaining.append(token)
        i += 1
    return port, remaining


def main(argv: list[str] | None = None) -> int:
    """Start the GUI server on localhost and open the system browser.

    argv accepts an optional ``--port N`` (default 8080, falls back to a
    random free port if busy) and ``--no-browser``.
    """
    port, argv = _parse_gui_args(argv)
    want_browser = "--no-browser" not in argv

    try:
        server = create_server(port)
    except OSError:  # default port busy — fall back to a random one
        server = create_server(0)
    host, bound_port = server.server_address[:2]
    url = f"http://{host}:{bound_port}/"
    if bound_port != port:
        print(f"porta {port} ocupada; usando {bound_port}")
    print(f"MiniSPICE GUI: {url}")
    print("deixe este terminal aberto — a GUI roda neste processo (Ctrl+C para sair).")
    if argv and argv[0]:
        print(f"netlist inicial: {argv[0]}")

    def launch() -> None:
        # sanity-check: prove the server answers before opening a tab
        try:
            import urllib.request

            urllib.request.urlopen(url, timeout=5).read(64)
            alive = True
        except Exception:
            alive = False
        opened, _hint = _open_quiet(url) if alive and want_browser else (None, "")
        if not opened:
            note = "" if alive else "  (aviso: o servidor não respondeu ao auto-teste!)"
            print(f"acesse manualmente no navegador: {url}{note}")

    threading.Timer(0.5, launch).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
