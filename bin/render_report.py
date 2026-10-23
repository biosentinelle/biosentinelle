#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
render_report.py -- BioSentinelle strain-identification HTML report renderer.

Reads one results JSON file and writes one self-contained interactive HTML
report. Standard library only (json, argparse, html, datetime); no templating
engine, no external packages, no CDN assets. The report is meant to be opened
from a file:// URL on an offline HPC cluster, so every byte of CSS, JS and data
lives inside the single output file.

Usage:
    python render_report.py --json results.json --out report.html [--title "..."]

PERFORMANCE / LAZY RENDERING STRATEGY
-------------------------------------
A panel may contain ~50000 SNP sites per strain. Building that many <tr>
elements eagerly (either in Python or at page load in JS) makes the file huge
and freezes the browser. Instead:

  * The whole results object is embedded ONCE as JSON text inside a
    <script type="application/json"> block and parsed with JSON.parse on the
    first interaction that needs it. No per-site markup is emitted by Python.
  * Each strain's per-site table is built LAZILY: nothing is created until the
    user expands that strain's detail panel.
  * Within a detail panel the table is PAGINATED (fixed page size, default 200
    rows). Only the rows of the current page exist in the DOM at any time; the
    previous page's rows are dropped before the new ones are appended. Filtering
    (all / mismatches only / uncovered only) recomputes an index array of
    integers (not DOM nodes) and resets to page 1. So the DOM cost is O(page
    size) and the JS cost per interaction is O(n_sites) over a plain array,
    which stays responsive at 50k rows.
  * Rows are built with document.createElement + textContent, never innerHTML,
    so no data-derived string is ever parsed as markup. Strings interpolated
    into markup on the Python side go through html.escape.
"""

import argparse
import datetime
import html
import json
import sys


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def esc(value):
    """html.escape() on the str() of any value, safe for None."""
    if value is None:
        return ""
    return html.escape(str(value), quote=True)


def num(value, default="n/a"):
    """Format a number for display, tolerating None / missing values."""
    if value is None:
        return default
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return "{:,}".format(value).replace(",", " ")
    if isinstance(value, float):
        return "{:.3f}".format(value)
    return str(value)


def pct(value, digits=1):
    """Fraction (0-1) as a percentage string."""
    if value is None:
        return "n/a"
    try:
        return "{:.{d}f}%".format(float(value) * 100.0, d=digits)
    except (TypeError, ValueError):
        return "n/a"


def clamp01(value):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return 0.0
    return max(0.0, min(1.0, v))


def colour_of(tool):
    """Normalise the colour status of a tool element.

    `colour` is taken from the JSON when present; otherwise it is derived from
    `score` using the documented thresholds (green >= 99, orange >= 90, red
    below). Anything unrecognised falls back to 'red' so nothing is silently
    reported as present.
    """
    c = tool.get("colour")
    if isinstance(c, str) and c.lower() in ("green", "orange", "red"):
        return c.lower()
    score = tool.get("score")
    try:
        score = float(score)
    except (TypeError, ValueError):
        return "red"
    if score >= 99.0:
        return "green"
    if score >= 90.0:
        return "orange"
    return "red"


# Status glyphs: colour is never the only cue -- each status also carries a
# distinct shape/glyph and a spelled-out text label.
COLOUR_GLYPH = {"green": "●", "orange": "▲", "red": "■"}
COLOUR_LABEL = {"green": "present", "orange": "partial", "red": "absent"}


def embed_json(obj):
    """Serialise data for a <script type="application/json"> block.

    JSON.parse(textContent) needs the text to survive HTML parsing: the only
    sequences that can terminate the block early are '<' (as in '</script>')
    and, for safety against legacy comment parsing, '>' and '&'. They are
    escaped as \\uXXXX, which is valid JSON and parses back to the same string.
    """
    text = json.dumps(obj, ensure_ascii=False, separators=(",", ":"))
    return (text.replace("&", "\\u0026")
                .replace("<", "\\u003c")
                .replace(">", "\\u003e")
                .replace(" ", "\\u2028")
                .replace(" ", "\\u2029"))


# --------------------------------------------------------------------------- #
# CSS
# --------------------------------------------------------------------------- #

CSS = """
:root {
  color-scheme: light dark;
  --bg: #ffffff;
  --bg-alt: #f4f6f8;
  --bg-raised: #ffffff;
  --fg: #16191d;
  --fg-muted: #5a6572;
  --border: #d2d8de;
  --border-strong: #9aa4b0;
  --accent: #1f5fa9;
  --accent-fg: #ffffff;
  --bar-bg: #e2e7ec;
  --bar-fill: #1f5fa9;
  --green: #1a7f37;
  --green-bg: #e6f4ea;
  --orange: #9a5b00;
  --orange-bg: #fdf1dd;
  --red: #a4232b;
  --red-bg: #fbe9ea;
  --warn-bg: #fdf1dd;
  --warn-border: #c88a1e;
  --alert-bg: #fbe9ea;
  --alert-border: #a4232b;
  --mono: ui-monospace, "Cascadia Mono", "Consolas", "DejaVu Sans Mono", monospace;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #14171b;
    --bg-alt: #1c2126;
    --bg-raised: #1f252b;
    --fg: #e7ecf1;
    --fg-muted: #a3aeba;
    --border: #39424c;
    --border-strong: #5c6773;
    --accent: #6aa9ef;
    --accent-fg: #0d1014;
    --bar-bg: #2b333b;
    --bar-fill: #6aa9ef;
    --green: #67d58a;
    --green-bg: #15301f;
    --orange: #e3ab51;
    --orange-bg: #33260f;
    --red: #f2868c;
    --red-bg: #361619;
    --warn-bg: #33260f;
    --warn-border: #e3ab51;
    --alert-bg: #361619;
    --alert-border: #f2868c;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0;
  padding: 0 0 4rem 0;
  background: var(--bg);
  color: var(--fg);
  font-family: system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  font-size: 15px;
  line-height: 1.5;
}
.wrap { max-width: 1180px; margin: 0 auto; padding: 0 16px; }
header.page {
  background: var(--bg-alt);
  border-bottom: 1px solid var(--border);
  padding: 20px 0 16px 0;
  margin-bottom: 20px;
}
header.page h1 { margin: 0 0 4px 0; font-size: 1.5rem; }
header.page .sample { font-family: var(--mono); font-size: 1.15rem; font-weight: 700; }
dl.meta {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
  gap: 6px 20px;
  margin: 12px 0 0 0;
}
dl.meta div { border-left: 3px solid var(--border-strong); padding-left: 8px; }
dl.meta dt { font-size: 0.78rem; text-transform: uppercase; letter-spacing: .04em; color: var(--fg-muted); }
dl.meta dd { margin: 0; font-weight: 600; word-break: break-word; }
section { margin: 28px 0; }
h2 { font-size: 1.2rem; border-bottom: 2px solid var(--border); padding-bottom: 4px; margin-bottom: 6px; }
h3 { font-size: 1.02rem; margin: 18px 0 6px 0; }
p.note { color: var(--fg-muted); font-size: 0.87rem; margin: 6px 0 12px 0; }
.banner {
  border: 2px solid var(--alert-border);
  background: var(--alert-bg);
  border-radius: 6px;
  padding: 12px 16px;
  margin: 0 0 18px 0;
}
.banner h2 { border: none; margin: 0 0 4px 0; color: var(--red); font-size: 1.12rem; }
.banner .glyph { font-weight: 700; margin-right: 6px; }
.warnings {
  border: 1px solid var(--warn-border);
  background: var(--warn-bg);
  border-radius: 6px;
  padding: 10px 16px 10px 16px;
  margin: 0 0 18px 0;
}
.warnings h2 { border: none; font-size: 1rem; margin: 0 0 4px 0; }
.warnings ul { margin: 0; padding-left: 22px; }
table { border-collapse: collapse; width: 100%; font-size: 0.92rem; }
th, td { text-align: left; padding: 5px 8px; border-bottom: 1px solid var(--border); vertical-align: top; }
thead th {
  background: var(--bg-alt);
  border-bottom: 2px solid var(--border-strong);
  position: sticky; top: 0; z-index: 1;
  font-size: 0.8rem; text-transform: uppercase; letter-spacing: .03em;
}
td.n, th.n { text-align: right; font-variant-numeric: tabular-nums; font-family: var(--mono); }
tbody tr:nth-child(even) { background: var(--bg-alt); }
.strain { border: 1px solid var(--border); border-radius: 6px; margin-bottom: 10px; background: var(--bg-raised); }
.strain.top { border-color: var(--accent); border-width: 2px; }
button.strain-head {
  display: block; width: 100%; text-align: left; cursor: pointer;
  background: transparent; border: 0; color: inherit; font: inherit;
  padding: 10px 12px;
}
button.strain-head:hover { background: var(--bg-alt); }
button.strain-head:focus-visible { outline: 3px solid var(--accent); outline-offset: -3px; }
.strain-row { display: flex; flex-wrap: wrap; align-items: baseline; gap: 10px; }
.rank {
  font-family: var(--mono); font-weight: 700;
  background: var(--accent); color: var(--accent-fg);
  border-radius: 4px; padding: 1px 8px; font-size: 0.85rem;
}
.sid { font-family: var(--mono); font-weight: 700; font-size: 1.05rem; }
.prob { font-variant-numeric: tabular-nums; font-weight: 700; min-width: 5.5em; }
.bar { flex: 1 1 180px; min-width: 140px; height: 12px; background: var(--bar-bg); border: 1px solid var(--border-strong); border-radius: 3px; overflow: hidden; }
.bar > i { display: block; height: 100%; background: var(--bar-fill); }
.chev { font-family: var(--mono); color: var(--fg-muted); }
.counts { margin-top: 6px; font-size: 0.85rem; color: var(--fg-muted); display: flex; flex-wrap: wrap; gap: 14px; }
.counts b { color: var(--fg); font-variant-numeric: tabular-nums; }
.detail { border-top: 1px solid var(--border); padding: 10px 12px; }
.detail[hidden] { display: none; }
.controls { display: flex; flex-wrap: wrap; gap: 12px; align-items: center; margin: 8px 0; }
.controls label { font-size: 0.88rem; }
select, .pager button {
  font: inherit; font-size: 0.88rem; padding: 3px 8px;
  background: var(--bg-raised); color: var(--fg);
  border: 1px solid var(--border-strong); border-radius: 4px;
}
.pager button { cursor: pointer; }
.pager button:disabled { opacity: .45; cursor: default; }
.pager { display: flex; gap: 8px; align-items: center; font-size: 0.88rem; }
.tablebox { max-height: 460px; overflow: auto; border: 1px solid var(--border); border-radius: 4px; }
.status { white-space: nowrap; font-weight: 600; }
.status .g { font-family: var(--mono); margin-right: 4px; }
.status.green { color: var(--green); }
.status.orange { color: var(--orange); }
.status.red { color: var(--red); }
.pill {
  display: inline-block; border-radius: 10px; padding: 0 8px;
  font-size: 0.8rem; font-weight: 700; border: 1px solid currentColor;
}
.pill.green { color: var(--green); background: var(--green-bg); }
.pill.orange { color: var(--orange); background: var(--orange-bg); }
.pill.red { color: var(--red); background: var(--red-bg); }
.tiles { display: flex; flex-wrap: wrap; gap: 12px; margin: 10px 0; }
.tile {
  border: 1px solid var(--border-strong); border-radius: 6px;
  padding: 8px 16px; min-width: 150px; background: var(--bg-raised);
}
.tile .v { font-size: 1.5rem; font-weight: 700; font-variant-numeric: tabular-nums; }
.tile .k { font-size: 0.82rem; color: var(--fg-muted); text-transform: uppercase; letter-spacing: .03em; }
.legend { border: 1px solid var(--border); border-radius: 6px; padding: 8px 14px; background: var(--bg-alt); font-size: 0.88rem; }
.legend ul { margin: 4px 0 0 0; padding-left: 20px; }
.short-note { font-size: 0.87rem; color: var(--fg-muted); margin: 4px 0 8px 0; }
.empty { color: var(--fg-muted); font-style: italic; padding: 6px 0; }
code { font-family: var(--mono); }
footer.page { margin-top: 36px; border-top: 1px solid var(--border); padding-top: 10px; color: var(--fg-muted); font-size: 0.82rem; }
"""


# --------------------------------------------------------------------------- #
# JavaScript
# --------------------------------------------------------------------------- #

JS = r"""
'use strict';
// The payload is parsed once, on first use. Per-site tables are built lazily
// (on expand) and paginated, so the DOM never holds more than PAGE_SIZE rows
// per strain even when a panel carries tens of thousands of sites.
var PAGE_SIZE = 200;
var DATA = null;

function data() {
  if (DATA === null) {
    var el = document.getElementById('report-data');
    DATA = JSON.parse(el.textContent);
  }
  return DATA;
}

var OBS_LABEL = {
  'R': 'R (reference allele)',
  'A': 'A (alternative allele)',
  'H': 'H (both alleles)',
  '?': '? (not covered)'
};

function obsLabel(v) {
  if (v === null || v === undefined || v === '') { return '-'; }
  var s = String(v);
  return OBS_LABEL[s] || s;
}

function fmtNum(v, digits) {
  if (v === null || v === undefined || v === '') { return '-'; }
  var n = Number(v);
  if (!isFinite(n)) { return String(v); }
  if (digits === 0) { return String(Math.round(n)); }
  return n.toFixed(digits);
}

function txt(tag, value, cls) {
  var el = document.createElement(tag);
  el.textContent = (value === null || value === undefined) ? '-' : String(value);
  if (cls) { el.className = cls; }
  return el;
}

// A site is "uncovered" when observed is '?' or depth is 0/absent.
function isUncovered(s) {
  if (s.observed === '?') { return true; }
  var d = s.depth;
  return (d === null || d === undefined || Number(d) === 0);
}

function isMismatch(s) {
  if (s.match === true) { return false; }
  if (s.match === false) { return true; }
  // No explicit flag: fall back to comparing observed with expected.
  if (isUncovered(s)) { return false; }
  return String(s.observed) !== String(s.expected);
}

function buildIndex(sites, mode) {
  var idx = [], i;
  for (i = 0; i < sites.length; i++) {
    var s = sites[i];
    if (mode === 'mismatch') {
      if (!isUncovered(s) && isMismatch(s)) { idx.push(i); }
    } else if (mode === 'uncovered') {
      if (isUncovered(s)) { idx.push(i); }
    } else {
      idx.push(i);
    }
  }
  return idx;
}

var STATE = {}; // strainIndex -> {idx: [], page: 0, mode: 'all', built: false}

function renderPage(si) {
  var st = STATE[si];
  var sites = data().strains[si].sites || [];
  var tbody = document.getElementById('sbody-' + si);
  var total = st.idx.length;
  var pages = Math.max(1, Math.ceil(total / PAGE_SIZE));
  if (st.page >= pages) { st.page = pages - 1; }
  if (st.page < 0) { st.page = 0; }

  // Drop the previous page's rows entirely, then append only this page's.
  while (tbody.firstChild) { tbody.removeChild(tbody.firstChild); }

  var start = st.page * PAGE_SIZE;
  var end = Math.min(total, start + PAGE_SIZE);
  var frag = document.createDocumentFragment();
  for (var j = start; j < end; j++) {
    var s = sites[st.idx[j]];
    var tr = document.createElement('tr');
    tr.appendChild(txt('td', s.site_id));
    tr.appendChild(txt('td', s.chrom));
    tr.appendChild(txt('td', s.pos, 'n'));
    tr.appendChild(txt('td', s.ref));
    tr.appendChild(txt('td', s.alt));
    tr.appendChild(txt('td', s.expected));
    tr.appendChild(txt('td', obsLabel(s.observed)));
    tr.appendChild(txt('td', fmtNum(s.depth, 0), 'n'));
    tr.appendChild(txt('td', fmtNum(s.alt_fraction, 3), 'n'));
    tr.appendChild(txt('td', fmtNum(s.weight, 2), 'n'));
    var state = isUncovered(s) ? 'not covered'
              : (isMismatch(s) ? 'mismatch' : 'match');
    var glyph = isUncovered(s) ? '□' : (isMismatch(s) ? '■' : '●');
    var td = document.createElement('td');
    td.className = 'status ' + (isUncovered(s) ? 'red' : (isMismatch(s) ? 'orange' : 'green'));
    var g = document.createElement('span');
    g.className = 'g';
    g.textContent = glyph;
    td.appendChild(g);
    td.appendChild(document.createTextNode(state));
    tr.appendChild(td);
    frag.appendChild(tr);
  }
  tbody.appendChild(frag);

  var info = document.getElementById('pinfo-' + si);
  if (total === 0) {
    info.textContent = 'no site matches this filter';
  } else {
    info.textContent = 'rows ' + (start + 1) + '–' + end + ' of ' + total +
                       '  (page ' + (st.page + 1) + ' of ' + pages + ')';
  }
  document.getElementById('pprev-' + si).disabled = (st.page <= 0);
  document.getElementById('pnext-' + si).disabled = (st.page >= pages - 1);
}

function setFilter(si, mode) {
  var st = STATE[si];
  st.mode = mode;
  st.idx = buildIndex(data().strains[si].sites || [], mode);
  st.page = 0;
  renderPage(si);
}

function toggleStrain(si) {
  var panel = document.getElementById('detail-' + si);
  var btn = document.getElementById('head-' + si);
  var open = btn.getAttribute('aria-expanded') === 'true';
  if (open) {
    panel.hidden = true;
    btn.setAttribute('aria-expanded', 'false');
    document.getElementById('chev-' + si).textContent = '▶';
    return;
  }
  panel.hidden = false;
  btn.setAttribute('aria-expanded', 'true');
  document.getElementById('chev-' + si).textContent = '▼';
  if (!STATE[si]) { STATE[si] = {mode: 'all', page: 0, idx: [], built: false}; }
  if (!STATE[si].built) {   // lazy: first expansion builds the table
    STATE[si].built = true;
    setFilter(si, 'all');
  }
}

function filterTools(cls) {
  var rows = document.querySelectorAll('tr[data-eclass]');
  for (var i = 0; i < rows.length; i++) {
    var keep = (cls === '__all__') || (rows[i].getAttribute('data-eclass') === cls);
    rows[i].hidden = !keep;
  }
}

document.addEventListener('DOMContentLoaded', function () {
  var heads = document.querySelectorAll('button.strain-head[data-si]');
  for (var i = 0; i < heads.length; i++) {
    (function (btn) {
      var si = Number(btn.getAttribute('data-si'));
      if (btn.getAttribute('data-sites') === '0') { return; }
      btn.addEventListener('click', function () { toggleStrain(si); });
    })(heads[i]);
  }
  var sels = document.querySelectorAll('select[data-filter-si]');
  for (var j = 0; j < sels.length; j++) {
    (function (sel) {
      var si = Number(sel.getAttribute('data-filter-si'));
      sel.addEventListener('change', function () { setFilter(si, sel.value); });
    })(sels[j]);
  }
  var btns = document.querySelectorAll('button[data-page-si]');
  for (var k = 0; k < btns.length; k++) {
    (function (b) {
      var si = Number(b.getAttribute('data-page-si'));
      var step = Number(b.getAttribute('data-step'));
      b.addEventListener('click', function () { STATE[si].page += step; renderPage(si); });
    })(btns[k]);
  }
  var tsel = document.getElementById('tool-class-filter');
  if (tsel) { tsel.addEventListener('change', function () { filterTools(tsel.value); }); }
});
"""


# --------------------------------------------------------------------------- #
# HTML sections
# --------------------------------------------------------------------------- #

def render_header(d, title):
    panel = d.get("panel") or {}
    rows = [
        ("Organism", d.get("organism")),
        ("Reference", d.get("reference")),
        ("Generated", d.get("generated")),
        ("Panel sites", num(panel.get("n_sites"))),
        ("k-mer size (k)", num(panel.get("k"))),
        ("Panel strains", ", ".join(str(s) for s in (panel.get("strains") or [])) or "n/a"),
    ]
    items = "\n".join(
        "      <div><dt>{}</dt><dd>{}</dd></div>".format(esc(k), esc(v) if v not in (None, "") else "n/a")
        for k, v in rows
    )
    return """<header class="page">
  <div class="wrap">
    <h1>{title}</h1>
    <div>Sample <span class="sample">{sample}</span></div>
    <dl class="meta">
{items}
    </dl>
  </div>
</header>""".format(title=esc(title), sample=esc(d.get("sample") or "unknown"), items=items)


def render_unassigned(d):
    if not d.get("unassigned"):
        return ""
    reason = d.get("unassigned_reason") or "No reason recorded."
    return """<div class="banner" role="alert">
  <h2><span class="glyph">&#9888;</span>UNASSIGNED &mdash; no strain in the panel matches this sample well enough</h2>
  <p>The sample could not be assigned to any strain of the panel. The ranking below is still
  shown for inspection, but the posterior probabilities are conditional on the panel and must
  <strong>not</strong> be read as an identification.</p>
  <p><strong>Reason:</strong> {reason}</p>
</div>""".format(reason=esc(reason))


def render_warnings(d):
    warns = d.get("warnings") or []
    if not warns:
        return ""
    items = "\n".join("    <li>{}</li>".format(esc(w)) for w in warns)
    return """<div class="warnings" role="note">
  <h2>Warnings ({n})</h2>
  <ul>
{items}
  </ul>
</div>""".format(n=len(warns), items=items)


def render_strains(d):
    strains = d.get("strains") or []
    if not strains:
        return """<section id="strains"><h2>Strain ranking</h2>
  <p class="empty">No strain was scored for this sample.</p></section>"""

    ordered = sorted(
        enumerate(strains),
        key=lambda t: (t[1].get("rank") if isinstance(t[1].get("rank"), (int, float)) else 10 ** 9, t[0]),
    )

    blocks = []
    for si, s in ordered:
        sites = s.get("sites") or []
        n_sites = len(sites)
        prob = clamp01(s.get("probability"))
        top = " top" if s.get("rank") == 1 else ""
        disabled = ' disabled aria-disabled="true"' if n_sites == 0 else ""
        chev = "▶" if n_sites else "—"

        head = """  <h3 style="margin:0">
    <button class="strain-head" id="head-{si}" type="button" data-si="{si}" data-sites="{n}"
            aria-expanded="false" aria-controls="detail-{si}"{dis}>
      <span class="strain-row">
        <span class="chev" id="chev-{si}" aria-hidden="true">{chev}</span>
        <span class="rank">rank {rank}</span>
        <span class="sid">{sid}</span>
        <span class="prob">{p}</span>
        <span class="bar" role="img" aria-label="posterior probability {p}"><i style="width:{w:.2f}%"></i></span>
        <span class="prob">concordance {conc}</span>
      </span>
      <span class="counts">
        <span>matching <b>{m}</b></span>
        <span>mismatching <b>{mm}</b></span>
        <span>uncovered <b>{un}</b></span>
        <span>heterozygous <b>{het}</b></span>
        <span>called <b>{cal}</b></span>
        <span>log-likelihood <b>{ll}</b></span>
        <span>{drill}</span>
      </span>
    </button>
  </h3>""".format(
            si=si, n=n_sites, dis=disabled, chev=chev,
            rank=esc(s.get("rank")), sid=esc(s.get("strain_id") or "?"),
            p=pct(s.get("probability")), w=prob * 100.0,
            conc=pct(s.get("concordance"), 1),
            m=num(s.get("sites_matching")), mm=num(s.get("sites_mismatching")),
            un=num(s.get("sites_uncovered")), het=num(s.get("sites_het")),
            cal=num(s.get("sites_called")), ll=num(s.get("log_likelihood")),
            drill=("per-site detail: {} rows (click to expand)".format(num(n_sites))
                   if n_sites else "per-site detail not available in this result file"),
        )

        if n_sites:
            detail = """  <div class="detail" id="detail-{si}" hidden>
    <p class="note"><strong>{sid}:</strong> {m} sites match the expectation, {mm} mismatch,
      {un} are not covered, {het} are heterozygous (out of {cal} called sites).</p>
    <div class="controls">
      <label for="filter-{si}">Show sites:</label>
      <select id="filter-{si}" data-filter-si="{si}">
        <option value="all">all sites</option>
        <option value="mismatch">mismatches only</option>
        <option value="uncovered">uncovered only</option>
      </select>
      <span class="pager">
        <button type="button" id="pprev-{si}" data-page-si="{si}" data-step="-1">&#8592; previous</button>
        <button type="button" id="pnext-{si}" data-page-si="{si}" data-step="1">next &#8594;</button>
        <span id="pinfo-{si}"></span>
      </span>
    </div>
    <div class="tablebox">
      <table>
        <thead><tr>
          <th>site_id</th><th>chrom</th><th class="n">pos</th><th>ref</th><th>alt</th>
          <th>expected</th><th>observed</th><th class="n">depth</th>
          <th class="n">alt_fraction</th><th class="n">weight</th><th>state</th>
        </tr></thead>
        <tbody id="sbody-{si}"></tbody>
      </table>
    </div>
  </div>""".format(
                si=si, sid=esc(s.get("strain_id") or "?"),
                m=num(s.get("sites_matching")), mm=num(s.get("sites_mismatching")),
                un=num(s.get("sites_uncovered")), het=num(s.get("sites_het")),
                cal=num(s.get("sites_called")),
            )
        else:
            detail = ""

        blocks.append('<div class="strain{top}">\n{head}\n{detail}\n</div>'.format(
            top=top, head=head, detail=detail))

    return """<section id="strains">
  <h2>Strain ranking</h2>
  <p class="note">The probability is a <strong>posterior over the strains present in the panel only</strong>.
  It is normalised across those strains and therefore always sums to 1: it does <em>not</em> say whether the
  true strain of the sample is in the panel at all. That question is answered separately by the
  <em>unassigned</em> flag shown as a banner at the top of this report when it is raised.</p>
  <p class="note">Click a strain row to expand its per-site detail table.</p>
{blocks}
</section>""".format(blocks="\n".join(blocks))


def tool_row(t):
    c = colour_of(t)
    return """      <tr data-eclass="{ec}">
        <td><code>{name}</code></td>
        <td>{ecls}</td>
        <td class="n">{ln}</td>
        <td class="n">{breadth}</td>
        <td class="n">{ident}</td>
        <td class="n">{depth}</td>
        <td class="n">{score}</td>
        <td class="status {c}"><span class="g" aria-hidden="true">{glyph}</span>{label}</td>
        <td><code>{eid}</code></td>
      </tr>""".format(
        ec=esc(t.get("element_class") or "unclassified"),
        name=esc(t.get("name") or "?"),
        ecls=esc(t.get("element_class") or "unclassified"),
        ln=num(t.get("length_bp")),
        breadth=pct(t.get("breadth"), 1),
        ident=pct(t.get("identity"), 1),
        depth=num(t.get("mean_depth")),
        score=num(t.get("score")),
        c=c, glyph=COLOUR_GLYPH[c], label=COLOUR_LABEL[c],
        eid=esc(t.get("element_id") or ""),
    )


TOOL_HEAD = """      <tr><th>name</th><th>class</th><th class="n">length (bp)</th>
        <th class="n">breadth</th><th class="n">identity</th><th class="n">mean depth</th>
        <th class="n">score</th><th>status</th><th>element_id</th></tr>"""


def render_tools(d):
    tools = d.get("tools") or []
    summary = d.get("tool_summary") or {}

    tiles = "\n".join(
        """    <div class="tile"><div class="k"><span class="pill {c}">{glyph} {label}</span></div>
      <div class="v">{v}</div></div>""".format(
            c=c, glyph=COLOUR_GLYPH[c], label=COLOUR_LABEL[c], v=num(summary.get(c, 0)))
        for c in ("green", "orange", "red")
    )

    legend = """  <div class="legend">
    <strong>Colour / status legend</strong> &mdash; the status is never signalled by hue alone:
    each level also carries a glyph and a written label.
    <ul>
      <li><span class="status green"><span class="g">&#9679;</span>present</span>
        &mdash; round glyph, green: score &ge; 99 (element present)</li>
      <li><span class="status orange"><span class="g">&#9650;</span>partial</span>
        &mdash; triangular glyph, orange: 90 &le; score &lt; 99 (partial or diverged)</li>
      <li><span class="status red"><span class="g">&#9632;</span>absent</span>
        &mdash; square glyph, red: score &lt; 90 (absent or too little evidence)</li>
    </ul>
  </div>"""

    if not tools:
        return """<section id="tools">
  <h2>Tool elements</h2>
  <div class="tiles">
{tiles}
  </div>
{legend}
  <p class="empty">No tool element was reported for this sample.</p>
</section>""".format(tiles=tiles, legend=legend)

    scorable = [t for t in tools if (t.get("scorable") or "yes") != "short"]
    short = [t for t in tools if (t.get("scorable") or "yes") == "short"]

    classes = sorted({str(t.get("element_class") or "unclassified") for t in scorable})
    options = "\n".join(
        '      <option value="{v}">{v}</option>'.format(v=esc(c)) for c in classes)

    main_tbl = """  <table>
    <thead>
{head}
    </thead>
    <tbody>
{rows}
    </tbody>
  </table>""".format(head=TOOL_HEAD,
                     rows="\n".join(tool_row(t) for t in sorted(
                         scorable,
                         key=lambda t: (str(t.get("element_class") or "~"),
                                        -(t.get("score") if isinstance(t.get("score"), (int, float)) else -1)))))

    if short:
        short_tbl = """  <h3>Short elements &mdash; shown for information only ({n})</h3>
  <p class="short-note">Elements shorter than 100 bp (<code>scorable = "short"</code>) cannot be scored
  reliably by breadth of coverage; they are listed separately and are <strong>not</strong> included in the
  green / orange / red counts above.</p>
  <table>
    <thead>
{head}
    </thead>
    <tbody>
{rows}
    </tbody>
  </table>""".format(n=len(short), head=TOOL_HEAD,
                     rows="\n".join(tool_row(t) for t in short))
    else:
        short_tbl = ""

    return """<section id="tools">
  <h2>Tool elements</h2>
  <p class="note">{ns} scorable element(s) and {nsh} short element(s) screened.</p>
  <div class="tiles">
{tiles}
  </div>
{legend}
  <div class="controls">
    <label for="tool-class-filter">Filter by element class:</label>
    <select id="tool-class-filter">
      <option value="__all__">all classes</option>
{options}
    </select>
  </div>
{main_tbl}
{short_tbl}
</section>""".format(ns=len(scorable), nsh=len(short), tiles=tiles, legend=legend,
                     options=options, main_tbl=main_tbl, short_tbl=short_tbl)


# --------------------------------------------------------------------------- #
# document
# --------------------------------------------------------------------------- #

def render(d, title):
    now = datetime.datetime.now().replace(microsecond=0).isoformat()
    parts = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        "<title>{}</title>".format(esc(title)),
        "<style>{}</style>".format(CSS),
        "</head>",
        "<body>",
        render_header(d, title),
        '<div class="wrap">',
        render_unassigned(d),
        render_warnings(d),
        render_strains(d),
        render_tools(d),
        '<footer class="page">BioSentinelle report &mdash; rendered {} by '
        "<code>render_report.py</code>. Self-contained: no network access required."
        "</footer>".format(esc(now)),
        "</div>",
        '<script type="application/json" id="report-data">{}</script>'.format(embed_json(d)),
        "<script>{}</script>".format(JS),
        "</body>",
        "</html>",
        "",
    ]
    return "\n".join(parts)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Render a self-contained interactive HTML report from a "
                    "BioSentinelle strain-identification results JSON file.")
    ap.add_argument("--json", required=True, metavar="FILE", help="input results JSON file")
    ap.add_argument("--out", required=True, metavar="FILE", help="output HTML file")
    ap.add_argument("--title", default=None, help="report title (default: BioSentinelle report -- <sample>)")
    args = ap.parse_args(argv)

    with open(args.json, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        sys.stderr.write("error: the JSON root must be an object\n")
        return 1

    title = args.title or "BioSentinelle report — {}".format(data.get("sample") or "unknown sample")

    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(render(data, title))

    sys.stderr.write("wrote {}\n".format(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
