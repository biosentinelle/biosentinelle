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
import re
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
# help popups
# --------------------------------------------------------------------------- #

# Text of the '?' popups, keyed by data-help. Authored here, never derived
# from the data, so the JS may insert it as markup.
HELP = {
    "organism": "Species of the strain database the sample is compared with.",
    "reference": "Reference genome the panel positions (CHROM, POS) and the REF alleles refer "
                 "to. For yeast: S288C, assembly R64-1-1.",
    "date": "Date and time the sample was scored against the panel.",
    "panel_sites": "Number of genome positions (SNVs and small indels) at which the panel strains "
                   "differ from one another. At each one the sample is genotyped by counting the "
                   "short k-mers that carry the REF allele and those that carry the ALT allele.",
    "panel_strains": "Number of known strains (laboratory and wild, from public assemblies and "
                     "read runs) whose genotype at the panel sites is known. The sample is "
                     "compared with every one of them; the 10 best are listed below.",
    "mating": "Measured on the sample reads, not taken from the strain: mean count of the k-mers "
              "specific to the MATa cassette against those specific to MAT&alpha;. Every genome "
              "carries one silent a (HMR) and one silent &alpha; (HML) copy, plus the active MAT "
              "copy, so the a fraction a / (a + &alpha;) is about 0.67 in a MATa haploid, 0.33 in "
              "a MAT&alpha; haploid and 0.5 in a MATa/&alpha; diploid (or in a mix of MATa and "
              "MAT&alpha; cells). Calls: MATa &ge; 0.58, MATalpha &le; 0.42, MATa/alpha in "
              "between.",
    "probability": "<b>Posterior probability</b> of this strain, among the panel strains only. "
                   "Each strain gets a log-likelihood (LL, see its own ?). The probability is "
                   "exp(LL &minus; LL<sub>best</sub>), divided by the sum over all panel strains "
                   "(equal prior for every strain), so the probabilities add up to 100 %.<br>"
                   "It is <b>not</b> the probability that the sample <i>is</i> this strain, nor a "
                   "measure of fit: one extra mismatch divides the odds by about 100, so a strain "
                   "differing by a few sites only already falls to ~0 %. How well the sample fits "
                   "is the concordance. Strains only a few sites away from the best one are listed in "
                   "the Decision, with how the mating type and the native markers settle them.",
    "concordance": "<b>Concordance = matching / called.</b> Fraction of the called sites where the "
                   "sample carries the allele expected for this strain. A sample that really is "
                   "this strain gives close to 100 % (only sequencing errors disagree). If the "
                   "best strain is below 90 % (min_concordance), the sample is reported as "
                   "<i>unassigned</i>: this is the closest strain, but not this strain.",
    "matching": "Called sites where the sample's genotype equals the strain's: REF where the "
                "strain carries REF, ALT where it carries ALT, H where it is heterozygous.",
    "mismatching": "Called sites where they differ: ALT where the strain carries REF or the "
                   "reverse, H where the strain is homozygous, or a homozygous call where the "
                   "strain is heterozygous (that last kind costs less in the log-likelihood, a "
                   "spore or a loss of heterozygosity gives it).",
    "uncovered": "Panel sites not used for this strain: the sample has fewer than 3 k-mer "
                 "observations there (kmer_min_depth, not sequenced deeply enough), or the "
                 "strain's own genotype is unknown there (gap in its assembly). "
                 "Called + uncovered = all panel sites.",
    "heterozygous": "Sites the sample itself calls H: both alleles seen, the ALT fraction being "
                    "between 20 % and 80 %. The same number for every strain. Many of them mean a "
                    "diploid, a cross or a mixed culture.",
    "called": "Sites where the sample has a call (REF, ALT or H) and the strain's genotype is "
              "known. Called = matching + mismatching.",
    "log_likelihood": "<b>LL = matching &times; ln(1 &minus; &epsilon;) + mismatching &times; "
                      "ln(&epsilon;)</b>, with &epsilon; = 0.01 the per-site error rate "
                      "(strain_error_rate). A homozygous call where the strain is heterozygous "
                      "counts ln(0.25) instead of ln(&epsilon;). At sites where the strain's "
                      "genotype is unknown, the sample's call is scored with its frequency among "
                      "the strain's clade (or the whole panel).<br>Closer to 0 is better; strains "
                      "are ranked by it. Each extra mismatch costs about 4.6.",
    "chrom": "Chromosome of the site, named after the reference genome: chrI to chrXVI, chrM for "
             "the mitochondrial genome. The RefSeq accession (NC_0011xx) of the reference "
             "sequence is shown when the mouse is on the name.",
    "expected": "Genotype of this strain at the site in the panel: REF (reference allele), ALT "
                "(alternative allele), H (both, heterozygous) or ? (unknown for this strain).",
    "observed": "Genotype called in the sample from the k-mer counts: ALT when the ALT fraction "
                "is &ge; 80 %, REF when it is &le; 20 %, H in between, ? when fewer than 3 k-mers "
                "were seen (not covered).",
    "depth": "Number of sample k-mers seen at the site, REF and ALT together.",
    "alt_fraction": "ALT k-mers / (REF + ALT k-mers) at the site: 0 = only REF, 1 = only ALT.",
    "weight": "Weight of the site in the score. All sites currently weigh 1.",
    "class": "Function of the element, from its name: selection marker (drug resistance or "
             "auxotrophy), reporter, epitope tag, degron, recombination site, promoter, "
             "terminator.",
    "length": "Length of the element's sequence in the database. Elements under 100 bp are listed "
              "apart: too short to be scored by coverage.",
    "breadth": "Fraction of the element's bases covered by at least one sample read, after "
               "aligning the reads on the tool database (bowtie2). 100 % = the whole element is "
               "there; partial = a truncated cassette, or only a shared promoter or terminator.",
    "identity": "Fraction of the aligned bases identical to the element (1 &minus; mismatches / "
                "aligned length, all reads together). Separates the element itself from a "
                "diverged relative of it.",
    "mean_depth": "Average number of reads covering each base of the element. About the genome "
                  "depth for a single-copy insertion, higher for several copies or a plasmid.",
    "score": "<b>Score = 100 &times; breadth &times; identity.</b> Both must be high: an element "
             "fully covered by a diverged sequence, or a perfect but partial one, is not the "
             "element of the database.",
    "status": "From the score: present (&ge; 99), partial (90 to 99), absent (&lt; 90).",
    "element_id": "Unique identifier of the element in the database: its name, two underscores, "
                  "then the first 6 characters of the SHA-256 fingerprint of its DNA sequence. "
                  "Two elements with the same name but different sequences get different "
                  "suffixes; the same sequence found in several plasmids is kept once. Click the "
                  "? of an element for its details.",
    "marker_status": "Status of the gene measured in the sample: active when the encoded protein "
                      "is intact (or a known functional form), inactive when it is lost "
                      "(deletion, frameshift, stop, or a known loss-of-function allele).",
    "conclusion": "Whether the status measured in the sample (status column) is the one expected "
                  "for the strain called: agrees, differs, or not compared when either status is "
                  "unknown.",
    "decided_by": "<b>The SNP panel first.</b> The strain called is the first of the SNP ranking, "
                  "unless other strains are too close to it: fewer than 200 called sites separate "
                  "them (min_discriminating_sites), or their log-likelihood is within ln 1000 "
                  "(about two net mismatches) of the best. Close relatives such as BY4741 and "
                  "BY4742 differ at a few panel sites or none, but by their mating type and by "
                  "whole-gene deletions (met15D0 against lys2D0).<br><b>Then, among those close "
                  "strains, the mating type and the native markers.</b> Each one is compared with "
                  "the mating type and marker status measured in the sample, from its documented "
                  "genotype or, without one, from the marker alleles of its assembly. A strain "
                  "contradicted at any point is ruled out; among the others, the one agreeing at "
                  "the most points is called, the SNP ranking breaking the ties. If every close "
                  "strain is contradicted, or none can be compared, the first of the SNP ranking "
                  "is kept.",
}

# Element classes, for the per-element popup.
CLASS_HELP = {
    "selection_marker_drug": "drug-resistance gene used to select transformed cells (kanMX: G418, "
                             "natMX: nourseothricin, hphMX: hygromycin).",
    "selection_marker_auxo": "gene restoring growth without a nutrient (HIS3, LEU2, URA3, TRP1, "
                             "...). Yeast has its own copy, so a hit can be the native gene.",
    "reporter": "fluorescent protein or enzyme reporting expression or localisation (GFP, "
                "mCherry, lacZ, ...).",
    "epitope_tag": "short peptide fused to a protein so that antibodies can detect or purify it "
                   "(FLAG, HA, Myc, V5, GST, TAP, ...). 5xFLAG = five FLAG tags in a row.",
    "degron": "sequence making a protein degradable on demand (auxin-inducible degron, ...).",
    "recombination_site": "target of a recombinase (loxP for Cre, FRT for Flp), used to remove "
                          "a marker.",
    "promoter": "promoter of a cassette; it may match the yeast's own copy of that promoter.",
    "terminator": "terminator of a cassette; it may match the yeast's own copy.",
    "other": "not classified.",
}


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
.muted { color: var(--fg-muted); font-size: 0.85rem; }
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
.strain-head {
  display: block; width: 100%; text-align: left; cursor: pointer;
  background: transparent; border: 0; color: inherit; font: inherit;
  padding: 10px 12px;
}
.strain-head[aria-disabled="true"] { cursor: default; }
.strain-head:hover { background: var(--bg-alt); }
.strain-head:focus-visible { outline: 3px solid var(--accent); outline-offset: -3px; }
button.help {
  display: inline-flex; align-items: center; justify-content: center;
  width: 1.15em; height: 1.15em; padding: 0; margin-left: 2px; vertical-align: 0.1em;
  border: 1px solid var(--border-strong); border-radius: 50%;
  background: var(--bg-raised); color: var(--accent);
  font: 700 0.72rem/1 system-ui, sans-serif; cursor: pointer;
  text-transform: none; letter-spacing: 0;
}
button.help:hover, button.help:focus-visible { background: var(--accent); color: var(--accent-fg); outline: none; }
.popover {
  position: absolute; z-index: 20; max-width: min(420px, calc(100vw - 16px));
  background: var(--bg-raised); color: var(--fg);
  border: 1px solid var(--border-strong); border-radius: 6px;
  box-shadow: 0 4px 16px rgba(0, 0, 0, .25);
  padding: 10px 30px 10px 12px; font-size: 0.86rem; line-height: 1.45;
  font-weight: 400; text-transform: none; letter-spacing: 0; text-align: left;
}
.popover p { margin: 0 0 6px 0; }
.popover .pop-close {
  position: absolute; top: 4px; right: 6px; border: 0; background: transparent;
  color: var(--fg-muted); font-size: 1.1rem; line-height: 1; cursor: pointer;
}
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
.scrollx { overflow-x: auto; }
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
.decision { margin: 6px 0; padding-left: 22px; }
.decision > li { margin: 3px 0; }
.decision ul { margin: 2px 0; padding-left: 20px; }
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

var OBS_LABEL = {'R': 'REF', 'A': 'ALT', 'H': 'H (REF + ALT)', '?': '? (not covered)'};
var EXP_LABEL = {'R': 'REF', 'A': 'ALT', 'H': 'H (REF + ALT)', '?': '? (unknown)'};

function alleleLabel(v, labels) {
  if (v === null || v === undefined || v === '') { return '-'; }
  var s = String(v);
  return labels[s] || s;
}

// RefSeq accessions of the S288C R64 chromosomes, shown by their usual names.
var CHROM_NAME = {
  NC_001133: 'chrI', NC_001134: 'chrII', NC_001135: 'chrIII', NC_001136: 'chrIV',
  NC_001137: 'chrV', NC_001138: 'chrVI', NC_001139: 'chrVII', NC_001140: 'chrVIII',
  NC_001141: 'chrIX', NC_001142: 'chrX', NC_001143: 'chrXI', NC_001144: 'chrXII',
  NC_001145: 'chrXIII', NC_001146: 'chrXIV', NC_001147: 'chrXV', NC_001148: 'chrXVI',
  NC_001224: 'chrM'
};

// 'ref|NC_001133|' -> {name: 'chrI', accession: 'NC_001133'}
function chromName(c) {
  var s = (c === null || c === undefined) ? '' : String(c);
  var m = /^(?:ref\|)?([A-Z]{2}_\d+)(?:\.\d+)?\|?$/.exec(s);
  if (!m) { return {name: s || '-', accession: ''}; }
  return {name: CHROM_NAME[m[1]] || m[1], accession: m[1]};
}

// '?' popups: one at a time, closed by its x, Escape or a click elsewhere.
var HELP = null;
var POP = null;

function helpText(key) {
  if (HELP === null) { HELP = JSON.parse(document.getElementById('help-data').textContent); }
  return HELP[key];
}

function closeHelp() {
  if (POP) { POP.box.parentNode.removeChild(POP.box); POP = null; }
}

function para(parent, label, value) {
  var p = document.createElement('p');
  if (label) {
    var b = document.createElement('b');
    b.textContent = label + ': ';
    p.appendChild(b);
  }
  p.appendChild(document.createTextNode(value));
  parent.appendChild(p);
}

// Details of one tool element, from the data (text only, never markup).
function fillElement(body, eid) {
  var tools = data().tools || [], t = null, i;
  for (i = 0; i < tools.length; i++) { if (tools[i].element_id === eid) { t = tools[i]; break; } }
  var cut = eid.lastIndexOf('__');
  var name = cut > 0 ? eid.slice(0, cut) : eid, hash = cut > 0 ? eid.slice(cut + 2) : '';
  para(body, eid, '"' + name + '" is the name of the element in the annotation of the plasmid it ' +
       'comes from' + (hash ? '; "' + hash + '" is the start of the SHA-256 fingerprint of its DNA ' +
       'sequence, which tells apart elements of the same name with different sequences.' : '.'));
  if (!t) { return; }
  var classes = helpText('__classes__') || {};
  var cls = t.element_class || 'other';
  para(body, 'Class', cls + (classes[cls] ? ': ' + classes[cls] : ''));
  para(body, 'Sequence', fmtNum(t.length_bp, 0) + ' bp' +
       (t.feature_type ? ', annotated as ' + t.feature_type : ''));
  if (t.source_plasmid) {
    var also = (t.seen_in || []).filter(function (p) { return p !== t.source_plasmid; });
    para(body, 'Taken from plasmid', t.source_plasmid +
         (also.length ? '; same sequence also in ' + also.join(', ') : ''));
  }
  if (t.origin === 'native') {
    para(body, 'Origin', 'native: this sequence is also in the yeast genome, so it can be ' +
         'found in any strain, with or without a cassette.');
  } else if (t.origin === 'foreign') {
    para(body, 'Origin', 'foreign: absent from the yeast genome, so a hit means a construct.');
  }
}

function openHelp(btn) {
  var again = POP && POP.btn === btn;
  closeHelp();
  if (again) { return; }
  var box = document.createElement('div');
  box.className = 'popover';
  box.setAttribute('role', 'dialog');
  var x = document.createElement('button');
  x.type = 'button';
  x.className = 'pop-close';
  x.setAttribute('aria-label', 'close');
  x.textContent = '×';
  x.addEventListener('click', function (e) { e.stopPropagation(); closeHelp(); btn.focus(); });
  box.appendChild(x);
  var body = document.createElement('div');
  var key = btn.getAttribute('data-help');
  if (key === 'element') {
    fillElement(body, btn.getAttribute('data-eid') || '');
  } else {
    body.innerHTML = helpText(key) || 'No explanation available.';  // authored text, not data
  }
  box.appendChild(body);
  document.body.appendChild(box);
  var r = btn.getBoundingClientRect();
  var width = document.documentElement.clientWidth;
  var left = Math.max(8, Math.min(r.left, width - box.offsetWidth - 8));
  box.style.left = (window.pageXOffset + left) + 'px';
  box.style.top = (window.pageYOffset + r.bottom + 6) + 'px';
  POP = {box: box, btn: btn};
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
    var chrom = chromName(s.chrom);
    var ctd = txt('td', chrom.name);
    if (chrom.accession) { ctd.title = 'RefSeq ' + chrom.accession; }
    tr.appendChild(ctd);
    tr.appendChild(txt('td', s.pos, 'n'));
    tr.appendChild(txt('td', s.ref));
    tr.appendChild(txt('td', s.alt));
    tr.appendChild(txt('td', alleleLabel(s.expected, EXP_LABEL)));
    tr.appendChild(txt('td', alleleLabel(s.observed, OBS_LABEL)));
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
  var heads = document.querySelectorAll('.strain-head[data-si]');
  for (var i = 0; i < heads.length; i++) {
    (function (btn) {
      var si = Number(btn.getAttribute('data-si'));
      if (btn.getAttribute('data-sites') === '0') { return; }
      btn.addEventListener('click', function () { toggleStrain(si); });
      btn.addEventListener('keydown', function (e) {
        if (e.target === btn && (e.key === 'Enter' || e.key === ' ')) {
          e.preventDefault();
          toggleStrain(si);
        }
      });
    })(heads[i]);
  }
  var helps = document.querySelectorAll('button.help');
  for (var h = 0; h < helps.length; h++) {
    (function (b) {
      // Not passed on to the strain row, which would open or close.
      b.addEventListener('click', function (e) { e.stopPropagation(); openHelp(b); });
    })(helps[h]);
  }
  document.addEventListener('click', function (e) {
    if (POP && !POP.box.contains(e.target)) { closeHelp(); }
  });
  document.addEventListener('keydown', function (e) {
    if (e.key === 'Escape' && POP) { var b = POP.btn; closeHelp(); b.focus(); }
  });
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

def fmt_date(value):
    """ISO timestamp as 'YYYY-MM-DD HH:MM'; anything else as given."""
    try:
        return datetime.datetime.fromisoformat(str(value)).strftime("%Y-%m-%d %H:%M")
    except (TypeError, ValueError, AttributeError):
        return value


# RefSeq accessions of the S288C R64 chromosomes, by their usual names (the
# same table as CHROM_NAME in the JS).
CHROM_NAME = dict(("NC_0011%02d" % (33 + i), "chr" + roman) for i, roman in enumerate(
    "I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI".split()))
CHROM_NAME["NC_001224"] = "chrM"


def chrom_names(text):
    """'ref|NC_001137|:115949' -> 'chrV:115949' in free text."""
    return re.sub(r"(?:ref\|)?\b([A-Z]{2}_\d+)(?:\.\d+)?\|?",
                  lambda m: CHROM_NAME.get(m.group(1), m.group(0)), str(text))


def help_btn(key, label):
    """Small '?' button opening the explanation HELP[key] in a popup."""
    return ('<button type="button" class="help" data-help="{k}" '
            'aria-label="Explain {l}" title="Explain {l}">?</button>').format(
        k=esc(key), l=esc(label))


def render_header(d, title):
    panel = d.get("panel") or {}
    rows = [
        ("Organism", d.get("organism"), "organism"),
        ("Reference genome", d.get("reference"), "reference"),
        ("Date", fmt_date(d.get("generated")), "date"),
        ("SNP sites in the panel", num(panel.get("n_sites")), "panel_sites"),
        ("Strains in the panel", num(panel.get("n_strains")), "panel_strains"),
    ]
    items = "\n".join(
        "      <div><dt>{} {}</dt><dd>{}</dd></div>".format(
            esc(k), help_btn(h, k.lower()), esc(v) if v not in (None, "") else "n/a")
        for k, v, h in rows
    )
    return """<header class="page">
  <div class="wrap">
    <h1>{title}</h1>
    <div>Submitted sample <span class="sample">{sample}</span></div>
    <dl class="meta">
{items}
    </dl>
  </div>
</header>""".format(title=esc(title), sample=esc(d.get("sample") or "unknown"), items=items)


VERDICT = {
    "agree": '<span class="status green"><span class="g">&#10003;</span>agree</span>',
    "disagree": '<span class="status red"><span class="g">&#10007;</span>disagree</span>',
}


def verdict(agreement):
    return VERDICT.get(agreement, '<span class="muted">not compared</span>')


OUTCOME = {
    "called": '<span class="status green"><span class="g">&#10003;</span>called</span>',
    "ruled out": '<span class="status red"><span class="g">&#10007;</span>ruled out</span>',
    "not excluded": '<span class="muted">not excluded</span>',
}


def render_decided_by(d):
    """How the strain was called: the SNP panel alone, or the SNP panel then
    the mating type and the native markers among the strains it does not
    separate, each listed with what agrees and what contradicts it."""
    res = d.get("resolution") or {}
    rows = res.get("candidates") or []
    if "method" not in res:   # result file older than the two-step call
        return ""
    if not rows:
        text = "SNP panel alone: no other strain close to {}".format(esc(d.get("call")))
        return '<li><strong>Decided by:</strong> {} {}</li>'.format(
            text, help_btn("decided_by", "how the strain is decided"))
    if res.get("method") == "snp+markers":
        text = ("SNP panel, which does not separate {n} strains, then the mating type and the "
                "native markers").format(n=len(rows))
    else:
        text = ("SNP panel, which does not separate {n} strains; the mating type and the native "
                "markers do not settle them, the first of the SNP ranking is kept").format(
            n=len(rows))
    items = []
    for row in rows:
        # A strain ruled out is shown with what rules it out only.
        points = []
        if row.get("disagree"):
            points.append("contradicted by " + esc("; ".join(row["disagree"])))
        elif row.get("agree"):
            points.append("agrees on " + esc(", ".join(row["agree"])))
        if not points:
            points.append("nothing to compare (no documented genotype)")
        snp = ("first of the SNP ranking" if not row.get("rank") or row["rank"] == 1 else
               "SNP rank {}, {} separating site(s)".format(row["rank"],
                                                          row.get("separating_sites")))
        items.append('    <li><span class="sid">{s}</span> {o} <span class="muted">({snp})'
                     '</span>: {p}</li>'.format(s=esc(row.get("strain")),
                                                o=OUTCOME.get(row.get("outcome"), ""),
                                                snp=esc(snp), p="; ".join(points)))
    return ('<li><strong>Decided by:</strong> {t} {h}<ul>\n{i}\n  </ul></li>').format(
        t=text, h=help_btn("decided_by", "how the strain is decided"), i="\n".join(items))


def render_decision(d):
    """The strain called, how it was decided (SNP panel, then mating type and
    markers among the strains the panel does not separate), then whether the
    mating type, each native marker and each foreign tool element agree with
    it. The evidence is in the sections below."""
    if d.get("unassigned") or not d.get("call"):
        return ""
    strain = esc(d.get("call"))
    # Sex: the mating type in yeast. Set against the strain called only when
    # its name carries one (BY4741 is MATa, W303 is either), else just shown.
    check = d.get("mating_check") or {}
    mating = ""
    measured = check.get("measured") or (d.get("mating_type") or {}).get("call")
    if measured:
        documented = check.get("agreement") in ("agree", "disagree")
        mating = '<li><strong>Sex:</strong> mating type {m}{v}</li>'.format(
            m=esc(measured), v=" " + verdict(check["agreement"]) if documented else "")
    genes = (d.get("markers_expected") or {}).get("genes") or []
    markers = ""
    if genes:
        markers = '<li><strong>Markers:</strong><ul>\n{}\n  </ul></li>'.format("\n".join(
            '    <li>{g} {o} {v}</li>'.format(
                g=esc(g.get("gene")), o=esc(g.get("observed") or ""),
                v=verdict("agree" if g.get("agree") is True else
                          "disagree" if g.get("agree") is False else None))
            for g in genes))
    tool_check = d.get("tool_check")
    if tool_check is None:   # result file older than tool_check
        elements = [{"name": t.get("name"), "status": COLOUR_LABEL[colour_of(t)],
                     "agreement": None}
                    for t in d.get("tools") or []
                    if t.get("origin") != "native" and (t.get("scorable") or "yes") == "yes"
                    and colour_of(t) == "green"]
        documented = False
    else:
        elements = tool_check.get("elements") or []
        documented = tool_check.get("documented")
    if elements:
        tools = '<li><strong>Tools:</strong><ul>\n{}\n  </ul></li>'.format("\n".join(
            '    <li>{n} {s} {v}</li>'.format(n=esc(e.get("name")), s=esc(e.get("status")),
                                          v=verdict(e.get("agreement")))
            for e in elements))
    else:
        tools = '<li><strong>Tools:</strong> no foreign element found {}</li>'.format(
            verdict("agree" if documented else None))
    return """<section id="decision">
  <h2>Decision</h2>
  <ul class="decision">
  <li><strong>Strain:</strong> <span class="sid">{strain}</span></li>
  {decided}
  {mating}
  {markers}
  {tools}
  </ul>
</section>""".format(strain=strain, decided=render_decided_by(d), mating=mating,
                    markers=markers, tools=tools)


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
    items = "\n".join("    <li>{}</li>".format(esc(chrom_names(w))) for w in warns)
    return """<div class="warnings" role="note">
  <h2>Warnings ({n})</h2>
  <ul>
{items}
  </ul>
</div>""".format(n=len(warns), items=items)


# Strains listed in the ranking, the first ones only.
N_STRAINS_SHOWN = 10


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
    for si, s in ordered[:N_STRAINS_SHOWN]:
        sites = s.get("sites") or []
        n_sites = len(sites)
        prob = clamp01(s.get("probability"))
        top = " top" if s.get("rank") == 1 else ""
        disabled = ' aria-disabled="true"' if n_sites == 0 else ""
        chev = "▶" if n_sites else "—"

        # A div acting as the toggle, not a <button>: it holds the '?' help
        # buttons, and a button cannot contain another one.
        head = """  <h3 style="margin:0">
    <div class="strain-head" id="head-{si}" role="button" tabindex="0" data-si="{si}" data-sites="{n}"
         aria-expanded="false" aria-controls="detail-{si}"{dis}>
      <span class="strain-row">
        <span class="chev" id="chev-{si}" aria-hidden="true">{chev}</span>
        <span class="rank">rank {rank}</span>
        <span class="sid">{sid}</span>
        <span class="prob">concordance {conc} {h_conc}</span>
        <span class="prob">PP {p} {h_prob}</span>
        <span class="bar" role="img" aria-label="posterior probability {p}"><i style="width:{w:.2f}%"></i></span>
      </span>
      <span class="counts">
        <span>matching <b>{m}</b> {h_m}</span>
        <span>mismatching <b>{mm}</b> {h_mm}</span>
        <span>uncovered <b>{un}</b> {h_un}</span>
        <span>heterozygous <b>{het}</b> {h_het}</span>
        <span>called <b>{cal}</b> {h_cal}</span>
        <span>log-likelihood <b>{ll}</b> {h_ll}</span>
        <span>{drill}</span>
      </span>
    </div>
  </h3>""".format(
            h_prob=help_btn("probability", "the probability"),
            h_conc=help_btn("concordance", "the concordance"),
            h_m=help_btn("matching", "matching"), h_mm=help_btn("mismatching", "mismatching"),
            h_un=help_btn("uncovered", "uncovered"), h_het=help_btn("heterozygous", "heterozygous"),
            h_cal=help_btn("called", "called"), h_ll=help_btn("log_likelihood", "the log-likelihood"),
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
          <th>site_id</th><th>chrom {h_chrom}</th><th class="n">pos</th><th>ref</th><th>alt</th>
          <th>expected {h_exp}</th><th>observed {h_obs}</th><th class="n">depth {h_depth}</th>
          <th class="n">alt_fraction {h_af}</th><th class="n">weight {h_w}</th><th>state</th>
        </tr></thead>
        <tbody id="sbody-{si}"></tbody>
      </table>
    </div>
  </div>""".format(
                si=si, sid=esc(s.get("strain_id") or "?"),
                m=num(s.get("sites_matching")), mm=num(s.get("sites_mismatching")),
                un=num(s.get("sites_uncovered")), het=num(s.get("sites_het")),
                cal=num(s.get("sites_called")),
                h_chrom=help_btn("chrom", "CHROM"), h_exp=help_btn("expected", "expected"),
                h_obs=help_btn("observed", "observed"), h_depth=help_btn("depth", "depth"),
                h_af=help_btn("alt_fraction", "alt_fraction"), h_w=help_btn("weight", "weight"),
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
  <p class="note">{shown}Click a strain row to expand its per-site detail table.</p>
{blocks}
</section>""".format(blocks="\n".join(blocks),
                     shown=("The {} first of the {} strains scored are shown. ".format(
                         N_STRAINS_SHOWN, len(strains)) if len(strains) > N_STRAINS_SHOWN else ""))


def render_sex(d):
    """Sex of the sample, measured on its reads. Yeast: the mating type, set
    beside the one documented for the strain called when there is one. Other
    organisms: nothing measured yet."""
    mating = d.get("mating_type")
    if not mating:
        return """<section id="sex">
  <h2>Sex</h2>
  <p class="empty">Not determined for this organism.</p>
</section>"""
    # (label, value), both as markup
    check = d.get("mating_check") or {}
    rows = [("Mating type of the sample " + help_btn("mating", "the mating type"),
             esc(mating.get("call") or "unknown"))]
    if mating.get("a_fraction") is not None:
        rows.append(("a fraction a / (a + &alpha;)", esc(num(mating["a_fraction"]))))
    if mating.get("reason"):
        rows.append(("Reason", esc(mating["reason"])))
    if check.get("agreement") in ("agree", "disagree"):
        rows.append(("Documented for " + esc(check.get("strain")),
                     esc(check.get("documented")) + " " + verdict(check["agreement"])))
    items = "\n".join("    <div><dt>{}</dt><dd>{}</dd></div>".format(k, v) for k, v in rows)
    return """<section id="sex">
  <h2>Sex</h2>
  <dl class="meta">
{items}
  </dl>
</section>""".format(items=items)


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
        <td style="white-space:nowrap"><code>{eid}</code> <button type="button" class="help" data-help="element" data-eid="{eid}"
            aria-label="Explain {eid}" title="Explain {eid}">?</button></td>
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


TOOL_HEAD = """      <tr><th>name</th><th>class {}</th><th class="n">length (bp) {}</th>
        <th class="n">breadth {}</th><th class="n">identity {}</th><th class="n">mean depth {}</th>
        <th class="n">score {}</th><th>status {}</th><th>element_id {}</th></tr>""".format(
    help_btn("class", "class"), help_btn("length", "length"), help_btn("breadth", "breadth"),
    help_btn("identity", "identity"), help_btn("mean_depth", "mean depth"),
    help_btn("score", "the score"), help_btn("status", "status"),
    help_btn("element_id", "element_id"))


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

    main_tbl = """  <div class="scrollx"><table>
    <thead>
{head}
    </thead>
    <tbody>
{rows}
    </tbody>
  </table></div>""".format(head=TOOL_HEAD,
                     rows="\n".join(tool_row(t) for t in sorted(
                         scorable,
                         key=lambda t: (str(t.get("element_class") or "~"),
                                        -(t.get("score") if isinstance(t.get("score"), (int, float)) else -1)))))

    if short:
        short_tbl = """  <h3>Short elements &mdash; shown for information only ({n})</h3>
  <p class="short-note">Elements shorter than 100 bp (<code>scorable = "short"</code>) cannot be scored
  reliably by breadth of coverage; they are listed separately and are <strong>not</strong> included in the
  green / orange / red counts above.</p>
  <div class="scrollx"><table>
    <thead>
{head}
    </thead>
    <tbody>
{rows}
    </tbody>
  </table></div>""".format(n=len(short), head=TOOL_HEAD,
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


MARKER_STATUS = {
    # status -> (colour class, glyph, label); never hue alone
    "active": ("green", "&#9679;", "active"),
    "inactive": ("red", "&#9632;", "inactive"),
    "uncertain": ("orange", "&#9650;", "uncertain"),
    "mixed": ("orange", "&#9650;", "mixed"),
    "no_data": ("orange", "&#9650;", "no data"),
}


def expected_cell(e):
    """Expected status of one gene in the first strain of the SNP ranking."""
    if not e or e.get("expected") == "unknown":
        return '<td class="muted">unknown</td>'
    colour, glyph, label = MARKER_STATUS.get(e["expected"], ("orange", "&#9650;", e["expected"]))
    return ('<td><span class="status {c}"><span class="g">{g}</span>{label}</span>'
            '<br><span class="muted">{allele}</span></td>').format(
        c=colour, g=glyph, label=esc(label), allele=esc(e.get("expected_allele")))


def conclusion_cell(e):
    """Whether the sample agrees with the expected status of one gene."""
    return {True: '<td><span class="status green"><span class="g">&#10003;</span>agrees</span></td>',
            False: '<td><span class="status red"><span class="g">&#10007;</span>differs</span></td>',
            None: '<td class="muted">not compared</td>'}[(e or {}).get("agree")]


def marker_row(m, expected=None):
    colour, glyph, label = MARKER_STATUS.get(m.get("status"), ("orange", "&#9650;", "?"))
    alleles = m.get("alleles") or []
    first = alleles[0] if alleles else {}
    changes = first.get("coding_events") or []
    change_text = ", ".join(changes[:6]) + (" (+%d more)" % (len(changes) - 6)
                                            if len(changes) > 6 else "")
    detail = [first.get("reason") or m.get("explanation") or ""]
    detail += m.get("structural") or []
    detail += m.get("notes") or []
    if m.get("heterozygous_calls"):
        detail.append("heterozygous: " + ", ".join(m["heterozygous_calls"][:6]))
    if len(alleles) > 1:
        detail.append("second allele: %s, %s %s" % (alleles[1].get("activity"),
                                                   alleles[1].get("effect"),
                                                   alleles[1].get("protein_change")))
    return """      <tr><td><strong>{gene}</strong><br><span class="muted">{sid}</span></td>
        {expected}
        <td>{effect}<br><span class="muted">{pchange}</span></td>
        <td><code>{changes}</code></td>
        <td>{detail}</td>
        <td><span class="status {c}"><span class="g">{g}</span>{label}</span></td>
        {conclusion}</tr>""".format(
        gene=esc(m.get("gene")), sid=esc(m.get("systematic_id")), c=colour, g=glyph,
        expected=expected_cell(expected) if expected is not None else "",
        conclusion=conclusion_cell(expected) if expected is not None else "",
        label=esc(label), effect=esc((first.get("effect") or "").replace("_", " ")),
        pchange=esc(first.get("protein_change") or ""), changes=esc(change_text or "none"),
        detail="<br>".join(esc(chrom_names(x)) for x in detail if x))


def render_markers(d):
    mk = d.get("markers") or {}
    markers = mk.get("markers") or []
    if not markers:
        return ""
    matches = [m for m in (mk.get("documented_genotype_matches") or [])][:4]
    match_items = "\n".join(
        "    <li><strong>{s}</strong>: {a} of {c} markers agree{d}</li>".format(
            s=esc(m["strain"]), a=m["agree"], c=m["compared"],
            d=esc(" (" + "; ".join(m["differences"]) + ")") if m["differences"] else "")
        for m in matches)
    limits = "\n".join("    <li>{}</li>".format(esc(x)) for x in mk.get("limits") or [])
    # Expected status per gene in the first strain of the SNP ranking. The two
    # readings are kept side by side, neither one corrects the other.
    ex = d.get("markers_expected") or {}
    expected_by_gene = {g["gene"]: g for g in ex.get("genes") or []}
    if ex:
        expected_head = "<th>expected in {}</th>".format(esc(ex.get("strain")))
        conclusion_head = "<th>conclusion {}</th>".format(help_btn("conclusion", "the conclusion"))
        expected_note = """  <p><strong>Expected from the SNP ranking:</strong> the status each marker should have in
  <strong>{s}</strong>, first strain of the SNP ranking ({src}). {summary}</p>""".format(
            s=esc(ex.get("strain")), src=esc(ex.get("source")),
            summary=("The sample agrees at <strong>{a} of {c}</strong> markers compared.".format(
                a=ex.get("agree"), c=ex.get("compared")) if ex.get("compared")
                     else "No marker could be compared."))
    else:
        expected_head, expected_note, conclusion_head = "", "", ""
    return """<section id="markers">
  <h2>Native marker genes</h2>
  <p class="note">The strain's own selection markers, judged by the protein each one encodes. A
  gene can be fully covered and still inactive (a point mutation such as <code>ura3-1</code>), or
  carry many SNPs and still be active (a natural form). Loss of function is called from the
  sequence; a missense form is called only when the allele catalogue of {n} assemblies knows it.</p>
  <p><strong>Genotype:</strong> <code>{genotype}</code></p>
{expected_note}
  <table>
    <thead><tr><th>gene</th>{expected_head}<th>effect</th><th>coding changes</th>
      <th>evidence</th><th>status {status_help}</th>{conclusion_head}</tr></thead>
    <tbody>
{rows}
    </tbody>
  </table>
  <h3>Documented genotypes</h3>
  <ul>
{matches}
  </ul>
  <h3>Limits</h3>
  <ul>
{limits}
  </ul>
</section>""".format(n=esc(mk.get("catalogue_assemblies")), genotype=esc(mk.get("genotype")),
                     expected_note=expected_note, expected_head=expected_head,
                     conclusion_head=conclusion_head,
                     status_help=help_btn("marker_status", "the marker status"),
                     rows="\n".join(marker_row(m, expected_by_gene.get(m.get("gene"), {})
                                               if ex else None) for m in markers),
                     matches=match_items or "    <li>none documented</li>", limits=limits)


# --------------------------------------------------------------------------- #
# document
# --------------------------------------------------------------------------- #

def render(d, title):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    help_data = dict(HELP, __classes__=CLASS_HELP)
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
        render_decision(d),
        render_warnings(d),
        render_strains(d),
        render_sex(d),
        render_markers(d),
        render_tools(d),
        '<footer class="page">BioSentinelle report &mdash; rendered {} by '
        "<code>render_report.py</code>. Self-contained: no network access required."
        "</footer>".format(esc(now)),
        "</div>",
        '<script type="application/json" id="report-data">{}</script>'.format(embed_json(d)),
        '<script type="application/json" id="help-data">{}</script>'.format(embed_json(help_data)),
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
    ap.add_argument("--title", default="BioSentinelle report", help="report title")
    args = ap.parse_args(argv)

    with open(args.json, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        sys.stderr.write("error: the JSON root must be an object\n")
        return 1

    title = args.title

    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(render(data, title))

    sys.stderr.write("wrote {}\n".format(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
