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
                            [--tree panel_tree.nwk --strains panel_strains.tsv]

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

TREE OF THE PANEL STRAINS
-------------------------
With --tree (and --strains), the "Strains in the panel" entry of the header
gets a "tree" button that opens the tree of every panel strain in a new page,
the strain called highlighted. The page is written from four blocks embedded
in the report (its style, body, script and data), so the report stays one
self-contained file; when the browser blocks the new window, the same page
opens over the report.
"""

import argparse
import csv
import datetime
import html
import json
import math
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
_SVG = ('<svg class="ico" viewBox="0 0 24 24" width="1em" height="1em" fill="none" '
        'stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round" '
        'aria-hidden="true">{}</svg>')
COLOUR_GLYPH = {
    "green": _SVG.format('<circle cx="12" cy="12" r="10"/><path d="M8 12.5l3 3 5-6"/>'),
    "orange": _SVG.format('<path d="M12 3L2 20h20L12 3z"/><path d="M12 10v5M12 17.5v.5"/>'),
    "red": _SVG.format('<circle cx="12" cy="12" r="10"/><path d="M9 9l6 6M15 9l-6 6"/>'),
}
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
    "a_fraction": "<b>Fraction a / (a + &alpha;)</b>: share of the sample's mating-type k-mers "
                  "that come from the a cassette. a = mean count, in the sample reads, of the "
                  "k-mers specific to the MATa cassette; &alpha; = same for the MAT&alpha; "
                  "cassette. Every genome carries one silent a (HMR) and one silent &alpha; (HML) "
                  "copy plus the active MAT copy, hence about 0.67 in a MATa haploid (2 a, "
                  "1 &alpha;), 0.33 in a MAT&alpha; haploid (1 a, 2 &alpha;) and 0.5 in a "
                  "MATa/&alpha; diploid or a mix of both. Calls: MATa &ge; 0.58, MATalpha &le; "
                  "0.42, MATa/alpha in between.",
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
    "mito": "Growth on glycerol needs respiration, hence a mitochondrial genome (mtDNA) that "
            "still encodes the respiratory chain. The reads aligned to the mtDNA of the reference "
            "are counted over the exons of the genes needed for respiration (COX1, COX2, COX3, "
            "COB, ATP6, ATP8, OLI1, VAR1, the two rRNAs, RPM1 and the tRNAs); the optional introns "
            "are left out, since strains differ in which they carry. <b>rho+</b>: every gene "
            "covered, growth on glycerol expected. <b>rho-</b>: whole genes missing, the mtDNA is "
            "reduced to a fragment, no growth. <b>rho0</b>: no mtDNA, no growth. "
            "<b>uncertain</b>: every gene there but some only in part. A point mutation in a "
            "mitochondrial gene (mit-) or a nuclear respiratory mutation (pet) is not seen.",
    "marker_status": "Status of the gene measured in the sample: active when the encoded protein "
                      "is intact (or a known functional form), inactive when it is lost "
                      "(deletion, frameshift, stop, or a known loss-of-function allele).",
    "conclusion": "Whether the status measured in the sample (status column) is the one expected "
                  "for the strain called: agrees, differs, or not compared when either status is "
                  "unknown.",
    "decided_by": "<b>The SNP panel first.</b> The strain called is the first of the SNP ranking, "
                  "unless other strains are too close to it for the SNPs alone: their "
                  "log-likelihood is within ln 1000 (about two net mismatches) of the best, or "
                  "fewer than 200 called sites separate them (min_discriminating_sites). A "
                  "separating site is a panel site called in the sample where the genotypes of "
                  "the two strains are known and differ. The second rule does not look at the "
                  "allele the sample carries there: it may be the allele of the best strain at "
                  "every separating site, as shown for each strain, but a lead built on so few "
                  "sites can come from panel errors, a misassembled region or the drift of a "
                  "laboratory stock. Such a strain can therefore have a PP of 0 and rank below "
                  "the 10 strains shown in Strain ranking: the PP treats the sites as "
                  "independent, so that a few of them already bring it to 0, and the "
                  "log-likelihood also pays at the sites where the genotype of the strain is "
                  "unknown. Close relatives such as BY4741 and "
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
.ico { vertical-align: -0.14em; }
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
.documented { list-style: none; padding-left: 0; }
.documented > li { margin: 0 0 10px 0; }
.marker-grid {
  display: grid; grid-template-columns: max-content max-content max-content;
  column-gap: 14px; row-gap: 2px; margin: 2px 0;
}
.legend { border: 1px solid var(--border); border-radius: 6px; padding: 8px 14px; background: var(--bg-alt); font-size: 0.88rem; }
.legend ul { margin: 4px 0 0 0; padding-left: 20px; }
.short-note { font-size: 0.87rem; color: var(--fg-muted); margin: 4px 0 8px 0; }
.empty { color: var(--fg-muted); font-style: italic; padding: 6px 0; }
code { font-family: var(--mono); }
footer.page { margin-top: 36px; border-top: 1px solid var(--border); padding-top: 10px; color: var(--fg-muted); font-size: 0.82rem; }
button.tree-btn {
  display: inline-flex; align-items: center; gap: 4px; margin-left: 8px; vertical-align: 0.08em;
  padding: 1px 8px; border: 1px solid var(--border-strong); border-radius: 10px;
  background: var(--bg-raised); color: var(--accent);
  font: 600 0.76rem/1.35 system-ui, -apple-system, "Segoe UI", sans-serif; cursor: pointer;
}
button.tree-btn:hover, button.tree-btn:focus-visible { background: var(--accent); color: var(--accent-fg); outline: none; }
.tree-overlay {
  position: fixed; top: 0; right: 0; bottom: 0; left: 0; z-index: 50;
  display: flex; flex-direction: column; background: var(--bg);
}
.tree-overlay[hidden] { display: none; }
.tree-overlay-bar {
  display: flex; flex-wrap: wrap; gap: 10px; align-items: center;
  padding: 6px 16px; border-bottom: 1px solid var(--border); font-size: 0.85rem; background: var(--bg-alt);
}
.tree-overlay-bar button {
  font: inherit; padding: 2px 10px; border: 1px solid var(--border-strong); border-radius: 4px;
  background: var(--bg-raised); color: var(--fg); cursor: pointer;
}
.tree-overlay iframe { flex: 1 1 auto; width: 100%; border: 0; }
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
  var mk = (data().markers || {}).markers || [];
  for (i = 0; i < mk.length; i++) {
    if (String(t.name || name).toUpperCase() === String(mk[i].gene).toUpperCase()) {
      para(body, 'Versus the marker gene', 'this element carries the name of the native marker ' +
           mk[i].gene + ', which is judged in the Native marker genes section from the genomic ' +
           'copy (reads on the chromosome locus). The element is the copy carried by the ' +
           'plasmid (' + (t.source_plasmid || 'unknown') + '), told apart by its own sequence ' +
           '(the fingerprint after the two underscores) and by its origin: ' +
           (t.origin === 'foreign' ? 'foreign here, its sequence differs from the genomic gene, ' +
            'so a hit means a construct, not the strain\'s own marker.' :
            'native, its sequence is the genomic one, so a hit alone does not prove a construct.'));
      break;
    }
  }
}

// Characteristics of one panel strain, with links to look it up (text and links only).
function fillStrain(body, si) {
  var s = (data().strains || [])[si] || {}, id = s.strain_id || '', i, a;
  para(body, 'Strain', (s.display_name && s.display_name !== id) ? s.display_name + ' (' + id + ')' : id);
  var rows = [['Clade', s.clade], ['Origin', s.origin], ['Ploidy', s.ploidy],
              ['Zygosity', s.zygosity], ['Sequence source', s.source_type]];
  for (i = 0; i < rows.length; i++) { if (rows[i][1]) { para(body, rows[i][0], rows[i][1]); } }
  var p = document.createElement('p'), links = [['SGD', 'https://www.yeastgenome.org/search?q=' +
        encodeURIComponent(id)], ['PubMed', 'https://pubmed.ncbi.nlm.nih.gov/?term=' +
        encodeURIComponent('"' + (s.display_name || id) + '" yeast')]];
  var runs = String(s.source || '').split(',').filter(function (r) { return /^[SED]RR\d+$/.test(r); });
  for (i = 0; i < runs.length; i++) {
    links.push([runs[i] + ' (ENA)', 'https://www.ebi.ac.uk/ena/browser/view/' + runs[i]]);
  }
  p.appendChild(document.createTextNode('Links: '));
  for (i = 0; i < links.length; i++) {
    a = document.createElement('a');
    a.href = links[i][1]; a.target = '_blank'; a.rel = 'noopener noreferrer';
    a.textContent = links[i][0];
    if (i) { p.appendChild(document.createTextNode(' · ')); }
    p.appendChild(a);
  }
  body.appendChild(p);
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
  } else if (key === 'strain') {
    fillStrain(body, +btn.getAttribute('data-si'));
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

// Tree of the panel strains, in a page of its own. The page is put together
// from the blocks embedded at the end of this report (style, body, data,
// script) and written into a new window, so that the report stays one
// self-contained file. If the browser blocks the new window, the same page
// opens over the report, in a frame.
function treeDocument() {
  function block(id) { var e = document.getElementById(id); return e ? e.textContent : ''; }
  return '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">' +
         '<meta name="viewport" content="width=device-width, initial-scale=1">' +
         '<title>Tree of the panel strains<\/title><style>' + block('tree-css') + '<\/style><\/head>' +
         '<body>' + block('tree-body') +
         '<script type="application/json" id="tree-data">' + block('tree-data') + '<\/script>' +
         '<script>' + block('tree-js') + '<\/script><\/body><\/html>';
}

function treeOverlay(html) {
  var box = document.getElementById('tree-overlay');
  if (!box) {
    box = document.createElement('div');
    box.id = 'tree-overlay';
    box.className = 'tree-overlay';
    box.setAttribute('role', 'dialog');
    box.setAttribute('aria-label', 'Tree of the panel strains');
    var bar = document.createElement('div');
    bar.className = 'tree-overlay-bar';
    var back = document.createElement('button');
    back.type = 'button';
    back.textContent = '← Back to the report';
    back.addEventListener('click', closeTreeOverlay);
    bar.appendChild(back);
    bar.appendChild(document.createTextNode(
      'The browser did not let the tree open in a new page, so it is shown here, over the report.'));
    box.appendChild(bar);
    var frame = document.createElement('iframe');
    frame.title = 'Tree of the panel strains';
    frame.srcdoc = html;
    box.appendChild(frame);
    document.body.appendChild(box);
  }
  box.hidden = false;
  box.querySelector('button').focus();
}

function closeTreeOverlay() {
  var box = document.getElementById('tree-overlay');
  if (box && !box.hidden) {
    box.hidden = true;
    var b = document.getElementById('open-tree');
    if (b) { b.focus(); }
  }
}

function openTree() {
  var html = treeDocument(), w = null;
  try { w = window.open('', '_blank'); } catch (e) { w = null; }
  if (w) {
    try {
      w.document.open();
      w.document.write(html);
      w.document.close();
      return;
    } catch (e) {
      try { w.close(); } catch (e2) { /* nothing to close */ }
    }
  }
  treeOverlay(html);
}

document.addEventListener('DOMContentLoaded', function () {
  var treeBtn = document.getElementById('open-tree');
  if (treeBtn) { treeBtn.addEventListener('click', openTree); }
  document.addEventListener('keydown', function (e) {
    if (e.key === 'Escape') { closeTreeOverlay(); }
  });
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
# tree page
# --------------------------------------------------------------------------- #

# The tree of the panel strains opens in a page of its own, written by the
# report's openTree() from four blocks embedded in the report: this style,
# this body, the tree data (tree_payload) and this script. None of them may
# hold '</script', which would end the block that carries it (checked in
# render()). Colour carries two marks only, the selected strain (blue) and the
# best ranked strains (orange), each with its own legend entry and label;
# clades are neutral bands named in text, never hues.

TREE_CSS = """
:root {
  color-scheme: light dark;
  --bg: #ffffff; --bg-alt: #f4f6f8; --fg: #16191d; --fg-muted: #5a6572;
  --border: #d2d8de; --border-strong: #9aa4b0; --accent: #1f5fa9; --accent-fg: #ffffff;
  --branch: #52514e; --leader: #e1e0d9; --node: #898781; --hover: #0b0b0b;
  --sel: #2a78d6; --rank: #eb6834; --band-a: #e1e0d9; --band-b: #c3c2b7;
  --row-sel: rgba(42, 120, 214, 0.16); --row-hov: rgba(11, 11, 11, 0.06);
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #14171b; --bg-alt: #1c2126; --fg: #e7ecf1; --fg-muted: #a3aeba;
    --border: #39424c; --border-strong: #5c6773; --accent: #6aa9ef; --accent-fg: #0d1014;
    --branch: #c3c2b7; --leader: #2c2c2a; --node: #898781; --hover: #f0efec;
    --sel: #3987e5; --rank: #d95926; --band-a: #2c2c2a; --band-b: #4a4a46;
    --row-sel: rgba(57, 135, 229, 0.24); --row-hov: rgba(255, 255, 255, 0.08);
  }
}
* { box-sizing: border-box; }
html, body { height: 100%; }
body {
  margin: 0; display: flex; flex-direction: column; background: var(--bg); color: var(--fg);
  font-family: system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  font-size: 14px; line-height: 1.45;
}
[hidden] { display: none !important; }
header.bar { flex: none; background: var(--bg-alt); border-bottom: 1px solid var(--border); padding: 8px 16px 6px 16px; }
header.bar h1 { font-size: 1.15rem; margin: 0; }
.sub { margin: 0; color: var(--fg-muted); font-size: 0.85rem; }
.tools { display: flex; flex-wrap: wrap; gap: 6px 16px; align-items: center; margin-top: 6px; }
.grp { display: inline-flex; flex-wrap: wrap; align-items: center; gap: 4px; }
.tools button, .tools input, .tools a, .tools label { font: inherit; font-size: 0.84rem; }
.tools button {
  padding: 2px 9px; border: 1px solid var(--border-strong); border-radius: 4px;
  background: var(--bg); color: var(--fg); cursor: pointer;
}
.tools button:hover:not(:disabled) { border-color: var(--accent); }
.tools button:disabled { opacity: .45; cursor: default; }
.tools button[aria-pressed="true"] { background: var(--accent); color: var(--accent-fg); border-color: var(--accent); }
.tools button:focus-visible, .tools input:focus-visible, .tools a:focus-visible, #c-ranked:focus-visible,
.info button:focus-visible, .info a:focus-visible { outline: 3px solid var(--accent); outline-offset: 1px; }
.tools input {
  width: 13em; padding: 2px 6px; border: 1px solid var(--border-strong); border-radius: 4px;
  background: var(--bg); color: var(--fg);
}
.tools a { color: var(--accent); }
.zl, .msg { font-size: 0.8rem; color: var(--fg-muted); }
.legend { display: flex; flex-wrap: wrap; gap: 4px 18px; align-items: center; margin-top: 6px; font-size: 0.82rem; color: var(--fg-muted); }
.key { display: inline-flex; align-items: center; gap: 5px; }
.key svg { flex: none; }
.scale i { display: inline-block; height: 0; border-top: 2px solid var(--branch); vertical-align: middle; }
details.about { margin-top: 4px; font-size: 0.82rem; color: var(--fg-muted); max-width: 75em; }
details.about summary { cursor: pointer; color: var(--accent); width: max-content; }
details.about p { margin: 4px 0; }
#stage { flex: 1 1 auto; min-height: 0; overflow: auto; position: relative; cursor: grab; outline: none; }
#stage.dragging { cursor: grabbing; }
#stage:focus-visible { box-shadow: inset 0 0 0 3px var(--accent); }
#stage .msg { display: block; padding: 16px; }
#tree { display: block; -webkit-user-select: none; user-select: none; }
.tip {
  position: fixed; z-index: 30; pointer-events: none; max-width: 340px;
  background: var(--bg); color: var(--fg); border: 1px solid var(--border-strong); border-radius: 6px;
  box-shadow: 0 4px 14px rgba(0, 0, 0, .22); padding: 7px 10px; font-size: 0.82rem; line-height: 1.4;
}
.tip b { font-size: 0.92rem; }
.tip div, .info p { margin: 0 0 2px 0; }
.muted { color: var(--fg-muted); }
.key-line { display: inline-block; width: 14px; height: 0; border-top: 3px solid; vertical-align: middle; margin-right: 6px; }
.info {
  position: fixed; z-index: 25; right: 16px; bottom: 16px; width: min(380px, calc(100vw - 32px));
  max-height: 60vh; overflow: auto; background: var(--bg); color: var(--fg);
  border: 1px solid var(--border-strong); border-radius: 8px; box-shadow: 0 6px 22px rgba(0, 0, 0, .25);
  padding: 10px 34px 10px 14px; font-size: 0.86rem;
}
.info h2 { font-size: 1rem; margin: 0 0 4px 0; }
.info a { color: var(--accent); }
.info .x {
  position: absolute; top: 4px; right: 8px; border: 0; background: transparent;
  color: var(--fg-muted); font-size: 1.2rem; line-height: 1; cursor: pointer;
}
.info .acts { margin-top: 8px; }
.info .acts button {
  font: inherit; font-size: 0.84rem; padding: 2px 9px; border: 1px solid var(--border-strong);
  border-radius: 4px; background: var(--bg); color: var(--fg); cursor: pointer;
}
"""

TREE_BODY = """
<header class="bar">
  <h1>Tree of the panel strains</h1>
  <p class="sub" id="t-sub"></p>
  <div class="tools" role="toolbar" aria-label="Tree controls">
    <div class="grp" role="group" aria-label="Layout">
      <button type="button" data-layout="rect" aria-pressed="true">Rectangular</button>
      <button type="button" data-layout="circ" aria-pressed="false">Circular</button>
    </div>
    <div class="grp" role="group" aria-label="Part of the tree shown">
      <button type="button" id="b-whole">Whole tree</button>
      <button type="button" id="b-around">Around the selected strain</button>
      <button type="button" id="b-up">Up one level</button>
    </div>
    <div class="grp" role="group" aria-label="Zoom">
      <span class="zl" id="z-label">Height</span>
      <button type="button" id="b-out" aria-label="Zoom out" title="Zoom out (also Ctrl + mouse wheel)">&minus;</button>
      <button type="button" id="b-in" aria-label="Zoom in" title="Zoom in (also Ctrl + mouse wheel)">+</button>
      <span class="zl rect-only">Branch length</span>
      <button type="button" class="rect-only" id="b-narrow" aria-label="Shorter branches" title="Shorter branches">&minus;</button>
      <button type="button" class="rect-only" id="b-wide" aria-label="Longer branches" title="Longer branches">+</button>
      <button type="button" id="b-fit" title="Fit the window">Fit</button>
    </div>
    <div class="grp">
      <label for="find">Find a strain</label>
      <input type="search" id="find" list="strain-names" placeholder="name or 1011 code" autocomplete="off" spellcheck="false">
      <datalist id="strain-names"></datalist>
      <span class="msg" id="find-msg" aria-live="polite"></span>
    </div>
    <div class="grp">
      <span class="zl">Download</span>
      <a id="dl-nwk" href="#" download="panel_tree.nwk">Newick</a>
      <a id="dl-svg" href="#" download="panel_tree.svg">SVG</a>
    </div>
  </div>
  <div class="legend">
    <span class="key"><svg width="14" height="14" aria-hidden="true"><circle cx="7" cy="7" r="5.5" class="mk-sel"/></svg><span id="key-sel"></span></span>
    <label class="key"><input type="checkbox" id="c-ranked" checked><svg width="14" height="14" aria-hidden="true"><circle cx="7" cy="7" r="4" class="mk-rank"/></svg><span id="key-rank"></span></label>
    <span class="key"><svg width="12" height="16" aria-hidden="true"><rect x="0" y="0" width="12" height="8" class="band-a"/><rect x="0" y="8" width="12" height="8" class="band-b"/></svg>clades of Peter et al. 2018, named next to their band</span>
    <span class="key scale"><i id="scale-bar"></i><span id="scale-text"></span></span>
    <span class="key" id="crumb"></span>
  </div>
  <details class="about"><summary>About this tree</summary><div id="about"></div></details>
</header>
<main id="stage" tabindex="0" aria-label="Tree of the panel strains: drag to move, Ctrl + wheel to zoom">
  <svg id="tree" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Tree of the panel strains"></svg>
</main>
<div class="tip" id="tip" role="tooltip" hidden></div>
<aside class="info" id="info" hidden></aside>
"""

TREE_JS = r"""
'use strict';
(function () {
var T = JSON.parse(document.getElementById('tree-data').textContent);
var PARENT = T.parent, LEN = T.length, NAME = T.name, N = PARENT.length;
var SVGNS = 'http://www.w3.org/2000/svg';
var stage = document.getElementById('stage'), svg = document.getElementById('tree');
var tip = document.getElementById('tip'), infoBox = document.getElementById('info');
var i;

function fail(msg) {
  var p = document.createElement('p');
  p.className = 'msg';
  p.textContent = msg;
  stage.textContent = '';
  stage.appendChild(p);
}
if (!N) { fail('The tree is empty.'); return; }
for (i = 1; i < N; i++) {
  if (!(PARENT[i] >= 0 && PARENT[i] < i)) { fail('The tree data are not in preorder: the tree cannot be drawn.'); return; }
}

// Nodes are in preorder: a node comes before its descendants, which follow it
// in one block, and children are in drawing order. Leaves carry strain names.
var KIDS = [];
for (i = 0; i < N; i++) { KIDS.push([]); }
for (i = 1; i < N; i++) { KIDS[PARENT[i]].push(i); }
var NNODES = new Int32Array(N), NLEAVES = new Int32Array(N);
for (i = N - 1; i >= 0; i--) {
  NNODES[i] += 1;
  if (!KIDS[i].length) { NLEAVES[i] += 1; }
  if (i > 0) { NNODES[PARENT[i]] += NNODES[i]; NLEAVES[PARENT[i]] += NLEAVES[i]; }
}
var LEAF = {}, FIELD = {};
for (i = 0; i < N; i++) { if (!KIDS[i].length) { LEAF[NAME[i]] = i; } }
T.fields.forEach(function (f, k) { FIELD[f] = k; });

function meta(v, field) {
  var row = T.strains[NAME[v]];
  return row && FIELD[field] !== undefined ? (row[FIELD[field]] || '') : '';
}
function label(v) { return meta(v, 'display_name') || NAME[v]; }
function cladeOf(v) { return meta(v, 'clade'); }
function shortClade(c) { var m = /^(M?\d+)\./.exec(c); return m ? m[1] : c.split(' ')[0]; }

var RANK = {};
(T.ranked || []).forEach(function (r) {
  var v = LEAF[r[0]];
  if (v !== undefined) { RANK[v] = {rank: r[1], pp: r[2], conc: r[3]}; }
});
var SEL = (T.selected && LEAF[T.selected] !== undefined) ? LEAF[T.selected] : -1;
var BEST = -1;
Object.keys(RANK).forEach(function (k) {
  var v = Number(k);
  if (BEST < 0 || RANK[v].rank < RANK[BEST].rank) { BEST = v; }
});

// drawing state; fit: the drawing follows the size of the window
var S = {layout: 'rect', sub: 0, row: 12, wide: 1, size: 1, ranked: true, fit: false};
var X = new Float64Array(N), Y = new Float64Array(N), DEP = new Float64Array(N);
var ANG = new Float64Array(N), RAD = new Float64Array(N);
var CX = 0, CY = 0;
var G = null, HOVER = null, PIN = -1;
var LEFT = 24, TOP = 16, BOTTOM = 18, GAP = 8, STRIP = 10, MARGIN = 14, ROW_MAX = 26, CALLOUT = 26;

// ---------------------------------------------------------------- helpers
var CTX = document.createElement('canvas').getContext('2d');
var FAMILY = window.getComputedStyle(document.body).fontFamily || 'sans-serif';
var WIDTH = {};
function textW(str, size, weight) {
  var key = (weight || 400) + '|' + str, w = WIDTH[key];
  if (w === undefined) {
    CTX.font = (weight || 400) + ' 12px ' + FAMILY;
    w = WIDTH[key] = CTX.measureText(str).width;
  }
  return w * size / 12;
}
function r1(v) { return Math.round(v * 10) / 10; }
function el(tag, attrs, parent) {
  var e = document.createElementNS(SVGNS, tag);
  for (var k in attrs) {
    if (Object.prototype.hasOwnProperty.call(attrs, k)) { e.setAttribute(k, attrs[k]); }
  }
  if (parent) { parent.appendChild(e); }
  return e;
}
function svgText(x, y, str, cls, size, parent) {
  var e = el('text', {x: r1(x), y: r1(y), dy: '0.35em', 'class': cls, 'font-size': size}, parent);
  e.textContent = str;
  return e;
}
function inSub(v, g) { return v >= g.s && v < g.end; }
function leafText(v) { return label(v) + (S.ranked && RANK[v] ? '  #' + RANK[v].rank : ''); }
function marked(g) {
  // the selected strain and the best ranked ones in the part shown, top first
  var out = [];
  if (SEL >= 0 && inSub(SEL, g)) { out.push(SEL); }
  if (S.ranked) {
    Object.keys(RANK).map(Number).sort(function (a, b) { return RANK[a].rank - RANK[b].rank; })
      .forEach(function (v) { if (v !== SEL && inSub(v, g)) { out.push(v); } });
  }
  return out;
}
function niceLength(max) {
  var p = Math.pow(10, Math.floor(Math.log(max) / Math.LN10)), m = max / p;
  return (m >= 5 ? 5 : m >= 2 ? 2 : 1) * p;
}
function fmtLength(v) {
  return String(+v.toPrecision(2));
}

// ---------------------------------------------------------------- geometry
function subtree(s) {
  var end = s + NNODES[s], leaves = [], inner = [], maxD = 0, v;
  DEP[s] = 0;
  for (v = s; v < end; v++) {
    if (v > s) { DEP[v] = DEP[PARENT[v]] + Math.max(0, LEN[v]); }
    if (KIDS[v].length) { inner.push(v); } else { leaves.push(v); if (DEP[v] > maxD) { maxD = DEP[v]; } }
  }
  return {s: s, end: end, leaves: leaves, inner: inner, maxD: maxD > 0 ? maxD : 1};
}

function cladeRuns(leaves) {
  // consecutive leaves of one clade; up to three leaves without a clade
  // between two leaves of the same clade stay in its band
  var runs = [], cur = null, k, c;
  for (k = 0; k < leaves.length; k++) {
    c = cladeOf(leaves[k]);
    if (!c) { continue; }
    if (cur && cur.clade === c && k - cur.k1 <= 4) { cur.k1 = k; }
    else { cur = {clade: c, k0: k, k1: k}; runs.push(cur); }
  }
  return runs;
}

function spread(items, gap, lo, hi) {
  // 1-D placement of the labels of the marked strains, in their order
  var k, n = items.length;
  for (k = 0; k < n; k++) { items[k].p = Math.max(items[k].p, k ? items[k - 1].p + gap : lo); }
  if (n && items[n - 1].p > hi) {
    items[n - 1].p = hi;
    for (k = n - 2; k >= 0; k--) { items[k].p = Math.min(items[k].p, items[k + 1].p - gap); }
  }
}

function pt(r, a) { return r1(CX + r * Math.sin(a)) + ',' + r1(CY - r * Math.cos(a)); }
function arcTo(r, a0, a1) {
  var sweep = a1 > a0 ? 1 : 0, rr = r1(r), head = 'A' + rr + ',' + rr + ' 0 0,' + sweep + ' ';
  if (Math.abs(a1 - a0) > Math.PI) { return head + pt(r, (a0 + a1) / 2) + head + pt(r, a1); }
  return head + pt(r, a1);
}
function sector(r0, r1_, a0, a1) {
  var m = (a0 + a1) / 2, o0 = r1(r0), o1 = r1(r1_);
  return 'M' + pt(r0, a0) + 'A' + o0 + ',' + o0 + ' 0 0,1 ' + pt(r0, m) + 'A' + o0 + ',' + o0 + ' 0 0,1 ' +
    pt(r0, a1) + 'L' + pt(r1_, a1) + 'A' + o1 + ',' + o1 + ' 0 0,0 ' + pt(r1_, m) + 'A' + o1 + ',' + o1 +
    ' 0 0,0 ' + pt(r1_, a0) + 'Z';
}
function dot(x, y, r) {
  return 'M' + r1(x - r) + ',' + r1(y) + 'a' + r + ',' + r + ' 0 1,0 ' + 2 * r + ',0a' + r + ',' + r +
    ' 0 1,0 ' + (-2 * r) + ',0';
}
// branch from u to its child v: the part drawn for v
function branch(u, v) {
  if (S.layout === 'circ') {
    return 'M' + pt(RAD[u], ANG[v]) + 'L' + pt(RAD[v], ANG[v]);
  }
  return 'M' + r1(X[u]) + ',' + r1(Y[v]) + 'H' + r1(X[v]);
}
// the joint of u: vertical line or arc over its children
function joint(u) {
  var c = KIDS[u];
  if (S.layout === 'circ') {
    return RAD[u] > 0.5 ? 'M' + pt(RAD[u], ANG[c[0]]) + arcTo(RAD[u], ANG[c[0]], ANG[c[c.length - 1]]) : '';
  }
  return 'M' + r1(X[u]) + ',' + r1(Y[c[0]]) + 'V' + r1(Y[c[c.length - 1]]);
}
// path from the top of the part shown down to v
function pathTo(v, g) {
  var d = [], u;
  while (v !== g.s && v > 0) {
    u = PARENT[v];
    d.push(branch(u, v));
    if (S.layout === 'circ') {
      if (RAD[u] > 0.5) { d.push('M' + pt(RAD[u], ANG[u]) + arcTo(RAD[u], ANG[u], ANG[v])); }
    } else {
      d.push('M' + r1(X[u]) + ',' + r1(Y[u]) + 'V' + r1(Y[v]));
    }
    v = u;
  }
  return d.join('');
}
function subtreePath(v) {
  var d = [], end = v + NNODES[v], u;
  for (u = v + 1; u < end; u++) { d.push(branch(PARENT[u], u)); }
  for (u = v; u < end; u++) { if (KIDS[u].length) { d.push(joint(u)); } }
  return d.join('');
}

function layers(parent) {
  var L = {};
  ['bands', 'rows', 'leaders', 'branches', 'hover', 'path', 'dots', 'marks', 'labels', 'clades', 'top']
    .forEach(function (name) { L[name] = el('g', {'class': 'layer-' + name}, parent); });
  return L;
}

function leafLabel(v, x, y, size, parent, a) {
  var e = svgText(x, y, label(v), v === SEL ? 'lb lb-sel' : 'lb', size, parent);
  if (S.ranked && RANK[v]) {
    var t = document.createElementNS(SVGNS, 'tspan');
    t.setAttribute('class', 'rk');
    t.textContent = '  #' + RANK[v].rank;
    e.appendChild(t);
  }
  e.setAttribute('data-v', v);
  if (a !== undefined) { rotate(e, x, y, a); }
  return e;
}
function rotate(e, x, y, a) {
  // radial text at angle a (from the top, clockwise), upright on the left half
  var deg = a * 180 / Math.PI - 90, flip = a > Math.PI;
  e.setAttribute('transform', 'rotate(' + r1(flip ? deg + 180 : deg) + ',' + r1(x) + ',' + r1(y) + ')');
  if (flip) { e.setAttribute('text-anchor', 'end'); }
}

// ---------------------------------------------------------------- rectangular
function drawRect(g, root) {
  var n = g.leaves.length, w0 = Math.max(stage.clientWidth, 320), h0 = Math.max(stage.clientHeight, 200);
  var rowMin = Math.min(Math.max(0.2, (h0 - TOP - BOTTOM) / n), ROW_MAX);
  if (S.fit) { S.row = rowMin; }
  var row = Math.max(rowMin, Math.min(ROW_MAX, S.row)), k, v, d;
  S.row = row;
  var font = Math.max(9, Math.min(13, row - 2)), showLabels = row >= 8, hl = marked(g);
  var labelW = 0;
  if (showLabels) {
    g.leaves.forEach(function (l) { labelW = Math.max(labelW, textW(leafText(l), font, l === SEL ? 700 : 400)); });
  } else {
    hl.forEach(function (l) { labelW = Math.max(labelW, textW(leafText(l), 11, l === SEL ? 700 : 400)); });
  }
  var cladeW = 0;
  g.runs.forEach(function (r) {
    if ((r.k1 - r.k0 + 1) * row >= 12) { cladeW = Math.max(cladeW, textW(r.clade, 11, 400)); }
  });
  var callout = showLabels ? 0 : (hl.length ? CALLOUT : 0);
  var right = GAP + STRIP + GAP + callout + (labelW ? labelW + GAP : 0) + (cladeW ? cladeW + GAP : 0) + MARGIN;
  var treeW = Math.max(160, w0 * S.wide - LEFT - right);
  g.w = Math.ceil(LEFT + treeW + right);
  g.h = Math.ceil(TOP + n * row + BOTTOM);
  g.row = row;
  g.unit = treeW / g.maxD;
  for (k = 0; k < n; k++) { Y[g.leaves[k]] = TOP + (k + 0.5) * row; }
  for (v = g.end - 1; v >= g.s; v--) {
    var c = KIDS[v];
    if (c.length) { Y[v] = (Y[c[0]] + Y[c[c.length - 1]]) / 2; }
    X[v] = LEFT + DEP[v] * g.unit;
  }
  var xStrip = LEFT + treeW + GAP, xLabel = xStrip + STRIP + GAP + callout;
  var xClade = xLabel + (labelW ? labelW + GAP : 0);
  g.labelEnd = xLabel + labelW;
  var L = layers(root);
  g.runs.forEach(function (r, j) {
    var y0 = TOP + r.k0 * row, h = (r.k1 - r.k0 + 1) * row;
    el('rect', {x: r1(xStrip), y: r1(y0), width: STRIP, height: r1(Math.max(h, 0.6)), 'class': j % 2 ? 'band-b' : 'band-a'}, L.bands);
    if (h >= 12) { svgText(xClade, y0 + h / 2, r.clade, 'cl', 11, L.clades); }
  });
  if (SEL >= 0 && inSub(SEL, g)) {
    var hs = Math.max(row, 4);
    el('rect', {x: 0, y: r1(Y[SEL] - hs / 2), width: g.w, height: r1(hs), 'class': 'row-sel'}, L.rows);
  }
  d = [];
  for (v = g.s + 1; v < g.end; v++) { d.push(branch(PARENT[v], v)); }
  g.inner.forEach(function (u) { d.push(joint(u)); });
  d.push('M' + r1(X[g.s] - 12) + ',' + r1(Y[g.s]) + 'H' + r1(X[g.s]));
  el('path', {d: d.join(''), 'class': 'br'}, L.branches);
  if (row >= 6) {
    d = [];
    g.inner.forEach(function (u) { d.push(dot(X[u], Y[u], 1.8)); });
    el('path', {d: d.join(''), 'class': 'nd'}, L.dots);
  }
  if (showLabels) {
    d = [];
    g.leaves.forEach(function (l) {
      if (xStrip - X[l] > 10) { d.push('M' + r1(X[l] + 4) + ',' + r1(Y[l]) + 'H' + r1(xStrip - 3)); }
      leafLabel(l, xLabel, Y[l], font, L.labels);
    });
    el('path', {d: d.join(''), 'class': 'ld'}, L.leaders);
  } else {
    var items = hl.map(function (l) { return {v: l, p: Y[l]}; }).sort(function (a, b) { return a.p - b.p; });
    spread(items, 14, TOP + 6, g.h - BOTTOM - 6);
    d = [];
    items.forEach(function (it) {
      d.push('M' + r1(xStrip + STRIP + 2) + ',' + r1(Y[it.v]) + 'L' + r1(xLabel - 3) + ',' + r1(it.p));
      leafLabel(it.v, xLabel, it.p, 11, L.labels);
    });
    el('path', {d: d.join(''), 'class': 'ld-hl'}, L.leaders);
  }
  if (SEL >= 0 && inSub(SEL, g)) { el('path', {d: pathTo(SEL, g), 'class': 'br-sel'}, L.path); }
  marks(hl, L.marks);
  return L;
}

// ---------------------------------------------------------------- circular
function drawCirc(g, root) {
  var n = g.leaves.length, w0 = Math.max(stage.clientWidth, 320), h0 = Math.max(stage.clientHeight, 200);
  if (S.fit) { S.size = 1; }
  var D = Math.max(360, Math.min(w0, h0) - 8) * S.size, step = 2 * Math.PI / n, hl = marked(g), k, v, d;
  var cladeW = 0, labelAll = 0, labelHl = 0;
  g.runs.forEach(function (r) { cladeW = Math.max(cladeW, textW(shortClade(r.clade), 11, 400)); });
  g.leaves.forEach(function (l) { labelAll = Math.max(labelAll, textW(leafText(l), 11, l === SEL ? 700 : 400)); });
  hl.forEach(function (l) { labelHl = Math.max(labelHl, textW(leafText(l), 11, l === SEL ? 700 : 400)); });
  var ring = GAP + STRIP + 4 + (cladeW ? cladeW + GAP : 0);
  var rAll = D / 2 - MARGIN - labelAll;
  var showLabels = rAll - ring > 60 && rAll * step >= 10, Rt, rL;
  if (showLabels) {
    // every label, radial, outside the ring of the clades
    Rt = Math.max(60, rAll - ring);
    rL = Rt + ring;
    g.w = Math.max(Math.ceil(2 * (rL + labelAll + MARGIN)), w0);
    g.h = Math.max(Math.ceil(2 * (rL + labelAll + MARGIN)), h0);
    g.outer = rL + labelAll;
  } else {
    // the marked strains only, in two columns beside the circle
    Rt = Math.max(60, D / 2 - MARGIN - ring);
    rL = Rt + ring;
    g.w = Math.max(Math.ceil(2 * (rL + (hl.length ? CALLOUT + labelHl : 0) + MARGIN)), w0);
    g.h = Math.max(Math.ceil(2 * (rL + MARGIN)), h0);
    g.outer = rL;
  }
  var rS0 = Rt + GAP, rS1 = rS0 + STRIP, rC = rS1 + 4;
  CX = g.w / 2;
  CY = g.h / 2;
  g.step = step;
  g.unit = Rt / g.maxD;
  for (k = 0; k < n; k++) { ANG[g.leaves[k]] = (k + 0.5) * step; }
  for (v = g.end - 1; v >= g.s; v--) {
    var c = KIDS[v];
    if (c.length) { ANG[v] = (ANG[c[0]] + ANG[c[c.length - 1]]) / 2; }
    RAD[v] = DEP[v] * g.unit;
    X[v] = CX + RAD[v] * Math.sin(ANG[v]);
    Y[v] = CY - RAD[v] * Math.cos(ANG[v]);
  }
  var L = layers(root);
  g.runs.forEach(function (r, j) {
    var a0 = ANG[g.leaves[r.k0]] - step / 2, a1 = ANG[g.leaves[r.k1]] + step / 2;
    el('path', {d: sector(rS0, rS1, a0, a1), 'class': j % 2 ? 'band-b' : 'band-a'}, L.bands);
    if ((a1 - a0) * rC >= 13) {
      var m = (a0 + a1) / 2, x = CX + rC * Math.sin(m), y = CY - rC * Math.cos(m);
      var e = svgText(x, y, shortClade(r.clade), 'cl', 11, L.clades);
      rotate(e, x, y, m);
      var t = document.createElementNS(SVGNS, 'title');
      t.textContent = r.clade;
      e.appendChild(t);
    }
  });
  d = [];
  for (v = g.s + 1; v < g.end; v++) { d.push(branch(PARENT[v], v)); }
  g.inner.forEach(function (u) { d.push(joint(u)); });
  el('path', {d: d.join(''), 'class': 'br'}, L.branches);
  if (rL * step >= 6) {
    d = [];
    g.inner.forEach(function (u) { d.push(dot(X[u], Y[u], 1.8)); });
    el('path', {d: d.join(''), 'class': 'nd'}, L.dots);
  }
  if (showLabels) {
    d = [];
    g.leaves.forEach(function (l) {
      var a = ANG[l];
      if (rS0 - RAD[l] > 10) { d.push('M' + pt(RAD[l] + 4, a) + 'L' + pt(rS0 - 3, a)); }
      leafLabel(l, CX + rL * Math.sin(a), CY - rL * Math.cos(a), 11, L.labels, a);
    });
    el('path', {d: d.join(''), 'class': 'ld'}, L.leaders);
  } else {
    d = [];
    [[0, Math.PI, 1], [Math.PI, 2 * Math.PI + 1, -1]].forEach(function (half) {
      var x = CX + half[2] * (rL + CALLOUT);
      var items = hl.filter(function (l) { return ANG[l] >= half[0] && ANG[l] < half[1]; })
        .map(function (l) { return {v: l, p: CY - (rS1 + 2) * Math.cos(ANG[l])}; })
        .sort(function (a, b) { return a.p - b.p; });
      spread(items, 14, MARGIN, g.h - MARGIN);
      items.forEach(function (it) {
        d.push('M' + pt(rS1 + 2, ANG[it.v]) + 'L' + r1(x - half[2] * 3) + ',' + r1(it.p));
        var e = leafLabel(it.v, x, it.p, 11, L.labels);
        if (half[2] < 0) { e.setAttribute('text-anchor', 'end'); }
      });
    });
    el('path', {d: d.join(''), 'class': 'ld-hl'}, L.leaders);
  }
  if (SEL >= 0 && inSub(SEL, g)) { el('path', {d: pathTo(SEL, g), 'class': 'br-sel'}, L.path); }
  marks(hl, L.marks);
  return L;
}

function marks(hl, parent) {
  hl.slice().reverse().forEach(function (v) {
    el('circle', {cx: r1(X[v]), cy: r1(Y[v]), r: v === SEL ? 5.5 : 4, 'class': v === SEL ? 'mk-sel' : 'mk-rank'}, parent);
  });
}

// ---------------------------------------------------------------- render
function render() {
  var size = [stage.clientWidth, stage.clientHeight], g = subtree(S.sub);
  g.runs = cladeRuns(g.leaves);
  var frag = document.createDocumentFragment();
  g.L = S.layout === 'circ' ? drawCirc(g, frag) : drawRect(g, frag);
  while (svg.firstChild) { svg.removeChild(svg.firstChild); }
  svg.appendChild(frag);
  svg.setAttribute('width', g.w);
  svg.setAttribute('height', g.h);
  svg.setAttribute('viewBox', '0 0 ' + g.w + ' ' + g.h);
  G = g;
  HOVER = null;
  if (PIN >= 0 && inSub(PIN, g)) { pinMark(PIN); }
  chrome();
  drawnFor = size;   // the size the drawing was laid out for
}

// keep the point under (ax, ay) of the stage in place while the drawing changes
function keepPoint(ax, ay, change) {
  var old = G, fx = 0.5, fy = 0.5, u = 0;
  if (old) {
    fx = (stage.scrollLeft + ax) / old.w;
    fy = (stage.scrollTop + ay) / old.h;
    if (S.layout === 'rect') { u = (stage.scrollTop + ay - TOP) / old.row; }
  }
  var oldCx = CX, oldCy = CY, oldUnit = old ? old.unit : 1;
  change();
  render();
  if (!old) { return; }
  if (S.layout === 'rect' && old.s === G.s) {
    stage.scrollTop = TOP + u * G.row - ay;
    stage.scrollLeft = fx * G.w - ax;
  } else if (S.layout === 'circ' && old.s === G.s) {
    var f = G.unit / oldUnit;
    stage.scrollLeft = CX + (stage.scrollLeft + ax - oldCx) * f - ax;
    stage.scrollTop = CY + (stage.scrollTop + ay - oldCy) * f - ay;
  } else {
    stage.scrollLeft = fx * G.w - ax;
    stage.scrollTop = fy * G.h - ay;
  }
}
function zoom(kind, factor, ax, ay) {
  if (ax === undefined) { ax = stage.clientWidth / 2; ay = stage.clientHeight / 2; }
  keepPoint(ax, ay, function () {
    S.fit = false;
    if (kind === 'wide') { S.wide = Math.max(1, Math.min(64, S.wide * factor)); }
    else if (S.layout === 'circ') { S.size = Math.max(1, Math.min(16, S.size * factor)); }
    else { S.row = Math.max(0.2, Math.min(ROW_MAX, S.row * factor)); }
  });
}

function centerOn(v) {
  if (!G || !inSub(v, G)) { return; }
  var cw = stage.clientWidth;
  stage.scrollLeft = S.layout === 'rect' ? Math.min(G.labelEnd + 20 - cw, X[v] - 40) : X[v] - cw / 2;
  stage.scrollTop = Y[v] - stage.clientHeight / 2;
}
function flash(v) {
  if (!G || !inSub(v, G)) { return; }
  var c = el('circle', {cx: r1(X[v]), cy: r1(Y[v]), r: 6, 'class': 'flash'}, G.L.top);
  setTimeout(function () { if (c.parentNode) { c.parentNode.removeChild(c); } }, 3000);
}
function show(sub, target) {
  // a part of the tree (the subtree under node sub), then a strain to center on
  S.sub = sub;
  if (S.layout === 'rect') {
    S.wide = 1;
    S.fit = sub !== 0;
    if (sub === 0) { S.row = 12; }
  } else {
    S.fit = true;
  }
  render();
  if (target !== undefined && target >= 0 && inSub(target, G)) { centerOn(target); flash(target); }
  else { stage.scrollLeft = 0; stage.scrollTop = 0; }
}
function around(v) {
  var a = v;
  while (a > 0 && NLEAVES[a] < 20) { a = PARENT[a]; }
  show(a, v);
}

// ---------------------------------------------------------------- chrome
function chrome() {
  var g = G;
  document.querySelectorAll('[data-layout]').forEach(function (b) {
    b.setAttribute('aria-pressed', String(b.getAttribute('data-layout') === S.layout));
  });
  document.querySelectorAll('.rect-only').forEach(function (e) { e.hidden = S.layout !== 'rect'; });
  document.getElementById('z-label').textContent = S.layout === 'rect' ? 'Height' : 'Size';
  document.getElementById('b-up').disabled = g.s === 0;
  document.getElementById('b-whole').disabled = g.s === 0;
  var near = document.getElementById('b-around');
  near.disabled = SEL < 0 && BEST < 0;
  near.textContent = SEL >= 0 ? 'Around the selected strain' : 'Around the best ranked strain';
  var len = niceLength(110 / g.unit);
  document.getElementById('scale-bar').style.width = r1(len * g.unit) + 'px';
  document.getElementById('scale-text').textContent = fmtLength(len) + (T.n_sites ?
    ' = about ' + Math.round(len * T.n_sites).toLocaleString('en-US') + ' differing panel sites' :
    ' of the panel sites differing');
  var crumb = g.s === 0 ? 'Whole tree, ' + NLEAVES[0] + ' strains' :
    'Part of the tree: ' + NLEAVES[g.s] + ' of the ' + NLEAVES[0] + ' strains';
  if (SEL >= 0 && !inSub(SEL, g)) { crumb += ', without the selected strain'; }
  document.getElementById('crumb').textContent = crumb;
  svg.setAttribute('aria-label', crumb + '. Drag to move, Ctrl and the mouse wheel to zoom, ' +
                   'a click on a strain for its details, on a node for its part of the tree.');
}

function fillLegend() {
  document.title = 'Tree of the panel strains' + (T.sample ? ', ' + T.sample : '');
  var sub = 'Sample ' + (T.sample || 'unknown') + ' · ' + NLEAVES[0] + ' strains' +
            (T.n_sites ? ', ' + Number(T.n_sites).toLocaleString('en-US') + ' panel sites' : '');
  document.getElementById('t-sub').textContent = sub;
  var ks = document.getElementById('key-sel');
  if (SEL >= 0) { ks.textContent = 'Selected strain: ' + label(SEL); }
  else if (T.selected) { ks.textContent = 'Selected strain ' + T.selected + ': not in this tree'; }
  else { ks.textContent = T.unassigned ? 'No strain selected: the sample is unassigned' : 'No strain selected'; }
  if (SEL < 0) { ks.previousElementSibling.style.display = 'none'; }
  var nr = Object.keys(RANK).length;
  document.getElementById('key-rank').textContent = nr ?
    'the ' + nr + ' best ranked strains of the sample (#rank)' : 'no ranked strain in this tree';
  document.getElementById('c-ranked').disabled = !nr;
  var about = document.getElementById('about');
  (T.about || []).forEach(function (para) {
    var p = document.createElement('p');
    p.textContent = para;
    about.appendChild(p);
  });
  var list = document.getElementById('strain-names'), names = [];
  Object.keys(LEAF).forEach(function (name) { names.push(label(LEAF[name])); });
  names.sort(function (a, b) { return a.toLowerCase() < b.toLowerCase() ? -1 : 1; });
  names.forEach(function (name) {
    var o = document.createElement('option');
    o.value = name;
    list.appendChild(o);
  });
}

// ---------------------------------------------------------------- hover, details
function line(parent, strong, rest, cls) {
  var p = document.createElement(parent === tip ? 'div' : 'p');
  if (cls) { p.className = cls; }
  if (strong) {
    var b = document.createElement('b');
    b.textContent = strong;
    p.appendChild(b);
  }
  if (rest) { p.appendChild(document.createTextNode(rest)); }
  parent.appendChild(p);
  return p;
}
function keyLine(parent, cls, txt) {
  var p = line(parent, '', '', ''), k = document.createElement('span');
  k.className = 'key-line';
  k.style.borderColor = cls === 'sel' ? 'var(--sel)' : 'var(--rank)';
  p.appendChild(k);
  p.appendChild(document.createTextNode(txt));
}
function facts(v) {
  return [['Clade', cladeOf(v) || 'none'], ['Origin', meta(v, 'origin')], ['Ploidy', meta(v, 'ploidy')],
          ['Zygosity', meta(v, 'zygosity')], ['Sequences', meta(v, 'source')]]
    .filter(function (f) { return f[1]; });
}
function strainLines(parent, v) {
  if (v === SEL) { keyLine(parent, 'sel', 'Selected strain: the strain called for the sample'); }
  if (RANK[v]) {
    keyLine(parent, 'rank', 'Rank ' + RANK[v].rank + ' of the sample ranking: concordance ' +
            RANK[v].conc + ', PP ' + RANK[v].pp);
  }
  facts(v).forEach(function (f) { line(parent, '', f[0] + ': ' + f[1]); });
}
function cladeCounts(v) {
  var count = {}, end = v + NNODES[v], u, c;
  for (u = v; u < end; u++) {
    if (!KIDS[u].length) { c = cladeOf(u) || 'no clade'; count[c] = (count[c] || 0) + 1; }
  }
  var keys = Object.keys(count).sort(function (a, b) { return count[b] - count[a]; });
  return keys.slice(0, 4).map(function (k) { return k + ' ' + count[k]; }).join(', ') +
         (keys.length > 4 ? ', ...' : '');
}

function locate(e) {
  var t = e.target, v;
  if (t && t.getAttribute) {
    var lab = t.closest ? t.closest('[data-v]') : null;
    if (lab) { return {leaf: Number(lab.getAttribute('data-v'))}; }
  }
  var g = G, box = svg.getBoundingClientRect(), px = e.clientX - box.left, py = e.clientY - box.top;
  // the strain whose row (or ray) the pointer is on, right of its tip
  var n = g.leaves.length, k, leaf = -1;
  if (S.layout === 'circ') {
    var dx = px - CX, dy = py - CY, r = Math.sqrt(dx * dx + dy * dy), a = Math.atan2(dx, -dy);
    if (a < 0) { a += 2 * Math.PI; }
    k = Math.floor(a / g.step);
    if (k >= 0 && k < n) { v = g.leaves[k]; if (r >= RAD[v] - 4 && r <= g.outer + 4) { leaf = v; } }
  } else {
    k = Math.floor((py - TOP) / g.row);
    if (k >= 0 && k < n) { v = g.leaves[k]; if (px >= X[v] - 4 && px <= g.labelEnd + 4) { leaf = v; } }
  }
  // the nearest node within 7 px; on a strain's row only when right on it
  // (a strain with a branch of length 0 has its tip on its parent node)
  var best = -1, bd = leaf >= 0 ? 9 : 49;
  g.inner.forEach(function (u) {
    var ux = X[u] - px, uy = Y[u] - py, dd = ux * ux + uy * uy;
    if (dd < bd) { bd = dd; best = u; }
  });
  if (best >= 0) { return {node: best}; }
  return leaf >= 0 ? {leaf: leaf} : null;
}

function clearHover() {
  if (G) { while (G.L.hover.firstChild) { G.L.hover.removeChild(G.L.hover.firstChild); } }
  HOVER = null;
  tip.hidden = true;
}
function hover(hit, e) {
  var key = hit ? (hit.leaf !== undefined ? 'l' + hit.leaf : 'n' + hit.node) : null;
  if (key !== HOVER) {
    clearHover();
    if (!hit) { return; }
    HOVER = key;
    var layer = G.L.hover, v = hit.leaf !== undefined ? hit.leaf : hit.node;
    tip.textContent = '';
    if (hit.leaf !== undefined) {
      if (S.layout === 'rect') {
        var hs = Math.max(G.row, 4);
        el('rect', {x: 0, y: r1(Y[v] - hs / 2), width: G.w, height: r1(hs), 'class': 'row-hov'}, layer);
      }
      el('path', {d: pathTo(v, G), 'class': 'br-hov'}, layer);
      el('circle', {cx: r1(X[v]), cy: r1(Y[v]), r: 6, 'class': 'mk-hov'}, layer);
      line(tip, label(v), label(v).indexOf(NAME[v]) < 0 ? '  ' + NAME[v] : '');
      strainLines(tip, v);
      line(tip, '', 'Click for its links', 'muted');
    } else {
      el('path', {d: subtreePath(v), 'class': 'br-hov'}, layer);
      el('circle', {cx: r1(X[v]), cy: r1(Y[v]), r: 5, 'class': 'mk-hov'}, layer);
      line(tip, NLEAVES[v] + ' strains', '');
      line(tip, '', cladeCounts(v));
      line(tip, '', 'Click to show this part of the tree alone', 'muted');
    }
    tip.hidden = false;
  }
  var w = tip.offsetWidth, h = tip.offsetHeight;
  var vw = document.documentElement.clientWidth, vh = document.documentElement.clientHeight;
  var x = e.clientX + 16, y = e.clientY + 16;
  if (x + w > vw - 8) { x = e.clientX - w - 16; }
  if (y + h > vh - 8) { y = e.clientY - h - 16; }
  tip.style.left = Math.max(8, x) + 'px';
  tip.style.top = Math.max(8, y) + 'px';
}

function pinMark(v) {
  var old = svg.querySelector('.mk-pin');
  if (old) { old.parentNode.removeChild(old); }
  if (v >= 0 && G && inSub(v, G)) {
    el('circle', {cx: r1(X[v]), cy: r1(Y[v]), r: 8, 'class': 'mk-pin'}, G.L.top);
  }
}
function closeInfo() {
  infoBox.hidden = true;
  PIN = -1;
  pinMark(-1);
}
function openInfo(v) {
  PIN = v;
  pinMark(v);
  infoBox.textContent = '';
  var x = document.createElement('button');
  x.type = 'button';
  x.className = 'x';
  x.setAttribute('aria-label', 'Close');
  x.textContent = '×';
  x.addEventListener('click', closeInfo);
  infoBox.appendChild(x);
  var h = document.createElement('h2');
  h.textContent = label(v);
  infoBox.appendChild(h);
  if (label(v) !== NAME[v]) { line(infoBox, '', 'Name in the panel: ' + NAME[v], 'muted'); }
  strainLines(infoBox, v);
  var p = line(infoBox, '', 'Links: '), links = [
    ['SGD', 'https://www.yeastgenome.org/search?q=' + encodeURIComponent(NAME[v])],
    ['PubMed', 'https://pubmed.ncbi.nlm.nih.gov/?term=' + encodeURIComponent('"' + NAME[v] + '" yeast')]];
  String(meta(v, 'runs')).split(',').forEach(function (r) {
    if (/^[SED]RR\d+$/.test(r)) { links.push([r + ' (ENA)', 'https://www.ebi.ac.uk/ena/browser/view/' + r]); }
  });
  links.forEach(function (l, k) {
    var a = document.createElement('a');
    a.href = l[1];
    a.target = '_blank';
    a.rel = 'noopener noreferrer';
    a.textContent = l[0];
    if (k) { p.appendChild(document.createTextNode(' · ')); }
    p.appendChild(a);
  });
  var acts = document.createElement('div');
  acts.className = 'acts';
  var b = document.createElement('button');
  b.type = 'button';
  b.textContent = 'Show the part of the tree around it';
  b.addEventListener('click', function () { around(v); });
  acts.appendChild(b);
  infoBox.appendChild(acts);
  infoBox.hidden = false;
}

// ---------------------------------------------------------------- find
var lastFind = {q: '', t: 0};
function findStrain(q) {
  var msg = document.getElementById('find-msg');
  q = String(q || '').trim().toLowerCase();
  msg.textContent = '';
  if (!q || (q === lastFind.q && Date.now() - lastFind.t < 600)) { return; }
  lastFind = {q: q, t: Date.now()};
  var hit = -1, part = -1;
  Object.keys(LEAF).some(function (name) {
    var v = LEAF[name], lab = label(v).toLowerCase(), id = name.toLowerCase();
    if (id === q || lab === q || lab.indexOf('(' + q + ')') >= 0) { hit = v; return true; }
    if (part < 0 && (id.indexOf(q) >= 0 || lab.indexOf(q) >= 0)) { part = v; }
    return false;
  });
  if (hit < 0) { hit = part; }
  if (hit < 0) { msg.textContent = 'not in the tree'; return; }
  if (!inSub(hit, G)) { S.sub = 0; }
  if (S.layout === 'rect' && S.row < 8) { S.fit = false; S.row = 12; }
  render();
  centerOn(hit);
  flash(hit);
  openInfo(hit);
}

// ---------------------------------------------------------------- export
var SVG_CSS = [
  '#tree text { font-family: system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }',
  '.br { fill: none; stroke: var(--branch); stroke-width: 1; }',
  '.br-sel { fill: none; stroke: var(--sel); stroke-width: 2.5; stroke-linecap: round; }',
  '.br-hov { fill: none; stroke: var(--hover); stroke-width: 2; }',
  '.ld { fill: none; stroke: var(--leader); stroke-width: 1; }',
  '.ld-hl { fill: none; stroke: var(--fg-muted); stroke-width: 1; }',
  '.nd { fill: var(--node); }',
  '.lb { fill: var(--fg); }',
  '.lb-sel { font-weight: 700; }',
  '.rk { fill: var(--fg-muted); font-weight: 400; }',
  '.cl { fill: var(--fg-muted); }',
  '.band-a { fill: var(--band-a); }',
  '.band-b { fill: var(--band-b); }',
  '.mk-sel { fill: var(--sel); stroke: var(--bg); stroke-width: 2; }',
  '.mk-rank { fill: var(--rank); stroke: var(--bg); stroke-width: 2; }',
  '.mk-hov { fill: none; stroke: var(--hover); stroke-width: 1.5; }',
  '.mk-pin { fill: none; stroke: var(--sel); stroke-width: 2; }',
  '.row-sel { fill: var(--row-sel); }',
  '.row-hov { fill: var(--row-hov); }',
  '.flash { fill: none; stroke: var(--sel); stroke-width: 2; transform-box: fill-box; transform-origin: center; animation: flash 0.9s ease-out 3; }',
  '@keyframes flash { from { transform: scale(1); opacity: 1; } to { transform: scale(3.5); opacity: 0; } }',
  '@media (prefers-reduced-motion: reduce) { .flash { animation: none; } }'
].join('\n');

function svgExport() {
  var copy = svg.cloneNode(true), cs = window.getComputedStyle(document.documentElement);
  copy.setAttribute('xmlns', SVGNS);
  ['.layer-hover', '.flash', '.mk-pin'].forEach(function (sel) {
    var e = copy.querySelectorAll(sel);
    for (var k = 0; k < e.length; k++) { e[k].parentNode.removeChild(e[k]); }
  });
  var css = SVG_CSS.split('\n').filter(function (r) { return r.indexOf('flash') < 0 && r.indexOf('@media') < 0; })
    .join('\n').replace(/var\((--[a-z0-9-]+)\)/g, function (m, name) { return cs.getPropertyValue(name).trim(); });
  var style = document.createElementNS(SVGNS, 'style');
  style.textContent = css.replace('#tree text', 'text');
  var bg = document.createElementNS(SVGNS, 'rect');
  bg.setAttribute('width', '100%');
  bg.setAttribute('height', '100%');
  bg.setAttribute('fill', cs.getPropertyValue('--bg').trim());
  copy.insertBefore(bg, copy.firstChild);
  copy.insertBefore(style, copy.firstChild);
  return new XMLSerializer().serializeToString(copy);
}

// ---------------------------------------------------------------- events
var drag = null, noClick = false;
stage.addEventListener('pointerdown', function (e) {
  var box = stage.getBoundingClientRect();
  if (e.button !== 0 || e.pointerType !== 'mouse' || e.clientX - box.left >= stage.clientWidth ||
      e.clientY - box.top >= stage.clientHeight) { return; }
  drag = {x: e.clientX, y: e.clientY, sl: stage.scrollLeft, st: stage.scrollTop, moved: false};
});
window.addEventListener('pointermove', function (e) {
  if (!drag) { return; }
  var dx = e.clientX - drag.x, dy = e.clientY - drag.y;
  if (!drag.moved && Math.abs(dx) + Math.abs(dy) > 4) { drag.moved = true; stage.classList.add('dragging'); clearHover(); }
  if (drag.moved) { stage.scrollLeft = drag.sl - dx; stage.scrollTop = drag.st - dy; }
});
window.addEventListener('pointerup', function () {
  if (drag && drag.moved) { noClick = true; setTimeout(function () { noClick = false; }, 0); }
  drag = null;
  stage.classList.remove('dragging');
});
svg.addEventListener('pointermove', function (e) { if (!drag || !drag.moved) { hover(locate(e), e); } });
svg.addEventListener('pointerleave', clearHover);
svg.addEventListener('click', function (e) {
  if (noClick) { return; }
  var hit = locate(e);
  if (!hit) { closeInfo(); return; }
  if (hit.leaf !== undefined) { openInfo(hit.leaf); }
  else { clearHover(); show(hit.node, PIN >= 0 && inSub(PIN, {s: hit.node, end: hit.node + NNODES[hit.node]}) ? PIN : SEL); }
});
stage.addEventListener('wheel', function (e) {
  if (!e.ctrlKey) { return; }
  e.preventDefault();
  var box = stage.getBoundingClientRect();
  zoom('main', e.deltaY < 0 ? 1.25 : 0.8, e.clientX - box.left, e.clientY - box.top);
}, {passive: false});
document.addEventListener('keydown', function (e) {
  if (e.key === 'Escape') { closeInfo(); clearHover(); }
});
document.querySelectorAll('[data-layout]').forEach(function (b) {
  b.addEventListener('click', function () {
    var layout = b.getAttribute('data-layout');
    if (layout === S.layout) { return; }
    S.layout = layout;
    S.fit = layout === 'circ';
    S.wide = 1;
    S.size = 1;
    if (layout === 'rect') { S.row = 12; }
    render();
    var target = PIN >= 0 ? PIN : (SEL >= 0 ? SEL : BEST);
    if (layout === 'rect' && target >= 0 && inSub(target, G)) { centerOn(target); }
    else { stage.scrollLeft = (G.w - stage.clientWidth) / 2; stage.scrollTop = (G.h - stage.clientHeight) / 2; }
  });
});
document.getElementById('b-whole').addEventListener('click', function () {
  var was = S.sub;
  show(0);
  centerOn(was > 0 ? was : (SEL >= 0 ? SEL : BEST));
});
document.getElementById('b-up').addEventListener('click', function () {
  if (S.sub > 0) { var was = S.sub; show(PARENT[S.sub]); centerOn(was); }
});
document.getElementById('b-around').addEventListener('click', function () {
  var v = SEL >= 0 ? SEL : BEST;
  if (v >= 0) { around(v); }
});
document.getElementById('b-in').addEventListener('click', function () { zoom('main', 1.5); });
document.getElementById('b-out').addEventListener('click', function () { zoom('main', 1 / 1.5); });
document.getElementById('b-wide').addEventListener('click', function () { zoom('wide', 1.5); });
document.getElementById('b-narrow').addEventListener('click', function () { zoom('wide', 1 / 1.5); });
document.getElementById('b-fit').addEventListener('click', function () {
  S.fit = true;
  S.wide = 1;
  S.size = 1;
  render();
  stage.scrollLeft = 0;
  stage.scrollTop = 0;
});
document.getElementById('c-ranked').addEventListener('change', function (e) {
  var top = stage.scrollTop, left = stage.scrollLeft;
  S.ranked = e.target.checked;
  render();
  stage.scrollTop = top;
  stage.scrollLeft = left;
});
var find = document.getElementById('find');
find.addEventListener('change', function () { findStrain(find.value); });
find.addEventListener('input', function () {
  var q = find.value.trim().toLowerCase();
  if (q && Object.keys(LEAF).some(function (name) { return label(LEAF[name]).toLowerCase() === q; })) {
    findStrain(find.value);
  }
});
find.addEventListener('keydown', function (e) { if (e.key === 'Enter') { e.preventDefault(); findStrain(find.value); } });
document.getElementById('dl-nwk').addEventListener('click', function () {
  this.href = 'data:text/plain;charset=utf-8,' + encodeURIComponent(T.newick || '');
});
document.getElementById('dl-svg').addEventListener('click', function () {
  this.href = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgExport());
});
// redrawn when the window or the header (which wraps) changes size
var resizing = 0, drawnFor = [0, 0];
function relayout() {
  clearTimeout(resizing);
  resizing = setTimeout(function () {
    if (Math.abs(stage.clientWidth - drawnFor[0]) < 2 && Math.abs(stage.clientHeight - drawnFor[1]) < 2) { return; }
    keepPoint(stage.clientWidth / 2, stage.clientHeight / 2, function () {});
  }, 120);
}
window.addEventListener('resize', relayout);
if (window.ResizeObserver) { new ResizeObserver(relayout).observe(document.querySelector('header.bar')); }

// ---------------------------------------------------------------- start
var style = document.createElement('style');
style.textContent = SVG_CSS;
document.head.appendChild(style);
fillLegend();
render();
var start = SEL >= 0 ? SEL : BEST;
if (start >= 0) { centerOn(start); flash(start); }
})();
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


# Button next to "Strains in the panel": opens the tree page (openTree() in JS).
TREE_BUTTON = (
    ' <button type="button" class="tree-btn" id="open-tree" '
    'title="Open the tree of the {n} panel strains in a new page, the strain called highlighted" '
    'aria-label="Open the tree of the panel strains in a new page">'
    '<svg class="ico" viewBox="0 0 24 24" width="1em" height="1em" fill="none" stroke="currentColor" '
    'stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">'
    '<path d="M2 12h5M7 6v12M7 6h6M7 18h6M13 3v6M13 3h8M13 9h8M13 15v6M13 15h8M13 21h8"/></svg>'
    'tree</button>')


def render_header(d, title, tree=None):
    panel = d.get("panel") or {}
    rows = [
        ("Organism", d.get("organism"), "organism"),
        ("Reference genome", d.get("reference"), "reference"),
        ("Date", fmt_date(d.get("generated")), "date"),
        ("SNP sites in the panel", num(panel.get("n_sites")), "panel_sites"),
        ("Strains in the panel", num(panel.get("n_strains")), "panel_strains"),
    ]
    extra = {"panel_strains": TREE_BUTTON.format(n=esc(tree["n_leaves"]))} if tree else {}
    items = "\n".join(
        "      <div><dt>{} {}</dt><dd>{}{}</dd></div>".format(
            esc(k), help_btn(h, k.lower()), esc(v) if v not in (None, "") else "n/a",
            extra.get(h, ""))
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


def snp_position(row, best, shown):
    """Where a strain close to the best one stands in the SNP panel: rank and
    PP as in Strain ranking (saying when it is not among the strains shown
    there), log-likelihood gap, separating sites and the allele the sample
    carries at them."""
    first = not row.get("rank") or row["rank"] == 1
    parts = ["first of the SNP ranking" if first else "SNP rank {}".format(row["rank"])]
    if row.get("strain") not in shown:
        parts.append("not in the ranking below")
    if row.get("probability") is not None:   # absent from older result files
        parts.append("PP " + pct(row["probability"]))
    if first:
        return ", ".join(parts)
    if row.get("log_likelihood_gap") is not None:
        parts.append("log-likelihood {:.1f} below".format(row["log_likelihood_gap"]))
    separating = row.get("separating_sites")
    if separating is None:
        return ", ".join(parts)
    if not separating:
        sites = "no separating site"
    else:
        sites = "{} separating site{}".format(num(separating), "s" if separating > 1 else "")
        with_best = row.get("separating_with_best")
        if with_best is not None:
            # The allele the sample carries at these sites: that of the best
            # strain, that of this one (or neither, where one is heterozygous).
            if with_best == separating:
                sites += ", sample with the {} allele {}".format(
                    best, "there" if separating == 1 else "at both" if separating == 2 else
                    "at all {}".format(num(separating)))
            else:
                sites += ", sample with the {} allele at {}".format(best, num(with_best))
                if row.get("separating_with_strain"):
                    sites += ", with the {} allele at {}".format(
                        row.get("strain"), num(row["separating_with_strain"]))
    return "{}; {}".format(", ".join(parts), sites)


def render_decided_by(d):
    """How the strain was called: the SNP panel alone, or the SNP panel then
    the mating type and the native markers among the strains too close to
    the best one for the SNPs alone, each listed with where it stands in the
    SNP panel, what agrees and what contradicts it."""
    res = d.get("resolution") or {}
    rows = res.get("candidates") or []
    if "method" not in res:   # result file older than the two-step call
        return ""
    if not rows:
        text = "SNP panel alone: no other strain close to {}".format(esc(d.get("call")))
        return '<li><strong>Decided by:</strong> {} {}</li>'.format(
            text, help_btn("decided_by", "how the strain is decided"))
    best = next((row.get("strain") for row in rows if row.get("rank") == 1), rows[0].get("strain"))
    # The two rules that make a strain close, with the thresholds of the run
    # (the defaults for a result file older than these fields).
    tie = res.get("tie_log_ratio", math.log(1000))
    rules = ("log-likelihood at most ln {:g} = {:.1f} below it, or fewer than {} separating "
             "sites".format(round(math.exp(tie), 2), tie, res.get("min_separating_sites", 200)))
    close = len(rows) - 1
    if res.get("method") == "snp+markers":
        text = ("SNP panel, then the mating type and the native markers among {b} and the {n} "
                "too close to it for the SNPs alone ({r})").format(
            b=esc(best), n="{} strains".format(close) if close > 1 else "strain", r=rules)
    else:
        text = ("SNP panel: {n} strain{s} too close to {b} for the SNPs alone ({r}), which the "
                "mating type and the native markers do not settle; the first of the SNP ranking "
                "is kept").format(b=esc(best), n=close,
                                  s=" is" if close == 1 else "s are", r=rules)
    shown = {s.get("strain_id") for _, s in ranked(d.get("strains") or [])[:N_STRAINS_SHOWN]}
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
        items.append('    <li><span class="sid">{s}</span> {o} <span class="muted">({snp})'
                     '</span>: {p}</li>'.format(s=esc(row.get("strain")),
                                                o=OUTCOME.get(row.get("outcome"), ""),
                                                snp=esc(snp_position(row, best, shown)),
                                                p="; ".join(points)))
    return ('<li><strong>Decided by:</strong> {t} {h}<ul>\n{i}\n  </ul></li>').format(
        t=text, h=help_btn("decided_by", "how the strain is decided"), i="\n".join(items))


def marker_token(record):
    """Genotype name of one marker as the call prints it: URA3 when active,
    ura3-1 (or ura3[stop gained]) when inactive."""
    if not record:
        return ""
    if record.get("status") == "no_data":
        return "{}(?)".format(record.get("gene"))
    alleles = record.get("alleles") or [{}]
    allele = alleles[0]
    if record.get("status") == "active":
        return str(record.get("gene")).upper()
    if record.get("status") == "inactive":
        known = allele.get("known_as")
        if known:
            return "/".join(sorted({k.split(" ")[0] for k in known}))
        return "{}[{}]".format(str(record.get("gene")).lower(),
                               (allele.get("effect") or "").replace("_", " "))
    return "{}({})".format(str(record.get("gene")).lower(), record.get("status"))


def marker_state(status, record):
    """Genotype of one marker and its status, in plain text."""
    label = MARKER_STATUS.get(status, (None, None, status or "?"))[2]
    return '<code>{t}</code> {l}'.format(t=esc(marker_token(record) or status or ""),
                                          l=esc(label))


def render_decision(d):
    """The strain called, how it was decided (SNP panel, then mating type and
    markers among the strains too close for the SNPs alone), then whether the
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
        records = {m.get("gene"): m for m in (d.get("markers") or {}).get("markers") or []}
        # One row per marker, "URA3:" then genotype and verdict in aligned columns.
        markers = ('<li><strong>Markers:</strong> {h}<div class="marker-grid">\n{}\n  </div></li>').format(
            "\n".join(
                '    <span>{g}:</span><span>{o}</span><span>{v}</span>'.format(
                    g=esc(g.get("gene")),
                    o=marker_state(g.get("observed"), records.get(g.get("gene"))),
                    v=verdict("agree" if g.get("agree") is True else
                              "disagree" if g.get("agree") is False else None))
                for g in genes),
            h=help_btn("markers_agree", "how the markers are compared"))
    mito = (d.get("markers") or {}).get("mitochondria")
    glycerol = ('<li><strong>Glycerol:</strong> {} {}</li>'.format(
        mito_state(mito), help_btn("mito", "the mtDNA status")) if mito else "")
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
            '    <li>{n} <span class="muted">{s}</span></li>'.format(
                n=esc(e.get("name")), s=esc(e.get("status")))
            for e in elements))
    else:
        tools = '<li><strong>Tools:</strong> no foreign element found</li>'
    return """<section id="decision">
  <h2>Decision</h2>
  <ul class="decision">
  <li><strong>Strain:</strong> <span class="sid">{strain}</span></li>
  {decided}
  {mating}
  {markers}
  {glycerol}
  {tools}
  </ul>
</section>""".format(strain=strain, decided=render_decided_by(d), mating=mating,
                    markers=markers, glycerol=glycerol, tools=tools)


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


def ranked(strains):
    """(index, strain) pairs in rank order, the strains without a rank last:
    the order of Strain ranking, whose first N_STRAINS_SHOWN are shown."""
    return sorted(
        enumerate(strains),
        key=lambda t: (t[1].get("rank") if isinstance(t[1].get("rank"), (int, float)) else 10 ** 9, t[0]),
    )


def render_strains(d):
    strains = d.get("strains") or []
    if not strains:
        return """<section id="strains"><h2>Strain ranking</h2>
  <p class="empty">No strain was scored for this sample.</p></section>"""

    ordered = ranked(strains)

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
        <span class="sid">{sid}</span> {h_strain}
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
            h_strain=('<button type="button" class="help" data-help="strain" data-si="{si}" '
                      'aria-label="Explain strain" title="Explain strain">?</button>').format(si=si),
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
        rows.append(("<span style=\"text-transform:none\">Fraction a / (a + &alpha;)</span> " +
                     help_btn("a_fraction", "the fraction a / (a + alpha)"),
                     esc(num(mating["a_fraction"]))))
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
        <td class="status {c}"><span class="g">{glyph}</span>{label}</td>
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
    # A native element named like a marker gene (URA3, TRP1...) is the genomic
    # gene itself, judged in the Native marker genes section: not listed here.
    marker_genes = {str(m.get("gene")).upper()
                    for m in (d.get("markers") or {}).get("markers") or []}
    kept = [t for t in tools
            if not (t.get("origin") == "native"
                    and str(t.get("name") or "").upper() in marker_genes)]
    if len(kept) != len(tools):
        summary = {c: sum(1 for t in kept if (t.get("scorable") or "yes") != "short"
                          and colour_of(t) == c) for c in ("green", "orange", "red")}
        tools = kept

    tiles = "\n".join(
        """    <div class="tile"><div class="k"><span class="pill {c}">{glyph} {label}</span></div>
      <div class="v">{v}</div></div>""".format(
            c=c, glyph=COLOUR_GLYPH[c], label=COLOUR_LABEL[c], v=num(summary.get(c, 0)))
        for c in ("green", "orange", "red")
    )

    legend = """  <div class="legend">
    <strong>Colour / status legend</strong>
    <ul>
      <li><span class="status green"><span class="g">{g}</span>present</span>
        &mdash; score &ge; 99 (element present)</li>
      <li><span class="status orange"><span class="g">{o}</span>partial</span>
        &mdash; 90 &le; score &lt; 99 (partial or diverged)</li>
      <li><span class="status red"><span class="g">{r}</span>absent</span>
        &mdash; score &lt; 90 (absent or too little evidence)</li>
    </ul>
  </div>""".format(g=COLOUR_GLYPH["green"], o=COLOUR_GLYPH["orange"], r=COLOUR_GLYPH["red"])

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


MITO_STATUS = {
    # status -> (colour class, glyph); never hue alone
    "rho+": ("green", "&#9679;"),
    "rho-": ("red", "&#9632;"),
    "rho0": ("red", "&#9632;"),
}


def mito_state(mito):
    """mtDNA status and the growth on glycerol it implies."""
    status = mito.get("status") or "?"
    colour, glyph = MITO_STATUS.get(status, ("orange", "&#9650;"))
    return ('<span class="status {c}"><span class="g">{g}</span>{gl}</span> '
            '<span class="muted">mtDNA {s}</span>').format(
        c=colour, g=glyph, gl=esc(mito.get("glycerol") or "unknown"),
        s=esc(status.replace("_", " ")))


def render_mito(d):
    """Mitochondrial genome: rho+, rho- or rho0, hence growth on glycerol."""
    mito = (d.get("markers") or {}).get("mitochondria")
    if not mito:
        return ""

    def depth(value):
        return "n/a" if value is None else "{:.1f}".format(float(value))
    rows = [("Growth on glycerol " + help_btn("mito", "the mtDNA status"), mito_state(mito)),
            ("Reason", esc(mito.get("reason") or ""))]
    if mito.get("nuclear_depth") is not None:
        rows.append(("Nuclear depth", esc(depth(mito["nuclear_depth"])) + "x"))
    if mito.get("mito_median_depth") is not None:
        rows.append(("mtDNA gene depth (median)", esc(depth(mito["mito_median_depth"])) + "x"))
    if mito.get("mtdna_copies_per_nuclear_genome") is not None:
        rows.append(("mtDNA copies per nuclear genome",
                     esc(depth(mito["mtdna_copies_per_nuclear_genome"]))))
    items = "\n".join("    <div><dt>{}</dt><dd>{}</dd></div>".format(k, v) for k, v in rows)
    genes = mito.get("genes") or []
    state_colour = {"present": ("green", "&#9679;"), "absent": ("red", "&#9632;"),
                    "partial": ("orange", "&#9650;")}
    table = ""
    if genes:
        table = """  <table>
    <thead><tr><th>gene</th><th>role</th><th>exon bp</th><th>covered</th><th>median depth</th>
      <th>state</th></tr></thead>
    <tbody>
{}
    </tbody>
  </table>""".format("\n".join(
            "      <tr><td><strong>{g}</strong></td><td>{r}</td><td>{bp}</td><td>{cov}</td>"
            "<td>{dp}</td><td><span class=\"status {c}\"><span class=\"g\">{gl}</span>{st}"
            "</span></td></tr>".format(
                g=esc(x.get("gene")), r=esc(x.get("role")), bp=esc(x.get("exon_bp")),
                cov=esc(pct(x.get("covered_fraction"), 0)), dp=esc(depth(x.get("median_depth"))),
                c=state_colour.get(x.get("state"), ("orange", ""))[0],
                gl=state_colour.get(x.get("state"), ("orange", "&#9650;"))[1],
                st=esc(x.get("state")))
            for x in genes))
    notes = (mito.get("notes") or []) + (mito.get("limits") or [])
    note_items = "\n".join("    <li>{}</li>".format(esc(x)) for x in notes)
    return """<section id="mito">
  <h2>Mitochondrial genome (growth on glycerol)</h2>
  <dl class="meta">
{items}
  </dl>
{table}
  <h3>Notes and limits</h3>
  <ul>
{notes}
  </ul>
</section>""".format(items=items, table=table, notes=note_items)


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


def documented_match(m, n_markers):
    """One documented strain against the sample: every marker with its
    documented allele, then the counts."""
    genes = m.get("genes")
    if genes is None:   # result file older than the per-gene list
        rows = "".join('    <span>{}</span><span></span><span></span>\n'.format(esc(x))
                       for x in m.get("differences") or [])
    else:
        cells = []
        for g in genes:
            if g.get("agree") is not None:
                value = "<code>{}</code>".format(esc(g.get("documented")))
                v = verdict("agree" if g["agree"] else "disagree")
            elif g.get("documented") is None:
                value, v = "heterozygous", '<span class="muted">not compared</span>'
            else:
                value, v = "not detected", ""
            cells.append('    <span>{g}:</span><span>{a}</span><span>{v}</span>\n'.format(
                g=esc(g.get("gene")), a=value, v=v))
        rows = "".join(cells)
    return """    <li><strong>{s}</strong>
  <div class="marker-grid">
{rows}  </div>
  {s}: {c} of {t} markers identified,<br>{a} of {c} markers agree.</li>""".format(
        s=esc(m["strain"]), rows=rows, a=m["agree"], c=m["compared"],
        t=m.get("total") or n_markers)


def render_markers(d):
    mk = d.get("markers") or {}
    markers = mk.get("markers") or []
    if not markers:
        return ""
    matches = [m for m in (mk.get("documented_genotype_matches") or [])][:4]
    match_items = "\n".join(documented_match(m, len(markers)) for m in matches)
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
  <p class="note">Loss of function is defined according to the detection of the variation in the
  EFFECT column.</p>
  <p><strong>Genotype identified:</strong> <code>{genotype}</code></p>
{expected_note}
  <table>
    <thead><tr><th>gene</th>{expected_head}<th>effect {effect_help}</th><th>coding changes</th>
      <th>evidence</th><th>status {status_help}</th>{conclusion_head}</tr></thead>
    <tbody>
{rows}
    </tbody>
  </table>
  <h3>Documented genotypes</h3>
  <ul class="documented">
{matches}
  </ul>
  <h3>Limits</h3>
  <ul>
{limits}
  </ul>
</section>""".format(genotype=esc(mk.get("genotype")),
                     expected_note=expected_note, expected_head=expected_head,
                     conclusion_head=conclusion_head,
                     status_help=help_btn("marker_status", "the marker status"),
                     effect_help=help_btn("marker_effect", "the effect"),
                     rows="\n".join(marker_row(m, expected_by_gene.get(m.get("gene"), {})
                                               if ex else None) for m in markers),
                     matches=match_items or "    <li>none documented</li>", limits=limits)


# --------------------------------------------------------------------------- #
# tree of the panel strains
# --------------------------------------------------------------------------- #

def parse_newick(text):
    """Nodes of a Newick tree in preorder, the root first: (parents, lengths,
    names), parent -1 for the root, length 0 where none is given, name '' for
    an unnamed node. Written without recursion: the tree can be deep. Quoted
    labels ('...', '' for a quote) and [comments] are understood."""
    parents, lengths, names = [], [], []
    stack, last, expect, i, n = [], None, True, 0, len(text)

    def new(parent):
        parents.append(parent)
        lengths.append(0.0)
        names.append("")
        return len(parents) - 1

    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
        elif c == "[":
            j = text.find("]", i)
            i = n if j < 0 else j + 1
        elif c == "(":
            node = new(stack[-1] if stack else -1)
            stack.append(node)
            last, expect = None, True
            i += 1
        elif c in ",)":
            if expect:   # an empty leaf, as in (,)
                new(stack[-1] if stack else -1)
            if c == ")":
                if not stack:
                    raise ValueError("unbalanced ')' at character %d" % i)
                last, expect = stack.pop(), False
            else:
                last, expect = None, True
            i += 1
        elif c == ":":
            if expect:
                last, expect = new(stack[-1] if stack else -1), False
            j = i + 1
            while j < n and text[j] not in ",);[" and not text[j].isspace():
                j += 1
            lengths[last] = float(text[i + 1:j])
            i = j
        elif c == ";":
            break
        else:
            if c == "'":
                j, label = i + 1, []
                while j < n:
                    if text[j] == "'":
                        if j + 1 < n and text[j + 1] == "'":
                            label.append("'")
                            j += 2
                            continue
                        break
                    label.append(text[j])
                    j += 1
                label, i = "".join(label), j + 1
            else:
                j = i
                while j < n and text[j] not in "(),:;[" and not text[j].isspace():
                    j += 1
                label, i = text[i:j], j
            if expect:
                last, expect = new(stack[-1] if stack else -1), False
            names[last] = label
    if stack:
        raise ValueError("unbalanced '(': %d clade(s) not closed" % len(stack))
    if not parents:
        raise ValueError("no node")
    return parents, lengths, names


def strain_source(row):
    """How the sequences of a panel strain were obtained, in words, and its
    read runs. The assembly paths of the database build machine are not
    shown, only the collection they come from."""
    kind, source = row.get("source_type") or "", row.get("source") or ""
    if kind == "reads":
        runs = [r for r in source.split(",") if re.match(r"^[SED]RR\d+$", r)]
        return ("reads, runs " + ", ".join(runs)) if runs else "reads", ",".join(runs)
    if kind == "matrix":
        code = (" (isolate " + source + ")") if re.match(r"^[A-Za-z0-9_]+$", source) else ""
        return "1011 Yeast Genomes matrix" + code, ""
    if kind == "assembly":
        low = source.lower()
        for key, collection in (("scrapdb", "ScRAPdb"), ("sgd", "SGD")):
            if key in low:
                return "genome assembly, " + collection, ""
        return "genome assembly", ""
    return kind, ""


TREE_FIELDS = ["display_name", "clade", "origin", "ploidy", "zygosity", "source", "runs"]


def tree_payload(d, newick_text, strain_rows):
    """Data of the tree page: the tree in preorder, the strains of the panel,
    the strain called and the best ranked ones."""
    parents, lengths, names = parse_newick(newick_text)
    children = set(parents)
    leaves = [names[v] for v in range(len(parents)) if v not in children]
    leaf_set = set(leaves)
    strains = {}
    for row in strain_rows:
        source, runs = strain_source(row)
        strains[row.get("strain")] = [row.get("display_name") or "", row.get("clade") or "",
                                      row.get("origin") or "", row.get("ploidy") or "",
                                      row.get("zygosity") or "", source, runs]
    panel = d.get("panel") or {}
    best = [[s.get("strain_id"), s.get("rank"), pct(s.get("probability")), pct(s.get("concordance"), 1)]
            for _, s in ranked(d.get("strains") or [])[:N_STRAINS_SHOWN]
            if isinstance(s.get("rank"), (int, float))]
    missing = sorted(set(strains) - leaf_set)
    about = [
        "Tree of the {n} strains of the panel, the strains the sample is compared with, built "
        "once by BIONJ (the method of Peter et al. 2018 for the 1011 isolates) from their "
        "genotypes at the {sites} panel sites: the distance between two strains is the fraction "
        "of the sites known in both where their genotypes differ, a heterozygous genotype "
        "against a homozygous one counting one half. The selected strain is the one called for "
        "the sample (Decision of the report); the best ranked strains are the first of its "
        "Strain ranking.".format(n=len(leaves), sites=num(panel.get("n_sites"))),
        "The tree shows which strains group together. Its branch lengths are distances over "
        "the panel sites, not genome-wide divergences: the panel keeps the rarest variants of "
        "each strain and a share of the common ones, so a small and very divergent lineage, "
        "such as the Taiwanese wild isolates, looks closer to the others than it is. The root, "
        "on the most divergent lineage known (the Taiwanese isolates, for yeast), only sets "
        "the drawing.",
        "Clades: Table S1 of Peter et al. 2018 (Nature 556:339), for the 1011 isolates and the "
        "S288C family; the strains taken from genome assemblies have none. Two strains at "
        "distance 0 are not separated by any panel site (in the yeast panel, BY4741 and BY4742).",
        "The tree is the panels/panel_tree.nwk file of the strain identity database "
        "(biosentinelle_yeast_db), made by its bin/build_panel_tree.py and described in its "
        "db/README.md; Download Newick saves it, for FigTree or iTOL.",
    ]
    if missing:
        about.append("{} strain(s) of the panel are not in this tree, which was built before "
                     "them: {}. Rebuild it with build_db.nf.".format(
                         len(missing), ", ".join(missing[:20]) + (", ..." if len(missing) > 20
                                                                  else "")))
    return {
        "parent": parents,
        "length": [round(x, 10) for x in lengths],
        "name": names,
        "fields": TREE_FIELDS,
        "strains": {k: v for k, v in strains.items() if k in leaf_set},
        "newick": newick_text.strip() + "\n",
        "sample": d.get("sample"),
        "selected": None if d.get("unassigned") else d.get("call"),
        "unassigned": bool(d.get("unassigned")),
        "ranked": best,
        "n_sites": panel.get("n_sites"),
        "n_leaves": len(leaves),
        "about": about,
    }


def read_tree(newick_path, strains_path, d):
    """Tree data for the report, or None (with a warning on stderr) when the
    tree file is absent or unreadable: the report is then made without it."""
    if not newick_path:
        return None
    try:
        with open(newick_path, "r", encoding="utf-8") as fh:
            text = fh.read()
        rows = []
        if strains_path:
            with open(strains_path, "r", encoding="utf-8", newline="") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
        return tree_payload(d, text, rows)
    except (OSError, ValueError) as err:
        sys.stderr.write("warning: no tree in the report, {} unreadable: {}\n".format(
            newick_path, err))
        return None


# --------------------------------------------------------------------------- #
# document
# --------------------------------------------------------------------------- #

def render(d, title, tree=None):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    help_data = dict(HELP, __classes__=CLASS_HELP)
    if tree:
        help_data["panel_strains"] += (" The <b>tree</b> button opens the tree of all these "
                                       "strains in a new page, the strain called highlighted.")
    help_data["marker_effect"] = (
        "<b>Effect of the sample's allele on the encoded protein</b>, from the changes found in "
        "the reads at the gene (coding changes column): whole ORF deleted, frameshift, stop "
        "gained, start lost, insertion, missense, or none. The protein the sample encodes is "
        "then looked up in an allele catalogue built from {n} genome assemblies of known "
        "strains: a form documented as a laboratory mutant allele is named (ura3-1, ...), a form "
        "carried by active strains is a neutral polymorphism, and a frameshift, stop, start loss "
        "or deletion is a loss of function even when not catalogued. The status column follows "
        "from this effect.").format(
        n=esc((d.get("markers") or {}).get("catalogue_assemblies") or "several"))
    ex = d.get("markers_expected") or {}
    if ex:
        strain = esc(ex.get("strain") or d.get("call"))
        alleles = ", ".join("{} ({})".format(esc(g.get("gene")),
                                             esc(g.get("expected_allele") or g.get("expected")))
                            for g in ex.get("genes") or [])
        help_data["markers_agree"] = (
            "<b>agree / disagree according to {s}.</b> The status measured in the sample "
            "(active or inactive) is compared with the status expected in {s} "
            "({src}). Expected: {a}. Not compared when either status is unknown.").format(
            s=strain, src=esc(ex.get("source")), a=alleles)
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
        render_header(d, title, tree),
        '<div class="wrap">',
        render_unassigned(d),
        render_decision(d),
        render_warnings(d),
        render_strains(d),
        render_sex(d),
        render_markers(d),
        render_mito(d),
        render_tools(d),
        '<footer class="page">BioSentinelle report &mdash; rendered {} by '
        "<code>render_report.py</code>. Self-contained: no network access required."
        "</footer>".format(esc(now)),
        "</div>",
        '<script type="application/json" id="report-data">{}</script>'.format(embed_json(d)),
        '<script type="application/json" id="help-data">{}</script>'.format(embed_json(help_data)),
        "<script>{}</script>".format(JS),
    ]
    if tree:
        # the four blocks of the tree page, read by openTree()
        for block in (TREE_CSS, TREE_BODY, TREE_JS):
            if "</script" in block.lower():
                raise ValueError("a block of the tree page holds '</script'")
        parts += [
            '<script type="text/plain" id="tree-css">{}</script>'.format(TREE_CSS),
            '<script type="text/plain" id="tree-body">{}</script>'.format(TREE_BODY),
            '<script type="application/json" id="tree-data">{}</script>'.format(embed_json(tree)),
            '<script type="text/plain" id="tree-js">{}</script>'.format(TREE_JS),
        ]
    parts += ["</body>", "</html>", ""]
    return "\n".join(parts)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Render a self-contained interactive HTML report from a "
                    "BioSentinelle strain-identification results JSON file.")
    ap.add_argument("--json", required=True, metavar="FILE", help="input results JSON file")
    ap.add_argument("--out", required=True, metavar="FILE", help="output HTML file")
    ap.add_argument("--title", default="BioSentinelle report", help="report title")
    ap.add_argument("--tree", metavar="FILE",
                    help="tree of the panel strains, Newick (panels/panel_tree.nwk of the "
                         "database): opened in a page of its own from the report")
    ap.add_argument("--strains", metavar="FILE",
                    help="panel_strains.tsv of the database: clade, origin and source of the "
                         "strains shown in the tree")
    args = ap.parse_args(argv)

    with open(args.json, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        sys.stderr.write("error: the JSON root must be an object\n")
        return 1

    title = args.title
    tree = read_tree(args.tree, args.strains, data)

    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(render(data, title, tree))

    sys.stderr.write("wrote {}\n".format(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
