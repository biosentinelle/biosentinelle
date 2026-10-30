#!/usr/bin/env python3
"""
Score a sample against the strain panel and the tool element database.

Takes the per-site allele calls from scan_kmers.py (or from a targeted pileup),
the site matrix from merge_snp_panels.py, and optionally the per-element
coverage of the tool database, and writes the single JSON that the HTML report
is rendered from.

Strain probability
------------------
Each called site votes. A site whose observed allele matches what a strain is
expected to carry supports that strain; a mismatch argues against it. Votes are
weighted by the site weight from the panel and combined as a log likelihood
under a per-site error rate, then normalised across the panel strains.

Two things this number is not:

  * It is not a probability that the sample *is* that strain. It is a posterior
    over the strains in the panel, which are assumed to be exhaustive. If the
    sample is a strain nobody put in the panel, one of the listed strains will
    still come out on top, with a high percentage.

  * It is not a measure of fit. Concordance is. A sample can be 99.9% "most
    like W303" while matching only 60% of the W303 expectation, which means it
    is not W303. That case is caught by --min-concordance and reported as
    unassigned, which is the honest answer.

Detecting a strain the panel has never seen
-------------------------------------------
Low concordance only catches an unknown strain that differs from a panel strain
*at the panel's own sites*. A strain whose variants lie somewhere else entirely
carries the reference allele at nearly every panel site, so it scores as a
near-perfect match to the reference strain. Concordance alone would report it
with full confidence, and would be wrong.

What gives it away is the sites that could not be called. Every panel k-mer is
31 bases wide, so a variant anywhere in those 31 bases stops both probes from
matching and the site reads as uncovered. A sample that really is a panel strain
breaks almost none of those windows, because the panel was built to exclude the
windows its own strains break. A strain carrying thousands of variants the panel
has never seen breaks them constantly.

That signal is only usable when depth is adequate, since low coverage produces
uncovered sites too. So the check is made conditional: an uncovered fraction
above --max-uncovered-fraction counts as evidence of an off-panel strain only
when the sites that *were* called have a median depth of at least
--uncovered-depth-floor. Below that, the uncovered sites are attributed to
shallow sequencing and reported as such instead.

Heterozygous sites are not counted as matches or mismatches for any strain;
they are counted separately, because a pile of them means a diploid, a cross or
a mixed culture, and none of those is a single-strain question.
"""

import argparse
import array
import csv
import datetime
import gzip
import json
import math
import os
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_tsv(path):
    with opener(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


if hasattr(int, "bit_count"):
    def popcount(value):
        return value.bit_count()
else:
    def popcount(value):
        return bin(value).count("1")


def bitset(indices, size):
    """Python int with the given bit positions set."""
    buffer = bytearray((size + 7) // 8)
    for index in indices:
        buffer[index >> 3] |= 1 << (index & 7)
    return int.from_bytes(bytes(buffer), "little")


def site_lists(site, key):
    value = site.get(key, "")
    return [int(x) for x in value.split(",")] if value else []


class Panel:
    """The site table as one bit set per strain and expected state.

    A thousand strains times hundreds of thousands of sites does not fit a
    per-site loop in Python. Each strain becomes three bit sets over the
    sites (alternative, heterozygous, unknown; reference is the rest), the
    sample three more (observed A, R, H), and every count below is a popcount
    of their intersections.
    """

    def __init__(self, sites, strains):
        self.sites, self.strains = sites, strains
        self.size = len(sites)
        self.full = (1 << self.size) - 1
        # Per clade and site, how many members carry the alternative allele,
        # are heterozygous, or are unknown, for scoring the unknown sites of
        # a strain against its relatives. A strain without a clade gets -1.
        names = sorted({strain.get("clade", "") for strain in strains} - {""})
        number = {name: i for i, name in enumerate(names)}
        self.clade_of = [number.get(strain.get("clade", ""), -1) for strain in strains]
        self.clade_size = [self.clade_of.count(i) for i in range(len(names))]
        self.clade_counts = {key: [array.array("H", [0]) * self.size for _ in names]
                             for key in ("alt_strains", "het_strains", "unknown_strains")}
        per = {"alt_strains": {}, "het_strains": {}, "unknown_strains": {}}
        for index, site in enumerate(sites):
            for key, table in per.items():
                counts = self.clade_counts[key]
                for strain in site_lists(site, key):
                    table.setdefault(strain, []).append(index)
                    clade = self.clade_of[strain]
                    if clade >= 0:
                        counts[clade][index] += 1
        self.alt, self.het, self.unknown = [], [], []
        for strain in range(len(strains)):
            self.alt.append(bitset(per["alt_strains"].get(strain, ()), self.size))
            self.het.append(bitset(per["het_strains"].get(strain, ()), self.size))
            self.unknown.append(bitset(per["unknown_strains"].get(strain, ()), self.size))

    def ref(self, strain):
        return self.full ^ (self.alt[strain] | self.het[strain] | self.unknown[strain])

    def expected(self, strain, index):
        bit = 1 << index
        if self.unknown[strain] & bit:
            return "?"
        if self.alt[strain] & bit:
            return "A"
        if self.het[strain] & bit:
            return "H"
        return "R"


def observed_sets(panel, calls):
    position = {site["site_id"]: i for i, site in enumerate(panel.sites)}
    found = {"A": [], "R": [], "H": []}
    for site_id, call in calls.items():
        index = position.get(site_id)
        if index is not None and call["call"] in found:
            found[call["call"]].append(index)
    return {key: bitset(value, panel.size) for key, value in found.items()}


def discriminating(panel, obs, a, b):
    """Sites called A or R where strains a and b expect different alleles."""
    known = (panel.full ^ panel.unknown[a]) & (panel.full ^ panel.unknown[b])
    differ = (panel.alt[a] ^ panel.alt[b]) | (panel.het[a] ^ panel.het[b])
    return popcount(differ & known & (obs["A"] | obs["R"]))


def score_strains(panel, calls, error_rate, max_detail_sites, detail_top):
    """Return one result block per strain, ranked; the top ones with site detail."""
    log_match = math.log(1.0 - error_rate)
    log_mismatch = math.log(error_rate)
    # A homozygous call where the strain is heterozygous is not the same
    # error as the wrong homozygote: a spore, or a stock that lost
    # heterozygosity, gives exactly that. It costs a quarter, not error_rate.
    log_partial = math.log(0.25)
    obs = observed_sets(panel, calls)
    called = obs["A"] | obs["R"] | obs["H"]

    # A site where a strain is unknown cannot be left out of its score:
    # a strain with a gappy assembly would then have fewer chances to
    # mismatch than any other, and win against everything. It is scored
    # instead with the probability a relative of the strain would give the
    # observed call there: the frequency within its clade, or within the
    # whole panel when no member of the clade is known at the site. The
    # panel-wide frequency penalises every gap of a laboratory strain:
    # where the S288C family carries the reference allele and most wild
    # strains do not, it costs the strain about a mismatch, and BY4742,
    # with three times the gaps of BY4741, lost to it while matching the
    # sample better. Even one panel-wide pseudo-strain added to a clade of
    # seven costs a tenth of a nat per gap, enough for the same result.
    n_strains = len(panel.strains)
    panel_freq = {}
    for index, site in enumerate(panel.sites):
        known = n_strains - int(site.get("n_unknown") or 0)
        if known > 0:
            panel_freq[index] = (int(site.get("n_alt") or 0) / float(known),
                                 int(site.get("n_het") or 0) / float(known))
    counts = panel.clade_counts

    def log_unknown(index, clade, call):
        if index not in panel_freq:
            return 0.0
        p_alt, p_het = panel_freq[index]
        known = (panel.clade_size[clade] - counts["unknown_strains"][clade][index]
                 if clade >= 0 else 0)
        if known > 0:
            p_alt = counts["alt_strains"][clade][index] / float(known)
            p_het = counts["het_strains"][clade][index] / float(known)
        p = {"A": p_alt, "H": p_het, "R": 1.0 - p_alt - p_het}[call]
        return math.log(min(max(p, error_rate), 1.0 - error_rate))

    call_of = {}
    for state in ("A", "R", "H"):
        bits = obs[state]
        while bits:
            low = bits & -bits
            call_of[low.bit_length() - 1] = state
            bits ^= low

    results = []
    for strain, info in enumerate(panel.strains):
        alt, het, unknown = panel.alt[strain], panel.het[strain], panel.unknown[strain]
        ref = panel.ref(strain)
        matching = (popcount(obs["A"] & alt) + popcount(obs["R"] & ref)
                    + popcount(obs["H"] & het))
        partial = popcount((obs["A"] | obs["R"]) & het)
        mismatching = (popcount(obs["A"] & ref) + popcount(obs["R"] & alt)
                       + popcount(obs["H"] & (alt | ref)))
        unknown_term = 0.0
        clade = panel.clade_of[strain]
        bits = unknown & called
        while bits:
            low = bits & -bits
            index = low.bit_length() - 1
            unknown_term += log_unknown(index, clade, call_of[index])
            bits ^= low
        mismatching += partial
        sites_called = matching + mismatching
        results.append({
            "strain_index": strain, "strain_id": info["strain"],
            "display_name": info.get("display_name", info["strain"]),
            "clade": info.get("clade", ""), "source_type": info.get("source_type", ""),
            "origin": info.get("origin", ""), "ploidy": info.get("ploidy", ""),
            "zygosity": info.get("zygosity", ""), "source": info.get("source", ""),
            "log_likelihood": (matching * log_match + (mismatching - partial) * log_mismatch
                               + partial * log_partial + unknown_term),
            "concordance": round(matching / float(sites_called), 4) if sites_called else 0.0,
            "sites_called": sites_called,
            "sites_matching": matching, "sites_mismatching": mismatching,
            "sites_uncovered": panel.size - popcount(called & (panel.full ^ unknown)),
            "sites_het": popcount(obs["H"]),
        })

    # Softmax over log likelihoods with a uniform prior, shifted by the
    # maximum so exp() stays in range.
    best = max(result["log_likelihood"] for result in results)
    weights = [math.exp(result["log_likelihood"] - best) for result in results]
    total = sum(weights)
    for result, weight in zip(results, weights):
        result["probability"] = round(weight / total, 6) if total else 0.0
    results.sort(key=lambda result: (-result["log_likelihood"], -result["concordance"]))
    for rank, result in enumerate(results, start=1):
        result["rank"] = rank

    # Site detail for the top strains only: with a thousand strains the full
    # detail would be a report of gigabytes. Mismatches are kept first,
    # everything else fills whatever room is left.
    for result in results[:detail_top]:
        strain = result["strain_index"]
        detail = []
        for index, site in enumerate(panel.sites):
            expected = panel.expected(strain, index)
            call = calls.get(site["site_id"])
            observed = call["call"] if call else "?"
            matched = None if expected == "?" or observed == "?" else observed == expected
            detail.append({
                "site_id": site["site_id"], "chrom": site["chrom"],
                "pos": int(site["pos"]), "ref": site["ref"], "alt": site["alt"],
                "expected": expected, "observed": observed,
                "depth": int(call["depth"]) if call else 0,
                "alt_fraction": float(call["alt_fraction"]) if call else 0.0,
                "match": matched, "weight": 1.0,
            })
        detail_total = len(detail)
        if detail_total > max_detail_sites:
            interesting = [row for row in detail if row["match"] is False]
            plain = [row for row in detail if row["match"] is not False]
            detail = interesting[:max_detail_sites]
            room = max_detail_sites - len(detail)
            if room > 0:
                step = max(len(plain) // room, 1)
                detail = detail + plain[::step][:room]
            detail.sort(key=lambda row: (row["chrom"], row["pos"]))
        result.update({"sites": detail, "sites_in_detail": len(detail),
                       "sites_total": detail_total,
                       "detail_truncated": detail_total > len(detail)})
    return results, obs


def mating_type(path, min_kmers_seen, seq_kind="reads"):
    """Mating type from the Ya and Yalpha k-mer counts.

    Per k-mer mean count, a / (a + alpha): about 2/3 for MATa, 1/3 for
    MATalpha and 1/2 for a MATa/MATalpha diploid, because the silent HMR
    and HML copies add one a and one alpha cassette to every genome.

    In an assembly the mean count is the number of cassette copies. A
    short-read assembly often collapses MAT, HML and HMR, which are nearly
    identical, into one a and one alpha contig: 1 : 1 would then read as a
    diploid. Fewer than the three cassettes of a haploid genome means the
    mating type cannot be read.
    """
    rows = {row["class"]: row for row in load_tsv(path)}
    a, alpha = rows.get("MATa"), rows.get("MATalpha")
    if not a or not alpha:
        return {"call": "unknown", "reason": "no mating-type k-mer counts"}
    seen = min(int(a["kmers_seen"]), int(alpha["kmers_seen"]))
    mean_a = int(a["total_count"]) / float(a["n_kmers"])
    mean_alpha = int(alpha["total_count"]) / float(alpha["n_kmers"])
    out = {"mean_count_a": round(mean_a, 2), "mean_count_alpha": round(mean_alpha, 2)}
    if seen < min_kmers_seen or mean_a + mean_alpha == 0:
        out.update({"call": "unknown", "reason": "too few mating-type k-mers seen"})
        return out
    copies_a, copies_alpha = int(round(mean_a)), int(round(mean_alpha))
    if seq_kind == "assembly" and copies_a + copies_alpha < 3:
        out.update({"call": "unknown",
                    "reason": "the assembly holds %d a and %d alpha cassette(s), fewer than the "
                              "three (MAT, HML, HMR) of a complete genome: the cassettes are "
                              "collapsed or missing" % (copies_a, copies_alpha)})
        return out
    ratio = mean_a / (mean_a + mean_alpha)
    out["a_fraction"] = round(ratio, 3)
    if ratio >= 0.58:
        out["call"] = "MATa"
    elif ratio <= 0.42:
        out["call"] = "MATalpha"
    else:
        out["call"] = "MATa/alpha"
    return out


def documented_mating(path):
    """strain -> MATa, MATalpha or MATa/alpha, from documented_genotypes.tsv."""
    out = {}
    for row in load_tsv(path):
        token = next((t for t in row["documented_genotype"].split()
                      if t.upper().startswith("MAT")), "")
        if not token:
            continue
        if "/" in token:
            out[row["strain"]] = "MATa/alpha"
        elif token.lower() in ("matalpha", "matα"):
            out[row["strain"]] = "MATalpha"
        elif token.lower() == "mata":
            out[row["strain"]] = "MATa"
    return out


def score_tools(elements, coverage, min_breadth_depth):
    """Combine breadth and identity into one 0-100 score and a colour."""
    by_id = {row["element_id"]: row for row in coverage} if coverage else {}
    tools = []
    for element in elements:
        measured = by_id.get(element["element_id"])
        breadth = float(measured["breadth"]) if measured else 0.0
        identity = float(measured.get("identity", 0.0)) if measured else 0.0
        mean_depth = float(measured.get("mean_depth", 0.0)) if measured else 0.0

        # Breadth answers "is the whole element there", identity answers "is it
        # the same sequence". An element can be fully covered by reads that
        # differ at every tenth base, which is a diverged homologue and not the
        # element in the database, so both have to be high for a green call.
        score = round(100.0 * breadth * identity, 2)
        if score >= 99.0:
            colour = "green"
        elif score >= 90.0:
            colour = "orange"
        else:
            colour = "red"
        tools.append({
            "element_id": element["element_id"], "name": element["name"],
            "element_class": element["element_class"],
            "length_bp": int(element["length_bp"]),
            "scorable": element.get("scorable", "yes"),
            "breadth": round(breadth, 4), "identity": round(identity, 4),
            "mean_depth": round(mean_depth, 2), "score": score,
            "colour": colour, "present": colour != "red",
            "origin": element.get("origin", "unknown"),
            "feature_type": element.get("feature_type", ""),
            "source_plasmid": element.get("source_plasmid", ""),
            "seen_in": [x for x in element.get("seen_in", "").split(";") if x],
            "host_kmer_fraction": element.get("host_kmer_fraction", ""),
        })
    tools.sort(key=lambda tool: (-tool["score"], tool["name"]))
    return tools


def documented_markers(path, genes):
    """strain -> {gene: allele} from documented_genotypes.tsv, read as
    genotype_markers.py reads it: the wild-type name for a gene not mentioned,
    None for a diploid heterozygous at it."""
    out = {}
    for row in load_tsv(path):
        tokens = row["documented_genotype"].split()
        alleles = {}
        for gene in genes:
            mutant = [t for t in tokens if t.lower().startswith(gene.lower())
                      and not t.startswith(gene.upper())]
            if mutant and gene.upper() in mutant[0].split("/"):
                alleles[gene] = None
            else:
                alleles[gene] = mutant[0].split("/")[0] if mutant else gene.upper()
        out[row["strain"]] = alleles
    return out


def expected_markers(strain, strain_rows, markers, documented_path, catalogue_path,
                     aliases_path):
    """What the native markers of one panel strain should be, gene by gene,
    next to what the sample shows. Kept apart from the SNP ranking: it says
    whether the two kinds of evidence agree, it does not change the ranking.

    The documented genotype is used when the strain has one. Otherwise the
    marker alleles of its own assembly, found in the catalogue under the
    strain name, an alias, its 1011 matrix code, or as a spore of that
    isolate (AGA_1a for AGA)."""
    genes = [m["gene"] for m in markers.get("markers", [])]
    if not genes:
        return {}
    out = {"strain": strain, "source": "", "genes": [], "agree": 0, "compared": 0}
    documented = documented_markers(documented_path, genes) if documented_path else {}
    if strain in documented:
        out["source"] = "documented genotype"
        expected = {gene: ("heterozygous" if allele is None else
                           ("active" if allele == gene.upper() else "inactive"), allele or "")
                    for gene, allele in documented[strain].items()}
    else:
        expected = {}
        if catalogue_path:
            row = next((r for r in strain_rows if r["strain"] == strain), {})
            aliases = ({r["alias"]: r["strain"] for r in load_tsv(aliases_path)}
                       if aliases_path else {})
            names = {strain, row.get("source", "")}
            names.update(alias for alias, name in aliases.items() if name == strain)
            by_name = {}
            for r in load_tsv(catalogue_path):
                by_name.setdefault(r["strain"], []).append(r)
            assembly = next((n for n in sorted(names) if n in by_name), None)
            if assembly is None:
                assembly = next((n for n in sorted(by_name)
                                 if n.split("_")[0] in names and "_" in n), None)
            if assembly is not None:
                out["source"] = "assembly %s in the marker catalogue" % assembly
                # A missense form is "uncertain" in the catalogue when no
                # documented strain settles it; it is shown, not compared.
                for r in by_name[assembly]:
                    expected[r["gene"]] = (r["activity"] or None,
                                           r.get("documented_allele") or
                                           r["effect"].replace("_", " "))
    if not expected:
        out["source"] = "none: no documented genotype and no assembly in the marker catalogue"
    for m in markers["markers"]:
        status, allele = expected.get(m["gene"], (None, ""))
        observed = m.get("status")
        comparable = status in ("active", "inactive") and observed in ("active", "inactive")
        out["genes"].append({"gene": m["gene"], "expected": status or "unknown",
                             "expected_allele": allele or "", "observed": observed,
                             "agree": (status == observed) if comparable else None})
        if comparable:
            out["compared"] += 1
            out["agree"] += status == observed
    return out


def resolve(panel, obs, results, mating, mating_doc, expect_markers, tie_log_ratio,
            min_separating):
    """The strain call: the SNP panel first, then, when it hesitates, the
    mating type and the native markers.

    Close relatives (S288C, FY4, BY4741, BY4742, ...) differ at a few dozen
    panel sites or none: BY4741 and BY4742 are sister segregants of one cross
    and differ mainly by MATa met15D0 against MATalpha lys2D0, whole-gene
    deletions and a mating type that no SNP site sees. A strain is close to
    the best one when fewer than min_separating called sites separate them,
    or when its likelihood is within tie_log_ratio of the best (about two net
    mismatches at ln 1000). With no close strain, the call is the first of
    the SNP ranking. Otherwise each candidate (the best and the close
    strains) is compared with the measured mating type and native markers,
    from its documented genotype or, without one, from the marker alleles of
    its assembly (expect_markers). A candidate contradicted by any of them is
    ruled out; among the others, the one agreeing at the most points is
    called, the SNP ranking breaking the ties. When every candidate is
    contradicted, or none can be compared, the call stays the first of the
    SNP ranking.
    """
    best = results[0]
    candidates = [(best, 0, 0.0)]
    for result in results[1:60]:
        separating = discriminating(panel, obs, best["strain_index"], result["strain_index"])
        gap = best["log_likelihood"] - result["log_likelihood"]
        if separating < min_separating or gap <= tie_log_ratio:
            candidates.append((result, separating, gap))
    if len(candidates) == 1:
        return {"call": best["strain_id"], "method": "snp", "candidates": [],
                "summary": "%s, first of the SNP ranking, with no other strain close to it"
                           % best["strain_id"],
                "explanation": []}

    rows = []
    for result, separating, gap in candidates:
        strain = result["strain_id"]
        mating_row = check_mating(mating, mating_doc, strain)
        expected = expect_markers(strain)
        agree, against = [], []
        if mating_row["agreement"] == "agree":
            agree.append("mating type %s" % mating_row["measured"])
        elif mating_row["agreement"] == "disagree":
            against.append("mating type %s, documented %s"
                           % (mating_row["measured"], mating_row["documented"]))
        for gene in expected.get("genes", []):
            if gene["agree"] is True:
                agree.append("%s %s" % (gene["gene"], gene["observed"]))
            elif gene["agree"] is False:
                allele = gene["expected_allele"]
                against.append("%s %s, expected %s%s" % (
                    gene["gene"], gene["observed"], gene["expected"],
                    " (%s)" % allele if allele and allele not in (gene["gene"].upper(),
                                                                  "reference") else ""))
        rows.append({"strain": strain, "rank": result["rank"],
                     "separating_sites": separating,
                     "mismatching": result["sites_mismatching"],
                     "log_likelihood_gap": round(gap, 2),
                     "documented_mating": mating_doc.get(strain, ""),
                     "marker_source": expected.get("source", ""),
                     "agree": agree, "disagree": against})

    kept = sorted((row for row in rows if not row["disagree"]),
                  key=lambda row: (-len(row["agree"]), row["rank"]))
    names = ", ".join(row["strain"] for row in rows)
    explanation = []
    if kept and kept[0]["agree"]:
        call, method = kept[0]["strain"], "snp+markers"
        tied = [row["strain"] for row in kept[1:] if len(row["agree"]) == len(kept[0]["agree"])]
        summary = ("the SNP panel does not separate %s; the mating type and the native "
                   "markers point to %s (%s)" % (names, call, ", ".join(kept[0]["agree"])))
        if tied:
            summary += ("; %s agree(s) as well, and %s is called as first of the SNP ranking"
                        % (", ".join(tied), call))
            explanation.append("%s and %s are separated neither by the SNP panel nor by the "
                               "mating type and the native markers" % (call, ", ".join(tied)))
    else:
        call, method = best["strain_id"], "snp"
        if kept:
            ruled_out = [row["strain"] for row in rows if row["disagree"]]
            summary = ("the SNP panel does not separate %s%s; %s ha%s no documented mating "
                       "type or markers to compare with, %s is kept as first of the SNP ranking"
                       % (names,
                          ", the mating type or the native markers rule out %s"
                          % ", ".join(ruled_out) if ruled_out else "",
                          ", ".join(row["strain"] for row in kept),
                          "s" if len(kept) == 1 else "ve", call))
        else:
            summary = ("the SNP panel does not separate %s, and every one of them contradicts "
                       "the mating type or the native markers; %s is kept as first of the SNP "
                       "ranking. The sample may be a derivative of it with another mating type "
                       "or other markers, or a strain outside the panel" % (names, call))
        explanation.append(summary)
    for row in rows:
        row["outcome"] = ("called" if row["strain"] == call else
                          "ruled out" if row["disagree"] else "not excluded")
    return {"call": call, "method": method, "candidates": rows, "summary": summary,
            "explanation": explanation}


def check_mating(mating, mating_doc, strain):
    """Measured mating type against the one documented for the strain called
    by the SNPs: agree, disagree, or not compared (one of them unknown)."""
    documented = mating_doc.get(strain, "")
    measured = mating.get("call", "unknown")
    out = {"strain": strain, "documented": documented or "not documented",
           "measured": measured}
    if not documented or measured in (None, "unknown"):
        out["agreement"] = "not compared"
        out["reason"] = ("no documented mating type for %s" % strain if not documented
                         else "mating type of the sample unknown")
    elif documented == measured:
        out["agreement"] = "agree"
    else:
        out["agreement"] = "disagree"
        if measured == "MATa/alpha":
            out["reason"] = ("the sample reads as a/alpha: a diploid (as BY4743 from BY4741 x "
                             "BY4742), a mix of both mating types, or another strain")
    return out


def check_tools(tools, documented_path, strain):
    """Foreign tool elements found in the sample against the documented
    genotype of the strain called: an element present (green) agrees when the
    genotype inserts it (kanMX in his3D1::kanMX), one inserted but not present
    disagrees, one present but not inserted disagrees. Only what follows '::'
    is an insertion: ura3D0 is a deletion, not a URA3 cassette. Native
    elements (a yeast promoter, TRP1) are left out: every strain carries
    them. Not compared when the strain has no documented genotype."""
    genotype = ""
    if documented_path:
        genotype = next((row["documented_genotype"] for row in load_tsv(documented_path)
                         if row["strain"] == strain), "")
    inserted = " ".join(token.split("::", 1)[1] for token in genotype.lower().split()
                        if "::" in token)
    out = {"strain": strain, "documented": bool(genotype), "elements": []}
    for tool in tools:
        if tool["origin"] == "native" or tool["scorable"] != "yes":
            continue
        found = tool["colour"] == "green"
        expected = bool(inserted) and tool["name"].lower() in inserted
        if not found and not expected:
            continue
        if not genotype:
            agreement = "not compared"
        else:
            agreement = "agree" if found == expected else "disagree"
        out["elements"].append({"name": tool["name"], "element_id": tool["element_id"],
                                "status": "present" if found else "absent",
                                "expected": expected, "agreement": agreement})
    return out


def marker_warnings(markers):
    """What the native marker genotype adds to the warnings."""
    out = []
    best = (markers.get("documented_genotype_matches") or [None])[0]
    if best and best["agree"] == best["compared"]:
        out.append("native marker genotype %s matches the documented genotype of %s at all "
                   "%d markers called" % (markers["genotype"], best["strain"], best["compared"]))
    unsure = [m["gene"] for m in markers["markers"]
              if m["status"] in ("uncertain", "mixed", "no_data")]
    if unsure:
        out.append("native marker(s) %s could not be called as active or inactive; see the "
                   "marker section for why" % ", ".join(unsure))
    out.extend("native markers: %s" % limit for limit in markers.get("limits", []))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", required=True, help="panel site table TSV (.gz)")
    parser.add_argument("--strains", required=True, help="panel strain table TSV")
    parser.add_argument("--mating", help="mating-type k-mer counts from scan_kmers.py")
    parser.add_argument("--documented", help="documented_genotypes.tsv")
    parser.add_argument("--marker-catalogue",
                        help="marker_alleles.tsv: the marker alleles of each assembly, "
                             "expected in a strain that has no documented genotype")
    parser.add_argument("--aliases",
                        help="aliases.tsv: catalogue names of strains the panel names "
                             "otherwise, used to find the expected marker status of the "
                             "best strain")
    parser.add_argument("--tie-log-ratio", type=float, default=math.log(1000),
                        help="strains whose log likelihood is within this of the best "
                             "are not separated by the SNP panel and are settled by the "
                             "mating type and the native markers (default ln 1000, about "
                             "two net mismatching sites)")
    parser.add_argument("--calls", required=True, help="per-site allele calls TSV")
    parser.add_argument("--tools", help="tool element manifest TSV")
    parser.add_argument("--tool-coverage", help="per-element coverage TSV")
    parser.add_argument("--markers", help="marker_genotype.json from genotype_markers.py")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--organism", default="Saccharomyces cerevisiae")
    parser.add_argument("--reference", default="")
    parser.add_argument("--out", required=True, help="output JSON")
    parser.add_argument("--error-rate", type=float, default=0.01,
                        help="per-site probability of an allele call disagreeing with "
                             "the true strain, covering sequencing error, culture "
                             "drift and panel error (default 0.01)")
    parser.add_argument("--min-concordance", type=float, default=0.90,
                        help="below this, the best strain is reported as unassigned")
    parser.add_argument("--min-called-fraction", type=float, default=0.25,
                        help="minimum fraction of the panel that must be called")
    parser.add_argument("--min-discriminating-sites", type=int, default=200,
                        help="if fewer called sites than this separate the top two "
                             "strains, the report says the two cannot be told apart "
                             "confidently, whatever the probability says. Closely "
                             "related strains such as an S288C derivative and S288C "
                             "itself routinely fall below it.")
    parser.add_argument("--max-uncovered-fraction", type=float, default=0.15,
                        help="above this fraction of uncovered panel sites, at adequate "
                             "depth, the sample is reported as a strain outside the "
                             "panel rather than matched to the closest member")
    parser.add_argument("--uncovered-depth-floor", type=float, default=10.0,
                        help="median depth over called sites needed before uncovered "
                             "sites are read as off-panel divergence instead of "
                             "shallow coverage")
    parser.add_argument("--max-detail-sites", type=int, default=5000,
                        help="per-strain site rows written into the report. Every "
                             "mismatch and uncovered site is kept regardless; this "
                             "caps the matching sites, which are the bulk and the "
                             "least interesting. At 5000 a five strain panel gives a "
                             "report of a few MB instead of tens.")
    parser.add_argument("--detail-top", type=int, default=5,
                        help="strains whose per-site detail goes into the report")
    parser.add_argument("--report-top", type=int, default=25,
                        help="strains listed in the report; the full ranking is in "
                             "--out-ranking")
    parser.add_argument("--out-ranking", help="TSV ranking of every panel strain")
    parser.add_argument("--min-mating-kmers", type=int, default=100,
                        help="k-mers of each mating type that must be seen to call it")
    parser.add_argument("--seq-kind", choices=("reads", "assembly"), default="reads",
                        help="reads or assembly: an assembly is checked for collapsed cassettes")
    args = parser.parse_args()

    if not 0.0 < args.error_rate < 0.5:
        parser.error("--error-rate must be between 0 and 0.5")

    sites = load_tsv(args.sites)
    calls = {row["site_id"]: row for row in load_tsv(args.calls)}
    if not sites:
        sys.exit("empty site matrix: %s" % args.sites)
    panel = Panel(sites, load_tsv(args.strains))

    results, obs = score_strains(panel, calls, args.error_rate, args.max_detail_sites,
                                 args.detail_top)

    tools, tool_summary = [], {"green": 0, "orange": 0, "red": 0}
    if args.tools:
        coverage = load_tsv(args.tool_coverage) if args.tool_coverage else []
        tools = score_tools(load_tsv(args.tools), coverage, 1)
        # Only foreign elements are counted. A native element covered at full
        # identity says the strain has its own genome, not that it carries a
        # cassette, so counting it as present would mean every wild-type sample
        # reports several tool elements present.
        for tool in tools:
            if tool["scorable"] == "yes" and tool["origin"] != "native":
                tool_summary[tool["colour"]] += 1

    # Depth over the sites that were called, used to tell an off-panel strain
    # apart from a shallow run. The median is taken rather than the mean so a
    # handful of repeat-collapsed sites cannot carry it.
    called_depths = sorted(int(row["depth"]) for row in calls.values()
                           if row["call"] != "?")
    median_depth = (called_depths[len(called_depths) // 2] if called_depths else 0)
    uncovered_fraction = 1.0 - (len(called_depths) / float(len(sites)))

    warnings = []
    unassigned, reason = False, ""
    best = results[0] if results else None
    if best:
        called_fraction = best["sites_called"] / float(len(sites))
        if (uncovered_fraction > args.max_uncovered_fraction and
                median_depth >= args.uncovered_depth_floor):
            unassigned = True
            reason = ("%.0f%% of the panel sites could not be called although the "
                      "called sites reach a median depth of %d. The sample carries "
                      "variation the panel does not describe, so it is probably a "
                      "strain that is not in the database, whatever the ranking "
                      "below suggests" %
                      (100 * uncovered_fraction, median_depth))
        elif called_fraction < args.min_called_fraction:
            unassigned = True
            reason = ("only %.0f%% of the %d panel sites were called; there is not "
                      "enough data to identify the strain" %
                      (100 * called_fraction, len(sites)))
        elif best["concordance"] < args.min_concordance:
            unassigned = True
            reason = ("the closest strain, %s, matches only %.1f%% of its expected "
                      "alleles; the sample is probably a strain that is not in the "
                      "panel, or a mixture" %
                      (best["strain_id"], 100 * best["concordance"]))
        if best["sites_het"] > 0.05 * max(best["sites_called"], 1):
            warnings.append(
                "%d sites show both alleles: the sample may be diploid, a cross or a "
                "mixed culture, and a single-strain answer may not apply"
                % best["sites_het"])
        if len(results) > 1:
            # How many called sites actually tell the top two apart. The
            # probability is computed over the whole panel, so two strains that
            # differ at 20 sites out of 36000 still produce a likelihood ratio
            # large enough to print as 1.000. That number is only as good as
            # the assumption that site errors are independent, and they are
            # not: a mis-assembled region, a panel built from one assembly, or
            # culture drift all produce correlated errors. Reporting the count
            # lets the reader judge the call instead of trusting the softmax.
            runner_up = results[1]
            separating = discriminating(panel, obs, best["strain_index"],
                                        runner_up["strain_index"])
            best["discriminating_sites_vs_next"] = separating
            best["runner_up"] = runner_up["strain_id"]
            # Fewer separating sites than the threshold is handled by
            # resolve(), which settles them by mating type and markers.
            if separating >= args.min_discriminating_sites and runner_up["probability"] > 0.01:
                warnings.append(
                    "%s and %s are not cleanly separated by this panel"
                    % (best["strain_id"], runner_up["strain_id"]))
    if tools:
        native_hits = [tool["name"] for tool in tools
                       if tool["origin"] == "native" and tool["colour"] != "red"]
        foreign_hits = [tool["name"] for tool in tools
                        if tool["origin"] != "native" and tool["colour"] != "red"
                        and tool["scorable"] == "yes"]
        if native_hits:
            warnings.append(
                "%d element(s) scored as covered but are native to the host "
                "genome (%s). An unmodified strain covers these end to end "
                "because it carries the native locus, so they are not evidence "
                "of a cassette and are excluded from the counts. Showing one of "
                "them as present would need reads spanning the junction with an "
                "insertion site, or a copy number above one"
                % (len(native_hits), ", ".join(sorted(set(native_hits)))))
        if not foreign_hits:
            warnings.append(
                "no foreign tool element was detected: no marker, tag, reporter "
                "or resistance cassette in the database is present in this sample")

    if best and best.get("detail_truncated"):
        warnings.append(
            "the per-site tables show %d of %d sites per strain. Every mismatch "
            "and uncovered site is included; matching sites are thinned. The "
            "counts above are over the whole panel, and the complete calls are "
            "in the site_calls TSV"
            % (best["sites_in_detail"], best["sites_total"]))
    markers = {}
    if args.markers and os.path.getsize(args.markers) > 2:
        with open(args.markers) as handle:
            markers = json.load(handle)
    if markers.get("markers"):
        warnings.extend(marker_warnings(markers))

    mating = (mating_type(args.mating, args.min_mating_kmers, args.seq_kind) if args.mating
              else {"call": "unknown", "reason": "not measured"})
    mating_doc = documented_mating(args.documented) if args.documented else {}

    def expect_markers(strain):
        return (expected_markers(strain, panel.strains, markers, args.documented,
                                 args.marker_catalogue, args.aliases)
                if markers.get("markers") else {})

    resolution = resolve(panel, obs, results, mating, mating_doc, expect_markers,
                         args.tie_log_ratio, args.min_discriminating_sites)
    warnings.extend(resolution["explanation"])

    # The mating type set beside the one documented for the strain called.
    mating_check = check_mating(mating, mating_doc, resolution["call"])
    if mating_check["agreement"] == "disagree":
        warnings.append("mating type disagrees with %s, the strain called: measured %s, "
                        "documented %s%s" % (
                            mating_check["strain"], mating_check["measured"],
                            mating_check["documented"],
                            ". " + mating_check["reason"] if mating_check.get("reason") else ""))

    # The markers expected from the strain called, set beside those measured.
    markers_expected = expect_markers(resolution["call"])
    if markers_expected.get("compared") and \
            markers_expected["agree"] < markers_expected["compared"]:
        warnings.append(
            "native markers differ from those expected of %s, the strain called "
            "(%s), at %d of %d genes: %s" % (
                markers_expected["strain"], markers_expected["source"],
                markers_expected["compared"] - markers_expected["agree"],
                markers_expected["compared"],
                ", ".join("%s %s, expected %s" % (g["gene"], g["observed"], g["expected"])
                          for g in markers_expected["genes"] if g["agree"] is False)))

    report = {
        "sample": args.sample,
        "generated": datetime.datetime.now().replace(microsecond=0).isoformat(),
        "organism": args.organism,
        "reference": args.reference,
        "panel": {"n_sites": len(sites), "n_strains": len(panel.strains),
                  "strains": [result["strain_id"] for result in results[:args.report_top]]},
        "call": resolution["call"],
        "resolution": resolution,
        "mating_type": mating,
        "mating_check": mating_check,
        "tool_check": check_tools(tools, args.documented, resolution["call"]),
        "unassigned": unassigned,
        "unassigned_reason": reason,
        "coverage": {"median_depth_called": median_depth,
                     "uncovered_fraction": round(uncovered_fraction, 4)},
        "strains": results[:args.report_top],
        "tools": tools,
        "tool_summary": tool_summary,
        "warnings": warnings,
        "markers": markers,
        "markers_expected": markers_expected,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as handle:
        json.dump(report, handle, indent=1)
    if args.out_ranking:
        with open(args.out_ranking, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            columns = ["rank", "strain_id", "display_name", "clade", "source_type",
                       "probability", "concordance", "sites_called", "sites_matching",
                       "sites_mismatching", "sites_uncovered"]
            writer.writerow(columns)
            for result in results:
                writer.writerow([result[column] for column in columns])

    if best:
        verdict = "UNASSIGNED" if unassigned else resolution["call"]
        print("%s: %s (by %s; SNP panel best %s, p=%.3f, concordance=%.1f%%, %d/%d sites "
              "called; mating type %s, %s with %s)" %
              (args.sample, verdict, resolution["method"], best["strain_id"],
               best["probability"],
               100 * best["concordance"], best["sites_called"], len(sites),
               mating.get("call", "unknown"), mating_check["agreement"],
               mating_check["documented"]))
    for warning in warnings:
        print("warning: %s" % warning)
    return 0


if __name__ == "__main__":
    sys.exit(main())
