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
import csv
import datetime
import json
import math
import os
import sys


def load_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def score_strains(sites, calls, error_rate, max_detail_sites):
    """Return one result block per strain, ranked, with per-site detail."""
    strains = [column[4:] for column in sites[0] if column.startswith("exp_")] if sites else []
    log_match = math.log(1.0 - error_rate)
    log_mismatch = math.log(error_rate)

    results = []
    for strain in strains:
        log_likelihood = 0.0
        matching = mismatching = uncovered = het = 0
        weight_matching = weight_total = 0.0
        detail = []
        for site in sites:
            expected = site["exp_" + strain]
            call = calls.get(site["site_id"])
            observed = call["call"] if call else "?"
            weight = float(site.get("weight", 1.0))

            if expected == "?" or observed == "?":
                uncovered += 1
                matched = None
            elif observed == "H":
                het += 1
                matched = None
            else:
                matched = (observed == expected)
                if matched:
                    matching += 1
                    weight_matching += weight
                    log_likelihood += weight * log_match
                else:
                    mismatching += 1
                    log_likelihood += weight * log_mismatch
                weight_total += weight

            detail.append({
                "site_id": site["site_id"], "chrom": site["chrom"],
                "pos": int(site["pos"]), "ref": site["ref"], "alt": site["alt"],
                "expected": expected, "observed": observed,
                "depth": int(call["depth"]) if call else 0,
                "alt_fraction": float(call["alt_fraction"]) if call else 0.0,
                "match": matched, "weight": weight,
            })

        # A panel of tens of thousands of sites makes the report unwieldy, so
        # the detail list is capped. Taking the first N in file order would
        # quietly drop the mismatches, which are the rows anyone opening the
        # drill-down is looking for. Everything that is not a plain match is
        # kept first, and matches fill whatever room is left.
        detail_total = len(detail)
        if detail_total > max_detail_sites:
            interesting = [row for row in detail if row["match"] is not True]
            plain = [row for row in detail if row["match"] is True]
            detail = interesting[:max_detail_sites]
            room = max_detail_sites - len(detail)
            if room > 0:
                step = max(len(plain) // room, 1)
                detail = detail + plain[::step][:room]
            detail.sort(key=lambda row: (row["chrom"], row["pos"]))

        concordance = weight_matching / weight_total if weight_total else 0.0
        results.append({
            "strain_id": strain, "log_likelihood": log_likelihood,
            "concordance": round(concordance, 4),
            "sites_called": matching + mismatching,
            "sites_matching": matching, "sites_mismatching": mismatching,
            "sites_uncovered": uncovered, "sites_het": het,
            "sites_in_detail": len(detail), "sites_total": detail_total,
            "detail_truncated": detail_total > len(detail),
            "sites": detail,
        })

    # Softmax over log likelihoods with a uniform prior. Shifting by the maximum
    # keeps exp() in range when the likelihoods are far apart, which they are as
    # soon as the panel has a few thousand sites.
    if results:
        best = max(result["log_likelihood"] for result in results)
        weights = [math.exp(result["log_likelihood"] - best) for result in results]
        total = sum(weights)
        for result, weight in zip(results, weights):
            result["probability"] = round(weight / total, 6) if total else 0.0
    results.sort(key=lambda result: (-result["probability"], -result["concordance"]))
    for rank, result in enumerate(results, start=1):
        result["rank"] = rank
    return results


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
        })
    tools.sort(key=lambda tool: (-tool["score"], tool["name"]))
    return tools


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", required=True, help="panel site matrix TSV")
    parser.add_argument("--calls", required=True, help="per-site allele calls TSV")
    parser.add_argument("--tools", help="tool element manifest TSV")
    parser.add_argument("--tool-coverage", help="per-element coverage TSV")
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
    args = parser.parse_args()

    if not 0.0 < args.error_rate < 0.5:
        parser.error("--error-rate must be between 0 and 0.5")

    sites = load_tsv(args.sites)
    calls = {row["site_id"]: row for row in load_tsv(args.calls)}
    if not sites:
        sys.exit("empty site matrix: %s" % args.sites)

    results = score_strains(sites, calls, args.error_rate, args.max_detail_sites)

    tools, tool_summary = [], {"green": 0, "orange": 0, "red": 0}
    if args.tools:
        coverage = load_tsv(args.tool_coverage) if args.tool_coverage else []
        tools = score_tools(load_tsv(args.tools), coverage, 1)
        for tool in tools:
            if tool["scorable"] == "yes":
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
            discriminating = sum(
                1 for site in sites
                if site["exp_" + best["strain_id"]] in ("A", "R")
                and site["exp_" + runner_up["strain_id"]] in ("A", "R")
                and site["exp_" + best["strain_id"]] != site["exp_" + runner_up["strain_id"]]
                and calls.get(site["site_id"], {}).get("call") in ("A", "R"))
            best["discriminating_sites_vs_next"] = discriminating
            best["runner_up"] = runner_up["strain_id"]
            if discriminating < args.min_discriminating_sites:
                warnings.append(
                    "only %d called sites distinguish %s from %s. The ranking "
                    "below separates them, but on so few sites the margin rests "
                    "on the assumption that errors are independent across sites, "
                    "which correlated error (assembly, panel construction, "
                    "culture drift) breaks. Treat this as %s or %s, not as a "
                    "confident call of %s"
                    % (discriminating, best["strain_id"], runner_up["strain_id"],
                       best["strain_id"], runner_up["strain_id"], best["strain_id"]))
            elif runner_up["probability"] > 0.01:
                warnings.append(
                    "%s and %s are not cleanly separated by this panel"
                    % (best["strain_id"], runner_up["strain_id"]))
    if best and best.get("detail_truncated"):
        warnings.append(
            "the per-site tables show %d of %d sites per strain. Every mismatch "
            "and uncovered site is included; matching sites are thinned. The "
            "counts above are over the whole panel, and the complete calls are "
            "in the site_calls TSV"
            % (best["sites_in_detail"], best["sites_total"]))
    if any(site.get("assumed_ref") == "yes" for site in sites):
        warnings.append(
            "some panel sites assume the reference allele where a strain's assembly "
            "gave no call; rebuild the panel with callable regions to remove this")

    report = {
        "sample": args.sample,
        "generated": datetime.datetime.now().replace(microsecond=0).isoformat(),
        "organism": args.organism,
        "reference": args.reference,
        "panel": {"k": len(next(iter(calls.values()), {}).get("site_id", "")) or None,
                  "n_sites": len(sites),
                  "strains": [result["strain_id"] for result in results]},
        "unassigned": unassigned,
        "unassigned_reason": reason,
        "coverage": {"median_depth_called": median_depth,
                     "uncovered_fraction": round(uncovered_fraction, 4)},
        "strains": results,
        "tools": tools,
        "tool_summary": tool_summary,
        "warnings": warnings,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as handle:
        json.dump(report, handle, indent=1)

    if best:
        verdict = "UNASSIGNED" if unassigned else best["strain_id"]
        print("%s: %s (p=%.3f, concordance=%.1f%%, %d/%d sites called)" %
              (args.sample, verdict, best["probability"],
               100 * best["concordance"], best["sites_called"], len(sites)))
    for warning in warnings:
        print("warning: %s" % warning)
    return 0


if __name__ == "__main__":
    sys.exit(main())
