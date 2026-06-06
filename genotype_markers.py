#!/usr/bin/env python3
# =============================================================================
# Génotypage complet des marqueurs génétiques de levure
# Gère : SNV, insertions, délétions partielles et complètes
# =============================================================================

import subprocess
import os
import sys
import csv
from collections import defaultdict
from pathlib import Path

# --- Paramètres ---
MARKERS_TSV = "/home/axel/Desktop/markers/markers_all.tsv"
SAMPLE_BAM  = "sample.bam"
REFERENCE   = "S288C.fsa"
OUTPUT_TSV  = "genotype_results.tsv"

# Ploïdie (1 = haploïde, 2 = diploïde)
PLOIDY = 1

# SNV / INS
MIN_DP_SNV  = 3    # profondeur min pour appeler un SNV/INS
MIN_AF_ALT  = 0.8  # AF min pour appeler ALT (haploïde)
MAX_AF_ALT  = 0.2  # AF max pour appeler REF

# Délétions : couverture moyenne région vs flancs
MIN_DP_DEL       = 2    # couverture min flancs pour avoir un appel
DEL_COVERAGE_MAX = 0.2  # ratio cov_del/cov_flank < 0.2  → délétion homozygote
DEL_HET_MIN      = 0.35 # ratio entre DEL_COVERAGE_MAX et DEL_HET_MIN → hémizygote/het
DEL_HET_MAX      = 0.65
DEL_FLANK_SIZE   = 500  # taille des régions flanquantes

# =============================================================================
# Utilitaires
# =============================================================================

def run(cmd, **kwargs):
    # On capture en bytes pour éviter les UnicodeDecodeError (sorties binaires
    # de samtools/bcftools peuvent contenir des octets non-UTF-8)
    kwargs.pop("text", None)
    kwargs.pop("capture_output", None)
    result = subprocess.run(cmd, shell=isinstance(cmd, str),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            **kwargs)
    result.stdout = result.stdout.decode("utf-8", errors="replace")
    result.stderr = result.stderr.decode("utf-8", errors="replace")
    return result

def ensure_bam_indexed(bam):
    bai = bam + ".bai"
    if not os.path.exists(bai):
        print(f"   Index BAM manquant, création de {bai}...")
        r = run(f"samtools index {bam}")
        if r.returncode != 0:
            sys.exit(f"Erreur indexation BAM :\n{r.stderr}")

def check_dependencies():
    for tool in ("samtools", "bcftools"):
        r = run(f"{tool} --version")
        if r.returncode != 0:
            sys.exit(f"Outil manquant : {tool}")

# =============================================================================
# Chargement des marqueurs
# =============================================================================

def load_markers(tsv):
    markers = []
    with open(tsv) as f:
        next(f)  # skip header
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            markers.append({
                "gene":   fields[0],
                "allele": fields[1],
                "type":   fields[2],   # SNV | INS | DEL
                "chrom":  fields[3],
                "start":  int(fields[4]),
                "end":    int(fields[5]),
                "ref":    fields[6],
                "alt":    fields[7],
                "desc":   fields[8],
            })
    return markers

# =============================================================================
# Génotypage SNV + INS par pileup (bcftools)
# =============================================================================

def genotype_snv_ins(markers, reference, bam, min_dp, min_af, max_af):
    """Génotype SNV et INS via bcftools mpileup."""
    target_markers = [m for m in markers if m["type"] in ("SNV", "INS")]
    if not target_markers:
        return {}

    # BED des positions cibles
    bed_file = "tmp_snv_ins_positions.bed"
    with open(bed_file, "w") as f:
        for m in target_markers:
            # Pour INS, on cible la base précédant l'insertion (convention VCF)
            pos = m["start"] - 1 if m["type"] == "INS" else m["start"] - 1
            f.write(f"{m['chrom']}\t{pos}\t{m['start']}\n")

    # Pileup avec annotation AD et DP
    pileup_vcf = "tmp_snv_ins_pileup.vcf.gz"
    cmd = (f"bcftools mpileup -f {reference} -R {bed_file} "
           f"-q 1 -Q 1 --annotate FORMAT/AD,FORMAT/DP {bam} "
           f"--indel-size 100 "        # autorise des indels jusqu'à 100 bp
           f"-O z -o {pileup_vcf} && bcftools index {pileup_vcf}")
    r = run(cmd)
    if r.returncode != 0:
        print(f"   [WARN] mpileup SNV/INS : {r.stderr[:200]}")

    cmd2 = ["bcftools", "query", "-f", "%CHROM\t%POS\t[%AD]\t[%DP]\n", pileup_vcf]
    result = run(cmd2)

    pileup_data = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        f = line.split("\t")
        chrom, pos = f[0], int(f[1])
        dp     = int(f[3]) if f[3] not in (".", "") else 0
        ad     = f[2].split(",")
        ad_ref = int(ad[0]) if ad[0] not in (".", "") else 0
        ad_alt = int(ad[1]) if len(ad) > 1 and ad[1] not in (".", "") else 0
        pileup_data[(chrom, pos)] = (ad_ref, ad_alt, dp)

    results = {}
    for m in target_markers:
        mid = marker_id(m)
        key = (m["chrom"], m["start"])

        if key not in pileup_data:
            results[mid] = {"call": "NC", "dp": 0, "af": 0.0,
                            "detail": "Non couvert"}
            continue

        ad_ref, ad_alt, dp = pileup_data[key]
        af = ad_alt / dp if dp > 0 else 0.0

        if dp < min_dp:
            call   = "NC"
            detail = f"DP trop faible ({dp}x)"
        elif af >= min_af:
            call   = "MUT" if m["type"] == "SNV" else "INS"
            detail = f"ALT détecté AF={af:.2f} DP={dp}x"
        elif af <= max_af:
            call   = "WT"
            detail = f"REF AF_alt={af:.2f} DP={dp}x"
        else:
            call   = "AMB"
            detail = f"Ambigu AF={af:.2f} DP={dp}x"

        results[mid] = {"call": call, "dp": dp, "af": round(af, 3),
                        "detail": detail}

    # Nettoyage fichiers temporaires
    for tmp in (bed_file, pileup_vcf, pileup_vcf + ".csi"):
        if os.path.exists(tmp):
            os.remove(tmp)

    return results

# =============================================================================
# Génotypage DEL par couverture
# =============================================================================

def get_mean_coverage(bam, chrom, start, end):
    if start >= end:
        return 0.0
    cmd    = ["samtools", "depth", "-a", "-r", f"{chrom}:{start}-{end}", bam]
    result = run(cmd)
    depths = []
    for l in result.stdout.strip().split("\n"):
        if l:
            parts = l.split("\t")
            if len(parts) >= 3:
                depths.append(int(parts[2]))
    return sum(depths) / len(depths) if depths else 0.0

def genotype_deletions(markers, bam, min_dp, del_cov_max, het_min, het_max, flank_size):
    del_markers = [m for m in markers if m["type"] == "DEL"]
    results = {}

    for m in del_markers:
        mid   = marker_id(m)
        chrom = m["chrom"]

        cov_del   = get_mean_coverage(bam, chrom, m["start"], m["end"])
        cov_left  = get_mean_coverage(bam, chrom,
                                      max(1, m["start"] - flank_size),
                                      m["start"] - 1)
        cov_right = get_mean_coverage(bam, chrom,
                                      m["end"] + 1,
                                      m["end"] + flank_size)
        cov_flank = (cov_left + cov_right) / 2 if (cov_left + cov_right) > 0 else 0.0
        ratio     = cov_del / cov_flank if cov_flank > 0 else 0.0

        if cov_flank < min_dp:
            call   = "NC"
            detail = f"Flancs non couverts ({cov_flank:.1f}x)"
        elif ratio <= del_cov_max:
            call   = "DEL"
            detail = (f"Délétion homozygote ratio={ratio:.2f} "
                      f"(del={cov_del:.1f}x flank={cov_flank:.1f}x)")
        elif het_min <= ratio <= het_max:
            # Pertinent surtout en diploïde — signalé aussi en haploïde
            call   = "DEL/WT" if PLOIDY == 2 else "AMB"
            detail = (f"Délétion hémizygote ? ratio={ratio:.2f} "
                      f"(del={cov_del:.1f}x flank={cov_flank:.1f}x)")
        elif ratio >= 0.7:
            call   = "WT"
            detail = (f"Région intacte ratio={ratio:.2f} "
                      f"(del={cov_del:.1f}x flank={cov_flank:.1f}x)")
        else:
            call   = "AMB"
            detail = (f"Ambigu ratio={ratio:.2f} "
                      f"(del={cov_del:.1f}x flank={cov_flank:.1f}x)")

        results[mid] = {"call": call,
                        "dp":   round(cov_del, 1),
                        "af":   round(ratio, 3),
                        "detail": detail}

    return results

# =============================================================================
# Clé unique par marqueur
# =============================================================================

def marker_id(m):
    if m["type"] == "DEL":
        return f"{m['gene']}|{m['allele']}|{m['start']}-{m['end']}"
    return f"{m['gene']}|{m['allele']}|{m['start']}"

def get_result(m, snv_ins_results, del_results):
    mid = marker_id(m)
    if m["type"] == "DEL":
        return del_results.get(mid, {"call": "NC", "dp": 0, "af": 0.0, "detail": ""})
    return snv_ins_results.get(mid, {"call": "NC", "dp": 0, "af": 0.0, "detail": ""})

# =============================================================================
# Export TSV
# =============================================================================

def export_tsv(markers, snv_ins_results, del_results, output_path):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Gene", "Allele", "Type", "Chrom", "Start", "End",
                         "Ref", "Alt", "Call", "DP", "AF_ratio", "Detail", "Description"])
        for m in markers:
            r = get_result(m, snv_ins_results, del_results)
            writer.writerow([
                m["gene"], m["allele"], m["type"],
                m["chrom"], m["start"], m["end"],
                m["ref"], m["alt"],
                r["call"], r["dp"], r["af"], r["detail"], m["desc"]
            ])
    print(f"\n   Résultats exportés → {output_path}")

# =============================================================================
# Affichage terminal
# =============================================================================

SYMBOLS = {
    "MUT":    "✓ MUT  ",
    "INS":    "✓ INS  ",
    "DEL":    "✓ DEL  ",
    "DEL/WT": "~ HET  ",
    "WT":     "  WT   ",
    "AMB":    "? AMB  ",
    "NC":     "  NC   ",
}

def print_results(markers, snv_ins_results, del_results):
    print("\n" + "="*100)
    print(f"{'Gene':<8} {'Allele':<14} {'Type':<4} {'Position':<28} "
          f"{'DP':>6} {'Ratio/AF':>8}  {'Call':<9} Description")
    print("-"*100)

    current_gene = ""
    for m in markers:
        if m["gene"] != current_gene:
            if current_gene:
                print()
            current_gene = m["gene"]

        r   = get_result(m, snv_ins_results, del_results)
        sym = SYMBOLS.get(r["call"], r["call"])

        if m["type"] == "DEL":
            pos = f"{m['chrom']}:{m['start']}-{m['end']}"
        else:
            pos = f"{m['chrom']}:{m['start']}"

        print(f"{m['gene']:<8} {m['allele']:<14} {m['type']:<4} {pos:<28} "
              f"{r['dp']:>6} {r['af']:>8.3f}  [{sym}] {m['desc']}")

    print("="*100)
    legend = ("Légende : ✓ MUT = SNV muté | ✓ INS = insertion | ✓ DEL = délétion homozygote | "
              "~ HET = hémizygote | WT = sauvage | AMB = ambigu | NC = non couvert")
    print(f"\n{legend}")

    # --- Résumé par allèle ---
    print("\n=== Résumé par allèle ===\n")
    allele_calls = defaultdict(lambda: {"mut": 0, "wt": 0, "nc": 0,
                                        "het": 0, "amb": 0, "total": 0})
    for m in markers:
        r   = get_result(m, snv_ins_results, del_results)
        key = f"{m['gene']} — {m['allele']}"
        allele_calls[key]["total"] += 1
        c = r["call"]
        if c in ("MUT", "INS", "DEL"):
            allele_calls[key]["mut"] += 1
        elif c == "DEL/WT":
            allele_calls[key]["het"] += 1
        elif c == "WT":
            allele_calls[key]["wt"] += 1
        elif c == "AMB":
            allele_calls[key]["amb"] += 1
        else:
            allele_calls[key]["nc"] += 1

    detected   = {k: v for k, v in allele_calls.items() if v["mut"] > 0}
    undetected = {k: v for k, v in allele_calls.items() if v["mut"] == 0}

    if detected:
        print("  Allèles détectés :")
        for allele, c in sorted(detected.items()):
            extras = []
            if c["het"] > 0: extras.append(f"{c['het']} hémizygote(s)")
            if c["nc"]  > 0: extras.append(f"{c['nc']} NC")
            extra_str = f"  [{', '.join(extras)}]" if extras else ""
            print(f"    ✓ {allele} : {c['mut']}/{c['total']} marqueurs confirmés{extra_str}")
    else:
        print("  Aucun allèle muté/délété détecté.")

    if undetected:
        print("\n  Allèles non détectés :")
        for allele, c in sorted(undetected.items()):
            status = "NC" if c["nc"] == c["total"] else "WT"
            print(f"    — {allele} : {status} ({c['wt']} WT, {c['nc']} NC, {c['amb']} AMB)")

# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    print("\n=== Génotypage des marqueurs génétiques de levure ===")
    print(f"    Ploïdie configurée : {PLOIDY}n\n")

    print("0. Vérification des dépendances...")
    check_dependencies()

    print("1. Indexation BAM si nécessaire...")
    ensure_bam_indexed(SAMPLE_BAM)

    print("2. Chargement des marqueurs...")
    markers = load_markers(MARKERS_TSV)
    n_snv = sum(1 for m in markers if m["type"] == "SNV")
    n_ins = sum(1 for m in markers if m["type"] == "INS")
    n_del = sum(1 for m in markers if m["type"] == "DEL")
    print(f"   {len(markers)} marqueurs : {n_snv} SNV, {n_ins} INS, {n_del} DEL")

    print("\n3. Génotypage SNV + INS par pileup...")
    snv_ins_results = genotype_snv_ins(
        markers, REFERENCE, SAMPLE_BAM, MIN_DP_SNV, MIN_AF_ALT, MAX_AF_ALT
    )

    print("\n4. Génotypage des délétions par couverture...")
    del_results = genotype_deletions(
        markers, SAMPLE_BAM, MIN_DP_DEL,
        DEL_COVERAGE_MAX, DEL_HET_MIN, DEL_HET_MAX, DEL_FLANK_SIZE
    )

    print_results(markers, snv_ins_results, del_results)

    export_tsv(markers, snv_ins_results, del_results, OUTPUT_TSV)
