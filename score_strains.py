#!/usr/bin/env python3
# =============================================================================
# Scoring de concordance pour identification de souche de levure haploïde
# Entrées : merged_snv.vcf.gz + experiment_pileup.vcf.gz
# =============================================================================

import subprocess
import pandas as pd

# --- Paramètres ---
MERGED_SNV        = "merged_snv.vcf.gz"
EXPERIMENT_PILEUP = "experiment_pileup.vcf.gz"
STRAINS           = ["S288C", "W303", "SK1"]  # ordre du header du VCF mergé
MIN_DP            = 5     # profondeur minimale pour considérer une position
MIN_AF_ALT        = 0.8   # seuil AF pour appeler ALT dans l'échantillon
MAX_AF_ALT        = 0.2   # seuil AF pour appeler REF dans l'échantillon

# =============================================================================
# Fonctions
# =============================================================================

def load_strain_genotypes(merged_vcf):
    """
    Charge la matrice de génotypes (0/1) des 3 souches depuis le VCF mergé.
    ./.  → 0 (REF)
    1/.  → 1 (ALT)
    """
    print(f"  Lecture de {merged_vcf}...")
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t[%GT\t]\n", merged_vcf]
    result = subprocess.run(cmd, capture_output=True, text=True)

    genotypes = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        chrom, pos = fields[0], int(fields[1])
        gts = [g.strip() for g in fields[2:] if g.strip()]

        s288c = 1 if gts[0].startswith("1") else 0
        w303  = 1 if gts[1].startswith("1") else 0
        sk1   = 1 if gts[2].startswith("1") else 0

        # Ignorer les positions monomorphes
        if s288c == w303 == sk1:
            continue

        genotypes[(chrom, pos)] = (s288c, w303, sk1)

    # Résumé des catégories de positions
    s288c_w303 = sum(1 for v in genotypes.values() if v[0]==v[1] and v[1]!=v[2])
    s288c_sk1  = sum(1 for v in genotypes.values() if v[0]==v[2] and v[0]!=v[1])
    w303_sk1   = sum(1 for v in genotypes.values() if v[1]==v[2] and v[0]!=v[1])
    print(f"  → {len(genotypes)} positions polymorphes chargées")
    print(f"     S288C=W303 vs SK1   : {s288c_w303}  (non discriminant S288C/W303)")
    print(f"     S288C=SK1  vs W303  : {s288c_sk1}  (discriminant pour W303)")
    print(f"     W303=SK1   vs S288C : {w303_sk1}   (discriminant pour S288C)")
    return genotypes


def load_pileup(pileup_vcf):
    """
    Charge les données AD/DP depuis le pileup de l'échantillon.
    Retourne un dict : (CHROM, POS) -> (ad_ref, ad_alt, dp)
    """
    print(f"  Lecture de {pileup_vcf}...")
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t[%AD]\t[%DP]\n", pileup_vcf]
    result = subprocess.run(cmd, capture_output=True, text=True)

    ad_data = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        chrom, pos = fields[0], int(fields[1])
        dp     = int(fields[3]) if fields[3] not in (".", "") else 0
        ad     = fields[2].split(",")
        ad_ref = int(ad[0]) if ad[0] not in (".", "") else 0
        ad_alt = int(ad[1]) if len(ad) > 1 and ad[1] not in (".", "") else 0
        ad_data[(chrom, pos)] = (ad_ref, ad_alt, dp)

    print(f"  → {len(ad_data)} positions chargées")
    return ad_data
    
def score_concordance(strain_gts, ad_data, strains, min_dp, min_af_alt, max_af_alt):
    """
    Chaque souche est scorée uniquement sur les positions où
    ELLE diffère des DEUX autres — ses positions vraiment privées.
    """
    # Définir les positions discriminantes pour chaque souche
    # S288C : positions où W303==SK1 mais S288C != W303
    # W303  : positions où S288C==SK1 mais W303 != S288C
    # SK1   : positions où S288C==W303 mais SK1  != S288C
    discriminant = {
        "S288C": {k for k, v in strain_gts.items() if v[1]==v[2] and v[0]!=v[1]},
        "W303":  {k for k, v in strain_gts.items() if v[0]==v[2] and v[1]!=v[0]},
        "SK1":   {k for k, v in strain_gts.items() if v[0]==v[1] and v[2]!=v[0]},
    }

    for s in strains:
        print(f"  Positions discriminantes {s} : {len(discriminant[s])}")

    counts = {s: {"concordant": 0, "discordant": 0, "covered": 0, "ambiguous": 0}
              for s in strains}

    for (chrom, pos), (s288c, w303, sk1) in strain_gts.items():
        if (chrom, pos) not in ad_data:
            continue

        ad_ref, ad_alt, dp = ad_data[(chrom, pos)]
        if dp < min_dp:
            continue

        af_alt = ad_alt / dp
        if af_alt >= min_af_alt:
            sample_gt = 1
        elif af_alt <= max_af_alt:
            sample_gt = 0
        else:
            for s in strains:
                if (chrom, pos) in discriminant[s]:
                    counts[s]["ambiguous"] += 1
            continue

        strain_gt_values = {"S288C": s288c, "W303": w303, "SK1": sk1}

        for strain in strains:
            if (chrom, pos) not in discriminant[strain]:
                continue  # position non discriminante pour cette souche
            counts[strain]["covered"] += 1
            if sample_gt == strain_gt_values[strain]:
                counts[strain]["concordant"] += 1
            else:
                counts[strain]["discordant"] += 1

    results = []
    for strain in strains:
        covered    = counts[strain]["covered"]
        concordant = counts[strain]["concordant"]
        discordant = counts[strain]["discordant"]
        ambiguous  = counts[strain]["ambiguous"]
        score      = concordant / covered if covered > 0 else 0
        results.append({
            "Strain":     strain,
            "Covered":    covered,
            "Concordant": concordant,
            "Discordant": discordant,
            "Ambiguous":  ambiguous,
            "Score":      round(score, 4)
        })

    return pd.DataFrame(results).sort_values("Score", ascending=False)

# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    print("\n=== Chargement des génotypes des souches ===")
    strain_gts = load_strain_genotypes(MERGED_SNV)

    print("\n=== Chargement du pileup de l'échantillon ===")
    ad_data = load_pileup(EXPERIMENT_PILEUP)

    print("\n=== Scoring de concordance ===")
    df = score_concordance(
        strain_gts, ad_data, STRAINS,
        MIN_DP, MIN_AF_ALT, MAX_AF_ALT
    )

    print("\n" + "="*60)
    print(df.to_string(index=False))
    print("="*60)

    best = df.iloc[0]
    print(f"\n→ Souche identifiée : {best['Strain']} "
          f"(score {best['Score']:.4f} sur {best['Covered']} positions)")
