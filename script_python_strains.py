import subprocess
import pandas as pd

# --- Paramètres ---
MERGED_SNV        = "merged_snv.vcf.gz"
EXPERIMENT_PILEUP = "experiment_pileup.vcf.gz"
STRAINS           = ["S288C", "W303", "SK1"]
MIN_DP            = 10
MIN_AF_ALT        = 0.8   # seuil pour appeler ALT dans l'échantillon
MAX_AF_ALT        = 0.2   # seuil pour appeler REF dans l'échantillon

# --- Charger la matrice de génotypes des 3 souches ---
def load_strain_genotypes(merged_vcf):
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

        # Garder uniquement les positions polymorphes
        if s288c == w303 == sk1:
            continue

        genotypes[(chrom, pos)] = (s288c, w303, sk1)

    print(f"  {len(genotypes)} positions polymorphes chargées")
    return genotypes

# --- Charger le pileup de l'échantillon ---
def load_pileup(pileup_vcf):
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t[%AD]\t[%DP]\n", pileup_vcf]
    result = subprocess.run(cmd, capture_output=True, text=True)

    ad_data = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        chrom, pos = fields[0], int(fields[1])
        dp  = int(fields[3]) if fields[3] not in (".", "") else 0
        ad  = fields[2].split(",")
        ad_ref = int(ad[0]) if ad[0] not in (".", "") else 0
        ad_alt = int(ad[1]) if len(ad) > 1 and ad[1] not in (".", "") else 0
        ad_data[(chrom, pos)] = (ad_ref, ad_alt, dp)

    print(f"  {len(ad_data)} positions chargées depuis le pileup")
    return ad_data

# --- Scoring de concordance ---
def score_concordance(strain_gts, ad_data, strains, min_dp, min_af_alt, max_af_alt):
    counts = {s: {"concordant": 0, "discordant": 0, "covered": 0} for s in strains}

    for (chrom, pos), (s288c, w303, sk1) in strain_gts.items():
        if (chrom, pos) not in ad_data:
            continue

        ad_ref, ad_alt, dp = ad_data[(chrom, pos)]
        if dp < min_dp:
            continue

        af_alt = ad_alt / dp

        # Appel du génotype de l'échantillon
        if af_alt >= min_af_alt:
            sample_gt = 1
        elif af_alt <= max_af_alt:
            sample_gt = 0
        else:
            continue  # ambigu, on ignore

        strain_gt_values = [s288c, w303, sk1]
        for i, strain in enumerate(strains):
            counts[strain]["covered"] += 1
            if sample_gt == strain_gt_values[i]:
                counts[strain]["concordant"] += 1
            else:
                counts[strain]["discordant"] += 1

    results = []
    for strain in strains:
        covered    = counts[strain]["covered"]
        concordant = counts[strain]["concordant"]
        score      = concordant / covered if covered > 0 else 0
        results.append({
            "Strain":      strain,
            "Covered":     covered,
            "Concordant":  concordant,
            "Discordant":  counts[strain]["discordant"],
            "Score":       round(score, 4)
        })

    return pd.DataFrame(results).sort_values("Score", ascending=False)

# --- Main ---
print("Chargement des génotypes des souches...")
strain_gts = load_strain_genotypes(MERGED_SNV)

print("\nChargement du pileup de l'échantillon...")
ad_data = load_pileup(EXPERIMENT_PILEUP)

print("\n=== Résultats ===\n")
df = score_concordance(strain_gts, ad_data, STRAINS, MIN_DP, MIN_AF_ALT, MAX_AF_ALT)
print(df.to_string(index=False))
