#!/bin/bash
# =============================================================================
# Pipeline d'identification de souche de levure haploïde
# Souches : S288C, W303, SK1 sur génome de référence SGDref
# =============================================================================

set -euo pipefail

# --- Paramètres à modifier ---
REFERENCE="SGDref.asm01.HP0.nuclear_genome.tidy.fa"  # based on S288C 
FASTQ="/home/axel/Desktop/alignment/alignment5/SRR38694169.bis_1.fastq"

bowtie2-build $REFERENCE ref_gen
echo "aligning "$FASTQ" on "$REFERENCE
bowtie2 --local --very-sensitive-local -p 10  $FASTQ  -x ref_gen -S sample.sam
samtools sort sample.sam  -o sample.bam


SAMPLE_BAM="sample.bam"

# from G. Liti work 
VCF_S288C="pav_S288C.asm01.HP0.vcf.gz"
VCF_W303="pav_W303.asm01.HP0.vcf.gz"
VCF_SK1="pav_SK1.asm01.HP0.vcf.gz"

# --- Couleurs pour les logs ---
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

log() { echo -e "${BLUE}[$(date '+%H:%M:%S')]${NC} $1"; }
ok()  { echo -e "${GREEN}[OK]${NC} $1"; }

# =============================================================================
# ETAPE 1 — Extraire les SNVs de chaque souche depuis les VCF PAV
# =============================================================================
log "Etape 1 : Extraction des SNVs depuis les VCF PAV..."

for strain in S288C W303 SK1; do
    eval VCF=\$VCF_${strain}
    bcftools view -i 'INFO/SVTYPE="SNV"' ${VCF} -O z -o ${strain}_snv.vcf.gz
    bcftools index ${strain}_snv.vcf.gz
    COUNT=$(bcftools stats ${strain}_snv.vcf.gz | grep "number of records" | tail -1 | awk '{print $NF}')
    ok "${strain} : ${COUNT} SNVs"
done

# =============================================================================
# ETAPE 2 — Merger les 3 VCF
# =============================================================================
log "Etape 2 : Merge des 3 VCF..."

bcftools merge S288C_snv.vcf.gz W303_snv.vcf.gz SK1_snv.vcf.gz \
    -O z -o merged_snv.vcf.gz
bcftools index merged_snv.vcf.gz

COUNT=$(bcftools stats merged_snv.vcf.gz | grep "number of records" | tail -1 | awk '{print $NF}')
ok "Merged : ${COUNT} SNVs totaux"

# Vérifier l'ordre des samples
log "Ordre des samples dans le VCF mergé :"
bcftools view -h merged_snv.vcf.gz | tail -1

# =============================================================================
# ETAPE 3 — Pileup de l'échantillon sur toutes les positions diagnostiques
# =============================================================================
log "Etape 3 : Pileup de l'échantillon sur les positions polymorphes..."

bcftools mpileup \
    -f ${REFERENCE} \
    -T merged_snv.vcf.gz \
    -q 20 -Q 20 \
    --annotate FORMAT/AD,FORMAT/DP \
    ${SAMPLE_BAM} \
    -O z -o experiment_pileup.vcf.gz

bcftools index experiment_pileup.vcf.gz

COUNT=$(bcftools query -f '%CHROM\t%POS\n' experiment_pileup.vcf.gz | wc -l)
ok "Pileup : ${COUNT} positions couvertes"

# =============================================================================
# ETAPE 4 — Scoring Python
# =============================================================================
log "Etape 4 : Scoring de concordance..."

python3 score_strains.py

ok "Pipeline terminé."
