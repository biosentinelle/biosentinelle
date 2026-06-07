process GenotypeMarkers {
    label 'python'
    cache 'true'

    input:
    tuple val(ref_name), path(ref), val(bam_name), path(bam)
    path markers_tsv

    output:
    path "genotype_results.tsv"
    path "genotype_markers_report.txt", emit: genotype_markers_log_ch

    script:
    """
    set -euo pipefail

    echo -e "\\n\\n<br /><br />\\n\\n### Genotype markers for ${ref_name}\\n\\n" > genotype_markers_report.txt
    echo "using existing BAM ${bam}" | tee -a genotype_markers_report.txt

    python3 - "${markers_tsv}" "${bam}" "${ref}" "genotype_results.tsv" <<'PY'
    import os
    import sys
    import subprocess
    from pathlib import Path

    markers_tsv, sample_bam, reference, output_tsv = sys.argv[1:5]

    MIN_DP_SNV = 3
    MIN_AF_ALT = 0.8
    MAX_AF_ALT = 0.2
    MIN_DP_DEL = 2
    DEL_COVERAGE_MAX = 0.2
    DEL_HET_MIN = 0.35
    DEL_HET_MAX = 0.65
    DEL_FLANK_SIZE = 500

    def run(cmd):
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise SystemExit(result.stderr)
        return result

    def load_markers(tsv):
        markers = []
        with open(tsv) as handle:
            next(handle, None)
            for line in handle:
                if not line.strip():
                    continue
                fields = line.rstrip("\\n").split("\\t")
                markers.append({
                    "gene": fields[0],
                    "allele": fields[1],
                    "type": fields[2],
                    "chrom": fields[3],
                    "start": int(fields[4]),
                    "end": int(fields[5]),
                    "ref": fields[6],
                    "alt": fields[7],
                    "desc": fields[8],
                })
        return markers

    def marker_id(marker):
        return f"{marker['gene']}:{marker['allele']}:{marker['type']}"

    def ensure_bam_indexed(bam):
        bai = bam + ".bai"
        if not os.path.exists(bai):
            run(["samtools", "index", bam])

    def get_mean_coverage(bam, chrom, start, end):
        if start >= end:
            return 0.0
        result = run(["samtools", "depth", "-a", "-r", f"{chrom}:{start}-{end}", bam])
        depths = []
        for line in result.stdout.strip().split("\\n"):
            if not line:
                continue
            parts = line.split("\\t")
            if len(parts) >= 3:
                depths.append(int(parts[2]))
        return sum(depths) / len(depths) if depths else 0.0

    def genotype_snv_ins(markers, reference, bam):
        target_markers = [m for m in markers if m["type"] in ("SNV", "INS")]
        results = {}
        if not target_markers:
            return results

        bed_file = "tmp_snv_ins_positions.bed"
        with open(bed_file, "w") as handle:
            for marker in target_markers:
                pos = marker["start"] - 1
                handle.write(f"{marker['chrom']}\\t{pos}\\t{marker['start']}\\n")

        pileup_vcf = "tmp_snv_ins_pileup.vcf.gz"
        run([
            "bcftools", "mpileup", "-f", reference, "-R", bed_file,
            "-q", "1", "-Q", "1", "--annotate", "FORMAT/AD,FORMAT/DP",
            bam, "--indel-size", "100", "-O", "z", "-o", pileup_vcf
        ])
        run(["bcftools", "index", pileup_vcf])

        result = run(["bcftools", "query", "-f", "%CHROM\\t%POS\\t[%AD]\\t[%DP]\\n", pileup_vcf])
        pileup_data = {}
        for line in result.stdout.strip().split("\\n"):
            if not line:
                continue
            chrom, pos, ad_field, dp_field = line.split("\\t")
            dp = int(dp_field) if dp_field not in (".", "") else 0
            ad = ad_field.split(",")
            ad_ref = int(ad[0]) if ad and ad[0] not in (".", "") else 0
            ad_alt = int(ad[1]) if len(ad) > 1 and ad[1] not in (".", "") else 0
            pileup_data[(chrom, int(pos))] = (ad_ref, ad_alt, dp)

        for marker in target_markers:
            mid = marker_id(marker)
            key = (marker["chrom"], marker["start"])
            if key not in pileup_data:
                results[mid] = {"call": "NC", "dp": 0, "af": 0.0, "detail": "Non couvert"}
                continue
            ad_ref, ad_alt, dp = pileup_data[key]
            af = ad_alt / dp if dp > 0 else 0.0
            if dp < MIN_DP_SNV:
                call = "NC"
                detail = f"DP trop faible ({dp}x)"
            elif af >= MIN_AF_ALT:
                call = "MUT" if marker["type"] == "SNV" else "INS"
                detail = f"ALT detecte AF={af:.2f} DP={dp}x"
            elif af <= MAX_AF_ALT:
                call = "WT"
                detail = f"REF AF_alt={af:.2f} DP={dp}x"
            else:
                call = "AMB"
                detail = f"Ambigu AF={af:.2f} DP={dp}x"
            results[mid] = {"call": call, "dp": dp, "af": round(af, 3), "detail": detail}

        for tmp in (bed_file, pileup_vcf, pileup_vcf + ".csi"):
            if os.path.exists(tmp):
                os.remove(tmp)
        return results

    def genotype_deletions(markers, bam):
        del_markers = [m for m in markers if m["type"] == "DEL"]
        results = {}
        for marker in del_markers:
            mid = marker_id(marker)
            cov_del = get_mean_coverage(bam, marker["chrom"], marker["start"], marker["end"])
            cov_left = get_mean_coverage(bam, marker["chrom"], max(1, marker["start"] - DEL_FLANK_SIZE), marker["start"] - 1)
            cov_right = get_mean_coverage(bam, marker["chrom"], marker["end"] + 1, marker["end"] + DEL_FLANK_SIZE)
            cov_flank = (cov_left + cov_right) / 2 if (cov_left + cov_right) > 0 else 0.0
            ratio = cov_del / cov_flank if cov_flank > 0 else 0.0

            if cov_flank < MIN_DP_DEL:
                call = "NC"
                detail = f"Flancs non couverts ({cov_flank:.1f}x)"
            elif ratio <= DEL_COVERAGE_MAX:
                call = "DEL"
                detail = f"Deletion homozygote ratio={ratio:.2f}"
            elif DEL_HET_MIN <= ratio <= DEL_HET_MAX:
                call = "AMB"
                detail = f"Deletion hemizygote ? ratio={ratio:.2f}"
            elif ratio >= 0.7:
                call = "WT"
                detail = f"Region intacte ratio={ratio:.2f}"
            else:
                call = "AMB"
                detail = f"Ratio intermediaire={ratio:.2f}"

            results[mid] = {
                "call": call,
                "dp": round(cov_flank, 2),
                "af": round(ratio, 3),
                "detail": detail,
            }
        return results

    def export_tsv(markers, snv_ins_results, del_results, output_path):
        with open(output_path, "w") as handle:
            handle.write("gene\\tallele\\ttype\\tchrom\\tstart\\tend\\tref\\talt\\tcall\\tdp\\taf\\tdescription\\n")
            for marker in markers:
                mid = marker_id(marker)
                result = snv_ins_results.get(mid) or del_results.get(mid) or {"call": "NC", "dp": 0, "af": 0.0, "detail": "Non traite"}
                handle.write(
                    f"{marker['gene']}\\t{marker['allele']}\\t{marker['type']}\\t{marker['chrom']}\\t"
                    f"{marker['start']}\\t{marker['end']}\\t{marker['ref']}\\t{marker['alt']}\\t"
                    f"{result['call']}\\t{result['dp']}\\t{result['af']}\\t{marker['desc']}\\n"
                )

    if not markers_tsv or not os.path.exists(markers_tsv):
        raise SystemExit("markers_tsv manquant")

    markers = load_markers(markers_tsv)
    ensure_bam_indexed(sample_bam)
    snv_ins_results = genotype_snv_ins(markers, reference, sample_bam)
    del_results = genotype_deletions(markers, sample_bam)
    export_tsv(markers, snv_ins_results, del_results, output_tsv)
    PY
    """
}
