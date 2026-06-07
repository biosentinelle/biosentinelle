# TODO - Genotype MATA/MATALPHA Pipeline

## Task: Create a pipeline to check if the genotype is MATA or MATALPHA using corresponding reference files

### Steps:
1. [x] Create directory `mating_type_refs/` containing MATA.fasta and MATALPHA.fasta
2. [x] Modify nextflow.config to set `ref_path` to the new directory
3. [x] Add genotype determination process that outputs "MATA" or "MATALPHA" based on alignment results
4. [x] Integrate the genotype output into the main workflow

### Implementation Details:
- Use existing Bowtie2_ref pipeline to align against both references
- Use SNP counting to determine genotype (lowest SNP count = correct genotype)
- Add specific output reporting "MATA" or "MATALPHA" as the genotype
