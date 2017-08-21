#!/bin/bash
#SBATCH --mem=60000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p vccc
#SBATCH -o step1_WES013_%j.out
#SBATCH -e step1_WES013_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pdiakumis@unimelb.edu.au
#SBATCH --time=600

set -eu -o pipefail

module load SAMtools/1.3.1-intel-2016.u3
# need older version of Perl or else outputs error
module load Perl/5.16.3-goolf-2015a

SAMPLE="WES013"
CHROM="21"
BAM_DIR="../data/${SAMPLE}"
BAM="${BAM_DIR}/WES013PFFR-sort_chr${CHROM}.bam"
GENOME="${BAM_DIR}/GRCh37.fa"

printf "[$(date)] Working on chr${CHROM} of ${SAMPLE}\n"
printf "[$(date)] Starting dna damage estimation analysis for ${BAM}\n"

printf "[$(date)] Step 1: splitting into R1 and R2\n"
# Step 1: split into R1 and R2
perl 0-split_mapped_reads.pl \
  --bam ${BAM} \
  --genome ${GENOME} \
  --mpileup1 ${BAM_DIR}/out/chr${CHROM}/R1.mpileup \
  --mpileup2 ${BAM_DIR}/out/chr${CHROM}/R2.mpileup

printf "[$(date)] Step 1 Finished!\n"
