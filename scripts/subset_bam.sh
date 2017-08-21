#!/bin/bash
#SBATCH --mem=10000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p vccc
#SBATCH -o subset_bam_%j.out
#SBATCH -e subset_bam_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pdiakumis@unimelb.edu.au
#SBATCH --time=600


set -eu -o pipefail
module load SAMtools/1.3.1-intel-2016.u3

BAM_DIR="../data/FFPE"
CHROM="chr21"

samtools view ${BAM_DIR}/WES005PVFE.bam -b ${CHROM} > ${BAM_DIR}/WES005PVFE_${CHROM}.bam
samtools view ${BAM_DIR}/WES005PVFR.bam -b ${CHROM} > ${BAM_DIR}/WES005PVFR_${CHROM}.bam
samtools view ${BAM_DIR}/WES006JCFE.bam -b ${CHROM} > ${BAM_DIR}/WES006JCFE_${CHROM}.bam
samtools view ${BAM_DIR}/WES006JCFR.bam -b ${CHROM} > ${BAM_DIR}/WES006JCFR_${CHROM}.bam
samtools view ${BAM_DIR}/WES011RCFE.bam -b ${CHROM} > ${BAM_DIR}/WES011RCFE_${CHROM}.bam
samtools view ${BAM_DIR}/WES011RCFR.bam -b ${CHROM} > ${BAM_DIR}/WES011RCFR_${CHROM}.bam
