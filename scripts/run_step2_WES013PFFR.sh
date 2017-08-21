#!/bin/bash
#SBATCH --mem=60000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p vccc
#SBATCH -o step2_WES013_%j.out
#SBATCH -e step2_WES013_%j.err
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
printf "[$(date)] Step 2: estimating basic damage\n"
# Step 2: Estimate basic damage
perl 1a-estimate_damage.pl \
    --mpileup1 ${BAM_DIR}/out/chr${CHROM}/R1.mpileup \
    --mpileup2 ${BAM_DIR}/out/chr${CHROM}/R2.mpileup \
    --id WES013_chr${CHROM} \
    > ${BAM_DIR}/out/chr${CHROM}/1-basic.damage

printf "[$(date)] Step 2 Finished!\n"
