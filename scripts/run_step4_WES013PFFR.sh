#!/bin/bash
#SBATCH --mem=20000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p vccc
#SBATCH -o step4_%j.out
#SBATCH -e step4_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pdiakumis@unimelb.edu.au
#SBATCH --time=600

set -eu -o pipefail

module load SAMtools/1.3.1-intel-2016.u3
# need older version of Perl or else outputs error
module load Perl/5.16.3-goolf-2015a

BAM_DIR="../data/MH17B001P004"
CHROM="21"
SAMPLE="FOO"
BAM="${BAM_DIR}/${SAMPLE}_chr${CHROM}.bam"

printf "[$(date)] Step 4: estimating damage relative to read position and context\n"
# Step 4: Estimate damage relative to read position and context
perl 3a-estimate_damage_location_context.pl \
    --mpileup1 ${BAM_DIR}/out/chr${CHROM}/${SAMPLE}_R1.mpileup \
    --mpileup2 ${BAM_DIR}/out/chr${CHROM}/${SAMPLE}_R2.mpileup \
    --id ${SAMPLE}_chr${CHROM} \
    --out ${BAM_DIR}/out/chr${CHROM}/3-${SAMPLE}_loc_cont.damage

printf "[$(date)] Step 4 Finished!\n"
