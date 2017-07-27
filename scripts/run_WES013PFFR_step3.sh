#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p vccc
#SBATCH -o WES013_step3_%j.out
#SBATCH -e WES013_step3_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pdiakumis@unimelb.edu.au
#SBATCH --time=1200

set -eu -o pipefail

module load SAMtools/1.3.1-intel-2016.u3
# need older version of Perl or else outputs error
module load Perl/5.16.3-goolf-2015a

SAMPLE="WES013"
BAM_DIR="../data/WES013"
BAM="${BAM_DIR}/WES013PFFR-sort_chr20.bam"
GENOME="${BAM_DIR}/GRCh37.fa"

printf "[$(date)] Starting analysis\n"
printf "[$(date)] Step 3: estimating damage relative to read position\n"

# Step 3: Estimate damage relative to read position
perl 2a-estimate_damage_location.pl \
    --mpileup1 ${BAM_DIR}/out/${SAMPLE}_R1.mpileup \
    --mpileup2 ${BAM_DIR}/out/${SAMPLE}_R2.mpileup \
    --id WES013 \
    --out ${BAM_DIR}/out/2-loc.damage

printf "[$(date)] Step 3 Finished!\n"
