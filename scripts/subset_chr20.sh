#!/bin/bash
#SBATCH --mem=8000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p vccc
#SBATCH -o WES013_%j.out
#SBATCH -e WES013_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pdiakumis@unimelb.edu.au
#SBATCH --time=360

set -eu -o pipefail
module load SAMtools/1.3.1-intel-2016.u3

BAM_DIR="../data/WES013"
BAM="${BAM_DIR}/WES013PFFR-sort.bam"
CHROM=20

[ ! -f $BAM ] && printf "BAM file not found!\n" && exit 0

printf "[$(date)] Subsetting chromosome ${CHROM} from ${BAM}\n"
samtools view ${BAM} -b ${CHROM} > "${BAM_DIR}/$(basename $BAM .bam)_chr${CHROM}.bam"
printf "[$(date)] Subsetting done successfully!"
