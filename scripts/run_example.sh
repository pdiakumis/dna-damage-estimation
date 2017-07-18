#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p vccc
#SBATCH --time=60
#SBATCH -o DamEst_%j.out
#SBATCH -e DamEst_%j.err

set -eu -o pipefail

module load SAMtools/1.3.1-intel-2016.u3
# need older version of Perl or else outputs error
module load Perl/5.16.3-goolf-2015a

DATA_DIR="../data/example"
SAMPLE="NA12891_CEU_sample.bam"
GENOME="../data/genome/human_g1k_v37.fasta.gz"

echo "[$(date)] Starting analysis"

echo "[$(date)] Step 1: splitting into R1 and R2"
# Step 1: split into R1 and R2
perl 0-split_mapped_reads.pl \
  --bam ${DATA_DIR}/${SAMPLE} \
  --genome ${GENOME} \
  --mpileup1 ${DATA_DIR}/out/out1.mpileup \
  --mpileup2 ${DATA_DIR}/out/out2.mpileup \

echo "Wait one minute... just so I can check stuff"
sleep 60

echo "[$(date)] Step 2: estimating basic damage"
# Step 2: Estimate basic damage 
perl 1a-estimate_damage.pl \
  --mpileup1 ${DATA_DIR}/out/out1.mpileup \
  --mpileup2 ${DATA_DIR}/out/out2.mpileup \
  --id foo \
  > ${DATA_DIR}/out/1-basic.damage

echo "[$(date)] Step 3: estimating damage relative to read position"
# Step 3: Estimate damage relative to read position
perl 2a-estimate_damage_location.pl \
  --mpileup1 ${DATA_DIR}/out/out1.mpileup \
  --mpileup2 ${DATA_DIR}/out/out2.mpileup \
  --id foo \
  --out ${DATA_DIR}/out/2-loc.damage

echo "[$(date)] Step 4: estimating damage relative to read position and context"
# Step 4: Estimate damage relative to read position and context
perl 3a-estimate_damage_location_context.pl \
  --mpileup1 ${DATA_DIR}/out/out1.mpileup \
  --mpileup2 ${DATA_DIR}/out/out2.mpileup \
  --id foo \
  --out ${DATA_DIR}/out/3-pos_loc.damage
