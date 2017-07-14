# DNA Damage Estimation
Peter Diakumis  
14 July 2017  



Here we'll run the DNA Damage Estimator scripts from the
[Ettwiller GitHub repo](https://github.com/Ettwiller/Damage-estimator).

# Example Run
An example BAM file was downloaded from
[bds](https://github.com/vsbuffalo/bds-files/tree/master/chapter-11-alignment):


```bash
ls -lLh ../data/example
```

```
## total 122680
## -rw-r--r--  1 diakumis  10908    60M 14 Jul 11:10 NA12891_CEU_sample.bam
## -rw-r--r--  1 diakumis  10908   106K 14 Jul 11:26 NA12891_CEU_sample.bam.bai
```

## Step 1: Split BAM into R1 and R2


```bash
perl ../scripts/split_mapped_reads.pl \
  --bam ../data/example/NA12891_CEU_sample.bam \
  --genome ../data/genome/human_g1k_v37.fasta \
  --mpileup1 ../data/out/out1.mpileup \
  --mpileup2 ../data/out/out2.mpileup
```

```
## [mpileup] 1 samples in 1 input files
## <mpileup> Set max per-file depth to 8000
## [mpileup] 1 samples in 1 input files
## <mpileup> Set max per-file depth to 8000
```

