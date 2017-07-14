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
total 122680
-rw-r--r--  1 diakumis  10908    60M 14 Jul 11:10 NA12891_CEU_sample.bam
-rw-r--r--  1 diakumis  10908   106K 14 Jul 11:26 NA12891_CEU_sample.bam.bai
```

## Step 1: Split BAM into R1 and R2


```bash
perl ../scripts/split_mapped_reads.pl \
  --bam ../data/example/NA12891_CEU_sample.bam \
  --genome ../data/genome/human_g1k_v37.fasta \
  --mpileup1 ../data/out/out1.mpileup \
  --mpileup2 ../data/out/out2.mpileup
```

Output:


```bash
ls -lLh ../data/out
```

```
total 221424
-rw-r--r--  1 diakumis  10908    54M 14 Jul 17:01 out1.mpileup
-rw-r--r--  1 diakumis  10908    54M 14 Jul 17:02 out2.mpileup
```


```bash
head -n5 ../data/out/*.mpileup
```

```
==> ../data/out/out1.mpileup <==
1	215622850	G	1	^].	;	]	1
1	215622851	G	1	.	=	]	2
1	215622852	A	1	.	<	]	3
1	215622853	A	1	.	>	]	4
1	215622854	T	1	.	=	]	5

==> ../data/out/out2.mpileup <==
1	215622860	T	1	^].	;	]	1
1	215622861	A	1	.	9	]	2
1	215622862	G	1	.	=	]	3
1	215622863	G	2	.^F,	;=	]F	4,17
1	215622864	A	2	.,	<<	]F	5,18
```

