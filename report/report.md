# DNA Damage Estimation
Peter Diakumis  
14 July 2017  



```r
library(ggplot2)
library(ggforce)
library(dplyr)
library(readr)
```

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

### Command Line


```bash
perl ../scripts/0-split_mapped_reads.pl \
  --bam ../data/example/NA12891_CEU_sample.bam \
  --genome ../data/genome/human_g1k_v37.fasta \
  --mpileup1 ../data/out/out1.mpileup \
  --mpileup2 ../data/out/out2.mpileup
```

### Output


```bash
ls -lLh ../data/out
```

```
total 221480
-rw-r--r--  1 diakumis  10908   1.2K 14 Jul 17:45 foo.damage
-rw-r--r--  1 diakumis  10908    23K 17 Jul 10:27 foo2.damage
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

The script uses `samtools mpileup` to output a summary of the read pileup in the
given BAM file. Options used are:

* `-O`: output base positions on reads 
* `-s`: output mapping quality
* `-q`: skip alignments with mapQ smaller than [10]
* `-Q`: skip bases with baseQ/BAQ smaller than [0]


## Step 2: Estimate basic damage 

### Command Line

```bash
perl ../scripts/1a-estimate_damage.pl \
  --mpileup1 ../data/out/out1.mpileup \
  --mpileup2 ../data/out/out2.mpileup \
  --id foo \
  > ../data/out/foo.damage
```

### Output


```bash
head -n5 ../data/out/foo.damage
```

```
33	T_C	foo	5.17234082904787e-05	A_G-T_C	0.380401132513564
29	T_G	foo	4.54539042552691e-05	A_C-T_G	2.17289737814566
12	A_C	foo	2.09185692395926e-05	A_C-T_G	0.460215015240801
9	C_+	foo	0.000178126113288208	C_+-G_+	3.20967422713059
8	G_C	foo	4.9330340626002e-05	C_G-G_C	0.113293854112244
```

Column description:

1. raw count of variant type
2. variant type (ex. G_T, G to T)
3. id (from the --id option)
4. frequency of variant
5. family (the variant type and reverse complement)
6. GIV-score

If you have followed the standard protocol for acoustic shearing during library preparation you should obtain a GIV score for G_T around 2.

### Plot


```r
type_clean <-c("G_T", "C_A", "C_T", "G_A", "T_A", "A_T",
               "A_G", "T_C", "C_G", "G_C", "T_G", "A_C")
mut <- readr::read_tsv("../data/out/foo.damage",
                       col_names =  c("abs", "type", "experiment", "count", "family", "damage"),
                       col_types = "iccdcd") %>% 
  filter(type %in% type_clean) %>% 
  mutate(type = factor(type, level = type_clean))

#coloring scheme (feel free to change)
local_color <- c("cornflowerblue", "royalblue4", paste0("grey", c(1, seq(10, 100, 10))))

g <- ggplot(mut, aes(x = reorder(type, damage), y = log2(damage), color = experiment))

g + geom_point(alpha = 0.6, size=1.5) +
  scale_colour_manual(values = local_color) +
  geom_hline(yintercept = log2(1.5), color = "#990000", linetype = "dashed") +
  annotate("text", x = 4, y = log2(1.6), color = "#990000",
           label = "Above this line 1/3 of variants is due to damage") +
  geom_hline(yintercept = 0, color = "grey") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
        legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1, size=11)) +
  ggtitle("GIV scores for variant types")
```

![](report_files/figure-html/example_plot1-1.png)<!-- -->


## Step 3: Estimate damage relative to read position

### Command Line

```bash
perl ../scripts/2a-estimate_damage_location.pl \
  --mpileup1 ../data/out/out1.mpileup \
  --mpileup2 ../data/out/out2.mpileup \
  --id foo2 \
  --out ../data/out/foo2.damage
```

### Output


```bash
head -n5 ../data/out/foo2.damage
```

```
foo2	G_+	R1	0.00142247510668563	2	36
foo2	G_+	R2	0.0103092783505155	1	36
foo2	A_G	R1	0.000260710871643348	3	4
foo2	A_G	R2	0.000272331154684096	1	4
foo2	A_G	R1	0.000270929287455974	3	6
```

Column description:

1. id (from the --id option)
2. variant type (ex. G_T, G to T)
3. R1 or R2
4. count (freq)
5. absolute counts
6. position on the read

### Plot


```r
mut <- readr::read_tsv("../data/out/foo2.damage",
                       col_names = c("experiment", "type", "read", "count", "abs", "loc"),
                       col_types = c("cccdii"))

ggplot(mut) +
  geom_point(aes(x = loc, y = count)) +
  theme_bw() +
  facet_grid(type~read, scales = "fixed")
```

![](report_files/figure-html/example_plot2-1.png)<!-- -->

