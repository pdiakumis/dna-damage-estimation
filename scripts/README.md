Description of Scripts
----------------------

Contents
========

<!-- vim-markdown-toc GFM -->
* [`split_mapped_reads.pl`](#split_mapped_readspl)
        * [Pileup format](#pileup-format)
            * [Examples](#examples)

<!-- vim-markdown-toc -->


`split_mapped_reads.pl`
======================
- Splits a `BAM` file into two files containing the first mate (R1 - flag 64) and the
second mate (R2 - flag 128) of the paired-end reads, using `samtools view -f <flag>`.
- Creates `pileup` files using `samtools mpileup -O -s -q <min_mapq> -Q <min_baseq>`, where:
    `-O`: output base positions on reads
    `-s`: output mapping quality

### Pileup format
In the pileup format (without -u or -g), each line represents
a genomic position, consisting of:

| column   | description                |
|----------|----------------------------|
| chr      | chromosome name            |
| pos      | 1-based coordinate         |
| ref      | reference base             |
| cov      | depth of coverage at site  |
| rbase    | read bases                 |
| baseq    | base qualities             |
| mapq     | mapping qualities          |

Information on match/mismatch, indel, strand, mapping quality and start/end of
a read are all encoded at the __read base__ column. In this column:

- dot: match ref on forward strand
- comma:  match ref on reverse strand
- '>' or '<': reference skip,
- [ACGTN]: mismatch on the forward strand
- [acgtn]: mismatch on the reverse strand.
- `\\+[0-9]+[ACGTNacgtn]+` pattern: insertion between this and next ref position
  Integer = length of insertion, followed by inserted sequence.
- `-[0-9]+[ACGTNacgtn]+` pattern : deletion from the reference.
  Deleted bases will be presented as `*` in the following lines.
- `^`: start of a read. The ASCII of the character following `^` minus 33 gives
  the mapping quality.
- `$': end of a read.

#### Examples

```
$> samtools mpileup --no-BAQ --region 1:215906528-215906567 --fasta-ref human_g1k_v37.fasta.gz NA12891_CEU_sample.bam
1	215906528	G	21	,,,,,,,,.,,,.,,..,.,,	;=?./:?>>;=7?>>@A?==:
[...]
1	215906534	G	19	,,,,.,,.,,..,.,,,,^].	=9?<>;;?>=@B>>??13>
[...]
1	215906539	C	14	,,$.,,,+1g,,....,.	6;244.76>15:6:
[...]
1	215906547	C	15	gGg$,GggGG,,....	<;80;><9=86=C>=
[...]
1	215906555	G	16	.$aaaaaA.AAAaAAA^:A	2@>?8?;<:335?:A>
[...]
```

