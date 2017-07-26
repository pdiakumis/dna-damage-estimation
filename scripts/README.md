<!-- vim-markdown-toc GFM -->
* [`split_mapped_reads.pl`](#split_mapped_readspl)
    * [Pileup format](#pileup-format)
    * [Examples](#examples)
* [`estimate_damage.pl`](#estimate_damagepl)

<!-- vim-markdown-toc -->

# `split_mapped_reads.pl`

* Splits a `BAM` file into two files containing the first mate (R1 - flag 64) and the
second mate (R2 - flag 128) of the paired-end reads, using `samtools view -f <flag>`.
* Creates `pileup` files using `samtools mpileup -O -s -q <min_mapq> -Q <min_baseq>`, where:
    - `-O`: output base positions on reads
    - `-s`: output mapping quality

## Pileup format

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

## Examples

* Row1: vanilla
* Row2: with `--no-BAQ`
* Row3: with `-O` (base positions on reads)
* Row4: with `-s` (mapping quality)

```
$> samtools mpileup --region 1:215906528-215906567 --fasta-ref human_g1k_v37.fasta.gz NA12891_CEU_sample.bam
1	215906528	G	21	,,,,,,,,.,,,.,,..,.,,	;=?./:?>>;=7?>>@A?==:
1	215906528	G	21	,,,,,,,,.,,,.,,..,.,,	;=?./:?>>;=7?>>@A?==:
1	215906528	G	21	,,,,,,,,.,,,.,,..,.,,	;=?./:?>>;=7?>>@A?==:	49,42,33,32,44,31,28,40,20,18,28,25,15,15,13,12,11,7,5,21,18
1	215906528	G	21	,,,,,,,,.,,,.,,..,.,,	;=?./:?>>;=7?>>@A?==:	]]]]F]]]F]FF]]]]:]]FF
[...]
1	215906534	G	19	,,,,.,,.,,..,.,,,,^].	=9?<>;;?>=@B>>??13>
1	215906534	G	19	,,,,.,,.,,..,.,,,,^].	=9?<>;;?>=@B>>??13>
1	215906534	G	19	,,,,.,,.,,..,.,,,,^].	=9?<>;;?>=@B>>??13>	48,37,34,46,26,24,34,21,21,19,18,17,13,11,27,24,21,2,1
1	215906534	G	19	,,,,.,,.,,..,.,,,,^].	=9?<>;;?>=@B>>??13>	]]]]F]F]]]]:]]FFF]]
[...]
1	215906539	C	14	,,$.,,,+1g,,....,.	6;244.76>15:6:
1	215906539	C	14	,,$.,,,+1g,,....,.	6;244.76>15:6:
1	215906539	C	14	,,$.,,,+1g,,....,.	6;244.76>15:6:	42,51,31,29,39,36,26,24,23,22,16,8,7,6
1	215906539	C	14	,,$.,,,+1g,,....,.	6;244.76>15:6:	]]F]FF]]]:]F]]
[...]
1	215906547	C	7	,,,....	086=C>=
1	215906547	C	15	gGg$,GggGG,,....	<;80;><9=86=C>=
1	215906547	C	7	,,,....	086=C>=	45,26,15,14,8,4,3
1	215906547	C	7	,,,....	086=C>=	F]]]F:S
[...]
1	215906555	G	12	aaaaAAAAaAAA	>?8?;:335?:A
1	215906555	G	16	.$aaaaaA.AAAaAAA^:A	2@>?8?;<:335?:A>
1	215906555	G	12	aaaaAAAAaAAA	>?8?;:335?:A	48,45,42,23,22,16,12,11,10,8,7,6
1	215906555	G	12	aaaaAAAAaAAA	>?8?;:335?:A	FFF]]F:SF]]:
[...]
```

So in summary, the `pileup` files generated from this script contain:

chr, pos, ref, cov, rbase, baseq, basepos, mapq:

```
[...]
1       215906546       G       6       ..,...  (!;;><  ]]]]:S  33,30,25,13,3,2
1       215906547       C       6       GG,...  !!8=>=  ]]]]:S  34,31,26,14,4,3
1       215906548       G       6       .C,...  !!=746  ]]]]:S  35,32,27,15,5,4
1       215906549       G       6       .$.,... !!>7?;  ]]]]:S  36,33,28,16,6,5
1       215906550       G       5       .,...   !>6?8   ]]]:S   34,29,17,7,6
1       215906551       G       5       .,...   !>7?5   ]]]:S   35,30,18,8,7
1       215906552       G       5       .$,...  !><87   ]]]:S   36,31,19,9,8
1       215906553       G       4       ,...    :<?9    ]]:S    32,20,10,9
1       215906554       C       4       ,...    @483    ]]:S    33,21,11,10
1       215906555       G       4       aAAA    $;31    ]]:S    34,22,12,11
[...]
```

# `estimate_damage.pl`
