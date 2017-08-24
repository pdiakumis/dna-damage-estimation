configfile: 'config.yaml'
shell.prefix("set -euo pipefail; ")

rule all:
    input:
        expand('{out_dir}/{sample}_tot_damage.tsv',
               out_dir = config['out_dir'],
               sample = config['samples']),
        expand('{out_dir}/{sample}_pos_damage.tsv',
               out_dir = config['out_dir'],
               sample = config['samples']),
        expand('{out_dir}/{sample}_con_damage.tsv',
               out_dir = config['out_dir'],
               sample = config['samples'])



rule samtools_filter:
    """Filter BAM based on flag"""
    input:
        bam = config['bam_dir'] + '{sample}.bam'
    output:
        bam = config['out_dir'] + '{sample}_f{flag}.bam'
    log:
        config['log_dir'] + '{sample}_f{flag}.bam.log'
    threads: 1
    shell:
        "samtools view -b -f {wildcards.flag} {input.bam} > {output.bam} 2> {log}"


rule samtools_mpileup:
    """Generate read pileup"""
    input:
        bam = config['out_dir'] + '{sample}_f{flag}.bam',
        fasta = config['fasta']
    output:
        pileup = config['out_dir'] + '{sample}_f{flag}.mpileup'
    params:
        mapq = config['map_qual'],
        baseq = config['base_qual']
    log:
        config['log_dir'] + '{sample}_f{flag}.mpileup.log'
    threads: 1
    shell:
        "samtools mpileup -O -s -q {params.mapq} -Q {params.baseq} "
        "-f {input.fasta} {input.bam} > {output.pileup} 2> {log}"

rule count_mutations:
    """Count mutations in read pileup"""
    input:
        pileup = config['out_dir'] + '{sample}_f{flag}.mpileup'
    output:
        counts_tot = config['out_dir'] + '{sample}_f{flag}_counts_tot.tsv',
        counts_pos = config['out_dir'] + '{sample}_f{flag}_counts_pos.tsv',
        counts_con = config['out_dir'] + '{sample}_f{flag}_counts_con.tsv'
    log:
        config['log_dir'] + '{sample}_f{flag}_counts.log'
    threads: 1
    shell:
        "perl ../scripts/count_mut.pl --in_mp {input.pileup} "
        "--out_ct {output.counts_tot} --out_cp {output.counts_pos} "
        "--out_cc {output.counts_con} 2> {log}"

rule estimate_damage:
    """Estimate damage scores"""
    input:
        ct1 = config['out_dir'] + '{sample}_f64_counts_tot.tsv',
        ct2 = config['out_dir'] + '{sample}_f128_counts_tot.tsv',
        cp1 = config['out_dir'] + '{sample}_f64_counts_pos.tsv',
        cp2 = config['out_dir'] + '{sample}_f128_counts_pos.tsv',
        cc1 = config['out_dir'] + '{sample}_f64_counts_con.tsv',
        cc2 = config['out_dir'] + '{sample}_f128_counts_con.tsv'
    output:
        tot = config['out_dir'] + '{sample}_tot_damage.tsv',
        pos = config['out_dir'] + '{sample}_pos_damage.tsv',
        con = config['out_dir'] + '{sample}_con_damage.tsv'
    log:
        config['log_dir'] + '{sample}_damest.log'
    threads: 1
    shell:
        "perl ../scripts/damest.pl "
        "--ct1 {input.ct1} --ct2 {input.ct2} "
        "--cp1 {input.cp1} --cp2 {input.cp2} "
        "--cc1 {input.cc1} --cc2 {input.cc2} "
        "--out_tot {output.tot} --out_pos {output.pos} "
        "--out_con {output.con} --id {wildcards.sample} 2> {log}"
