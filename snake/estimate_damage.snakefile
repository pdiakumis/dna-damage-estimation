configfile: 'config.yaml'
shell.prefix("set -euo pipefail; ")

rule all:
    input:
        expand('{out_dir}/{sample}_giv_scores.tsv',
               out_dir = config['out_dir'],
               sample = config['samples2'])


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
        fasta = config['fasta2']
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
        counts = config['out_dir'] + '{sample}_f{flag}_counts.tsv'
    log:
        config['log_dir'] + '{sample}_f{flag}_counts.log'
    threads: 1
    shell:
        "perl ../scripts/count_mutations.pl {input.pileup} > {output.counts} 2> {log}"

rule estimate_giv:
    """Estimate GIV score"""
    input:
        counts1 = config['out_dir'] + '{sample}_f64_counts.tsv',
        counts2 = config['out_dir'] + '{sample}_f128_counts.tsv'
    output:
        giv = config['out_dir'] + '{sample}_giv_scores.tsv'
    log:
        config['log_dir'] + '{sample}_giv_scores.log'
    threads: 1
    shell:
        "perl ../scripts/damest.pl --c1 {input.counts1} --c2 {input.counts2} "
        "--id {wildcards.sample} > {output.giv} 2> {log}"
