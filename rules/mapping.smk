# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: mapping.smk
# time: 2021/05/24
import os


def find_trimed_fastq(wildcards):
    sample, condition = find_library(wildcards)
    num = len(config['samples'][sample][condition][wildcards.library])
    fastq = [
        config['workspace'] + f'/preprocess/{wildcards.library}/{wildcards.library}_r{i}.fq.gz'
        for i in range(1, num + 1)
    ]
    for file in fastq:
        if not os.path.isfile(file):
            raise ValueError(
                f'file {file} not found, you may want to use `--forcerun trim` to regenerate!'
            )
    else:
        return fastq


rule mapping:
    output:
        config['workspace'] + '/mapping/{library}/{library}.sorted.bam'
    input:
        rules.trim.output.fq1
    priority:
        10
    log:
       bwa=config['workspace'] + '/log/mapping/{library}_bwa.log',
       samblaster=config['workspace'] + '/log/mapping/{library}_samblaster.log'
    params:
        rg='\'@RG\\tID:{library}\\tSM:{library}\\tLB:{library}\\tPL:ILLUMINA\'',
        index=config['genome']['bwa_index'],
        tmp=config['workspace'] + '/mapping/{library}/{library}.tmp',
    threads:
        8 if workflow.cores > 8 else workflow.cores
    run:
        fastq = find_trimed_fastq(wildcards)
        bwa = f'bwa mem -t {threads} -R {params.rg} {params.index} {" ".join(fastq)} 2>{log.bwa}'
        if is_pe(wildcards):
            samblaster = f'samblaster 2>{log.samblaster}'
        else:
            samblaster = f'samblaster --ignoreUnmated 2>{log.samblaster}'
        samtools = f'samtools sort -T {params.tmp} -o {output} -'
        shell(' | '.join((bwa, samblaster, samtools)))


rule index:
    output:
        config['workspace'] + '/mapping/{library}/{library}.sorted.bam.bai'
    input:
        rules.mapping.output
    shell:
        'samtools index {input}'


rule chrom_sizes:
    output:
        config['workspace'] + '/mapping/{library}/{library}.chrom_sizes'
    input:
        bam=rules.mapping.output,
        index=rules.index.output
    shell:
        'samtools idxstats {input.bam}'
        ' | grep \'^chr[0-9XY]\{{1,2\}}[[:space:]]\''
        ' | awk \'BEGIN{{OFS="\\t"}}{{print $1, $2}}\''
        ' | sort -k1,1V > {output}'


rule clean_bed:
    output:
        config['workspace'] + '/mapping/{library}/{library}_clean.bed'
    input:
        rules.chrom_sizes.output
    params:
        blacklist=config['genome']['blacklist']
    shell:
        'bedtools complement -g {input} -i {params.blacklist} > {output}'


rule fragment:
    output:
        temp(config['workspace'] + '/mapping/{library}/{library}_tag.bed')
    input:
        bam=rules.mapping.output,
        index=rules.index.output,
        bed=rules.clean_bed.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/bam2tag.py'
    run:
        if is_pe(wildcards):
            shell(f'python {BASE_DIR}/tools/bam2tag.py -p -L {input.bed} {input.bam} {output}')
        else:
            shell(f'python {BASE_DIR}/tools/bam2tag.py -L {input.bed} {input.bam} {output}')


rule sort_fragment:
    output:
        config['workspace'] + '/mapping/{library}/{library}_tag.sorted.bed'
    input:
        rules.fragment.output
    shell:
        'sort -k1,1V -k2,2n -k3,3n {input} > {output}'
    
    
