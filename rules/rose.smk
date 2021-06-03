# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: rose.smk
# time: 2021/05/31
import os


rule bed2gff:
    output:
        config['workspace'] + '/callpeak/{sample}/{sample}_peaks.gff',
    input:
        rules.callpeak.output.narrowPeak
    shell:
        'awk \'BEGIN {{OFS="\\t"}}{{print $1,$4,".",$2,$3,".",".",".",$4}}\' {input} > {output}'


def find_treat_bams(wildcards):
    return [
       config['workspace'] + f'/mapping/{library}/{library}.sorted.bam'
       for library in config['samples'][wildcards.sample]['treat']
    ]


def find_control_bams(wildcards):
    if 'input' in config['samples'][wildcards.sample]:
        return [
            config['workspace'] + f'/mapping/{library}/{library}.sorted.bam'
            for library in config['samples'][wildcards.sample]['input']
        ]
    return list()


rule rose:
    output:
        config['workspace'] + '/ROSE/{sample}/{sample}_peaks_AllEnhancers.table.txt'
    input:
        gff=rules.bed2gff.output,
        treat=find_treat_bams,
        control=find_control_bams
    log:
        config['workspace'] + '/log/rose/{sample}_rose.log'
    run:
        outdir = os.path.dirname(str(output))
        rose = (
            f'python {BASE_DIR}/tools/ROSE/ROSE_main.py'
            f' -g {config["genome"]["name"].upper()} -s 12500 -t {config["tss"]}'
            f' -r {",".join(input.treat)} -i {input.gff} -o {outdir}'
        )
        if input.control:
            rose += f' -c {",".join(input.control)}'
        rose += f' >{log} 2>&1'
        shell(rose)


rule gene_mapper:
    output:
        gene=config['workspace'] + '/ROSE/{sample}/{sample}_peaks_AllEnhancers_GENE_TO_ENHANCER.txt',
        enhancer=config['workspace'] + '/ROSE/{sample}/{sample}_peaks_AllEnhancers_ENHANCER_TO_GENE.txt'
    input:
        rules.rose.output
    log:
        config['workspace'] + '/log/rose/{sample}_rose_annotate.log'
    run:
        mapper = (
            f'python {BASE_DIR}/tools/ROSE/ROSE_geneMapper.py'
            f' -g {config["genome"]["name"].upper()} -i {input}'
            f' >{log} 2>&1'
        )
        shell(mapper)
