# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: preprocess.smk
# time: 2021/05/24
import os


def find_raw_fastq(wildcards):
    sample, condition = find_library(wildcards)
    return config['samples'][sample][condition][wildcards.library]


rule trim:
    output:
        fq1=config['workspace'] + '/preprocess/{library}/{library}_r1.fq.gz',
        json=config['workspace'] + '/qc/{library}_fastp.json',
        html=config['workspace'] + '/qc/{library}_fastp.html'
    input:
        find_raw_fastq
    log:
        config['workspace'] + '/log/preprocess/{library}_fastp.log'
    threads:
        8 if workflow.cores > 8 else workflow.cores
    run:
        if len(input) == 1:
            fqin = f'-i {input[0]}'
            fqout = f'-o {output.fq1}'
        elif len(input) == 2:
            fqin = f'-i {input[0]} -I {input[1]}'
            fqout = f'-o {output.fq1} -O {output.fq1.replace("r1.fq.gz", "r2.fq.gz")}'
        else:
            raise ValueError(f'Fastq file number for library {wildcards.library} must be 1 or 2')
        shell(f'fastp -w {threads} {fqin} {fqout} -j {output.json} -h {output.html} 2>{log}')