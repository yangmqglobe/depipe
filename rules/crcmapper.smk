# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: crcmapper.smk
# time: 2021/06/08
import os


def find_treat_bams(wildcards):
    return [
       config['workspace'] + f'/mapping/{library}/{library}.sorted.bam'
       for library in config['samples'][wildcards.sample]['treat']
    ]


rule crcmapper:
    output:
        crc=config['workspace'] + '/CRC/crcmapper/{sample}/{sample}_CRC_SCORES.txt',
        enhancer=config['workspace'] + '/CRC/crcmapper/{sample}/{sample}_SE_ASSIGNMENT_GENE.txt'
    input:
        enhancer=rules.rose.output,
        treat=find_treat_bams,
        peak=rules.callpeak.output.narrowPeak
    log:
        config['workspace'] + '/log/crcmapper/{sample}_crcmapper.log'
    run:
        outdir = os.path.dirname(str(output))
        crcmapper = (
            f'python {BASE_DIR}/tools/CRCmapper/CRCmapper.py'
            f' -e {input.enhancer} -g {config["genome"]["name"]} -f {config["genome"]["fasta"]}'
            f' -b {",".join(input.treat)} -s {input.peak} -n {wildcards.sample} -o {outdir}'
            f' >{log} 2>&1'
        )
        shell(crcmapper)


rule crcmapper_downstream:
    output:
        config['workspace'] + '/CRC/crcmapper/{sample}/downstream/{sample}_crc_regulate.txt'
    input:
        crc=rules.crcmapper.output.crc,
        enhancer=rules.crcmapper.output.enhancer,
        peak=rules.callpeak.output.narrowPeak
    threads:
        workflow.cores
    run:
        subpeaks = f'{config["workspace"]}/CRC/crcmapper/{wildcards.sample}/downstream/{wildcards.sample}_subpeaks.txt'
        bedtools = (
            'awk \'BEGIN{{OFS="\\t"}}{{if($3 < $2){{temp=$2;$2=$3;$3=temp}};print $0}}\' {input.enhancer}'
            f' | bedtools intersect -wao -a stdin -b {input.peak} > {subpeaks}'
        )
        shell(bedtools)
        downstream = (
            f'python {BASE_DIR}/tools/CRCmapper/downstream.py -c {input.crc}'
            f' -o {config["workspace"]}/CRC/crcmapper/{wildcards.sample}/downstream'
            f' -n {wildcards.sample} -p {subpeaks} -f {config["genome"]["fasta"]}'
        )
        shell(downstream)

