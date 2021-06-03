# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: callpeak.smk
# time: 2021/04/16


def find_treat_tags(wildcards):
    return [
       config['workspace'] + f'/mapping/{library}/{library}_tag.sorted.bed'
       for library in config['samples'][wildcards.sample]['treat']
    ]


def find_control_tags(wildcards):
    if 'input' in config['samples'][wildcards.sample]:
        return [
            config['workspace'] + f'/mapping/{library}/{library}_tag.sorted.bed'
            for library in config['samples'][wildcards.sample]['input']
        ]
    return list()


def find_chrom_sizes(wildcards):
    for library in config['samples'][wildcards.sample]['treat']:
        return config['workspace'] + f'/mapping/{library}/{library}.chrom_sizes'


rule callpeak:
    output:
        narrowPeak=config['workspace'] + '/callpeak/{sample}/{sample}_peaks.narrowPeak',
        treat_bdg=config['workspace'] + '/callpeak/{sample}/{sample}_treat_pileup.bdg'
    input:
        treat=find_treat_tags,
        control=find_control_tags
    log:
        config['workspace'] + '/log/callpeak/{sample}_macs2.log'
    params:
        genome_size=config['genome']['size'],
        outdir=config['workspace'] + '/callpeak/{sample}',
        name='{sample}'
    run:
        macs2 = (f'macs2 callpeak -B --SPMR --keep-dup all -q {config["qvalue"]}'
                f' -g {params.genome_size} --outdir {params.outdir} -n {params.name}'
                f' -t {" ".join(input.treat)}')
        if input.control:
            macs2 += f' -c {" ".join(input.control)}'
        if is_pe(wildcards):
            macs2 += ' -f BEDPE'
        else:
            macs2 += ' -f BED'
        macs2 += f' 2>{log}'
        shell(macs2)


rule bdg_clip:
    output:
        temp(config['workspace'] + '/callpeak/{sample}/{sample}_treat_pileup_clip.bdg')
    input:
        rules.callpeak.output.treat_bdg,
        find_chrom_sizes
    shell:
        'bedClip {input} {output}'


rule bdg_sort:
    output:
        temp(config['workspace'] + '/callpeak/{sample}/{sample}_treat_pileup_sort.bdg')
    input:
        rules.bdg_clip.output
    shell:
        'bedSort {input} {output}'


rule bdg2bw:
    output:
        config['workspace'] + '/callpeak/{sample}/{sample}_treat_pileup.bw'
    input:
        rules.bdg_sort.output,
        find_chrom_sizes
    shell:
        'bedGraphToBigWig {input} {output}'