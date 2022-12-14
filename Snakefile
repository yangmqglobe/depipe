# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: Snakefile
# time: 2021/05/24
include: 'rules/utils.smk'
include: 'rules/preprocess.smk'
include: 'rules/mapping.smk'
include: 'rules/callpeak.smk'
include: 'rules/aggregate.smk'
include: 'rules/comparison.smk'
include: 'rules/rose.smk'
include: 'rules/crcmapper.smk'

import os


BASE_DIR = os.path.dirname(workflow.snakefile)


rule basic_all:
    input:
        expand(
            config['workspace'] + '/callpeak/{sample}/{sample}_peaks.bed',
            sample=config['samples']
        ),
        expand(
            config['workspace'] + '/callpeak/{sample}/{sample}_treat_pileup.bw',
            sample=config['samples']
        ),
        expand(
            config['workspace'] + '/callpeak/{sample}/{sample}_control_lambda.bw',
            sample=config['samples']
        ),
        expand(
            config['workspace'] + '/callpeak/{sample}/{sample}_peaks_annotated.txt',
            sample=config['samples']
        ),
        config['workspace'] + '/aggregate/all_uniq_peaks_annotated.bed',
        config['workspace'] + '/aggregate/all_sample_fpkm_qnorm.txt'


rule plot_all:
    input:
        expand(
            config['workspace'] + '/aggregate/all_sample_pcaplot.{fmt}',
            fmt=config['plot_formats']
        )


rule compare_all:
    input:
        expand(
            config['workspace'] + '/comparisons/{comparison}/{comparison}_result.txt',
            comparison=config['comparisons'] if 'comparisons' in config else list()
        )


rule rose_all:
    input:
        expand(
            config['workspace'] + '/ROSE/{sample}/{sample}_peaks_AllEnhancers_GENE_TO_ENHANCER.txt',
            sample=config['samples']
        )


rule crcmapper_all:
    input:
        expand(
            config['workspace'] + '/CRC/crcmapper/{sample}/{sample}_CRC_SCORES.txt',
            sample=config['samples']
        ),
        expand(
            config['workspace'] + '/CRC/crcmapper/{sample}/downstream/{sample}_crc_regulate.txt',
            sample=config['samples']
        )


rule all:
    input:
        rules.basic_all.input,
        rules.plot_all.input,
        rules.compare_all.input,
        rules.rose_all.input,
        rules.crcmapper_all.input