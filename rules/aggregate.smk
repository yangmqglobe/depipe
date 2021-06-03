# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: aggregate.smk
# time: 2021/05/26
rule merge_peaks:
    output:
        config['workspace'] + '/aggregate/all_uniq_peaks.bed'
    input:
        expand(
            config['workspace'] + '/callpeak/{sample}/{sample}_peaks.narrowPeak',
            sample=config['samples']
        )
    shell:
        'cat {input} | sort -k1,1V -k2,2n -k3,3n'
        ' | bedtools merge -i stdin'
        ' | awk \'BEGIN{{OFS="\\t"}}{{print $1, $2, $3, "PEAK_"NR, ".", "."}}\''
        ' > {output}'


rule annotate:
    output:
        config['workspace'] + '/aggregate/all_uniq_peaks_annotated.bed'
    input:
        rules.merge_peaks.output,
        config['genome']['txdb']
    params:
        script=lambda wildcards: BASE_DIR + '/tools/annotate.R'
    shell:
        'Rscript {params.script} {input} {output}'


def find_chrom_sizes(wildcards):
    return [
        config['workspace'] + f'/mapping/{library}/{library}.chrom_sizes'
        for library in config['samples'][wildcards.sample]['treat']
    ]


rule merge_chrom_sizes:
    output:
        temp(config['workspace'] + '/aggregate/{sample}.chrom_sizes')
    input:
        find_chrom_sizes
    shell:
        'cat {input} | sort -k1,1V | uniq > {output}'


def find_fragments(wildcards):
    return [
        config['workspace'] + f'/mapping/{library}/{library}_tag.sorted.bed'
        for library in config['samples'][wildcards.sample]['treat']
    ]


rule count:
    output:
        temp(config['workspace'] + '/aggregate/{sample}_tag_counts.bed')
    input:
        peaks=rules.merge_peaks.output,
        tag=find_fragments,
        chrom_sizes=rules.merge_chrom_sizes.output
    shell:
        'bedtools intersect -wa -c -a {input.peaks} -b {input.tag} -sorted -g {input.chrom_sizes} > {output}'


def find_all_fragment(wildcards):
    return [
        config['workspace'] + f'/mapping/{library}/{library}_tag.sorted.bed'
        for sample in config['samples'].values() for library in sample['treat']
    ]


rule total:
    output:
        temp(config['workspace'] + '/aggregate/all_sample_tag_total.txt')
    input:
        find_all_fragment
    run:
        from snakemake.utils import linecount
        import pandas as pd
        data = pd.Series(list(input), name='path').to_frame()
        samples = pd.Series({
            config['workspace'] + f'/mapping/{library}/{library}_tag.sorted.bed': sample_name
            for sample_name, sample in config['samples'].items() for library in sample['treat']
        }, name='sample')
        data = data.join(samples, on='path')
        data['total'] = data['path'].apply(linecount)
        total = data.groupby('sample')['total'].sum()
        total.to_csv(output[0], sep='\t')



rule merge_counts:
    output:
        raw=config['workspace'] + '/aggregate/all_sample_raw_counts.txt',
        rpkm=config['workspace'] + '/aggregate/all_sample_fpkm.txt',
        qnorm=config['workspace'] + '/aggregate/all_sample_fpkm_qnorm.txt'
    input:
        peaks=rules.merge_peaks.output,
        counts=expand(
            config['workspace'] + '/aggregate/{sample}_tag_counts.bed',
            sample=config['samples']
        ),
        total=rules.total.output
    params:
        names=expand('{sample}', sample=config['samples'])
    run:
        from snakemake.utils import linecount
        import pandas as pd
        import qnorm

        counts = [
            pd.read_table(file, usecols=[3, 6], names=['peak', name], index_col=0, squeeze=True)
            for name, file in zip(params.names, input.counts)
        ]
        counts = pd.concat(counts, axis=1)
        counts.to_csv(output.raw, sep='\t', index_label='')

        peaks = pd.read_table(
            input.peaks[0], usecols=[1, 2, 3], names=['start', 'end', 'peak'], index_col=2
        )
        length = peaks['end'] - peaks['start']
        rpk = counts.divide(length / 1000, axis=0)

        total = pd.read_table(str(input.total), index_col=0, squeeze=True)
        rpkm = rpk.divide(total / 1000000, axis=1)
        rpkm.to_csv(output.rpkm, sep='\t', index_label='')

        norm = qnorm.quantile_normalize(rpkm)
        norm.to_csv(output.qnorm, sep='\t', index_label='')


rule comparisons:
    output:
        config['workspace'] + '/aggregate/comparisons.txt'
    input:
        rules.merge_counts.output
    run:
        import pandas as pd

        meta = {sample: meta['meta'] for sample, meta in config['samples'].items()}
        meta = pd.DataFrame.from_dict(meta, orient='index')

        meta_cols = list({
            comparison['condition'] for comparison in config['comparisons'].values()
        })
        if len(meta_cols) == 0:
            raise ValueError('comparisons not defind!')
        meta = meta[meta_cols]

        meta.to_csv(output[0], sep='\t')


rule pcaplot:
    output:
        config['workspace'] + '/aggregate/all_sample_pcaplot.{fmt}',
    input:
        rules.comparisons.output,
        rules.merge_counts.output.qnorm
    params:
        script=lambda wildcards: BASE_DIR + '/tools/pcaplot.R'
    shell:
        'Rscript {params.script} {input} {output}'
