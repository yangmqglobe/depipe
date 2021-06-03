# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: comparison.smk
# time: 2021/05/27
import os


rule plan:
    output:
        config['workspace'] + '/comparisons/{comparison}/{comparison}_plan.txt'
    run:
        import pandas as pd
        # construct metadata
        meta = {sample: meta['meta'] for sample, meta in config['samples'].items()}
        meta = pd.DataFrame.from_dict(meta, orient='index')
        # filter out unuse samples
        use = pd.Series(True, index=meta.index)
        if 'exclude' in config['comparisons'][wildcards.comparison]:
            for key, value in config['comparisons'][wildcards.comparison]['exclude'].items():
                if isinstance(value, list):
                    use = use & (~meta[key].isin(value))
                else:
                    use = use & (meta[key] != value)
        # construct comparison
        condition = config['comparisons'][wildcards.comparison]['condition']
        numerator = config['comparisons'][wildcards.comparison]['numerator']
        denominator = config['comparisons'][wildcards.comparison]['denominator']
        comparison = pd.Series('unuse', index=meta.index, name='condition')
        comparison[(meta[condition] == numerator) & use] = 'numerator'
        comparison[(meta[condition] == denominator) & use] = 'denominator'
        comparison.to_csv(output[0], sep='\t')


rule compare:
    output:
        config['workspace'] + '/comparisons/{comparison}/{comparison}_result.txt'
    input:
        rules.plan.output,
        rules.merge_counts.output.raw
    params:
        script=lambda wildcards: BASE_DIR + '/tools/compare.R'
    shell:
        'Rscript {params.script} {input} {output}'
        