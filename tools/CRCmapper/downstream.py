# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: downstream.py
# time: 2021/06/10
from multiprocessing import TimeoutError
from multiprocessing.dummy import Pool
from collections import defaultdict
from functools import partial, reduce
from argparse import ArgumentParser
from pyfaidx import Fasta
import pandas as pd
import subprocess
import sys
import os


def call_fimo(bgfile, pwm, fasta, motif, fimo='fimo'):
    args = [
        fimo, '--motif', motif, '-verbosity', '1', '--bgfile',
        bgfile, '-text', pwm, fasta
    ]
    p = subprocess.check_output(args)
    return p.decode()


def run(crc, motif, genome, peak, out, name, processes, bin, pwm):
    crc = pd.read_table(crc, names=['crc', 'score', 'count'])
    crc['crc'] = crc['crc'].apply(eval)

    crc_genes = reduce(set.union, crc['crc'].apply(set))

    motif = pd.read_table(motif, names=['motif', 'tf'], index_col=0)

    crc_motifs = motif.loc[motif['tf'].isin(crc_genes)].index.values

    genome = Fasta(genome)

    peak = pd.read_table(
        peak, usecols=(0, 1, 2, 3, 5, 6), names=[
            'chrom', 'se_start', 'se_end', 'se_gene', 'peak_start', 'peak_end'
        ]
    )

    fasta = os.path.join(out, f'{name}_peaks.fa')
    with open(fasta, 'w') as f:
        for _, chrom, se_start, se_end, se_gene, peak_start, peak_end in peak.itertuples():
            seq = genome[chrom][peak_start - 200:peak_end + 200].seq.upper()
            print(
                f'>{chrom}|{se_start}|{se_end}|{se_gene}|200|{peak_start}|{peak_end}\n{seq}', file=f)

    bgfile = os.path.join(out, f'{name}_peaks_bg.meme')
    subprocess.check_call([
        f'{bin}/fasta-get-markov', '-m', '1', fasta, bgfile
    ])

    fimo_process = partial(call_fimo, bgfile, pwm, fasta, fimo=f'{bin}/fimo')
    binding = os.path.join(out, f'{name}_fimo.txt')
    with open(binding, 'w') as f:
        with Pool(processes=processes) as pool:
            results = [
                pool.apply_async(fimo_process, (motif,)) for motif in crc_motifs
            ]
            while results:
                result = results.pop()
                try:
                    output = result.get(1)
                    for line in output.splitlines()[1:]:
                        print(line, file=f)
                except TimeoutError:
                    results.append(result)

    fimo = pd.read_table(binding, names=[
        'motif_id', 'gene', 'sequence_name', 'start', 'stop', 'strand',
        'score', 'p-value', 'q-value', 'na', 'matched_sequence'
    ])
    fimo = fimo[fimo['score'] > 10].copy()
    fimo = fimo.join(motif['tf'], on='motif_id')
    origin = pd.DataFrame(
        fimo['sequence_name'].str.split('|').values.tolist(),
        columns=[
            'chrom', 'se_start', 'se_end', 'se_gene', 'extend', 'peak_start', 'peak_end'
        ]
    )
    origin = origin.astype({
        'se_start': int, 'se_end': int, 'extend': int, 'peak_start': int, 'peak_end': int
    })
    fimo = fimo.merge(origin, left_index=True, right_index=True)
    fimo['abs_start'] = fimo['peak_start'] - fimo['extend'] + fimo['start']
    fimo['abs_end'] = fimo['peak_start'] - fimo['extend'] + fimo['stop']
    fimo = fimo[
        (fimo['abs_start'] >= fimo['peak_start'])
        & (fimo['abs_end'] <= fimo['peak_end'])
    ].copy()
    fimo = fimo.drop_duplicates(['chrom', 'abs_start', 'abs_end', 'tf'])
    count = fimo.groupby(['se_gene', 'tf']).size()
    count = count[count >= 3]

    se_binding = count.reset_index().groupby('se_gene')['tf'].apply(set)
    crc_bind_to = defaultdict(set)
    for se, tf in se_binding.items():
        for i, genes, *_ in crc.itertuples():
            if len(set(genes) - tf) == 0:
                crc_bind_to[i].add(se)
    crc_bind_to = pd.Series({
        i: ','.join(genes) for i, genes in crc_bind_to.items()
    })

    crc = crc.join(crc_bind_to.rename('regulate'))
    crc['crc_genes'] = crc['crc'].apply(','.join)

    crc[['crc_genes', 'score', 'count', 'regulate']].to_csv(
        os.path.join(out, f'{name}_crc_regulate.txt'), index=False, sep='\t'
    )


def main():
    parser = ArgumentParser(
        description='analysis CRC bind SE-associate genes'
    )
    parser.add_argument(
        '-c', dest='crc', help='crc score output from crcmapper', required=True
    )
    parser.add_argument(
        '-f', dest='genome', help='samtools indexed genome fasta file', required=True
    )
    parser.add_argument(
        '-p', dest='peak', help='narrow peak output from macs2', required=True
    )
    parser.add_argument(
        '-o', dest='out', help='output dir', required=True
    )
    parser.add_argument(
        '-n', dest='name', help='output file prefix', required=True
    )
    parser.add_argument(
        '-b', dest='bin', help='bin dictionary which contains fimo', default=None
    )
    parser.add_argument(
        '-t', dest='processes', help='maximum process to use', type=int, default=0
    )
    args = vars(parser.parse_args())

    cwd = os.path.dirname(os.path.abspath(__file__))

    args['motif'] = os.path.join(cwd, 'MotifDictionary.txt')
    args['pwm'] = os.path.join(cwd, 'VertebratePWMs.txt')

    if args['bin'] is None:
        args['bin'] = os.path.dirname(sys.executable)

    if args['processes'] == 0:
        args['processes'] = os.cpu_count()
    else:
        args['processes'] = min(os.cpu_count(), args['process'])

    run(**args)


if __name__ == '__main__':
    main()
