# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: bam2tag.py
# time: 2021/05/26
from argparse import ArgumentParser
import pysam


def read_regions(filename):
    with open(filename) as f:
        for line in f:
            chrom, start, end = line.strip('\n').split('\t')[:3]
            start, end = int(start), int(end)
            yield chrom, start, end


def fetch_reads(bam, regions=None):
    if regions is None:
        yield from bam
    else:
        regions = read_regions(regions)
        for region in regions:
            for read in bam.fetch(*region):
                yield read


def run_se(in_bam, regions=None, include=0, exclude=1028, mapq=10):
    with pysam.AlignmentFile(in_bam, 'rb') as bam:
        for read in fetch_reads(bam, regions):
            if include != 0 and not read.flag & include:
                continue
            if read.flag & exclude:
                continue
            if read.mapq < mapq:
                continue
            record = (
                read.reference_name, read.reference_start, read.reference_end,
                read.query_name, read.mapq, '-' if read.is_reverse else '+'
            )
            yield record


def run_pe(in_bam, regions=None, include=0, exclude=1028, mapq=10):
    with pysam.AlignmentFile(in_bam, 'rb') as bam:
        for read in fetch_reads(bam, regions):
            if include != 0 and not read.flag & include:
                continue
            if read.flag & exclude:
                continue
            if read.mapq < mapq:
                continue
            if not read.is_proper_pair:
                continue
            if read.template_length <= 0:
                continue
            record = (
                read.reference_name, read.reference_start,
                read.reference_start + read.template_length,
                read.query_name, read.mapq, '.'
            )
            yield record


def main():
    parser = ArgumentParser(
        description='Extract fragments from bam file'
    )
    parser.add_argument(
        '-p', dest='paired', help='use parired end model', action='store_true'
    )
    parser.add_argument(
        '-L', dest='regions', help='only include reads overlapping this BED FILE', default=None
    )
    parser.add_argument(
        '-f', dest='include', help='only include reads with all  of the FLAGs in INT present',
        default=0
    )
    parser.add_argument(
        '-F', dest='exclude', help='only include reads with none of the FLAGS in INT present',
        default=1028
    )
    parser.add_argument(
        '-q', dest='mapq', help='only include reads with mapping quality >= INT', default=10, type=int
    )
    parser.add_argument(
        '-l', dest='length', help='maximum length of fragment', default=1500, type=int
    )
    parser.add_argument('<in_bam>', help='input bam file')
    parser.add_argument('<out_bed>', help='output bed file')
    args = vars(parser.parse_args())
    if args['paired']:
        records = run_pe(
            args['<in_bam>'], args['regions'], args['include'], args['exclude'], args['mapq']
        )
    else:
        records = run_se(
            args['<in_bam>'], args['regions'], args['include'], args['exclude'], args['mapq']
        )
    with open(args['<out_bed>'], 'w') as bed:
        for record in records:
            if record[2] - record[1] < args['length']:
                bed.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(*record))


if __name__ == '__main__':
    main()
