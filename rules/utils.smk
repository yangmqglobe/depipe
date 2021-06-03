# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: utils.smk
# time: 2021/05/24


def find_library(wildcards):
    for sample_name, sample in config['samples'].items():
        for library in sample['treat'].keys():
            if library == wildcards.library:
                return sample_name, 'treat'
        try:
            for library in sample['input'].keys():
                if library == wildcards.library:
                    return sample_name, 'input'
        except KeyError:
            pass


def is_pe(wildcards):
    if getattr(wildcards, 'library', None) is not None:
        sample, condition = find_library(wildcards)
        num = len(config['samples'][sample][condition][wildcards.library])
        return num == 2
    elif getattr(wildcards, 'sample', None) is not None:
        for library in config['samples'][wildcards.sample]['treat'].values():
            return len(library) == 2
