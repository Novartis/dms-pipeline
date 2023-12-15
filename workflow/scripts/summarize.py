#!/usr/bin/env python
# coding: utf-8

# In[4]:



'''
Creates a summary table of the number of reads that passed or failed each step of the DMS pipeline for all samples.


'''

import os
import sys
import csv
import numpy as np


def count_seq_fastqgz(file):
    cmd =  'echo $(zcat %s|wc -l)/4|bc' % (file)
    return (int(os.popen(cmd).read().split()[0]))


def formatLine(l):
    return('\t'.join(['%s' % (x) for x in l]) + '\n')

log_file = open(snakemake.log.log1,"w")
sys.stdout = log_file

samples = snakemake.params.samples
OUTPUT_DIR = snakemake.params.output_dir
expected_mutants = False

summary = []

out_log = ('Number of reads remaining after each bioinformatic step\n'
           '-----------------------------------\n\n'
           'Output are:\n'
           'absolute number, percentage wrt total number of sequenced reads, percentage wrt absolute number in previous step\n\n')

for i, s in enumerate(samples):
    
    # 1.1) count number of reads
    total_reads = min(count_seq_fastqgz('%s/%s/file_R1.fastq.gz' % (OUTPUT_DIR, s)),
                      count_seq_fastqgz('%s/%s/file_R2.fastq.gz' % (OUTPUT_DIR, s)))

    # 1.2) count number of stitched reads
    stitched_reads= count_seq_fastqgz('%s/%s/stitched.fastq.gz' % (OUTPUT_DIR, s))
    stitched_discarded = total_reads - stitched_reads

    # 1.3) count number of mapped stitched reads
    # 1.4) count number of expected reads
    exp_reads = 0
    mapped_reads = 0
    with open(snakemake.input.counts[i]) as csv_file:
        reader = csv.reader(csv_file)
        header = np.array(next(reader))
        if 'expected' in header:
            expected_mutants = True
        for row in reader:
            mapped_reads = mapped_reads + int(row[int(np.where(header == 'count')[0])])
            if expected_mutants:
                exp_reads = exp_reads + int(row[int(np.where(header == 'count')[0])]) *                 int(row[int(np.where(header == 'expected')[0])] == 'True')
        csv_file.close()
    
    unmapped_reads = stitched_reads - mapped_reads
    unexpected_reads = mapped_reads - exp_reads
    
    summary.append([s,
                    total_reads, 
                    stitched_reads,
                    stitched_discarded,
                    round(stitched_reads / total_reads * 100, 2), 
                    mapped_reads,
                    unmapped_reads,
                    round(mapped_reads / total_reads * 100, 2),
                    round(mapped_reads / stitched_reads * 100, 2),
                    exp_reads,
                    unexpected_reads,
                    round(exp_reads / total_reads * 100, 2),
                    round(exp_reads / mapped_reads * 100, 2)])
    
    out_log = out_log +     '\nsample: %s\n' % (s) +     'sequenced: %s\n' % (total_reads) +     'stitched: %s\t%s%%\n' % (stitched_reads, 
                              round(stitched_reads / total_reads *100, 2)) + \
    '  discarded in stitch: %s\t%s%%\n' % (stitched_discarded, 
                                           round(stitched_discarded / total_reads *100, 2)) + \
    'mapped to wt sequence: %s\t%s%%\t%s%%\n' % (mapped_reads,
                                                round(mapped_reads / total_reads *100, 2),
                                                round(mapped_reads / stitched_reads *100, 2)) + \
    '  stitched but unmapped: %s\t%s%%\t%s%%\n' % (unmapped_reads,
                                                  round(unmapped_reads / total_reads *100, 2),
                                                  round(unmapped_reads / stitched_reads *100, 2))
    
    if expected_mutants:
        out_log = out_log +         'expected amplicon: %s\t%s%%\t%s%%\n' % (exp_reads, 
                                                round(exp_reads / total_reads *100, 2),
                                                round(exp_reads / mapped_reads *100, 2)) + \
        '  unexpected amplicon: %s\t%s%%\t%s%%\n' % (unexpected_reads, 
                                                    round(unexpected_reads / total_reads *100, 2),
                                                    round(unexpected_reads / mapped_reads *100, 2))

# 1.5) create summary table
with open(snakemake.output.summary, 'w') as fOut:
    header = ['sample',
              'total_reads', 
              'stitched_reads',
              'discarded_in_stitch',
              'perc_stitched_of_total',
              'mapped_reads', 
              'unmapped_reads', 
              'perc_mapped_of_total', 
              'perc_mapped_of_stitched',
              'expected_reads', 
              'unexpected_reads', 
              'perc_expected_of_total', 
              'perc_expected_of_mapped']
    if expected_mutants:
        fOut.write(formatLine(header))
        for l in summary:
            fOut.write(formatLine(l))
    else:
        fOut.write(formatLine(header[:-4]))
        for l in summary:
            fOut.write(formatLine(l[:-4]))
    fOut.close()
            
print(out_log)
log_file.close()

