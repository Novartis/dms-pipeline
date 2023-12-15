#!/usr/bin/env python
# coding: utf-8

import os
import subprocess 
import csv

"""
def get_gene_length(path_to_gene_fasta):
    first_record = SeqIO.read(path_to_gene_fasta, "fasta")
    return(int(len(first_record.seq)))
""" 

def get_gene_length(path_to_gene_fasta):
    with open(path_to_gene_fasta) as csv_file:
        reader = csv.reader(csv_file)
        try:
            seq_name = next(reader)[0]
        except:
            print("""Your file with the reference sequence is not a valid fasta file.
            Please provide a valid fasta file as input.""")
        else:
            if seq_name[0] != '>':
                raise Exception("""Your file with the reference sequence is not a valid fasta file.
                First line should start with \'>\' followed by the name of the sequence.""")
            try:
                seq = next(reader)[0]
                csv_file.close()
            except:
                print("""Your file with the reference sequence is not a valid fasta file.
                Please provide a valid fasta file as input.""")
            else:
                return(int(len(seq)))

with open(snakemake.log[0], 'w') as log_file:

    gene_length = get_gene_length(snakemake.params.gene_fasta)

    callString = "NGmerge -1 %s -2 %s -o %s -f %s -l %s -n %s -d -e %s" % (
        snakemake.input.R1, 
        snakemake.input.R2, 
        snakemake.output.stitched, 
        snakemake.params.unstitched,
        snakemake.params.ngmerge_log, 
        snakemake.threads,
        gene_length
    )

    log_file.write('gene length: %s\n\n' % (gene_length))
    log_file.write('call string: %s\n\n' % (callString))

    ret = subprocess.run(callString, shell=True)
    log_file.write('%s \n' % (ret))


