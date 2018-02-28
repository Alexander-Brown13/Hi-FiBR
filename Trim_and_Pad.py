'''
Author: Alexander Brown
Date: 2/27/2018
Purpose: Please refer to "High-throughput analysis of DNA break-induced chromosome rearrangements by amplicon sequencing" in Methods in Enzymology by Brown et al. (2018)
'''

import csv
import sys
import os
from os import listdir
from os.path import isfile, join

def increase_field_size():
    maxInt = sys.maxsize
    decrease = True
    while decrease:
        decrease = False
        try:
            csv.field_size_limit(maxInt)
        except OverflowError:
            maxInt = int(maxInt/10)
            decrease = True
            
def get_files(directory):
    files = [directory + '/' + f for f in listdir(directory) if isfile(join(directory, f))]
    return files

def index_or_error(seq, line_2):
    try:
        seq_crd = line_2.index(seq)
        return seq_crd
    except ValueError:
        seq_crd = -1
        return seq_crd

def rindex_or_error(seq, line_2):
    try:
        seq_crd = line_2.rindex(seq)
        return seq_crd
    except ValueError:
        seq_crd = -1
        return seq_crd

def main():
    increase_field_size()
    directory = input("Enter directory of fastq files: ")
    input_files = get_files(directory)
    
    seq1 = input('Enter first sequence (5\'): ')
    left_pad = input('Enter \'padding\' sequence for this side: ')
    seq2 = input('Enter second sequence (3\'): ')
    right_pad = input('Enter \'padding\' sequence for this side: ')
    lenseq2 = len(seq2)

    len_left_pad = len(left_pad)
    left_pad_score = "G" * len_left_pad
    len_right_pad = len(right_pad)
    right_pad_score = "G" * len_right_pad

    for f in input_files:
        reader = csv.reader(open(f,"r"), delimiter = '\t')
        match_file = open(f[:-6] + "Match.fastq", "w")
        nomatch_file = open(f[:-6] + "NoMatch.fastq", "w")

        for line in reader:
            if line[0][0] != "@":
                continue
            line_2 = next(reader)[0]
            line_3 = next(reader)[0]
            line_4 = next(reader)[0]
            seq1_crd = index_or_error(seq1, line_2)
            seq2_crd = rindex_or_error(seq2, line_2)

            if seq2_crd > seq1_crd and seq1_crd != -1 and seq2_crd != -1:
                line_2_slice = line_2[seq1_crd:seq2_crd+lenseq2]
                line_2_pad = left_pad + line_2_slice + right_pad
                line_4_slice = line_4[seq1_crd:seq2_crd+lenseq2]
                line_4_pad = left_pad_score + line_4_slice + right_pad_score
                line_list = [line[0], line_2_pad, line_3, line_4_pad]
                for out_line in line_list:
                    match_file.write(out_line + '\n')

            else:
                line_list = [line[0], line_2, line_3, line_4]
                for out_line in line_list:
                    nomatch_file.write(out_line + '\n')

main()
