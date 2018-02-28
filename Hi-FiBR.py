'''
Author: Alexander Brown
Date: 2/27/2018
Purpose: Please refer to "High-throughput analysis of DNA break-induced chromosome rearrangements by amplicon sequencing" in Methods in Enzymology by Brown et al. (2018)
'''

import csv
import sys
import os
from os import listdir, getcwd
from os.path import isfile, join, isdir

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

def retrieve_ref_seq(reference):
    reader = csv.reader(open(reference, "r"), delimiter = '\t')
    header = next(reader)
    ref_seq = next(reader)[0].strip()
    return ref_seq

def make_cigar_substring_list(cigar_string):
    MDI_list = ["M", "D", "I"]
    cigar_substring_list = []
    cigar_len = len(cigar_string)
    i = 0
    while i < cigar_len:
        sub_cigar = ""
        while cigar_string[i] not in MDI_list:
            sub_cigar = sub_cigar + cigar_string[i]
            i += 1
        sub_cigar = sub_cigar + cigar_string[i]
        cigar_substring_list.append(sub_cigar)
        i += 1
    return cigar_substring_list

def report_deletion_insertion_info(cigar_substring_list, breakpoint_from_left, breakpoint_from_right,
                                   read_length, read_seq, ref_seq):

    M_left = int(cigar_substring_list[0][:-1])
    M_right = int(cigar_substring_list[-1][:-1])

    dist_from_brk_left = M_left - breakpoint_from_left
    dist_from_brk_right = M_right - breakpoint_from_right

    deletion_left = dist_from_brk_left if dist_from_brk_left < 0 else 0
    deletion_right = dist_from_brk_right if dist_from_brk_right < 0 else 0
    deletion_total = deletion_left + deletion_right

    if deletion_left < 0:
        insert_start = M_left
    else:
        insert_start = breakpoint_from_left

    if deletion_right < 0:
        insert_end = M_right
    else:
        insert_end = breakpoint_from_right

    insert_length = read_length - insert_start - insert_end
    insert_length = 0 if insert_length < 0 else insert_length
    insert_seq = read_seq[insert_start:-insert_end]

    return M_left, M_right, dist_from_brk_left, dist_from_brk_right, deletion_left, deletion_right, deletion_total, insert_start, insert_end, insert_length, insert_seq

def check_cigar_for_fake_insertions(cigar_substring_list, breakpoint_from_left, breakpoint_from_right,
                                   read_length, read_seq, ref_seq):
    M_left, M_right, dist_from_brk_left, dist_from_brk_right, deletion_left, deletion_right, deletion_total, insert_start, insert_end, insert_length, insert_seq = report_deletion_insertion_info(cigar_substring_list, breakpoint_from_left, breakpoint_from_right, read_length, read_seq, ref_seq)
    m_left_sub_cigar = cigar_substring_list[0]
    m_right_sub_cigar = cigar_substring_list[-1]
    if m_left_sub_cigar.find("M") != -1 and m_right_sub_cigar.find("M") != -1 and deletion_total < 0 and insert_length > 0:
        del_start_index = breakpoint_from_left + deletion_left
        del_end_index = breakpoint_from_right + deletion_right
        del_seq = ref_seq[del_start_index:-del_end_index]

        if len(cigar_substring_list) != 4:
            if 0 > deletion_left and 0 > deletion_right:
                i = 0
                try:
                    while del_seq[i] == insert_seq[i]:
                        i += 1
                except IndexError:
                    pass
                M_left = str(int(m_left_sub_cigar[:-1]) + i) + "M"
                i = 1
                try:
                    while del_seq[-i] == insert_seq[-i]:
                        i += 1
                except IndexError:
                    pass
                i -= 1
                M_right = str(int(m_right_sub_cigar[:-1]) + i) + "M"
                new_cigar_substring_list = [M_left] + cigar_substring_list[1:-1] + [M_right]
                return new_cigar_substring_list

            if 0 > deletion_left:
                i = 0
                try:
                    while del_seq[i] == insert_seq[i]:
                        i += 1
                except IndexError:
                    pass
                M_left = str(int(m_left_sub_cigar[:-1]) + i) + "M"
                new_cigar_substring_list = [M_left] + cigar_substring_list[1:-1] + [m_right_sub_cigar]
                return new_cigar_substring_list

            if 0 > deletion_right:
                i = 1
                try:
                    while del_seq[-i] == insert_seq[-i]:
                        i += 1
                except IndexError:
                    pass
                i -= 1
                M_right = str(int(m_right_sub_cigar[:-1]) + i) + "M"
                new_cigar_substring_list = [m_left_sub_cigar] + cigar_substring_list[1:-1] + [M_right]
                return new_cigar_substring_list

        substring_index = 0
        for cigar_substring in cigar_substring_list:
            if cigar_substring.find("D") != -1:
                del_substring_index = substring_index
                del_sub_cigar = cigar_substring_list[del_substring_index]
            if cigar_substring.find("I") != -1:
                insert_substring_index = substring_index
                insert_sub_cigar = cigar_substring_list[insert_substring_index]
            substring_index += 1

        if insert_substring_index > del_substring_index:
            i = 0
            try:
                while del_seq[i] == insert_seq[i]:
                    i += 1
            except IndexError:
                pass
            M_left = str(int(m_left_sub_cigar[:-1]) + i) + "M"
            new_cigar_substring_list = [M_left, del_sub_cigar, insert_sub_cigar, m_right_sub_cigar]
            return new_cigar_substring_list

        if insert_substring_index < del_substring_index:
            i = 1
            try:
                while del_seq[-i] == insert_seq[-i]:
                    i += 1
            except IndexError:
                pass
            i -= 1
            M_right = str(int(m_right_sub_cigar[:-1]) + i) + "M"
            new_cigar_substring_list = [m_left_sub_cigar, insert_sub_cigar, del_sub_cigar,  M_right]
            return new_cigar_substring_list
    else:
        return cigar_substring_list

def check_cigar_for_non_temp_insertions(cigar_substring_list, breakpoint_from_left, breakpoint_from_right, left_count_threshold, right_count_threshold, read_seq, ref_seq):
    m_left_sub_cigar = cigar_substring_list[0]
    m_right_sub_cigar = cigar_substring_list[-1]
    M_left = m_left_sub_cigar
    M_right = m_right_sub_cigar
    exact_ticker = 0
    if len(cigar_substring_list) == 1:
        if len(ref_seq) > len(read_seq):
            return cigar_substring_list
        exact_ticker = 1
        m_left_sub_cigar = str(breakpoint_from_left) + "M"
        m_right_sub_cigar = str(breakpoint_from_right) + "M"
        M_left = m_left_sub_cigar
        M_right = m_right_sub_cigar
        
    if m_left_sub_cigar.find("M") != -1 and m_right_sub_cigar.find("M") != -1:
        left_match_index = int(m_left_sub_cigar[:-1])
        left_count = 1
        left_mismatch_list = []
        while left_count < left_count_threshold:
            if read_seq[left_match_index - left_count] != ref_seq[left_match_index - left_count]:
                mismatch_index = left_match_index - left_count + 1
                left_mismatch_list.append(mismatch_index)
            left_count += 1
        if len(left_mismatch_list) > 0:
            min_mismatch_index = min(left_mismatch_list)
            M_left = str(min_mismatch_index - 1) + "M"

        right_match_index = int(m_right_sub_cigar[:-1])
        right_count = 0
        right_mismatch_list = []
        while right_count < right_count_threshold:
            if read_seq[-(right_match_index - right_count)] != ref_seq[-(right_match_index - right_count)]:
                mismatch_index = right_match_index - right_count
                right_mismatch_list.append(mismatch_index)
            right_count += 1
        if len(right_mismatch_list) > 0:
            min_mismatch_index = min(right_mismatch_list)
            M_right = str(min_mismatch_index - 1) + "M"

        if int(M_left[:-1]) + int(M_right[:-1]) == len(ref_seq) and exact_ticker == 1:
            return cigar_substring_list

        new_cigar_substring_list = [M_left] + cigar_substring_list[1:-1] + [M_right]
        return new_cigar_substring_list
    else:
        return cigar_substring_list

def report_class(deletions, insert_length):
    if deletions == 0:
        if insert_length == 0:
            return 'exact'
        else:
            return 'insertion'
    else:
        if insert_length == 0:
            return 'deletion'
        else:
            return 'complex'

def get_working_sequence(read_length, ref_seq_length, sequence_class, reconstructed_read, read_seq, ref_seq):
    if read_length != ref_seq_length and sequence_class != "exact":
        working_sequence = reconstructed_read
    if read_length != ref_seq_length and sequence_class == "exact":
        working_sequence = read_seq
    if read_length == ref_seq_length and sequence_class != "exact":
        working_sequence = reconstructed_read
    if read_length == ref_seq_length and sequence_class == "exact":
        working_sequence = ref_seq
    return working_sequence

def prepare_line(new_line_list):
    new_line_list = [str(s) for s in new_line_list]
    new_line = '\t'.join(new_line_list)
    return new_line

def remove_duplicates(alist):
    return list(set(alist))
        
def report_seq_micro(seq, read_micro_index):
    try:
        read_micro = seq[read_micro_index]
        return read_micro
    except IndexError:
        read_micro = ''
        return read_micro

def report_microhomologies(item, ref_seq, breakpoint_from_right, breakpoint_from_left):
    side_list = ["left", "right"]
    for side in side_list:
        if side == "left":
            read_micro_index = int(item[4]) - 1 #left
            read_micro_index_2 = read_micro_index #left
            ref_micro_index = -(int(item[5]) + 1) #right
            ref_micro_index_2 = ref_micro_index #right
            if_statement = ref_micro_index_2 < -(breakpoint_from_right)
            increment = -1
        if side == "right":
            read_micro_index = -int(item[5]) #right
            read_micro_index_2 = read_micro_index #right
            ref_micro_index = int(item[4]) #left
            ref_micro_index_2 = ref_micro_index #left
            if_statement = ref_micro_index_2 > breakpoint_from_left-1
            increment = 1

        read_micro = report_seq_micro(item[-1], read_micro_index)
        ref_micro = report_seq_micro(ref_seq, ref_micro_index)

        while read_micro == ref_micro:
            if if_statement:
                break
            read_micro_index_2 += increment
            ref_micro_index_2 += increment
            read_micro = report_seq_micro(item[-1], read_micro_index_2)
            ref_micro = report_seq_micro(ref_seq, ref_micro_index_2)
        if side == "left":
            read_micro_homology_left = item[-1][read_micro_index_2+1:read_micro_index+1]
            ref_micro_homology_right = ref_seq[ref_micro_index_2+1:ref_micro_index+1]
        if side == "right":
            read_micro_homology_right = item[-1][read_micro_index:read_micro_index_2]
            ref_micro_homology_left = ref_seq[ref_micro_index:ref_micro_index_2]
    return read_micro_homology_left, ref_micro_homology_right, read_micro_homology_right, ref_micro_homology_left

def write_fasta_line(name_list, out_sequence, out_file):
    fasta_name = '>' + '_'.join(name_list)
    out_file.write(fasta_name + '\n')
    out_file.write(out_sequence + '\n')

def main():
    increase_field_size()
    reference = input("Enter reference sequence file(fasta format): ")
    breakpoint_from_left = int(input("Enter distance to breakpoint from 5' (left) side: "))
    breakpoint_from_right = int(input("Enter distance to breakpoint from 3' (right) side: "))

    ref_seq = retrieve_ref_seq(reference)
    ref_seq_length = len(ref_seq)
    directory = input("Enter directory of SAM files (make sure file names end with '.sam'): ")
    input_files = get_files(directory)

    left_count_threshold = int(input("Enter mismatch window from break 5' (left) side: ")) + 1
    right_count_threshold = int(input("Enter mismatch window from break 3' (right) side: "))
    
    for f in input_files:
        if f.find(".sam") == -1:
            continue
        reader = csv.reader(open(f,"r"), dialect="excel-tab")
        out = open(f[:-4] + "_Extended.sam", "w")
        recon_list = []
        no_duplicate_list = []
        read_error_dict = {}
        read_error_distribution_dict = {}
        for line in reader:
            if line[0][0] == '@':
                continue
            cigar_string = line[5]
            blacklisted = ["S", "*", "H"] #S appears in CLC, * in Bowtie2, H in BWA
            blacklisted_ticker = 0
            for item in blacklisted:
                if cigar_string.find(item) != -1:
                    blacklisted_ticker = 1
                    break
            if blacklisted_ticker != 0:
                continue
            read_seq = line[9]
            read_length = len(read_seq)
            
            cigar_substring_list = make_cigar_substring_list(cigar_string)
            cigar_substring_list = check_cigar_for_fake_insertions(cigar_substring_list, breakpoint_from_left, breakpoint_from_right, read_length, read_seq, ref_seq)
            cigar_substring_list = check_cigar_for_non_temp_insertions(cigar_substring_list, breakpoint_from_left, breakpoint_from_right, left_count_threshold, right_count_threshold, read_seq, ref_seq)

            M_left, M_right, dist_from_brk_left, dist_from_brk_right, deletion_left, deletion_right, deletion_total, insert_start, insert_end, insert_length, insert_seq = report_deletion_insertion_info(cigar_substring_list, breakpoint_from_left, breakpoint_from_right, read_length, read_seq, ref_seq)
            sequence_class = report_class(deletion_total, insert_length)

            left_end_index = breakpoint_from_left + deletion_left
            right_start_index = breakpoint_from_right + deletion_right

            left_half_ref = ref_seq[:left_end_index]
            right_half_ref = ref_seq[-right_start_index:]

            reconstructed_read = left_half_ref + insert_seq + right_half_ref
            working_sequence = get_working_sequence(read_length, ref_seq_length, sequence_class, reconstructed_read, read_seq, ref_seq)
            if len(read_seq) != len(working_sequence):
                continue#This should not occur often; it is an additional precaution for junk reads

            recon_list.append(working_sequence)
            if working_sequence not in read_error_dict:
                read_error_dict[working_sequence] = 0
                read_error_distribution_dict[working_sequence] = {}
                bp = 0
                while bp < len(working_sequence):
                    read_error_distribution_dict[working_sequence][bp + 1] = 0
                    bp += 1
            if read_seq != working_sequence:
                read_error_dict[working_sequence] += 1
                bp_position = 0
                while bp_position < len(working_sequence):
                    read_bp = read_seq[bp_position]
                    recon_bp = working_sequence[bp_position]
                    if read_bp != recon_bp:
                        read_error_distribution_dict[working_sequence][bp_position + 1] += 1
                    bp_position += 1

            additional_line_list = [read_length, cigar_string, M_left, M_right, dist_from_brk_left,
                                    dist_from_brk_right, deletion_left, deletion_right, deletion_total, insert_start,
                                    insert_end, insert_length, insert_seq, sequence_class, working_sequence]
            new_line_list = line.copy()
            new_line_list.extend(additional_line_list)
            new_line = prepare_line(new_line_list)
            out.write(new_line + '\n')

            no_duplicate_line = [line[2], line[5]]
            no_duplicate_line.extend(additional_line_list)
            no_duplicate_line = prepare_line(no_duplicate_line)
            no_duplicate_list.append(no_duplicate_line)

        no_duplicate_list = remove_duplicates(no_duplicate_list)
        split_list = []
        for item in no_duplicate_list:
            split_list.append(item.split("\t"))
        no_duplicate_list = split_list

        micro_out = open(f[:-4] + "_Count_Columns.sam", "w")
        control_list = []
        sorting_list = []
        for item in no_duplicate_list:
            w = 1
            line_count = recon_list.count(item[-1])
            control_list.append(item[-1])

            if item[-2] == "deletion":
                read_micro_homology_left, ref_micro_homology_right, read_micro_homology_right, ref_micro_homology_left = report_microhomologies(item, ref_seq, breakpoint_from_right, breakpoint_from_left)
                read_micro_homology_total = read_micro_homology_left + read_micro_homology_right
                len_micro_total = len(read_micro_homology_total)
                end_line_list = [line_count, read_micro_homology_total, len_micro_total]
                full_line_list = item.copy()
                full_line_list.extend(end_line_list)

            if item[-2] != "deletion":
                end_line_list = [line_count, "", ""]
                full_line_list = item.copy()
                full_line_list.extend(end_line_list)

            sorting_list.append(full_line_list)
        sorting_list.sort(key=lambda x:x[-3], reverse=True)
        for sorted_item in sorting_list:
            full_line = prepare_line(sorted_item)
            micro_out.write(full_line + '\n')
        micro_out.close()

        micro_in_reader = csv.reader(open(f[:-4] + "_Count_Columns.sam", "r"), delimiter = '\t')
        complex_out = open(f[:-4] + "_Complex_Seqs.fasta", "w")
        insert_out = open(f[:-4] + "_Insert_Seqs.fasta", "w")
        fasta_out = open(f[:-4] + "_Fasta_Reads.fasta", "w")
        final_out = open(f[:-4] + "_Final.sam", "w")
        read_error_dist_out = open(f[:-4] + "_read_error_distribution.txt", "w")

        for line in micro_in_reader:
            insert_sequence = line[-6]
            recon_read = line[-4]
            read_count = str(line[-3])
            cigar_str = line[1]
            read_class = line[-5]
            name_list = [read_count, cigar_str]
            if read_class == "complex":
                write_fasta_line(name_list, insert_sequence, complex_out)
            if read_class == "insertion":
                write_fasta_line(name_list, insert_sequence, insert_out)
            write_fasta_line(name_list, recon_read, fasta_out)
            
            control_count = str(control_list.count(line[-4]))
            read_error = str((read_error_dict[recon_read]/int(read_count)) * 100)
            final_line = '\t'.join(line) + '\t' + control_count + '\t' + read_error
            final_out.write(final_line + '\n')

            read_error_dist_line = [recon_read, read_count, read_error]
            error_dist_key = 1
            for key in read_error_distribution_dict[recon_read]:
                bp_count = read_error_distribution_dict[recon_read][error_dist_key]
                read_error_dist_line.append(bp_count)
                error_dist_key += 1
            read_error_out_line = [str(s) for s in read_error_dist_line]
            read_error_out_line = '\t'.join(read_error_out_line)
            read_error_dist_out.write(read_error_out_line + '\n')

main()
