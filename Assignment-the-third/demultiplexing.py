#!usr/bin/env python

import gzip
import itertools
import argparse
import os
import re

def get_args():
    '''collects user inputs during runtime'''
    parser = argparse.ArgumentParser(description='\
    This program demultiplexes sequencing data generated with parallel runs. Indices for records\
    are checked for congruence. Records with mismatched indices are set aside. Records with low quality\
    indices (based on qscores) or unknown indices (with N base calls) are set aside in the "unknown" output files.\
    Records with indices not in the expected list of indices are set in the unknown outputs. Must specify which\
    read each FASTQ file is. Read quality scores are assumed to be Ascii base 33. Assumed to be paired end reads with dual indexing')
    parser.add_argument('-r', '--read_one', type=str, help='file containing first read. FASTQ gzipped format.')
    parser.add_argument('-R', '--read_two', type=str, help='file containing second read. FASTQ gzipped format.')
    parser.add_argument('-i', '--index_one', type=str, help='file contaiing first index. FASTQ gzipped format.')
    parser.add_argument('-I', '--index_two', type=str, help='file containing second index. FASTQ gzipped format.')
    parser.add_argument('-o', '--output_directory', type=str, help='Specify output directory name')
    parser.add_argument('-q', '--qscore_cutoff', type=int, help='Define qscore cutoff. Indices with a mean qscore below the threshold. Defaults to 26',default=26)
    parser.add_argument('-x', '--indices', type=str, help='file with all expected indices. 4 column TSV format with indices in 4th column')
    parser.add_argument('-f', '--output_format', type=str, help='Filename base for demultiplexed FASTQ files. Each file will have this information + index information')
    args = parser.parse_args()
    return args.output_directory, args.read_one, args.read_two, args.index_one, args.index_two, args.qscore_cutoff, args.indices, args.output_format
output_dir, read_one_file, read_two_file, index_one_file, index_two_file, qscore_cutoff, indices_file, output_format = get_args()


def open_inputs(read_one_file, read_two_file, index_one_file, index_two_file):
    '''
    input: filepaths for 4 zipped FASTQ input files 
    output: open file pointers dictionary
            input file name format, used later for output format'''
    files = {}
    locations = [read_one_file, read_two_file, index_one_file, index_two_file]
    file_names = ['read1', 'read2', 'index1', 'index2']
    count = 0
    for file in locations:
        #store open file objects in dictionary with short file names as keys
        files[file_names[count]] = gzip.open(file, 'rt')
        count += 1
    #grab the first part of an input file name, to reference later in output file format
    return files
input_files = open_inputs(read_one_file, read_two_file, index_one_file, index_two_file)


def get_indices(indices):
    '''
    input: indices file
    output: list of indices
            dictionary of indices with their substitutions'''
    indexes = []
    index_subs = {}
    with open(indices, 'r') as fh:
        for line in fh:
            #removing non-indices
            if ' ' not in line:
                #remove newline
                line = line.strip()
                #grab the 4th column of index file
                index = line.split('\t')[4]
                indexes.append(index)
                #store the index info for concise output filenaming
                sub = line.split('\t')[:3]
                index_subs[index] = '_'.join(sub)
    return indexes, index_subs
indices, index_subs = get_indices(indices_file)


def index_combiner(indices):
    '''
    Input: list of indeces
    Output: dictionary with all possible index combinations as the keys and a running counter as the value'''
    #gives tuples of combos
    combos = itertools.product(indices, indices)
    formated_combos = []
    for index1, index2 in combos:
        #going from tuple to single str
        formated_combos.append(str(index1) + '-' + str(index2))
    return {index:0 for index in combos}
index_combos = index_combiner(indices)


def rev_comp(sequence):
    '''
    input: nucleotide sequence
    output: the reverser complement of that seqeunce'''
    #nucleotide dict
    sequence = sequence.strip()
    crickpairs = {'A':'T', 'T':'A', 'G':'C', 'C':'G',}
    reverse_comp = ''
    #reverse the input string
    for letter in sequence[::-1]:
        reverse_comp += crickpairs[letter]
    return reverse_comp


def open_outputs(file_name_format, output_dir, subs):
    '''
    input: file name format from opened fastq files
           output directory
           index id information
    output: 2 open files for each index to store paired reads with correct indices
            2 files for reads with low quality indices
            2 files for reads with unknown indices'''
    pwd = os.getcwd()
    os.chdir(output_dir)
    #a dictionary for each output file (forward and reverse) for each index pair
    output_files_for = {}
    output_files_rev = {}
    #make a list of paired indices
    index_subs = []
    for index, sub in subs.items():
        index_subs.append(sub)
    #open 2 files for each pair
    for index in index_subs:
        #key is substituted index
        file_name = (str(index)+file_name_format+'_R1_001.fastq')
        output_files_for[index] = open(file_name, 'w') 
        #open two files, forward and reverse
        file_name = (str(index)+file_name_format+'_R2_001.fastq')
        output_files_rev[index] = open(file_name, 'w')
    #add some outputs for mismatches and low quality reads
    file_name = 'misindexed'+file_name_format+'_R1_001.fastq'
    output_files_for['misindexed'] = open(file_name, 'w')
    file_name = 'misindexed'+file_name_format+'_R2_001.fastq'
    output_files_rev['misindexed'] = open(file_name, 'w')
    file_name = 'undetermined'+file_name_format+'_R1_001.fastq'
    output_files_for['undetermined'] = open(file_name, 'w')
    file_name = 'undetermined'+file_name_format+'_R2_001.fastq'
    output_files_rev['undetermined'] = open(file_name, 'w')
    os.chdir(pwd)
    return output_files_for, output_files_rev
output_files_for, output_files_rev = open_outputs(output_format, output_dir, index_subs)


def evaluate_record(record, qscore_cutoff):
    '''
    input: full record from 4 files (4x4 list)
           qscore threshold for cutoff
    output: whether the read is paired, index-hopped, or unknown'''
    read1 = record[0]
    index1 = record[1]
    index2 = record[2]
    read2 = record[3]
    #detecting unknown indices
    for call in (index1[1] + index2[1]):
        if call == 'N':
            return 'undetermined'
    #grab quality scores of indices
    raw_qscore_index1 = index1[3]
    raw_qscore_index2 = index2[3]
    raw_qscores = [raw_qscore_index1, raw_qscore_index2]
    #keeping converted qscores separate
    phred_scores = [0,0]
    count = 0
    #checking for low quality indices
    for raw_score in raw_qscores:
        phred_scores
        for pos in raw_score:
            #phred conversion assumes ascii base 33
            phred = ord(pos)-33
            #keep running total of phred scores for index
            phred_scores[count] += phred
        #finding average 
        avg_phred = phred_scores[count] / len(raw_qscore_index1)
        #detecting low quality indices based on average qscore
        if avg_phred < qscore_cutoff:
            return 'undetermined'
        count = 1
    #gets reverse complement of second index
    index2_rev_comp = rev_comp(index2[1])
    #check if the indices are actually in the expected indices file
    if index1[1] not in indices or index2_rev_comp not in indices:
        return 'undetermined'
    #if the qsores were acceptable, proceed to check for mismatches in indices
    for pos in range(len(index1[1])):
        #go through bp by bp, comparing to make sure they're the same
        if index1[1][pos] != index2_rev_comp[pos]:
            #the indices are not matched
            return 'misindexed'
    #if all other tests passed with no return, must be matched and high quality
    return 'paired'


#dictionary to count occurances of index pairs
index_counts = {}
index_hoping_counts = {'undetermined':0,'misindexed':0,'paired':0}


def write_output(record, evaluation, output_files_for, output_files_rev, index_subs):
    '''
    input: 4x4 list of records
        evaluation of index pairs
        files to write out to
        indices information
    output:
        write the record to the corresponding output files (2)'''
    read1 = record[0]
    index1 = record[1]
    index2 = record[2]
    read2 = record[3]
    #string of both indices for reads
    indices = index1[1] + '-' + index2[1]
    #add this to header of reads
    read1[0] += ' ' + indices
    read2[0] += ' ' + indices
    #count records per index pair
    if evaluation != 'undetermined':
        if indices in index_counts.keys():
            index_counts[indices] += 1
        else:
            index_counts[indices] = 1
    index_hoping_counts[evaluation] += 1
    count = 0
    #adding newlines back
    for line in read1:
        read1[count] = line + '\n'
        count += 1
    count = 0
    for line in read2:
        read2[count] = line + '\n'
        count += 1
    #sending to output file based on evaluation
    if evaluation == 'undetermined':
        (output_files_for['undetermined']).writelines(read1)
        (output_files_rev['undetermined']).writelines(read2)
    elif evaluation == 'misindexed':
        (output_files_for['misindexed']).writelines(read1)
        (output_files_rev['misindexed']).writelines(read2)
    elif evaluation == 'paired':
        #if the index is expected from the index log
        if index1[1] in index_subs.keys():
            index_sub = index_subs[index1[1]]
            (output_files_for[index_sub]).writelines(read1)
            (output_files_rev[index_sub]).writelines(read2)
        else:
            #this record has paired indices, but they aren't expected
            (output_files_for['undetermined']).writelines(read1)
            (output_files_rev['undetermined']).writelines(read2)

def file_reader(input_files, index_combos, qscore_cutoff, output_files_for, output_files_rev, index_subs):
    '''
    Input: 4 FASTQ files. Two indeces and two reads
    Output: a dictionary counter for the occurances of each index pair (including mismatched indeces)'''
    inputs = []
    file_names = ['read1', 'index1', 'index2', 'read2']
    for file in file_names:
        inputs.append(input_files[file])
    ln = 0
    records = []
    #4x4 list
    records = [['' for x in range(4)] for y in range(4)]
    #read through all 4 files one line at a time
    for r1, r2, r3, r4 in zip(inputs[0],inputs[1],inputs[2],inputs[3]):
        lines = [r1.strip(), r2.strip(), r3.strip(), r4.strip()]
        #store lines in 4x4 list
        for x in range(4):
            records[x][ln%4] = lines[x]
        if (ln%4 == 3):
            #reached the end of a record for all 4 files
            evaluation = evaluate_record(records, qscore_cutoff)
            write_output(records, evaluation, output_files_for, output_files_rev, index_subs)
        ln += 1
file_reader(input_files, index_combos, qscore_cutoff, output_files_for, output_files_rev, index_subs)


def file_closer(input_files, output_files_for, output_files_rev):
    '''close all open input and output files'''
    for key in input_files:
        input_files[key].close()
    for file in output_files_for.keys():
        output_files_for[file].close()
    for file in output_files_rev.keys():
        output_files_rev[file].close()
file_closer(input_files, output_files_for, output_files_rev)


def stats(output_dir):
    '''report run statistics. A full record of index pair counts is stored in a stats.txt file'''
    os.chdir(output_dir)
    with open('stats.txt', 'w') as fh:
        fh.write('Run information:\nInputs:\n')
        fh.write(str('Read 1 file:\t' + read_one_file + '\n'))
        fh.write(str('Read 2 file:\t' + read_two_file + '\n'))
        fh.write(str('Index 1 file:\t' + index_one_file + '\n'))
        fh.write(str('Index 2 file:\t' + index_two_file + '\n'))
        fh.write(str('Expected indices file' + indices_file + '\n'))
        fh.write(str('Qscore threshold:\t' + qscore_cutoff + '\n'))
        fh.write(str('Total records for index pair category:\n'))
        print('Total records for index pair category:')
        for item in index_hoping_counts.keys():
            percentage = (index_hoping_counts[item] / sum(index_hoping_counts.values())) * 100
            line = item+'\t'+str(index_hoping_counts[item])+'\t%'+str(percentage)+'\n'
            print(line, end='')
            fh.write(line)
        print()
        fh.write('\n')
        line = 'Total number of records\t' + str(sum(index_hoping_counts.values()))+'\n'
        print(line)
        fh.write(line)
        line = 'Full FASTQ output files can be found on Talapas at\n' + output_dir
        fh.write()
        print('See stats.txt in the output directory for full run stats')
        fh.write('Frequencies of each index pair:\n')
        sorted_counts = {k:v for k,v in sorted(index_counts.items(),key=lambda item: item[1], reverse=True)}
        for item in sorted_counts.keys():
            percentage = (index_counts[item] / sum(index_counts.values())) * 100
            line = item+'\t'+str(index_counts[item])+'\t%'+str(percentage)+'\n'
            fh.write(line)          
stats(output_dir)