#!/usr/ bin/env python

import gzip 
import argparse
import matplotlib.pyplot as plt


def arguments():
    '''captures runtime arguments'''
    parser = argparse.ArgumentParser(description='Program to calculate average qscore values at each position for every record in a FASTQ file')
    parser.add_argument('-i', '--input', type=str, help='Specify input file, should be FASTQ')
    parser.add_argument('-o', '--output', type=str, help='Specify output file name. Will be in .png format')
    return parser.parse_args()


def initialize_qscore_list(length):
    '''
    Arguments:
        length of qscore line
    return:
        empty list of provided length'''
    qscores_list =  []
    for x in range(length):
        qscores_list.append(0)
    return qscores_list


def read_file(inp):
    '''
    Arguments:
        input file (FASTQ)
    return:
        sum of qscore values at each nucleotide position in record
        length of records
        number of records'''
    ln = 0
    with gzip.open(inp, 'tr') as fh:
        for line in fh:
            if ln == 3:
                #initializing qscore list based on length of record, assuming all lengths to be equal
                record_length = len(line.strip())
                qscores_list = initialize_qscore_list(record_length)
            if ln % 4 == 3:
                line = line.strip()
                #quality line
                nucleotide_pos = 0
                for letter in line:
                    #convert qscore ascii to int
                    phred_score = ord(letter)-33
                    #running phred score sum for each position
                    qscores_list[nucleotide_pos] += phred_score
                    nucleotide_pos += 1
            ln += 1
    #number of records in file
    recordcount = ln / 4
    return qscores_list, recordcount, record_length


def calculate_means(qscores, recordcount):
    '''caclulates the mean qscore at each position
    arguments:
        list of sums of qscores for each position
        number of records
    return:
        list of mean qscores for each position'''
    means = []
    for item in qscores:
        means.append(item / recordcount)
    return means


def plotting(means, out, record_length):
    '''Outputs the mean qscores for each position for the entire file in a bar graph
    arguments:
        list of mean qscores
        length of mean qscore list
        output file name'''
    x = [x for x in range(int(record_length))]
    plt.xlabel('Nucleotide position')
    plt.ylabel('Qscore mean')
    plt.title('Averge qscores at each position in reads')
    plt.bar(x, means)
    plt.savefig(out + '.png')

#store arguments
args = arguments()
inp = args.input
out = args.output
#gather qscores from file and get record length and count
qscores, recordcount, record_length = read_file(inp)
#calculate mean qscore at each nucleotide position
means = calculate_means(qscores, recordcount)
#plot means
plotting(means, out, record_length)
