#!/usr/bin/env python

import gzip.open 
import arparse
import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Program to calculate average qscore values at each position for every record in a FASTQ file')
parser.add_argument('-i', '--input', type=str, help='Specify input file, should be FASTQ')
parser.add_argument('-o', '--output', type=str, help='Specify output file name. Will be in .png format')
args = parser.parse_args()


num_reads = 363246735
qscores_array =  []
for x in range(101):
    qscores_array.append(0)


ln = 0
with open(input_file, 'tr') as fh:
    for line in fh:
        if ln % 4 == 3:
            line = line.strip()
            #quality line
            nucleotide_pos = 0
            for letter in line:
                #convert qscore ascii to int
                phred_score = ord(letter)-33
                #which record we're on
                #record = ((ln+1) // 4) -1
                #add phredscore to 2d array
                qscores_array[nucleotide_pos] = phred_score
                nucleotide_pos += 1
        ln += 1


#calculating mean
means = []
for item in qscores_array:
    means.append(item / num_reads)


#plotting
x = [x for x in range(101)]
plt.xlabel('Nucleotide position')
plt.ylabel('Qscore mean')
plt.title('Averge qscores at each position in reads')
plt.bar(x, means)
plt.savefig(args.output + '.png')

