#!/usr/bin/python

# process fasta file containing all reads and extract sequences according to originating well
# output is a directory with fasta files containing T cell receptor (TCR) or B cell receptor (BCR)
# reads sorted by well
#
# process reads: 01_extract_wells.py <fasta file>
#
# the fasta file can be generated for example with PEAR

import subprocess
import re
import os
import sys
import argparse

#command line arguments
parser = argparse.ArgumentParser(description='Analyze TCR and BCR reads and extract')
parser.add_argument('inputfile', type=str, help='Filename of fasta input file')
parser.add_argument('-output', type=str, help='output directory')
args = parser.parse_args()
sequences=args.inputfile

# create output directorires if necessary
if args.output:
	os.mkdir(args.output)
	os.mkdir(args.output+'/TCR/')
	os.mkdir(args.output+'/BCR/')

# plate, row, and column  barcodes in numerical/alphabetical order
TCR_plate_barcodes =[ 'GCAGA', 'TCGAA', 'AACAA', 'GGTGC', 'TTGGT', 'CATTC', 'ATTGG', 'CGGTT', 'ATCCT', 'ATGTC', 'TCACG', 'AGACC', 'CCCCA', 'GCGCT', 'TCCTT', 'TATAT', 'CGTAA', 'AAGGT', 'AGCTC', 'CTTGC', 'GTATC', 'TATGA', 'CACAC', 'ACACT', 'ACTAC', 'GTTAC']
TCR_row_barcodes = [ 'TAAGC', 'TGCAC', 'CTCAG', 'GGAAT', 'CGAGG', 'AGGAG', 'TGTTG', 'CAACT']
TCR_column_barcodes = ['TGAAC', 'TCCTG', 'TATAA', 'ACAGG', 'GCGGT', 'TAAGT', 'CTAGC', 'ACGTC', 'TAGCC', 'CATTC', 'GTTGG', 'GTCTC']

BCR_plate_barcodes = ['ATGATACA', 'AACATACA', 'ATGTAACA', 'ATGATTGA', 'AAGATAGA', 'ATCATTCA', 'ATCTTACA', 'ATGAATCA', 'ATCTATCA', 'TTCATACA', 'TTGATACT', 'ATGATTCT', 'TGAGCCGC', 'CCTGCCAG', 'CCACCGTC', 'CCCATACA', 'CCGTTACA', 'CCGTGGCA', 'CCGTATCA', 'GGTATACA']
BCR_column_barcodes = [ 'ATGATACA', 'AACATACA', 'ATGTAACA', 'ATGATTGA', 'AAGATAGA', 'ATCATTCA', 'ATCTTACA', 'ATGAATCA', 'ATCTATCA', 'TTCATACA', 'TTGATACT', 'ATGATTCT']
BCR_row_barcodes = [ 'TGTATCAT', 'TGTATGTT', 'TGTTACAT', 'TCAATCAT', 'TCTATCTT', 'TGAATGAT', 'TGTAAGAT', 'TGATTCAT']

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A', 'N':'N'}[B] for B in x][::-1])

t_ass_seqs={}
b_ass_seqs={}
unassigned=0
assigned=0
read_number=0
tcr_count=0
bcr_count=0

# open fasta file
seq=open(sequences, mode='rU')
DNA=set('ATGC\n')

# parse fasta file line by line
for line in seq:
    if set(line) <= DNA:
	read_number+=1
        # print out statistics of progress
	print '\rAssigned', assigned, '/', read_number, 'TCR: ', tcr_count, 'BCR: ', bcr_count,
        line=line.replace('\n','')
        # T cell receptor read found
	if line[2:7] in TCR_plate_barcodes and line[9:14] in TCR_row_barcodes and line[len(line)-7:-2] in TCR_column_barcodes and not 'N' in line and not 'AAAAAAAAAAAA' in line and not 'GGGGGGGGGGGG' in line and not 'GGGGGGGGGGGG' in line:
		location=str(TCR_plate_barcodes.index(line[2:7])+1) + str(chr(ord('A')+TCR_row_barcodes.index(line[9:14]))) + str(TCR_column_barcodes.index(line[len(line)-7:-2])+1)
		if location in t_ass_seqs:
			t_ass_seqs[location].append(revcompl(line[14:].replace('\n','')))
		else:
			t_ass_seqs[location]=[revcompl(line[14:].replace('\n', ''))]
		tcr_count+=1
		assigned+=1
        # B cell receptor read found
	elif line[4:12] in BCR_plate_barcodes and line[14:22] in BCR_column_barcodes and line[len(line)-12:-4] in BCR_row_barcodes and not 'N' in line and not 'AAAAAAAAAAAA' in line and not 'GGGGGGGGGGGG' in line and not 'GGGGGGGGGGGG' in line:
		location=str(BCR_plate_barcodes.index(line[4:12])+1) + str(chr(ord('A')+BCR_row_barcodes.index(line[len(line)-12:-4]))) + str(BCR_column_barcodes.index(line[14:22])+1)
		if location in b_ass_seqs:
			b_ass_seqs[location].append(revcompl(line[22:].replace('\n','')))
		else:
                	b_ass_seqs[location]=[revcompl(line[22:].replace('\n',''))]
		bcr_count+=1
		assigned+=1
        # sequence not assigned
        else:
                unassigned+=1

# statistics of assigned and unassigned sequences
print 'Unassigned: ', unassigned, ' Assigned: ', assigned
seq.close()

# write sequences to files
# T cell receptor reads
for key in t_ass_seqs.keys():
    counter=1
    outfile=open(args.output+'/TCR/' + key + '.fasta', mode='w')
    for line in t_ass_seqs[key]:
        outfile.write('>' + str(counter) + '\n')
        outfile.write(line + '\n')
        counter=counter+1
    outfile.close()
# B cell receptor reads
for key in b_ass_seqs.keys():
    counter=1
    outfile=open(args.output+'/BCR/' + key + '.fasta', mode='w')
    for line in b_ass_seqs[key]:
        outfile.write('>' + str(counter) + '\n')
        outfile.write(line + '\n')
        counter=counter+1
    outfile.close()

print('well files written')
