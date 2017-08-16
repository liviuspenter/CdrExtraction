#!/usr/bin/python

# process output from IMGT

import argparse

def Process_IMGT_Junction_Output(filename):
	a1_info = {}
	a2_info = {}
	b_info = {}

	# open IMGT Junction file
	f = open(filename, 'rU')
	# collect headline information
	header = f.readline().strip().split('\t')

	# go through file and process each line
	for line in f:
		splitted_line = line.strip().split('\t')
		# discard lines without results or unknown CDR3
		if splitted_line[header.index('V-DOMAIN Functionality')] != 'No results' and splitted_line[header.index('V-DOMAIN Functionality')] != 'unknown' and len(splitted_line) > 5:
			# decode information contained in Sequence ID
			well, seq_id, reads = splitted_line[header.index('Sequence ID')].split(':')
			# gather other information
			CDR3 = splitted_line[header.index('CDR3-IMGT (AA)')]
			CDR3_nt = splitted_line[header.index('CDR3-IMGT')]
			if (splitted_line[header.index('JUNCTION frame')] == 'in-frame'):
				CDR3_frameshift = 'no'
			else:
				CDR3_frameshift = 'yes'
			V_Region = splitted_line[header.index('V-GENE and allele')].split(' ')[1]
			J_Region = splitted_line[header.index('J-GENE and allele')].split(' ')[1]	
			# assemble information to store
			assembled_information = [reads, CDR3, CDR3_nt, CDR3_frameshift, V_Region, J_Region, splitted_line[header.index('Sequence ID')]]
			# store information obtained
			if V_Region[2] == 'A':
				if well not in a1_info.keys():
					a1_info[well] = assembled_information
				elif well in a1_info.keys() and well not in a2_info.keys():
					a2_info[well] = assembled_information
			if V_Region[2] == 'B' and well not in b_info.keys():
					b_info[well] = assembled_information

	return a1_info, a2_info, b_info
	
def Process_IMGT_Summary_Output(filename):
	c_sequences = {}
	imgt_summary = open(args.imgt_output_summary, 'rU')
	imgt_summary_header = imgt_summary.readline().strip().split('\t')
	for line in imgt_summary:
		if len(line.strip().split('\t')) > 5:
			c_sequences[line.strip().split('\t')[imgt_summary_header.index('Sequence ID')]] = line.strip().split('\t')[imgt_summary_header.index('Sequence')]
	imgt_summary.close()

	return c_sequences

def Read_Cytokine_Data(filename):
	c_info = {}
	c_file = open(filename, 'rU')
	cytokine_header = c_file.readline().strip().split('\t')
	for line in c_file:
		c_info[line.strip().split('\t')[0]] = line.strip().split('\t')[1:]
	c_file.close()

	return c_info, cytokine_header

# parsing arguments
parser = argparse.ArgumentParser(description='Process output from IMGT and merge with cytokine data.')
parser.add_argument('--imgt_output_summary', required=True, help='1_Summary.txt from IMGT High/V-Quest')
parser.add_argument('--imgt_output_junction', required=True, help='6_Junction.txt from IMGT High/V-Quest')
parser.add_argument('--cytokine_output', required=True, help='File that contains output from cytokine analysis')
parser.add_argument('--output', required=True, help='File that will contain single cell sequencing data plus cytokine data. ')
args = parser.parse_args()

# read IMGT file 6_Junction.txt
a1_info, a2_info, b_info = Process_IMGT_Junction_Output(args.imgt_output_junction)

# read IMGT file 1_Summary.txt to obtain consensus sequences
c_sequences = Process_IMGT_Summary_Output(args.imgt_output_summary)

# read cytokine list
c_info, cytokine_header = Read_Cytokine_Data(args.cytokine_output)

# write file containing single cell sequencing data plus cytokine data
out_file = open(args.output, 'w')
out_file.write ('Well\tReads alpha\tCDR3 alpha\tCDR3 alpha nt\tCDR3 alpha frameshift\tV-Region alpha\tJ-Region alpha\tReads alpha2\tCDR3 alpha2\t')
out_file.write ('CDR3 alpha2 nt\tCDR3 alpha2 frameshift\tV-Region alpha2\tJ-Region alpha2\tReads beta\tCDR3 beta\tCDR3 beta nt\tCDR3 beta frameshift\tV-Region beta\tJ-Region beta\t')
out_file.write ('\t'.join(cytokine_header[1:]))
out_file.write ('\tConsensus Alpha1\tConsensus Alpha2\tConsensus Beta\n')

# go through all wells and print out information
for well in sorted(c_info.keys()):
	consensus_a1 = 'XXXXXXXXXXXX'
	consensus_a2 = 'XXXXXXXXXXXX'
	consensus_b = 'XXXXXXXXXXXX'
	out_file.write(well + '\t')
	if well in a1_info.keys():
		out_file.write('\t'.join(a1_info[well][0:6]) + '\t')
		consensus_a1 = c_sequences[a1_info[well][6]]
	else:
		out_file.write('0\t' + '\t'.join(5*['XXXXXXXXXXXX']) + '\t')

	if well in a2_info.keys():
		out_file.write('\t'.join(a2_info[well][0:6]) + '\t')
		consensus_a2 = c_sequences[a2_info[well][6]]
	else:
		out_file.write('0\t' + '\t'.join(5*['XXXXXXXXXXXX']) + '\t')

	if well in b_info.keys():
		out_file.write('\t'.join(b_info[well][0:6]) + '\t')
		consensus_b = c_sequences[b_info[well][6]]
	else:
		out_file.write('0\t' + '\t'.join(5*['XXXXXXXXXXXX']) + '\t')
	
	out_file.write('\t'.join(c_info[well]) + '\t' + consensus_a1 + '\t' + consensus_a2 + '\t' + consensus_b + '\n')

out_file.close()
