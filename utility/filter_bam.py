import sys
import pysam
import numpy as np
import argparse
import subprocess
import os
import random
import string


def get_average_identity(alns):
	"""
	usage: get average identity if there are many alignment blocks of one target 

	input: synteny[query][target] from function filter()

	return: the average identity
	"""

	tmp = []
	for a in alns:
		tmp.append(a[-1])
	average = sum(tmp) / len(alns)
	return average


def merge_alns_properties(alns, x, y):
	"""
	usage: merge overlapped aligned blocks of either query or target based on the inputted parameters x and y

	input: ynteny[query][target] from function filter(),
		   query (1, 2) or target (3, 4)

	return: mapped_length of the query (inputting 1, 2), the lowest and highest position of target (inputting 3, 4)
	"""

	# extract bed list form alns
	bed_list = []
	for a in alns:
		bed_list.append([a[x], a[y]])

	# merge bed to generate non-overlap blocks
	mapped_length = 0
	sort_bed_list = sorted(bed_list)

	low_est = sort_bed_list[0][0]
	high_est = sort_bed_list[0][1]

	for index, block in enumerate(sort_bed_list):
		low, high = block
		if high_est >= low:
			if high_est < high:
				high_est = high
		else:
			mapped_length += (high_est - low_est)
			low_est, high_est = sort_bed_list[index]
	mapped_length += (high_est - low_est)
	return mapped_length, sort_bed_list[0][0], high_est


def filter(paf_files=[], bam_files=[], prefix='GCI', map_qual=30, iden_percent=0.9, clip_percent=0.1, ovlp_percent=0.9, directory='.', force=False):
	"""
	usage: filter the paf and bam file(s) based on many metrics and output the filtered bam file(s)

	input: paf file(s),
		   bam file(s),
		   the filtered mapping quality,
		   the filtered identity percentage,
		   the filtered clipped percentage (S / (M + I + S)),
		   the filtered overlapping percentage (overlap/length),
		   the path to output,
		   whether to rewrite the existing files (force)
	
	output: the filtered bam file(s)

	return: the whole-genome depth dictionary,
			one dictionary keyed by the targets with the length value
	"""
	
	output_files = []
	#output_files_name = []
	samfile = pysam.AlignmentFile(bam_files[0], 'rb')
	for i, file in enumerate(bam_files):
		if isinstance(prefix, str):
			if os.path.exists(f'{directory}/{".".join(file.split(".")[:-1])}.filter.bam') and force == False:
				print(f'ERROR!!! The file "{directory}/{".".join(file.split(".")[:-1])}.filter.bam" exists\nPlease using "-f" or "--force" to rewrite', file=sys.stderr)
				raise SystemExit
			#output_files_name.append(f'{directory}/{".".join(file.split(".")[:-1])}.filter.bam')
			f = pysam.AlignmentFile(f'{directory}/{".".join(file.split(".")[:-1])}.filter.bam', 'wb', template=samfile)
		else:
			if os.path.exists(f'{directory}/{prefix[i]}.bam') and force == False:
				print(f'ERROR!!! The file "{directory}/{prefix[i]}.bam" exists\nPlease using "-f" or "--force" to rewrite', file=sys.stderr)
				raise SystemExit
			#output_files_name.append(f'{directory}/{prefix[i]}.bam')
			f = pysam.AlignmentFile(f'{directory}/{prefix[i]}.bam', 'wb', template=samfile)
		output_files.append(f)
	samfile.close()


	paf_lines = [{} for i in range(len(paf_files))]
	if len(paf_files) != 0:
		synteny = {}
		for i, file in enumerate(paf_files):
			with open(file, 'r') as f:
				for line in f:
					paf = line.strip().split("\t")
					query = paf[0]
					query_length = int(paf[1])
					query_start = int(paf[2])
					query_end = int(paf[3])
					target = paf[5]
					target_start = int(paf[7])
					target_end = int(paf[8])
					num_match_res = int(paf[9])
					len_aln = int(paf[10])
					mapq = int(paf[11])
					
					identity = num_match_res/len_aln
					if (mapq >= map_qual) and (identity >= iden_percent):
						if query not in synteny.keys():
							synteny[query] = {}
						if target not in synteny[query].keys():
							synteny[query][target] = []
						synteny[query][target].append((query_length, query_start, query_end, target_start, target_end, identity))
				for query in synteny.keys():
					mapping_results = {}
					for target in synteny[query].keys():
						alns = synteny[query][target]
						non_overlap_qry_aligned, _, _ = merge_alns_properties(alns, 1, 2)
						query_length = alns[0][0]
						alignrate = non_overlap_qry_aligned / query_length
						#? if alignrate >= xxx 
						average_identity = get_average_identity(alns)
						score = average_identity * alignrate
						_, start, end = merge_alns_properties(alns, 3, 4)
						mapping_results[target] = (score, start, end, query_length)
					primary_target = sorted(mapping_results, key=lambda k: (mapping_results[k][0], k), reverse=True)[0]
					primary_target_result = mapping_results[primary_target]
					paf_lines[i][query] = (primary_target, primary_target_result[1], primary_target_result[2], primary_target_result[-1])


	samfile_dicts = [{} for i in range(len(bam_files))]
	for i, file in enumerate(bam_files):
		samfile = pysam.AlignmentFile(file, 'rb')
		for segment in samfile.fetch():
			if (segment.is_mapped == True) and (segment.is_secondary == False) and (segment.is_supplementary == False) and (segment.mapping_quality >= map_qual):
				M = segment.get_cigar_stats()[0][0]
				I = segment.get_cigar_stats()[0][1]
				D = segment.get_cigar_stats()[0][2]
				S = segment.get_cigar_stats()[0][4]
				NM = segment.get_tag('NM')
				mm = NM - (I + D)
				if (S/(M+I+S) <= clip_percent) and ((M-mm)/(M+I+D) >= iden_percent):
					samfile_dicts[i][segment.query_name] = (segment.reference_name, segment.reference_start, segment.reference_end, segment.query_length)
		samfile.close()


	files = paf_lines + samfile_dicts
	if len(files) > 1:
		files_sets = []
		for file in files:
			files_sets.append(set(file.keys()))
		querys = set.intersection(*files_sets)

		file1 = {query:segment for query, segment in files[0].items() if query in querys}
		for file in files[1:]:
			for query, segment in file.items():
				if query in file1.keys():
					segment1 = file1[query]
					if segment[0] ==  segment1[0]:
						start1 = segment[1]
						end1 = segment[2]
						start2 = segment1[1]
						end2 = segment1[2]

						ovlp = min(end1, end2) - max(start1, start2)
						if ovlp/segment[-1] < ovlp_percent:
							del file1[query]
						else:
							file1[query] = (segment1[0], max(start1, start2), min(end1, end2))
	else:
		file1 = files[0]
	for i, bamfile in enumerate(bam_files):
		samfile = pysam.AlignmentFile(bamfile, 'rb')
		for query, segment in file1.items():
			for segment1 in samfile.fetch(segment[0], segment[1], segment[2]):
				if (segment1.query_name == query) and (segment1.is_mapped == True) and (segment1.is_secondary == False) and (segment1.is_supplementary == False) and (segment1.mapping_quality >= map_qual):
					output_files[i].write(segment1)
					break
		samfile.close()

	"""for file in output_files_name:
		random_string = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
		r1 = subprocess.run(f'samtools sort {file} -o {directory}/{random_string}.bam', shell=True, capture_output=True, check=True)
		print(r1.stdout)
		print(r1.stderr)
		subprocess.run(f'mv {directory}/{random_string}.bam {file}', shell=True, capture_output=True)"""


def preprocessing(files=[], directory='.', prefix='GCI', threads=1, map_qual=30, iden_percent=0.9, ovlp_percent=0.9, clip_percent=0.1, plot=False, region='', regions_file='', force=False):
	if directory.endswith('/'):
		directory = directory.split('/')[0]
	if os.path.exists(directory):
		if not os.access(directory, os.R_OK):
			print(f'ERROR!!! The path "{directory}" is unable to read', file=sys.stderr)
			raise SystemExit
		if not os.access(directory, os.W_OK):
			print(f'ERROR!!! The path "{directory}" is unable to write', file=sys.stderr)
			raise SystemExit
	else:
		os.makedirs(directory)


	bam_files = []
	paf_files = []
	for file in files:
		CompletedProcess = subprocess.run(f'file -r {file}', shell=True, check=True, capture_output=True)
		if 'gzip' in str(CompletedProcess.stdout):
			bam_files.append(file)
		else:
			paf_files.append(file)
	
	filter(paf_files, bam_files, prefix, map_qual, iden_percent, clip_percent, ovlp_percent, directory, force)


if __name__=='__main__':
	###! this version can't plot and use multiple threads (-p, -t)

	parser = argparse.ArgumentParser(prog=sys.argv[0], add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description='This is the core part of the main program GCI.py with additional functions\nYou can provide some bam and paf files (like GCI.py) and set the filtered parameters\nand then you will get the filtered bam file which can be used for subsequent analysis', epilog='Examples:\npython filter_bam.py bam1 paf1 ...')

	group_io = parser.add_argument_group("Input/Output")
	group_io.add_argument('files', nargs='+', metavar='ALIGNMENT-FILE', help='Long reads alignment files (at least one bam file)')
	group_io.add_argument('-d', dest='directory', metavar='PATH', help='The directory of output files [.]', default='.')
	group_io.add_argument('-o', '--output', nargs='*', dest='prefix', metavar='STR', help='Prefix of output files; one prefix corresponds to one bam file in order and if provide the parameter "-p", the last one was used as the prefix for bamsnap outputs [[$input.filter] for filtered bam files and [GCI] for bamsnap outputs]', default='GCI')
	group_io.add_argument('-t', '--threads', metavar='INT', type=int, help='Number of threads [1]', default=1)

	group_fo = parser.add_argument_group("Filter Options")
	group_fo.add_argument('-mq', '--map-qual', metavar='INT', type=int, help='Minium mapping quality for alignments [30]', default=30)
	group_fo.add_argument('-ip', '--iden-percent', metavar='FLOAT', type=float, help='Minimum identity (num_match_res/len_aln) of the reads [0.9]', default=0.9)
	group_fo.add_argument('-op', '--ovlp-percent', metavar='FLOAT', type=float, help='Minimum overlapping percentage of the reads if inputting more than one alignment files [0.9]', default=0.9)
	group_fo.add_argument('-cp', '--clip-percent', metavar='FLOAT', type=float, help='Maximum clipped percentage of the reads [0.1]', default=0.1)

	group_po = parser.add_argument_group("Plot Options")
	group_po.add_argument('-p', '--plot', action='store_const', help='Visualize the filtered bam files', const=True, default=False)
	group_po.add_argument('-r', '--region', metavar='STR', help='The region to plot in chr:pos or chr:start-end format')
	group_po.add_argument('-R', '--regions-file', metavar='FILE', help='Bed file contains the regions to plot')

	group_op = parser.add_argument_group("Other Options")
	group_op.add_argument('-f', '--force', action='store_const', help='Force rewriting of existing files', const=True, default=False)
	group_op.add_argument('-h', '--help', action="help", help="Show this help message and exit")
	
	args = vars(parser.parse_args())
	print(f'Used arguments:{args}')
	
	bam_num = 0
	for file in args['files']:
		if os.path.exists(file) and os.access(file, os.R_OK):
			CompletedProcess = subprocess.run(f'file -r {file}', shell=True, check=True, capture_output=True)
			if 'gzip' in str(CompletedProcess.stdout):
				bam_num += 1
		else:
			print(f'ERROR!!! "{file}" is not an available file', file=sys.stderr)
			raise SystemExit
	if bam_num == 0:
			print('ERROR!!! Please input at least one bam file\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit
	
	if isinstance(args['prefix'], str):
		pass
	elif len(args['prefix']) ==0:
		print('ERROR!!! Please input at least one prefix\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
		raise SystemExit
	elif (args['plot'] == False):
		if len(args['prefix']) != bam_num:
			print(f'ERROR!!! The number of prefixes and bam files is inconsistent\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit
	elif (args['plot'] == True):
		if len(args['prefix']) != (bam_num + 1):
			print(f'ERROR!!! Expect {bam_num + 1} prefixes but only provide {len(args["prefix"])}\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit

	preprocessing(**args)