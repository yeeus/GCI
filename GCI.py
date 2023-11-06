import sys
import pysam
import numpy as np
import argparse
import subprocess
import os
import math

def read_paf(paf_files=[], map_qual=30, iden_percent=0.9):
	"""
	input: paf files,
		   the filtered mapping quality,
		   the filtered identity fraction

	usage: read paf file and simply filter the mapping quality and identity

	return: one dictionary with the queries as the keys and valued by a tuple containing the target, start, and end;
			one dictionary keyed by the targets with the length value 
	"""
	
	paf_lines = {}
	targets_length = {}
	if len(paf_files) == 1:
		with open(paf_files[0], 'r') as f:
			for line in f:
				paf = line.strip().split("\t")
				query = paf[0]
				target = paf[5]
				target_length = int(paf[6])
				start = int(paf[7])
				end = int(paf[8])
				num_match_res = int(paf[9])
				len_aln = int(paf[10])
				mapq = int(paf[11])

				if (mapq >= map_qual) and (num_match_res/len_aln >= iden_percent):
					paf_lines[query] = (target, start, end)

				if target not in targets_length.keys():
					targets_length[target] = target_length

	return paf_lines, targets_length


def filter_bam(paf_lines={}, targets_length={}, bam_files=[], prefix='GCI', map_qual=30, clip_percent=0.1, ovlp_percent=0.9, flank_len=10, dictionary='.', force=False):
	"""
	input: paf_lines and targets_length generated by the function read_paf(),
		   bam files,
		   the prefix of output depth.bed file, 
		   the filtered mapping quality,
		   the filtered clipped fraction (S / (M + I + S)),
		   the filtered overlapping fraction (overlap/length),
		   the length of flanking bases ([start + num, end - num + 1]),
		   the path to output,
		   whether to rewrite the existing files (force)

	usage: filter the bam (and paf) file based on many metrics, and finally generate the depth file 
	
	output: the whole-genome depth file

	return: the whole-genome depth dictionary
	"""
	if len(bam_files) == 1:
		samfile = pysam.AlignmentFile(bam_files[0], 'rb')
	depths = {target:np.zeros(length, dtype=int) for target, length in targets_length.items()}
	
	for segment in samfile.fetch():
		if segment.query_name in paf_lines:
			if (segment.is_mapped == True) and (segment.is_secondary == False) and (segment.is_supplementary == False) and (segment.mapping_quality >= map_qual) and (segment.reference_name == paf_lines[segment.query_name][0]):
				M = segment.get_cigar_stats()[0][0]
				I = segment.get_cigar_stats()[0][1]
				S = segment.get_cigar_stats()[0][4]
				if S / (M + I + S) <= clip_percent:
					paf_start = paf_lines[segment.query_name][1]
					paf_end = paf_lines[segment.query_name][2]
					sam_start = segment.reference_start
					sam_end = segment.reference_end

					ovlp = min(paf_end, sam_end) - max(paf_start, sam_start)
					if ovlp/segment.query_length >= ovlp_percent:
						start = max(paf_start, sam_start) + flank_len
						end = min(paf_end, sam_end) - flank_len
						depths[segment.reference_name][start:end+1] += 1
			del paf_lines[segment.query_name]


	if os.path.exists(f'{dictionary}/{prefix}.depth') and force == False:
		print(f'ERROR!!! The file "{dictionary}/{prefix}.depth" exists', file=sys.stderr)
		raise SystemExit
	else:
		with open(f'{dictionary}/{prefix}.depth', 'w') as f:
			for target, depth_list in depths.items():
				for i, depth in enumerate(depth_list):
					f.write(f'{target}\t{i}\t{depth}\n')

	return depths


def merge_depth(depths={}, prefix='GCI', threshold=0, flank_len=10, dictionary='.', force=False):
	"""
	input: the whole-genome depth dictionary generated by filter_bam(),
		   the prefix of output threshold.depth.bed file,
		   the threshold of depth,
		   the length of flanking bases,
		   the path to output,
		   whether to rewrite the existing files (force)
	
	usage: merge positions with depth lower than the threshold
	
	output: the merged depth file in bed format

	return: the dictionary containing the merged depth bed file
	"""

	merged_depths_bed = {target:[] for target in depths.keys()}
	if os.path.exists(f'{dictionary}/{prefix}.{threshold}.depth.bed') and force == False:
		print(f'ERROR!!! The file "{dictionary}/{prefix}.{threshold}.depth.bed" exists', file=sys.stderr)
		raise SystemExit
	else:
		with open(f'{dictionary}/{prefix}.{threshold}.depth.bed', 'w') as f:
			for target, depth_list in depths.items():
				start_flag = 0
				end_flag = 1
				chr_len = len(depth_list)
				for i, depth in enumerate(depth_list[flank_len:chr_len-flank_len]):
					if depth <= threshold:
						if start_flag == 0:
							start = i + flank_len
							start_flag = 1
							end_flag = 0
						if i == (chr_len - flank_len*2 - 1): 
							end = i + flank_len + 1
							f.write(f'{target}\t{start}\t{end}\n')
							merged_depths_bed[target].append((start, end))
					else:
						if end_flag == 0:
							end = i + flank_len
							f.write(f'{target}\t{start}\t{end}\n')
							merged_depths_bed[target].append((start, end))
							end_flag = 1
							start_flag = 0

	return merged_depths_bed

def compute_index(targets_length={}, prefix='GCI', dictionary='.', force=False, merged_depths_bed={}, flank_len=10):
	###########################
	### version = 1.0
	### compute the radio of N50 generated by the candidate gaps
	### and that of the original genome
	###########################
	"""
	input: targets_length generated by the function read_paf(),
		   the prefix of the output gci file,
		   the path to output,
		   whether to rewrite the existing files (force),
		   merged_depths_bed generated by merge_depth(),
		   the length of flanking bases

	usage: remove the regions with depth lower than the threshold and compute the index
		   
	output: an index file containing the original N50, new N50, genome continuity index
	"""

	origin_lengths = sorted([length for length in targets_length.values()], reverse=True)
	origin_cum = np.cumsum(origin_lengths)
	for i, number in enumerate(origin_cum):
		if number >= origin_cum[-1] / 2:
			origin_n50 = origin_lengths[i]
			break
	

	new_lengths = []
	for target, length in targets_length.items():
		last = flank_len
		n = len(merged_depths_bed[target])
		for i, segment in enumerate(merged_depths_bed[target]):
			if i != n-1:
				new_lengths.append(segment[0] - last)
				last = segment[1]
			else:
				new_lengths.append(segment[0] - last)
				new_lengths.append(length - flank_len - segment[1])
	
	new_lengths = sorted(new_lengths, reverse=True)
	new_cum = np.cumsum(new_lengths)
	for i, number in enumerate(new_cum):
		if number >= new_cum[-1] / 2:
			new_n50 = new_lengths[i]
			break
	
	if os.path.exists(f'{dictionary}/{prefix}.gci') and force == False:
		print(f'ERROR!!! The file "{dictionary}/{prefix}.gci" exists', file=sys.stderr)
		raise SystemExit
	else:
		with open(f'{dictionary}/{prefix}.gci', 'w') as f:
			f.write('Original N50\tNew N50\tGenome Continuity Index\n')
			f.write(f'{origin_n50}\t{new_n50}\t{100 * math.log2(new_n50/origin_n50 + 1)}\n')


def GCI(hifi=[], nano=[], dictionary='.', prefix='GCI', threads=1, map_qual=30, iden_percent=0.9, ovlp_percent=0.9, clip_percent=0.1, flank_len=10, threshold=0, plot=False, force=False):
	if dictionary.endswith('/'):
		dictionary = dictionary.split('/')[0]

	if os.path.exists(dictionary):
		if not os.access(dictionary, os.R_OK):
			print(f'ERROR!!! The path "{dictionary}" is unable to read', file=sys.stderr)
			raise SystemExit
		if not os.access(dictionary, os.W_OK):
			print(f'ERROR!!! The path "{dictionary}" is unable to write', file=sys.stderr)
			raise SystemExit
	else:
		os.makedirs(dictionary)


	if prefix.endswith('/'):
		print(f'ERROR!!! The prefix "{prefix}" is not allowed', file=sys.stderr)
		raise SystemExit
	

	if len(hifi) > 0:
		hifi_paf = []
		hifi_bam = []
		for file in hifi:
			CompletedProcess = subprocess.run(f'file -r {file}', shell=True, check=True, capture_output=True)
			if 'gzip' in str(CompletedProcess.stdout):
				hifi_bam.append(file)
			else:
				hifi_paf.append(file)
		
		paf_lines, targets_length = read_paf(hifi_paf, map_qual, iden_percent)
		depths = filter_bam(paf_lines, targets_length, hifi_bam, prefix, map_qual, clip_percent, ovlp_percent, flank_len, dictionary, force)
		merged_depth_bed = merge_depth(depths, prefix, threshold, flank_len, dictionary, force)
		compute_index(targets_length, prefix, dictionary, force, merged_depth_bed, flank_len)
	

	if len(nano) > 0:
		nano_paf = []
		nano_bam = []
		for file in nano:
			CompletedProcess = subprocess.run(f'file -r {file}', shell=True, check=True, capture_output=True)
			if 'gzip' in str(CompletedProcess.stdout):
				nano_bam.append(file)
			else:
				nano_paf.append(file)
		
		paf_lines, targets_length = read_paf(nano_paf, map_qual, iden_percent)
		depths = filter_bam(paf_lines, targets_length, nano_bam, prefix, map_qual, clip_percent, ovlp_percent, flank_len, dictionary, force)
		merged_depth_bed = merge_depth(depths, prefix, threshold, flank_len, dictionary, force)
		compute_index(targets_length, prefix, dictionary, force, merged_depth_bed, flank_len)


if __name__=='__main__':
	######################
	### version = 1.0
	### ideally at least one bam file for one type of reads (HiFi and/or ONT) alignment
	### but this initial version can only analyze the case of inputting one paf and one bam
	### and can only parse one type of reads alignment 
	######################
	version = '1.0'

	parser = argparse.ArgumentParser(prog=sys.argv[0], add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description='A program for assessing the T2T genome', epilog='Examples:\npython GCI.py --hifi hifi.bam hifi.paf --nano ont.bam ont.paf')

	group_io = parser.add_argument_group("Input/Output")
	group_io.add_argument('--hifi', nargs='+', metavar='', help='PacBio HiFi reads alignment files (.bam and/or .paf) including at least one bam file')
	group_io.add_argument('--nano', nargs='+', metavar='', help='Oxford Nanopore long reads alignment files (.bam and/or .paf) including at least one bam file')
	group_io.add_argument('-d', nargs='?', dest='dictionary', metavar='PATH', help='The dictionary of output files [.]', const='.', default='.')
	group_io.add_argument('-o', '--output', nargs='?', dest='prefix', metavar='STR', help='Prefix of output files [GCI]', const='GCI', default='GCI')
	group_io.add_argument('-t', '--threads', nargs='?', metavar='INT', type=int, help='Number of threads [1]', const=1, default=1)

	group_fo = parser.add_argument_group("Filter Options")
	group_fo.add_argument('-mq', '--map-qual', nargs='?', metavar='INT', type=int, help='Minium mapping quality for alignments in both bam and paf files [30]', const=30, default=30)
	group_fo.add_argument('-ip', '--iden-percent', nargs='?', metavar='FLOAT', type=float, help='Minimum identity (num_match_res/len_aln) of the reads in paf files [0.9]', const=0.9, default=0.9)
	group_fo.add_argument('-op', '--ovlp-percent', nargs='?', metavar='FLOAT', type=float, help='Minimum overlapping percentage of the reads if inputting more than one file [0.9]', const=0.9, default=0.9)
	group_fo.add_argument('-cp', '--clip-percent', nargs='?', metavar='FLOAT', type=float, help='Maximum clipped percentage of the reads in bam files [0.1]', const=0.1, default=0.1)
	group_fo.add_argument('-fl', '--flank-len', nargs='?', metavar='INT', type=int, help='The flanking length of the clipped bases [10]', const=10, default=10)
	group_fo.add_argument('-ts', '--threshold', nargs='?', metavar='INT', type=int, help='The threshold of depth in the final bed file [0]', const=0, default=0)

	group_op = parser.add_argument_group("Other Options")
	group_op.add_argument('-f', '--force', action='store_const', help='Force rewriting of existing files', const=True, default=False)
	group_op.add_argument('-p', '--plot', action='store_const', help='Visualize the final result', const=True, default=False)
	group_op.add_argument('-h', '--help', action="help", help="Show this help message and exit")
	group_op.add_argument('-v', '--version', action="version", version=version, help="Show program's version number and exit")
	
	args = vars(parser.parse_args())
	#print(args)


	if (args['hifi'] == None) and (args['nano'] == None):
		print('ERROR!!! Please input at least one type of TGS reads alignment files (PacBio HiFi and/or Oxford Nanopore long reads)\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
		raise SystemExit
	
	if args['hifi'] != None:
		bam_num = 0
		for file in args['hifi']:
			if os.path.exists(file) and os.access(file, os.R_OK):
				CompletedProcess = subprocess.run(f'file -r {file}', shell=True, check=True, capture_output=True)
				if 'gzip' in str(CompletedProcess.stdout):
					bam_num += 1
			else:
				print(f'ERROR!!! "{file}" is not an available file', file=sys.stderr)
				raise SystemExit
		if bam_num == 0:
			print('ERROR!!! Please input at least one PacBio HiFi reads bam file\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit
		
	if args['nano'] != None:
		bam_num = 0
		for file in args['nano']:
			if os.path.exists(file) and os.access(file, os.R_OK):
				CompletedProcess = subprocess.run(f'file -r {file}', shell=True, check=True, capture_output=True)
				if 'gzip' in str(CompletedProcess.stdout):
					bam_num += 1
			else:
				print(f'ERROR!!! "{file}" is not an available file', file=sys.stderr)
				raise SystemExit
		if bam_num == 0:
			print('ERROR!!! Please input at least one Oxford Nanopore long reads bam file\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit
	
	GCI(**args)

### benchmark
### clip_percent: 0.05, 0.1, 0.2
### flank_len: 10, 20, 50
### performance: distribution, `index`(candidate N50 / original N50), visualization: karyotypeR
