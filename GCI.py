import sys
import pysam
import numpy as np
import argparse
import subprocess
import os
from math import log2


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


def filter(paf_files=[], bam_files=[], prefix='GCI', map_qual=30, iden_percent=0.9, clip_percent=0.1, ovlp_percent=0.9, flank_len=10, directory='.', force=False, generate=False):
	"""
	usage: filter the paf and bam file(s) based on many metrics(, and finally generate the depth file)

	input: paf file(s),
		   bam file(s),
		   the prefix of output depth file,
		   the filtered mapping quality,
		   the filtered identity percentage,
		   the filtered clipped percentage (S / (M + I + S)),
		   the filtered overlapping percentage (overlap/length),
		   the length of flanking bases ([start + num, end - num + 1]),
		   the path to output,
		   whether to rewrite the existing files (force),
		   whether to generate the depth file

	(output: the whole-genome depth file)

	return: the whole-genome depth dictionary,
			one dictionary keyed by the targets with the length value
	"""
	if generate == True:
		if os.path.exists(f'{directory}/{prefix}.depth') and force == False:
			print(f'ERROR!!! The file "{directory}/{prefix}.depth" exists\nPlease using "-f" or "--force" to rewrite', file=sys.stderr)
			raise SystemExit
	
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


	samfile = pysam.AlignmentFile(bam_files[0], 'rb')
	depths = {reference:np.zeros(length, dtype=int) for (reference, length) in zip(samfile.references, samfile.lengths)}
	targets_length = {reference:length for (reference, length) in zip(samfile.references, samfile.lengths)}
	samfile.close()

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
	if len(files) == 1:
		for segment in files[0].values():
			target = segment[0]
			start = segment[1] + flank_len
			end = segment[2] - flank_len
			depths[target][start:end+1] += 1
	else:
		querys = []
		for file in files:
			querys += list(file.keys())
		querys = set(querys)

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
		for segment in file1.values():
			target = segment[0]
			start = segment[1] + flank_len
			end = segment[2] - flank_len
			depths[target][start:end+1] += 1


	if generate == True:
		with open(f'{directory}/{prefix}.depth', 'w') as f:
			for target, depth_list in depths.items():
				f.write(f'>{target}\n')
				for i, depth in enumerate(depth_list):
					f.write(f'{i}\t{depth}\n')
	
	return depths, targets_length


def merge_two_type_depth(hifi_depths={}, nano_depths={}, prefix='GCI', directory='.', force=False, generate=False):
	"""
	usage: merge the depths dictionary generated by two types of long reads,
		   the prefix of output depth file,
		   the path to output,
		   whether to rewrite the existing files (force),
		   whether to generate the depth file

	input: the whole-genome depth dictionaries of two types of long reads generated by filter()

	(output: the whole-genome depth file)

	return: the merged whole-genome depth dictionary
	"""
	if generate == True:
		if os.path.exists(f'{directory}/{prefix}.depth') and force == False:
			print(f'ERROR!!! The file "{directory}/{prefix}.depth" exists\nPlease using "-f" or "--force" to rewrite', file=sys.stderr)
			raise SystemExit


	merged_two_type_depths = {target:[] for target in hifi_depths.keys()}
	for target, hifi_depth_list in hifi_depths.items():
		nano_depth_list = nano_depths[target]
		for i, depth in enumerate(hifi_depth_list):
			merged_two_type_depths[target].append(max(depth, nano_depth_list[i]))
	
	if generate == True:
		with open(f'{directory}/{prefix}.depth', 'w') as f:
			for target, depth_list in merged_two_type_depths.items():
				f.write(f'>{target}\n')
				for i, depth in enumerate(depth_list):
					f.write(f'{i}\t{depth}\n')

	return merged_two_type_depths


def merge_depth(depths={}, prefix='GCI', threshold=0, flank_len=10, directory='.', force=False):
	"""
	usage: merge positions with depth lower than the threshold

	input: the whole-genome depth dictionary generated by filter() and merge_two_type_depth(),
		   the prefix of output threshold.depth.bed file,
		   the threshold of depth,
		   the length of flanking bases,
		   the path to output,
		   whether to rewrite the existing files (force)
	
	output: the merged depth file in bed format

	return: the dictionary containing the merged depth bed file
	"""

	if os.path.exists(f'{directory}/{prefix}.{threshold}.depth.bed') and force == False:
		print(f'ERROR!!! The file "{directory}/{prefix}.{threshold}.depth.bed" exists\nPlease using "-f" or "--force" to rewrite', file=sys.stderr)
		raise SystemExit
	
	merged_depths_bed = {target:[] for target in depths.keys()}
	with open(f'{directory}/{prefix}.{threshold}.depth.bed', 'w') as f:
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


def complement_merged_depth(merged_depths_bed={}, targets_length={}, flank_len=10):
	"""
	usage: generate the complement of the merged_depth

	input: merged_depths_bed generated by the function merge_depth(),
		   targets_length generated by the function filter(),
		   the length of flanking bases
		     
	return: a list containing the sorted lengths of the complement
	"""

	lengths_com_merged_depth = []
	for target, length in targets_length.items():
		last = flank_len
		n = len(merged_depths_bed[target])
		if n > 0:
			for i, segment in enumerate(merged_depths_bed[target]):
				if i != n-1:
					if segment[0] > last:
						lengths_com_merged_depth.append(segment[0] - last)
					last = segment[1]
				else:
					if segment[0] > last:
						lengths_com_merged_depth.append(segment[0] - last)
					if (length - flank_len) > segment[1]:
						lengths_com_merged_depth.append(length - flank_len - segment[1])
		else:
			lengths_com_merged_depth.append(length - 2*flank_len)
	lengths_com_merged_depth = sorted(lengths_com_merged_depth, reverse=True)
	return lengths_com_merged_depth


def compute_index(targets_length={}, prefix='GCI', directory='.', force=False, merged_depths_bed_list=[], type_list=[], flank_len=10, dist_percent=0.001):
	"""
	usage: remove the regions with depth lower than the threshold and compute the index

	input: targets_length generated by the function filter(),
		   the prefix of the output gci file,
		   the path to output,
		   whether to rewrite the existing files (force),
		   a list of merged_depths_bed generated by merge_depth(),
		   a list of the type of reads,
		   the length of flanking bases,
		   the percentage of the distance between the gap intervals in the chromosome (contig)
		   
	output: an index file containing the reads type, expected N50, observed N50, expected number of contigs, observed number of contigs, genome continuity index
	"""

	if os.path.exists(f'{directory}/{prefix}.gci') and force == False:
		print(f'ERROR!!! The file "{directory}/{prefix}.gci" exists\nPlease using "-f" or "--force" to rewrite', file=sys.stderr)
		raise SystemExit
	with open(f'{directory}/{prefix}.gci', 'w') as f:
		f.write('Type of reads\tExpected N50\tObserved N50\tExpected number of contigs\tObserved number of contigs\tGenome Continuity Index\n')
	

	exp_lengths = sorted([length for length in targets_length.values()], reverse=True)
	exp_cum = np.cumsum(exp_lengths)
	for i, number in enumerate(exp_cum):
		if number >= exp_cum[-1] / 2:
			exp_n50 = exp_lengths[i]
			break
	exp_num_ctg = len(exp_lengths)
	

	for j, merged_depths_bed in enumerate(merged_depths_bed_list):
		obs_lengths = complement_merged_depth(merged_depths_bed, targets_length, flank_len)
		obs_cum = np.cumsum(obs_lengths)
		for i, number in enumerate(obs_cum):
			if number >= obs_cum[-1] / 2:
				obs_n50 = obs_lengths[i]
				break
		
		new_merged_depths_bed = {}
		for target, length in targets_length.items():
			new_merged_depths_bed[target] = []
			dist = length * dist_percent
			current_segment = (flank_len, flank_len)
			for segment in merged_depths_bed[target]:
				if (segment[0] - current_segment[1]) <= dist:
					current_segment = (current_segment[0], segment[1])
				else:
					new_merged_depths_bed[target].append(current_segment)
					current_segment = segment
			if (length - flank_len - current_segment[1]) <= dist:
				current_segment = (current_segment[0], length - flank_len)
			new_merged_depths_bed[target].append(current_segment)
		new_obs_lengths = complement_merged_depth(new_merged_depths_bed, targets_length, flank_len)
		obs_num_ctg = len(new_obs_lengths)


		with open(f'{directory}/{prefix}.gci', 'a') as f:
			f.write(f'{type_list[j]}\t{exp_n50}\t{obs_n50}\t{exp_num_ctg}\t{obs_num_ctg}\t{100 * log2(obs_n50/exp_n50 + 1) / log2(obs_num_ctg/exp_num_ctg + 1)}\n')


def GCI(hifi=[], nano=[], directory='.', prefix='GCI', threads=1, map_qual=30, iden_percent=0.9, ovlp_percent=0.9, clip_percent=0.1, flank_len=10, threshold=0, plot=False, force=False, generate=False, dist_percent=0.001):
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


	if prefix.endswith('/'):
		print(f'ERROR!!! The prefix "{prefix}" is not allowed', file=sys.stderr)
		raise SystemExit
	

	hifi_bam = []
	hifi_paf = []
	nano_bam = []
	nano_paf = []
	if hifi != None:
		for file in hifi:
			CompletedProcess = subprocess.run(f'file -r {file}', shell=True, check=True, capture_output=True)
			if 'gzip' in str(CompletedProcess.stdout):
				hifi_bam.append(file)
			else:
				hifi_paf.append(file)
	if nano != None:
		for file in nano:
			CompletedProcess = subprocess.run(f'file -r {file}', shell=True, check=True, capture_output=True)
			if 'gzip' in str(CompletedProcess.stdout):
				nano_bam.append(file)
			else:
				nano_paf.append(file)
	

	if nano == None:
		depths, targets_length = filter(hifi_paf, hifi_bam, prefix, map_qual, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, generate)
		merged_depth_bed = merge_depth(depths, prefix, threshold, flank_len, directory, force)
		compute_index(targets_length, prefix, directory, force, [merged_depth_bed], ['HiFi'], flank_len, dist_percent)

	elif hifi == None:
		depths, targets_length = filter(nano_paf, nano_bam, prefix, map_qual, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, generate)
		merged_depth_bed = merge_depth(depths, prefix, threshold, flank_len, directory, force)
		compute_index(targets_length, prefix, directory, force, [merged_depth_bed], ['Nano'], flank_len, dist_percent)
	
	else:
		hifi_refs = []
		nano_refs = []
		for samfile in hifi_bam:
			hifi_samfile = pysam.AlignmentFile(samfile, 'rb')
			hifi_refs += hifi_samfile.references
			hifi_samfile.close()
		for samfile in nano_bam:
			nano_samfile = pysam.AlignmentFile(samfile, 'rb')
			nano_refs += nano_samfile.references
			nano_samfile.close()
		
		if set(hifi_refs) != set(nano_refs):
			print(f'ERROR!!! The targets in hifi and nano alignment files are inconsistent\nPlease check the reference used in mapping both hifi and ont reads', file=sys.stderr)
			raise SystemExit
		

		hifi_depths, targets_length = filter(hifi_paf, hifi_bam, prefix, map_qual, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, generate)
		nano_depths, targets_length = filter(nano_paf, nano_bam, prefix, map_qual, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, generate)
		merged_two_type_depths = merge_two_type_depth(hifi_depths, nano_depths, prefix+'_two_type', directory, force, generate)

		hifi_merged_depth_bed = merge_depth(hifi_depths, prefix+'_hifi', threshold, flank_len, directory, force)
		nano_merged_depth_bed = merge_depth(nano_depths, prefix+'_nano', threshold, flank_len, directory, force)
		two_type_merged_depth_bed = merge_depth(merged_two_type_depths, prefix+'_two_type', threshold, flank_len, directory, force)
		compute_index(targets_length, prefix, directory, force, [hifi_merged_depth_bed, nano_merged_depth_bed, two_type_merged_depth_bed], ['HiFi', 'Nano', 'HiFi + Nano'], flank_len, dist_percent)


if __name__=='__main__':
	###########################
	### version = 1.0
	### this version have some limits: -t, -p
	###########################
	version = '1.0'

	parser = argparse.ArgumentParser(prog=sys.argv[0], add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description='A program for assessing the T2T genome', epilog='Examples:\npython GCI.py --hifi hifi.bam hifi.paf ... --nano nano.bam nano.paf ...')

	group_io = parser.add_argument_group("Input/Output")
	group_io.add_argument('--hifi', nargs='+', metavar='', help='PacBio HiFi reads alignment files (at least one bam file)')
	group_io.add_argument('--nano', nargs='+', metavar='', help='Oxford Nanopore long reads alignment files (at least one bam file)')
	group_io.add_argument('-dp', '--dist-percent', nargs='?', metavar='FLOAT', type=float, help='The percentage of the distance between the candidate gap intervals in the whole chromosome (contig) [0.001]', const=0.001, default=0.001)
	group_io.add_argument('-d', nargs='?', dest='directory', metavar='PATH', help='The directory of output files [.]', const='.', default='.')
	group_io.add_argument('-o', '--output', nargs='?', dest='prefix', metavar='STR', help='Prefix of output files [GCI]', const='GCI', default='GCI')
	group_io.add_argument('-t', '--threads', nargs='?', metavar='INT', type=int, help='Number of threads [1]', const=1, default=1)

	group_fo = parser.add_argument_group("Filter Options")
	group_fo.add_argument('-mq', '--map-qual', nargs='?', metavar='INT', type=int, help='Minium mapping quality for alignments [30]', const=30, default=30)
	group_fo.add_argument('-ip', '--iden-percent', nargs='?', metavar='FLOAT', type=float, help='Minimum identity (num_match_res/len_aln) of the reads [0.9]', const=0.9, default=0.9)
	group_fo.add_argument('-op', '--ovlp-percent', nargs='?', metavar='FLOAT', type=float, help='Minimum overlapping percentage of the reads if inputting more than one alignment files [0.9]', const=0.9, default=0.9)
	group_fo.add_argument('-cp', '--clip-percent', nargs='?', metavar='FLOAT', type=float, help='Maximum clipped percentage of the reads [0.1]', const=0.1, default=0.1)
	group_fo.add_argument('-fl', '--flank-len', nargs='?', metavar='INT', type=int, help='The flanking length of the clipped bases [10]', const=10, default=10)
	group_fo.add_argument('-ts', '--threshold', nargs='?', metavar='INT', type=int, help='The threshold of depth in the final bed file [0]', const=0, default=0)

	group_op = parser.add_argument_group("Other Options")
	group_op.add_argument('-g', '--generate', action='store_const', help='Generate the depth files', const=True, default=False)
	group_op.add_argument('-f', '--force', action='store_const', help='Force rewriting of existing files', const=True, default=False)
	group_op.add_argument('-p', '--plot', action='store_const', help='Visualize the final result', const=True, default=False)
	group_op.add_argument('-h', '--help', action="help", help="Show this help message and exit")
	group_op.add_argument('-v', '--version', action="version", version=version, help="Show program's version number and exit")
	
	args = vars(parser.parse_args())
	print(f'Used arguments:{args}')

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
