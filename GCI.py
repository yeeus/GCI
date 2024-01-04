import sys
import pysam
import numpy as np
import argparse
import os
from math import log2
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import ScalarFormatter


def get_average_identity(alns):
	"""
	usage: get average identity if there are many alignment blocks of one target 

	input: synteny[query][target] from function filter()

	return: the average identity
	"""
	###! just average without weights
	tmp = []
	for a in alns:
		tmp.append(a[-1])
	average = sum(tmp) / len(alns)
	return average


def merge_alns_properties(alns, x, y):
	"""
	usage: merge overlapped aligned blocks of either query or target based on the inputted parameters x and y

	input: synteny[query][target] from function filter(),
		   query (1, 2) or target (3, 4)

	return: mapped_length of the query (inputting 1, 2), the leftmost and rightmost position of target (inputting 3, 4)
	"""

	# extract bed list form alns
	bed_list = []
	for a in alns:
		bed_list.append([a[x], a[y]])
	sort_bed_list = sorted(bed_list)
	# merge bed to generate non-overlap blocks
	
	target_length = []
	mapped_length = 0
	low_est = sort_bed_list[0][0]
	high_est = sort_bed_list[0][1]

	for index, block in enumerate(sort_bed_list):
		low, high = block
		if high_est >= low:
			if high_est < high:
				high_est = high
		else:
			target_length.append((high_est-low_est, low_est, high_est))
			mapped_length += (high_est - low_est)
			low_est, high_est = sort_bed_list[index]
	target_length.append((high_est-low_est, low_est, high_est))
	mapped_length += (high_est - low_est)

	target_length = sorted(target_length, key=lambda x:x[0], reverse=True)
	return mapped_length, target_length[0][1], target_length[0][2]


def filter(paf_files=[], bam_files=[], prefix='GCI', map_qual=30, mq_cutoff=50, iden_percent=0.9, clip_percent=0.1, ovlp_percent=0.9, flank_len=15, directory='.', force=False, generate=False):
	"""
	usage: filter the paf and bam file(s) based on many metrics(, and finally generate the depth file)

	input: paf file(s),
		   bam file(s),
		   the prefix of output depth file,
		   the filtered mapping quality,
		   the cutoff of mapping quality,
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
	
	high_qual_querys = set()
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
						if mapq >= mq_cutoff:
							high_qual_querys.add(query)
				for query in synteny.keys():
					mapping_results = {}
					for target in synteny[query].keys():
						alns = synteny[query][target]
						non_overlap_qry_aligned, _, _ = merge_alns_properties(alns, 1, 2)
						query_length = alns[0][0]
						alignrate = non_overlap_qry_aligned / query_length
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
					if segment.mapping_quality >= mq_cutoff:
						high_qual_querys.add(segment.query_name)
		samfile.close()


	files = paf_lines + samfile_dicts
	if len(files) > 1:
		files_sets = []
		for file in files:
			files_sets.append(set(file.keys()))
		comm_querys = set.intersection(*files_sets)

		final_querys = high_qual_querys | comm_querys
		file1 = {query:segment for query, segment in files[0].items() if query in final_querys}
		for file in files[1:]:
			for query, segment in file.items():
				if query in file1.keys():
					segment1 = file1[query]
					if segment[0] == segment1[0]:
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
						del file1[query]
				elif query in high_qual_querys:
					file1.update({query:(segment[0], segment[1], segment[2])})
	else:
		file1 = files[0]
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


def collapse_depth_range(depths={}, leftmost=-1, rightmost=0, flank_len=15):
	"""
	usage: collapse positions with depth in the range (leftmost, rightmost]

	input: the whole-genome depth dictionary generated by filter() and merge_two_type_depth(),
		   the leftmost threshold of depth,
		   the rightmost threshold of depth,
		   the length of flanking bases

	return: the dictionary containing the merged depth bed file
	"""

	merged_depths_bed = {target:[] for target in depths.keys()}
	for target, depth_list in depths.items():
		start_flag = 0
		end_flag = 1
		chr_len = len(depth_list)
		for i, depth in enumerate(depth_list[flank_len:chr_len-flank_len]):
			if leftmost < depth <= rightmost:
				if start_flag == 0:
					start = i + flank_len
					start_flag = 1
					end_flag = 0
				if i == (chr_len - flank_len*2 - 1):
					end = i + flank_len + 1
					merged_depths_bed[target].append((start, end))
			else:
				if end_flag == 0:
					if i > flank_len: #! look better
						end = i + flank_len
						merged_depths_bed[target].append((start, end))
					end_flag = 1
					start_flag = 0
	return merged_depths_bed


def merge_depth(depths={}, prefix='GCI', threshold=0, flank_len=15, directory='.', force=False):
	"""
	usage: merge positions with depth lower than the threshold (used in the main function and based on the function collapse_depth_range)

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
	
	merged_depths_bed = collapse_depth_range(depths, -1, threshold, flank_len)
	with open(f'{directory}/{prefix}.{threshold}.depth.bed', 'w') as f:
		for target, segments in merged_depths_bed.items():
			for segment in segments:
				f.write(f'{target}\t{segment[0]}\t{segment[1]}\n')
	
	return merged_depths_bed


def complement_merged_depth(merged_depths_bed={}, targets_length={}, flank_len=15):
	"""
	usage: generate the complement of the merged_depth

	input: merged_depths_bed generated by the function merge_depth(),
		   targets_length generated by the function filter(),
		   the length of flanking bases

	return: a dict containing a list with the content of the sorted lengths of the complement
	"""
	
	lengths_com_merged_depth_dict = {}
	for target, length in targets_length.items():
		lengths_com_merged_depth = []
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
		lengths_com_merged_depth_dict.update({target:lengths_com_merged_depth})
	return lengths_com_merged_depth_dict


def compute_n50(lengths=[]):
	"""
	usage: compute n50

	input: a list of the lengths

	return: n50
	"""
	n50 = 0
	lengths = sorted(lengths, reverse=True)
	cum = np.cumsum(lengths)
	for i, number in enumerate(cum):
		if number >= cum[-1] / 2:
			n50 = lengths[i]
			break
	return n50


def compute_index(targets_length={}, prefix='GCI', directory='.', force=False, merged_depths_bed_list=[], type_list=[], flank_len=15, dist_percent=0.005):
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
		pass

	
	exp_lengths = [length for length in targets_length.values()]
	exp_n50 = compute_n50(exp_lengths)
	exp_num_ctg = len(exp_lengths)
	exp_n50_dict = dict(targets_length)
	exp_n50_dict.update({'Genome':exp_n50})
	exp_num_ctg_dict = {target:1 for target in targets_length.keys()}
	exp_num_ctg_dict.update({'Genome':exp_num_ctg})


	for i, merged_depths_bed in enumerate(merged_depths_bed_list):
		obs_lengths_dict = complement_merged_depth(merged_depths_bed, targets_length, flank_len)
		obs_lengths = [item for value in obs_lengths_dict.values() for item in value]
		obs_n50 = compute_n50(obs_lengths)
		obs_n50_dict = {target:compute_n50(lengths) for target, lengths in obs_lengths_dict.items()}
		obs_n50_dict.update({'Genome':obs_n50})
		
		
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
		new_obs_lengths_dict = complement_merged_depth(new_merged_depths_bed, targets_length, flank_len)
		new_obs_lengths = [item for value in new_obs_lengths_dict.values() for item in value]
		obs_num_ctg = len(new_obs_lengths)
		obs_num_ctg_dict = {target:len(lengths) for target, lengths in new_obs_lengths_dict.items()}
		obs_num_ctg_dict.update({'Genome':obs_num_ctg})


		with open(f'{directory}/{prefix}.gci', 'a') as f:
			f.write(f'{type_list[i]}:\n')
			f.write('Chromosome\tExpected N50\tObserved N50\tExpected number of contigs\tObserved number of contigs\tGenome Continuity Index\n')
			for target in exp_n50_dict.keys():
				exp_n50 = exp_n50_dict[target]
				obs_n50 = obs_n50_dict[target]
				exp_num_ctg = exp_num_ctg_dict[target]
				obs_num_ctg = obs_num_ctg_dict[target]
				f.write(f'{target}\t{exp_n50}\t{obs_n50}\t{exp_num_ctg}\t{obs_num_ctg}\t{100 * log2(obs_n50/exp_n50 + 1) / log2(obs_num_ctg/exp_num_ctg + 1):.4f}\n')
			f.write('----------------------------------------------------------------------------------------------------------------------------------------\n\n\n')


def sliding_window_average_depth(depths=[], window_size=1, max_depth=None):
	"""
	usage: get the averaged depths via sliding window

	input: the single chromosome depth list from the whole-genome depth dictionary generated by filter(),
		   the window size in bytes,
		   the max depth to plot

	return: the positions and averaged depths list
	"""
	averaged_positions = []
	averaged_depths = []
	window_depths = []
	for i, depth in enumerate(depths):
		if depth == 0:
			if len(window_depths) > 0:
				average_depth = sum(window_depths) / len(window_depths)
				if average_depth > max_depth:
					average_depth = max_depth
				averaged_depths.append(average_depth)
				averaged_positions.append(i-1)
				window_depths = []
			averaged_depths.append(0)
			averaged_positions.append(i)
		else:
			window_depths.append(depth)
			if len(window_depths) == window_size:
				average_depth = sum(window_depths) / window_size
				if average_depth > max_depth:
					average_depth = max_depth
				averaged_depths.append(average_depth)
				averaged_positions.append(i)
				window_depths = []
	if len(window_depths) > 0:
		average_depth = sum(window_depths) / len(window_depths)
		if average_depth > max_depth:
			average_depth = max_depth
		averaged_depths.append(average_depth)
		averaged_positions.append(i)
	return averaged_positions, averaged_depths


def plot_depth(depths_list=[], depth_min=0.1, depth_max=4.0, window_size=0.001, image_type='png', directory='.', prefix='GCI', force=False):
	"""
	usage: plot whole genome depth

	input: a list of the whole-genome depth dictionary generated by filter(),
		   the cutoff in folds of depth to plot,
		   the max folds of depth to plot,
		   the window size in chromosome units (0-1) when plotting,
		   the format of output images,
		   the path to output,
		   the prefix of the output gci file,
		   whether to rewrite the existing files (force)

	output: the whole genome depth plots
	"""
	if image_type == 'pdf' or image_type == 'png':
		for target in depths_list[0].keys():
			if os.path.exists(f'{directory}/images/{prefix}.{target}.{image_type}') and force == False:
				print(f'ERROR!!! The file "{directory}/images/{prefix}.{target}.{image_type}" exists\nPlease using "-f" or "--force" to rewrite', file=sys.stderr)
				raise SystemExit
	else:
		print(f'ERROR!!! The format of output images only supports pdf and png', file=sys.stderr)
		raise SystemExit
	
	if len(depths_list) == 1:
		sum_depths = []
		for depths in depths_list[0].values():
			sum_depths = np.concatenate((sum_depths, depths))
		mean_depth = np.mean(sum_depths)
	elif len(depths_list)  == 2:
		sum_depths = []
		for depths in depths_list[0].values():
			sum_depths = np.concatenate((sum_depths, depths))
		mean_depth1 = np.mean(sum_depths)
		sum_depths = []
		for depths in depths_list[-1].values():
			sum_depths = np.concatenate((sum_depths, depths))
		mean_depth2 = np.mean(sum_depths)

	if len(depths_list) == 1:
		max_depth = mean_depth * depth_max
		averaged_dict = {}
		max_averaged_depths = []
		for target in depths_list[0].keys():
			depths = depths_list[0][target]
			averaged_positions, averaged_depths = sliding_window_average_depth(depths, max(1, round(len(depths) * window_size)), max_depth)
			averaged_dict.update({target:(averaged_positions, averaged_depths)})
			max_averaged_depths.append(max(averaged_depths))
		y_high = max(max_averaged_depths) + 10
	elif len(depths_list)  == 2:
		depth_max1 = mean_depth1 * depth_max
		depth_max2 = mean_depth2 * depth_max
		averaged_dict1 = {}
		averaged_dict2 = {}
		max_averaged_depths1 = []
		max_averaged_depths2 = []
		for target in depths_list[0].keys():
			depths1 = depths_list[0][target]
			averaged_positions1, averaged_depths1 = sliding_window_average_depth(depths1, max(1, round(len(depths1) * window_size)), depth_max1)
			averaged_dict1.update({target:(averaged_positions1, averaged_depths1)})
			max_averaged_depths1.append(max(averaged_depths1))
			depths2 = depths_list[-1][target]
			averaged_positions2, averaged_depths2 = sliding_window_average_depth(depths2, max(1, round(len(depths2) * window_size)), depth_max2)
			averaged_dict2.update({target:(averaged_positions2, averaged_depths2)})
			max_averaged_depths2.append(max(averaged_depths2))
		y_max = max(max_averaged_depths1) + 10
		y_min = max(max_averaged_depths2) + 10


	for target in depths_list[0].keys():
		if len(depths_list) == 1:
			fig, ax = plt.subplots(figsize=(20, 4))

			depths = depths_list[0][target]
			averaged_positions, averaged_depths = averaged_dict[target]
			ax.stackplot(averaged_positions, averaged_depths, lw=0.8, color='#2ca25f', zorder=4)
			ax.axhline(mean_depth, color="r", ls='--', dash_capstyle='butt', lw=1, zorder=5)
			
			ax.set_ylim(top=y_high)
			merged_min_bed = collapse_depth_range({target:depths}, 0, mean_depth * depth_min, 0)
			for segment in merged_min_bed[target]:
				ax.axvspan(segment[0], segment[1], facecolor='#B7DBEA')
			merged_0_bed = collapse_depth_range({target:depths}, -1, 0, 0)
			for segment in merged_0_bed[target]:
				ax.axvspan(segment[0], segment[1], facecolor='#FAD7DD')
			
			ax.xaxis.set_minor_locator(AutoMinorLocator())
			ax.xaxis.set_minor_formatter(ScalarFormatter())
			ax.yaxis.set_minor_locator(AutoMinorLocator())

		elif len(depths_list)  == 2:
			fig, ax = plt.subplots(figsize=(20, 8))

			depths1 = depths_list[0][target]
			averaged_positions1, averaged_depths1 = averaged_dict1[target]
			ax.stackplot(averaged_positions1, averaged_depths1, lw=0.8, color='#2ca25f', zorder=4)
			ax.axhline(mean_depth1, color="r", ls='-.', dash_capstyle='butt', lw=1, zorder=5)
			
			ax.axhline(0, color="black")
			depths2 = depths_list[-1][target]
			averaged_positions2, averaged_depths2 = averaged_dict2[target]
			ax.stackplot(averaged_positions2, -np.array(averaged_depths2), lw=0.8, color='#3C5488', zorder=4)
			ax.axhline(-mean_depth2, color="r", ls='-.', dash_capstyle='butt', lw=1, zorder=5)
			
			ax.set_ylim(bottom=-y_min, top=y_max)
			y_frac = y_min / (y_max + y_min)
			merged_min_bed = collapse_depth_range({target:depths1}, 0, mean_depth1 * depth_min, 0)
			for segment in merged_min_bed[target]:
				ax.axvspan(segment[0], segment[1], y_frac, 1, facecolor='#B7DBEA')
			merged_0_bed = collapse_depth_range({target:depths1}, -1, 0, 0)
			for segment in merged_0_bed[target]:
				ax.axvspan(segment[0], segment[1], y_frac, 1, facecolor='#FAD7DD')
			merged_min_bed = collapse_depth_range({target:depths2}, 0, mean_depth2 * depth_min, 0)
			for segment in merged_min_bed[target]:
				ax.axvspan(segment[0], segment[1], 0, y_frac, facecolor='#B7DBEA')
			merged_0_bed = collapse_depth_range({target:depths2}, -1, 0, 0)
			for segment in merged_0_bed[target]:
				ax.axvspan(segment[0], segment[1], 0, y_frac, facecolor='#FAD7DD')
			
			hifi_line = mlines.Line2D([], [], color='#2ca25f', label='HiFi', lw=0.8)
			nano_line = mlines.Line2D([], [], color='#3C5488', label='Nano', lw=0.8)
			legend1 = plt.legend(handles=[hifi_line, nano_line], loc='upper left')
			plt.gca().add_artist(legend1)
			ax.xaxis.set_minor_locator(AutoMinorLocator())
			ax.xaxis.set_minor_formatter(ScalarFormatter())
			ax.yaxis.set_minor_locator(AutoMinorLocator())


		merged_min_line = mlines.Line2D([], [], color='#B7DBEA', label=f'The region with the depth in the range of (0, {depth_min}*mean_depth]')
		merged_0_line = mlines.Line2D([], [], color='#FAD7DD', label='The region of zero depth')
		mean_line = mlines.Line2D([], [], color="r", ls='-.', dash_capstyle='butt', lw=1, label='Mean Coverage')
		legend2 = plt.legend(handles=[mean_line, merged_min_line, merged_0_line], loc='lower center', bbox_to_anchor=(0.5, 1), ncols=3)
		plt.gca().add_artist(legend2)
		plt.title(f'Filtered depth across the whole genome:{target}', fontsize=18, pad=30)
		plt.xlabel('Genomic Position (bp)', fontsize=14)
		plt.ylabel('Depth', fontsize=14)
		plt.xticks(fontsize=12)
		plt.yticks(fontsize=12)
		plt.tight_layout()
		plt.savefig(f'{directory}/images/{prefix}.{target}.{image_type}', dpi=200)
		plt.close()


def GCI(hifi=[], nano=[], directory='.', prefix='GCI', threads=1, map_qual=30, mq_cutoff=50, iden_percent=0.9, ovlp_percent=0.9, clip_percent=0.1, flank_len=15, threshold=0, plot=False, depth_min=0.1, depth_max=4.0, window_size=0.001, image_type='png', force=False, generate=False, dist_percent=0.005):
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
	
	if plot == True:
		if os.path.exists(f'{directory}/images'):
			if not os.access(f'{directory}/images', os.R_OK):
				print(f'ERROR!!! The path "{directory}/images" is unable to read', file=sys.stderr)
				raise SystemExit
			if not os.access(f'{directory}/images', os.W_OK):
				print(f'ERROR!!! The path "{directory}/images" is unable to write', file=sys.stderr)
				raise SystemExit
		else:
			os.makedirs(f'{directory}/images')
		
		image_type = image_type.lower()
	

	hifi_bam = []
	hifi_paf = []
	nano_bam = []
	nano_paf = []
	if hifi != None:
		for file in hifi:
			if file.endswith('.bam'):
				hifi_bam.append(file)
			else:
				hifi_paf.append(file)
	if nano != None:
		for file in nano:
			if file.endswith('.bam'):
				nano_bam.append(file)
			else:
				nano_paf.append(file)
	

	if nano == None:
		depths, targets_length = filter(hifi_paf, hifi_bam, prefix, map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, generate)
		merged_depth_bed = merge_depth(depths, prefix, threshold, flank_len, directory, force)
		compute_index(targets_length, prefix, directory, force, [merged_depth_bed], ['HiFi'], flank_len, dist_percent)
		if plot == True:
			plot_depth([depths], depth_min, depth_max, window_size, image_type, directory, prefix, force)

	elif hifi == None:
		depths, targets_length = filter(nano_paf, nano_bam, prefix, map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, generate)
		merged_depth_bed = merge_depth(depths, prefix, threshold, flank_len, directory, force)
		compute_index(targets_length, prefix, directory, force, [merged_depth_bed], ['Nano'], flank_len, dist_percent)
		if plot == True:
			plot_depth([depths], depth_min, depth_max, window_size, image_type, directory, prefix, force)
	
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
		

		hifi_depths, targets_length = filter(hifi_paf, hifi_bam, prefix+'_hifi', map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, generate)
		nano_depths, targets_length = filter(nano_paf, nano_bam, prefix+'_nano', map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, generate)
		merged_two_type_depths = merge_two_type_depth(hifi_depths, nano_depths, prefix+'_two_type', directory, force, generate)

		hifi_merged_depth_bed = merge_depth(hifi_depths, prefix+'_hifi', threshold, flank_len, directory, force)
		nano_merged_depth_bed = merge_depth(nano_depths, prefix+'_nano', threshold, flank_len, directory, force)
		two_type_merged_depth_bed = merge_depth(merged_two_type_depths, prefix+'_two_type', threshold, flank_len, directory, force)
		compute_index(targets_length, prefix, directory, force, [hifi_merged_depth_bed, nano_merged_depth_bed, two_type_merged_depth_bed], ['HiFi', 'Nano', 'HiFi + Nano'], flank_len, dist_percent)
		if plot == True:
			plot_depth([hifi_depths, nano_depths], depth_min, depth_max, window_size, image_type, directory, prefix, force)


if __name__=='__main__':
	###########################
	### version = 0.1
	### this version have some limits: -t
	###########################
	version = 'GCI version 0.1'

	parser = argparse.ArgumentParser(prog=sys.argv[0], add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description='A program for assessing the T2T genome', epilog='Examples:\npython GCI.py --hifi hifi.bam hifi.paf ... --nano nano.bam nano.paf ...')

	group_io = parser.add_argument_group("Input/Output")
	group_io.add_argument('--hifi', nargs='+', metavar='', help='PacBio HiFi reads alignment files (at least one bam file)')
	group_io.add_argument('--nano', nargs='+', metavar='', help='Oxford Nanopore long reads alignment files (at least one bam file)')
	group_io.add_argument('-ts', '--threshold', metavar='INT', type=int, help='The threshold of depth in the final bed file [0]', default=0)
	group_io.add_argument('-dp', '--dist-percent', metavar='FLOAT', type=float, help='The percentage of the distance between the candidate gap intervals in the whole chromosome (contig) [0.005]', default=0.005)
	group_io.add_argument('-d', dest='directory', metavar='PATH', help='The directory of output files [.]', default='.')
	group_io.add_argument('-o', '--output', dest='prefix', metavar='STR', help='Prefix of output files [GCI]', default='GCI')
	#group_io.add_argument('-t', '--threads', metavar='INT', type=int, help='Number of threads [1]', default=1)

	group_fo = parser.add_argument_group("Filter Options")
	group_fo.add_argument('-mq', '--map-qual', metavar='INT', type=int, help='Minium mapping quality for alignments [30]', default=30)
	group_fo.add_argument('--mq-cutoff', metavar='INT', type=int, help='The cutoff of mapping quality for keeping the alignment [50]', default=50)
	group_fo.add_argument('-ip', '--iden-percent', metavar='FLOAT', type=float, help='Minimum identity (num_match_res/len_aln) of the reads [0.9]', default=0.9)
	group_fo.add_argument('-op', '--ovlp-percent', metavar='FLOAT', type=float, help='Minimum overlapping percentage of the reads if inputting more than one alignment files [0.9]', default=0.9)
	group_fo.add_argument('-cp', '--clip-percent', metavar='FLOAT', type=float, help='Maximum clipped percentage of the reads [0.1]', default=0.1)
	group_fo.add_argument('-fl', '--flank-len', metavar='INT', type=int, help='The flanking length of the clipped bases [15]', default=15)

	group_po = parser.add_argument_group("Plot Options")
	group_po.add_argument('-p', '--plot', action='store_const', help='Visualize the finally filtered whole genome depth', const=True, default=False)
	group_po.add_argument('-dmin', '--depth-min', metavar='FLOAT', type=float, help='Minimum depth in folds of mean coverage for plotting [0.1]', default=0.1)
	group_po.add_argument('-dmax', '--depth-max', metavar='FLOAT', type=float, help='Maximum depth in folds of mean coverage for plotting [4.0]', default=4.0)
	group_po.add_argument('-ws', '--window-size', metavar='FLOAT', type=float, help='The window size in chromosome units (0-1) when plotting [0.001]', default=0.001)
	group_po.add_argument('-it', '--image-type', metavar='STR', help='The format of the output images: png or pdf [png]', default='png')

	group_op = parser.add_argument_group("Other Options")
	group_op.add_argument('-g', '--generate', action='store_const', help='Generate the depth files', const=True, default=False)
	group_op.add_argument('-f', '--force', action='store_const', help='Force rewriting of existing files', const=True, default=False)
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
				if file.endswith('.bam'):
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
				if file.endswith('.bam'):
					bam_num += 1
			else:
				print(f'ERROR!!! "{file}" is not an available file', file=sys.stderr)
				raise SystemExit
		if bam_num == 0:
			print('ERROR!!! Please input at least one Oxford Nanopore long reads bam file\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit
	
	if args['map_qual'] > args['mq_cutoff']:
		print(f'WARNING!!! The minium mapping quality is {args["map_qual"]} and higher than the cutoff {args["mq_cutoff"]}, which means that wouldn\'t filter any reads\nPlease read the help message using "-h" or "--help"', file=sys.stdout)

	GCI(**args)
	
