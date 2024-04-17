import sys
import pysam
import numpy as np
import argparse
import os
import subprocess
import re
import gzip
from Bio import SeqIO
from math import log2
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import ScalarFormatter


def get_Ns_ref(reference=None, prefix='GCI', directory='.', force=False):
    """
    usage: get Ns (gaps) of reference and return the bed file (if have)

    input: the reference file

    output: the gaps bed file

    return: the dictionary containing the gaps bed file and the file path
    """
    Ns_bed = {}
    pattern = re.compile(r'(?i)N+')
    for record in SeqIO.parse(reference, 'fasta'):
        for match in pattern.finditer(str(record.seq)):
            target = record.id
            if target not in Ns_bed.keys():
                Ns_bed[target] = []
            Ns_bed[target].append((match.start(), match.end()))
    
    if len(Ns_bed) > 0:
        if os.path.exists(f'{directory}/{prefix}.gaps.bed') and force == False:
            sys.exit(f'ERROR!!! The file "{directory}/{prefix}.gaps.bed" exists\nPlease use "-f" or "--force" to rewrite')
        with open(f'{directory}/{prefix}.gaps.bed', 'w') as f:
            for target, segments in Ns_bed.items():
                for segment in segments:
                    f.write(f'{target}\t{segment[0]}\t{segment[1]}\n')
        return Ns_bed, f'{directory}/{prefix}.gaps.bed'
    else:
        return None, None
    

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


def filter(paf_files=[], bam_files=[], prefix='GCI', map_qual=30, mq_cutoff=50, iden_percent=0.9, clip_percent=0.1, ovlp_percent=0.9, flank_len=15, directory='.', force=False, log_reads_type=''):
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
           whether to rewrite the existing files (force)

    output: the gzipped whole-genome depth file

    return: the whole-genome depth dictionary,
            one dictionary keyed by the targets with the length value
    """
    if os.path.exists(f'{directory}/{prefix}.depth.gz') and force == False:
        sys.exit(f'ERROR!!! The file "{directory}/{prefix}.depth.gz" exists\nPlease use "-f" or "--force" to rewrite')
    print(f'Filtering {log_reads_type} alignment files ...')
    
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

    print(f'Filtering {log_reads_type} alignment files done!!!')
    print(f'Writing depths into "{directory}/{prefix}.depth.gz" ...')
    with gzip.open(f'{directory}/{prefix}.depth.gz', 'wb') as f:
        for target, depth_list in depths.items():
            content = f'>{target}\n'
            f.write(content.encode('utf-8'))
            for i, depth in enumerate(depth_list):
                content = f'{depth}\n'
                f.write(content.encode('utf-8'))
    print(f'Writing depths done!!!\n\n')
    return depths, targets_length


def merge_gaps_depths(depths={}, Ns_bed=None):
    """
    usage: merge gaps and issues detected by filter()

    input: depths generated by filter(),
           Ns_bed generated by get_Ns_ref()

    return: the merged depths
    """
    if Ns_bed != None:
        for target, segments in Ns_bed.items():
            for segment in segments:
                depths[target][segment[0]:segment[1]] = 0
    return depths


def merge_two_type_depth(hifi_depths={}, nano_depths={}, prefix='GCI', directory='.', force=False):
    """
    usage: merge the depths dictionary generated by two types of long reads,
           the prefix of output depth file,
           the path to output,
           whether to rewrite the existing files (force)

    input: the whole-genome depth dictionaries of two types of long reads generated by filter()

    output: the gzipped whole-genome depth file

    return: the merged whole-genome depth dictionary
    """
    print('Merging HiFi and ONT depth file ...')
    if os.path.exists(f'{directory}/{prefix}.depth.gz') and force == False:
        sys.exit(f'ERROR!!! The file "{directory}/{prefix}.depth.gz" exists\nPlease use "-f" or "--force" to rewrite')

    merged_two_type_depths = {target:[] for target in hifi_depths.keys()}
    for target, hifi_depth_list in hifi_depths.items():
        nano_depth_list = nano_depths[target]
        for i, depth in enumerate(hifi_depth_list):
            merged_two_type_depths[target].append(max(depth, nano_depth_list[i]))
        merged_two_type_depths[target] = np.array(merged_two_type_depths[target])
    

    with gzip.open(f'{directory}/{prefix}.depth.gz', 'wb') as f:
        for target, depth_list in merged_two_type_depths.items():
            content = f'>{target}\n'
            f.write(content.encode('utf-8'))
            for i, depth in enumerate(depth_list):
                content = f'{depth}\n'
                f.write(content.encode('utf-8'))
    print('Merging HiFi and ONT depth file done!!!\n\n')
    return merged_two_type_depths


def collapse_depth_range(depths={}, leftmost=-1, rightmost=0, flank_len=15, start_pos=0):
    """
    usage: collapse positions with depth in the range (leftmost, rightmost]

    input: the whole-genome depth dictionary generated by filter() and merge_two_type_depth(),
           the leftmost threshold of depth,
           the rightmost threshold of depth,
           the length of flanking bases,
           the position of start

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
                    merged_depths_bed[target].append((start+start_pos, end+start_pos))
            else:
                if end_flag == 0:
                    if i > flank_len: #! look better
                        end = i + flank_len
                        merged_depths_bed[target].append((start+start_pos, end+start_pos))
                    end_flag = 1
                    start_flag = 0
    return merged_depths_bed


def merge_depth(depths={}, prefix='GCI', threshold=0, flank_len=15, directory='.', force=False, log_reads_type=''):
    """
    usage: merge positions with depth lower than the threshold (used in the main function and based on the function collapse_depth_range)

    input: the whole-genome depth dictionary generated by filter(), merge_two_type_depth(), and merge_gaps_depths()
           the prefix of output threshold.depth.bed, gaps.bed and final.bed file,
           the threshold of depth,
           the length of flanking bases,
           the path to output,
           whether to rewrite the existing files (force),
           reads type for logging
    
    output: the merged depth file

    return: the dictionary containing the merged depth bed file
    """
    print(f'Getting {log_reads_type} issues bed file detected by GCI ...')
    if os.path.exists(f'{directory}/{prefix}.{threshold}.depth.bed') and force == False:
        sys.exit(f'ERROR!!! The file "{directory}/{prefix}.{threshold}.depth.bed" exists\nPlease use "-f" or "--force" to rewrite')
    
    merged_depths_bed = collapse_depth_range(depths, -1, threshold, flank_len, 0)
    with open(f'{directory}/{prefix}.{threshold}.depth.bed', 'w') as f:
        for target, segments in merged_depths_bed.items():
            for segment in segments:
                f.write(f'{target}\t{segment[0]}\t{segment[1]}\n')
    print(f'Getting {log_reads_type} issues bed file done!!!\n\n')
    return merged_depths_bed


def complement_merged_depth(merged_depths_bed={}, targets_length={}, flank_len=15, start=None, end=None):
    """
    usage: generate the complement of the merged_depth

    input: merged_depths_bed generated by the function merge_depth(),
           targets_length generated by the function filter(),
           the length of flanking bases,
           the position of start,
           the position of end

    return: a dict containing a list with the content of the sorted lengths of the complement
    """
    start_flag = False
    end_flag = False
    if start != None and end != None:
        start_flag = True
        end_flag = True

    lengths_com_merged_depth_dict = {}
    for target, length in targets_length.items():
        if start_flag == False and end_flag == False:
            start = flank_len
            end = length - flank_len
        lengths_com_merged_depth = []
        last = start
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
                    if end > segment[1]:
                        lengths_com_merged_depth.append(end - segment[1])
        else:
            lengths_com_merged_depth.append(end - start)
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


def merge_merged_depth_bed(merged_depths_bed={}, targets_length={}, dist_percent=0.005, flank_len=15, start=None, end=None):
    """
    usage: merge the adjacent intervals with the distance lower than chr_length * dist_percent 

    input: merged_depths_bed generated by merge_depth(),
           targets_length generated by the function filter(),
           the percentage of the distance between the gap intervals in the chromosome,
           the length of flanking bases,
           the position of start,
           the position of end
    
    return: the merged merged_depths_bed
    """
    start_flag = False
    end_flag = False
    if start != None and end != None:
        start_flag = True
        end_flag = True

    new_merged_depths_bed = {}
    for target, length in targets_length.items():
        new_merged_depths_bed[target] = []
        dist = length * dist_percent
        if start_flag == False and end_flag == False:
            start = flank_len
            end = length - flank_len
        current_segment = (start, start)
        for segment in merged_depths_bed[target]:
            if (segment[0] - current_segment[1]) <= dist:
                current_segment = (current_segment[0], segment[1])
            else:
                new_merged_depths_bed[target].append(current_segment)
                current_segment = segment
        if (end - current_segment[1]) <= dist:
            current_segment = (current_segment[0], end)
        new_merged_depths_bed[target].append(current_segment)
    return new_merged_depths_bed


def compute_index(targets_length={}, prefix='GCI', directory='.', force=False, merged_depths_bed_list=[], type_list=[], flank_len=15, dist_percent=0.005, regions_bed={}, depths_list=[], threshold=0):
    """
    usage: remove the regions with depth lower than the threshold and compute the index

    input: targets_length generated by the function filter(),
           the prefix of the output gci file,
           the path to output,
           whether to rewrite the existing files (force),
           a list of merged_depths_bed generated by merge_depth(),
           a list of the type of reads,
           the length of flanking bases,
           the percentage of the distance between the gap intervals in the chromosome,
           the regions bed file,
           a list of depths generated by filter(),
           the threshold of depth
           
    output: an index file containing the reads type, Theoretical maximum N50, Corrected N50, Theoretical minimum contigs number, Corrected contigs number, GCI score,
            and regions gci file containing chromosome, start, end, GCI score for each types reads
    """
    if os.path.exists(f'{directory}/{prefix}.gci') and force == False:
        sys.exit(f'ERROR!!! The file "{directory}/{prefix}.gci" exists\nPlease use "-f" or "--force" to rewrite')
    with open(f'{directory}/{prefix}.gci', 'w') as f:
        pass
    if len(regions_bed) > 0:
        if os.path.exists(f'{directory}/{prefix}.regions.gci') and force == False:
            sys.exit(f'ERROR!!! The file "{directory}/{prefix}.regions.gci" exists\nPlease use "-f" or "--force" to rewrite')
        with open(f'{directory}/{prefix}.regions.gci', 'w') as f:
            f.write('Chromosome\tStart\tEnd\t' + '\t'.join(type_list) + '\n')
    
    print('Computing Theoretical minimum N50 and contigs number ...')
    exp_lengths = [length for length in targets_length.values()]
    exp_n50 = compute_n50(exp_lengths)
    exp_num_ctg = len(exp_lengths)
    exp_n50_dict = dict(targets_length)
    exp_n50_dict.update({'Genome':exp_n50})
    exp_num_ctg_dict = {target:1 for target in targets_length.keys()}
    exp_num_ctg_dict.update({'Genome':exp_num_ctg})
    print('Computing Theoretical minimum N50 and contigs number done!!!')

    for i, merged_depths_bed in enumerate(merged_depths_bed_list):
        print(f'Computing Curated N50 and contigs number for {type_list[i]} ...')
        obs_lengths_dict = complement_merged_depth(merged_depths_bed, targets_length, flank_len)
        obs_lengths = [item for value in obs_lengths_dict.values() for item in value]
        obs_n50 = compute_n50(obs_lengths)
        obs_n50_dict = {target:compute_n50(lengths) for target, lengths in obs_lengths_dict.items()}
        obs_n50_dict.update({'Genome':obs_n50})
        
        
        new_merged_depths_bed = merge_merged_depth_bed(merged_depths_bed, targets_length, dist_percent, flank_len)
        new_obs_lengths_dict = complement_merged_depth(new_merged_depths_bed, targets_length, flank_len)
        new_obs_lengths = [item for value in new_obs_lengths_dict.values() for item in value]
        obs_num_ctg = len(new_obs_lengths)
        obs_num_ctg_dict = {target:len(lengths) for target, lengths in new_obs_lengths_dict.items()}
        obs_num_ctg_dict.update({'Genome':obs_num_ctg})

        print(f'Computing Curated N50 and contigs number for {type_list[i]} done!!!')
        print(f'Writing results to {directory}/{prefix}.gci ...')
        with open(f'{directory}/{prefix}.gci', 'a') as f:
            f.write(f'{type_list[i]}:\n')
            f.write('Chromosome\tTheoretical maximum N50\tCurated N50\tTheoretical minimum contigs number\tCurated contigs number\tGCI score\n')
            for target in exp_n50_dict.keys():
                exp_n50 = exp_n50_dict[target]
                obs_n50 = obs_n50_dict[target]
                exp_num_ctg = exp_num_ctg_dict[target]
                obs_num_ctg = obs_num_ctg_dict[target]
                
                if obs_num_ctg == 0:
                    gci = 0
                else:
                    gci = 100 * log2(obs_n50/exp_n50 + 1) / log2(obs_num_ctg/exp_num_ctg + 1)
                f.write(f'{target}\t{exp_n50}\t{obs_n50}\t{exp_num_ctg}\t{obs_num_ctg}\t{gci:.4f}\n')
            f.write('----------------------------------------------------------------------------------------------------------------------------------------\n\n\n')
        print(f'Writing results to {directory}/{prefix}.gci done!!!\n\n')

    if len(regions_bed) > 0:
        print('Computing GCI scores for regions ...')
        region_all_lengths = []
        region_all_obs_length = [[] for i in range(len(depths_list))]
        region_all_obs_num_ctg = [0 for i in range(len(depths_list))]
        for target, segments in regions_bed.items():
            for segment in segments:
                start = segment[0]
                end = segment[1]
                exp_n50 = end - start
                region_all_lengths.append(exp_n50)
                exp_num_ctg = 1
                gci = []
                for i, depthss in enumerate(depths_list):
                    depths = depthss[target][start:end]
                    merged_depths_bed = collapse_depth_range({target:depths}, -1, threshold, 0, start)
                    obs_lengths_dict = complement_merged_depth(merged_depths_bed, {target:exp_n50}, start, start, end)
                    obs_n50 = compute_n50(obs_lengths_dict[target])
                    region_all_obs_length[i] += obs_lengths_dict[target]
                    new_merged_depths_bed = merge_merged_depth_bed(merged_depths_bed, {target:exp_n50}, dist_percent, start, start, end)
                    new_obs_lengths_dict = complement_merged_depth(new_merged_depths_bed, {target:exp_n50}, start, start, end)
                    obs_num_ctg = len(new_obs_lengths_dict[target])
                    region_all_obs_num_ctg[i] += obs_num_ctg

                    gci.append(round(100 * log2(obs_n50/exp_n50 + 1) / log2(obs_num_ctg/exp_num_ctg + 1), 4))
                with open(f'{directory}/{prefix}.regions.gci', 'a') as f:
                    f.write(f'{target}\t{segment[0]}\t{segment[1]}\t' + '\t'.join(map(str, gci)) + '\n')
        region_all_exp_n50 = compute_n50(region_all_lengths)
        region_all_exp_num_ctg = len(region_all_lengths)
        region_all_gci = []
        for i in range(len(depths_list)):
            region_all_obs_n50 = compute_n50(region_all_obs_length[i])
            region_all_gci.append(round(100 * log2(region_all_obs_n50/region_all_exp_n50 + 1) / log2(region_all_obs_num_ctg[i]/region_all_exp_num_ctg + 1), 4))
        with open(f'{directory}/{prefix}.regions.gci', 'a') as f:
            f.write('----------------------------------------------------------------------------------------------------------------------------------------\n\n\n')
            f.write(f'All_regions\t*\t*\t' + '\t'.join(map(str, region_all_gci)) + '\n')
        print('Computing GCI scores for regions done!!!\n\n')


def sliding_window_average_depth(depths=[], window_size=50000, max_depth=None, start=0, target=None):
    """
    usage: get the averaged depths via sliding window

    input: the single chromosome depth list from the whole-genome depth dictionary generated by filter(),
           the window size in bytes,
           the max depth to plot,
           the position of start,
           the target for plotting

    return: the positions and averaged depths list
    """
    averaged_positions = []
    averaged_depths = []
    window_depths = []
    if len(depths) < window_size:
        print(f'Warning!!! The length ({len(depths)}) of plotting region ({target}:{start}-{start + len(depths)}) is less than the window size ({window_size}), and therefore the window size will be 1 bp', file=sys.stderr)
        window_size = 1
    
    for i, depth in enumerate(depths):
        if depth == 0:
            if len(window_depths) > 0:
                average_depth = sum(window_depths) / len(window_depths)
                if average_depth > max_depth:
                    average_depth = max_depth
                averaged_depths.append(average_depth)
                averaged_positions.append((i+start-1)/1e6)
                window_depths = []
            averaged_depths.append(0)
            averaged_positions.append((i+start)/1e6)
        else:
            window_depths.append(depth)
            if len(window_depths) == window_size:
                average_depth = sum(window_depths) / window_size
                if average_depth > max_depth:
                    average_depth = max_depth
                averaged_depths.append(average_depth)
                averaged_positions.append((i+start)/1e6)
                window_depths = []
    if len(window_depths) > 0:
        average_depth = sum(window_depths) / len(window_depths)
        if average_depth > max_depth:
            average_depth = max_depth
        averaged_depths.append(average_depth)
        averaged_positions.append((i+start)/1e6)
    return averaged_positions, np.array(averaged_depths)


def pre_plot_base(depths_list=[], max_depths=[], window_size=50000, start=0):
    """
    usage: get some prerequisite objects

    input: a list of the whole-genome depth dictionary generated by filter(),
           a list of max depths for the depths_list,
           the window size in chromosome units (0-1) when plotting,
           the position of start

    return: a dictionary keyed by the target with the value positions and averaged depths list,
            fractions of y axis,
            minimum y value,
            max y value
    """
    averaged_dicts = [{} for i in range(len(depths_list))]
    max_averaged_depths_list = [[] for i in range(len(depths_list))]
    
    for target in depths_list[0].keys():
        for i, depthss in enumerate(depths_list):
            depths = depthss[target]
            averaged_positions, averaged_depths = sliding_window_average_depth(depths, window_size, max_depths[i], start, target)
            averaged_dicts[i].update({target:(averaged_positions, averaged_depths)})
            max_averaged_depths_list[i].append(max(averaged_depths))

    y_max = max(max_averaged_depths_list[0]) + 10
    if len(depths_list) == 1:
        y_min = 0
    elif len(depths_list)  == 2:
        y_min = max(max_averaged_depths_list[1]) + 10
    y_frac = y_min / (y_max + y_min)

    return averaged_dicts, y_frac, y_min, y_max


def plot_base(depths_list=[], target=None, averaged_dicts=[], mean_depths=[], y_frac=0, start=0, depth_min=0.1, dist_percent=0.005, y_min=0, y_max=None, image_type='png', directory='.', prefix='GCI', end=None, regions_flag=False, threshold=0):
    """
    usage: the core of plot_depth

    input: a list of the whole-genome depth dictionary generated by filter(),
           the target,
           a list of averaged_dict generated by pre_plot_base(),
           a list of mean depths for the depths_list,
           fractions of y axis,
           the position of start,
           the cutoff in folds of depth to plot,
           the percentage of the distance between the gap intervals in the chromosome,
           minimum y value,
           max y value,
           the format of output images,
           the path to output,
           the prefix of the output gci file,
           the position of end,
           the flag presenting plotting regions or whole genome,
           the threshold of depth

    output: the depth plots for whole genome or specific regions

    """
    depth_colors = ['#2ca25f', '#3C5488']
    flags = [1, -1]

    if len(depths_list) == 1:
        fig, ax = plt.subplots(figsize=(20, 4))
    elif len(depths_list)  == 2:
        fig, ax = plt.subplots(figsize=(20, 8))
        ax.axhline(0, color="black")
        hifi_line = mlines.Line2D([], [], color='#2ca25f', label='HiFi', lw=0.8)
        nano_line = mlines.Line2D([], [], color='#3C5488', label='Nano', lw=0.8)
        legend1 = plt.legend(handles=[hifi_line, nano_line], loc='upper left')
        plt.gca().add_artist(legend1)
    

    blue_flag = False
    red_flag = False
    for i, depthss in enumerate(depths_list):
        depths = depthss[target]
        arguments = (y_frac, 1) if i == 0 else (0, y_frac)
        merged_min_bed = collapse_depth_range({target:depths}, threshold, mean_depths[i] * depth_min, 0, start)
        if len(merged_min_bed[target]) > 0:
            merged_min_bed = merge_merged_depth_bed(merged_min_bed, {target: (end-start)}, dist_percent, start, start, end)
            for segment in merged_min_bed[target]:
                ax.axvspan(segment[0]/1e6, segment[1]/1e6, *arguments, facecolor='#B7DBEA')
            blue_flag = True
        
        merged_0_bed = collapse_depth_range({target:depths}, -1, threshold, 0, start)
        if len(merged_0_bed[target]) > 0:
            merged_0_bed = merge_merged_depth_bed(merged_0_bed, {target: (end-start)}, dist_percent, start, start, end)
            for segment in merged_0_bed[target]:
                ax.axvspan(segment[0]/1e6, segment[1]/1e6, *arguments, facecolor='#FAD7DD')
            red_flag = True
        
        averaged_positions, averaged_depths = averaged_dicts[i][target]
        ax.stackplot(averaged_positions, flags[i] * averaged_depths, lw=0.8, color=depth_colors[i], zorder=4)
        ax.axhline(flags[i] * mean_depths[i], color="r", ls='-.', dash_capstyle='butt', lw=1, zorder=5)


    ax.set_ylim(bottom=-y_min, top=y_max)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    lines = []
    if blue_flag == True:
        merged_min_line = mlines.Line2D([], [], color='#B7DBEA', label=f'The region with the depth in the range of (0, {depth_min}*mean_depth]')
        lines.append(merged_min_line)
    if red_flag == True:
        merged_0_line = mlines.Line2D([], [], color='#FAD7DD', label='The region of zero depth')
        lines.append(merged_0_line)
    mean_line = mlines.Line2D([], [], color="r", ls='-.', dash_capstyle='butt', lw=1, label='Mean Coverage')
    lines.append(mean_line)

    legend2 = plt.legend(handles=lines, loc='lower center', bbox_to_anchor=(0.5, 1), ncols=len(lines))
    plt.gca().add_artist(legend2)

    
    plt.xlabel('Genomic Position (Mb)', fontsize=14)
    plt.ylabel('Depth', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    if regions_flag == False:
        plt.title(f'Filtered depth across the whole genome:{target}', fontsize=18, pad=30)
        plt.tight_layout()
        plt.savefig(f'{directory}/images/{prefix}.{target}.{image_type}', dpi=200)
    else:
        plt.title(f'Filtered depth across the region:{target}:{start}-{end}', fontsize=18, pad=30)
        plt.tight_layout()
        plt.savefig(f'{directory}/images/{prefix}.{target}:{start}-{end}.{image_type}', dpi=200)
    plt.close()


def plot_depth(depths_list=[], depth_min=0.1, depth_max=4.0, window_size=50000, image_type='png', directory='.', prefix='GCI', force=False, targets_length={}, dist_percent=0.005, regions_bed={}, threshold=0):
    """
    usage: plot whole genome depth

    input: a list of the whole-genome depth dictionary generated by filter(),
           the cutoff in folds of depth to plot,
           the max folds of depth to plot,
           the window size in chromosome units (0-1) when plotting,
           the format of output images,
           the path to output,
           the prefix of the output gci file,
           whether to rewrite the existing files (force),
           targets_length generated by the function filter(),
           the percentage of the distance between the gap intervals in the chromosome,
           the regions bed file,
           the threshold of depth

    output: the depth plots for whole genome and specific regions
    """
    if image_type == 'pdf' or image_type == 'png':
        pass
    else:
        sys.exit(f'ERROR!!! The format of output images only supports pdf and png')
    
    
    mean_depths = []
    for depthss in depths_list:
        sum_depths = []
        for depths in depthss.values():
            sum_depths = np.concatenate((sum_depths, depths))
        mean_depths.append(np.mean(sum_depths))
    max_depths = [mean_depth * depth_max for mean_depth in mean_depths]
    
    for target in depths_list[0].keys():
        if os.path.exists(f'{directory}/images/{prefix}.{target}.{image_type}') and force == False:
            sys.exit(f'ERROR!!! The file "{directory}/images/{prefix}.{target}.{image_type}" exists\nPlease use "-f" or "--force" to rewrite')
    print(f'Plotting whole genome depth ...')
    averaged_dicts, y_frac, y_min, y_max = pre_plot_base(depths_list, max_depths, window_size, 0)
    for target in depths_list[0].keys():
        plot_base(depths_list, target, averaged_dicts, mean_depths, y_frac, 0, depth_min, dist_percent, y_min, y_max, image_type, directory, prefix, targets_length[target], False, threshold)
    print(f'Plotting whole genome depth done!!!\n\n')
    
    if len(regions_bed) > 0:
        print(f'Plotting depth for regions ...')
        for target, segments in regions_bed.items():
            for segment in segments:
                start = segment[0]
                end = segment[1]
                if os.path.exists(f'{directory}/images/{prefix}.{target}:{start}-{end}.{image_type}') and force == False:
                    sys.exit(f'ERROR!!! The file "{directory}/images/{prefix}.{target}:{start}-{end}.{image_type}" exists\nPlease use "-f" or "--force" to rewrite')
            
                regions_depths_list = []
                for depthss in depths_list:
                    depths = depthss[target]
                    regions_depths_list.append({target:depths[start:end]})
                averaged_dicts, y_frac, y_min, y_max = pre_plot_base(regions_depths_list, max_depths, window_size, start)
                plot_base(regions_depths_list, target, averaged_dicts, mean_depths, y_frac, start, depth_min, dist_percent, y_min, y_max, image_type, directory, prefix, end, True, threshold)
        print(f'Plotting depth for regions done!!!\n\n')


def GCI(hifi=[], nano=[], directory='.', prefix='GCI', map_qual=30, mq_cutoff=50, iden_percent=0.9, ovlp_percent=0.9, clip_percent=0.1, flank_len=15, threshold=0, plot=False, depth_min=0.1, depth_max=4.0, window_size=50000, image_type='png', force=False, dist_percent=0.005, reference=None, regions=None):
    if directory.endswith('/'):
        directory = '/'.join(directory.split('/')[:-1])
    if os.path.exists(directory):
        if not os.access(directory, os.R_OK):
            sys.exit(f'ERROR!!! The path "{directory}" is unable to read')
        if not os.access(directory, os.W_OK):
            sys.exit(f'ERROR!!! The path "{directory}" is unable to write')
    else:
        os.makedirs(directory)

    if prefix.endswith('/'):
        sys.exit(f'ERROR!!! The prefix "{prefix}" is not allowed')
    
    if plot == True:
        if os.path.exists(f'{directory}/images'):
            if not os.access(f'{directory}/images', os.R_OK):
                sys.exit(f'ERROR!!! The path "{directory}/images" is unable to read')
            if not os.access(f'{directory}/images', os.W_OK):
                sys.exit(f'ERROR!!! The path "{directory}/images" is unable to write')
        else:
            os.makedirs(f'{directory}/images')
        
        image_type = image_type.lower()
    
    ref_refs = []
    for record in SeqIO.parse(reference, 'fasta'):
        ref_refs.append(record.id)
    hifi_bam = []
    hifi_paf = []
    nano_bam = []
    nano_paf = []
    hifi_refs = []
    nano_refs = []
    if hifi != None:
        for file in hifi:
            if file.endswith('.bam'):
                hifi_bam.append(file)
                hifi_samfile = pysam.AlignmentFile(file, 'rb')
                hifi_refs += hifi_samfile.references
                hifi_samfile.close()
            else:
                hifi_paf.append(file)
        if set(hifi_refs) != set(ref_refs):
            sys.exit('ERROR!!! The targets in hifi alignment files are inconsistent with the reference file\nPlease check both hifi alignment files and the reference')
    if nano != None:
        for file in nano:
            if file.endswith('.bam'):
                nano_bam.append(file)
                nano_samfile = pysam.AlignmentFile(file, 'rb')
                nano_refs += nano_samfile.references
                nano_samfile.close()
            else:
                nano_paf.append(file)
        if set(nano_refs) != set(ref_refs):
            sys.exit('ERROR!!! The targets in ont alignment files are inconsistent with the reference file\nPlease check both ont alignment files and the reference')
        

    print('Finding gaps ...')
    Ns_bed, Ns_bed_file = get_Ns_ref(reference, prefix, directory, force)
    if Ns_bed_file != None:
        print(f'Finding gaps done!!! The gaps are in {Ns_bed_file}\n\n')
    else:
        print('Finding gaps done!!! Awesome! No gaps were found!\n\n')

    regions_bed = {}
    if regions != None:
        if os.path.exists(regions) and os.access(regions, os.R_OK):
            with open(regions, 'r') as f:
                for line in f:
                    target, start, end = line.strip().split('\t')
                    if target not in regions_bed.keys():
                        regions_bed[target] = []
                    regions_bed[target].append((int(start), int(end)))
        else:
            sys.exit(f'ERROR!!! "{regions}" is not an available file')
    

    if nano == None:
        depths, targets_length = filter(hifi_paf, hifi_bam, prefix, map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, 'HiFi')
        depths = merge_gaps_depths(depths, Ns_bed)
        merged_depth_bed = merge_depth(depths, prefix, threshold, flank_len, directory, force, 'HiFi')
        compute_index(targets_length, prefix, directory, force, [merged_depth_bed], ['HiFi'], flank_len, dist_percent, regions_bed, [depths], threshold)
        if plot == True:
            plot_depth([depths], depth_min, depth_max, window_size, image_type, directory, prefix, force, targets_length, dist_percent, regions_bed, threshold)

    elif hifi == None:
        depths, targets_length = filter(nano_paf, nano_bam, prefix, map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, 'ONT')
        depths = merge_gaps_depths(depths, Ns_bed)
        merged_depth_bed = merge_depth(depths, prefix, threshold, flank_len, directory, force, 'ONT')
        compute_index(targets_length, prefix, directory, force, [merged_depth_bed], ['Nano'], flank_len, dist_percent, regions_bed, [depths], threshold)
        if plot == True:
            plot_depth([depths], depth_min, depth_max, window_size, image_type, directory, prefix, force, targets_length, dist_percent, regions_bed, threshold)
    
    else:
        if set(hifi_refs) != set(nano_refs):
            sys.exit(f'ERROR!!! The targets in hifi and nano alignment files are inconsistent\nPlease check the reference used in mapping both hifi and ont reads')

        hifi_depths, targets_length = filter(hifi_paf, hifi_bam, prefix+'_hifi', map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, 'HiFi')
        hifi_depths = merge_gaps_depths(hifi_depths, Ns_bed)
        nano_depths, targets_length = filter(nano_paf, nano_bam, prefix+'_nano', map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, flank_len, directory, force, 'ONT')
        nano_depths = merge_gaps_depths(nano_depths, Ns_bed)
        merged_two_type_depths = merge_two_type_depth(hifi_depths, nano_depths, prefix+'_two_type', directory, force)
        merged_two_type_depths = merge_gaps_depths(merged_two_type_depths, Ns_bed)

        hifi_merged_depth_bed = merge_depth(hifi_depths, prefix+'_hifi', threshold, flank_len, directory, force, 'HiFi')
        nano_merged_depth_bed = merge_depth(nano_depths, prefix+'_nano', threshold, flank_len, directory, force, 'ONT')
        two_type_merged_depth_bed = merge_depth(merged_two_type_depths, prefix+'_two_type', threshold, flank_len, directory, force, 'two_types')
        compute_index(targets_length, prefix, directory, force, [hifi_merged_depth_bed, nano_merged_depth_bed, two_type_merged_depth_bed], ['HiFi', 'Nano', 'HiFi + Nano'], flank_len, dist_percent, regions_bed, [hifi_depths, nano_depths, merged_two_type_depths], threshold)
        if plot == True:
            plot_depth([hifi_depths, nano_depths], depth_min, depth_max, window_size, image_type, directory, prefix, force, targets_length, dist_percent, regions_bed, threshold)

    print('GCI finished!!!\nBye!!!')


if __name__=='__main__':
    ###########################
    ### version = 0.4
    ###########################
    version = 'GCI version 0.4'

    parser = argparse.ArgumentParser(prog=sys.argv[0], add_help=False, formatter_class=argparse.RawTextHelpFormatter, description='A program for assessing the T2T genome', epilog='Examples:\npython GCI.py -r ref.fa --hifi hifi.bam hifi.paf ... --nano nano.bam nano.paf ...')

    group_io = parser.add_argument_group("Input/Output")
    group_io.add_argument('-r', '--reference', metavar='FILE', help='The reference file')
    group_io.add_argument('--hifi', nargs='+', metavar='', help='PacBio HiFi reads alignment files (at least one bam file)')
    group_io.add_argument('--nano', nargs='+', metavar='', help='Oxford Nanopore long reads alignment files (at least one bam file)')
    group_io.add_argument('-R', '--regions', metavar='FILE', help='Bed file containing regions')
    group_io.add_argument('-ts', '--threshold', metavar='INT', type=int, help='The threshold of depth to be reported as issues [0]', default=0)
    group_io.add_argument('-dp', '--dist-percent', metavar='FLOAT', type=float, help='The distance between the candidate gap intervals for combining in chromosome units [0.005]', default=0.005)
    group_io.add_argument('-d', dest='directory', metavar='PATH', help='The directory of output files [.]', default='.')
    group_io.add_argument('-o', '--output', dest='prefix', metavar='STR', help='Prefix of output files [GCI]', default='GCI')

    group_fo = parser.add_argument_group("Filter Options")
    group_fo.add_argument('-mq', '--map-qual', metavar='INT', type=int, help='Minium mapping quality for alignments [30]', default=30)
    group_fo.add_argument('--mq-cutoff', metavar='INT', type=int, help='The cutoff of mapping quality for keeping the alignment [50]\n(only used when inputting more than one alignment files)', default=50)
    group_fo.add_argument('-ip', '--iden-percent', metavar='FLOAT', type=float, help='Minimum identity (num_match_res/len_aln) of alignments [0.9]', default=0.9)
    group_fo.add_argument('-op', '--ovlp-percent', metavar='FLOAT', type=float, help='Minimum overlapping percentage of the same read alignment if inputting more than one alignment files [0.9]', default=0.9)
    group_fo.add_argument('-cp', '--clip-percent', metavar='FLOAT', type=float, help='Maximum clipped percentage of the alignment [0.1]', default=0.1)
    group_fo.add_argument('-fl', '--flank-len', metavar='INT', type=int, help='The flanking length of the clipped bases [15]', default=15)

    group_po = parser.add_argument_group("Plot Options")
    group_po.add_argument('-p', '--plot', action='store_const', help='Visualize the finally filtered whole genome (and regions if providing the option `-R`) depth [False]', const=True, default=False)
    group_po.add_argument('-dmin', '--depth-min', metavar='FLOAT', type=float, help='Minimum depth in folds of mean coverage for plotting [0.1]', default=0.1)
    group_po.add_argument('-dmax', '--depth-max', metavar='FLOAT', type=float, help='Maximum depth in folds of mean coverage for plotting [4.0]', default=4.0)
    group_po.add_argument('-ws', '--window-size', metavar='INT', type=int, help='The window size when plotting [50000]', default=50000)
    group_po.add_argument('-it', '--image-type', metavar='STR', help='The format of the output images: png or pdf [png]', default='png')

    group_op = parser.add_argument_group("Other Options")
    group_op.add_argument('-f', '--force', action='store_const', help='Force rewriting of existing files [false]', const=True, default=False)
    group_op.add_argument('-h', '--help', action="help", help="Show this help message and exit")
    group_op.add_argument('-v', '--version', action="version", version=version, help="Show program's version number and exit")

    args = vars(parser.parse_args())
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    if (args['hifi'] == None) and (args['nano'] == None):
        sys.exit('ERROR!!! Please input at least one type of TGS reads alignment files (PacBio HiFi and/or Oxford Nanopore long reads)\nPlease read the help message use "-h" or "--help"')
    
    if args['hifi'] != None:
        bam_num = 0
        for file in args['hifi']:
            if os.path.exists(file) and os.access(file, os.R_OK):
                if file.endswith('.bam'):
                    bam_num += 1
            else:
                sys.exit(f'ERROR!!! "{file}" is not an available file')
        if bam_num == 0:
            sys.exit('ERROR!!! Please input at least one PacBio HiFi reads bam file\nPlease read the help message use "-h" or "--help"')
        
    if args['nano'] != None:
        bam_num = 0
        for file in args['nano']:
            if os.path.exists(file) and os.access(file, os.R_OK):
                if file.endswith('.bam'):
                    bam_num += 1
            else:
                sys.exit(f'ERROR!!! "{file}" is not an available file')
        if bam_num == 0:
            sys.exit('ERROR!!! Please input at least one Oxford Nanopore long reads bam file\nPlease read the help message use "-h" or "--help"')
    
    if args['reference'] == None:
        sys.exit('ERROR!!! Please input the reference file\nPlease read the help message use "-h" or "--help"')
    else:
        if os.path.exists(args['reference']) and os.access(args['reference'], os.R_OK):
            pass
        else:
            sys.exit(f'ERROR!!! \"{args["reference"]}\" is not an available file')

    if args['map_qual'] > args['mq_cutoff']:
        print(f'WARNING!!! The minium mapping quality ({args["map_qual"]}) is higher than the cutoff ({args["mq_cutoff"]}), which means that wouldn\'t filter any reads\nPlease read the help message use "-h" or "--help"', file=sys.stdout)

    print(f'Used arguments:{args}')
    GCI(**args)
    
