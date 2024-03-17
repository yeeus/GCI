import sys
import numpy as np
import argparse
import os
import re
import gzip
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import ScalarFormatter


def get_Ns_ref(reference=None, prefix='GCI', directory='.', force=False):
    """
    usage: get Ns (gaps) of reference and return the bed file (if have)

    input: the reference file

    output: the gaps bed file

    return: the dictionary containing the gaps bed file and the file
    """
    Ns_bed = {}
    pattern = re.compile(r'(?i)N+')
    for record in SeqIO.parse(reference, 'fasta'):
        for match in pattern.finditer(str(record.seq)):
            target = record.id
            if target not in Ns_bed.keys():
                Ns_bed[target] = []
            Ns_bed[target].append((match.start(), match.end()))
    
    if len(Ns_bed) != 0:
        if os.path.exists(f'{directory}/{prefix}.gaps.bed') and force == False:
            sys.exit(f'ERROR!!! The file "{directory}/{prefix}.gaps.bed" exists\nPlease using "-f" or "--force" to rewrite')
        with open(f'{directory}/{prefix}.gaps.bed', 'w') as f:
            for target, segments in Ns_bed.items():
                for segment in segments:
                    f.write(f'{target}\t{segment[0]}\t{segment[1]}\n')
        return Ns_bed, f'{directory}/{prefix}.gaps.bed'
    else:
        return None, None


def parse_depth(depth_file=None, log_reads_type=''):
    """
    usage: parse the depth file and get the whole-genome depth dictionary

    input: the gzipped depth file,
           reads type for logging

    return: the dictionary containing lengths of each chromosome
            one dictionary keyed by the targets with the length value
    """
    print(f'Parsing {log_reads_type} depth file ...')
    depths = {}
    targets_length = {}
    target = ''
    with gzip.open(depth_file, 'rb') as f:
        for line in f:
            item = line.decode('utf-8').strip()
            if item.startswith('>'):
                if target != '':
                    depths[target] = np.array(depths[target])
                    targets_length.update({target:len(depths[target])})
                target = item.split('>')[-1]
                depths[target] = []
            else:
                depths[target].append(int(item))
        depths[target] = np.array(depths[target])
        targets_length.update({target:len(depths[target])})
    print(f'Parsing {log_reads_type} depth file ... done!!!\n\n')
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

    new_merged_depths_bed = {}
    for target, length in targets_length.items():
        new_merged_depths_bed[target] = []
        dist = length * dist_percent
        if start == None and end == None:
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
        plt.savefig(f'{directory}/{prefix}.{target}.{image_type}', dpi=200)
    else:
        plt.title(f'Filtered depth across the region:{target}:{start}-{end}', fontsize=18, pad=30)
        plt.tight_layout()
        plt.savefig(f'{directory}/{prefix}.{target}:{start}-{end}.{image_type}', dpi=200)
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
    
    
    if len(regions_bed) == 0:
        for target in depths_list[0].keys():
            if os.path.exists(f'{directory}/{prefix}.{target}.{image_type}') and force == False:
                sys.exit(f'ERROR!!! The file "{directory}/{prefix}.{target}.{image_type}" exists\nPlease use "-f" or "--force" to rewrite')
        print(f'Plotting whole genome depth ...')
        averaged_dicts, y_frac, y_min, y_max = pre_plot_base(depths_list, max_depths, window_size, 0)
        for target in depths_list[0].keys():
            plot_base(depths_list, target, averaged_dicts, mean_depths, y_frac, 0, depth_min, dist_percent, y_min, y_max, image_type, directory, prefix, targets_length[target], False, threshold)
        print(f'Plotting whole genome depth done!!!\n\n')
    
    else:
        print(f'Plotting depth for regions ...')
        for target, segments in regions_bed.items():
            for segment in segments:
                start = segment[0]
                end = segment[1]
                if os.path.exists(f'{directory}/{prefix}.{target}:{start}-{end}.{image_type}') and force == False:
                    sys.exit(f'ERROR!!! The file "{directory}/{prefix}.{target}:{start}-{end}.{image_type}" exists\nPlease use "-f" or "--force" to rewrite')
            
                regions_depths_list = []
                for depthss in depths_list:
                    depths = depthss[target]
                    regions_depths_list.append({target:depths[start:end]})
                averaged_dicts, y_frac, y_min, y_max = pre_plot_base(regions_depths_list, max_depths, window_size, start)
                plot_base(regions_depths_list, target, averaged_dicts, mean_depths, y_frac, start, depth_min, dist_percent, y_min, y_max, image_type, directory, prefix, end, True, threshold)
        print(f'Plotting depth for regions done!!!\n\n')


def preprocessing(reference=None, hifi=None, nano=None, directory='.', prefix='GCI', depth_min=0.1, depth_max=4.0, window_size=50000, image_type='png', force=False, regions=None, dist_percent=0.005, threshold=0):
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
    
    image_type = image_type.lower()

    print('Finding gaps ...')
    Ns_bed, Ns_bed_file = get_Ns_ref(reference, prefix, directory, force)
    if Ns_bed_file != None:
        print(f'Finding gaps done!!! The gaps are in {Ns_bed_file}\n\n')
    else:
        print('Finding gaps done!!! Awesome! No gaps were found!\n\n')
    
    ref_refs = []
    for record in SeqIO.parse(reference, 'fasta'):
        ref_refs.append(record.id)
    
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
        depths, targets_length = parse_depth(hifi)
        if set(targets_length.keys()) != set(ref_refs):
            sys.exit('ERROR!!! The targets in hifi depth file are inconsistent with the reference file\nPlease check both hifi depth file and the reference')
        depths = merge_gaps_depths(depths, Ns_bed)
        plot_depth([depths], depth_min, depth_max, window_size, image_type, directory, prefix, force, targets_length, dist_percent, regions_bed, threshold)
    elif hifi == None:
        depths, targets_length = parse_depth(nano)
        if set(targets_length.keys()) != set(ref_refs):
            sys.exit('ERROR!!! The targets in ont depth file are inconsistent with the reference file\nPlease check both ont depth file and the reference')
        depths = merge_gaps_depths(depths, Ns_bed)
        plot_depth([depths], depth_min, depth_max, window_size, image_type, directory, prefix, force, Ns_bed)
    else:
        hifi_depths, targets_length = parse_depth(hifi)
        if set(targets_length.keys()) != set(ref_refs):
            sys.exit('ERROR!!! The targets in hifi depth file are inconsistent with the reference file\nPlease check both hifi depth file and the reference')
        hifi_depths = merge_gaps_depths(hifi_depths, Ns_bed)
        nano_depths, targets_length = parse_depth(nano)
        if set(targets_length.keys()) != set(ref_refs):
            sys.exit('ERROR!!! The targets in ont depth file are inconsistent with the reference file\nPlease check both ont depth file and the reference')
        nano_depths = merge_gaps_depths(nano_depths, Ns_bed)
        plot_depth([hifi_depths, nano_depths], depth_min, depth_max, window_size, image_type, directory, prefix, force, Ns_bed, targets_length, dist_percent, regions)


if __name__=='__main__':
    parser = argparse.ArgumentParser(prog=sys.argv[0], add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description='This is the plot function in GCI', epilog='Examples:\npython plot_depth.py -r ref.fa --hifi hifi.depth.gz --nano nano.depth.gz')

    group_io = parser.add_argument_group("Input/Output")
    group_io.add_argument('-r', '--reference', metavar='FILE', help='The reference file')
    group_io.add_argument('--hifi', metavar='FILE', help='The gzipped whole-genome depth file generated by the hifi alignment file')
    group_io.add_argument('--nano', metavar='FILE', help='The gzipped whole-genome depth file generated by the ont alignment file')
    group_io.add_argument('-d', dest='directory', metavar='PATH', help='The directory of output files [.]', default='.')
    group_io.add_argument('-o', '--output', dest='prefix', metavar='STR', help='Prefix of output files [GCI]', default='GCI')

    group_po = parser.add_argument_group("Plot Options")
    group_po.add_argument('-R', '--regions', metavar='FILE', help='Bed file containing regions to plot')
    group_po.add_argument('-ts', '--threshold', metavar='INT', type=int, help='The threshold of depth used in GCI.py [0]', default=0)
    group_po.add_argument('-dmin', '--depth-min', metavar='FLOAT', type=float, help='Minimum depth in folds of mean coverage for plotting [0.1]', default=0.1)
    group_po.add_argument('-dmax', '--depth-max', metavar='FLOAT', type=float, help='Maximum depth in folds of mean coverage for plotting [4.0]', default=4.0)
    group_po.add_argument('-ws', '--window-size', metavar='INT', type=int, help='The window size when plotting [50000]', default=50000)
    group_po.add_argument('-it', '--image-type', metavar='STR', help='The format of the output images: png or pdf [png]', default='png')


    group_op = parser.add_argument_group("Other Options")
    group_op.add_argument('-f', '--force', action='store_const', help='Force rewriting of existing files', const=True, default=False)
    group_op.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    args = vars(parser.parse_args())
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    if (args['hifi'] == None) and (args['nano'] == None):
        sys.exit('ERROR!!! Please input at least one depth file\nPlease read the help message using "-h" or "--help"')
    
    if (args['hifi'] != None):
        if os.path.exists(args['hifi']) and os.access(args['hifi'], os.R_OK):
            pass
        else:
            sys.exit(f'ERROR!!! \"{args["hifi"]}\" is not an available file')
    if (args['nano'] != None):
        if os.path.exists(args['nano']) and os.access(args['nano'], os.R_OK):
            pass
        else:
            sys.exit(f'ERROR!!! \"{args["nano"]}\" is not an available file')
    

    if args['reference'] == None:
        sys.exit('ERROR!!! Please input the reference file\nPlease read the help message using "-h" or "--help"')
    else:
        if os.path.exists(args['reference']) and os.access(args['reference'], os.R_OK):
            pass
        else:
            sys.exit(f'ERROR!!! \"{args["reference"]}\" is not an available file')
    
    print(f'Used arguments:{args}')
    preprocessing(**args)
