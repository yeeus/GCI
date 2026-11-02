import sys
import pysam
import numpy as np
import argparse
import subprocess
import os
import random
import string

# Reuse the partitioning and alignment-scoring primitives from GCI v1.1.
# The repository root must contain gci_partitioned.py and gci_depth_events.py.
from pathlib import Path
from multiprocessing import get_context
import sqlite3
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gci_partitioned import _BucketWriter, _bam_task, _paf_best_target




def _selected_bucket(task):
    """Resolve one Phase 2C read-name partition and retain query identifiers."""
    bucket, root, n_paf, n_bam, targets, map_qual, mq_cutoff, identity_min, ovlp_min = task
    root = Path(root)
    high_quality = set()
    synteny = {}
    paf_maps = [{} for _ in range(n_paf)]
    paf_path = root / f'paf.{bucket:04d}.tsv'
    current_file = 0
    if paf_path.exists():
        with open(paf_path) as src:
            for line in src:
                parts = line.rstrip('\n').split('\t')
                file_id = int(parts[0])
                while current_file < file_id:
                    paf_maps[current_file] = _paf_best_target(synteny)
                    current_file += 1
                paf = parts[1:]
                target = paf[5]
                if target not in targets:
                    continue
                qname = paf[0]
                qlen, qs, qe = int(paf[1]), int(paf[2]), int(paf[3])
                ts, te = int(paf[7]), int(paf[8])
                identity = int(paf[9]) / int(paf[10])
                mapq = int(paf[11])
                if mapq >= map_qual and identity >= identity_min:
                    synteny.setdefault(qname, {}).setdefault(target, []).append(
                        (qlen, qs, qe, ts, te, identity)
                    )
                    if mapq >= mq_cutoff:
                        high_quality.add(qname)
    for file_id in range(current_file, n_paf):
        paf_maps[file_id] = _paf_best_target(synteny)

    del synteny
    bam_maps = []
    for file_id in range(n_bam):
        reads = {}
        bam_path = root / f'bam{file_id}.{bucket:04d}.tsv'
        if bam_path.exists():
            with open(bam_path) as src:
                for line in src:
                    qname, chrom, start, end, qlen, high = line.rstrip('\n').split('\t')
                    reads[qname] = (chrom, int(start), int(end), int(qlen))
                    if high == '1':
                        high_quality.add(qname)
        bam_maps.append(reads)

    files = paf_maps + bam_maps
    if len(files) > 1:
        common = set.intersection(*(set(f) for f in files))
        final = high_quality | common
        selected = {q: v for q, v in files[0].items() if q in final}
        for mapping in files[1:]:
            for q, segment in mapping.items():
                if q in selected:
                    previous = selected[q]
                    if segment[0] == previous[0]:
                        overlap = min(segment[2], previous[2]) - max(segment[1], previous[1])
                        if overlap / segment[-1] < ovlp_min:
                            del selected[q]
                        else:
                            selected[q] = (previous[0], max(segment[1], previous[1]), min(segment[2], previous[2]))
                    else:
                        del selected[q]
                elif q in high_quality:
                    selected[q] = (segment[0], segment[1], segment[2])
    else:
        selected = files[0]

    out_path = root / f'selected.{bucket:04d}.tsv'
    with open(out_path, 'w') as out:
        for qname, segment in selected.items():
            out.write(f'{qname}\t{segment[0]}\t{segment[1]}\t{segment[2]}\n')
    return str(out_path)


def filter(paf_files=[], bam_files=[], prefix=None, map_qual=30, mq_cutoff=50,
           iden_percent=0.9, clip_percent=0.1, ovlp_percent=0.9,
           directory='.', force=False, threads=1, partitions=256):
    """Write filtered BAMs using GCI v1.1 Phase 2C read selection."""
    if threads < 1 or partitions < 1:
        raise ValueError('threads and partitions must both be >= 1')
    output_files_name = [
        os.path.join(directory, f'{prefix[i]}.bam' if prefix is not None
                     else f'{os.path.basename(path)[:-4]}.filter.bam')
        for i, path in enumerate(bam_files)
    ]
    for path in output_files_name:
        if os.path.exists(path) and not force:
            raise SystemExit(f'ERROR!!! The file "{path}" exists\nPlease use "-f" or "--force" to rewrite')

    with pysam.AlignmentFile(bam_files[0], 'rb') as bam:
        targets = dict(zip(bam.references, bam.lengths))

    with tempfile.TemporaryDirectory(prefix='gci-part-', dir=directory) as root:
        with _BucketWriter(root, 'paf', partitions) as writer:
            for i, path in enumerate(paf_files):
                with open(path) as source:
                    for line in source:
                        writer.write(line.split('\t', 1)[0], str(i) + '\t' + line)

        for file_id, bam_path in enumerate(bam_files):
            tasks = []
            for target, length in targets.items():
                step = 1 + (length - 1) // threads
                for left in range(0, length, step):
                    tasks.append((len(tasks), target, left, min(left + step, length),
                                  bam_path, map_qual, clip_percent,
                                  iden_percent, mq_cutoff, root))
            with _BucketWriter(root, f'bam{file_id}', partitions) as writer:
                if tasks:
                    with get_context('fork').Pool(processes=min(threads, len(tasks))) as pool:
                        for path, count in pool.imap(_bam_task, tasks, chunksize=1):
                            with open(path) as source:
                                for line in source:
                                    writer.write(line.split('\t', 1)[0], line)
                            os.remove(path)

        worker_tasks = [
            (bucket, root, len(paf_files), len(bam_files), targets,
             map_qual, mq_cutoff, iden_percent, ovlp_percent)
            for bucket in range(partitions)
        ]
        with sqlite3.connect(Path(root) / 'selected.sqlite') as db:
            db.execute('CREATE TABLE selected (query TEXT PRIMARY KEY, chrom TEXT, '
                       'start INTEGER, end INTEGER, emitted INTEGER DEFAULT -1)')
            with get_context('fork').Pool(processes=min(threads, partitions)) as pool:
                for path in pool.imap_unordered(_selected_bucket, worker_tasks, chunksize=1):
                    with open(path) as source:
                        db.executemany(
                            'INSERT INTO selected (query, chrom, start, end) VALUES (?, ?, ?, ?)',
                            (line.rstrip('\n').split('\t') for line in source),
                        )
                    os.remove(path)
            db.commit()

            for file_id, bam_path in enumerate(bam_files):
                with pysam.AlignmentFile(bam_path, 'rb') as bam:
                    with pysam.AlignmentFile(output_files_name[file_id], 'wb', template=bam) as out:
                        for aln in bam.fetch(until_eof=True):
                            if (aln.is_unmapped or aln.is_secondary or aln.is_supplementary
                                    or aln.mapping_quality < map_qual):
                                continue
                            hit = db.execute(
                                'SELECT chrom, start, end FROM selected WHERE query=? AND emitted!=?',
                                (aln.query_name, file_id),
                            ).fetchone()
                            if (hit is not None and aln.reference_name == hit[0]
                                    and aln.reference_start < hit[2]
                                    and aln.reference_end > hit[1]):
                                out.write(aln)
                                db.execute('UPDATE selected SET emitted=? WHERE query=?',
                                           (file_id, aln.query_name))
                db.commit()

    for path in output_files_name:
        random_string = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
        sorted_path = os.path.join(directory, f'{random_string}.bam')
        subprocess.run(['samtools', 'sort', path, '-o', sorted_path], capture_output=True, check=True)
        subprocess.run(['mv', sorted_path, path], capture_output=True, check=True)
        subprocess.run(['samtools', 'index', path], capture_output=True, check=True)
    return output_files_name



def bamsnap(files=[], output_files=[], reference=None, region=None, regions_file=None, directory='.', prefix='bamsnap', force=False):
	"""
	usage: visualize the alignments

	input: bam file(s),
		   the file names of the filtered bam file(s),
		   the reference file,
		   the region to plot,
		   the regions file,
		   the path to output,
		   the prefix of the output,
		   whether to rewrite the existing files (force)

	output: the plots
	"""

	final_files = []
	for (file, output_file) in zip(files, output_files):
		final_files.append(file)
		final_files.append(output_file)
	
	bam_argument = ' '.join(final_files)
	title_arguments = ' '.join([f'"{file}"' for file in final_files])
	if region != None:
		if os.path.exists(f'{directory}/{prefix}.png') and force == False:
			print(f'ERROR!!! The file "{directory}/{prefix}.png" exists\nPlease using "-f" or "--force" to rewrite', file=sys.stderr)
			raise SystemExit
		subprocess.run(f'bamsnap -bam {bam_argument} -title {title_arguments} -pos {region} -margin 100 -out {directory}/{prefix}.png -ref {reference} -draw coordinates bamplot base -bamplot coverage read -read_color_by strand', shell=True, check=True)
	
	elif regions_file != None:
		if prefix.endswith('/'):
			prefix = prefix.split('/')[0]
		if os.path.exists(f'{directory}/{prefix}'):
			if not os.access(f'{directory}/{prefix}', os.R_OK):
				print(f'ERROR!!! The path "{directory}/{prefix}" is unable to read', file=sys.stderr)
				raise SystemExit
			if not os.access(f'{directory}/{prefix}', os.W_OK):
				print(f'ERROR!!! The path "{directory}/{prefix}" is unable to write', file=sys.stderr)
				raise SystemExit
		else:
			os.makedirs(f'{directory}/{prefix}')
		subprocess.run(f'bamsnap -bam {bam_argument} -title {title_arguments} -bed {regions_file} -margin 100 -out {directory}/{prefix} -save_image_only -ref {reference} -draw coordinates bamplot base -bamplot coverage read -read_color_by strand', shell=True, check=True)


def preprocessing(files=[], directory='.', prefix='bamsnap', threads=1, partitions=256, map_qual=30, mq_cutoff=50, iden_percent=0.9, ovlp_percent=0.9, clip_percent=0.1, plot=False, reference=None, region=None, regions_file=None, force=False):
	if directory.endswith('/'):
		directory = directory.rstrip('/') or '.'
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
		if file.endswith('.bam'):
			bam_files.append(file)
		else:
			paf_files.append(file)
	

	if plot == True:
		if isinstance(prefix, str):
			prefix = [prefix]
		if len(prefix) == 1:
			output_files_name = filter(paf_files, bam_files, None, map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, directory, force, threads, partitions)
		else:
			output_files_name = filter(paf_files, bam_files, prefix[:-1], map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, directory, force, threads, partitions)
		bamsnap(bam_files, output_files_name, reference, region, regions_file, directory, prefix[-1], force)
	else:
		if isinstance(prefix, str):
			output_files_name = filter(paf_files, bam_files, None, map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, directory, force, threads, partitions)
		else:
			output_files_name = filter(paf_files, bam_files, prefix, map_qual, mq_cutoff, iden_percent, clip_percent, ovlp_percent, directory, force, threads, partitions)


if __name__=='__main__':
	parser = argparse.ArgumentParser(prog=sys.argv[0], add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description='This is the core part of the main program GCI.py with additional functions\nYou can provide some bam and paf files (like GCI.py) and set the filtered parameters\nand then you will get the filtered bam file which can be used for subsequent analysis', epilog='Examples:\npython filter_bam.py bam1 paf1 ...')

	group_io = parser.add_argument_group("Input/Output")
	group_io.add_argument('files', nargs='+', metavar='ALIGNMENT-FILE', help='Long reads alignment files (at least one bam file)')
	group_io.add_argument('-d', dest='directory', metavar='PATH', help='The directory of output files [.]', default='.')
	group_io.add_argument('-o', '--output', nargs='*', dest='prefix', metavar='STR', help='Prefix of output files; one prefix corresponds to one bam file in order and if provide the parameter "-p", the last one is used as the prefix for bamsnap outputs (output directory if inputting the regions file) [[$input.filter] for filtered bam files and [bamsnap] for bamsnap outputs]', default='bamsnap')
	group_io.add_argument('-t', '--threads', metavar='INT', type=int, help='Number of worker processes [1]', default=1)
	group_io.add_argument('--partitions', metavar='INT', type=int, help='Number of read-name partitions [256]', default=256)

	group_fo = parser.add_argument_group("Filter Options")
	group_fo.add_argument('-mq', '--map-qual', metavar='INT', type=int, help='Minium mapping quality for alignments [30]', default=30)
	group_fo.add_argument('--mq-cutoff', metavar='INT', type=int, help='The cutoff of mapping quality for keeping the alignment [50]', default=50)
	group_fo.add_argument('-ip', '--iden-percent', metavar='FLOAT', type=float, help='Minimum identity (num_match_res/len_aln) of the reads [0.9]', default=0.9)
	group_fo.add_argument('-op', '--ovlp-percent', metavar='FLOAT', type=float, help='Minimum overlapping percentage of the reads if inputting more than one alignment files [0.9]', default=0.9)
	group_fo.add_argument('-cp', '--clip-percent', metavar='FLOAT', type=float, help='Maximum clipped percentage of the reads [0.1]', default=0.1)

	group_po = parser.add_argument_group("Plot Options")
	group_po.add_argument('-p', '--plot', action='store_const', help='Visualize the filtered bam files', const=True, default=False)
	group_po.add_argument('-ref', '--reference', metavar='FILE', help='The reference file')
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
			if file.endswith('.bam'):
				bam_num += 1
		else:
			print(f'ERROR!!! "{file}" is not an available file', file=sys.stderr)
			raise SystemExit
	if bam_num == 0:
			print('ERROR!!! Please input at least one bam file\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit
	
	
	if args['map_qual'] > args['mq_cutoff']:
		print(f'WARNING!!! The minium mapping quality is {args["map_qual"]} and higher than the cutoff {args["mq_cutoff"]}, which means that wouldn\'t filter any reads\nPlease read the help message using "-h" or "--help"', file=sys.stdout)
	

	if isinstance(args['prefix'], str):
		pass
	elif len(args['prefix']) == 0:
		print('ERROR!!! Please input at least one prefix\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
		raise SystemExit
	elif (args['plot'] == False):
		if len(args['prefix']) != bam_num:
			print(f'ERROR!!! The number of prefixes and bam files is inconsistent\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit
	elif (args['plot'] == True):
		if len(args['prefix']) == 1:
			pass
		elif len(args['prefix']) != (bam_num + 1):
			print(f'ERROR!!! Expect {bam_num + 1} prefixes but provide {len(args["prefix"])}\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit


	if args['plot'] == True :
		if args['reference'] == None:
			print(f'ERROR!!! Please input the reference file\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit
		if (args['region'] == None) and (args['regions_file'] == None):
			print(f'ERROR!!! Please provide the genomic positions (or in bed format)\nPlease read the help message using "-h" or "--help"', file=sys.stderr)
			raise SystemExit
	
	preprocessing(**args)
