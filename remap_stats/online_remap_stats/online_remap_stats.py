"""
This provides toil job functions for output from --remap_stats option in the aligner.
"""
from toil.job import Job
from toil.common import Toil #todo: necessary? added later

from cigar import Cigar
from Bio import SeqIO
import statistics as stat
import collections as col

def save_input_assembly_stats(job, all_contig_lengths, remap_stats_internal_file, options):
    with open(job.fileStore.readGlobalFile(remap_stats_internal_file), "a+") as f:
        f.write("sequences included in each assembly file: \n")
        for assembly, seq_lens in all_contig_lengths.items():
            total_seq_len = int()
            for seq_len in seq_lens.values():
                total_seq_len += seq_len
            f.write(assembly + "\t" + str(total_seq_len) + "\n")
        
        if options.remap_stats_raw:
            print("dictionary of all input sequence lengths: " + str(total_seq_len) + "\n")

        f.write("\n")

def save_sequence_remapped_stats(job, assembly_to_align_file, poor_mapping_sequence_file, remap_stats_internal_file, options):
    # with open(job.fileStore.readGlobalFile(poor_mapping_sequence_file), "w+") as inf:
    with open(job.fileStore.readGlobalFile(remap_stats_internal_file), "a+") as outf:
        outf.write("stats for assembly file " + assembly_to_align_file + "\n")
        total_seq_len = int()
        remap_seq_lengths = dict()
        remap_seqs = SeqIO.index(job.fileStore.readGlobalFile(poor_mapping_sequence_file), "fasta")
        for remap_seq_name, seq in remap_seqs.items():
            total_seq_len += len(seq)
            remap_seq_lengths[remap_seq_name] = len(seq)
        outf.write("number of sequences sent to remapping: " + str(len(remap_seqs)) + "\n")
        outf.write("total length of sequences sent to remapping: " + str(total_seq_len) + "\n")
        outf.write("average length of sequences sent to remapping: " + str(total_seq_len/len(remap_seqs)) + "\n")
        outf.write("mode length of sequences sent to remapping: " + str(stat.mode(remap_seq_lengths.values())) + "\n")
        if options.remap_stats_raw:
            outf.write("dict of lengths of sequences sent to remapping:\n")
            job.addFollowOnJobFn(save_raw_data, remap_seq_lengths)
        outf.write("\n")

def save_all_to_all_mappings_stats(job, all_to_all_mapping_file, remap_stats_internal_file, options):
    """
    
    #todo: add options.remap_stats raw output
    #todo: make fxn generalized to any given mapping file? Or make a near-duplicate, with dif annotation.
    """
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ", all_to_all_mapping_file)
    
    with open(job.fileStore.readGlobalFile(remap_stats_internal_file), "a+") as outf:
        outf.write("all_to_all mapping stats:\n")
        
    mapping_coverage_points_job = job.addFollowOnJobFn(get_all_to_all_mapping_coverage_points, all_to_all_mapping_file, options)
    mapping_coverage_points = mapping_coverage_points_job.rv()
    
    mapping_coverage_coords_job = mapping_coverage_points_job.addFollowOnJobFn(get_mapping_coverage_coords, mapping_coverage_points)
    mapping_coverage_coords = mapping_coverage_coords_job.rv()

    mapping_coverage_lengths_job = mapping_coverage_points_job.addFollowOnJobFn(get_mapping_coverage_lengths, mapping_coverage_coords)
    mapping_coverage_lengths = mapping_coverage_lengths_job.rv()

    get_mapping_coverage_stats_job = mapping_coverage_coords_job.addFollowOnJobFn(get_mapping_coverage_stats, mapping_coverage_coords)
    mapping_coverage_stats = get_mapping_coverage_stats_job.rv()

    get_mapping_coverage_stats_job.addFollowOnJobFn(write_mapping_stats_to_file, mapping_coverage_stats, remap_stats_internal_file, options)

def write_mapping_stats_to_file(job, mapping_coverage_stats, remap_stats_internal_file, options):
    [num_mappings, bases_mapped, avg, mode] = mapping_coverage_stats
    with open(job.fileStore.readGlobalFile(remap_stats_internal_file), "a+") as outf:
        # outf.write("stats for mapping file " + mapping_file + ":\n")
        outf.write("total number of mappings: " + num_mappings + ":\n")
        outf.write("total bases mapped: " + bases_mapped + ":\n")
        outf.write("average mapping length: " + avg + ":\n")
        outf.write("mode mapping length: " + mode + ":\n")
        
        if options.remap_stats_raw:
            outf.write("dict of lengths of sequences sent to remapping:\n")
            #todo: check to see if opening filestream in two places is valid.
            job.addFollowOnJobFn(save_raw_data, remap_stats_internal_file, mapping_coverage_coords)

        outf.write("\n")

def get_all_to_all_mapping_coverage_points(job, map_to_ref_file, options):
    """
    Returns:
    all start and stop points of mappings, sorted by contigs.
        key: tuple(fasta_file, contig_id), value: list[regions in tuple(point_value, start_bool) format].
        where start_bool is true if the point is a start of a region, and false if the point is a stop of the region.
    """
    # record all regions of each contig that map well, broken down into the 
    # start-points and stop-points of that region. 
    # key: tuple(fasta_file, contig_id), value: list[regions in tuple(point_value, start_bool) format].
    mapping_coverage_points = col.defaultdict(list)

    # add all start and end points for regions that map well 
    with open(job.fileStore.readGlobalFile(map_to_ref_file)) as f:
        for mapping in f:
            # skip the header files:
            if mapping[0] == "@":
                continue

            # parse line in map_file:
            mapping = mapping.split("\t")
            
            cig_str = mapping[5]

            query_contig_id = mapping[0]
            query_start = get_query_start(cig_str)
            query_stop = get_query_stop(cig_str)

            ref_contig_id = mapping[2]
            ref_start = mapping[3]
            ref_stop = get_target_alignment_length(cig_str)

            # add these coordinates to mapping_coverage_points
            mapping_coverage_points[query_contig_id].append((query_start, True))
            mapping_coverage_points[query_contig_id].append((query_stop, False))
            mapping_coverage_points[ref_contig_id].append((ref_start, True))
            mapping_coverage_points[ref_contig_id].append((ref_stop, False))

    return mapping_coverage_points

def get_query_start(cigar_string):
    """
    Reads the cigar string, returns how many bases from index 0 the alignment begins.
    """
    cig = Cigar(cigar_string)

    cig_list = list(cig.items())

    if cig_list[0][1] in ["H", "S"]:
        # then the mapping has a clipping. Return clipping length
        # print("start", cig_list[0][0])
        return str(cig_list[0][0])
    else:
        # print("start", 0)
        return str(0)
    

def get_query_stop(cigar_string):
    """
    Reads the cigar string, returns how many bases from index 0 of the original read
        the alignment's stop is.
    """
    cig = Cigar(cigar_string)
    
    cig_list = list(cig.items())

    # modify len(cig) to get the length of the contig up to the point of the stop:
    stop_position = len(cig)
    if cig_list[0][1] == "H":
        # then the length of the hard clipping at the beginning won't be included in the
        # contig. Manually include it:
        stop_position += cig_list[0][0]
    if cig_list[-1][1] == "S":
        # then the length of the soft clipping after the mapping will be included in the
        # contig. Manually remove it:
        stop_position -= cig_list[-1][0]
    return str(stop_position)

def get_target_alignment_length(cigar_string):
    """
    Determines how much of the target was aligned to the query, based on the query's 
    cig_str. This is equal to the sum of matches/substitutions ("M") and deletions ("D")
    in the cig_str. Insertions ("I") don't count because these are bases contained in the 
    query but not the target. 
    """
    cig = Cigar(cigar_string)
    cig_list = list(cig.items())

    target_alignment_length = 0
    for i in cig_list:
        if i[1] in ["M", "D"]:
            target_alignment_length += i[0]
    return target_alignment_length

def get_mapping_coverage_coords(job, mapping_coverage_points):
    """
    Returns all the coords (defined by tuple(start,stop)) that are covered by at least one mapping in 
    mapping_coverage_points.
    """
    # mapping_coverage_coords is key: contig_id, value: list of coords: [(start, stop)]
    mapping_coverage_coords = col.defaultdict(list)
    for contig_id in mapping_coverage_points:
        contig_coverage_points = sorted(mapping_coverage_points[contig_id], key=operator.itemgetter(0, 1))
        # contig_coverage_points = sorted(mapping_coverage_points[contig_id], key=lambda point: point[0])
        open_points = 0
        current_region = [0, 0] # format (start, stop)
        for i in range(len(contig_coverage_points)):
            if open_points:
                # then we have at least one read overlapping this region.
                # expand the stop point of current_region
                current_region[1] = contig_coverage_points[i][0]
            if contig_coverage_points[i][1]:
                # if start_bool is true, the point represents a start of mapping
                open_points += 1
                if open_points == 1:
                    # that is, if we've found the starting point of a new current_region,
                    # so we should set the start of the current_region.
                    current_region[0] = contig_coverage_points[i][0]
            else:
                # if start_bool is not true, the point represents the end of a mapping.
                open_points -= 1
                if not open_points:
                    # if there's no more open_points in this region, then this is the 
                    # end of the current_region. Save current_region.
                    mapping_coverage_coords[contig_id].append(current_region.copy())
    return mapping_coverage_coords

def get_mapping_coverage_lengths(job, mapping_coverage_coords):
    lengths = list()
    for coord in mapping_coverage_coords.values():
        lengths.append(coord[1] - coord[0])
    return lengths

    
def get_mapping_coverage_stats(job, mapping_coverage_lengths):
    num_mappings = len(mapping_coverage_lengths)
    total = sum(mapping_coverage_lengths)
    avg = total/num_mappings
    mode = stat.mode(mapping_coverage_lengths)

    return num_mappings, total, avg, mode

def save_raw_data(job, remap_stats_internal_file, raw_data):
    with open(job.fileStore.readGlobalFile(remap_stats_internal_file), "a+") as outf:
        outf.write(raw_data)
        outf.write("\n")