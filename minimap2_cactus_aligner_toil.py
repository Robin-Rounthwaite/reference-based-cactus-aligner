
"""
Given a directory of asseblies (.fa) for alignment to a reference genome fasta,
aligns all the assemblies to the genome, and then excludes all sequences of size 
>=minimum_size_remap that either don't have alignments or have an alignment with mapq 
<mapq_cutoff. These sequences that fall below the mapQ cutoff are then realigned to 
all non-reference sequences that aren't from their input assembly set of contigs.

Note: "_segment_" keyword in the contig ids is used internally. Adding that to your contig
 ids before running this pipeline may result in problems.

#TODO: use mappy to streamline calls to minimap2, to avoid unnecessary io.
https://pypi.org/project/mappy/2.2rc0/
"""

from argparse import ArgumentParser

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

import os
import subprocess
from Bio import SeqIO
import cigar
import collections as col
import operator

import shutil

from sam_to_lastz_cigar import make_lastz_output

def empty(job):
    """
    An empty job, for easier toil job organization.
    """
    return

def rename_duplicate_contig_ids(job, assembly_files):
    """
    Sometimes, when combining assemblies from multiple sources, multiple contigs get the 
    same name. This function slightly modifies all but one of the contigs with the same
    name to ensure that there are no duplicates. Renamed contigs are in a format that
    should be easy to reverse.

    Given a list of assembly files (in job.fileStore.writeGlobalFile id output),
    outputs the list of edited assembly files, with the only 
    difference being that all the contigs have been given unique names. Unique names 
    follow this formula:
    x = original contig id
    y = unique_integer
    new id = x_renamed_y
    """
    
    contig_ids = set()
    unique_id = int()

    for assembly in assembly_files:
        assembly_file = job.fileStore.readGlobalFile(assembly)
        assembly_contigs = SeqIO.parse(assembly_file, "fasta")
        output_contigs = list()
        
        for contig in assembly_contigs:
            
            if contig.id in contig_ids:
                old_id = contig.id
                
                while contig.id in contig_ids:
                    # then there is a duplicate contig_id. edit this one.
                    # keep changing the contig_id until we get a completely unique id.
                    contig.id = old_id + "_renamed_" + str(unique_id)
                    contig.description = old_id + "_renamed_" + str(unique_id)
                    unique_id += 1
                    
                #record the new contig id as an observed id.
                contig_ids.add(contig.id)
                
            else:
                # this isn't a duplicate contig_id. record it.
                contig_ids.add(contig.id)

            output_contigs.append(contig)

        # write the altered assembly to original file.    
        SeqIO.write(output_contigs, assembly_file, "fasta")

    # return the file ids of the original, overwritten assembly files.
    return assembly_files

def align_all_assemblies(job, reference_file, assembly_files, options):
    """
    The head job for the main function of this workflow. Iterates through each assembly in
    assembly_files, and calls align_assembly on it. This aligns it to the reference, 
    then aligns poorly-mapping regions to the rest of assembly_files.
    """
    mapping_files = list()
    
    head_job = job.addChildJobFn(empty)
    head_job.rv()

    for assembly_to_align_file in assembly_files:
        assembly_mapping_file = head_job.addFollowOnJobFn(align_assembly, reference_file, assembly_files, assembly_to_align_file, options).rv()
        mapping_files.append(assembly_mapping_file)

    return job.addFollowOnJobFn(consolidate_mapping_files, mapping_files).rv()


    
def align_assembly(job, reference_file, assembly_files, assembly_to_align_file, options):
    """
    Aligns assembly_to_align_file to the reference, then aligns poorly-mapping regions to the rest of assembly_files.
    """
    # mapping_files compiles all the mappings to reference and mappings between assemblies.
    mapping_files = list()
    
    ## map to reference phase:
    # map assembly to reference. Get the id of the map-to-ref file.
    map_to_ref_job = job.addChildJobFn(map_assembly_to_ref, assembly_to_align_file, reference_file)
    map_to_ref_file = map_to_ref_job.rv()
    mapping_files.append(map_to_ref_file)
    # map_to_ref_job.addChildJobFn(debug_get_contig_mappings, map_to_ref_file, 20)

    ## re-mapping phase, from regions of the assemblies that map poorly to the reference, 
    ## to all the assembly sequence.
    # extract all the start & stop regions of good mapping regions in the assemblies
    # (this doesn't yet include options.sequence_context, but does include mapq cutoff).
    mapping_coverage_points_job = map_to_ref_job.addFollowOnJobFn(get_mapping_coverage_points, map_to_ref_file, options)
    mapping_coverage_points = mapping_coverage_points_job.rv()

    # assemble mapping_coverage_points into a list of proper coverage coordinates in 
    # (start_of_coverage, stop_of_coverage) format.
    mapping_coverage_coords_job = mapping_coverage_points_job.addFollowOnJobFn(get_mapping_coverage_coordinates, mapping_coverage_points)
    mapping_coverage_coords = mapping_coverage_coords_job.rv()

    # find the coordinates of regions that don't have mapping coverage. These are 
    # organized in a dictionary with
    # dict{contig_id: list[tuple(start_of_poor_coverage, stop_of_poor_coverage)]} format.
    # this includes options.sequence_context, and options.minimum_size_remap.
    poor_mapping_coverage_coords_job = mapping_coverage_coords_job.addFollowOnJobFn(get_poor_mapping_coverage_coordinates, assembly_to_align_file, mapping_coverage_coords, options)
    poor_mapping_coverage_coords = poor_mapping_coverage_coords_job.rv()

    # extract the actual sequence that has poor mapping coverage from the contigs.
    # saved in a fasta file in the filestore.
    poor_mapping_sequence_file_job = poor_mapping_coverage_coords_job.addFollowOnJobFn(get_poor_mapping_sequences, assembly_to_align_file, poor_mapping_coverage_coords)
    poor_mapping_sequence_file = poor_mapping_sequence_file_job.rv()

    # # map the poor mapping sequence to all the other assemblies!
    map_to_assemblies_file_job = poor_mapping_sequence_file_job.addFollowOnJobFn(remap_poor_mapping_sequences, poor_mapping_sequence_file, assembly_to_align_file, assembly_files)
    mapping_files.append(map_to_assemblies_file_job.rv())
    # map_to_assemblies_file_job.addChildJobFn(debug_get_contig_mappings, mapping_files[-1], 20)

    # consolidate all the mapping_files to become a single file.
    consolidate_mapping_files_job = map_to_assemblies_file_job.addFollowOnJobFn(consolidate_mapping_files, mapping_files)
    consolidated_mapping_files = consolidate_mapping_files_job.rv()

    return consolidate_mapping_files_job.addFollowOnJobFn(relocate_remapped_fragments_to_source_contigs, consolidated_mapping_files, assembly_to_align_file).rv()

    # return poor_mapping_sequence_file_job.addFollowOnJobFn(consolidate_mapping_files, mapping_files).rv()
    # return map_to_assemblies_file_job.addFollowOnJobFn(consolidate_mapping_files, mapping_files).rv()

def map_assembly_to_ref(job, assembly_to_align_file, reference_file):
    map_to_ref_file = job.fileStore.getLocalTempFile()
    print("*************************************about to call minimap to ref")
    subprocess.call(["minimap2", "-ax", "asm5", "-o",
                    map_to_ref_file, job.fileStore.readGlobalFile(reference_file), job.fileStore.readGlobalFile(assembly_to_align_file)])
    print("*************************************about to call minimap to ref")
    return job.fileStore.writeGlobalFile(map_to_ref_file)

def get_mapping_coverage_points(job, map_to_ref_file, options):
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
            
            # note: to keep contig_ids unique between fasta files, pair with the
            # original fasta filename of the assembly to the contig name.
            # contig_id = (fasta_file_name, mapping[0])
            contig_id = mapping[0]
            cig_str = mapping[5]
            if int(mapping[4]) >= options.mapq_cutoff:
                # then we have an mapping that counts as a high mapq.

                # find out start coord for high mapq coverage:
                start = get_start_of_cigar(cig_str)
                # find out stop coord for high mapq coverage:
                stop = get_stop_of_cigar(cig_str)

                # add these coordinates to mapping_coverage_points
                mapping_coverage_points[contig_id].append((start, True))
                mapping_coverage_points[contig_id].append((stop, False))
    # job.log("------------------------------------------------------len(mapping_coverage_points): " + str(len(mapping_coverage_points["624"])))

    return mapping_coverage_points

def get_start_of_cigar(cigar_string):
    """
    Reads the cigar string, returns how many bases from index 0 the alignment begins.
    """
    cig = cigar.Cigar(cigar_string)

    #todo: make it so that we can check for bad cigs:
    # if cig == ["not a real cig"]:
    #     #then we've been given an empty cigar string, which shouldn't ever happen.
    #     print("WARNING: get_start_of_cigar has been requested of an empty cigar string. "
    #             + "This shouldn't happen.")
    #     return
    # print("cigar_string", cigar_string)
    cig_list = list(cig.items())

    if cig_list[0][1] in ["H", "S"]:
        # then the mapping has a clipping. Return clipping length
        # print("start", cig_list[0][0])
        return cig_list[0][0]
    else:
        # print("start", 0)
        return 0
    

def get_stop_of_cigar(cigar_string):
    """
    Reads the cigar string, returns how many bases from index 0 of the original read
        the alignment's stop is.
    """
    # print("cigar_string", cigar_string)

    cig = cigar.Cigar(cigar_string)
    
    #todo: make it so that we can check for bad cigs:
    # if cig == ["not a real cig"]:
        # #then we've been given an empty cigar string, which shouldn't ever happen.
        # print("WARNING: get_start_of_cigar has been requested of an empty cigar string. "
        #         + "This shouldn't happen.")
        # return
    
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
    # print("stop", stop_position)
    return stop_position

def get_mapping_coverage_coordinates(job, mapping_coverage_points):
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

def get_poor_mapping_coverage_coordinates(job, assembly_file, mapping_coverage_coords, options):
    """
    mapping_coverage_coords is a dictionary of lists of coords in (start, stop) format.
    This function returns poor mapping coords, which is essentially the gaps between 
        those coords.
    example: mapping_coverage_coords{contig_1:[(3,5), (7, 9)]} would result in
                mapping_coverage_coords{contig_1:[(0,3), (5,7), (9, 11)]}, if contig_1 had a
                length of 11.
    variables:
        contig_lengths: A dictionary of the length of all the contigs in 
            {key: contig_id value: len(contig)} format.
        mapping_coverage_coords: a dictionary of lists of coords in 
            {key: contig_id, value:[(start, stop)]}
        sequence_context: an integer, representing the amount of sequence you would 
            want to expand each of the poor_mapping_coords by, to include context
            sequence for the poor mapping sequence. 
    """
    #TODO: moving contig_lenghts calculation to outer fxn so I can distribute across two fxns
    # get the length of the contigs in assembly_file:
    contig_lengths = directly_calculate_contig_lengths(job.fileStore.readGlobalFile(assembly_file))
    
    # poor_mapping_coords has key: contig_id, value list(tuple_of_positions(start, stop))
    poor_mapping_coords = col.defaultdict(list)
    for contig_id in contig_lengths:
        if contig_id in mapping_coverage_coords:
            if mapping_coverage_coords[contig_id][0][0] > 0:
                # if the first mapping region for the contig doesn't start at the start of
                # the contig, the first region is between the start of the contig and the 
                # start of the good_mapping_region.
                poor_mapping_stop = mapping_coverage_coords[contig_id][0][0] + options.sequence_context
                if poor_mapping_stop > contig_lengths[contig_id]:
                    poor_mapping_stop = contig_lengths[contig_id]
                if poor_mapping_stop - 0 >= options.minimum_size_remap: # implement size threshold.
                    poor_mapping_coords[contig_id].append((0, poor_mapping_stop))
                    # print("#################################################included coords for remap:", 0, poor_mapping_stop, poor_mapping_stop - 0, (poor_mapping_stop - 0)<100)
                else:
                    pass
                    # print("1________________Blocked:", 0, poor_mapping_stop)
            for i in range(len(mapping_coverage_coords[contig_id]) - 1):
                # for every pair of mapping coords i and i + 1,
                # make a pair of (stop_from_ith_region, start_from_i+1th_region) to
                # represent the poor_mapping_coords. Include sequence_context as necessary.
                poor_mapping_start = mapping_coverage_coords[contig_id][i][1] - options.sequence_context
                if poor_mapping_start < 0:
                    poor_mapping_start = 0
                    
                poor_mapping_stop = mapping_coverage_coords[contig_id][i + 1][0] + options.sequence_context
                if poor_mapping_stop > contig_lengths[contig_id]:
                    poor_mapping_stop = contig_lengths[contig_id]

                if poor_mapping_stop - poor_mapping_start >= options.minimum_size_remap: # implement size threshold.
                    poor_mapping_coords[contig_id].append((poor_mapping_start, poor_mapping_stop))
                else:
                    pass
                    # print("2________________Blocked:", poor_mapping_start, poor_mapping_stop)
            if mapping_coverage_coords[contig_id][-1][1] < contig_lengths[contig_id]:
                # if the last mapping region for the contig stops before the end of
                # the contig, the last region is between the end of the mapping and the 
                # end of the contig.
                poor_mapping_start = mapping_coverage_coords[contig_id][-1][1] - options.sequence_context
                if poor_mapping_start < 0:
                    poor_mapping_start = 0
                if contig_lengths[contig_id] - poor_mapping_start >= options.minimum_size_remap: # implement size threshold.
                    poor_mapping_coords[contig_id].append((poor_mapping_start, contig_lengths[contig_id]))
                else:
                    pass
                    # print("3________________Blocked:", poor_mapping_start, contig_lengths[contig_id])

        else:
            # there isn't a good_mapping region for this contig. The full length of 
            # the contig belongs in poor_mapping_coords.
            poor_mapping_coords[contig_id].append((0, contig_lengths[contig_id]))
    return poor_mapping_coords

def directly_calculate_contig_lengths(assembly_file):
    contig_lengths = dict()
    contigs = SeqIO.index(assembly_file, "fasta")
    for contig_name, seq in contigs.items():
        contig_lengths[contig_name] = len(seq)
    return contig_lengths

def get_poor_mapping_sequences(job, assembly_file, poor_mapping_coords):
    """
    ---Sequence extraction:---
    Read in the entire fasta file.
    
    for each contig-subdivided list of (start, stop) coordinates in 
    self.files.loc["fasta_files"], extract the sequence associated with 
    each coordinate in that contig.
    
    """
    
    # make fasta file for later remapping all_to_all.
    poor_mapping_sequence_file = job.fileStore.getLocalTempFile()

    contigs = SeqIO.index(job.fileStore.readGlobalFile(assembly_file), "fasta")
    with open(poor_mapping_sequence_file, "w+") as outf:
        for contig_name, contig_record in contigs.items():
            for coord in poor_mapping_coords[contig_name]:
                # for each coord in the low_mapq_coords corresponding to a specific contig:
                # extract the sequence for that contig.
                low_mapq_sequence = contig_record.seq[coord[0]: coord[1]]
                outf.write(
                    ">" + contig_name + "_segment_start_" + str(coord[0]) + "_stop_" + str(coord[1]) + "\n")
                outf.write(str(low_mapq_sequence) + "\n")
            
    return job.fileStore.writeGlobalFile(poor_mapping_sequence_file)

def remap_poor_mapping_sequences(job, poor_mapping_sequence_file, assembly_to_align_file, assembly_files):
    """
    ---Minimap2 all-to-all alignment:---
    input: poor-mapQ only fasta files,
    output: minimap2 all-to-all alignments of all low_mapq segments 
        to all the fasta_files.
    """
    minimap_calls = int()
    all_assemblies_but_to_align = assembly_files.copy()
    #todo: remove below line!! We want mappings between contigs in the same file.
    # all_assemblies_but_to_align.remove(assembly_to_align_file)

    remapping_files = list()
    for target_mapping_file in all_assemblies_but_to_align:
        output_file = job.fileStore.getLocalTempFile()
        output_file_global = job.fileStore.writeGlobalFile(output_file)
        # for every target mapping file (but check to make sure the target isn't 
        # the same original fasta as the poor_mapping_sequence_file),
        
        # map low_mapq_file to target_fasta_file.
        subprocess.call(["minimap2", "-ax", "asm5", "-o",
                            job.fileStore.readGlobalFile(output_file_global), job.fileStore.readGlobalFile(target_mapping_file), job.fileStore.readGlobalFile(poor_mapping_sequence_file)])
        minimap_calls += 1
        remapping_files.append(output_file_global)

    # job.log("minimap call count: " + str(minimap_calls))
    return job.addChildJobFn(consolidate_mapping_files, remapping_files).rv()

def consolidate_mapping_files(job, mapping_files):
    """
    Warning: discards headers of all mapping files.
    Given a list of mapping files, consolidates the contents (not counting headers) into a
    single file.
    """
    consolidated_mappings = job.fileStore.getLocalTempFile()
    with open(consolidated_mappings,"w") as outfile:
        for x in mapping_files:
            with open(job.fileStore.readGlobalFile(x)) as f1:
                line_cnt = 0
                for line in f1:
                    if not line.startswith("@"):
                        outfile.write(line)
                    line_cnt += 1
    return job.fileStore.writeGlobalFile(consolidated_mappings)

def relocate_remapped_fragments_to_source_contigs(job, mapping_file, fasta_file):
    """
    renames contig fragments from the remapping step to the original
    contig name. Changes the cigar clipping based on the fragment's coordinates too. 
    Arguments:
        mapping_file {[type]} -- [description]
        fasta_file {[type]} -- [description]
    """
    
    modified_mapping_file = job.fileStore.getLocalTempFile()
    contig_lengths  = directly_calculate_contig_lengths(job.fileStore.readGlobalFile(fasta_file))

    with open(job.fileStore.readGlobalFile(mapping_file)) as inf:
        with open(modified_mapping_file, "w") as outf:
            for line in inf:
                parsed = line.split()
                if (int(parsed[1])%8)//4==1:
                    # if the line is flagged as unmapped, just write it to the outfile.
                    outf.write(line)
                elif "segment" in parsed[0]:
                    # construct the new name by dropping all the parts of the old name from "_segment_" onwards.
                    name_parsed = parsed[0].split("_segment_")
                    print("---------------------------------------------name_parsed:", name_parsed)
                    new_name = name_parsed[0]

                    # extract the start and stop of the 
                    seg_start = name_parsed[1].split("_")[-3]
                    seg_stop = name_parsed[1].split("_")[-1]
                    cig_str = parsed[5]
                    cig = cigar.Cigar(cig_str)
                    cig_list = list(cig.items())
                    #cig_temp_start is used to determine where the alignment starts *with relation
                    # to the start of the contig fragment*. This is useful for calculating
                    # cig_stop, below.

                    alignment_start_pos_in_seg = int()
                    if cig_list[0][1] in ["H", "S"]:
                        alignment_start_pos_in_seg = int(cig_list[0][0])
                    else:
                        alignment_start_pos_in_seg = 0

                    alignment_length = len(cig)
                    if cig_list[0][1] == "S":
                        alignment_length = alignment_length - cig_list[0][0]
                    if cig_list[-1][1] == "S":
                        alignment_length = alignment_length - cig_list[-1][0]

                    alignment_end_pos_in_seg = alignment_start_pos_in_seg + alignment_length

                    ## modify cig_start (the clipping at the beginning of the cigar)
                    if cig_list[0][1] in ["H", "S"]:
                        # then the mapping has a clipping at start. Modify cig_str clipping.
                        cig_list[0] = (int(cig_list[0][0]) + int(seg_start), cig_list[0][1])
                    elif seg_start:
                        # then there is a nonzero seg_start, but there isn't a clipping for cig_str. Add clipping.
                        cig_list.insert(0, (seg_start, "H"))

                    ## modify clipping at the end of cig
                    # print("alignment_end_pos_in_seg", alignment_end_pos_in_seg, "seg_start", seg_start, "contig_lengths[new_name]", contig_lengths[new_name])
                    if alignment_end_pos_in_seg + int(seg_start) < contig_lengths[new_name]:
                        # then there is additional clipping that needs to be added to the end of the cig.
                        end_clipping_len = contig_lengths[new_name] - (alignment_end_pos_in_seg + int(seg_start))
                        if cig_list[-1][1] in ["H", "S"]:
                            # then the mapping has a clipping at the end. Modify the clipping.
                            #TODO: make it so you properly remove end clipping, not first two chars!
                            cig_list[-1] = (end_clipping_len, cig_list[-1][1])
                        else:
                            #the mapping doesn't have a clipping. Add one.
                            cig_list.append((end_clipping_len, "H"))
                    
                    # compose the new cig
                    new_cig = ""
                    for tup in cig_list:
                        new_cig += str(tup[0]) + tup[1]

                    #now, alter the line
                    new_line = new_name + "\t" + "\t".join(parsed[1:5]) + "\t" + new_cig + "\t" + "\t".join(parsed[6:])
                    # print("line after modification: ", new_line)

                    #add it to the outfile
                    outf.write(new_line)
                else:
                    outf.write(line)

    return job.fileStore.writeGlobalFile(modified_mapping_file)


def main(options=None):
    # if not options:
    #     # deal with command line arguments
    #     parser = ArgumentParser()
    #     Job.Runner.addToilOptions(parser)
    #     parser.add_argument(
    #         '--ref_file', default="/home/robin/paten_lab/kube_toil_minimap2_cactus/hg38_chr21.fa", help='replace_me', type=str)
    #         # '--ref_file', default="/home/robin/paten_lab/kube_toil_minimap2_cactus/chr21/hg38_chr21.fa", help='replace_me', type=str)
    #     parser.add_argument(
    #         '--assemblies_dir', default="/home/robin/paten_lab/kube_toil_minimap2_cactus/small_chr21/assemblies", help='replace_me', type=str)
    #     parser.add_argument(
    #         '--output_file', default="/home/robin/paten_lab/kube_toil_minimap2_cactus/small_chr21/output_alignment_2", help='replace_me', type=str)
    #     parser.add_argument('--mapq_cutoff', default=20,
    #                         help='replace_me', type=int)
    #     parser.add_argument('--minimum_size_remap', default=100, help='replace_me', type=int)
    #     parser.add_argument('--sequence_context',
    #                         default=10000, help='replace_me', type=int)

    #     options = parser.parse_args()
    #     options.clean = "always"

    # Now we are ready to run
    with Toil(options) as workflow:
        if not workflow.options.restart:
            # reference file
            print(os.path.abspath(options.ref_file))
            ref_file_url = 'file://' + os.path.abspath(options.ref_file)
            print(ref_file_url)
            ref_id = workflow.importFile(ref_file_url)

            # list of assembly files
            assembly_file_names = os.listdir(options.assemblies_dir)
            print(assembly_file_names)
            assembly_files = list()
            for assembly_file_name in assembly_file_names:
                assembly_file_url = 'file://' + os.path.abspath(options.assemblies_dir) + "/" + assembly_file_name
                assembly_files.append(workflow.importFile(assembly_file_url))

            if options.no_duplicate_contig_ids == False:
                # if the user hasn't guaranteed that the contig ids are all unique, then rename 
                # any that are duplicates.
                edited_assembly_files = workflow.start(Job.wrapJobFn(rename_duplicate_contig_ids, assembly_files))
                
                if options.overwrite_assemblies == True:
                    for i in range(len(edited_assembly_files)):
                        workflow.exportFile(edited_assembly_files[i], 'file://' + os.path.abspath(options.assemblies_dir) + "/" + assembly_file_names[i])
                else:
                    # make sure that the folder we want to save the assembly files in exists:
                    edited_assemblies_save_folder = os.path.abspath(options.assemblies_dir) + "_edited_for_duplicate_contig_ids/"
                    if not os.path.isdir(edited_assemblies_save_folder):
                        os.mkdir(edited_assemblies_save_folder)
                    for i in range(len(edited_assembly_files)):
                        # rename_duplicated_contig_ids outputs the edited_assembly_files in the same order as the input list of original assembly files.
                        old_assembly_name = assembly_file_names[i]
                        edited_assembly = edited_assembly_files[i]
                        # save the new assembly files:
                        workflow.exportFile(edited_assembly, 'file://' + edited_assemblies_save_folder + old_assembly_name)

                # now, replace assembly_files with the edited assembly_files:
                assembly_files = edited_assembly_files

            # perform the alignments:
            alignments = workflow.start(Job.wrapJobFn(
                align_all_assemblies,  ref_id, assembly_files, options=options))

            # reformat the alignments as lastz cigars:
            (lastz_cigar_primary_alignments, lastz_cigar_secondary_alignments) = workflow.start(Job.wrapJobFn(
                make_lastz_output, alignments))
                
            workflow.exportFile(lastz_cigar_primary_alignments, 'file://' + os.path.abspath(options.primary_output_file))
            workflow.exportFile(lastz_cigar_secondary_alignments, 'file://' + os.path.abspath(options.secondary_output_file))

        else:
            output = workflow.restart()

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
    
if __name__ == "__main__":
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    # small chr21 example:
    parser.add_argument(
        '--ref_file', default="chr21/hg38_chr21.fa", help='replace_me', type=str)
    parser.add_argument(
        '--assemblies_dir', default="small_chr21/assemblies", help='replace_me', type=str)
    parser.add_argument(
        '--primary_output_file', default="small_chr21_test_primary.out", help='replace_me', type=str)
    parser.add_argument(
        '--secondary_output_file', default="small_chr21_test_secondary.out", help='replace_me', type=str)
    parser.add_argument('--minimum_size_remap', default=100, help='replace_me', type=int)
    parser.add_argument('--no_duplicate_contig_ids', type=str2bool, nargs='?', const=True, default=False,
                        help="replace_me")
    parser.add_argument('--overwrite_assemblies', type=str2bool, nargs='?', const=True, default=False,
                        help="When cleaning the assembly files to make sure there are no duplicate contig ids, don't overwrite the assembly files. Copy them to a neigboring folder with the affix '_edited_for_duplicate_contig_ids' instead.")
    parser.add_argument('--mapq_cutoff', default=20,
                        help='replace_me', type=int)
    parser.add_argument('--sequence_context', default=10000,
                        help='replace_me', type=int)

    options = parser.parse_args()
    main(options=options)
