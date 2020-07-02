from toil.job import Job

import subprocess
from Bio import SeqIO
import cigar
import operator
import collections as col



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

def get_poor_mapping_coverage_coordinates(job, contig_lengths, assembly_file, mapping_coverage_coords, options):
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

def get_poor_mapping_sequences(job, assembly_file, poor_mapping_coords, options):
    """
    ---Sequence extraction:---
    Read in the entire fasta file.
    
    for each contig-subdivided list of (start, stop) coordinates in 
    self.files.loc["fasta_files"], extract the sequence associated with 
    each coordinate in that contig.
    
    If there are multiple identical copies of the same region, this only includes one of 
    them.
    """
    
    # make fasta file for later remapping all_to_all.
    poor_mapping_sequence_file = job.fileStore.getLocalTempFile()

    contigs = SeqIO.index(job.fileStore.readGlobalFile(assembly_file), "fasta")

    sequences_written = set()
    with open(poor_mapping_sequence_file, "w+") as outf:
        for contig_name, contig_record in contigs.items():
            for coord in poor_mapping_coords[contig_name]:
                # for each coord in the low_mapq_coords corresponding to a specific contig:
                # extract the sequence for that contig.
                sequence_name = ">" + contig_name + "_segment_start_" + str(coord[0]) + "_stop_" + str(coord[1])
                if sequence_name not in sequences_written:
                    low_mapq_sequence = contig_record.seq[coord[0]: coord[1]]
                    outf.write(sequence_name + "\n")
                    outf.write(str(low_mapq_sequence) + "\n")
                    sequences_written.add(sequence_name)
            
    return job.fileStore.writeGlobalFile(poor_mapping_sequence_file)

def remap_poor_mapping_sequences(job, poor_mapping_sequence_file, assembly_to_align_file, assembly_files, options):
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

    # if options.export_all_to_all_files:
    #     for remap_file in remapping_files:
    #         options.all_to_all_mapping_files.append(remap_file)

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

def relocate_remapped_fragments_to_source_contigs(job, contig_lengths, mapping_file, fasta_file):
    """
    renames contig fragments from the remapping step to the original
    contig name. Changes the cigar clipping based on the fragment's coordinates too. 
    Arguments:
        mapping_file {[type]} -- [description]
        fasta_file {[type]} -- [description]
    """
    
    modified_mapping_file = job.fileStore.getLocalTempFile()

    debug_line_count = int()
    debug_first_if = int()
    debug_second_if = int()
    debug_third_if = int()

    with open(job.fileStore.readGlobalFile(mapping_file)) as inf:
        with open(modified_mapping_file, "w") as outf:
            for line in inf:
                debug_line_count += 1
                parsed = line.split()
                if (int(parsed[1])%8)//4==1:
                    # if the line is flagged as unmapped, just write it to the outfile.
                    outf.write(line)
                    debug_first_if += 1
                elif "segment" in parsed[0]:
                    debug_second_if += 1
                    # construct the new name by dropping all the parts of the old name from "_segment_" onwards.
                    name_parsed = parsed[0].split("_segment_")
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
                    new_line = new_name + "\t" + "\t".join(parsed[1:5]) + "\t" + new_cig + "\t" + "\t".join(parsed[6:]) + "\n"
                    # print("line after modification: ", new_line)

                    #add it to the outfile
                    outf.write(new_line)
                else:
                    debug_third_if += 1
                    outf.write(line)
            # outf.write("debug_line_count:" + str(debug_line_count) + " 0 0 0 0 10M 0 0\n")
            # outf.write("debug_first_if:" + str(debug_first_if) + " 0 0 0 0 10M 0 0\n")
            # outf.write("debug_second_if:" + str(debug_second_if) + " 0 0 0 0 10M 0 0\n")
            # outf.write("debug_third_if:" + str(debug_third_if) + " 0 0 0 0 10M 0 0\n")

    return job.fileStore.writeGlobalFile(modified_mapping_file)
