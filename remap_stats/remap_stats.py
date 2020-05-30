"""
Get statistics on coverage on all assembly contigs (not counting the reference.)
#todo: add in reference base counting? 

Main question I want to answer is: 
    What percentage of bases get passed along to the all-to-all alignment phase?

Secondary questions: 
    What percentage of bases are involved in mappings from the primary phase? 
    What about from the secondary phase? 
    From both combined?
"""

from argparse import ArgumentParser
import collections as col
import operator
from Bio import SeqIO
import os

def get_chroms_in_ref(ref_fasta):
    """
    Extracts names of sequences from a fasta, and returns them in a set.
    """
    chroms_in_ref = set()
    chroms = SeqIO.index(ref_fasta, "fasta")
    for chrom in chroms:
        chroms_in_ref.add(chrom)
    return chroms_in_ref

def get_mapping_coords_from_mapping_files(mapping_files, chroms_in_ref):
    """
    "Mapping coords" are the start, stop pairs representing regions involved in a mapping.
    This function extracts mapping_coords from mapping files in lastz cigar format.
    """
    ref_based_mapping_coords = col.defaultdict(list)
    all_to_all_mapping_coords = col.defaultdict(list)


    # first, get the points of coverage on all assembly contigs (not counting the reference.)
    for mapping_file in mapping_files:
        with open(mapping_file) as f:
            for line in f:
                parsed = line.split()
                if parsed[5] in chroms_in_ref:
                    # at some point in the file, we'll reach lines that are mapped between contigs
                    #  from the assemblies, without the reference being the target sequence. These
                    #  should be filtered into a different dictionary.
                    #order of points is flipped depending on strand:
                    if parsed[4] == "+":
                        ref_based_mapping_coords[parsed[1]].append((int(parsed[2]), int(parsed[3])))
                    elif parsed[4] == "-":
                        ref_based_mapping_coords[parsed[1]].append((int(parsed[3]), int(parsed[2])))
                else:
                    #order of points is flipped depending on strand:
                    if parsed[4] == "+":
                        all_to_all_mapping_coords[parsed[1]].append((int(parsed[2]), int(parsed[3])))
                    elif parsed[4] == "-":
                        all_to_all_mapping_coords[parsed[1]].append((int(parsed[3]), int(parsed[2])))
                    
    return ref_based_mapping_coords, all_to_all_mapping_coords

def get_coverage_points_from_mapping_coords(mapping_coords):
    """
    mapping_coords is a dictionary with key: seq_name and value: list((start_coordinate, stop_coordinate)).
    These may overlap, which prevents easy counting of seq bases covered by mappings.
    This function breaks mapping_coords down to mapping_coverage_points, with each point
    containing a number representing its position in the sequence and a bool, representing
    whether the point is start(=True) or stop (=False). These can be used to easily
    calculate mapping_coverage_coords, which have no overlaps.
    """
    mapping_coverage_points = col.defaultdict(list)
    for seq_id, coords_list in mapping_coords.items():
        for coord in coords_list:
            mapping_coverage_points[seq_id].append((coord[0], True))
            mapping_coverage_points[seq_id].append((coord[1], False))
    return mapping_coverage_points

def get_coverage_coords_from_coverage_points(mapping_coverage_points):
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
                    # that is, if we've found the starting point of a all_to_all current_region,
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

def get_mapping_coverage_lengths_from_coverage_coords(mapping_coverage_coords):
    """
    Returns Dictionary of {contig_id: total_bases_passed_to_involved_in_a_mapping} 
    """
    mapping_coverage_lengths = col.defaultdict(int)
    for contig_id, coords in mapping_coverage_coords.items():
        for coord in coords:
            mapping_coverage_lengths[contig_id] += coord[1] - coord[0]
    return mapping_coverage_lengths


def get_all_contig_lengths(assembly_files):
    """
    Get the lengths of each contig. For this to be compatible with the mapping output,
     needs the contig files with updated names.
    In addition, this will calculate incorrect lengths for contigs with duplicate names
     (not an issue with the updated names).
    """
    contig_lengths = dict()
    for assembly_file in assembly_files:
        contigs = SeqIO.index(assembly_file, "fasta")
        for contig_name, seq in contigs.items():
            contig_lengths[contig_name] = len(seq)
    return contig_lengths

def get_poor_mapping_coverage_coordinates(contig_lengths, mapping_coverage_coords, sequence_context, minimum_size_remap):
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
                poor_mapping_stop = mapping_coverage_coords[contig_id][0][0] + sequence_context
                if poor_mapping_stop > contig_lengths[contig_id]:
                    poor_mapping_stop = contig_lengths[contig_id]
                if poor_mapping_stop - 0 >= minimum_size_remap: # implement size threshold.
                    poor_mapping_coords[contig_id].append((0, poor_mapping_stop))
                    # print("#################################################included coords for remap:", 0, poor_mapping_stop, poor_mapping_stop - 0, (poor_mapping_stop - 0)<100)
                else:
                    pass
                    # print("1________________Blocked:", 0, poor_mapping_stop)
            for i in range(len(mapping_coverage_coords[contig_id]) - 1):
                # for every pair of mapping coords i and i + 1,
                # make a pair of (stop_from_ith_region, start_from_i+1th_region) to
                # represent the poor_mapping_coords. Include sequence_context as necessary.
                poor_mapping_start = mapping_coverage_coords[contig_id][i][1] - sequence_context
                if poor_mapping_start < 0:
                    poor_mapping_start = 0
                    
                poor_mapping_stop = mapping_coverage_coords[contig_id][i + 1][0] + sequence_context
                if poor_mapping_stop > contig_lengths[contig_id]:
                    poor_mapping_stop = contig_lengths[contig_id]

                if poor_mapping_stop - poor_mapping_start >= minimum_size_remap: # implement size threshold.
                    poor_mapping_coords[contig_id].append((poor_mapping_start, poor_mapping_stop))
                else:
                    pass
                    # print("2________________Blocked:", poor_mapping_start, poor_mapping_stop)
            if mapping_coverage_coords[contig_id][-1][1] < contig_lengths[contig_id]:
                # if the last mapping region for the contig stops before the end of
                # the contig, the last region is between the end of the mapping and the 
                # end of the contig.
                poor_mapping_start = mapping_coverage_coords[contig_id][-1][1] - sequence_context
                if poor_mapping_start < 0:
                    poor_mapping_start = 0
                if contig_lengths[contig_id] - poor_mapping_start >= minimum_size_remap: # implement size threshold.
                    poor_mapping_coords[contig_id].append((poor_mapping_start, contig_lengths[contig_id]))
                else:
                    pass
                    # print("3________________Blocked:", poor_mapping_start, contig_lengths[contig_id])

        else:
            # there isn't a good_mapping region for this contig. The full length of 
            # the contig belongs in poor_mapping_coords.
            poor_mapping_coords[contig_id].append((0, contig_lengths[contig_id]))
    return poor_mapping_coords

def get_sequence_lengths_remapped(remapped_sequence_coords):
    """
    Near-identical code to get_mapping_coverage_lengths_from_coverage_coords.
    Returns Dictionary of {contig_id: total_bases_passed_to_all_to_all_phase} 
    """
    sequence_lengths_remapped = col.defaultdict(int)
    for contig_id, coords in remapped_sequence_coords.items():
        for coord in coords:
            sequence_lengths_remapped[contig_id] += coord[1] - coord[0]
    return sequence_lengths_remapped

def get_sequence_coverage(mapping_files, assembly_files, chroms_in_ref, contig_lengths, sequence_context, minimum_size_remap):
    """
    Gets several dictionaries useful for calculating stats for sequence coverage.
    This is passed to print_stats
    """
    ref_based_mapping_coords, all_to_all_mapping_coords = get_mapping_coords_from_mapping_files(mapping_files, chroms_in_ref)
    
    ref_based_mapping_coverage_points = get_coverage_points_from_mapping_coords(ref_based_mapping_coords)
    ref_based_mapping_coverage_coords = get_coverage_coords_from_coverage_points(ref_based_mapping_coverage_points)
    ref_based_mapping_coverage_lengths = get_mapping_coverage_lengths_from_coverage_coords(ref_based_mapping_coverage_coords)

    all_to_all_mapping_coverage_points = get_coverage_points_from_mapping_coords(all_to_all_mapping_coords)
    all_to_all_mapping_coverage_coords = get_coverage_coords_from_coverage_points(all_to_all_mapping_coverage_points)
    all_to_all_mapping_coverage_lengths = get_mapping_coverage_lengths_from_coverage_coords(all_to_all_mapping_coverage_coords)

    remapped_sequence_coords = get_poor_mapping_coverage_coordinates(contig_lengths, ref_based_mapping_coverage_coords, sequence_context, minimum_size_remap)
    sequence_lengths_remapped = get_sequence_lengths_remapped(remapped_sequence_coords)

    return ref_based_mapping_coverage_lengths, all_to_all_mapping_coverage_lengths, sequence_lengths_remapped

def print_stats(ref_based_mapping_coverage_lengths, all_to_all_mapping_coverage_lengths, sequence_lengths_remapped, contig_lengths, contig_stats, raw_data):
    """
    For playing with following questions:
        What percentage of bases are involved in mappings from the primary phase? 
        What about from the secondary phase? 
        From both combined?
        What percentage of bases get passed along to the all-to-all alignment phase?
    """
    #first: how many bases are there involved in the assemblies?
    bases_total = int()
    for length in contig_lengths.values():
        bases_total += length
    print("total number of bases in all input assemblies:", bases_total) 
    print()
    #what is the percentage of bases involved in mappings from each phase?
    print("***analysis of ref_based mappings:***")
    ref_based_bases_covered = int()
    for length in ref_based_mapping_coverage_lengths.values():
        ref_based_bases_covered += length
    print("bases covered in ref_based mappings:", ref_based_bases_covered)
    print("percentage of bases covered by mappings from ref_based:", (ref_based_bases_covered/bases_total)*100, "%")

    print()
    print("***analysis of all_to_all mappings:***")
    all_to_all_bases_covered = int()
    for length in all_to_all_mapping_coverage_lengths.values():
        all_to_all_bases_covered += length
    print("bases covered in ref_based mappings:", all_to_all_bases_covered)
    print("percentage of bases covered by mappings from all_to_all:", (all_to_all_bases_covered/bases_total)*100, "%")

    print()
    print("total percentage of bases covered by mappings:", ((ref_based_bases_covered + all_to_all_bases_covered)/bases_total)*100, "%")

    print()
    print("***analysis of quantity of sequence passed to the all-to-all alignment phase:***")
    bases_remapped = int()
    for length in sequence_lengths_remapped.values():
        bases_remapped += length
    print("bases passed to all-to-all phase:", bases_remapped)
    print("percentage of bases sent to all-to-all phase:", (bases_remapped/bases_total)*100, "%")
    print()
    print("ratio of (bases covered by mappings in all-to-all)/(bases sent to all-to-all-phase) :", ((all_to_all_bases_covered/bases_remapped)))

    if contig_stats:
        print()
        print("***analysis of coverage levels on individual contigs in assemblies:***")
        for contig_id in contig_lengths:
            print()
            print("stats_for " + contig_id + ":")
            print("contig_length:", contig_lengths[contig_id])
            print("bases_mapped_in_ref_based_phase/contig_length", ref_based_mapping_coverage_lengths[contig_id]/contig_lengths[contig_id])
            print("bases_passed_to_all_to_all_phase/contig_length", sequence_lengths_remapped[contig_id]/contig_lengths[contig_id])
            print("bases_mapped_in_all_to_all_phase/contig_length", all_to_all_mapping_coverage_lengths[contig_id]/contig_lengths[contig_id])

    if raw_data:
        print()
        print("***printing out raw data as python dictionaries:***")
        print("ref_based_mapping_coverage_lengths:\n", ref_based_mapping_coverage_lengths) 
        print("all_to_all_mapping_coverage_lengths:\n", all_to_all_mapping_coverage_lengths) 
        print("sequence_lengths_remapped:\n", sequence_lengths_remapped) 
        print("contig_lengths:\n", contig_lengths)
            

def main(options=None):
    if not options:
        #deal with command line arguments: NOTE: THIS IS NOT COMPATIBLE WITH JUPYTER, must
        #  always pass args in through options=options. (which works because I use
        #  parser.parse_known_args)
        parser = ArgumentParser()
        parser.add_argument(
            'ref_file', help='The reference fasta file for the initial alignment phase of the reference-based-cactus-aligner.', type=str)
        parser.add_argument('assemblies_dir', help='replace_me', type=str)
        parser.add_argument('mappings_dir', help='replace_me', type=str)
        parser.add_argument('--minimum_size_remap', default=100, help='replace_me', type=int)
        parser.add_argument('--sequence_context', default=10000, help='replace_me', type=int)
        parser.add_argument('--contig_stats', action='store_true')
        parser.add_argument('--raw_data', action='store_true')
        options = parser.parse_args()

    ## parsing options:s
    # list of mapping files
    mapping_files = os.listdir(options.mappings_dir)
    mapping_files = [options.mappings_dir + mapping_files[i] for i in range(len(mapping_files))]
    # list of assembly files
    assembly_files = os.listdir(options.assemblies_dir)
    assembly_files = [options.assemblies_dir + assembly_files[i] for i in range(len(assembly_files))]
    # list of chromosomes in ref:
    chroms_in_ref = get_chroms_in_ref(options.ref_file)

    contig_lengths = get_all_contig_lengths(assembly_files)
    ref_based_mapping_coverage_lengths, all_to_all_mapping_coverage_lengths, sequence_lengths_remapped = get_sequence_coverage(mapping_files, assembly_files, chroms_in_ref, contig_lengths, options.minimum_size_remap, options.sequence_context)
    print_stats(ref_based_mapping_coverage_lengths, all_to_all_mapping_coverage_lengths, sequence_lengths_remapped, contig_lengths, options.contig_stats, options.raw_data)
    # #todo: return here is for debug purposes in jupyter.
    # return ref_based_mapping_coverage_lengths, all_to_all_mapping_coverage_lengths, sequence_lengths_remapped, contig_lengths


if __name__ == "__main__":
    ### example command line call to test the small_chr21 test set with normal main():
    ### python remap_stats.py ../chr21/hg38_chr21.fa ../small_chr21/assemblies_edited_for_duplicate_contig_ids/ small_chr21_output/
    main()
    #for testing the small_chr21 test set:
    # parser = ArgumentParser()
    # options, unknown = parser.parse_known_args()
    # options.mappings_dir = "small_chr21_output/"
    # options.assemblies_dir = "../small_chr21/assemblies_edited_for_duplicate_contig_ids/"
    # options.ref_file = "../chr21/hg38_chr21.fa"
    # options.minimum_size_remap = 100
    # options.sequence_context = 10000
    # options.contig_stats = True
    # options.raw_data = True
    # ref_based_mapping_coverage_lengths, all_to_all_mapping_coverage_lengths, sequence_lengths_remapped, contig_lengths = main(options=options)
