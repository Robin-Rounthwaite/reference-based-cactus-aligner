import cigar
from Bio import SeqIO
from argparse import ArgumentParser
import os

"""
Test cigar length is same length as mapping length, both in terms of query and target:
"""

def test_mapping_len_is_cigar_len(mapping_files):
    """
    Test for making sure that the length of the mappings equal the length of the cigars  
    """
    line_cnt = int()
    for mappings in mapping_files:
        with open(mappings) as inf:
            for line in inf:
                parsed = line.split()
                query_length = int(parsed[3])-int(parsed[2])
                ref_length = int(parsed[7])-int(parsed[6])
                cig = cigar.Cigar(lastz_to_sam_cig(parsed[10:]))
                if query_length != len(cig):
                    if query_length != -len(cig): # additional check just for checking my bad-formatted cigars.
                        print("query length off:")
                        print(cig)
                        print(query_length, len(cig))

                if ref_length != cig.reference_length():
                    print("ref length off:")
                    print(cig)
                    print(ref_length, cig.reference_length())
                # print(len(cig))
                # print(cig.reference_length())
                # break
                line_cnt+=1 
    print(line_cnt, "lines tested.")

def lastz_to_sam_cig(lastz_cig_list):
    # lastz_cig_list = cig_str.split()
    cig = str()
    for i in range(0, len(lastz_cig_list), 2):
        cig += lastz_cig_list[i+1] + lastz_cig_list[i]
    return cig


"""
Generate contig_lengths:
"""
def get_merged_contig_lengths(reference_file, assembly_files):
    # first, get contig_lengths separated by assembly file:
    contig_lengths = get_contig_lengths(reference_file, assembly_files)

    # Check that all contig length are unique (essential for the proper function of 
    # mapping pipeline, as well as for proper merging of all contig lengths into a single
    # dictionary)
    check_all_contig_ids_unique(contig_lengths)

    # Consolidate contig_lengths so each contig_id is a key in the dictionary, and length is value.
    all_contig_lengths = dict()
    for assemb_file, length_dict in contig_lengths.items():
        all_contig_lengths.update(length_dict)

    print("found", len(all_contig_lengths), "contig lengths.")
    return all_contig_lengths

def get_contig_lengths(reference_file, assembly_files):
    # calculate lengths of all input contigs. + reference:
    all_contig_lengths = dict()
    all_contig_lengths[reference_file] = count_contig_lengths(reference_file)
    for assembly in assembly_files:
        all_contig_lengths[assembly] = count_contig_lengths(assembly)
    return all_contig_lengths

    
def count_contig_lengths(assembly_file):
    contig_lengths = dict()
    contigs = SeqIO.index(assembly_file, "fasta")
    for contig_name, seq in contigs.items():
        contig_lengths[contig_name] = len(seq)
    return contig_lengths

def check_all_contig_ids_unique(contig_lengths):
    """
    Is every id in contig_lengths unique?
    """
    all_ids = set()
    for assemb, contig_lens in contig_lengths.items():
        for contig_id in contig_lens:
            if contig_id in all_ids:
                print("WARNING: duplicate contig id found (so, broken pipeline):", contig_id)
            all_ids.add(contig_id)
            
"""
Test for making sure that the length of the mappings doesn't extend past the length of the
contigs.
"""
def test_mapping_coords_within_contig_lengths(mapping_files, contig_lengths):
    line_count = int()
    for mappings in mapping_files:
        with open(mappings) as inf:
            for line in inf:
                parsed = line.split()
                # See if the query mapping coords match total query length
                if int(parsed[2]) < 0:
                    print("WARNING: query mapping start<0.")
                if int(parsed[2]) >= contig_lengths[parsed[1]]:
                    print("WARNING: query mapping start extends past length of contig.")
                if int(parsed[3]) >= contig_lengths[parsed[1]]:
                    print("WARNING: query mapping end extends past length of contig.")
                if int(parsed[3]) < 0:
                    print("WARNING: query mapping end<0.")

                # See if the target mapping coords match total target length
                if int(parsed[6]) < 0:
                    print("WARNING: target mapping start<0.")
                if int(parsed[6]) >= contig_lengths[parsed[5]]:
                    print("WARNING: target mapping start extends past length of contig.")
                if int(parsed[7]) >= contig_lengths[parsed[5]]:
                    print("WARNING: target mapping end extends past length of contig.")
                if int(parsed[7]) < 0:
                    print("WARNING: target mapping end<0.")

                line_count += 1

    print(line_count, "lines tested")

"""
Main file and command line interface:
"""
def get_options():
    parser = ArgumentParser()
    parser.add_argument('ref_file', help='replace_me', type=str)
    parser.add_argument('assemblies_dir', help='The assemblies used in the pipeline, with all duplicatee contig_ids deduplicated.', type=str)
    parser.add_argument('mappings_dir', help='replace_me', type=str)
    return parser.parse_args()

def main():
    # parse options
    options = get_options()
    
    # reference file
    ref_file = os.path.abspath(options.ref_file)

    # list of assembly files
    assembly_file_names = os.listdir(options.assemblies_dir)
    assembly_files = list()
    for assembly in assembly_file_names:
        assembly_files.append(os.path.abspath(options.assemblies_dir) + "/" + assembly)
    
    # list of mapping_files
    mapping_file_names = os.listdir(options.mappings_dir)
    mapping_files = list()
    for mapping in mapping_file_names:
        mapping_files.append(os.path.abspath(options.mappings_dir) + "/" + mapping)
    
    # mapping_files = ["../primary.cigar", "../secondary.cigar"]
    
    # ref_file = "../chr21/hg38_chr21.fa"
    # assembly_files = ["../small_chr21/assemblies_edited_for_duplicate_contig_ids/HG03098_paf_chr21.fa", "../small_chr21/assemblies_edited_for_duplicate_contig_ids/HG03492_paf_chr21.fa"]
    
    print("Testing that all mapping lengths are equal to cigar lengths.")
    test_mapping_len_is_cigar_len(mapping_files)
    
    print("\nCalculating contig_lengths (and seeing if there are duplicate contig_ids)")
    contig_lengths = get_merged_contig_lengths(ref_file, assembly_files)

    print("\nTesting that all mapping start and end-points are above 0 and <= length of contig.")
    test_mapping_coords_within_contig_lengths(mapping_files, contig_lengths)

if __name__ == "__main__":
    main()
