"""
A pytest battery for verifying that the sequences extracted and manipulated by
minimap2_cactus_aligner_toil.py act as expected.
"""
import minimap2_cactus_aligner_toil

from toil.common import Toil
from toil.job import Job
from argparse import ArgumentParser


#get options

def get_options():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    options.test_fasta = "my_test_fasta.fa"
    options.test_sam = "my_test_sam.fa"
    options.mapq_cutoff = 20
    options.sequence_context = 5 # currently normally 10000
    options = parser.parse_args()
    
    return options

# test get_mapping_coverage_points
def test_get_mapping_coverage_points():
    options = get_options()
    with Toil(options) as workflow:
        mapping_coverage_points = workflow.start(Job.wrapJobFn(get_mapping_coverage_points, options.test_sam, options))
        



# test get_mapping_coverage_coordinates

# test get_poor_mapping_coverage_coordinates

# test get_poor_mapping_sequences

# test relocate_remapped_fragments_to_source_contigs

# test sam_to_lastz_cigar.py 

##Test draft for making sure that cigars have proper coordinates.
#%%
from Bio import SeqIO
from cigar import Cigar
import os

def directly_calculate_contig_lengths(assembly_file):
    contig_lengths = dict()
    contigs = SeqIO.index(assembly_file, "fasta")
    for contig_name, seq in contigs.items():
        contig_lengths[contig_name] = len(seq)
    return contig_lengths

# cig_str = "S 20 M 374 D 2 M 252 D 3 M 63 D 2 M 139 D 2 M 194 D 2 M 11 D 1 M 75 S 30"
def lastz_to_sam_cig(lastz_cig_str):
    brkdwn = lastz_cig_str.split()
    lets = list()
    nums = list()
    for i in range(len(brkdwn)):
        if i%2 == 0:
            #it's a letter:
            lets.append(brkdwn[i])
        else:
            #it's a num:
            nums.append(brkdwn[i])
    orig_cig = ""
    for i in range(len(lets)):
        orig_cig += nums[i] + lets[i]

    return orig_cig

def get_len_target_consumed(cig_str):
    """M + D length is length of target consumed."""
    cig = Cigar(cig_str)
    target_consumed = 0
    for i in cig.items():
        if i[1] in ["M", "D"]:
            target_consumed += int(i[0])
    return target_consumed

def get_len_query_consumed(cig_str):
    """M + I length is length of query consumed."""
    cig = Cigar(cig_str)
    query_consumed = 0
    for i in cig.items():
        if i[1] in ["M", "I"]:
            query_consumed += int(i[0])
    return query_consumed


# %%
###This should all be moved into my unit tests after I'm done with it here. (including above cell)
# goal: identify any size mismatches (that appear to be causing segfaults in cactus) in the Cigars in the sam file/lastz files
#       - if there is a well-defined mismatch, see if the coordinate size mismatch is preserved between the query and the target
###put this in tests!!


# get size of contigs:
fasta_files = ["./small_chr21/assemblies_edited_for_duplicate_contig_ids/HG03098_paf_chr21.fa", "./small_chr21/assemblies_edited_for_duplicate_contig_ids/HG03492_paf_chr21.fa", "./chr21/hg38_chr21.fa"]
contig_lengths = dict()
for fa in fasta_files:
    contig_lengths.update(directly_calculate_contig_lengths(fa))
print(contig_lengths)

#%%
lastz_file = "./small_chr21_test_primary.out"
with open(lastz_file) as inf:
    for line in inf:
        line = "cigar: 6348 26932 0 - chr21 6906547 6933497 + 45 M 652 D 6 M 530 I 1 M 404 D 2 M 139 D 1 M 149 I 1 M 1920 I 1 M 206 I 1 M 446 D 1 M 234 D 1 M 281 D 1 M 9 I 1 M 354 D 3 M 692 D 2 M 485 I 1 M 137 I 1 M 43 D 2 M 195 I 1 M 1786 D 1 M 2538 I 1 M 54 D 4 M 1515 D 2 M 327 I 5 M 1219 D 1 M 748 I 11 M 673 D 14 M 46 I 12 M 984 D 1 M 508 I 1 M 8 I 1 M 1277 D 1 M 411 D 1 M 2533 D 2 M 616 D 1 M 59 D 1 M 244 D 1 M 7 D 1 M 531 D 3 M 505 I 1 M 139 D 2 M 485 D 1 M 153 I 1 M 1579 I 1 M 300 D 1 M 203 D 1 M 37 D 2 M 529 " # line 590
        # line = "cigar: 11146 10607 8584 - 10208 1 2022 + 60 M 52 I 1 M 291 I 1 M 337 I 2 M 119 I 2 M 44 I 1 M 113 I 1 M 306 D 1 M 406 D 3 M 19 D 2 M 25 I 1 M 30 D 1 M 272" 
        # print(line)
        parsed = line.split()
        # print(parsed)
        query_start = parsed[2]
        query_stop = parsed[3]
        query_strand = parsed[4]
        target_start = parsed[6]
        target_stop = parsed[7]
        lastz_cig_str = " ".join(parsed[10:])
        # print(lastz_cig_str)
        cig_str = lastz_to_sam_cig(lastz_cig_str)
        cig = Cigar(cig_str)
        cig_len = len(cig)

        target_consumed = get_len_target_consumed(cig_str)
        query_consumed = get_len_query_consumed(cig_str)
        query_length = contig_lengths[parsed[1]]
        target_length = contig_lengths[parsed[5]]

        target_length_involved_in_mapping = int(target_stop) - int(target_start)
        query_length_involved_in_mapping = int(query_stop) - int(query_start)
        #todo: make "target_length", see if makes sense with my second cigar.
        print("cig_len", cig_len, "target_length_involved_in_mapping", target_length_involved_in_mapping, "target_consumed", target_consumed, "query_length_involved_in_mapping", query_length_involved_in_mapping, "query_consumed", query_consumed, "query_length", query_length, "target_length", target_length)

        # # check for starts/stops that are past the allowed borders of the actual contigs.
        # if int(query_start) < 0:
        #     print("line with query start (", query_start, ") is less than zero:", line)
        # if int(query_start) > query_length:  # covers case with query strand "-"
        #     print("line with query start (", query_start, ") past query length (", query_length, "):", line)

        # if int(query_stop) > query_length:
        #     print("line with query stop (", query_stop, ") past query length (", query_length, "):", line)
        # if int(query_stop) < 0: # covers case with query strand "-"
        #     print("line with query stop (", query_stop, ") is less than zero:", line)

        # if int(target_start) > target_length: # covers case with target strand "-"
        #     print("line with target start (", target_start, ") past target length (", target_length, "):", line)
        # if int(target_start) < 0:
        #     print("line with target start (", target_start, ") is less than zero:", line)

        # if int(target_stop) > target_length:
        #     print("line with target stop (", target_stop, ") past target length (", target_length, "):", line)
        # if int(target_stop) < 0:
        #     print("line with target stop (", target_stop, ") is less than zero:", line)
        break

# %%
#problem with line 359 of small_chr21 output: the target length extends past the maximum size of the target.
#first check .sam.
contig_lengths["10208"]


print(lastz_to_sam_cig("M 52 I 1 M 291 I 1 M 337 I 2 M 119 I 2 M 44 I 1 M 113 I 1 M 306 D 1 M 406 D 3 M 19 D 2 M 25 I 1 M 30 D 1 M 272"))
# 52M1I291M1I337M2I119M2I44M1I113M1I306M1D406M3D19M2D25M1I30M1D272M

print(lastz_to_sam_cig("M 652 D 6 M 530 I 1 M 404 D 2 M 139 D 1 M 149 I 1 M 1920 I 1 M 206 I 1 M 446 D 1 M 234 D 1 M 281 D 1 M 9 I 1 M 354 D 3 M 692 D 2 M 485 I 1 M 137 I 1 M 43 D 2 M 195 I 1 M 1786 D 1 M 2538 I 1 M 54 D 4 M 1515 D 2 M 327 I 5 M 1219 D 1 M 748 I 11 M 673 D 14 M 46 I 12 M 984 D 1 M 508 I 1 M 8 I 1 M 1277 D 1 M 411 D 1 M 2533 D 2 M 616 D 1 M 59 D 1 M 244 D 1 M 7 D 1 M 531 D 3 M 505 I 1 M 139 D 2 M 485 D 1 M 153 I 1 M 1579 I 1 M 300 D 1 M 203 D 1 M 37 D 2 M 529"))
#%%

sam_file = "./small_chr21/toilified_output_10k_context_20_mapq_cutoff_2_remap_thresh_100_frag_renamed.sam"
with open(sam_file) as inf:
    for line in inf:
        # print(line)
        parsed = line.split()
        target_start_pos = parsed[3]

        if target_start_pos == "1" and parsed[1] != "4":
            print(parsed[0:9])
# %%
