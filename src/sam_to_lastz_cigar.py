#Plan is to be able to convert minimap cigars in a sam file into lastz cigars.
# This should help with joining the output of my minimap2 cactus alignment algorithm
# into Cactus itself. 

from argparse import ArgumentParser

from cigar import Cigar
from toil.common import Toil
from toil.job import Job
import os
import subprocess

def make_lastz_output(job, paf_file):
    """
    Makes lastz output using paftools.js. Also splits the input paf_file into two files
    in the output, one for the primary and the other for secondary.
    """
    primary = list()
    secondary = list()
    other = list()
    
    with open(job.fileStore.readGlobalFile(paf_file)) as inf:
        for line in inf:
            if "tp:A:P" in line or "tp:A:I" in line:
                #then the line is a primary output file.
                primary.append(line)
            # elif "tp:A:S" in line:
            else:
                #then the line is a secondary output file.
                secondary.append(line)

    # write output to files; convert to lastz:
    lines = [primary, secondary]
    sort_files = [job.fileStore.getLocalTempFile() for i in range(len(lines))]
    paftool_files = [job.fileStore.getLocalTempFile() for i in range(len(lines))]
    out_files = [job.fileStore.writeGlobalFile(job.fileStore.getLocalTempFile()) for i in range(len(lines))]
    for i in range(len(lines)):
        with open(sort_files[i], "w") as sortf:
            sortf.writelines(lines[i])
        with open(paftool_files[i], "w") as outf:
            subprocess.run(["paftools.js", "view", "-f", "lastz-cigar", sort_files[i]], stdout=outf)
        fix_negative_strand_mappings(paftool_files[i], job.fileStore.readGlobalFile(out_files[i]))
    return out_files

def fix_negative_strand_mappings(infile, outfile):
    """
    paftools.js outputs lastz files with a problem. On "-" strand mappings, the start, stop
    coordinates of the mapping are shown with [smaller coord], [larger coord], when in fact
    the larger coord should be first. This fixes that.
    Note: requires that outfile != infile.
    """
    with open(infile) as inf:
        with open(outfile, "w") as outf:
            for line in inf:
                parsed = line.split()
                if parsed[4] == "-" or parsed[8] == "-":
                    if parsed[4] == "-":
                        first = parsed[3]
                        parsed[3] = parsed[2]
                        parsed[2] = first
                    if parsed[8] == "-":
                        first = parsed[6]
                        parsed[6] = parsed[5]
                        parsed[5] = first
                    outf.write(" ".join(parsed) + "\n")
                else:
                    outf.write(line)
            

# def get_query_start(cigar_string):
#     """
#     Reads the cigar string, returns how many bases from index 0 the alignment begins.
#     """
#     cig = Cigar(cigar_string)

#     #todo: make it so that we can check for bad cigs:
#     # if cig == ["not a real cig"]:
#     #     #then we've been given an empty cigar string, which shouldn't ever happen.
#     #     print("WARNING: _get_start_ has been requested of an empty cigar string. "
#     #             + "This shouldn't happen.")
#     #     return
    
#     cig_list = list(cig.items())

#     if cig_list[0][1] in ["H", "S"]:
#         # then the mapping has a clipping. Return clipping length
#         # print("start", cig_list[0][0])
#         return str(cig_list[0][0])
#     else:
#         # print("start", 0)
#         return str(0)
    

# def get_query_stop(cigar_string):
#     """
#     Reads the cigar string, returns how many bases from index 0 of the original read
#         the alignment's stop is.
#     """
#     cig = Cigar(cigar_string)
    
#     #todo: make it so that we can check for bad cigs:
#     # if cig == ["not a real cig"]:
#         # #then we've been given an empty cigar string, which shouldn't ever happen.
#         # print("WARNING: _get_start_ has been requested of an empty cigar string. "
#         #         + "This shouldn't happen.")
#         # return
#     cig_list = list(cig.items())

#     # modify len(cig) to get the length of the contig up to the point of the stop:
#     stop_position = len(cig)
#     if cig_list[0][1] == "H":
#         # then the length of the hard clipping at the beginning won't be included in the
#         # contig. Manually include it:
#         stop_position += cig_list[0][0]
#     if cig_list[-1][1] == "S":
#         # then the length of the soft clipping after the mapping will be included in the
#         # contig. Manually remove it:
#         stop_position -= cig_list[-1][0]
#     return str(stop_position)


# def get_target_alignment_length(cigar_string):
#     """
#     Determines how much of the target was aligned to the query, based on the query's 
#     cig_str. This is equal to the sum of matches/substitutions ("M") and deletions ("D")
#     in the cig_str. Insertions ("I") don't count because these are bases contained in the 
#     query but not the target. 
#     """
#     cig = Cigar(cigar_string)
#     cig_list = list(cig.items())

#     target_alignment_length = 0
#     for i in cig_list:
#         if i[1] in ["M", "D"]:
#             target_alignment_length += i[0]
#     return target_alignment_length

# def get_lastz_cig_str(cigar_string):
#     cig = Cigar(cigar_string)
#     lastz_cig = ''
#     for item in cig.items():
#         if item[1] not in ["H", "S"]:
#             lastz_cig += item[1] + " " + str(item[0]) + " "
#     return lastz_cig[:-1]

# def make_lastz_output(job, sam_file):
#     output_lines = list()
#     secondary_output_lines = list()
#     debug_deleted_lines = int()
#     with open(job.fileStore.readGlobalFile(sam_file)) as sam:
#         for line in sam:
#             parsed = line.split()
#             flag = int(parsed[1])
#             #check to make sure line is mapped: is there the 4 flag? If so, it's unmapped.
#             if (flag%8)//4 == 1:
#                 # then this line is marked as unmapped. Skip it.
#                 debug_deleted_lines += 1
#                 continue
#             cig_str = parsed[5]
#             query_id = parsed[0]
#             query_start = get_query_start(cig_str)
#             query_stop = get_query_stop(cig_str)
#             target_id = parsed[2]
#             #note: target_start is 1-based in SAM and 0-based in LASTZ. Hence the -1 here.
#             target_start = str(int(parsed[3]) - 1)
#             target_stop = str(int(target_start) + get_target_alignment_length(cig_str))
#             # target_start = str(int(parsed[3]) - 1)
#             # if int(target_start) < 0:
#             #     print("------------------------------------------------------------------------target_start less than 0!")
#             # #todo: I'm currently test the -1 in the line below. Is that correct?
#             # target_stop = str(int(target_start) + get_target_alignment_length(cig_str) - 1)

#             #TODO: make sure you want both flags: 64 and 32. Main question: do I support flag 64, which is for paired end reads?
#             # check to see if query strand is + or -. if flag 16 or 32 is flipped, then we have query strand -.
#             if (flag%64)//32 == 1 or (flag%32)//16 == 1:
#                 query_strand = "-"
#             else:
#                 query_strand = "+"

#             #TODO: use raw score instead of minimap's mapq?
#             #NOTE: currently using mapq for output, even though original lastz supposedly uses raw alignment score.
#             score = parsed[4]

#             #get the "cigar alignment" part of the lastz_cigar. e.g. the "M 23 D 4 M 16" part.
#             lastz_cig_alignment = get_lastz_cig_str(cig_str)

#             # NOTE: as far as I can tell, the strand for the reference/target is always +.
#             # This is supported by the rough description of the lastz algorithm in the wiki. 
#             target_strand = "+"

#             # #NOTE: currently hardcoding in id=0| prepends to id names, in the hopes that it will be compatible with cactus-align:
#             # query_id= "id=0|" + query_id
#             # target_id= "id=0|" + target_id

#             if query_strand == "+":
#                 full_lastz_cigar = " ".join(["cigar:", query_id, query_start, query_stop, query_strand, target_id, target_start, target_stop, target_strand, score, lastz_cig_alignment, "\n"])
#             else: # query_strand == "-", swap the query coordinates, because we're walking along the strand backwards.
#                 full_lastz_cigar = " ".join(["cigar:", query_id, query_stop, query_start, query_strand, target_id, target_start, target_stop, target_strand, score, lastz_cig_alignment, "\n"])

#             #if it's a secondary alignment (flag=256), then put it in a separate file.
#             if (flag%512)//256 == 1:
#                 secondary_output_lines.append(full_lastz_cigar)
#             else:
#                 output_lines.append(full_lastz_cigar)
#         # #todo: !!!remove!!!
#         # output_lines.append("lines deleted: " + str(debug_deleted_lines))


#     output_file = job.fileStore.getLocalTempFile()
#     secondary_output_file = job.fileStore.getLocalTempFile()

#     with open(output_file, mode="w") as output:
#         output.writelines(output_lines)
#     with open(secondary_output_file, mode="w") as secondary_output:
#         secondary_output.writelines(secondary_output_lines)

#     return (job.fileStore.writeGlobalFile(output_file), job.fileStore.writeGlobalFile(secondary_output_file))
    

# def main():
#     sam_file = "./small_chr21/toilified_output_10k_context_20_mapq_cutoff_2_remap_thresh_100_frag_renamed.sam"
#     output_file = "./test.out"
#     secondary_output_file = "./test_secondary.out"
    
#     parser = ArgumentParser()
#     Job.Runner.addToilOptions(parser)
#     options = parser.parse_args()
#     options.clean = "always"
#     with Toil(options) as workflow:
#         sam_file_url = 'file://' + os.path.abspath(sam_file)
#         sam_id = workflow.importFile(sam_file_url)

#         (lastz_cigar_primary_alignments, lastz_cigar_secondary_alignments) = workflow.start(Job.wrapJobFn(
#                 make_lastz_output, sam_id))

                
#         workflow.exportFile(lastz_cigar_primary_alignments, 'file://' + os.path.abspath(output_file))
#         workflow.exportFile(lastz_cigar_secondary_alignments, 'file://' + os.path.abspath(secondary_output_file))

    
    
#     # make_lastz_output(sam_file, output_file, secondary_output_file)
    
    

# if __name__ == "__main__":
#     main()

