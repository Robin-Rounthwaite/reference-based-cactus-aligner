#Plan is to be able to convert minimap cigars in a sam file into lastz cigars.
# This should help with joining the output of my minimap2 cactus alignment algorithm
# into Cactus itself. 

from cigar import Cigar
from toil.job import Job

def get_query_start(cigar_string):
    """
    Reads the cigar string, returns how many bases from index 0 the alignment begins.
    """
    cig = Cigar(cigar_string)

    #todo: make it so that we can check for bad cigs:
    # if cig == ["not a real cig"]:
    #     #then we've been given an empty cigar string, which shouldn't ever happen.
    #     print("WARNING: _get_start_ has been requested of an empty cigar string. "
    #             + "This shouldn't happen.")
    #     return
    
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
    
    #todo: make it so that we can check for bad cigs:
    # if cig == ["not a real cig"]:
        # #then we've been given an empty cigar string, which shouldn't ever happen.
        # print("WARNING: _get_start_ has been requested of an empty cigar string. "
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

def get_lastz_cig_str(cigar_string):
    cig = Cigar(cigar_string)
    lastz_cig = ''
    for item in cig.items():
        lastz_cig += item[1] + " " + str(item[0]) + " "
    return lastz_cig[:-1]

# def make_lastz_output(job, sam_file, output_file, secondary_output_file):
def make_lastz_output(job, sam_file):
    output_lines = list()
    secondary_output_lines = list()
    with open(job.fileStore.readGlobalFile(sam_file)) as sam:
        for line in sam:
            parsed = line.split()
            flag = int(parsed[1])
            #check to make sure line is mapped: is there the 4 flag? If so, it's unmapped.
            if (flag%8)//4 == 1:
                # then this line is marked as unmapped. Skip it.
                continue
            cig_str = parsed[5]
            #todo: if lastz needs rev_comped cigar, rev_comp the cig_str here (if flag=16 or 32). I'm leaning towards "no, rev_comp is not needed."
            query_id = parsed[0]
            query_start = get_query_start(cig_str)
            query_stop = get_query_stop(cig_str)
            target_id = parsed[2]
            target_start = parsed[3]
            target_stop = str(int(target_start) + get_target_alignment_length(cig_str))
            # check to see if strand is + or -. if flag 16 or 32 is flipped, then we have strand -. 
            if (flag%64)//32 == 1 or (flag%32)//16 == 1:
                strand = "-"
            else:
                strand = "+"
            # todo: change score to raw score? Or does lastz use mapq after all?
            score = parsed[4]
            lastz_cig = get_lastz_cig_str(cig_str)

            #todo: if it's a secondary alignment (flag=256), then put it in a separate file.
            if (flag%512)//256 == 1:
                secondary_output_lines.append(" ".join(["cigar:", query_id, query_start, query_stop, target_id, target_start, target_stop, strand, score, lastz_cig, "\n"]))
            
            #todo: if it's a "supplementary alignment," then. . . What does that mean?
            output_lines.append(" ".join(["cigar:", query_id, query_start, query_stop, target_id, target_start, target_stop, strand, score, lastz_cig, "\n"]))

    output_file = job.fileStore.getLocalTempFile()
    secondary_output_file = job.fileStore.getLocalTempFile()

    with open(output_file, mode="w") as output:
        output.writelines(output_lines)
    with open(secondary_output_file, mode="w") as secondary_output:
        secondary_output.writelines(secondary_output_lines)

    return (job.fileStore.writeGlobalFile(output_file), job.fileStore.writeGlobalFile(secondary_output_file))
    

# def main():
#     sam_file = "/home/robin/paten_lab/kube_toil_minimap2_cactus/small_chr21/toilified_output_10k_context_20_mapq_cutoff_2_remap_thresh_100.sam"
#     output_file = "/home/robin/paten_lab/convert_cigar_minimap_to_lastz/test.out"
#     secondary_output_file = "/home/robin/paten_lab/convert_cigar_minimap_to_lastz/test_secondary.out"
#     make_lastz_output(sam_file, output_file, secondary_output_file)

# if __name__ == "__main__":
#     main()