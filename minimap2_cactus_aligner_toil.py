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
from toil.common import Toil
from toil.job import Job

from argparse import ArgumentParser
import os
from Bio import SeqIO

from sam_to_lastz_cigar import make_lastz_output
from remap_stats.online_remap_stats import online_remap_stats
from src import mapping_functions
from src import fasta_preprocessing


def align_all_assemblies(job, reference_file, assembly_files, remap_stats_internal_file, options):
    """
    The head job for the main function of this workflow. Iterates through each assembly in
    assembly_files, and calls align_assembly on it. This aligns it to the reference, 
    then aligns poorly-mapping regions to the rest of assembly_files.
    """
    mapping_files = list()

    # get lengths of all contigs
    head_job = job.addChildJobFn(empty)
    all_contig_lengths_job = head_job.addChildJobFn(calc_all_contig_lengths, reference_file, assembly_files)
    all_contig_lengths = all_contig_lengths_job.rv()

    if options.remap_stats:
        all_contig_lengths_job.addFollowOnJobFn(online_remap_stats.save_input_assembly_stats, all_contig_lengths, remap_stats_internal_file, options)

    if options.export_all_to_all_files:
        all_to_all_fastas = list()
        all_to_all_sams = list()

    for assembly_to_align_file in assembly_files:
        assembly_mapping_file_job = all_contig_lengths_job.addFollowOnJobFn(align_assembly, reference_file, assembly_files, assembly_to_align_file, all_contig_lengths, remap_stats_internal_file, options)
        assembly_mapping_file = assembly_mapping_file_job.rv()
        mapping_files.append(assembly_mapping_file)
        
        if options.export_all_to_all_files:
            # Then assembly_mapping_file is actually a triplet, including useful data on all_to_all.
            assembly_mapping_file = assembly_mapping_file_job.addFollowOnJobFn(unpack_promise, triplet, 0).rv()
            poor_mapping_sequence_file = assembly_mapping_file_job.addFollowOnJobFn(unpack_promise, triplet, 1).rv()
            all_to_all_mapping_file = assembly_mapping_file_job.addFollowOnJobFn(unpack_promise, triplet, 2).rv()
            
            all_to_all_fastas.append(poor_mapping_sequence_file)
            all_to_all_sams.append(all_to_all_mapping_file)

    if options.export_all_to_all_files:
        return job.addFollowOnJobFn(mapping_functions.consolidate_mapping_files, mapping_files).rv(), all_to_all_fastas, all_to_all_sams
    else:
        return job.addFollowOnJobFn(mapping_functions.consolidate_mapping_files, mapping_files).rv()


def align_assembly(job, reference_file, assembly_files, assembly_to_align_file, all_contig_lengths, remap_stats_internal_file, options):
    """
    Aligns assembly_to_align_file to the reference, then aligns poorly-mapping regions to the rest of assembly_files.
    """
    # mapping_files compiles all the mappings to reference and mappings between assemblies.
    mapping_files = list()
    
    contig_lengths = all_contig_lengths[assembly_to_align_file]

    ## map to reference phase:
    # map assembly to reference. Get the id of the map-to-ref file.
    map_to_ref_job = job.addChildJobFn(mapping_functions.map_assembly_to_ref, assembly_to_align_file, reference_file)
    map_to_ref_file = map_to_ref_job.rv()
    mapping_files.append(map_to_ref_file)

    if not options.all_to_ref_only:
        ## re-mapping phase, from regions of the assemblies that map poorly to the reference, 
        ## to all the assembly sequence.
        # extract all the start & stop regions of good mapping regions in the assemblies
        # (this doesn't yet include options.sequence_context, but does include mapq cutoff).
        mapping_coverage_points_job = map_to_ref_job.addFollowOnJobFn(mapping_functions.get_mapping_coverage_points, map_to_ref_file, options)
        mapping_coverage_points = mapping_coverage_points_job.rv()

        # assemble mapping_coverage_points into a list of proper coverage coordinates in 
        # (start_of_coverage, stop_of_coverage) format.
        mapping_coverage_coords_job = mapping_coverage_points_job.addFollowOnJobFn(mapping_functions.get_mapping_coverage_coordinates, mapping_coverage_points)
        mapping_coverage_coords = mapping_coverage_coords_job.rv()

        # find the coordinates of regions that don't have mapping coverage. These are 
        # organized in a dictionary with
        # dict{contig_id: list[tuple(start_of_poor_coverage, stop_of_poor_coverage)]} format.
        # this includes options.sequence_context, and options.minimum_size_remap.
        poor_mapping_coverage_coords_job = mapping_coverage_coords_job.addFollowOnJobFn(mapping_functions.get_poor_mapping_coverage_coordinates, contig_lengths, assembly_to_align_file, mapping_coverage_coords, options)
        poor_mapping_coverage_coords = poor_mapping_coverage_coords_job.rv()

        # extract the actual sequence that has poor mapping coverage from the contigs.
        # saved in a fasta file in the filestore.
        poor_mapping_sequence_file_job = poor_mapping_coverage_coords_job.addFollowOnJobFn(mapping_functions.get_poor_mapping_sequences, assembly_to_align_file, poor_mapping_coverage_coords, options)
        poor_mapping_sequence_file = poor_mapping_sequence_file_job.rv()

        if options.remap_stats:
            poor_mapping_sequence_file_job.addFollowOnJobFn(online_remap_stats.save_sequence_remapped_stats, assembly_to_align_file, poor_mapping_sequence_file, remap_stats_internal_file, options).rv()

        # # map the poor mapping sequence to all the other assemblies!
        map_to_assemblies_file_job = poor_mapping_sequence_file_job.addFollowOnJobFn(mapping_functions.remap_poor_mapping_sequences, poor_mapping_sequence_file, assembly_to_align_file, assembly_files, options)
        all_to_all_mapping_file = map_to_assemblies_file_job.rv()
        mapping_files.append(all_to_all_mapping_file)

        if options.remap_stats:
            # count bases involved in all_to_all mappings
            map_to_assemblies_file_job.addFollowOnJobFn(debug_print, all_to_all_mapping_file)
            map_to_assemblies_file_job.addFollowOnJobFn(online_remap_stats.save_all_to_all_mappings_stats, all_to_all_mapping_file, remap_stats_internal_file, options)

        # consolidate all the mapping_files to become a single file.
        consolidate_mapping_files_job = map_to_assemblies_file_job.addFollowOnJobFn(mapping_functions.consolidate_mapping_files, mapping_files)
        consolidated_mapping_files = consolidate_mapping_files_job.rv()


        if options.export_all_to_all_files:
            # return (consolidate_mapping_files_job.addFollowOnJobFn(relocate_remapped_fragments_to_source_contigs, contig_lengths, consolidated_mapping_files, assembly_to_align_file).rv(), poor_mapping_sequence_file, all_to_all_mapping_file)
            return (consolidate_mapping_files_job.addFollowOnJobFn(mapping_functions.relocate_remapped_fragments_to_source_contigs, contig_lengths, consolidated_mapping_files, assembly_to_align_file).rv(), poor_mapping_sequence_file, consolidate_mapping_files_job.addFollowOnJobFn(mapping_functions.relocate_remapped_fragments_to_source_contigs, contig_lengths, all_to_all_mapping_file, assembly_to_align_file).rv())
        else:
            return consolidate_mapping_files_job.addFollowOnJobFn(mapping_functions.relocate_remapped_fragments_to_source_contigs, contig_lengths, consolidated_mapping_files, assembly_to_align_file).rv()
    else:
        return map_to_ref_file

def empty(job):
    """
    An empty job, for easier toil job organization.
    """
    return

def debug_print(job, all_to_all_mapping_files):
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ", job.fileStore.readGlobalFile(all_to_all_mapping_files))

def unpack_promise(job, iterable, i):
    """
    passed an iterable and a location i, returns ith item.
    """
    return iterable[i]

def calc_all_contig_lengths(job, reference_file, assembly_files):
    # calculate lengths of all input contigs. + reference:
    all_contig_lengths = dict()
    all_contig_lengths[reference_file] = job.addChildJobFn(directly_calculate_contig_lengths, reference_file).rv()
    for assembly in assembly_files:
        all_contig_lengths[assembly] = job.addChildJobFn(directly_calculate_contig_lengths, assembly).rv()
    return all_contig_lengths
    
def directly_calculate_contig_lengths(job, assembly_file):
    contig_lengths = dict()
    contigs = SeqIO.index(job.fileStore.readGlobalFile(assembly_file), "fasta")
    for contig_name, seq in contigs.items():
        contig_lengths[contig_name] = len(seq)
    return contig_lengths


def make_remap_stats_internal_file(job):
    f = job.fileStore.getLocalTempFile()
    return job.fileStore.writeGlobalFile(f)

def get_options():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    
    parser.add_argument(
        'ref_file', help='The reference fasta file for the initial alignment phase of the reference-based-cactus-aligner.', type=str)
    parser.add_argument(
        'assemblies_dir', help='replace_me', type=str)
   
    parser.add_argument(
        '--primary_output_file', default="primary.cigar", help='replace_me', type=str)
    parser.add_argument(
        '--secondary_output_file', default="secondary.cigar", help='replace_me', type=str)
    parser.add_argument('--minimum_size_remap', default=100, help='replace_me', type=int)
    parser.add_argument('--no_duplicate_contig_ids', action='store_true', 
                        help="Don't clean the assembly files; the user promises that they don't contain any duplicate contig ids.")
    parser.add_argument('--overwrite_assemblies', action='store_true',
                        help="When cleaning the assembly files to make sure there are no duplicate contig ids, don't overwrite the assembly files. Copy them to a neigboring folder with the affix '_edited_for_duplicate_contig_ids' instead.")
    parser.add_argument('--mapq_cutoff', default=20,
                        help='replace_me', type=int)
    parser.add_argument('--sequence_context', default=10000,
                        help='replace_me', type=int)
    parser.add_argument('--all_to_ref_only', help='replace_me', action='store_true')
    
    # Mostly for debug:
    parser.add_argument('--remap_stats', action='store_true', 
                        help='directly counts the total number of bases in the reference and input assemblies; the number of bases mapped in the all-to-ref phase; the bases sent to the all-to-all phase; and the bases mapped in the all-to-all phase.')
    parser.add_argument('--remap_stats_raw', action='store_true', 
                        help='requires --remap_stats to be set. Additionally prints the dictionaries of relevant data.')
    parser.add_argument('--remap_stats_output_file', type=str, default='remap_stats_in_pipeline.txt', 
                        help='Defines where to save the remap_stats, if --remap_stats is called.')
    parser.add_argument('--export_all_to_all_files', action='store_true', 
                        help='Exports both input fasta and output sam files for all-to-all.')
    options = parser.parse_args()
    return options

def main(options=None):
    if not options:
        options = get_options()

    # Now we are ready to run
    with Toil(options) as workflow:
        if not workflow.options.restart:
            # below variables only used when specific flags are called.
            remap_stats_internal_file = workflow.start(Job.wrapJobFn(make_remap_stats_internal_file))
            # if options.export_all_to_all_files:
            #     options.all_to_all_fastas = list()
            #     options.all_to_all_mapping_files = list()
            #     print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~both lists:", options.all_to_all_fastas, options.all_to_all_mapping_files)
            
            # reference file
            ref_file_url = 'file://' + os.path.abspath(options.ref_file)
            ref_id = workflow.importFile(ref_file_url)

            # list of assembly files
            assembly_file_names = os.listdir(options.assemblies_dir)

            assembly_files = list()
            for assembly_file_name in assembly_file_names:
                assembly_file_url = 'file://' + os.path.abspath(options.assemblies_dir) + "/" + assembly_file_name
                assembly_files.append(workflow.importFile(assembly_file_url))

            if options.no_duplicate_contig_ids == False:
                # if the user hasn't guaranteed that the contig ids are all unique, then rename 
                # any that are duplicates.
                #todo: make sure that the reference genome sequence ids are also included 
                #todo:      in the renaming. Ideally, make it so that the reference is 
                #todo:      evaluated first, which should ensure that the reference never
                #todo:      includes renamed sequence ids (assuming that the reference 
                #todo:      doesn't contain duplicate sequence ids internally)
                edited_assembly_files = workflow.start(Job.wrapJobFn(fasta_preprocessing.rename_duplicate_contig_ids, assembly_files, ref_id))
                
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
            if options.export_all_to_all_files:
                alignments, all_to_all_fastas, all_to_all_sams = workflow.start(Job.wrapJobFn(
                    align_all_assemblies,  ref_id, assembly_files, remap_stats_internal_file, options=options)) 
                print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~both lists:", all_to_all_fastas, all_to_all_sams)
            else:
                alignments = workflow.start(Job.wrapJobFn(
                    align_all_assemblies,  ref_id, assembly_files, remap_stats_internal_file, options=options))


            # reformat the alignments as lastz cigars:
            (lastz_cigar_primary_alignments, lastz_cigar_secondary_alignments) = workflow.start(Job.wrapJobFn(
                make_lastz_output, alignments))
                
            workflow.exportFile(lastz_cigar_primary_alignments, 'file://' + os.path.abspath(options.primary_output_file))
            workflow.exportFile(lastz_cigar_secondary_alignments, 'file://' + os.path.abspath(options.secondary_output_file))

            if options.remap_stats:
                workflow.exportFile(remap_stats_internal_file, 'file://' + os.path.abspath(options.remap_stats_output_file))

            if options.export_all_to_all_files:
                file_count = 0
                for fasta in all_to_all_fastas:
                    file_count += 1
                    workflow.exportFile(fasta, 'file://' + os.path.abspath(".") + "/all_to_all_fasta_" + str(file_count) + ".fa")

                file_count = 0
                for mapping_file in all_to_all_sams:
                    file_count += 1
                    (lastz_cigar_primary_alignments, lastz_cigar_secondary_alignments) = workflow.start(Job.wrapJobFn(
                make_lastz_output, mapping_file))
                    workflow.exportFile(lastz_cigar_primary_alignments, 'file://' + os.path.abspath(".") + "/all_to_all_mapping_file_primary" + str(file_count) + ".cigar")
                    workflow.exportFile(lastz_cigar_secondary_alignments, 'file://' + os.path.abspath(".") + "/all_to_all_mapping_file_secondary" + str(file_count) + ".cigar")
                    workflow.exportFile(mapping_file, 'file://' + os.path.abspath(".") + "/all_to_all_mapping_file_" + str(file_count) + ".sam")
        else:
            output = workflow.restart()
    
if __name__ == "__main__":
    main()
