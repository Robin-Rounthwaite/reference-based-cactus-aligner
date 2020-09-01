from toil.common import Toil
from toil.job import Job

import subprocess
import os
from argparse import ArgumentParser

from src import paf_to_lastz
from src import fasta_preprocessing


## utilitary fxns:

def unpack_promise(job, iterable, i):
    """
    passed an iterable and a location i, returns ith item.
    """
    return iterable[i]

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
                for line in f1:
                    if not line.startswith("@"):
                        outfile.write(line)
    return job.fileStore.writeGlobalFile(consolidated_mappings)

def get_asms_from_seqfile(workflow, seqfile):
    asm_files = dict()

    with open(seqfile) as inf:
        #skip the Newick tree on the first line:
        next(inf)
        
        for line in inf:
            parsed = line.split()
            if len(parsed) >= 2:
                #note: fastas can be in a directory containing single consecutive spaces. Two spaces in a row breaks my file parsing system. 
                fasta_url = 'file://' + os.path.abspath(" ".join(parsed[1:]))
                asm_files[parsed[0]] = workflow.importFile(fasta_url)
    
    return asm_files

def import_asms(job, options, workflow):
    """Import asms; deduplicating contig ids if not --all_unique_ids

    Args:
        seqfile ([type]): [description]
        workflow ([type]): [description]

    Returns:
        [type]: [description]
    """
    # asms is dictionary of all asms (not counting reference) with key: asm_name, value: imported global toil file.
    asms = get_asms_from_seqfile(workflow, options.seqFile)
    refFile_location = asms.pop(options.refFile)
    print(refFile_location)
    if not options.all_unique_ids:
        # deduplicate contig id names, if user hasn't guaranteed unique contig ids.
        # new_fastas is the location of the asms with unique ids.
        if options.overwrite_assemblies:
            # overwrite the original assemblies. Note that the reference is never overwritten, as it is never altered. 
            # (Code assumes that the reference is internally free of duplicate ids, and just ensures other asms don't use reference ids.)
            asms = fasta_preprocessing.rename_duplicate_contig_ids(job, asms, refFile_location, asms, workflow)
        else:
            # don't overwrite the original assemblies.
            # first, determine the new asm save locations.
            if not os.path.isdir(options.assembly_save_dir):
                os.mkdir(options.assembly_save_dir)

            new_asms = dict()
            for asm_id, asm in asms.items():
                new_asms[asm_id] = options.assembly_save_dir + asm.split("/")[-1]
            
            asms = fasta_preprocessing.rename_duplicate_contig_ids(job, asms, refFile_location, new_asms, workflow)

    # Import asms.
    for asm_id, asm in asms:
        asms[asm_id] = workflow.importFile('file://' + os.path.abspath(asm))

    return asms

def empty(job):
    """
    An empty job, for easier toil job organization.
    """
    return

## mapping fxns:

def map_all_to_ref(job, assembly_files, reference_file):
    """
    Primarily for use with option_all_to_ref_only. Otherwise, use map_all_to_ref_and_get_poor_mappings.
    """
    lead_job = job.addChildJobFn(empty)

    ref_mappings = dict()
    for assembly_file in assembly_files:
        ref_mappings[assembly_file] = lead_job.addChildJobFn(mapping_functions.map_a_to_b, assembly_file, reference_file).rv()
    
    consolidate_job = lead_job.addFollowOnJobFn(consolidate_mappings, ref_mappings)
    paf_mappings = consolidate_job.rv()

    conversion_job = consolidate_job.addFollowOnJobFn(paf_to_lastz.paf_to_lastz, paf_mappings)
    lastz_mappings = conversion_job.rv()

    primary_mappings = conversion_job.addChildJobFn(unpack_promise, lastz_mappings, 0)
    secondary_mappings = conversion_job.addChildJobFn(unpack_promise, lastz_mappings, 1)

    return (primary_mappings, secondary_mappings)

def map_a_to_b(job, a, b):
    """Maps fasta a to fasta b.

    Args:
        a (global file): fasta file a. In map_all_to_ref, a is an assembly fasta.
        b (global file): fasta file b. In map_all_to_ref, b is the reference.

    Returns:
        [type]: [description]
    """
    
    map_to_ref_paf = job.fileStore.writeGlobalFile(job.fileStore.getLocalTempFile())

    subprocess.call(["minimap2", "-cx", "asm5", "-o", job.fileStore.readGlobalFile(map_to_ref_paf),
                    job.fileStore.readGlobalFile(reference_file), job.fileStore.readGlobalFile(assembly_to_align_file)])
     
    return map_to_ref_paf


## main fxn and interface:

def get_options():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    
    # options for basic input/output
    parser.add_argument('seqFile', type=str,
                        help='A file containing all the information specified by cactus in construction. This aligner ignores the newick tree.')
    parser.add_argument('refFile', type=str, 
                        help='Specifies which asm in seqFile should be treated as the reference.')
    parser.add_argument('--primary', default="primary.cigar", type=str, 
                        help='Filename for where to write lastz cigar output for primary mappings.')
    parser.add_argument('--secondary', default="secondary.cigar", type=str, 
                        help='Filename for where to write lastz cigar output for secondary mappings.')
                        
    # options for importing assemblies:
    parser.add_argument('--all_unique_ids', action='store_true', 
                        help="Don't clean the assembly files; the user promises that they don't contain any duplicate contig ids.")
    parser.add_argument('--overwrite_assemblies', action='store_true', 
                        help="When cleaning the assembly files to make sure there are no duplicate contig ids, don't overwrite the assembly files. Copy them to a neigboring folder with the affix '_edited_for_duplicate_contig_ids' instead.")
    parser.add_argument('--assembly_save_dir', type=str, default='./unique_id_assemblies/',
                        help='While deduplicating contig ids in the input fastas, save the assemblies in this directory. Ignored when used in conjunction with --overwrite_assemblies.')
                        
    options = parser.parse_args()
    return options

def main():
    options = get_options()

    with Toil(options) as workflow:
        ## Preprocessing:
        # Import asms; deduplicating contig ids if not --all_unique_ids
        asms = workflow.start(Job.wrapJobFn(import_asms, options, workflow))
            
        # Import reference:
        reference = workflow.importFile('file://' + os.path.abspath(options.refFile))

        ## Perform alignments:
        if not workflow.options.restart:
            alignments = workflow.start(Job.wrapJobFn(map_all_to_ref, asms, reference))

        else:
            alignments = workflow.restart()

        ## Save alignments:
        workflow.exportFile(alignments[0], 'file://' + os.path.abspath(options.primary))
        workflow.exportFile(alignments[1], 'file://' + os.path.abspath(options.secondary))

if __name__ == "__main__":
    main()