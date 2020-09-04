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

def consolidate_mappings(job, mapping_files):
    """
    Warning: discards headers of all mapping files.
    Given a list of mapping files, consolidates the contents (not counting headers) into a
    single file.
    """
    consolidated_mappings = job.fileStore.getLocalTempFile()
    with open(consolidated_mappings,"w") as outfile:
        for mapping_file in mapping_files.values():
            with open(job.fileStore.readGlobalFile(mapping_file)) as inf:
                for line in inf:
                    if not line.startswith("@"):
                        outfile.write(line)
    return job.fileStore.writeGlobalFile(consolidated_mappings)

def get_asms_from_seqfile(seqfile):
    asm_files = dict()

    with open(seqfile) as inf:
        #skip the Newick tree on the first line:
        next(inf)
        
        for line in inf:
            parsed = line.split()
            if len(parsed) >= 2:
                #note: fastas can be in a directory containing single consecutive spaces. Two spaces in a row breaks my file parsing system. 
                if parsed[0][0] == "*":
                    asm_files[parsed[0][1:]] = " ".join(parsed[1:])
                else:
                    asm_files[parsed[0]] = " ".join(parsed[1:])

    return asm_files

def import_asms(options, workflow):
    """Import asms; deduplicating contig ids if not --all_unique_ids

    Args:
        seqfile ([type]): [description]
        workflow ([type]): [description]

    Returns:
        [type]: [description]
    """
    # asms is dictionary of all asms (not counting reference) with key: asm_name, value: imported global toil file.
    asms = get_asms_from_seqfile(options.seqFile)

    if not options.all_unique_ids:
        # deduplicate contig id names, if user hasn't guaranteed unique contig ids.
        # new_fastas is the location of the asms with unique ids.
        if options.overwrite_assemblies:
            # overwrite the original assemblies. Note that the reference is never overwritten, as it is never altered. 
            # (Code assumes that the reference is internally free of duplicate ids, and just ensures other asms don't use reference ids.)
            asms = fasta_preprocessing.rename_duplicate_contig_ids(asms, options.refID, asms)
        else:
            # don't overwrite the original assemblies.
            # first, determine the new asm save locations.
            if not os.path.isdir(options.assembly_save_dir):
                os.mkdir(options.assembly_save_dir)

            new_asms = dict()
            for asm_id, asm in asms.items():
                if asm_id != options.refID:
                    new_asms[asm_id] = options.assembly_save_dir + asm.split("/")[-1]
                else:
                    # reference file never needs deduplication of contig ids, since ref contig ids are counted before all other asms.
                    new_asms[asm_id] = asm
                    
            asms = fasta_preprocessing.rename_duplicate_contig_ids(asms, options.refID, new_asms)

    # Import asms.
    for asm_id, asm in asms.items():
        asms[asm_id] = workflow.importFile('file://' + os.path.abspath(asm))

    return asms

def empty(job):
    """
    An empty job, for easier toil job organization.
    """
    return

## mapping fxns:

def map_all_to_ref(job, assembly_files, reference, debug_export):
    """
    Primarily for use with option_all_to_ref_only. Otherwise, use map_all_to_ref_and_get_poor_mappings.
    """
    lead_job = job.addChildJobFn(empty)

    # map all assemblies to the reference. Don't map reference to reference, though.
    ref_mappings = dict()
    for assembly, assembly_file in assembly_files.items():
        if assembly != reference:
            ref_mappings[assembly] = lead_job.addChildJobFn(map_a_to_b, assembly_file, assembly_files[reference]).rv()
    
    consolidate_job = lead_job.addFollowOnJobFn(consolidate_mappings, ref_mappings)
    paf_mappings = consolidate_job.rv()

    conversion_job = consolidate_job.addFollowOnJobFn(paf_to_lastz.paf_to_lastz, paf_mappings)
    lastz_mappings = conversion_job.rv()

    primary_mappings = conversion_job.addChildJobFn(unpack_promise, lastz_mappings, 0).rv()
    secondary_mappings = conversion_job.addChildJobFn(unpack_promise, lastz_mappings, 1).rv()

    if debug_export:
        return (primary_mappings, secondary_mappings, ref_mappings, paf_mappings )
    else:
        return (primary_mappings, secondary_mappings)

def map_a_to_b(job, a, b):
    """Maps fasta a to fasta b.

    Args:
        a (global file): fasta file a. In map_all_to_ref, a is an assembly fasta.
        b (global file): fasta file b. In map_all_to_ref, b is the reference.

    Returns:
        [type]: [description]
    """
    
    # map_to_ref_paf = job.fileStore.writeGlobalFile(job.fileStore.getLocalTempFile())
    tmp = job.fileStore.getLocalTempFile()
    map_to_ref_paf = job.fileStore.writeGlobalFile(tmp)


    subprocess.call(["minimap2", "-cx", "asm5", "-o", job.fileStore.readGlobalFile(map_to_ref_paf),
                    job.fileStore.readGlobalFile(b), job.fileStore.readGlobalFile(a)])
    
    return map_to_ref_paf


## main fxn and interface:

def get_options():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    
    # options for basic input/output
    parser.add_argument('seqFile', type=str,
                        help='A file containing all the information specified by cactus in construction. This aligner ignores the newick tree.')
    parser.add_argument('refID', type=str, 
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
                        
    # for debugging:
    parser.add_argument('--debug_export', action='store_true',
                        help='Export several other files for debugging inspection.')


    options = parser.parse_args()
    return options

def main():
    options = get_options()

    with Toil(options) as workflow:
        ## Preprocessing:
        # Import asms; deduplicating contig ids if not --all_unique_ids
        asms = import_asms(options, workflow)
            
        ## Perform alignments:
        if not workflow.options.restart:
            alignments = workflow.start(Job.wrapJobFn(map_all_to_ref, asms, options.refID, options.debug_export))

        else:
            alignments = workflow.restart()

        if options.debug_export:
            print(alignments)
            # Then return value is: (primary_mappings, secondary_mappings, ref_mappings, paf_mappings )
            for asm, mapping_file in alignments[2].items():
                workflow.exportFile(mapping_file, 'file://' + os.path.abspath("debug_" + asm + "_mapping_to_ref.txt"))
                break
            workflow.exportFile(alignments[3], 'file://' + os.path.abspath("debug_paf_mappings.txt"))

        ## Save alignments:
        workflow.exportFile(alignments[0], 'file://' + os.path.abspath(options.primary))
        workflow.exportFile(alignments[1], 'file://' + os.path.abspath(options.secondary))

if __name__ == "__main__":
    main()