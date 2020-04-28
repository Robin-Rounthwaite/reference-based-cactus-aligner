import pytest

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '..')
import minimap2_cactus_aligner_toil as cactus_aligner

from toil.common import Toil
from toil.job import Job
from argparse import ArgumentParser
import os
import shutil
from Bio import SeqIO
import collections as col


def get_options():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    # options.test_fasta = "my_test_fasta.fa"
    # options.test_sam = "my_test_sam.fa"
    # options.mapq_cutoff = 20
    # options.sequence_context = 5 # currently normally 10000
    options = parser.parse_args()
    
    return options

def make_assembly_files_for_test_rename_duplicate_contig_ids(test_dir):
    test_fa1 = test_dir + "fa1.fa"
    with open(test_fa1, "w") as f:
        f.write(
""">1
AATTAACC
>1
AATTCCGG
>2
AATTCCGA
>3
AATTGCGA"""
                )

    test_fa2 = test_dir + "fa2.fa"
    with open(test_fa2, "w") as f:
        f.write(
""">1
AATTAACC
>2
AATTCCGA
>4
AATTGCGA"""
                )

    return [os.path.abspath(test_fa1), os.path.abspath(test_fa2)]

def test_rename_duplicate_contig_ids_overwrite():
    """
    Note: test_tmp is a directory reserved by this function. If you have that directory in
    the test directory, it will be deleted.
    """
    
    options = get_options()

    # make the directories for testing:
    if os.path.isdir("./test_tmp/"):
        shutil.rmtree("./test_tmp/")
    os.mkdir("./test_tmp/")
    os.mkdir("./test_tmp/test_rename_duplicate_contig_ids/")
    file_dir = "./test_tmp/test_rename_duplicate_contig_ids/overwrite_fastas/"
    os.mkdir(file_dir)

    # make the files for testing
    assembly_files = make_assembly_files_for_test_rename_duplicate_contig_ids(file_dir)
        
    with Toil(options) as workflow:
        # import the files into the workflow
        assembly_ids = list()
        for assembly_file in assembly_files:
            assembly_ids.append(workflow.importFile('file://' + assembly_file))

        # run the function being tested
        new_assembly_ids = workflow.start(Job.wrapJobFn(cactus_aligner.rename_duplicate_contig_ids, assembly_ids))

        # export the files from the workflow
        for i in range(len(new_assembly_ids)):
            # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",assembly_files[i])
            # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++", new_assembly_ids[i])
            workflow.exportFile(new_assembly_ids[i], 'file://' + assembly_files[i])
            # output_file = assembly_files[i]
            # workflow.exportFile(new_assembly_ids[i], "file://" + output_file)
    
    # print("new_assembly_files", new_assembly_files)
    id_counts = col.Counter()
    for assembly in assembly_files:
        # print(os.path.abspath(assembly))
        assembly_contigs = SeqIO.parse(assembly, "fasta")
        for contig in assembly_contigs:
            # print("contig.id", contig.id)
            id_counts[contig.id] += 1
    
    for contig_id, count in id_counts.items():
        assert count == 1, "contig id %r exists still exists multiple (%r) times in fasta files" % (contig_id, count)
        
    shutil.rmtree("./test_tmp/")

def test_rename_duplicate_contig_ids_no_overwrite():
    """
    Note: test_tmp is a directory reserved by this function. If you have that directory in
    the test directory, it will be deleted.
    """
    
    options = get_options()

    # make the directories for testing:
    if os.path.isdir("./test_tmp/"):
        shutil.rmtree("./test_tmp/")
    os.mkdir("./test_tmp/")
    os.mkdir("./test_tmp/test_rename_duplicate_contig_ids/")
    input_file_dir = "./test_tmp/test_rename_duplicate_contig_ids/input_fastas/"
    os.mkdir(input_file_dir)
    output_file_dir = "./test_tmp/test_rename_duplicate_contig_ids/output_fastas/"
    os.mkdir(output_file_dir)

    # make the files for testing
    assembly_files = make_assembly_files_for_test_rename_duplicate_contig_ids(input_file_dir)
        
    # new_assembly_files will be the function output.
    new_assembly_files = list()

    with Toil(options) as workflow:
        # import the files into the workflow
        assembly_ids = list()
        for assembly_file in assembly_files:
            assembly_ids.append(workflow.importFile('file://' + assembly_file))

        # run the function being tested
        new_assembly_ids = workflow.start(Job.wrapJobFn(cactus_aligner.rename_duplicate_contig_ids, assembly_ids, overwrite=False))

        # export the files from the workflow
        for i in range(len(new_assembly_ids)):
            output_file = os.path.abspath(output_file_dir + "fa_" + str(i) + ".fa")
            workflow.exportFile(new_assembly_ids[i], "file://" + output_file)
            new_assembly_files.append(output_file)
    
    # print("new_assembly_files", new_assembly_files)
    id_counts = col.Counter()
    for new_assembly in new_assembly_files:
        # print(os.path.abspath(new_assembly))
        assembly_contigs = SeqIO.parse(new_assembly, "fasta")
        for contig in assembly_contigs:
            # print("contig.id", contig.id)
            id_counts[contig.id] += 1
    
    for contig_id, count in id_counts.items():
        assert count == 1, "contig id %r exists still exists multiple (%r) times in fasta files" % (contig_id, count)
        
    shutil.rmtree("./test_tmp/")




def main():
        test_rename_duplicate_contig_ids_overwrite()

if __name__ == "__main__":
    main()
    