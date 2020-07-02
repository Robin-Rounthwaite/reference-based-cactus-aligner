from Bio import SeqIO

def rename_duplicate_contig_ids(job, assembly_files, reference_file):
    """
    Sometimes, when combining assemblies from multiple sources, multiple contigs get the 
    same name. This function slightly modifies all but one of the contigs with the same
    name to ensure that there are no duplicates. Renamed contigs are in a format that
    should be easy to reverse.

    Given a list of assembly files (in job.fileStore.writeGlobalFile id output),
    outputs the list of edited assembly files, with the only 
    difference being that all the contigs have been given unique names. Unique names 
    follow this formula:
    x = original contig id
    y = unique_integer
    new id = x_renamed_y
    """
    
    contig_ids = set()
    unique_id = int()

    #first, record the sequence ids in reference. (It is assumed that the reference 
    # doesn't contain duplicate ids internally)
    reference_file = job.fileStore.readGlobalFile(reference_file)
    reference_contigs = SeqIO.parse(reference_file, "fasta")
    for seq in reference_contigs:
        contig_ids.add(seq.id)

    for assembly in assembly_files:
        assembly_file = job.fileStore.readGlobalFile(assembly)
        assembly_contigs = SeqIO.parse(assembly_file, "fasta")
        output_contigs = list()
        
        for contig in assembly_contigs:
            
            if contig.id in contig_ids:
                old_id = contig.id
                
                while contig.id in contig_ids:
                    # then there is a duplicate contig_id. edit this one.
                    # keep changing the contig_id until we get a completely unique id.
                    contig.id = old_id + "_renamed_" + str(unique_id)
                    contig.description = old_id + "_renamed_" + str(unique_id)
                    unique_id += 1
                    
                #record the new contig id as an observed id.
                contig_ids.add(contig.id)
                
            else:
                # this isn't a duplicate contig_id. record it.
                contig_ids.add(contig.id)

            output_contigs.append(contig)

        # write the altered assembly to original file.    
        SeqIO.write(output_contigs, assembly_file, "fasta")

    # return the file ids of the original, overwritten assembly files.
    return assembly_files
