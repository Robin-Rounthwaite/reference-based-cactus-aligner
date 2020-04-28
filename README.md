# reference-based-cactus-constructor
Faster than an all-to-all alignment for constructing cactus graphs. Starts by aligning all sequences to the reference. All the sequences that align poorly (or not at all) to the reference are then re-aligned to all the non-reference sequence.

## Install
### To clone:
git clone https://github.com/Robin-Rounthwaite/reference-based-cactus-constructor.git
### Python prerequisites:
pip install toil[aws,google,htcondor,encryption,cwl,wdl] biopython argparse cigar 
