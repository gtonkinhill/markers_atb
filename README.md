# Markers

Some simple scripts for extracting markers for the All the Bacteria project.

```
usage: aln_markers.py [-h] -m MODELS (-p PROTEINS | -a ASSEMBLY) -o OUTPUT
                      [-c] [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

Search HMM models against protein sequences and return best-hit alignments.
Each model's best hit (or gaps if no hit) is extracted and optionally
concatenated into a single FASTA output.

options:
  -h, --help            show this help message and exit
  -m MODELS, --models MODELS
                        Path to the HMM models file (HMMER3 format). Can
                        contain multiple HMMs.
  -p PROTEINS, --proteins PROTEINS
                        Path to the protein sequences FASTA file.
  -a ASSEMBLY, --assembly ASSEMBLY
                        Path to the nucleotide assembly FASTA file. Proteins
                        will be predicted using Pyrodigal.
  -o OUTPUT, --output OUTPUT
                        Path to the output FASTA file for the best-hit
                        alignments.
  -c, --concatenate     If set, concatenate all alignments into a single FASTA
                        entry. Default: False.
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging level. Default is WARNING.
```