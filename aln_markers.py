import argparse
import sys
import os
import io
import collections
import re
from pathlib import Path
import pyhmmer
import pyrodigal
from pyhmmer.plan7 import HMMFile, TopHits
from pyhmmer.easel import Alphabet, SequenceFile, MSAFile, SequenceBlock, TextSequence


def get_basename(file_path):
    """Gets base name of file"""
    base_name = os.path.basename(file_path)
    file_name_without_extension, _ = os.path.splitext(base_name)
    return file_name_without_extension

def main():
    """
    Main function to parse arguments, run HMMER searches, and generate output.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Search HMM models against protein sequences and return best-hit alignments.\n"
            "Each model's best hit (or gaps if no hit) is extracted and optionally concatenated "
            "into a single FASTA output."
        )
    )

    parser.add_argument(
        "-m", "--models",
        type=Path,
        required=True,
        help="Path to the HMM models file (HMMER3 format). Can contain multiple HMMs."
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-p", "--proteins",
        type=Path,
        help="Path to the protein sequences FASTA file."
    )
    group.add_argument(
        "-a", "--assembly",
        type=Path,
        help="Path to the nucleotide assembly FASTA file. Proteins will be predicted using Pyrodigal."
    )

    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Path to the output FASTA file for the best-hit alignments."
    )

    parser.add_argument(
        "-c", "--concatenate",
        action="store_true",
        help="If set, concatenate all alignments into a single FASTA entry. Default: False."
    )

    args = parser.parse_args()

    hmm_filepath = args.models
    output_filepath = args.output
    print("here")
    # --- Input Validation ---
    # Check if the provided HMM models file exists.
    if not os.path.exists(hmm_filepath):
        print(f"Error: HMM models file not found at '{hmm_filepath}'", file=sys.stderr)
        sys.exit(1)

    # List to store the best alignment (or hyphen string) for each HMM model.
    all_best_hit_alignments = []

    try:
        # --- Load HMM models ---
        # `HMMFile` is used to parse HMMER3 format HMM files.
        # `pyhmmer` typically infers the alphabet (e.g., ProteinAlphabet) automatically.
        print(f"Loading HMM models from {hmm_filepath}...")
        hmms = []
        with HMMFile(hmm_filepath) as hmm_file:
            for hmm in hmm_file:
                hmms.append(hmm)
        
        # Handle the case where no HMM models are found in the input file.
        if not hmms:
            print("Warning: No HMM models found in the provided file. Output will be an empty FASTA entry.", file=sys.stderr)
            sys.exit(0)

        # --- Load protein sequences ---
        if args.proteins:
            seq_filepath = args.proteins
            if not os.path.exists(seq_filepath):
                print(f"Error: Protein sequences file not found at '{seq_filepath}'", file=sys.stderr)
                sys.exit(1)
            print(f"Loading protein sequences from {seq_filepath}...")
            with SequenceFile(seq_filepath, digital=True) as seqs_file:
                sequences = list(seqs_file)

        elif args.assembly:
            seq_filepath = args.assembly
            if not os.path.exists(seq_filepath):
                print(f"Error: Assembly file not found at '{seq_filepath}'", file=sys.stderr)
                sys.exit(1)
            print(f"Predicting proteins from contigs using Pyrodigal: {seq_filepath}...")
            orf_finder = pyrodigal.GeneFinder(meta=True, closed=False)
            
            gene_count = 0
            with SequenceFile(seq_filepath, digital=False) as seqs_file:
                contigs = list(seqs_file)

            sequences = []
            for contig in contigs:
                orfs = orf_finder.find_genes(contig.sequence)
                for gene in orfs:
                    sequences.append(pyhmmer.easel.TextSequence(
                        name = str(gene_count).encode(), 
                        sequence = gene.translate()).digitize(
                            pyhmmer.easel.Alphabet.amino()))
                    gene_count += 1

            if not sequences:
                print("Warning: No proteins predicted from assembly.", file=sys.stderr)
                sys.exit(1)

        # Determine the header for the output FASTA file based on the input sequence file's name.
        output_header_prefix = get_basename(seq_filepath)

        # index sequences to speed things up later
        sequence_dict = {seq.name.decode(): i for i,seq in enumerate(sequences)}
        
        # Warn if no protein sequences were loaded.
        if not sequences:
            print("Warning: No protein sequences found in the provided FASTA file!!!", file=sys.stderr)
            sys.exit(1)
            
        # --- Perform searches for each HMM ---
        print("Performing HMM searches for each model...")
        concatenated_hits = ""
        for i, hmm in enumerate(hmms):
            # Decode HMM name from bytes to string for display purposes.
            # Use a fallback name if HMM name is not available.
            hmm_name = hmm.name.decode() if hmm.name else f"unnamed_HMM_{i+1}"
            # print(f"  Searching with HMM '{hmm_name}' ({i+1}/{len(hmms)})...")

            # Search current hmm against all sequences
            hits = next(pyhmmer.hmmsearch(hmm, sequences, bit_cutoffs="trusted", cpus=5))

            if len(hits) < 1:
                # No matches - output a string oh hypens
                current_hit = "-" * hmm.consensus_length
            else:
                # Obtain hit with the best score
                argmax, best_hit = max(enumerate(hits), key=lambda x: x[1].score)

                # Generate full alignment
                aln = pyhmmer.hmmer.hmmalign(hmm, [sequences[sequence_dict[best_hit.name.decode()]]], 
                                           trim=True, digitize=False, all_consensus_cols=True)
                    
                # Remove inserts to maintain the same length as the HMM model.
                current_hit = re.sub(r'[a-z]', '', aln.alignment[0])
                
                # print("HMM_length", best_hit.best_domain.alignment.hmm_length)
                # print("Aln length", len(current_hit))

                if len(current_hit) != best_hit.best_domain.alignment.hmm_length:
                    print("Warning: alignment length does not match HMM length!!!")
                    sys.exit(1)

                if args.concatenate:
                   concatenated_hits += current_hit
                else:
                    print(f">{hmm_name}_{output_header_prefix}")
                    print(current_hit)

        if args.concatenate:
            print(f">{output_header_prefix}")
            print(concatenated_hits)

        print("Script finished successfully.")

    except FileNotFoundError as e:
        # Catch and report specific file not found errors.
        print(f"Error: A required file was not found: {e}. Please verify the paths.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        # Catch any other unexpected errors during script execution.
        print(f"An unexpected error occurred during execution: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    # Ensure the main function is called when the script is executed.
    main()