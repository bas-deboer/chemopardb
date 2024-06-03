from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import Protein
from structure.models import Structure

import os
import pandas as pd
import re
import urllib.request
import csv
import os
import urllib
import json
from Bio import SeqIO
from Bio import PDB
from Bio.PDB import Polypeptide, PDBParser, Superimposer
from Bio.PDB.PDBIO import PDBIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
from Bio import PDB
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser, PDBIO, PDBList
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1
from Bio import Align
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def parse_pdb(file_path):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", file_path)
    pdb_residue_positions = []
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_position = residue.id[1]
                pdb_residue_positions.append(residue_position)
    return pdb_residue_positions

def parse_sequence(sequence_files, target_sequence_name):
    sequence_positions = {}
    sequence_file_path = os.path.join(sequence_files, f"{target_sequence_name}_aligned.fasta")
    sequence = SeqIO.read(sequence_file_path, "fasta")
    # Extract positions of non-gap residues in the target sequence
    sequence_positions = [i + 1 for i, residue in enumerate(sequence.seq) if residue != "-"]
    return sequence_positions

def extract_sequence_from_chain(chain):
    sequence = ""
    for residue in chain:
        # Check if the residue has an amino acid code (to exclude water and ligands)
        if is_aa(residue):
            # Get the three-letter amino acid code
            three_letter_code = residue.get_resname()
            # Convert the three-letter code to one-letter code using the dictionary
            one_letter_code = protein_letters_3to1.get(three_letter_code, 'X')  # 'X' for unknown
            # Append the one-letter amino acid code to the sequence
            sequence += one_letter_code
    return sequence

def find_closest(ref_seqs, query_id, query_seq):
    best_score = 0
    best_id = None
    best_alignment = None
    for ref_id, ref_seq in ref_seqs.items():
        pw = pairwise2.align.localms(ref_seq, query_seq, 3, 1, -3, -.1)
        score = pw[0][2]
        if score > best_score:
            best_score = score
            best_id = ref_id
            best_alignment = pw
            
    return best_id, best_score, best_alignment

def align_to_ref(ref_seq, query_seq, ident_score=4, sim_score=2, gap_open=-2, gap_ext=-0.5, verbose=False):
    pw = pairwise2.align.localms(ref_seq, query_seq, ident_score, sim_score, gap_open, gap_ext)
    score = pw[0][2]
    if verbose:
        print(format_alignment(*pw[0]))
        print(score)
        print(len(ref_seq), len(pw[0][1]))  # Use len(ref_seq) to get the length of the reference sequence
        
    return pw[0][1]  # Return the aligned sequence


class Command(BaseBuild):

    def add_arguments(self, parser):
        parser.add_argument("--debug", default=False, action="store_true", help="Debug mode")
        parser.add_argument("--purge", default=False, action="store_true", help="Purge data")

    def handle(self, *args, **options):
        #if options['purge']:
        #    PDBData.objects.all().delete()
        #    Publication.objects.all().delete()
        #    StructureDomain.objects.all().delete()
        #    Chain.objects.all().delete()
        #    Structure.objects.all().delete()
        #    StructureType.objects.all().delete()

        # Align sequence to MSA
        self.align_model()
        
        # Creating mapping
        self.create_mapping()

    def align_model(self):
        input_dir = os.sep.join([settings.DATA_DIR, 'all_chains'])
        pdb_files = [os.path.join(input_dir, filename) for filename in os.listdir(input_dir) if filename.endswith(".pdb")]

        for pdb_file in pdb_files:
            parser = PDB.PDBParser(QUIET=True)

            # Parse the PDB file
            structure = parser.get_structure("structure", pdb_file)
            pdb_base_name = os.path.splitext(os.path.basename(pdb_file))[0]  # Extract only the file name

            # Create a directory to store sequence files for this PDB file
            output_dir = os.sep.join([settings.DATA_DIR, 'sequences'])
            os.makedirs(output_dir, exist_ok=True)

            for chain in structure.get_chains():
                protein_sequence = extract_sequence_from_chain(chain)
                # Construct the header for the FASTA entry
                header = f">{pdb_base_name}"

                # Construct the output file path (inside the output directory)
                output_file_name = f"{pdb_base_name}_{chain.id}.fasta"  # Include chain ID in the filename
                output_file_path = os.path.join(output_dir, output_file_name)

                # Save the protein sequence to a FASTA file
                with open(output_file_path, "w") as output_file:
                    output_file.write(f"{header}\n{protein_sequence}\n")

        aligner = Align.PairwiseAligner()
        msa = AlignIO.read(os.sep.join([settings.DATA_DIR, 'chemokines_ref_MSA.aln']), "clustal")

        # Directory where your sequence files are located
        sequence_directory = os.sep.join([settings.DATA_DIR, 'sequences'])

        # Create a dictionary to store the reference sequences
        ref_seqs = {}

        # Populate the reference sequences dictionary from the MSA
        for record in msa:
            ref_seqs[record.id] = str(record.seq)

        # Create a list to store aligned sequences
        aligned_sequences = []

        # Output directory where aligned sequences will be saved
        output_directory = os.sep.join([settings.DATA_DIR, 'aligned_sequences'])
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        # Loop over your query sequences
        for file in os.listdir(sequence_directory):
            if file.endswith(".fasta"):
                file_path = os.path.join(sequence_directory, file)
                for record in SeqIO.parse(file_path, "fasta"):
                    query_id = record.id
                    query_seq = str(record.seq)
                    header = record.id
        
                    # Construct the output file path
                    output_filename = os.path.join(output_directory, f"{header}_aligned.fasta")
                    
                    # Check if the output file already exists
                    if os.path.exists(output_filename):
                        print(f"Sequence '{header}' already exists in the output directory. Skipping.")
                    else:
                        # Find best match from the reference MSA
                        best_match_id, best_match_score, best_match_alignment = find_closest(ref_seqs, query_id, query_seq)
                        
                        # Align the query sequence to the best match
                        aligned_seq = align_to_ref(ref_seqs[best_match_id], query_seq)
                        
                        # Create a SeqRecord for the aligned sequence
                        aligned_record = SeqRecord(Seq(aligned_seq), id=query_id, description=f"Aligned to {best_match_id}")
                        
                        # Append the aligned record to the list
                        aligned_sequences.append(aligned_record)
                        print(f"Sequence '{header}' aligned to {best_match_id} with score of {best_match_score}")

                        # Save the aligned sequence to a new fasta file
                        header = aligned_record.id
                        output_filename = os.path.join(output_directory, f"{header}_aligned.fasta")

                        # Check if the output file already exists (this is just an additional check)
                        if os.path.exists(output_filename):
                            print(f"Sequence '{header}' already exists in the output directory. Skipping.")
                        else:
                            SeqIO.write([aligned_record], output_filename, "fasta")
                            print(f"Aligned sequence '{header}' saved")


    def create_mapping(self):
        pdb_directory = os.sep.join([settings.DATA_DIR, 'all_chains'])
        sequence_files = os.sep.join([settings.DATA_DIR, 'aligned_sequences'])

        # Create a directory to save mapping files
        output_directory = os.sep.join([settings.DATA_DIR, 'MSA_mapping'])
        os.makedirs(output_directory, exist_ok=True)

        # List all PDB files in the directory
        pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith(".pdb")]
        for pdb_file in pdb_files:
            pdb_file_path = os.path.join(pdb_directory, pdb_file)
            residue_positions = parse_pdb(pdb_file_path)
            target_sequence_name = pdb_file.split('.')[0]

            sequence_positions = parse_sequence(sequence_files, target_sequence_name)

            # Create the mapping between residue_positions and sequence_positions as a dictionary
            mapping = {key: value for key, value in zip(residue_positions, sequence_positions)}

            # Define the output file path and save mapping as json file
            output_file_path = os.path.join(output_directory, f"{target_sequence_name}_mapping.json")
            with open(output_file_path, "w") as output_file:
                json.dump(mapping, output_file)