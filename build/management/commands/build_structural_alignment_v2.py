from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import Protein
from structure.models import Structure

import os
import numpy as np
import pandas as pd
import re
import urllib.request
import csv
import os
import pprint
import urllib
import json
from datetime import datetime
from Bio import PDB
from Bio.PDB import Polypeptide, PDBParser, Superimposer
from Bio.PDB.PDBIO import PDBIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser, PDBIO, PDBList
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1
from Bio import Align
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment















def apply_transformation_to_structure(input_pdb_path, transformation_data, output_pdb_path):
    # Extract rotation and translation matrices from the transformation_data dictionary
    rotation_matrix = transformation_data["rotation"]
    translation_vector = transformation_data["translation"]

    # Load the structure you want to apply the transformation to
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("TransformedStructure", input_pdb_path)

    # Apply the rotation matrix and translation vector to all atoms in the structure
    for atom in structure.get_atoms():
        atom.transform(rotation_matrix, translation_vector)

    # Save the superimposed structure to the output PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_path)

def superimpose_to_pdb_raw_directory(pdb_raw_directory, json_directory, output_pdb_directory):
    # List all JSON mapping files in the JSON directory
    json_files = [f for f in os.listdir(json_directory) if f.endswith(".json")]

    # Iterate over each JSON mapping file
    for json_file in json_files:
        # Extract the base name without the file extension and split it at underscores
        json_base_name = os.path.splitext(json_file)[0].split("_")[0]
        json_name = os.path.splitext(json_file)[0].split(".")[0]

        # Construct the corresponding PDB Raw file path
        pdb_raw_path = os.path.join(pdb_raw_directory, f"{json_base_name}.pdb")

        # Check if the PDB Raw file exists
        if os.path.exists(pdb_raw_path):
            # Construct the output PDB file path based on the JSON mapping file name
            output_pdb_path = os.path.join(output_pdb_directory, f"{json_name}_superimposed.pdb")

            # Load the JSON data
            with open(os.path.join(json_directory, json_file), "r") as f:
                transformation_data = json.load(f)

            # Apply the transformation and save the superimposed structure using your function
            apply_transformation_to_structure(pdb_raw_path, transformation_data, output_pdb_path)
        else:
            print(f'PDB file not found: {pdb_raw_path}')


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

        # superimpose model to reference structure
        self.superimpose()
        
        # superimpose original pdb
        self.superimpose_original()

    def superimpose(self):
        # function to superimpose model to reference
        pdb_directory = os.sep.join([settings.DATA_DIR, 'all_chains'])
        mapping_directory = os.sep.join([settings.DATA_DIR, 'MSA_mapping'])
        output_directory = os.sep.join([settings.DATA_DIR, 'aligned_models'])
        transformation_dir = os.sep.join([settings.DATA_DIR, 'transformation_matrices'])
        os.makedirs(output_directory, exist_ok=True)
        os.makedirs(transformation_dir, exist_ok=True)

        template_structure = os.sep.join([settings.DATA_DIR, 'structural_ref.pdb'])
        core_residues = [65, 77, 78, 79, 83, 84, 85, 86, 87, 101, 102, 103, 104, 105, 106, 115, 116, 117, 118]

        for pdb_file in os.listdir(pdb_directory):
            if pdb_file.endswith(".pdb"):
                pdb_path = os.path.join(pdb_directory, pdb_file)

                parser = PDB.PDBParser(QUIET=True)
                ref_model = parser.get_structure("structure1", template_structure)
                sample_model = parser.get_structure("structure2", pdb_path)

                # Load the corresponding JSON file for the current PDB file
                json_filename = os.path.splitext(pdb_file)[0] + "_mapping.json"
                json_path = os.path.join(mapping_directory, json_filename)

                # Check if the JSON file exists
                if os.path.exists(json_path):
                    # Load the JSON data from the file into a Python dictionary
                    with open(json_path, 'r') as json_file:
                        residue_mapping = json.load(json_file)

                    # Filter the residue_mapping dictionary for specific key values
                    residue_mapping = {key: value for key, value in residue_mapping.items() if value in core_residues}

                    # Extract template and target selections
                    reference_residues = list(residue_mapping.values())
                    target_residues = list(residue_mapping.keys())

                    target_residues = [int(x) for x in target_residues]

                    # Create lists of atoms to match
                    fixed_atoms = []
                    moving_atoms = []

                    # Iterate through the residues and atoms in the reference structure
                    for model in ref_model:
                        for chain in model:
                            for residue in chain:
                                if residue.get_id()[1] in reference_residues:  # Check residue number
                                    for atom in residue:
                                        if atom.get_name() in ['CA']:
                                            # Add the atom to the list of backbone atoms
                                            fixed_atoms.append(atom)

                    # Iterate through the residues and atoms in the sample structure
                    for model in sample_model:
                        for chain in model:
                            for residue in chain:
                                if residue.get_id()[1] in target_residues:  # Check residue number
                                    for atom in residue:
                                        if atom.get_name() in ['CA']:
                                            # Add the atom to the list of backbone atoms
                                            moving_atoms.append(atom)

                    # Check if either of the lists is empty
                    if not fixed_atoms or not moving_atoms:
                        print(f"Skipping alignment for {pdb_file} due to missing atoms.")
                    else:
                        # Create a Superimposer object
                        sup = Superimposer()
    
                        # Set the atom lists in the Superimposer
                        sup.set_atoms(fixed_atoms, moving_atoms)
    
                        # Perform superposition and obtain rotation/translation matrix
                        sup.apply(sample_model.get_atoms())
    
                        # Create a dictionary to store the transformation matrix
                        transformation_matrix = {
                            "rotation": sup.rotran[0].tolist(),
                            "translation": sup.rotran[1].tolist()
                        }
    
                        print(f"RMSD for {pdb_file}:")
                        print(sup.rms)
    
                        # Save the aligned structure with a unique filename
                        output_filename = os.path.join(output_directory, f"{pdb_file}")
                        io = PDBIO()
                        io.set_structure(sample_model)
                        io.save(output_filename)
                        
                        # Save the transformation matrix to a JSON file
                        matrix_filename = os.path.join(transformation_dir, f"{pdb_file}.json")
                        with open(matrix_filename, 'w') as matrix_file:
                            json.dump(transformation_matrix, matrix_file)
                    

    def superimpose_original(self):
        # Specify the input directory containing PDB files, JSON directory, and
        # the output directory to save the superimposed PDB files.
        pdb_raw_directory = os.sep.join([settings.DATA_DIR, 'raw_pdbs'])
        json_directory = os.sep.join([settings.DATA_DIR, 'transformation_matrices'])
        output_pdb_directory = os.sep.join([settings.DATA_DIR, 'aligned_fullstructures'])

        # Create the output directory if it doesn't exist
        os.makedirs(output_pdb_directory, exist_ok=True)

        # Perform superimposition for each JSON file using the corresponding PDB files
        superimpose_to_pdb_raw_directory(pdb_raw_directory, json_directory, output_pdb_directory)
