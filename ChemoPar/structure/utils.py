import os
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser


def csv_to_matrix(file_path):
    """
    Reads a CSV file and returns a binary matrix representing the presence of interactions.
    """
    data = pd.read_csv(file_path[0], header=None)
    interaction_details = data.iloc[0]
    all_residues = set(range(1, 71))
    binary_matrix = np.zeros((1, len(all_residues)))  # Adjust based on your data if needed

    for i, detail in enumerate(interaction_details):
        parts = detail.split('.')
        for part in parts:
            if part.startswith('WAT'):
                continue
            if part[:3].isalpha() and part[3:].isdigit():
                residue_num = int(part[3:])
                if residue_num in all_residues:
                    col_index = sorted(list(all_residues)).index(residue_num)
                    binary_matrix[0, col_index] = 1
    return binary_matrix

def calc_distance_matrix(pdb_file):
    # Initialize the parser
    parser = PDBParser()
    
    # Read the structure
    structure = parser.get_structure('structure_name', pdb_file)
    ca_atoms = []
    
    # Extract all C alpha atoms
    ca_atoms = [atom for atom in structure.get_atoms() if atom.get_name() == 'CA']
    
    # Number of C alpha atoms
    n = len(ca_atoms)
    
    # Initialize the distance matrix
    distance_matrix = np.zeros((n, n))
    
    # Calculate distances and fill the matrix
    for i, atom_i in enumerate(ca_atoms):
        for j, atom_j in enumerate(ca_atoms):
            # Calculate the distance
            distance = atom_i - atom_j
            # Fill the matrix
            distance_matrix[i, j] = distance
    
    return distance_matrix