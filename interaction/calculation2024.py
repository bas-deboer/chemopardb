from django.conf import settings
from django.utils.text import slugify
from django.core.cache import cache

import os
import prolif as plf
from rdkit import Chem
import MDAnalysis as mda
from io import StringIO


def prolif_calculation(structure, chain_ids):
    pdbdata = StringIO(structure.pdb_data.pdb)
    
    # Load the structure using MDAnalysis
    u = mda.Universe(pdbdata, format="pdb")
    print("Now guessing bonds...")
    u.atoms.guess_bonds(vdwradii={"H": 0.4, "O": 1.48})
    
    # Select multiple chains and concatenate the AtomGroups
    ref_chain = u.select_atoms("chainID {}".format(chain_ids[0]))
    for chain_id in chain_ids[1:]:
        ref_chain += u.select_atoms("chainID {}".format(chain_id))
    
    target_chain = u.select_atoms(f"all and not resname HOH")
    
    # Generate an RDKit molecule for the reference chain
    ref_mol = plf.Molecule.from_mda(ref_chain, force=True)
    target_mol = plf.Molecule.from_mda(target_chain, force=True)
    
    # Calculate interactions using ProLIF
    fp = plf.Fingerprint()
    fp.run_from_iterable([ref_mol], target_mol)
    interactions_dict = fp.ifp
    
    return interactions_dict




def save_interactions(interactions_dict):
    # Dictionary to keep track of chemokine partner pairs to avoid duplicates
    pair_instances = {}

    for key, value in interactions_dict.items():
        for residue_pair, interactions in value.items():
            # Define identifiers for chemokine and partner
            chemokine_id = f"{residue_pair[0].name} {residue_pair[0].number}"
            partner_id = f"{residue_pair[1].name} {residue_pair[1].number}"
            pair_key = (chemokine_id, partner_id)

            # Check if this pair is already processed
            if pair_key not in pair_instances:
                # Create or get ChemokinePartnerPair instance
                pair_instance, created = ChemokinePartnerPair.objects.get_or_create(
                    structure=structure_instance,
                    partner=partner_instance,
                    defaults={'chemokine_id': chemokine_id, 'partner_id': partner_id}
                )
                pair_instances[pair_key] = pair_instance
            else:
                pair_instance = pair_instances[pair_key]

            # Create interactions
            for interaction_type, details in interactions.items():
                ChemokinePartnerInteraction.objects.create(
                    chemokine_partner_pair=pair_instance,
                    chemokine_residue=chemokine_id,
                    partner_residue=partner_id,
                    interaction_type=interaction_type,
                    # You might want to save more detailed interaction details here
                )
