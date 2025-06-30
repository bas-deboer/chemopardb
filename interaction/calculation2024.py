from django.conf import settings
from django.utils.text import slugify
from django.core.cache import cache
import warnings
from io import StringIO
import prolif as plf
from rdkit import Chem
import MDAnalysis as mda
import os


IGNORED_RESIDUES = ['HOH', 'SO4', 'AOP', 'PCA']


# proLIF - redefining HB-acceptors (slightly wider angle)
class HBAcceptor(plf.interactions.HBAcceptor):
    def __init__(self): 
        super().__init__(DHA_angle=(120, 180))


# proLIF - redefining HB-donors (slightly wider angle)
class HBDonor(plf.interactions.HBDonor):
    def __init__(self): 
        super().__init__(DHA_angle=(120, 180))
        

def prolif_calculation(chemokine_partner_complex, chemokine_chain, chemokine_residues=None, partner_chains=None, partner_residues=None):
    """
    Calculate interactions between a chemokine and its partner(s) using ProLIF.

    Parameters:
        structure: Structure object with a .pdb_data.pdb field.
        chemokine_chain: str, chain ID of the chemokine.
        chemokine_residues: list of int/str (optional), residue numbers in the chemokine chain.
        partner_chains: list of str (optional), chain IDs of the binding partner(s).
        partner_residues: list or set of int/str (optional), residue numbers to filter in partner chains.

    Returns:
        ProLIF Fingerprint object.
    """
    pdbdata = StringIO(chemokine_partner_complex.pdb_data.pdb)
    u = mda.Universe(pdbdata, format="pdb")
    u.atoms.guess_bonds(vdwradii={"H": 0.4, "O": 1.48})

    # Build chemokine selection
    if chemokine_residues:
        chem_clause = " or ".join(f"resnum {r}" for r in chemokine_residues)
        chemokine_selection = f"chainID {chemokine_chain} and ({chem_clause})"
    else:
        chemokine_selection = f"chainID {chemokine_chain}"
    print(chemokine_selection)
    # Build partner selection (based on chain and optional residues)
    if partner_chains:
        partner_clauses = []
        for chain in partner_chains:
            if partner_residues:
                res_clause = " or ".join(f"resnum {r}" for r in partner_residues)
                sel = f"chainID {chain} and ({res_clause})"
            else:
                sel = f"chainID {chain}"
            partner_clauses.append(sel)
        partner_selection = " or ".join(partner_clauses)
    
    print(partner_selection)
    
    # Select atoms
    ref_chain = u.select_atoms(chemokine_selection)
    target_chain = u.select_atoms(partner_selection)

    # Convert to ProLIF molecules
    ref_mol = plf.Molecule.from_mda(ref_chain, force=True)
    target_mol = plf.Molecule.from_mda(target_chain, force=True)

    # Calculate interactions
    fp = plf.Fingerprint()
    fp.run_from_iterable([ref_mol], target_mol)

    return fp