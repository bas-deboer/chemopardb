import csv
import gc
import json
import logging
import os
import threading
import pickle
from io import StringIO
from urllib.request import urlopen
import json
from Bio import Align

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm
from urllib.parse import urljoin
from urllib.request import urlopen

import django.apps
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import connection, IntegrityError, transaction
from django.utils.text import slugify

# Django models
from build.management.commands.base_build import Command as BaseBuild
from common.models import WebLink, WebResource, Publication, ResiduePosition
from protein.models import Protein, ProteinSegment
from structure.models import (  Structure,
                                PdbData,
                                StructureType,
                                Rotamer,
                                Entity,
                                EntityType,
                                EntityInstance,
                                ChemokineBindingPartner,
                            )
from residue.models import Residue
from interaction.models import *

# Custom APIs and calculations
from common.protossapi import *
from interaction.calculation2024 import *

# BioPython libraries
from Bio import PDB, SeqIO, pairwise2, Align
from Bio.PDB import (
    PDBParser,
    PDBIO,
    Select,
    PPBuilder,
    MMCIFParser,
    Polypeptide,
    Selection,
    Superimposer,
    parse_pdb_header,
)
from Bio.SeqUtils import seq1

# Third-party libraries
import prolif as plf
import MDAnalysis as mda

# =============================================================================
# Constants
# =============================================================================

STRUCTURE_METHOD_MAPPING = {
    "x-ray diffraction": "X-ray",
    "electron microscopy": "Cryo-EM",
    "solution-nmr": "NMR",
    "Solution nmr": "NMR",
    "solution nmr": "NMR",
    "solution nmr; solution scattering": "NMR",
    "unknown": "Unknown",
}

HARD_CODED_STRUCTURE_MAP = {
    "7XBX": {
        "protein": "P78423",        
        "chemokine_residues": list(range(0, 100))
    },
    "7F1R": {
        "protein": "P13501",
        "chemokine_residues": list(range(100)),
        "generic_mapping": {
            **{i: i + 23 for i in range(1, 15)},   # 1->24, 2->25, ... 14->37
            **{i: i + 33 for i in range(15, 32)}, # 15->48, 16->49, ... 31->64
            **{i: i + 38 for i in range(32, 51)}, # 32->70, 33->71, ... 50->88
        },
    },
    #"6EHZ": {
    #    "protein": "",
    #    "chemokine_residues": list(range(1, 68))
    #},    
}

# =============================================================================
# Helper Classes and Functions
# =============================================================================
# Custom Select class to extract a specific chain.
class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id


class ExcelChainSelect(Select):
    def __init__(self, allowed_chains):
        """
        allowed_chains: dict with chain IDs as keys.
          - If the value is None, the entire chain is kept.
          - If the value is a list of residue numbers, only those residues will be kept.
        """
        self.allowed_chains = allowed_chains

    def accept_chain(self, chain):
        return chain.id in self.allowed_chains

    def accept_residue(self, residue):
        allowed = self.allowed_chains.get(residue.parent.id, None)
        if allowed is None:
            return True
        # residue.id is typically a tuple: (hetfield, resseq, icode); we use resseq here.
        return residue.id[1] in allowed


class ChainResidueSelect(Select):
    """Selection class for saving only a specified chain and (optionally) a set of residues."""

    def __init__(self, chain_id, residues=None):
        """
        Initialize with a chain ID and an optional list of residue IDs.
        Example: chain_id="A", residues=["1", "2", "3"]
        """
        self.chain_id = chain_id
        self.residues = set(residues) if residues else None

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        if self.residues is None:
            return True
        return str(residue.id[1]) in self.residues


def fetch_entity_info(pdb_id, entity_id, entity_type):
    """
    Fetch entity info from RCSB PDB API based on type.
    For polymer entities, retrieves additional residue mapping and annotations.
    """
    url = f"https://data.rcsb.org/rest/v1/core/{entity_type}/{pdb_id}/{entity_id}"
    response = requests.get(url)

    if response.status_code != 200:
        return None

    entity_data = response.json()
    if entity_type == "polymer_entity":
        description = (
            entity_data.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
        )
        # Look for the first Pfam annotation
        pfam_accession = None
        for annotation in entity_data.get("rcsb_polymer_entity_annotation", []):
            if annotation.get("type") == "Pfam":
                pfam_accession = annotation.get("annotation_id")
                break

        organism_data = entity_data.get("rcsb_entity_source_organism", [])
        ncbi_scientific_name = (
            organism_data[0].get("ncbi_scientific_name") if organism_data else None
        )
        chain_ids = entity_data.get("rcsb_polymer_entity_container_identifiers", {}).get(
            "auth_asym_ids", []
        )
        residues = fetch_residue_mapping(chain_ids)

        return {
            "id": entity_data.get("rcsb_id"),
            "name": description,
            "type": entity_data.get("entity_poly", {}).get(
                "rcsb_entity_polymer_type", "Unknown"
            ),
            "chain": ", ".join(chain_ids),
            "residues": residues,
            "unp_accession": ", ".join(
                entity_data.get("rcsb_polymer_entity_container_identifiers", {}).get(
                    "uniprot_ids", []
                )
            ),
            "pfam_accession": pfam_accession,
            "pubchem_id": None,
            "chembl_id": None,
            "comp_id": None,
            "smiles": None,
            "inchikey": None,
            "organism": ncbi_scientific_name,
        }

    elif entity_type == "branched_entity":
        name = entity_data.get("rcsb_branched_entity_name_com", {}).get("name")
        comp_id = entity_data.get("rcsb_branched_entity_container_identifiers", {}).get(
            "chem_ref_def_id"
        )
        pubchem_id, chembl_id, smiles, inchikey = fetch_related_ids(comp_id)
        return {
            "id": entity_data.get("rcsb_id"),
            "name": name,
            "type": entity_data.get("entity_poly", {}).get(
                "rcsb_entity_polymer_type", "Unknown"
            ),
            "chain": ", ".join(
                entity_data.get("rcsb_branched_entity_container_identifiers", {}).get(
                    "auth_asym_ids", []
                )
            ),
            "residues": None,
            "unp_accession": None,
            "pfam_accession": None,
            "pubchem_id": None,
            "chembl_id": None,
            "comp_id": comp_id,
            "smiles": smiles,
            "inchikey": inchikey,
            "organism": None,
        }

    elif entity_type == "nonpolymer_entity":
        description = (
            entity_data.get("rcsb_nonpolymer_entity", {})
            .get("pdbx_description", "Unknown")
        )
        comp_id = entity_data.get("pdbx_entity_nonpoly", {}).get("comp_id")
        pubchem_id, chembl_id, smiles, inchikey = fetch_related_ids(comp_id)
        return {
            "id": entity_data.get("rcsb_id"),
            "name": description,
            "type": "Non-polymer",
            "chain": ", ".join(
                entity_data.get("rcsb_nonpolymer_entity_container_identifiers", {}).get(
                    "auth_asym_ids", []
                )
            ),
            "residues": None,
            "unp_accession": None,
            "pfam_accession": None,
            "pubchem_id": pubchem_id,
            "chembl_id": chembl_id,
            "comp_id": comp_id,
            "smiles": smiles,
            "inchikey": inchikey,
            "organism": None,
        }
    return None

def get_ifp_ccn_numbers():
    try:
        start_pos = ResiduePosition.objects.get(ccn_number="NTc.Cm50").position
        end_pos = ResiduePosition.objects.get(ccn_number="CT.50").position
        return list(
            ResiduePosition.objects.filter(
                position__gte=start_pos, position__lte=end_pos
            ).order_by("position").values_list("ccn_number", flat=True)
        )
    except ResiduePosition.DoesNotExist:
        return []


def fetch_residue_mapping(chain_ids):
    """
    Fetch the residue mapping for each chain specified in chain_ids using a GraphQL call.

    Args:
        chain_ids (list of str): List of chain IDs (e.g. ["1A15.A", "1A15.B"]).

    Returns:
        dict: Mapping from chain ID to list of residue indices.
    """
    query = """
    query getResidueMapping($instance_ids: [String!]!) {
      polymer_entity_instances(instance_ids: $instance_ids) {
        rcsb_polymer_entity_instance_container_identifiers {
          auth_to_entity_poly_seq_mapping
        }
      }
    }
    """
    variables = {"instance_ids": chain_ids}
    url = "https://data.rcsb.org/graphql"
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, headers=headers, json={"query": query, "variables": variables})

    if response.status_code != 200:
        print(f"Error: Received status code {response.status_code}")
        print("Response content:", response.content)
        return None

    try:
        data = response.json().get("data", {}).get("polymer_entity_instances", [])
        residue_mapping = {}
        for idx, entry in enumerate(data):
            residue_mapping[chain_ids[idx]] = entry["rcsb_polymer_entity_instance_container_identifiers"]["auth_to_entity_poly_seq_mapping"]
        return residue_mapping
    except (ValueError, AttributeError, IndexError) as e:
        print("Error parsing JSON response:", e)
        print("Raw response content:", response.content)
        return None


def get_pdb_entities(pdb_id):
    """
    Fetch the full PDB entity data from RCSB.

    Args:
        pdb_id (str): PDB identifier.

    Returns:
        dict or None: Parsed JSON data if successful, else None.
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    print(f"Failed to retrieve data for PDB ID: {pdb_id}")
    return None


def fetch_related_ids(comp_id):
    """
    Fetch related chemical IDs (PubChem, ChEMBL) and descriptors for a given chem comp ID.

    Args:
        comp_id (str): Chemical component identifier.

    Returns:
        tuple: (pubchem_id, chembl_id, smiles, inchikey)
    """
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
    response = requests.get(url)
    if response.status_code != 200:
        return None, None, None, None

    chemcomp_data = response.json()
    pubchem_id = None
    chembl_id = None
    smiles = None
    inchikey = None

    for entry in chemcomp_data.get("rcsb_chem_comp_related", []):
        if entry.get("resource_name") == "PubChem":
            pubchem_id = entry.get("resource_accession_code")
        elif entry.get("resource_name") == "ChEMBL":
            chembl_id = entry.get("resource_accession_code")

    descriptors = chemcomp_data.get("rcsb_chem_comp_descriptor", {})
    smiles = descriptors.get("smiles")
    inchikey = descriptors.get("in_ch_ikey")
    return pubchem_id, chembl_id, smiles, inchikey


class ParseStructureCSV:
    """
    Parse structure annotation CSV file and build a dictionary of structure data.
    """

    def __init__(self):
        self.pdb_ids = []
        self.structures = {}
        csv_path = os.path.join(
            settings.DATA_DIR, "structure_data", "annotation", "chemokine_structures.csv"
        )
        with open(csv_path, newline="") as csvfile:
            reader = csv.reader(csvfile, delimiter=";")
            next(reader, None)
            for row in reader:
                pdb_id = row[0]
                self.pdb_ids.append(pdb_id)
                self.structures[pdb_id] = {
                    "pdb": row[0],
                    "protein": row[1],
                    "name": row[0].lower(),
                    "method_from_file": row[2],
                    "resolution": row[3],
                    "state": row[4],
                    "preferred_chain": row[5],
                }

    def __str__(self):
        return f"<ParsedStructures: {len(self.pdb_ids)} entries>"


def create_binding_partner_pairs(self, structure, interactions_dict, chemokine_chain):
    """
    Group the interactions by partner chain and create binding partner pair records.
    The interactions_dict is assumed to be structured as:
    
        { interaction_type: { (chemokine_residue, partner_residue): details, ... }, ... }
        
    where each residue object has at least a 'chain' attribute.
    """
    # Group interactions by partner chain
    pairs = {}
    for interaction_type, pairs_dict in interactions_dict.items():
        for residue_pair, details in pairs_dict.items():
            chemokine_residue, partner_residue = residue_pair
            partner_chain = partner_residue.chain
            key = (chemokine_chain, partner_chain)
            pairs.setdefault(key, set()).add(interaction_type)
    
    # For each unique chemokine–partner chain pair, fetch the corresponding entities
    # and create the binding partner record.
    for (chem_chain, partner_chain), interaction_types in pairs.items():
        try:
            chemokine_entity = Entity.objects.filter(structure=structure, chain=chem_chain).first()
            partner_entity = Entity.objects.filter(structure=structure, chain=partner_chain).first()
            if not chemokine_entity or not partner_entity:
                self.logger.warning(f"Entities not found for chains {chem_chain} and {partner_chain} in structure {structure}")
                continue
            self.create_chemokine_binding_partner(
                structure=structure,
                chemokine_entity=chemokine_entity,
                partner_entity=partner_entity,
                chemokine_chain=chem_chain,
                partner_chain=partner_chain,
            )
        except Exception as e:
            self.logger.error(f"Error creating binding partner pair for chains {chem_chain} and {partner_chain}: {e}")


# =============================================================================
# Main Command
# =============================================================================




class Command(BaseBuild):
    help = "Reads source data and creates pdb structure objects"

    # Shared attributes
    lock = threading.Lock()
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    pdb_data_dir = os.path.join(settings.DATA_DIR, "structure_data", "pdbs")

    # Preload ProteinSegment objects by slug
    segments = {segment.slug: segment for segment in ProteinSegment.objects.all()}

    def add_arguments(self, parser):
        parser.add_argument(
            "-p",
            "--proc",
            type=int,
            action="store",
            dest="proc",
            default=1,
            help="Number of processes to run",
        )
        parser.add_argument(
            "-s",
            "--structure",
            dest="structure",
            nargs="+",
            help="Structure to import (PDB ID)",
        )
        parser.add_argument(
            "-u",
            "--purge",
            action="store_true",
            dest="purge",
            default=False,
            help="Purge existing records",
        )
        parser.add_argument(
            "--skip_cn",
            action="store_true",
            default=False,
            help="Skip building contact network for test build",
        )
        parser.add_argument(
            "-i",
            "--incremental",
            action="store_true",
            dest="incremental",
            default=False,
            help="Incremental update to structures for small live update",
        )
        parser.add_argument(
            "--debug",
            action="store_true",
            dest="debug",
            default=False,
            help="Print info for debugging",
        )

    def handle(self, *args, **options):
        self.stdout.write("Build command started...")

        if options["purge"]:
            try:
                self.purge_structures()
                self.purge_structuremethods()
                self.tracker = {}
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

        self.construct_errors = []
        self.rotamer_errors = []
        self.contactnetwork_errors = []
        self.interaction_errors = []
        self.missing_uniprot_ids = []

        # Step 1: Read structure list from CSV
        self.parsed_structures = ParseStructureCSV()
        if options["structure"]:
            self.parsed_structures.pdb_ids = [
                i
                for i in self.parsed_structures.pdb_ids
                if i in options["structure"] or i.lower() in options["structure"]
            ]

        self.incremental_mode = options["incremental"]
        self.debug = options["debug"]

        try:
            self.logger.info("CREATING STRUCTURES")
            # Step 2: Run main processing function
            self.main_func()
            self.logger.info("COMPLETED CREATING STRUCTURES")
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        #print("Construct errors:")
        #print(self.construct_errors)
        #print("Rotamer errors:")
        #print(self.rotamer_errors)
        #print("Interaction errors:")
        #print(self.interaction_errors)
        if self.missing_uniprot_ids:
            print("Missing UniProt IDs:")
            print(self.missing_uniprot_ids)
        else:
            print("No missing UniProt IDs.")

    # -------------------------------------------------------------------------
    # Helper Methods
    # -------------------------------------------------------------------------

    def purge_structures(self):
        Structure.objects.all().delete()
        PdbData.objects.all().delete()

    def purge_structuremethods(self):
        StructureType.objects.all().delete()

    def get_pdb_sequence(self, pdb_data, chain_id):
        pdb_io = StringIO(pdb_data)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("PDB", pdb_io)

        sequence = ""
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    ppb = PPBuilder()
                    for pp in ppb.build_peptides(chain):
                        sequence += str(pp.get_sequence())
        return sequence

    def align_sequences(self, seq1, seq2):
        """
        Align two sequences using a global pairwise alignment with custom gap penalties.
        """
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        alignment = aligner.align(seq1, seq2)
        return alignment[0]

    def create_mapping(self, alignment):
        """
        Create a mapping from PDB sequence positions to database sequence positions.
        """
        target_seq, query_seq = alignment[0], alignment[1]
        target_pos = 1
        query_pos = 1
        mapping = {}
        for t, q in zip(target_seq, query_seq):
            if t != "-":
                if q != "-":
                    mapping[target_pos] = query_pos
                target_pos += 1
            if q != "-":
                query_pos += 1
        return mapping
    
    def get_segment_from_ccn_number(self, ccn_number):
       """
       Given a CCN (chemokine common numbering) string (e.g. 'NTc.Cm9', 'B1.Cm14', ...),
       return the segment name.
       """
       if ccn_number:
           if ccn_number.startswith('NTc.Cm'):
               return 'N-term'
           if ccn_number.startswith('CX.'):
               return 'CX'
           if ccn_number.startswith('cxb1.'):
               return 'N-loop'
           if ccn_number.startswith('B1.'):
               return 'B1'
           if ccn_number.startswith('b1b2.'):
               return '30s-loop'
           if ccn_number.startswith('B2.'):
               return 'B2'
           if ccn_number.startswith('b2b3.'):
               return '40s-loop'
           if ccn_number.startswith('B3.'):
               return 'B3'
           if ccn_number.startswith('b3h.'):
               return '50s-loop'
           if ccn_number.startswith('H.'):
               return 'Helix'
           if ccn_number.startswith('CT.'):
               return 'C-term'
       return None

    
    

    def create_rotamers(self, structure, sd):
        # If the chain_id field contains multiple chains separated by commas,
        # split them into a list and strip whitespace.
        if isinstance(structure.chain_id, str):
            preferred_chains = [c.strip() for c in structure.chain_id.split(",") if c.strip()]
        else:
            preferred_chains = structure.chain_id

        if not structure.protein or not structure.protein.sequence:
            self.logger.error(f"Protein data or sequence not available for structure {structure}")
            return

        pdb_code = sd["pdb"].upper()
        allowed_residues = None
        manual_mapping = None

        if pdb_code in HARD_CODED_STRUCTURE_MAP:
            hard_coded_info = HARD_CODED_STRUCTURE_MAP[pdb_code]
            allowed_residues = set(str(r) for r in hard_coded_info.get("chemokine_residues", []))
            manual_mapping = hard_coded_info.get("generic_mapping")

        # Iterate over each valid chain.
        for chain_id in preferred_chains:
            pdb_sequence = self.get_pdb_sequence(structure.pdb_data.pdb, chain_id)
            db_sequence = structure.protein.sequence

            if manual_mapping:
                mapping = manual_mapping
                self.logger.info(f"Using manual generic mapping for structure {pdb_code} on chain {chain_id}")
            else:
                alignment = self.align_sequences(pdb_sequence, db_sequence)
                mapping = self.create_mapping(alignment)

            protein = structure.protein
            pdb_io = StringIO(structure.pdb_data.pdb)
            parser = PDBParser(QUIET=True)
            pdb_structure = parser.get_structure("structure", pdb_io)
            model = pdb_structure[0]
            try:
                chain = model[chain_id]
            except KeyError:
                self.logger.error(f"Chain {chain_id} not found in structure {structure.pdb_code.index}")
                continue

            residues = list(chain.get_residues())
            with transaction.atomic():
                for seq_pos, residue in enumerate(residues, start=1):
                    if residue.id[0] != " ":
                        continue
                    if seq_pos not in mapping:
                        self.logger.warning(f"Sequence position {seq_pos} not found in mapping.")
                        continue
                    db_res_id = mapping[seq_pos]
                    if allowed_residues and str(db_res_id) not in allowed_residues:
                        continue
                    try:
                        db_residue = Residue.objects.get(protein=protein, sequence_number=db_res_id)
                        residue_type = (
                            "canonical"
                            if db_residue.amino_acid == seq1(residue.resname)
                            else "mutated"
                        )
                        Rotamer.objects.create(
                            residue=db_residue,
                            structure=structure,
                            pdbdata=structure.pdb_data,
                            pdbseq_number=residue.id[1],
                            sequence_number=db_res_id,
                            #generic_number=db_residue.generic_number,
                            ccn_number=db_residue.ccn_number,
                            amino_acid=seq1(residue.resname),
                            amino_acid_three_letter=residue.resname,
                            residue_type=residue_type,
                            segment=db_residue.segment,
                            chain=chain_id,
                        )
                    except Residue.DoesNotExist:
                        self.logger.error(
                            f"Database residue not found for sequence number {db_res_id} in protein {structure.protein}"
                        )

    def parse_residue_range(self, rangestr):
        residues = []
        if rangestr is None or pd.isna(rangestr):
            return residues
        if not isinstance(rangestr, str):
            try:
                val = int(float(rangestr))
                residues.append(val)
                return residues
            except Exception:
                return residues
        for part in rangestr.split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                start, end = map(float, part.split("-"))
                residues.extend(range(int(start), int(end) + 1))
            else:
                try:
                    residues.append(int(float(part)))
                except Exception:
                    continue
        return residues




    def superimpose_structures(self, template_structure, target_structure, chains):
        """
        Superimpose the target structure onto the template structure for specified chains.
        """
        parser = PDBParser(QUIET=True)
        super_imposer = Superimposer()
        template = parser.get_structure("template", StringIO(template_structure))
        target = parser.get_structure("target", StringIO(target_structure.pdb_data.pdb))

        for chain_id in chains:
            try:
                template_chain = template[0][chain_id]
                target_chain = target[0][chain_id]
                template_atoms = list(template_chain.get_atoms())
                target_atoms = list(target_chain.get_atoms())
                super_imposer.set_atoms(template_atoms, target_atoms)
                super_imposer.apply(target.get_atoms())
                self.logger.info(
                    f"Successfully superimposed chain {chain_id} of target structure {target_structure.pdb_code.index} onto template"
                )
            except KeyError:
                self.logger.error(f"Chain {chain_id} not found in one of the structures.")
        return target

    def save_structure_with_contacts(self, structure, chains, file_path):
        """
        Save a PDB file containing only the specified chains and any chain that makes contact.
        """

        class ChainSelect(Select):
            def accept_chain(self, chain):
                if chain.id in chains:
                    return True
                # Return chains that contain any residue with a CA atom
                return any(residue.has_id("CA") for residue in chain.get_residues())

        parser = PDBParser(QUIET=True)
        structure_obj = parser.get_structure(structure.pdb_code.index, StringIO(structure.pdb_data.pdb))
        io = PDBIO()
        io.set_structure(structure_obj)
        io.save(file_path, ChainSelect())

    def calculate_ligand_interactions(self, protein_file, ligand_file, output_html, structure_id):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(protein_file)
            u.atoms.guess_bonds(vdwradii={"H": 0.4, "O": 1.48})
            protein_mol = plf.Molecule.from_mda(u)
            pose_iterable = plf.sdf_supplier(ligand_file)
            pose_list = list(pose_iterable)
            if not pose_list:
                raise ValueError("No ligands found in the ligand file.")
            fp = plf.Fingerprint()
            fp.run_from_iterable(pose_list, protein_mol)
            structure_instance = Structure.objects.get(id=structure_id)

            for i, lig_mol in enumerate(pose_list):
                try:
                    network = fp.plot_lignetwork(lig_mol)
                    html_content = (
                        network._repr_html_() if hasattr(network, "_repr_html_")
                        else network if isinstance(network, str)
                        else None
                    )
                    if html_content is None:
                        raise TypeError("The returned object from LigNetwork does not have HTML content")
                    with open(f"{output_html}_{i}.html", "w") as f:
                        f.write(html_content)
                    LigandNetworkHTML.objects.create(
                        ligand_name=f"Ligand {i}",
                        html_content=html_content,
                        structure=structure_instance,
                    )
                except Exception as e:
                    print(f"Error processing ligand at index {i}: {e}")
                    continue

    def create_binding_partner_pairs(self, structure, interactions_dict, chemokine_chain):
        """
        Group interactions by partner chain and create binding partner pair records.
        Assumes interactions_dict is structured as:
        
            { interaction_type: { (chemokine_residue, partner_residue): details, ... }, ... }
        
        where each residue object has at least a 'chain' attribute.
        """
        # Group interactions by partner chain
        pairs = {}
        for interaction_type, pairs_dict in interactions_dict.items():
            for residue_pair, details in pairs_dict.items():
                chemokine_residue, partner_residue = residue_pair
                partner_chain = partner_residue.chain  # Use the chain value from the partner residue
                key = (chemokine_chain, partner_chain)
                pairs.setdefault(key, set()).add(interaction_type)
        
        # For each unique chemokine–partner chain pair, look up the corresponding entities
        # Using an icontains lookup in case the stored chain string is comma–separated.
        for (chem_chain, partner_chain), interaction_types in pairs.items():
            try:
                chemokine_entity = Entity.objects.filter(structure=structure, chain__icontains=chem_chain).first()
                partner_entity = Entity.objects.filter(structure=structure, chain__icontains=partner_chain).first()
                if not chemokine_entity or not partner_entity:
                    self.logger.warning(
                        f"Entities not found for chains {chem_chain} and {partner_chain} in structure {structure}"
                    )
                    continue
                self.create_chemokine_binding_partner(
                    structure=structure,
                    chemokine_entity=chemokine_entity,
                    partner_entity=partner_entity,
                    chemokine_chain=chem_chain,
                    partner_chain=partner_chain,
                )
            except Exception as e:
                self.logger.error(f"Error creating binding partner pair for chains {chem_chain} and {partner_chain}: {e}")


    def fetch_and_create_protein_with_ccn(
        self,
        accession,
        allpara_sequences,
        alignment_ccn_map,
        alignment_ccn_columns,
        logger,
        create_residues_fn,
    ):
        """
        Creates a Protein entry for a missing UniProt accession,
        aligns to closest homolog in ALL_para_df_mod.csv, and creates residues with mapped CCN numbers.
        
        Args:
            accession (str): UniProt accession to fetch and create.
            allpara_sequences (dict): {uniprot: sequence, ...} from ALL_para_df_mod.csv.
            alignment_ccn_map (dict): {uniprot: [CCN_label or None, ...]}.
            alignment_ccn_columns (list): Alignment column headers for CCN numbers.
            logger (logging.Logger): Logger for info/warning/error.
            create_residues_fn (function): Function to create residues (protein, sequence, ccn_map).

        Returns:
            Protein object, or None on error.
        """
        # 1. Fetch UniProt Data
        url = f"https://rest.uniprot.org/uniprotkb/{accession}"
        try:
            with urlopen(url) as uf:
                data = json.load(uf)
            sequence = data.get("sequence", {}).get("value", "")
            if not sequence:
                logger.error(f"UniProt sequence not found for {accession}")
                return None
            full_name = (
                data.get("proteinDescription", {})
                .get("recommendedName", {})
                .get("fullName", {})
                .get("value", accession)
            )
            species_latin = data.get("organism", {}).get("scientificName", "Unknown species")
        except Exception as e:
            logger.error(f"Failed to fetch UniProt entry {accession}: {e}")
            return None

        # 2. Find closest homolog in ALL_para_df_mod
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        best_score = -float("inf")
        best_acc = None
        for acc, ref_seq in allpara_sequences.items():
            score = aligner.score(sequence, ref_seq)
            if score > best_score:
                best_score = score
                best_acc = acc
        if not best_acc:
            logger.error(f"No homolog found for {accession} in ALL_para_df_mod.csv")
            return None
        logger.info(f"Closest homolog for {accession} is {best_acc} (score={best_score:.1f})")

        ref_seq = allpara_sequences[best_acc]
        ref_ccn_list = alignment_ccn_map[best_acc]
        # ref_ccn_columns is available if needed

        # 3. Map CCN numbers using alignment
        alignment = aligner.align(sequence, ref_seq)[0]
        aligned_seq, aligned_ref = str(alignment[0]), str(alignment[1])
        pos_seq = 1
        pos_ref = 1
        ccn_map = {}
        for i in range(len(aligned_seq)):
            aa_seq = aligned_seq[i]
            aa_ref = aligned_ref[i]
            if aa_seq != "-" and aa_ref != "-":
                label = ref_ccn_list[pos_ref - 1] if pos_ref - 1 < len(ref_ccn_list) else None
                if label:
                    ccn_map[pos_seq] = label
            if aa_seq != "-":
                pos_seq += 1
            if aa_ref != "-":
                pos_ref += 1

        # 4. Create Protein and Species objects
        from protein.models import Protein, Species  # or your project import
        protein_name = f"{full_name} ({accession})"
        species, _ = Species.objects.get_or_create(latin_name=species_latin)
        protein = Protein.objects.create(
            name=protein_name,
            gene_name=accession,
            uniprot_id=accession,
            sequence=sequence,
            species=species,
            # Optionally set default family, type, etc.
        )
        logger.info(f"Auto-created Protein {accession} ({protein_name}) with {len(sequence)} residues.")

        # 5. Create Residues with CCN numbers
        create_residues_fn(protein, sequence, ccn_map)

        # Optionally mark as "auto-created" for later curation (add a flag or annotation)
        return protein

    import os











    # Main Function
    # -------------------------------------------------------------------------
    def main_func(self):
        """
        Main processing function that iterates over parsed structures, downloads and processes PDB files,
        creates Structure and Rotamer objects, fetches entity annotations, and calculates interactions.
        Additionally, for each chemokine-partner combination from the Excel file, a new PDB file is created
        containing only the specified chemokine and partner chains (with optional residue filtering).
        Finally, the chemokine chain is extracted, superimposed onto a template via rotamer generic numbers,
        the transformation is applied to the chemokine–partner complex, and the final aligned complex is saved
        in the database.
        """
        print(">>> main_func called")
        
        # Ensure required entity types exist in the database.
        for entity_name in ["Protein", "DNA", "RNA", "Ligand", "Non-polymer"]:
            EntityType.objects.get_or_create(name=entity_name)

        # Set up necessary directories.
        interactions_data_dir = os.path.join(settings.DATA_DIR, "structure_data", "interactions")
        os.makedirs(interactions_data_dir, exist_ok=True)
        protonated_dir = os.path.join(settings.DATA_DIR, "structure_data", "output_protonated_pdbs")
        complexes_dir = os.path.join(settings.DATA_DIR, "structure_data", "chemokine_partner_complexes")
        os.makedirs(complexes_dir, exist_ok=True)
        aligned_dir = os.path.join(settings.DATA_DIR, "structure_data", "aligned_complexes")
        os.makedirs(aligned_dir, exist_ok=True)

        # Load the Excel file (ChemoPar_masterlist) with binding partner information.
        excel_file_path = os.path.join(settings.DATA_DIR, "ChemoPar_masterlist_2025.xlsx")
        df_excel = pd.read_excel(excel_file_path, sheet_name="Partners")

        # Helper functions used for filtering PDB lines.
        def is_valid_line(line, chain):
            """
            Check if a PDB line is valid for a given chain and excludes water or unwanted molecules.
            """
            return (line.startswith("ATOM") or line.startswith("HETATM")) and (
                len(line) > 21 and line[21] == chain and
                line[17:20].strip() not in {"HOH", "WAT", "DOD", "TIP", "SOL", "OH2", "D20", "SO4"}
            )


        # Needed for superpositioning
        parser = PDBParser(QUIET=True)
        core_ccn_numbers = ["CX.5",                                                     # Conserved cysteine
                            "B1.1", "B1.2", "B1.3", "B1.4", "B1.5", "B1.6", "B1.7",     # B1 sheet
                            "B2.1", "B2.2", "B2.3", "B2.4", "B2.5", "B2.6", "B3.1",     # B2 sheet
                            "B3.2", "B3.3", "B3.4"                                      # B3 sheet
                            ]

        template_seq_ccn_dict = {
            "4XDX": {
                "CX.5": 9,
                "B1.1": 21, "B1.2": 22, "B1.3": 23, "B1.4": 24, "B1.5": 25, "B1.6": 26, "B1.7": 27,
                "B2.1": 38, "B2.2": 39, "B2.3": 40, "B2.4": 41, "B2.5": 42, "B2.6": 43,
                "B3.1": 48, "B3.2": 49, "B3.3": 50, "B3.4": 51
            }
        }

        template_ref_path = os.path.join(settings.DATA_DIR, "structure_data", "template_4xdx.pdb")
        if not os.path.exists(template_ref_path):
            self.stdout.write(self.style.ERROR(f"Template reference file not found at {template_ref_path}."))
            return        
        template_chain_id = "A"
        template_ccn_map = template_seq_ccn_dict["4XDX"]

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("template", template_ref_path)
        model = structure[0]
        chain = model[template_chain_id]

        template_backbone_atoms = {}

        for ccn_number, pdb_resnum in template_ccn_map.items():
            residue = None
            for res in chain:
                if res.id[1] == pdb_resnum:
                    residue = res
                    break
            if residue and all(atom in residue for atom in ["N", "CA", "C"]):
                template_backbone_atoms[ccn_number] = [residue["N"], residue["CA"], residue["C"]]
            else:
                print(f"Warning: Could not find all backbone atoms for {ccn_number} at residue {pdb_resnum}")

        print(f"Loaded {len(template_backbone_atoms)} reference backbone CCN positions from 4XDX chain {template_chain_id}")



        # Create a PDB parser and load the template structure.
        parser = PDBParser(QUIET=True)
        template_ref_structure = parser.get_structure("template_ref", template_ref_path)
        template_ref_model = template_ref_structure[0]

        def get_backbone_atoms(model, rotamers, offset=0):
            """
            Given a Bio.PDB structure (model[0]) and a list of rotamers
            (with keys 'chain', 'pdbseq_number', 'ccn_number'),
            return a dictionary mapping each ccn_number to a list of backbone atoms (N, CA, C).
            """
            backbone_atoms = {}
            for rotamer in rotamers:
                chain_id = rotamer['chain']
                resnum = rotamer['pdbseq_number'] + offset
                ccn_number = rotamer['ccn_number']
                if not ccn_number:
                    continue  # skip rotamers without a CCN number
                try:
                    chain = model[0][chain_id]
                except KeyError:
                    continue
                found_res = None
                for residue in chain:
                    if residue.id[1] == resnum:
                        found_res = residue
                        break
                if not found_res:
                    continue
                try:
                    backbone_atoms[ccn_number] = [
                        found_res['N'],
                        found_res['CA'],
                        found_res['C']
                    ]
                except KeyError:
                    continue
            return backbone_atoms


        # Loop over each parsed structure.
        for pdb_id in self.parsed_structures.pdb_ids:
            sd = self.parsed_structures.structures[pdb_id]
            pdb_code = sd["pdb"].upper()

            # SKIP IF STRUCTURE ALREADY EXISTS
            if Structure.objects.filter(pdb_code__index=pdb_code).exists():
                self.logger.info(f"Skipping {pdb_code}, structure already exists in the database.")
                continue

            # Determine protein ID and any hard-coded chemokine residue overrides.
            if pdb_code in HARD_CODED_STRUCTURE_MAP:
                hard_coded_info = HARD_CODED_STRUCTURE_MAP[pdb_code]
                protein_id = hard_coded_info["protein"]
                chemokine_residues = [str(r) for r in hard_coded_info["chemokine_residues"]]
            else:
                protein_id = sd["protein"]
                chemokine_residues = None

            # Fetch the Protein object from the database.
            try:
                con = Protein.objects.get(uniprot_id=protein_id)
            except Protein.DoesNotExist:
                self.logger.error(f"Protein {protein_id} does not exist for structure {pdb_code}. Skipping!")
                # record the missing UniProt ID
                self.missing_uniprot_ids.append(protein_id)
                continue

            # Fetch or download the PDB file.
            pdb_path = os.path.join(self.pdb_data_dir, f"{sd['pdb']}.pdb")
            if not os.path.isfile(pdb_path):
                self.logger.info(f"Downloading PDB file for {sd['pdb']}")
                try:
                    url = f"http://www.rcsb.org/pdb/files/{sd['pdb']}.pdb"
                    pdbdata_raw = urlopen(url).read().decode("utf-8")
                    with open(pdb_path, "w") as f:
                        f.write(pdbdata_raw)
                except Exception as e:
                    self.logger.error(f"Failed to download PDB file for {sd['pdb']}: {e}")
                    continue
            else:
                with open(pdb_path, "r") as pdb_file:
                    pdbdata_raw = pdb_file.read()

            # Parse annotations from the raw PDB file.
            for line in pdbdata_raw.splitlines():
                if line.startswith("EXPDTA"):
                    sd["structure_method"] = line[10:].strip()
                if line.startswith("REVDAT   1"):
                    sd["publication_date"] = line[13:22]
                if line.startswith("JRNL        PMID"):
                    sd["pubmed_id"] = line[19:].strip()
                if line.startswith("JRNL        DOI"):
                    sd["doi_id"] = line[19:].strip()

            header_dict = parse_pdb_header(StringIO(pdbdata_raw))
            sd["publication_date"] = header_dict.get("release_date", None)
            sd["resolution"] = header_dict.get("resolution", None)

            # Map the full structure method to a shortened name and ensure a StructureType exists.
            structure_method_full = sd.get("structure_method", "unknown").lower()
            structure_method_short = STRUCTURE_METHOD_MAPPING.get(structure_method_full, structure_method_full).capitalize()
            structure_type_slug = slugify(structure_method_short)
            try:
                structure_type, created = StructureType.objects.get_or_create(
                    slug=structure_type_slug, defaults={"name": structure_method_short}
                )
                if created:
                    print(f"Created new structure type: {structure_method_short}")
            except IntegrityError:
                structure_type = StructureType.objects.get(slug=structure_type_slug)

            # Get protonated PDB file.
            prot_pdb__file = os.path.join(protonated_dir, f"{pdb_code}_protein.pdb")
            if not os.path.isfile(prot_pdb__file):
                self.logger.warning(f"Protoss PDB not found: {prot_pdb__file}")
                continue

            with open(prot_pdb__file, "r") as f:
                prot_pdb_file_content = f.read()

            # Save the raw PDB data to the database.
            pdbdata_obj, _ = PdbData.objects.get_or_create(pdb=prot_pdb_file_content)

            # Save the base Structure object.
            try:
                s = Structure()
                s.protein = con
                s.pdb_code, _ = WebLink.objects.get_or_create(
                    index=pdb_code,
                    web_resource=WebResource.objects.get(slug="pdb"),
                )
                s.chain_id = sd["preferred_chain"]
                s.structure_type = structure_type
                s.resolution = float(sd["resolution"]) if sd.get("resolution") else None
                s.state = sd["state"]
                s.publication_date = sd.get("publication_date")
                s.pdb_data = pdbdata_obj
                try:
                    if "doi_id" in sd:
                        s.publication = Publication.get_or_create_from_doi(sd["doi_id"])
                    elif "pubmed_id" in sd:
                        s.publication = Publication.get_or_create_from_pubmed(sd["pubmed_id"])
                except Exception as e:
                    self.logger.error(f"Error saving publication for {pdb_code}: {e}")
                s.save()
                print("Saved Structure with pdb_code__index =", s.pdb_code.index)
                self.logger.info(f"Saved base Structure object for PDB {pdb_code}")
            except Exception as e:
                self.logger.error(f"Error creating base Structure object for PDB {pdb_code}: {e}")
                continue

            # Create rotamer objects by aligning the PDB sequence with the protein sequence.
            try:
                self.create_rotamers(structure=s, sd=sd)
            except Exception as e:
                self.logger.error(f"Error creating rotamers for PDB {pdb_code}: {e}")

            # Fetch additional entity annotations from the RCSB PDB API.
            response = get_pdb_entities(s.pdb_code.index)
            if not response:
                self.logger.warning(f"Failed to fetch entity annotations for {s.pdb_code.index}")
                continue

            pdb_id_resp = response["entry"]["id"]
            entity_ids = response["rcsb_entry_container_identifiers"]["entity_ids"]
            polymer_entity_ids = response.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
            nonpolymer_entity_ids = response.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids", [])
            branched_entity_ids = response.get("rcsb_entry_container_identifiers", {}).get("branched_entity_ids", [])

            for entity_id in entity_ids:
                if entity_id in polymer_entity_ids:
                    entity_type = "polymer_entity"
                elif entity_id in nonpolymer_entity_ids:
                    entity_type = "nonpolymer_entity"
                elif entity_id in branched_entity_ids:
                    entity_type = "branched_entity"
                else:
                    continue

                entity_info = fetch_entity_info(pdb_id_resp, entity_id, entity_type)
                if entity_info:
                    entity_type_instance, _ = EntityType.objects.get_or_create(name=entity_info["type"])
                    Entity.objects.create(
                        structure=s,
                        entity_type=entity_type_instance,
                        name=entity_info["name"],
                        chain=entity_info.get("chain", ""),
                        unp_accession=entity_info.get("unp_accession", ""),
                        pfam_accession=entity_info.get("pfam_accession", ""),
                        pubchem_id=entity_info.get("pubchem_id", None),
                        comp_id=entity_info.get("comp_id", None),
                        chembl_id=entity_info.get("chembl_id", None),
                        smiles=entity_info.get("smiles", None),
                        inchikey=entity_info.get("inchikey", None),
                        organism=entity_info.get("organism", None),
                        residues=entity_info.get("residues", None),
                    )
                    
            # Filter Excel rows for the current PDB code.
            df_pdb = df_excel[df_excel["PDB"].str.upper() == pdb_code]
            if df_pdb.empty:
                self.logger.info(f"No binding partner data found for PDB {pdb_code}. Skipping Excel-driven filtering.")
                continue
            
            # Process each chemokine–partner combination from the Excel file.
            for idx, row in df_pdb.iterrows():
                chem_chain = str(row["ChemokineChain"]).strip()
                chem_name = str(row["ChemokineName"]).strip()
                partner_chain_val = str(row["PartnerChain"]).strip() if not pd.isna(row["PartnerChain"]) else None
                partner_name_val = str(row["PartnerName"]).strip() if not pd.isna(row["PartnerName"]) else None
                partner_type_val = str(row["PartnerType"]).strip() if not pd.isna(row["PartnerType"]) else None
                combination_id = f"{pdb_code}_{chem_chain}"
                partner_residues_raw = row.get("PartnerResidues", None)

                # Parse partner_residues from Excel to a list of ints
                partner_residues_raw = row.get("PartnerResidues", None)
                partner_residues_clean = self.parse_residue_range(partner_residues_raw)
                residue_ids = self.parse_residue_range(partner_residues_raw)


                if partner_chain_val:
                    combination_id += f"_{partner_chain_val}"

                if partner_chain_val is None:
                    continue

                # -----------------------------
                # Step A: Create the chemokine–partner complex file first.
                # -----------------------------
                pdb_protein_file = os.path.join(protonated_dir, f"{pdb_code}_protein.pdb")
                if not os.path.isfile(pdb_protein_file):
                    self.logger.warning(f"Protoss PDB not found: {pdb_protein_file}")
                    continue
                
                with open(pdb_protein_file, "r") as f:
                    pdbdata_prot_lines = f.readlines()

                # Get chemokine lines from the protonated PDB.
                chem_lines_row = [l for l in pdbdata_prot_lines if is_valid_line(l, chem_chain)]

                partner_chain_raw = row.get("PartnerChain", "")
                partner_chains = [c.strip() for c in str(partner_chain_raw).split(",") if c.strip()]


                partner_lines = []
                for pchain in partner_chains:
                    for line in pdbdata_prot_lines:
                        if is_valid_line(line, pchain):
                            try:
                                res_id = int(line[22:26].strip())
                                if not residue_ids or res_id in residue_ids:
                                    partner_lines.append(line)
                            except ValueError:
                                if not residue_ids:
                                    partner_lines.append(line)
                if not partner_lines:
                    for pchain in partner_chains:
                        partner_lines.extend([l for l in pdbdata_prot_lines if is_valid_line(l, pchain)])

                # Combine chemokine and partner lines.
                combined_lines = chem_lines_row + partner_lines
                output_filename = f"{pdb_code}_{chem_chain}_{'_'.join(partner_chains)}_complex.pdb"
                partner_complex_path = os.path.join(complexes_dir, output_filename)
                try:
                    with open(partner_complex_path, "w") as f:
                        f.write("".join(combined_lines))
                    self.stdout.write(self.style.SUCCESS(f"Complex saved: {partner_complex_path}"))
                except Exception as e:
                    self.stdout.write(self.style.ERROR(f"Could not save {output_filename}: {e}"))
                    continue
                
                # -----------------------------
                # Step B: Superposition of the chemokine chain.
                # -----------------------------
                if not chem_chain:
                    self.stdout.write(self.style.WARNING(f"No chemokine chain found for {combination_id}. Skipping."))
                    continue
                
                # Extract the chemokine chain from the protonated PDB.
                chem_output_dir = os.path.join(settings.DATA_DIR, "structure_data", "temp")
                os.makedirs(chem_output_dir, exist_ok=True)
                chem_output_path = os.path.join(chem_output_dir, f"{combination_id}_chem.pdb")
                io = PDBIO()
                full_structure = parser.get_structure("full", prot_pdb__file)
                io.set_structure(full_structure)
                io.save(chem_output_path, ChainSelect(chem_chain))
                self.stdout.write(self.style.SUCCESS(f"Extracted chemokine chain {chem_chain} to {chem_output_path}."))

                # Load the extracted chemokine structure and its rotamers.
                chem_structure = parser.get_structure("chem", chem_output_path)
                chem_model = chem_structure[0]
                try:
                    chem_chain_struct = chem_model[chem_chain]
                except KeyError:
                    self.stdout.write(self.style.WARNING(f"Chain {chem_chain} not found in extracted file for {combination_id}."))
                    continue
                
                print("Trying to fetch Structure with pdb_code:", pdb_code)
                print("All Structure pdb_codes in DB:", list(Structure.objects.values_list('pdb_code__index', flat=True)))
                try:
                    sample_structure = Structure.objects.get(pdb_code__index=pdb_code)
                except Structure.DoesNotExist:
                    self.stdout.write(self.style.WARNING(f"Structure {pdb_code} not found in database."))
                    continue
                
                chem_rotamers = Rotamer.objects.filter(
                    structure=sample_structure,
                    ccn_number__in=core_ccn_numbers,
                    chain=chem_chain
                ).values('pdbseq_number', 'ccn_number', 'chain')

                # Get backbone atoms from the template and the extracted chem structure.
                template_atoms = template_backbone_atoms
                sample_atoms = get_backbone_atoms(chem_structure, chem_rotamers)

                common_ccn_numbers = set(template_atoms.keys()) & set(sample_atoms.keys())
                if not common_ccn_numbers:
                    self.stdout.write(self.style.WARNING(f"No common generic numbers found for {combination_id}."))
                    continue
                
                fixed_atoms = [atom for gn in common_ccn_numbers for atom in template_atoms[gn]]
                moving_atoms = [atom for gn in common_ccn_numbers for atom in sample_atoms[gn]]

                super_imposer = Superimposer()
                super_imposer.set_atoms(fixed_atoms, moving_atoms)
                super_imposer.apply(chem_model.get_atoms())
                rmsd = super_imposer.rms
                rotation_matrix, translation_vector = super_imposer.rotran
                self.stdout.write(self.style.SUCCESS(f"Superposition RMSD for {combination_id}: {rmsd:.3f} Å"))

                # -----------------------------
                # Step C: Apply the transformation to the chemokine-partner complex.
                # -----------------------------
                if not os.path.exists(partner_complex_path):
                    self.stdout.write(self.style.WARNING(f"Partner complex file missing for {combination_id}."))
                    continue
                complex_structure = parser.get_structure("complex", partner_complex_path)
                super_imposer.apply(complex_structure.get_atoms())

                final_complex_path = os.path.join(aligned_dir, f"{combination_id}_final.pdb")
                io.set_structure(complex_structure)
                io.save(final_complex_path)
                self.stdout.write(self.style.SUCCESS(f"Final aligned complex for {combination_id} saved to {final_complex_path}."))

                # -----------------------------
                # Step D: Save transformation data and update the database.
                # -----------------------------
                EntityInstance.objects.create(
                    structure=sample_structure,
                    chain_id=chem_chain,
                    rmsd=rmsd,
                    rotation_matrix=rotation_matrix.tolist(),
                    translation_vector=translation_vector.tolist()
                )
                with open(final_complex_path, "r") as f:
                    final_pdb_content = f.read()
                final_pdb_obj, created = PdbData.objects.get_or_create(pdb=final_pdb_content)
                self.stdout.write(self.style.SUCCESS(f"Database updated for combination {combination_id}."))


                # Create or link Partner object if PartnerType != "Chemokine"
                if partner_name_val and partner_type_val.lower() != "chemokine":
                    from partner.models import Partner  # Import at the top if needed

                    # Create or get the Partner object
                    partner_obj, created = Partner.objects.get_or_create(
                        name=partner_name_val,
                        defaults={"type": partner_type_val}
                    )

                    # Link structure to the partner if not already linked
                    if hasattr(partner_obj, "structures"):  # if many-to-many
                        partner_obj.structures.add(s)
                    elif hasattr(s, "partner"):  # if structure has FK to Partner
                        s.partner = partner_obj
                        s.save()

                # Continue with creating and saving the chemokine-partner binding record...
                try:
                    chemokine_partner = ChemokineBindingPartner.objects.create(
                        structure=s,
                        chemokine_chain=chem_chain,
                        chemokine_name=chem_name,
                        partner_chain=partner_chain_val,
                        partner_name=partner_name_val,
                        partner_type=partner_type_val,
                        partner_residues=partner_residues_clean,
                        pdb_data=final_pdb_obj,
                        partner=partner_obj,
                    )
                    print(f"Created chemokine-partner pair record: {combination_id}")
                except Exception as e:
                    print(f"Error saving chemokine-partner pair for {combination_id}: {e}")
                    continue
                
                # -----------------------------
                # Step E: Calculate interactions using ProLIF.
                # -----------------------------
                interactions_file_path = os.path.join(interactions_data_dir, f"{pdb_code}_{combination_id}_interactions.pkl")
                try:
                    if os.path.exists(interactions_file_path):
                        self.logger.info(f"Interactions file {interactions_file_path} already exists. Loading interactions data.")
                        with open(interactions_file_path, 'rb') as f:
                            fp = pickle.load(f)
                    else:
                        fp = prolif_calculation(
                            chemokine_partner_complex=chemokine_partner,
                            chemokine_chain=chem_chain,
                            chemokine_residues=chemokine_residues,
                            partner_chains=partner_chains,
                            partner_residues=partner_residues_clean,
                        )
                    
                    fp.to_pickle(interactions_file_path)

                    # Iterate through the interactions and create ChemokinePartnerInteraction records.
                    interactions_dict = fp.ifp
                    for key, value in interactions_dict.items():
                        for residue_pair, interactions in value.items():
                            chemokine_residue_id = f"{residue_pair[0].name} {residue_pair[0].number}"
                            partner_residue_id = f"{residue_pair[1].name} {residue_pair[1].number}"
                            try:
                                rotamer_object = Rotamer.objects.get(
                                    structure=s,
                                    pdbseq_number=residue_pair[0].number,
                                    chain=chem_chain,
                                )
                                for interaction_type, details in interactions.items():
                                    ChemokinePartnerInteraction.objects.create(
                                        chemokine_residue=rotamer_object,
                                        partner_residue=partner_residue_id,
                                        partner_chain=residue_pair[1].chain,
                                        interaction_type=interaction_type,
                                        structure=s,
                                        chemokine_binding_partner=chemokine_partner,  # Link here
                                    )
                            except Rotamer.DoesNotExist:
                                self.logger.warning(
                                    f"Rotamer not found for residue {residue_pair[0]} in {s.pdb_code.index}"
                                )
                                continue
                                             
                    print(f"Interaction data processed for {s.pdb_code.index} combination {combination_id}")
                        
                except Exception as e:
                    self.logger.error(f"Error during interaction calculation for {s.pdb_code.index} combination {combination_id}: {e}")
                    continue
                
                # Generation of interaction fingerprints (IFP)
                try:
                    interactions = ChemokinePartnerInteraction.objects.filter(
                        structure=s,
                        chemokine_binding_partner=chemokine_partner,
                    ).select_related('chemokine_residue')

                    interaction_dict = {}
                    for interaction in interactions:
                        ccn_number = interaction.chemokine_residue.ccn_number
                        if ccn_number not in interaction_dict:
                            interaction_dict[ccn_number] = set()
                        interaction_dict[ccn_number].add(interaction.interaction_type)

                    # Build the fingerprint bitstring
                    POSSIBLE_INTERACTIONS = [
                        "Hydrophobic", "Pistacking", "HBDonor", "HBAcceptor",
                        "Anionic", "Cationic", "CationPi", "PiCation", "VdWContact"
                    ]

                    ifp_ccn_numbers = get_ifp_ccn_numbers()
                    
                    print(f"IFP CCN numbers {ifp_ccn_numbers}")

                    bitstring_parts = []
                    for ccn_number in ifp_ccn_numbers:
                        interactions = interaction_dict.get(ccn_number, set())
                        bitstring = ["0"] * len(POSSIBLE_INTERACTIONS)
                        for interaction in interactions:
                            if interaction in POSSIBLE_INTERACTIONS:
                                index = POSSIBLE_INTERACTIONS.index(interaction)
                                bitstring[index] = "1"
                        bitstring_parts.append("".join(bitstring))
                    ifp_string = "".join(bitstring_parts)
                    print(ifp_string)

                    ChemokinePartnerIFP.objects.create(
                        structure=s,
                        binding_pair=chemokine_partner,
                        ifp_string=ifp_string
                    )
                    print(f"✅ IFP generated for {combination_id}")
                except Exception as e:
                    self.logger.error(f"Failed to generate IFP for {combination_id}: {e}")


                finally:
                    if 'u' in locals() and u is not None:
                        del u
                    if 'protein_mol' in locals() and protein_mol is not None:
                        del protein_mol
                    gc.collect()

                self.logger.info(f"Processed Structure for PDB {pdb_code} combination {combination_id}")
                                         
                s.save()
