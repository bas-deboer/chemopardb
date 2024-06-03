from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.db import transaction

from protein.models import Protein
from structure.models import Structure, PDBData, Structure_PDB, StructureType, Chain, Sequence, Residue


import os
import requests
import pandas as pd
import requests


import csv
import urllib
from datetime import datetime
from Bio import PDB
from Bio import SeqIO
from Bio.PDB import Polypeptide, PDBParser, Superimposer
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBIO import PDBIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqUtils import seq1

from copy import deepcopy
import urllib.request
import urllib.error
import json
from rcsbsearchapi.search import TextQuery
from rcsbsearchapi import rcsb_attributes as attrs
from tqdm import tqdm


def ensure_list(value):
    if not isinstance(value, list):
        return [value]
    return value
 
def fetch_pdb_ids_by_pfam(pfam_accession):
    """
    Fetches PDB IDs associated with a given PFAM accession number.

    Parameters:
    - pfam_accession: PFAM accession number as a string.

    Returns:
    - A list of PDB IDs.
    """
    url = "https://www.ebi.ac.uk/pdbe/search/pdb/select"
    query = f"pfam_accession:{pfam_accession}"
    params = {
        'q': query,
        'wt': 'json',
        'rows': 1000
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        data = response.json()
        return [doc.get('pdb_id') for doc in data['response']['docs'] if doc.get('pdb_id')]
    else:
        print("Failed to fetch data from PDBe API")
        return []

def download_mmcif_files(pdb_ids, output_dir):
    """
    Downloads mmCIF files for a list of PDB IDs.

    Parameters:
    - pdb_ids: List of PDB IDs to download.
    - output_dir: Directory where the mmCIF files will be saved.
    """
    base_url = "http://www.ebi.ac.uk/pdbe/entry-files/download/"
    for pdb_id in pdb_ids:
        url = f"{base_url}{pdb_id}_updated.cif"
        response = requests.get(url)
        if response.status_code == 200:
            file_path = os.path.join(output_dir, f"{pdb_id}.cif")
            with open(file_path, 'wb') as file:
                file.write(response.content)
            print(f"Downloaded {pdb_id}.cif")
        else:
            print(f"Failed to download {pdb_id}.cif. Status code: {response.status_code}")

def get_entry_from_api(pdb_id):
    api_url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/"

    response = requests.get(f"{api_url}{pdb_id}")
    if response.status_code == 200:
        data = response.json()
        # Adjusted to navigate the JSON structure correctly
        if pdb_id in data:
            entry_data = data[pdb_id][0]
            
            # Extracting the experimental_method_class
            method = entry_data.get('experimental_method_class')
            date = entry_data.get('release_date')

    return date

def get_comname_from_api(unp_id):
    api_url = "https://rest.uniprot.org/uniprotkb/"
    response = requests.get(f"{api_url}{unp_id}")
    if response.status_code == 200:
        data = response.json()
        # Assuming we are extracting the protein's recommended name
        if 'proteinDescription' in data and 'recommendedName' in data['proteinDescription'] and 'fullName' in data['proteinDescription']['recommendedName']:
            fullName_value = data['proteinDescription']['recommendedName']['fullName']['value']
            return fullName_value
    return None 

def find_closest(ref_seqs, query_seq):

    best_score = 0
    best_id = None
    best_alignment = None

    for ref_id, ref_seq in ref_seqs.items():
        # Perform pairwise alignment
        alignments = pairwise2.align.localms(ref_seq, query_seq, 3, 1, -3, -.1, one_alignment_only=True)
        
        if alignments:
            score = alignments[0][2]
            if score > best_score:
                best_score = score
                best_id = ref_id
                best_alignment = alignments[0]  # Store the best alignment

    return best_id, best_score, best_alignment


def align_to_ref(ref_seq, query_seq, ident_score=4, sim_score=2, gap_open=-2, gap_ext=-0.5):
    pw = pairwise2.align.localms(ref_seq, query_seq, ident_score, sim_score, gap_open, gap_ext, one_alignment_only=True)
    if pw:
        return pw[0][1]  # Return the aligned sequence
    return None


def map_sequence_to_aligned(sequence, aligned_sequence):
    """
    Map sequence positions to their corresponding positions in the aligned sequence.
    This function accounts for gaps in the aligned sequence.
    
    Parameters:
    - sequence: The original sequence (string) without gaps.
    - aligned_sequence: The aligned sequence (string) that may contain gaps ('-').
    
    Returns:
    A dictionary where keys are sequence positions (1-based) and values are the positions
    in the aligned sequence (1-based), ignoring gaps.
    """
    seq_index = 1
    aligned_index = 1
    seq_to_aligned_index_map = {}
    
    for aligned_char in aligned_sequence:
        if aligned_char != '-':
            if seq_index <= len(sequence) and sequence[seq_index-1] == aligned_char:
                seq_to_aligned_index_map[seq_index] = aligned_index
                seq_index += 1
        aligned_index += 1
    
    return seq_to_aligned_index_map



class Command(BaseBuild):
    help = "Downloads PDB files in mmCIF format and processes them."
  
    pdbs = []
    chains = {}
    seqnums = {}
    protein_seqs_aligned = {}
    gns = []
    pdbs_path = os.sep.join([settings.DATA_DIR, 'raw_pdbs'])
    af_path = os.sep.join([settings.DATA_DIR, 'af_structures'])

    def add_arguments(self, parser):
        parser.add_argument("--debug", default=False, action="store_true", help="Debug mode")
        parser.add_argument("--purge", default=False, action="store_true", help="Purge data")

    def handle(self, *args, **options):
        self.stdout.write("Command started")
        if options['purge']:
            PDBData.objects.all().delete()
            #Publication.objects.all().delete()
            #StructureDomain.objects.all().delete()
            Chain.objects.all().delete()
            Structure.objects.all().delete()
            StructureType.objects.all().delete()

        #self.download_structures()

        #self.parse_mmcif()
        
        self.parse_and_store_sequences()
        
        self.align_to_MSA()
        
        #self.parse_and_store_sequences_and_residues()


    def download_structures(self):
        """ 
        Searches PDBe for structures annotated with the chemokine Pfam family (PF00048)
        and downloads them in specified folder.
        """
        
        pfam_accession = 'PF00048'  
        print(f"Fetching PDB IDs for PFAM accession: {pfam_accession}")
        pdb_ids = fetch_pdb_ids_by_pfam(pfam_accession)
        print(f"Found {len(pdb_ids)} PDB IDs.")
    
        if pdb_ids:
            output_dir = os.path.join(settings.DATA_DIR, 'raw_cifs')
            os.makedirs(output_dir, exist_ok=True)
            download_mmcif_files(pdb_ids, output_dir)
        else:
            print("No PDB IDs found or failed to fetch.")

        #def download_protonated_structure(pdb_id, download_dir):
        #    # Ensure the directory exists
        #    os.makedirs(download_dir, exist_ok=True)
        #    response = requests.get(f'https://www.ebi.ac.uk/pdbe/model-server/v1/{pdb_id}/full?model_nums=1%2C2&encoding=cif&copy_all_categories=true&download=true&data_source=pdb-h')
        #    cif_path = os.path.join(download_dir, f'{pdb_id}.cif')
#
        #    with open(cif_path, 'wb') as fp:
        #        fp.write(response.content)
#
        #    return cif_path
#
        ## Specify your desired directory here
        #download_dir = os.path.join(settings.DATA_DIR, 'protonated_cifs_modelserver')
#
        ## Download a protonated quaternary structure to the specified directory
        #for pdb_id in tqdm(pdb_ids, desc="Downloading structures", leave=True):
        #    download_protonated_structure(pdb_id, download_dir)


    def parse_mmcif(self):
        """ 
        Reads directory of downloaded structures in mmCIF format to parse information of the structures.
        Gets general information on the structure and annotation of all entities in the structure.
        """
        
        def extract_general_info(mmcif_dict):
            method = mmcif_dict.get('_exptl.method', [None])[0]
            doi = mmcif_dict.get('_citation.pdbx_database_id_DOI', [None])[0]
            state = mmcif_dict.get('_pdbx_struct_assembly.oligomeric_details', [None])[0]
            resolution_value = mmcif_dict.get('_refine.ls_d_res_high', [None])[0]
            try:
                resolution = float(resolution_value) if resolution_value else None
            except ValueError:
                resolution = None
            return method, doi, state, resolution

        def extract_full_chain_entity_info(mmcif_file):
            mmcif_dict = MMCIF2Dict(mmcif_file)

            return {
                'asym_ids': ensure_list(mmcif_dict.get("_struct_asym.id", [])),
                'chain_to_entity_mapping': dict(zip(ensure_list(mmcif_dict.get("_struct_asym.id", [])), ensure_list(mmcif_dict.get("_struct_asym.entity_id", [])))),
                'entity_descriptions': {eid: (desc, etype) for eid, desc, etype in zip(ensure_list(mmcif_dict.get("_entity.id", [])), ensure_list(mmcif_dict.get("_entity.pdbx_description", [])), ensure_list(mmcif_dict.get("_entity.type", [])))},
                'entity_to_db_accession': dict(zip(ensure_list(mmcif_dict.get("_struct_ref.entity_id", [])), ensure_list(mmcif_dict.get("_struct_ref.pdbx_db_accession", [])))),
                'chain_to_pfam': {asym_id: acc for asym_id, xref_db, acc in zip(ensure_list(mmcif_dict.get("_pdbx_sifts_xref_db_segments.asym_id", [])), ensure_list(mmcif_dict.get("_pdbx_sifts_xref_db_segments.xref_db", [])), ensure_list(mmcif_dict.get("_pdbx_sifts_xref_db_segments.xref_db_acc", []))) if xref_db == "Pfam"},
            }
            
        # Define a function to determine if two chains are in contact
        def are_chains_in_contact(chain1, chain2, threshold=5.0):
            for atom1 in chain1.get_atoms():
                for atom2 in chain2.get_atoms():
                    if atom1 - atom2 < threshold:
                        return True
            return False

        def update_contacts_chemokine(structure, chains_queryset):
            # Identify chemokine chains
            chemokine_chains = [chain for chain in structure.get_chains() if chains_queryset.get(chain=chain.id).is_chemokine]
            # For all other chains, check contact with any chemokine chain
            for chain in structure.get_chains():
                chain_db = chains_queryset.get(chain=chain.id)
                if not chain_db.is_chemokine:  # Skip chemokine chains
                    for chemokine_chain in chemokine_chains:
                        if are_chains_in_contact(chain, chemokine_chain):
                            chain_db.contacts_chemokine = True
                            chain_db.save()
                            break
                            
        Structure.objects.all().delete()
        Chain.objects.all().delete()

        pdb_directory = os.path.join(settings.DATA_DIR, 'raw_cifs')
        parser = MMCIFParser(auth_chains=False, QUIET=True)

        for pdb_file in tqdm([f for f in os.listdir(pdb_directory) if f.endswith(".cif")], desc="Processing mmCIF files"):
            pdb_id = pdb_file[:-4]
            pdb_file_path = os.path.join(pdb_directory, pdb_file)
            structure = parser.get_structure(pdb_id, pdb_file_path)
            mmcif_dict = MMCIF2Dict(pdb_file_path)

            # Extract general structure information
            method, doi, state, resolution = extract_general_info(mmcif_dict)
            date_str = get_entry_from_api(pdb_id)  # Make sure to implement or replace `get_entry_from_api` as necessary
            date = datetime.strptime(date_str, "%Y%m%d").date()

            # Initialize a Structure instance
            pdb_info = Structure(
                pdb_id=pdb_id,
                method=method,
                doi=doi,
                resolution=resolution,
                date=date,
                state=state,
                protein=None,
            )
            pdb_info.save()

            chain_info = extract_full_chain_entity_info(pdb_file_path)
            accessions_checked = set()

            # Fetch all Protein objects once to reduce database queries
            all_proteins = Protein.objects.all()
            # Create a set of uniprot_ids for quick lookup
            protein_uniprot_ids = set(all_proteins.values_list('uniprot_id', flat=True))

            for chain_id in chain_info['asym_ids']:
                entity_id = chain_info['chain_to_entity_mapping'].get(chain_id, None)
                description, entity_type = chain_info['entity_descriptions'].get(entity_id, ("No description", "Unknown"))
                pdbx_accession = chain_info['entity_to_db_accession'].get(entity_id, None)
                pfam_accession = chain_info['chain_to_pfam'].get(chain_id, None)

                # Determine if the chain is a chemokine based on pdbx_accession
                is_chemokine = pdbx_accession in protein_uniprot_ids

                # Keep track of pdbx_accessions checked to avoid duplicate database queries
                accessions_checked.add(pdbx_accession)

                # Create and save the Chain instance with chemokine field set based on pdbx_accession
                Chain.objects.create(
                    structure=pdb_info,
                    chain=chain_id,
                    type=entity_type,
                    pdbx_accession=pdbx_accession,
                    pfam_accession=pfam_accession,
                    name=description,
                    is_chemokine=is_chemokine,
                )
                
            # Iterate over BioPython Chain objects for residue processing
            for chain in structure.get_chains():
                chain_id = chain.id
                db_chain = Chain.objects.get(structure=pdb_info, chain=chain_id)  # Get the corresponding Chain object from DB
                residues_to_create = []  # Prepare a list to hold Residue instances for bulk creation

                for residue in chain.get_residues():
                    if residue.id[0] == ' ':
                        residue_id = f"{residue.id[1]}{residue.id[2].strip()}"
                        amino_acid = seq1(residue.resname, custom_map={"UNK": "X"})  # Handle unknown amino acids as 'X'
                        # Prepare Residue instance (don't save yet)
                        residues_to_create.append(Residue(
                            chain=db_chain,
                            residue_id=residue_id,
                            amino_acid=amino_acid
                        ))

                # Bulk create Residue instances for current chain
                Residue.objects.bulk_create(residues_to_create)
                

            # After all chains are processed, try to link the structure to a Protein
            if not pdb_info.protein:
                for accession in accessions_checked:
                    if accession:
                        matching_protein = all_proteins.filter(uniprot_id=accession).first()
                        if matching_protein:
                            pdb_info.protein = matching_protein
                            pdb_info.save()
                            break  # Stop searching once a matching protein is found and assigned

            # BioPython structure parsing
            structure = parser.get_structure(pdb_info.pdb_id, pdb_file_path)

            # Fetch the queryset of all chains for this structure to update contacts
            chains_queryset = Chain.objects.filter(structure=pdb_info)

            # Update the contacts_chemokine field for chains based on contact calculation
            update_contacts_chemokine(structure, chains_queryset)

        print("Data saved to the database")


    def parse_and_store_sequences(self):
        """ 
        Reads all downloaded .cif files, extracts and stores the sequences for all chemokine chains
        """
        
        Sequence.objects.all().delete()
        d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
                    }
        
        cif_directory = os.path.join(settings.DATA_DIR, 'raw_cifs')
        parser = MMCIFParser(auth_chains=False, QUIET=True)

        for filename in os.listdir(cif_directory):
            if filename.endswith(".cif"):
                cif_path = os.path.join(cif_directory, filename)
                try:
                    structure = parser.get_structure(filename[:-4], cif_path)
                    pdb_id = filename[:-4].lower()

                    model = structure[0]
                    for chain in model:
                        chain_id = chain.id
                        sequence = []
                        sequence = ''.join(d3to1[residue.resname] for residue in chain.get_residues() if residue.id[0] == ' ' )

                        try:
                            chain_obj = Chain.objects.get(structure__pdb_id=pdb_id, chain=chain_id)

                            if chain_obj.is_chemokine:
                                seq_obj, created = Sequence.objects.get_or_create(
                                    chain=chain_obj,
                                    defaults={'sequence': sequence, 'aligned_sequence': ''}
                                )
                                if not created:
                                    seq_obj.sequence = sequence
                                    seq_obj.save()
                                    print(f"Updated sequence for chain {chain_id} in structure {pdb_id}.")
                                else:
                                    print(f"Stored sequence for chain {chain_id} in structure {pdb_id}.")
                        except Chain.DoesNotExist:
                            print(f"Chain {chain_id} in structure {pdb_id} not found in database.")

                except Exception as e:
                    print(f"Error processing {filename}: {e}")


    def align_to_MSA(self):
        """ 
        Aligns all chemokine sequences to the reference MSA.
        """
        
        def load_reference_sequences(msa_path):
            """Load reference sequences from an MSA file into a dictionary."""
            ref_seqs = {}
            for record in SeqIO.parse(msa_path, "clustal"):
                ref_seqs[record.id] = str(record.seq)
            return ref_seqs

        def update_sequences_with_alignment(cif_directory, msa_path):
            ref_seqs = load_reference_sequences(msa_path)

            for chain in Chain.objects.filter(is_chemokine=True):
                if not chain.sequence_set.exists():
                    print(f"No sequence for chain {chain.chain}. Skipping...")
                    continue
                
                query_seq = chain.sequence_set.first().sequence
                best_id, _, best_alignment = find_closest(ref_seqs, query_seq)

                if best_alignment and best_id:
                    # Ensure you're using the reference sequence corresponding to the best ID
                    ref_seq = ref_seqs[best_id]
                    aligned_seq = align_to_ref(ref_seq, query_seq)
                    if aligned_seq:
                        # Update the aligned_sequence field in the Sequence model
                        chain.sequence_set.update_or_create(
                            defaults={'aligned_sequence': aligned_seq}
                        )
                        print(f"Updated aligned sequence for chain {chain.chain}.")
                    else:
                        print(f"Alignment failed for chain {chain.chain}.")
                else:
                    print(f"No close reference found for chain {chain.chain}.")

        # Example usage
        cif_directory = os.path.join(settings.DATA_DIR, 'raw_cifs')
        msa_path = os.path.join(settings.DATA_DIR, 'chemokines_ref_MSA.aln')
        update_sequences_with_alignment(cif_directory, msa_path)


    def parse_and_store_sequences_and_residues(self):
        cif_directory = os.path.join(settings.DATA_DIR, 'raw_cifs')
        parser = MMCIFParser(QUIET=True, auth_chains=False)

        d3to1 = {
            'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
        }

        for filename in os.listdir(cif_directory):
            if filename.endswith(".cif"):
                cif_path = os.path.join(cif_directory, filename)
                try:
                    structure = parser.get_structure(filename[:-4], cif_path)
                    pdb_id = filename[:-4].lower()

                    try:
                        structure_obj = Structure.objects.get(pdb_id=pdb_id)
                    except Structure.DoesNotExist:
                        print(f"Structure {pdb_id} not found in the database.")
                        continue
                    
                    for model in structure:
                        for chain in model:
                            chain_id = chain.id

                            try:
                                chain_obj = Chain.objects.get(structure=structure_obj, chain=chain_id)
                                sequence_obj = Sequence.objects.filter(chain=chain_obj).first()

                                if chain_obj.is_chemokine and sequence_obj:
                                    aligned_sequence = sequence_obj.aligned_sequence
                                    sequence_str = ''.join(d3to1[residue.resname] for residue in chain.get_residues() if residue.id[0] == ' ')
                                    seq_to_aligned_index_map = map_sequence_to_aligned(sequence_str, aligned_sequence)

                                    with transaction.atomic():
                                        for residue in chain.get_residues():
                                            if residue.id[0] == ' ' and residue.resname in d3to1:
                                                sequence_number = residue.id[1]
                                                generic_number = seq_to_aligned_index_map.get(sequence_number, 0)
                                                Residue.objects.create(
                                                    structure=structure_obj,
                                                    chain=chain_obj,
                                                    amino_acid=d3to1[residue.resname],
                                                    amino_acid_three_letter=residue.resname,
                                                    sequence_number=sequence_number,
                                                    generic_number=generic_number
                                                )
                            except Chain.DoesNotExist:
                                print(f"Chain {chain_id} in structure {pdb_id} not found in the database.")
                except Exception as e:
                    print(f"Error processing {filename}: {e}")

        


# from protwis
    def build_interactions(self, pdb_code):
        try:
            # interacting_pairs, distances  = compute_interactions(pdb_code, save_to_db=True)
            compute_interactions(pdb_code, protein=None, lig=None, do_interactions=True, do_complexes=False, do_peptide_ligand=True, save_to_db=True, file_input=False)
            # compute_interactions(pdb_code, do_interactions=True, do_peptide_ligand=True, save_to_db=True)
        except:
            self.logger.error('Error with computing interactions (%s)' % (pdb_code))
            return













#####################################################################################################


    def download_alphafold(self):
        """
        Load UniProt IDs from a specified CSV file.
        Args: 
            uniprot_file (str): Path to the file containing UniProt IDs.
        Returns: 
            list: A list of UniProt IDs.
        """
        
        print("Reading Pfam protein data")
        uniprot_list = []
        with open(os.sep.join([settings.DATA_DIR, 'pfam_protein_data.csv']), newline='') as csvfile:
            uniprot_reader = csv.reader(csvfile, delimiter=',')
            for row in uniprot_reader:
                if row[3] not in uniprot_list:
                    uniprot_list.append(row[3])
        print(uniprot_list)
            
        af_directory = os.sep.join([settings.DATA_DIR, 'af_structures'])
        os.makedirs(af_directory, exist_ok=True)

        for uniprot_id in uniprot_list:
            alphafold_id = f'AF-{uniprot_id}-F1'
            model_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_id}-model_v4.pdb'
            file_path = os.path.join(af_directory, f'{alphafold_id}.pdb')
            print(f"downloading {alphafold_id}")
            os.system(f'curl {model_url} -o "{file_path}"')


    def build_alphafold(self):
        af_type, created = StructureType.objects.get_or_create(slug='AF', name='Alphafold')
        month_dict = {'JAN':'01', 'FEB':'02', 'MAR':'03', 'APR':'04', 'MAY':'05', 'JUN':'06', 'JUL':'07', 'AUG':'08', 'SEP':'09', 'OCT':'10', 'NOV':'11', 'DEC':'12'}

        for f in os.listdir(self.af_path):
            this_file = os.sep.join([self.af_path, f])
            accession = f.split('-')[1]
            print(accession)
            protein = Protein.objects.get(uniprot_id=accession)
            with open(this_file, 'r') as f_handle:
                lines = f_handle.readlines()
            day, month, year = lines[0].split('HEADER')[1].strip().split('-')
            model_date = '20{}-{}-{}'.format(year, month_dict[month], day)
            structure, created = Structure.objects.get_or_create(pdb_id=None, date=model_date, doi=None, resolution=None,
                                                                 structure_type=af_type, protein=protein)
            chain, created = Chain.objects.get_or_create(chain_ID='A', structure=structure)
            
            # Store the PDB data
            pdbdata = self.build_pdbdata(this_file)
            name = this_file.replace(".pdb", "")
            structure_pdb, created = Structure_PDB.objects.get_or_create(name=name, pdbdata=pdbdata)
   
                
    def build_structure_residues(self):
        def get_structure_id(pdb_id):
            if pdb_id != None:
                try:
                    structure = Structure.objects.get(pdb_id=pdb_id)
                    return structure.id
                except Structure.DoesNotExist:
                    return None        
        
        # Mapping from three-letter to one-letter amino acid codes
        AMINO_ACID_MAPPING = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        
        directory = os.path.join(settings.DATA_DIR, 'raw_pdbs')
        
        parser = PDB.PDBParser()
        for filename in os.listdir(directory):
            if filename.endswith(".pdb"):
                pdb_id = filename[:4]  # Assumes the file name starts with a 4-character PDB ID
                structure_id = get_structure_id(pdb_id)
                if structure_id is None:
                    continue  # Skip if no matching structure is found
                
                filepath = os.path.join(directory, filename)
                structure = parser.get_structure('PDB_structure', filepath)
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.id[0] == ' ' and residue.resname in AMINO_ACID_MAPPING:  # Filter out hetero-atoms and unknown residues
                                amino_acid = AMINO_ACID_MAPPING[residue.resname]
                                sequence_number = residue.id[1]
                                Residue.objects.create(
                                    amino_acid=amino_acid,
                                    structure_id=structure_id,
                                    chain=chain.id,
                                    sequence_number=sequence_number
                                )


    def build_pdbdata(self, filename):
        with open(filename, 'r') as f:
            pdbdata = f.read()
            
        pdbdata, created = PDBData.objects.get_or_create(pdb=pdbdata)
        return pdbdata


    def build_models(self):
        p = PDBParser()

        for filename in os.listdir(self.pdbs_path):
            if filename.endswith(".pdb"):
                model_filename = os.path.join(self.pdbs_path, filename)
                
                pdbdata = self.build_pdbdata(model_filename)
                name = filename.replace(".pdb", "")
                
                # Store the PDB data
                structdomain, created = Structure_PDB.objects.get_or_create(name=name, pdbdata=pdbdata)


    def structure_state(self):
        file_path = os.path.join(settings.DATA_DIR, 'ChemoPar_partner_list.xlsx')

        # Read the Excel file
        df = pd.read_excel(file_path)

        # Iterate through the DataFrame
        for index, row in df.iterrows():
            pdb_id = row['PDB']  # Replace 'pdb_id' with the actual column name in your Excel file
            state = row['State']    # Replace 'state' with the actual column name for the state

            # Update the Structure instance
            try:
                structure = Structure.objects.get(pdb_id=pdb_id)
                structure.state = state  # Replace 'state' with the actual field name in your Structure model
                structure.save()
                self.stdout.write(self.style.SUCCESS(f'Successfully updated {pdb_id}'))
            except Structure.DoesNotExist:
                self.stdout.write(self.style.WARNING(f'Structure with pdb_id {pdb_id} not found'))

## FROM SH2DB function:   
#    def pdb_request_by_pdb(self, pdb):
#        data = {}
#        try:
#            response = urllib.request.urlopen('https://data.rcsb.org/rest/v1/core/entry/{}'.format(pdb))
#        except urllib.error.HTTPError:
#            print('Error: {} PDB ID not found on RCSB.'.format(pdb))
#            return data
#        json_data = json.loads(response.read())
#        response.close()
#        # pprint.pprint(json_data)
#        data['method'] = json_data['exptl'][0]['method']
#        data['journal'] = json_data['citation'][0]['rcsb_journal_abbrev']
#        if 'pdbx_database_id_doi' in json_data['citation'][0]:
#            data['doi'] = json_data['citation'][0]['pdbx_database_id_doi']
#        else:
#            data['doi'] = None
#        data['authors'] = json_data['citation'][0]['rcsb_authors']
#        data['title'] = json_data['citation'][0]['title']
#        if 'year' in json_data['citation'][0]:
#            data['year'] = json_data['citation'][0]['year']
#        else:
#            data['year'] = None
#        if 'resolution_combined' in json_data['rcsb_entry_info']:
#            data['resolution'] = json_data['rcsb_entry_info']['resolution_combined'][0]
#        else:
#            data['resolution'] = None
#        data['publication_date'] = json_data['rcsb_accession_info']['initial_release_date'][:10]
#        return data