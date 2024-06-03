from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import transaction

from protein.models import Protein
from structure.models import Structure, Chain, ChemokineInteraction
from model.models import Model

import os
import io
import pandas as pd
from Bio import PDB
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from pathlib import Path
import requests
import sys
import time
from urllib.parse import urljoin
import warnings
import glob
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio import pairwise2
from Bio.SeqUtils import seq1
from Bio.PDB import MMCIFParser, MMCIFIO
from Bio.SeqUtils import seq1
import os
from tqdm import tqdm
from Bio.PDB import MMCIFParser, MMCIFIO
from Bio.PDB import FastMMCIFParser, MMCIFIO
from Bio.SeqUtils import seq1
from io import StringIO
from Bio import AlignIO
from Bio.Align import PairwiseAligner
from tqdm import tqdm
import os
from Bio.PDB import MMCIFParser, PDBIO
import pdbfixer
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from moleculekit.molecule import Molecule
from moleculekit.tools.preparation import systemPrepare
import pdbfixer
from pdbfixer import PDBFixer
from openmm.app import PDBFile

import os
from django.db.models import Q
from Bio.PDB import MMCIFParser, PDBIO, NeighborSearch
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import Select


def are_chains_in_contact(chain1, chain2, threshold=5.0):
    for atom1 in chain1.get_atoms():
        for atom2 in chain2.get_atoms():
            if atom1 - atom2 < threshold:
                return True
    return False

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


class ChemokineContactSelect(Select):
    """
    This class allows selective writing of chains.
    Only the chemokine chain and chains in contact with it are written.
    """
    def __init__(self, chemokine_chain_id, contacting_chain_ids):
        self.chemokine_chain_id = chemokine_chain_id
        self.contacting_chain_ids = contacting_chain_ids

    def accept_chain(self, chain):
        # Accept the chain if it's the chemokine chain or one of the contacting chains
        return chain.id == self.chemokine_chain_id or chain.id in self.contacting_chain_ids

class Command(BaseBuild):

    def add_arguments(self, parser):
        parser.add_argument("--debug", default=False, action="store_true", help="Debug mode")
        parser.add_argument("--purge", default=False, action="store_true", help="Purge data")

    def handle(self, *args, **options):
        #Model.objects.all().delete()
        #if options['purge']:
        #    PDBData.objects.all().delete()
        #    Publication.objects.all().delete()
        #    StructureDomain.objects.all().delete()
        #    Chain.objects.all().delete()
        #    Structure.objects.all().delete()
        #    StructureType.objects.all().delete()

        # Download protonated models from ModelServer
        #self.cif_to_pdb()

        #self.generate_chemokine_pdb_files()
        #self.fix_cif_files()
        self.protoss()
        
        # Split chains
        #self.split_chains()
        
        # Calculate chemokine similarity
        #self.chemokine_similarity()
        
        # Protonation
        #self.protoss()
        
        # Split chains
        #self.split_model_to_chains()

        # Filter out chemokine chains
        #self.filter_chemokines()
        


    def generate_chemokine_pdb_files(self):
        cif_directory = os.path.join(settings.DATA_DIR, 'raw_cifs')
        output_directory = os.path.join(settings.DATA_DIR, 'output_pdb2')

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        for cif_file in os.listdir(cif_directory):
            if cif_file.endswith(".cif"):
                pdb_id = cif_file[:-4]
                structure_records = Structure.objects.filter(pdb_id=pdb_id)

                for structure_record in structure_records:
                    chemokine_chains = Chain.objects.filter(structure=structure_record, is_chemokine=True)
                    if chemokine_chains.exists():
                        parser = MMCIFParser(QUIET=True, auth_chains=False)
                        structure = parser.get_structure(pdb_id, os.path.join(cif_directory, cif_file))

                        for chemokine_chain_record in chemokine_chains:
                            chemokine_chain = structure[0][chemokine_chain_record.chain]

                            contacting_chain_ids = []
                            for chain in structure[0]:
                                if chain.id != chemokine_chain.id and are_chains_in_contact(chain, chemokine_chain):
                                    contacting_chain_ids.append(chain.id)

                            io = MMCIFIO()
                            select = ChemokineContactSelect(chemokine_chain.id, contacting_chain_ids)
                            output_path = os.path.join(output_directory, f"{pdb_id}_{chemokine_chain_record.chain}.cif")
                            io.set_structure(structure)
                            io.save(output_path, select=select)

                            # Create and save the ChemokineInteraction object
                            chemokine_interaction = ChemokineInteraction(
                                structure=structure_record,
                                chemokine_chain=chemokine_chain_record,
                                output_pdb_path=output_path
                            )
                            chemokine_interaction.save()

                            # Add the contacting chains
                            contacting_chains = Chain.objects.filter(structure=structure_record, chain__in=contacting_chain_ids)
                            for chain in contacting_chains:
                                chemokine_interaction.contacting_chains.add(chain)
        
        
    def fix_cif_files(self):
        # Define the source and target directories
        source_dir = os.path.join(settings.DATA_DIR, 'output_pdb2')
        target_dir = os.path.join(settings.DATA_DIR, 'output_pdb2_fixed')

        # Ensure the target directory exists
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        # Iterate over all .cif files in the source directory
        for filename in os.listdir(source_dir):
            if filename.endswith(".cif"):
                try:
                    # Construct the full path of the source file
                    pdb_file_path = os.path.join(source_dir, filename)

                    # Initialize PDBFixer with the source file
                    fixer = pdbfixer.PDBFixer(filename=pdb_file_path)
                    fixer.findMissingResidues()
                    fixer.missingResidues = {}
                    fixer.findMissingAtoms()
                    fixer.addMissingAtoms()

                    # Construct the output file path, replacing .cif with _output.pdb
                    output_file_path = os.path.join(target_dir, filename.replace(".cif", "_output.pdb"))

                    # Write the fixed structure to the output file
                    with open(output_file_path, 'w') as f:
                        PDBFile.writeFile(fixer.topology, fixer.positions, f)

                    print(f"Processed {filename} to {output_file_path}")
                except Exception as e:
                    print(f"Error processing file {filename}: {e}")
                    
                    
    def protoss(self):
        
        def poll_job(job_id, poll_url, poll_interval=1, max_polls=10):
            """Poll the progress of a job

            Continuosly polls the server in regular intervals and updates the job information, especially the status.

            :param job_id: UUID of the job to poll
            :type job_id: str
            :param poll_url: URl to send the polling request to
            :type poll_url: str
            :param poll_interval: time interval between polls in seconds
            :type poll_interval: int
            :param max_polls: maximum number of times to poll before exiting
            :type max_polls: int
            :return: polled job
            :rtype: dict
            """
            job = requests.get(poll_url + job_id + '/').json()
            status = job['status']
            current_poll = 0
            while status == 'pending' or status == 'running':
                print(f'Job {job_id} is { status }')
                current_poll += 1
                if current_poll >= max_polls:
                    print(f'Job {job_id} has not completed after {max_polls} polling requests' \
                          f' and {poll_interval * max_polls} seconds')
                    return job
                time.sleep(poll_interval)
                job = requests.get(poll_url + job_id + '/').json()
                status = job['status']
            print(f'Job {job_id} completed with { status }')
            return job

        # constants
        PROTEINS_PLUS_URL = 'https://proteins.plus/api/v2/'
        UPLOAD = urljoin(PROTEINS_PLUS_URL, 'molecule_handler/upload/')
        UPLOAD_JOBS = urljoin(PROTEINS_PLUS_URL, 'molecule_handler/upload/jobs/')
        PROTEINS = urljoin(PROTEINS_PLUS_URL, 'molecule_handler/proteins/')
        LIGANDS = urljoin(PROTEINS_PLUS_URL, 'molecule_handler/ligands/')
        PROTOSS = urljoin(PROTEINS_PLUS_URL, 'protoss/')
        PROTOSS_JOBS = urljoin(PROTEINS_PLUS_URL, 'protoss/jobs/')

        pdb_directory = os.path.join(settings.DATA_DIR, 'output_pdb2_fixed')
        output_directory = os.path.join(settings.DATA_DIR, 'output_pdb2_fixed_protoss')

        # Create the output directory if it doesn't exist
        os.makedirs(output_directory, exist_ok=True)

        # List all files in the directory
        pdb_files = os.listdir(pdb_directory)

        # Loop over each PDB file in the directory
        for pdb_file in pdb_files:
            # Construct the full path to the PDB file
            pdb_file_path = os.path.join(pdb_directory, pdb_file)

            # Check if the file is a regular file (not a directory)
            if os.path.isfile(pdb_file_path):
                try:
                    # Your existing code goes here
                    with open(pdb_file_path, 'rb') as upload_file:
                        query = {'protein_file': upload_file}
                        job_submission = requests.post(PROTOSS, files=query).json()
                    protoss_job = poll_job(job_submission['job_id'], PROTOSS_JOBS)
                    protossed_protein = requests.get(PROTEINS + protoss_job['output_protein'] + '/').json()

                    # Save the new PDB file
                    new_pdb_file_path = os.path.join(output_directory, f"{pdb_file.split('.')[0]}.pdb")
                    protein_file = io.StringIO(protossed_protein['file_string'])
                    protein_structure = PDBParser().get_structure(protossed_protein['name'], protein_file)

                    # Create a PDBIO instance and save the structure
                    pdb_io = PDBIO()
                    pdb_io.set_structure(protein_structure)
                    pdb_io.save(new_pdb_file_path)

                except KeyError as e:
                    # Handle the case where 'job_id' is missing in the response
                    print(f"Error processing '{pdb_file}': {e}")
                except Exception as e:
                    # Handle other exceptions that may occur
                    print(f"Error processing '{pdb_file}': {e}")
        
        
        
        
        
        
        
        
        
        
        
        
        
        

    def split_cif_models(self):
        input_directory = os.path.join(settings.DATA_DIR, 'raw_cifs')
        output_directory = os.path.join(settings.DATA_DIR, 'all_models')
        os.makedirs(output_directory, exist_ok=True)
    
        parser = FastMMCIFParser(QUIET=True)
    
        # List all CIF files in the input directory
        cif_files = glob.glob(os.path.join(input_directory, "*.cif"))
        print(cif_files)
    
        for cif_file in cif_files:
            try:
                structure = parser.get_structure("structure", cif_file)
    
                for model in structure:
                    base_filename = os.path.basename(cif_file)
                    output_filename = f"{base_filename.rsplit('.', 1)[0]}_model_{model.id}.cif"
                    
                    io = MMCIFIO()
                    io.set_structure(model)
    
                    output_path = os.path.join(output_directory, output_filename)
                    io.save(output_path)
                    print(f"Saved {output_path}")
            except Exception as e:
                print(f"Error processing file {cif_file}: {e}")


    def cif_to_prepared_pdb(self, input_dir, output_dir):
        """
        Converts all .cif files in the input directory to .pdb files in the output directory.
        """
        os.makedirs(output_dir, exist_ok=True)
        error_files = []  # Initialize the list to track files with errors

        for file in os.listdir(input_dir):
            if file.endswith(".cif"):
                # Construct the final output file path in the output directory
                output_file_path = os.path.join(output_dir, file[:-4] + '_prepared.pdb')
            
                # Check if the output file already exists
                if os.path.exists(output_file_path):
                    print(f"Output file {output_file_path} already exists. Skipping conversion.")
                    continue
                try:
                    file_path = os.path.join(input_dir, file)

                    mol = Molecule(file_path)
                    # Prepare the molecule
                    prepared = systemPrepare(mol, verbose=False, ignore_ns_errors=True)

                    # Handle the systemPrepare output
                    if isinstance(prepared, tuple):
                        mol_prepared, df = prepared
                        # Construct report CSV path and save it in the output directory
                        report_csv_path = os.path.join(output_dir, file[:-4] + '-report.csv')
                        df.to_csv(report_csv_path)
                    else:
                        # If systemPrepare returns only a Molecule object
                        mol_prepared = prepared
                    # Construct the output file path for the prepared PDB file in the output directory
                    prepared_pdb_path = os.path.join(output_dir, file[:-4] + '_prepared.pdb')
                    # Write the prepared PDB file
                    mol_prepared.write(prepared_pdb_path)

                    # Run PDB through PDBFixer
                    #pdb_file = prepared_pdb_path
#
                    #fixer = pdbfixer.PDBFixer(pdb_file)
                    #fixer.removeHeterogens()   
                    #fixer.findMissingResidues()
                    #fixer.missingResidues = {}
                    #fixer.findNonstandardResidues()  
                    #fixer.replaceNonstandardResidues()  
                    #fixer.findMissingAtoms()  
                    #fixer.addMissingAtoms()  
                    #fixer.addMissingHydrogens(7.0)  
                    #
                    ## Construct the final output file path in the output directory
                    #output_file_path = os.path.join(output_dir, file[:-4] + '_prepared.pdb')
                    #with open(output_file_path, 'w') as output_file:
                    #    PDBFile.writeFile(fixer.topology, fixer.positions, output_file)
                    print(f"Converted {file} successfully")
                    
                except Exception as e:
                    print(f"Error converting {file}: {e}")
                    error_files.append({'file': file, 'error': str(e)})  # Append file and error message to the list
        
        # After processing all files, check if there are any errors to report
        if error_files:
            df_errors = pd.DataFrame(error_files)
            errors_csv_path = os.path.join(output_dir, 'conversion_errors.csv')
            df_errors.to_csv(errors_csv_path, index=False)
            print(f"Errors reported in {errors_csv_path}")


    def cif_to_pdb(self):
        input_directory = os.path.join(settings.DATA_DIR, 'raw_cifs')
        output_directory = os.path.join(settings.DATA_DIR, 'prepared_pdbs')
        self.cif_to_prepared_pdb(input_directory, output_directory)
        
        
    def split_pdb_by_chains(self, pdb_path, output_dir):
        parser = PDB.PDBParser()
        pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]
        structure = parser.get_structure(pdb_id, pdb_path)

        # Create a subdirectory for the current structure
        structure_dir = os.path.join(output_dir, pdb_id)
        os.makedirs(structure_dir, exist_ok=True)

        io = PDB.PDBIO()

        for model in structure:
            for chain in model:
                # Create a new structure to hold the current chain
                new_struct = PDB.Structure.Structure("chain_" + chain.id)
                new_model = PDB.Model.Model(0)
                new_model.add(chain)
                new_struct.add(new_model)

                output_file_path = os.path.join(structure_dir, f"chain_{chain.id}.pdb")
                io.set_structure(new_struct)
                io.save(output_file_path)
                self.stdout.write(self.style.SUCCESS(f'Successfully written {output_file_path}'))


    def split_chains(self):
        prepared_pdbs_dir = os.path.join(settings.DATA_DIR, 'prepared_pdbs')
        if not os.path.exists(prepared_pdbs_dir):
            raise CommandError(f"Directory '{prepared_pdbs_dir}' does not exist")

        for pdb_file in os.listdir(prepared_pdbs_dir):
            if pdb_file.endswith(".pdb"):
                pdb_path = os.path.join(prepared_pdbs_dir, pdb_file)
                self.split_pdb_by_chains(pdb_path, prepared_pdbs_dir)
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

    def parse_mmcif(self, mmcif_file):
        parser = PDB.MMCIFParser(auth_chains=False, QUIET=True)
        return parser.get_structure("structure", mmcif_file)

    #def split_chains(self):
    #    downloaded_mmcifs = os.path.join(settings.DATA_DIR, 'raw_cifs')
    #    output_dir = os.path.join(settings.DATA_DIR, 'splitted_chains')
    #    os.makedirs(output_dir, exist_ok=True)
#
    #    mmcif_directory = downloaded_mmcifs
#
    #    for filename in tqdm(os.listdir(mmcif_directory), desc="Parsing mmCIFs", leave=True):
    #        if filename.endswith(".cif"):
    #           mmcif_file_path = os.path.join(mmcif_directory, filename)
    #           structure = self.parse_mmcif(mmcif_file_path)
    #           pdb_id = os.path.splitext(filename)[0]

    #           for model_id, model in enumerate(structure):
    #               for chain in model:
    #                   chain_id = chain.id
    #                   sequence = seq1("".join([residue.resname for residue in chain.get_residues() if residue.id[0] == ' ']))  # Only standard residues
    #                   
    #                   # Prepare MMCIF content for this chain
    #                   string_io = StringIO()
    #                   io = MMCIFIO()
    #                   io.set_structure(chain)
    #                   io.save(string_io)
    #                   mmcif_chain_content = string_io.getvalue()
    #                   string_io.close()

    #                   model_instance = Model(
    #                       name=f"{pdb_id}_{model_id}_{chain_id}",
    #                       pdb_id=pdb_id,
    #                       model_id=str(model_id),
    #                       #chain_id=chain_id,
    #                       sequence=sequence,
    #                       mmcif_text=mmcif_chain_content,  # Store the MMCIF content for this chain
    #                   )
    #                   model_instance.save()

    #                   # Optional: Write chain to a separate file
    #                   file_output_dir = os.path.join(output_dir, pdb_id, f"model_{model_id}")
    #                   os.makedirs(file_output_dir, exist_ok=True)
    #                   chain_file = os.path.join(file_output_dir, f"chain_{chain_id}.cif")
    #                   with open(chain_file, 'w') as file:
    #                        file.write(mmcif_chain_content)
#
    #    self.stdout.write(self.style.SUCCESS('Successfully processed chains and stored their sequences and MMCIF contents'))























    def split_model_to_chains(self):
        input_dir = os.sep.join([settings.DATA_DIR, 'protonated_models'])
        output_dir = os.sep.join([settings.DATA_DIR, 'all_chains'])
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parser = PDB.PDBParser(QUIET=True)

        with transaction.atomic():
            for pdb_file in os.listdir(input_dir):
                if pdb_file.endswith('.pdb'):
                    filename_parts = pdb_file.split('_')  # Split the filename by underscores
                    if len(filename_parts) == 2:
                        pdb_id, model_str = filename_parts
                        model_number = int(model_str.replace("model", "").split('.')[0])  # Extract model number
                    else:
                        print("Incorrect name detected")
                        continue

                    structure = parser.get_structure(pdb_id, os.path.join(input_dir, pdb_file))

                    # Iterate through each model (conformation) in the structure
                    for model in structure:
                        for chain in model:
                            # Create a new structure with only the current chain
                            chain_structure = PDB.Structure.Structure(pdb_id)
                            chain_model = PDB.Model.Model(model.id)
                            chain_model.add(chain)
                            chain_structure.add(chain_model)

                            # Write the new structure to a PDB file in the output directory
                            output_filename = f"{pdb_id}_model{model_number}_chain{chain.id}.pdb"
                            output_path = os.path.join(output_dir, output_filename)
                            io = PDB.PDBIO()
                            io.set_structure(chain_structure)
                            io.save(output_path)

                            # Create and save a new Model instance for this chain
                            model_instance = Model(
                                name=f"{pdb_id}_model{model_number}_chain{chain.id}",
                                pdb_id=pdb_id,
                                model_id=model_number,
                                chain_id=chain.id,
                            )
                            model_instance.save()
    
    def filter_chemokines(self):
        # Define the directory containing PDB files and the path to the Excel sheet
        pdb_directory = os.sep.join([settings.DATA_DIR, 'all_chains'])
        excel_file = os.sep.join([settings.DATA_DIR, 'chemokine_filters.xlsx'])

        # Read the Excel sheet into a DataFrame
        df = pd.read_excel(excel_file)

        # Iterate through the rows of the DataFrame and delete the corresponding PDB files
        for index, row in df.iterrows():
            pdb_id = row['pdb_id']
            chain_nums = str(row['chain#']).split()  # Split the chain numbers into a list
            matching_files = []

            # Find matching files for each chain number
            for chain_num in chain_nums:
                file_pattern = f"{pdb_id}_model*_chain{chain_num}.pdb"
                matching_files.extend(glob.glob(os.path.join(pdb_directory, file_pattern)))

            # Delete all matching files
            for file_path in matching_files:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    print(f"Deleted file: {file_path}")
                else:
                    print(f"File not found: {file_path}")

        print("Finished filtering PDB files.")