from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from django.db import IntegrityError
from django.db import transaction

from common.models import WebLink, WebResource, Publication
from partner.models import PartnerProteinStructure
from protein.models import Protein, ProteinConformation, ProteinState, ProteinSegment
from structure.models import Structure, PdbData, StructureType, Chain, Rotamer, Fragment, Entity, EntityType
from residue.models import Residue
from interaction.models import *

from common.protossapi import *
from interaction.calculation2024 import *

import ast
import os
import re
import time
import django.apps
import requests
import csv
import gc
import json
import logging
import threading

from urllib.parse import urljoin
from collections import OrderedDict
from copy import deepcopy
from datetime import datetime
from urllib.request import urlopen
from urllib.error import urllib
from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder
from io import StringIO
from Bio import PDB, SeqIO, pairwise2
from Bio.PDB import (MMCIFParser, PDBIO, PDBParser, PPBuilder, Polypeptide,
                     Selection, Superimposer, parse_pdb_header)
from Bio.SeqUtils import seq1
from tqdm import tqdm
from rcsbsearchapi import rcsb_attributes as attrs
from rcsbsearchapi.search import TextQuery
from Bio import Align


def fetch_entity_info(pdb_id, entity_id, entity_type):
    """Fetch entity info from RCSB PDB API based on type."""
    url = f"https://data.rcsb.org/rest/v1/core/{entity_type}/{pdb_id}/{entity_id}"
    response = requests.get(url)
    if response.status_code == 200:
        entity_data = response.json()
        if entity_type == "polymer_entity":
            description = entity_data.get("rcsb_polymer_entity", {}).get("pdbx_description", "Unknown")
            return {
                "id": entity_data.get("rcsb_id"),
                "name": description,
                "type": entity_data.get("entity_poly", {}).get("rcsb_entity_polymer_type", "Unknown"),
                "chain": entity_data.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids", "Unknown"),
                "unp_accession": entity_data.get("rcsb_polymer_entity_container_identifiers", {}).get("uniprot_ids", "Unknown"),
            }
        elif entity_type == "nonpolymer_entity":
            description = entity_data.get("rcsb_nonpolymer_entity", {}).get("pdbx_description", "Unknown")
            return {
                "id": entity_data.get("rcsb_id"),
                "name": description,
                "type": "Non-polymer",
                "chain": entity_data.get("rcsb_nonpolymer_entity_container_identifiers", {}).get("auth_asym_ids", "Unknown"),
                "unp_accession": "Unknown",
            }
    return None

def get_pdb_entities(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        print(f"Failed to retrieve data for PDB ID: {pdb_id}")
        return None


class ParseStructureCSV:
    def __init__(self):
        self.pdb_ids = []
        self.structures = {}
        with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'chemokine_structures.csv']), newline='') as csvfile:
            structures = csv.reader(csvfile)
            next(structures, None)  # Skip header
            for s in structures:
                self.pdb_ids.append(s[0])
                self.structures[s[0]] = {
                    'protein': s[1],
                    'name': s[0].lower(),
                    'method_from_file': s[2],
                    'resolution': s[3],
                    'state': s[4],
                    'preferred_chain': s[5]
                }

    def __str__(self):
        return '<ParsedStructures: {} entries>'.format(len(self.pdb_ids))


class Command(BaseBuild):
    help = 'Reads source data and creates pdb structure objects'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
                            type=int,
                            action='store',
                            dest='proc',
                            default=1,
                            help='Number of processes to run')
        parser.add_argument('-s', '--structure',
                            dest='structure',
                            help='Structure to import (PDB ID)',
                            nargs='+')
        parser.add_argument('-u', '--purge',
                            action='store_true',
                            dest='purge',
                            default=False,
                            help='Purge existing records')
        parser.add_argument('--skip_cn',
                            action='store_true',
                            default=False,
                            help='Skip building contact network for test build')
        parser.add_argument('-i', '--incremental',
                            action='store_true',
                            dest='incremental',
                            default=False,
                            help='Incremental update to structures for small live update')
        parser.add_argument('--debug',
                            action='store_true',
                            dest='debug',
                            default=False,
                            help='Print info for debugging')

    lock = threading.Lock()
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]

    pdb_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])

    s = ProteinSegment.objects.all()
    segments = {}
    for segment in s:
        segments[segment.slug] = segment

    exp_method_dict = {'X-ray': 'X-ray diffraction', 'NMR': 'Nuclear magnetic resonance',
                       'cryo-EM': 'Electron microscopy', 'Electron crystallography': 'Electron crystallography'}

    def handle(self, *args, **options):
        self.stdout.write("Build command started...")

        if options['purge']:
            try:
                self.purge_structures()
                self.tracker = {}
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

        if options['skip_cn']:
            self.run_contactnetwork = False
        else:
            self.run_contactnetwork = True

        self.construct_errors, self.rotamer_errors, self.contactnetwork_errors, self.interaction_errors = [], [], [], []

        self.parsed_structures = ParseStructureCSV()

        if options['structure']:
            self.parsed_structures.pdb_ids = [i for i in self.parsed_structures.pdb_ids if i in options['structure'] or i.lower() in options['structure']]

        self.incremental_mode = options['incremental']
        self.debug = options['debug']

        try:
            self.logger.info('CREATING STRUCTURES')
            self.main_func()
            self.logger.info('COMPLETED CREATING STRUCTURES')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        print('Construct errors:')
        print(self.construct_errors)
        print('Rotamer errors:')
        print(self.rotamer_errors)
        print('Contact network errors')
        print(self.contactnetwork_errors)
        print('Interaction errors')
        print(self.interaction_errors)

    def purge_structures(self):
        Structure.objects.all().delete()
        PdbData.objects.all().delete()

    def get_pdb_sequence(self, pdb_data, chain_id):
        """Extract amino acid sequence from a PDB file's specified chain, handling gaps."""
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
        """Align two sequences using pairwise alignment with custom gap penalties."""
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'

        # Set gap penalties
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        alignment = aligner.align(seq1, seq2)
        
        print(alignment[0])
        
        return alignment[0]

    def create_mapping(self, alignment):
        """Create a mapping from PDB sequence positions to database sequence positions, including gaps."""
        # Extract the sequences
        target_seq = alignment[0]
        query_seq = alignment[1]

        # Initialize counters
        target_pos = 1
        query_pos = 1

        # Initialize mapping dictionary
        mapping = {}

        # Iterate through the sequences
        for t, q in zip(target_seq, query_seq):
            if t != '-':
                if q != '-':
                    mapping[target_pos] = query_pos
                target_pos += 1
            if q != '-':
                query_pos += 1

        print(mapping)
        return mapping


    def get_segment_from_generic_number(self, generic_number):
        if generic_number is not None:
            if 1 <= generic_number <= 62:
                return 'N-term'
            if 62 <= generic_number <= 67:
                return 'CX'
            if 67 <= generic_number <= 82:
                return 'N-loop'
            if 82 <= generic_number <= 89:
                return 'B1'
            if 89 <= generic_number <= 100:
                return '30s-loop'
            if 100 <= generic_number <= 105:
                return 'B2'
            if 105 <= generic_number <= 114:
                return '40s-loop'
            if 114 <= generic_number <= 118:
                return 'B3'
            if 118 <= generic_number <= 122:
                return '50s-loop'
            if 122 <= generic_number <= 135:
                return 'Helix'
            elif 135 <= generic_number <= 500:
                return 'C-term'
        return None

    def create_rotamers(self, structure, sd):
        preferred_chains = ast.literal_eval(sd["preferred_chain"])
        print(preferred_chains)

        if not structure.protein or not structure.protein.sequence:
            self.logger.error(f"Protein data or sequence not available for structure {structure}")
            return

        for chain_id in preferred_chains:
            print(chain_id)
            pdb_sequence = self.get_pdb_sequence(structure.pdb_data.pdb, chain_id)
            db_sequence = structure.protein.sequence
            alignment = self.align_sequences(pdb_sequence, db_sequence)
            mapping = self.create_mapping(alignment)
            protein = structure.protein

            pdb = StringIO(structure.pdb_data.pdb)
            parser = PDBParser(QUIET=True)
            pdb_structure = parser.get_structure("structure", pdb)
            model = pdb_structure[0]
            chain = model[chain_id]

            # Get a list of residues in the chain to use their positions
            residues = list(chain.get_residues())

            with transaction.atomic():
                for seq_pos, residue in enumerate(residues, start=1):  # sequence position starts at 1
                    if residue.id[0] == ' ':
                        if seq_pos in mapping:
                            db_res_id = mapping[seq_pos]
                            try:
                                db_residue = Residue.objects.get(protein=protein, sequence_number=db_res_id)
                                residue_type = 'canonical' if db_residue.amino_acid == seq1(residue.resname) else 'mutated'
                                structure_amino_acid = seq1(residue.resname)
                                structure_amino_acid_three_letter = residue.resname

                                Rotamer.objects.create(
                                    residue=db_residue,
                                    structure=structure,
                                    pdbdata=structure.pdb_data,
                                    pdbseq_number=residue.id[1],
                                    sequence_number=db_res_id,
                                    generic_number=db_residue.generic_number,
                                    amino_acid=structure_amino_acid,
                                    amino_acid_three_letter=structure_amino_acid_three_letter,
                                    residue_type=residue_type,
                                    segment=db_residue.segment,
                                    chain=chain_id
                                )
                            except Residue.DoesNotExist:
                                self.logger.error(f"Database residue not found for sequence number {db_res_id} in protein {structure.protein}")
                        else:
                            self.logger.warning(f"Sequence position {seq_pos} not found in mapping.")
                    else:
                        self.logger.warning(f"Skipping non-standard residue: {residue}")

    def main_func(self):
        print(self.parsed_structures)
        Structure.objects.filter().delete()
        Entity.objects.filter().delete()

        entity_types = ["Protein", "DNA", "RNA", "Ligand", "Non-polymer"]

        for entity_name in entity_types:
            EntityType.objects.get_or_create(name=entity_name)

        for pdb_id in self.parsed_structures.pdb_ids:
            sd = self.parsed_structures.structures[pdb_id]
            self.logger.info('Building structure {}'.format(sd['name']))

            try:
                con = Protein.objects.get(accession=sd['protein'])
                print(con)
            except Protein.DoesNotExist:
                print('BIG ERROR Construct {} does not exist, skipping!'.format(sd['name'].lower()))
                self.logger.error('Construct {} does not exist, skipping!'.format(sd['name'].lower()))
                continue

            s = Structure()

            if 'state' not in sd:
                self.logger.warning('State not defined, using default state {}'.format(settings.DEFAULT_PROTEIN_STATE))
                state = settings.DEFAULT_STATE.title()
            else:
                state = sd['state']
            state_slug = slugify(state)
            try:
                ps, created = ProteinState.objects.get_or_create(slug=state_slug, defaults={'name': state})
                if created:
                    self.logger.info('Created protein state {}'.format(ps.name))
            except IntegrityError:
                ps = ProteinState.objects.get(slug=state_slug)
            s.state = ps
            s.author_state = ps

            sd['pdb'] = sd['name'].upper()
            if not os.path.exists(self.pdb_data_dir):
                os.makedirs(self.pdb_data_dir)

            pdb_path = os.sep.join([self.pdb_data_dir, sd['pdb'] + '.pdb'])
            if not os.path.isfile(pdb_path):
                self.logger.info('Fetching PDB file {}'.format(sd['pdb']))
                url = 'http://www.rcsb.org/pdb/files/%s.pdb' % sd['pdb']
                pdbdata_raw = urlopen(url).read().decode('utf-8')
                with open(pdb_path, 'w') as f:
                    f.write(pdbdata_raw)
            else:
                with open(pdb_path, 'r') as pdb_file:
                    pdbdata_raw = pdb_file.read()

            output_directory = os.path.join(settings.DATA_DIR, 'structure_data', 'output_protonated_pdbs')
            os.makedirs(output_directory, exist_ok=True)

            # Check if the protonated PDB already exists
            protonated_pdb_path = os.path.join(settings.DATA_DIR, 'structure_data', 'output_protonated_pdbs', f"{sd['pdb']}_protein.pdb")
            if os.path.exists(protonated_pdb_path):
                self.logger.info(f"Protonated PDB file already exists: {protonated_pdb_path}")
                with open(protonated_pdb_path, 'r') as pdb_file:
                    pdbdata_prot = pdb_file.read()
                    pdbdata, created = PdbData.objects.get_or_create(pdb=pdbdata_prot)
                    s.pdb_data = pdbdata
            else:
                protonate_pdb = handle_protoss_request(sd['pdb'], output_directory)
                print(protonate_pdb)
                if os.path.exists(protonated_pdb_path):
                    with open(protonated_pdb_path, 'r') as pdb_file:
                        pdbdata_prot = pdb_file.read()
                        pdbdata, created = PdbData.objects.get_or_create(pdb=pdbdata_prot)
                        s.pdb_data = pdbdata
                else:
                    self.logger.error('Failed to protonate PDB file {}'.format(pdb_path))

            hetsyn = {}
            hetsyn_reverse = {}
            for line in pdbdata_raw.splitlines():
                if line.startswith('HETSYN'):
                    m = re.match("HETSYN[\s]+([\w]{3})[\s]+(.+)", line)
                    if m:
                        hetsyn[m.group(2).strip()] = m.group(1).upper()
                        hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                if line.startswith('HETNAM'):
                    m = re.match("HETNAM[\s]+([\w]{3})[\s]+(.+)", line)
                    if m:
                        hetsyn[m.group(2).strip()] = m.group(1).upper()
                        hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                if line.startswith('REVDAT   1'):
                    sd['publication_date'] = line[13:22]
                if line.startswith('JRNL        PMID'):
                    sd['pubmed_id'] = line[19:].strip()
                if line.startswith('JRNL        DOI'):
                    sd['doi_id'] = line[19:].strip()

            if len(hetsyn) == 0:
                self.logger.info("PDB file contained NO hetsyn")

            with open(pdb_path, 'r') as header:
                header_dict = parse_pdb_header(header)
            sd['publication_date'] = header_dict['release_date']
            if str(header_dict['resolution']).strip() != 'None':
                sd['resolution'] = str(header_dict['resolution']).strip()
            sd['structure_method'] = header_dict['structure_method']

            if 'structure_method' in sd and sd['structure_method']:
                if sd['structure_method'] == 'unknown':
                    sd['structure_method'] = self.exp_method_dict[sd['method_from_file']]

                structure_type = sd['structure_method'].capitalize()
                structure_type_slug = slugify(sd['structure_method'])

                try:
                    st, created = StructureType.objects.get_or_create(slug=structure_type_slug,
                                                                      defaults={'name': structure_type})
                    if created:
                        self.logger.info('Created structure type {}'.format(st))
                except IntegrityError:
                    st = StructureType.objects.get(slug=structure_type_slug)
                s.structure_type = st
            else:
                self.logger.warning('No structure type specified in PDB file {}'.format(sd['pdb']))

            matched = 0
            if 'ligand' in sd and sd['ligand']:
                if isinstance(sd['ligand'], list):
                    ligands = sd['ligand']
                else:
                    ligands = [sd['ligand']]
                for ligand in ligands:
                    if 'name' in ligand:
                        if ligand['name'].upper() in hetsyn:
                            self.logger.info('Ligand {} matched to PDB records'.format(ligand['name']))
                            matched = 1
                            ligand['name'] = hetsyn[ligand['name'].upper()]
                        elif ligand['name'].upper() in hetsyn_reverse:
                            matched = 1

            if matched == 0 and len(hetsyn) > 0:
                self.logger.info('No ligand names found in HET in structure {}'.format(sd['pdb']))

            if 'pdb' in sd:
                web_resource = WebResource.objects.get(slug='pdb')
                s.pdb_code, created = WebLink.objects.get_or_create(index=sd['pdb'], web_resource=web_resource)
            else:
                self.logger.error('PDB code not specified for structure {}, skipping!'.format(sd['pdb']))
                continue

            if 'preferred_chain' in sd:
                chain_ids = sd['preferred_chain'].split(',')
                s.preferred_chain = ','.join(chain_ids)
            else:
                self.logger.warning('Preferred chain not specified for structure {}'.format(sd['pdb']))
                chain_ids = []

            if 'resolution' in sd:
                s.resolution = float(sd['resolution'])
            else:
                self.logger.warning('Resolution not specified for structure {}'.format(sd['pdb']))

            if 'publication_date' in sd:
                s.publication_date = sd['publication_date']
                if int(s.publication_date[:4]) < 1990:
                    s.publication_date = sd['date_from_file']
                    print('WARNING: publication date for {} is incorrect ({}), switched to ({}) from structures.csv'.format(s, sd['publication_date'], sd['date_from_file']))
            else:
                self.logger.warning('Publication date not specified for structure {}'.format(sd['pdb']))

            try:
                if 'doi_id' in sd:
                    s.publication = Publication.get_or_create_from_doi(sd['doi_id'])
                elif 'pubmed_id' in sd:
                    s.publication = Publication.get_or_create_from_pubmed(sd['pubmed_id'])
            except:
                self.logger.error('Error saving publication'.format(sd['pdb']))

            s.annotated = True
            s.refined = False
            s.stats_text = None
            s.protein = con
            s.save()

            response = get_pdb_entities(s.pdb_code.index)
            if not response:
                return

            pdb_id = response["entry"]["id"]

            entity_ids = response["rcsb_entry_container_identifiers"]["entity_ids"]
            polymer_entity_ids = response["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
            nonpolymer_entity_ids = response.get("nonpolymer_entity_ids", [])

            for entity_id in entity_ids:
                entity_type = "polymer_entity" if entity_id in polymer_entity_ids else "nonpolymer_entity"
                entity_info = fetch_entity_info(pdb_id, entity_id, entity_type)
                if entity_info:
                    print(entity_info)
                    entity_type_instance, _ = EntityType.objects.get_or_create(name=entity_info["type"])
                    description = entity_info["name"]

                    Entity.objects.create(
                        structure=s,
                        entity_type=entity_type_instance,
                        name=description,
                        chain=entity_info["chain"],
                        unp_accession=entity_info["unp_accession"],
                    )

            ResidueFragmentInteraction.objects.filter(structure_partner_pair__structure=s).delete()
            Fragment.objects.filter(structure=s).delete()
            Rotamer.objects.filter(structure=s).delete()

            self.logger.error('Creating rotamers for {}'.format(sd['pdb']))
            try:
                current = time.time()
                self.create_rotamers(structure=s, sd=sd)
                end = time.time()
                diff = round(end - current, 1)
                self.logger.info('Create residues done for {}. {} seconds.'.format(
                    s.protein.protein.entry_name, diff))
            except Exception as msg:
                print(msg)
                print('ERROR WITH ROTAMERS {}'.format(sd['pdb']))
                self.logger.error('Error with rotamers for {}'.format(sd['pdb']))
                self.rotamer_errors.append(s)

            d = {}

            print("NOW CALCULATING INTERACTIONS FOR {}".format(sd['pdb']))

            try:
                current = time.time()

                # Calculate interactions for the structure considering all preferred chains
                interactions_dict = prolif_calculation(structure=s, chain_ids=chain_ids)

                chain_interactions = {}  # To store unique chemokine-partner chain pairs

                for key, value in interactions_dict.items():
                    for residue_pair, interactions in value.items():
                        chemokine_residue_id = f"{residue_pair[0].name} {residue_pair[0].number}"
                        partner_residue_id = f"{residue_pair[1].name} {residue_pair[1].number}"

                        try:
                            # Extract chain IDs from residue_pair
                            chain_id_ligand = residue_pair[0].chain
                            chain_id_protein = residue_pair[1].chain

                            if (chain_id_ligand, chain_id_protein) not in chain_interactions:
                                chain_interactions[(chain_id_ligand, chain_id_protein)] = ChemokinePartnerPair.objects.create(
                                    structure=s,
                                    chemokine_chain=chain_id_ligand,
                                    partner_chain=chain_id_protein,
                                )

                            rotamer_object = Rotamer.objects.get(
                                structure=s,
                                pdbseq_number=residue_pair[0].number,
                                chain=chain_id_ligand
                            )
                        except Rotamer.DoesNotExist:
                            print(f"Rotamer not found for PDB sequence number {residue_pair[0].number} and amino acid {residue_pair[0].name}")
                            continue

                        for interaction_type, details in interactions.items():
                            print(f"Pair Instance, {chemokine_residue_id}, {partner_residue_id}, {interaction_type}")

                            ChemokinePartnerInteraction.objects.create(
                                chemokine_residue=rotamer_object,
                                partner_residue=partner_residue_id,
                                partner_chain=chain_id_protein,
                                interaction_type=interaction_type,
                                chemokine_partner_pair=chain_interactions[(chain_id_ligand, chain_id_protein)]
                            )

                end = time.time()
                diff = round(end - current, 1)
                self.logger.info('Create interaction residue pairs done for {}. {} seconds.'.format(
                    s.protein.protein.entry_name, diff))

            except Exception as msg:
                print(msg)
                print('ERROR WITH CONTACTNETWORK {}'.format(sd['pdb']))
                self.logger.error('Error with contactnetwork for {}'.format(sd['pdb']))
                self.contactnetwork_errors.append(s)

            s.build_check = True
            s.save()
        pass
