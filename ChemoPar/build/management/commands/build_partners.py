from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import Protein
from structure.models import Structure
from partner.models import Partner, Partner_PDB, PDBData

import os
import pandas as pd
import re
import urllib.request
import csv
import pprint
import urllib
import json
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select
from datetime import datetime
from Bio.PDB import Polypeptide, PDBParser, Superimposer
from Bio.PDB.PDBIO import PDBIO
from copy import deepcopy
import urllib.error
from Bio.PDB import PDBParser, PDBIO, Select
from io import StringIO

class ChainSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.get_id() in self.chain_ids

class Command(BaseCommand):
    help = 'Builds PDB data from partner models.'

    pdbs_path = os.sep.join([settings.DATA_DIR, 'raw_pdbs'])

    def add_arguments(self, parser):
        parser.add_argument("--debug", default=False, action="store_true", help="Debug mode")
        parser.add_argument("--purge", default=False, action="store_true", help="Purge data")

    def handle(self, *args, **options):
        if options['purge']:
            Partner.objects.all().delete()

        # Read partner table and store in database
        self.store_partner_list()
        
        # Store PDB chains of partners
        #self.store_partner_pdbdata()
        
        self.stdout.write(self.style.SUCCESS('Successfully built PDB data.'))


    def store_partner_list(self):
        file_path = os.path.join(settings.DATA_DIR, 'partner_list.xlsx')
        data = pd.read_excel(file_path)

        for _, row in data.iterrows():
            # Create or get the structure
            pdb_id = row['PDB']
            structure, _ = Structure.objects.get_or_create(pdb_id=pdb_id)

            # If partner information is present, process it
            if not pd.isna(row['Partner']):
                self.store_row_partners(structure, row)

        self.stdout.write(self.style.SUCCESS('Data imported successfully.'))

    def store_row_partners(self, structure, row):
        partner_name = row['Partner']
        partner_type = row['Partner_type']

        try:
            # Get or create the partner
            partner, created = Partner.objects.get_or_create(
                name=partner_name, 
                defaults={'type': partner_type}
            )

            # Update partner's type if it already existed
            if not created and partner.type != partner_type:
                partner.type = partner_type
                partner.save()

            # Add structure to partner's pdb list
            partner.structures.add(structure)

        except IntegrityError as e:
            self.stdout.write(self.style.ERROR(f"Failed to save structure with PDB ID {structure.pdb_id} due to integrity error: {e}"))


    def store_partner_pdbdata(self):
        # Set up directories
        directory_path = os.sep.join([settings.DATA_DIR, 'raw_pdbs'])
        output_directory = os.sep.join([settings.DATA_DIR, 'partner_pdbs'])
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        # Initialize PDB parser
        parser = PDB.PDBParser(QUIET=True)

        # Process each partner entry
        for partner in Partner.objects.all():
            pdb_id = partner.pdb_id
            chain_ids = partner.chain

            # Construct the filename for the corresponding PDB file
            pdb_filename = f"{pdb_id}.pdb"

            try:
                # Parse the PDB file
                structure = parser.get_structure('structure', os.path.join(directory_path, pdb_filename))

                for model in structure:
                    for chain in model:
                        if chain.id in chain_ids:
                            # Create output filename and save the chain
                            output_filename = f"{pdb_id}_{chain.id}.pdb"
                            io = PDB.PDBIO()
                            io.set_structure(chain)
                            io.save(os.path.join(output_directory, output_filename))
            except Exception as e:
                # Handle exceptions
                print(f"Error processing {pdb_id}: {e}")










    #def store_partner_pdbdata(self):
    #    directory_path = os.sep.join([settings.DATA_DIR, 'raw_pdbs'])
    #    output_directory = os.sep.join([settings.DATA_DIR, 'partner_pdbs'])
    #    if not os.path.exists(output_directory):
    #        os.makedirs(output_directory)
#
    #    parser = PDB.PDBParser(QUIET=True)
#
    #    for filename in os.listdir(directory_path):
    #        if filename.endswith(".pdb"):
    #            try:
    #                structure = parser.get_structure('structure', os.path.join(directory_path, filename))
    #                pdb_id = filename[:4]
    #                partners = Partner.objects.filter(pdb_id=pdb_id)
    #                chain_ids = [partner.chain for partner in partners if partner.chain]
#
    #                for model in structure:
    #                    for chain in model:
    #                        if chain.id in chain_ids:
    #                            output_filename = f"{filename[:-4]}_{chain.id}.pdb"
    #                            io = PDB.PDBIO()
    #                            io.set_structure(chain)
    #                            io.save(os.path.join(output_directory, output_filename))
    #            except Exception as e:
    #                self.stdout.write(self.style.ERROR(f"Error processing {filename}: {e}"))
#
    #def build_pdb_data(self):
    #        for filename in os.listdir(self.pdbs_path):
    #            if filename.endswith(".pdb"):
    #                model_filename = os.path.join(self.pdbs_path, filename)
    #                pdbdata = self.build_pdbdata(model_filename)
    #                name = filename.replace(".pdb", "")
    #                Partner_PDB.objects.get_or_create(name=name, pdbdata=pdbdata)
#