from django.shortcuts import render
from django.http import HttpResponse, FileResponse
from django.template.loader import render_to_string
from django.views.generic import View, TemplateView
from django.db.models import Count, Q, Prefetch, TextField
from django.conf import settings
from django.utils.html import format_html
from django.urls import reverse

from common.diagrams_chemokine import *
from structure.utils import *

from structure.models import Structure, PdbData, Chain, Rotamer, Entity, EntityType
from protein.models import Protein
from residue.models import Residue
from interaction.models import ChemokinePartnerInteraction
from .forms import ChainSelectionForm

import requests
from io import StringIO, BytesIO
from Bio.PDB import PDBIO, PDBParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import zipfile
import pandas as pd
import glob
import networkx as nx
import numpy as np
from matplotlib.colors import ListedColormap
import mpld3
import re
import os
import subprocess
from random import randint
from IPython.display import HTML
import plotly.express as px


class StructureBrowser(TemplateView):
    """
    Fetching Structure data for browser
    """
    template_name = "structure/structure_browser.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        structures = Structure.objects.all().select_related('protein')
        context['structures'] = structures
        return context


class StructureDetails(TemplateView):
    template_name = 'structure/structure.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        pdb_id = self.kwargs.get('pdb_id', '')
    
        try:
            structure = Structure.objects.get(pdb_code__index=pdb_id)
            protein = structure.protein
            canonical_sequence = protein.sequence

            # Get all residues for alignment table
            # Retrieve protein residues
            protein_residues = Residue.objects.filter(protein_id=protein).order_by('sequence_number')
            protein_dict = {res.generic_number: res for res in protein_residues if res.generic_number is not None}

            # Retrieve structure rotamers grouped by chain
            structure_rotamers = Rotamer.objects.filter(structure_id=structure).order_by('sequence_number')
            rotamers_by_chain = {}
            for rotamer in structure_rotamers:
                chain_id = rotamer.chain
                if chain_id not in rotamers_by_chain:
                    rotamers_by_chain[chain_id] = {}
                rotamers_by_chain[chain_id][rotamer.generic_number] = rotamer

            # Determine the range of generic numbers
            all_generic_numbers = sorted(set(protein_dict.keys()).union(*(chain_dict.keys() for chain_dict in rotamers_by_chain.values())))

            aligned_sequences = []
            for gn in all_generic_numbers:
                protein_residue = protein_dict.get(gn, None)
                structure_rotamers_for_gn = {chain_id: rotamers_by_chain[chain_id].get(gn, None) for chain_id in rotamers_by_chain}
                aligned_sequences.append((gn, protein_residue, structure_rotamers_for_gn))

            context['aligned_sequences'] = aligned_sequences
            context['generic_numbers'] = all_generic_numbers
            context['rotamers_by_chain'] = rotamers_by_chain
            context['residues_aligned'] = [protein_dict.get(gn) for gn in all_generic_numbers]
            context['chains'] = rotamers_by_chain.keys()

            # Prepare residues grouped by segment
            residues = Residue.objects.filter(protein=protein).order_by('sequence_number')
            residues_by_segment = {}
            for residue in residues:
                segment = residue.segment
                if segment not in residues_by_segment:
                    residues_by_segment[segment] = []
                residues_by_segment[segment].append(residue)

            context['structure'] = structure
            context['protein'] = protein
            context['canonical_sequence'] = canonical_sequence
            context['residues_by_segment'] = residues_by_segment

            # Entities
            entity_types = EntityType.objects.prefetch_related(
                Prefetch(
                    'entity_set',
                    queryset=Entity.objects.filter(structure=structure).order_by('name'),
                    to_attr='filtered_entities'
                )
            )
            context['entity_types'] = entity_types

            # Prepare residues data for the context
            context['residues_data'] = structure_rotamers

            # Get the IDs of these rotamers
            rotamer_ids = structure_rotamers.values_list('id', flat=True)
            # Fetch interactions that are linked to the fetched rotamers
            interactions = ChemokinePartnerInteraction.objects.filter(chemokine_residue_id__in=rotamer_ids)
            context['interactions'] = interactions

            # Manually construct interaction URLs for each chain and count interactions per chain
            interaction_urls = {}
            interactions_count = {}
            for chain_id in rotamers_by_chain.keys():
                interaction_urls[chain_id] = f"/interaction/{structure.pdb_code.index}/{chain_id}/"
                interactions_count[chain_id] = interactions.filter(chemokine_residue__chain=chain_id).exclude(
                partner_chain=chain_id
            ).count()

            context['interaction_urls'] = interaction_urls
            context['interactions_count'] = interactions_count

        except Structure.DoesNotExist:
            return render(request, 'error.html')
        
        return context



        
    def get(self, request, *args, **kwargs):
        context = self.get_context_data(**kwargs)
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        form = ChainSelectionForm(request.POST)
        context = self.get_context_data(**kwargs)  # Reuse get_context_data to set up initial context
        pdb_id = self.kwargs.get('pdb_id', '')

        if form.is_valid():
            selected_chain = '.' + form.cleaned_data['chain']
            file_path = os.path.join(settings.DATA_DIR, 'prepared_pdbs', f"{pdb_id}_prepared", f"{pdb_id}_all_interactions.csv")

            if os.path.exists(file_path):
                df = pd.read_csv(file_path, header=None)
                first_filter = df.columns[df.iloc[0].str.endswith(selected_chain)]
                first_filtered_df = df[first_filter]
                second_filter = first_filtered_df.columns[~first_filtered_df.iloc[1].str.endswith(selected_chain)]
                final_filtered_df = first_filtered_df[second_filter]
                interactions = final_filtered_df.to_html(classes=["table", "table-bordered", "table-striped"], index=False)
            else:
                interactions = "File not found."

            context['interactions'] = interactions  # Add filtered interactions to the context

        context['form'] = form  # Include the form to maintain form state or display errors
        return render(request, self.template_name, context)
    
    
    
    
class StructureInteractions(TemplateView):
    """
    Show chemokine interactions from selected structure and chain
    """
    template_name = 'structure/structure_interactions.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        
        pdb_id = self.kwargs.get('pdb_id', '')
        chain_id = self.kwargs.get('chain', '')
        
        try:
            structure = Structure.objects.get(pdb_id=pdb_id)
            context['structure'] = structure
            context['protein'] = structure.protein
            context['partner'] = structure.partner.all()
            #context['entities'] = structure.chain_set.all()
            context['chain_id'] = chain_id
            

            try:
                chain = Chain.objects.get(structure=structure, chain=chain_id)
            except Chain.DoesNotExist:
                residues = []   


            # Prepare residues data for the context
            rotamers = Rotamer.objects.filter(structure=structure).select_related('residue')
            residues_data = [{
                        'residue_id': rotamer.residue.sequence_number,
                        'amino_acid': rotamer.residue.amino_acid_three_letter,
                        'chain_id': rotamer.residue.chain.chain_id if hasattr(rotamer.residue, 'chain') else None
                    } for rotamer in rotamers]
            context['residues_data'] = residues_data
            
            
            ## If focusing on a specific chemokine chain, filter interactions accordingly; otherwise, get all interactions for the structure
            #if chain_id:
            #    chemokine_interactions = ChemokineInteraction.objects.filter(structure=structure, chemokine_chain__chain=chain_id)
            #else:
            #    chemokine_interactions = ChemokineInteraction.objects.filter(structure=structure)
#
            #context['chemokine_interactions'] = chemokine_interactions
            
            
            file_path = os.path.join(settings.DATA_DIR, 'prepared_pdbs', f"{pdb_id}_prepared", f"{pdb_id}_all_interactions.csv")
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, header=None).transpose()

                # Adding column headers for readability and ease of filtering
                df.columns = ['ChemokineResidue', 'PartnerResidue', 'InteractionType', 'InteractionPresence']

                # Filter: Select rows where 'ChemokineResidue' ends with the selected chain ID
                df_filtered_by_chain = df[df['ChemokineResidue'].str.endswith('.' + chain_id)]

                # Further filter: Exclude rows where 'PartnerResidue' also ends with the selected chain ID
                df_filtered_excluding_partner = df_filtered_by_chain[~df_filtered_by_chain['PartnerResidue'].str.endswith('.' + chain_id)]

                # Additional filter for "VdWContact" in 'InteractionType'
                final_filtered_df = df_filtered_excluding_partner[~df_filtered_excluding_partner['InteractionType'].str.endswith("VdWContact")]
                
                interaction_data = []

                for index, row in final_filtered_df.iterrows():
                    interaction_row = {
                        'ChemokineResidue': row[0],
                        'PartnerResidue': row[1],
                        'InteractionType': row[2],
                        'InteractionPresent': row[3]
                    }
                    interaction_data.append(interaction_row)
    
            else:
                interaction_data = None
    
            context['interactions'] = interaction_data
                        


        except Structure.DoesNotExist:
            context['error'] = "Structure does not exist."

        return context