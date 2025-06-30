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

from common.models import ResiduePosition
from structure.models import Structure, PdbData, Rotamer, Entity, EntityType, ChemokineBindingPartner
from protein.models import Protein
from residue.models import Residue
from interaction.models import ChemokinePartnerInteraction
from .forms import ChainSelectionForm

import requests
from io import StringIO, BytesIO
from Bio.PDB import PDBIO, PDBParser
import zipfile
import pandas as pd
import glob
import networkx as nx
import numpy as np
import re
import os
import subprocess
from random import randint
import plotly.express as px
from io import BytesIO


from django.shortcuts import render, get_object_or_404
from django.http import HttpResponse
from structure.models import Structure, Rotamer, Residue
from protein.models import Protein

def download_alignment(request, pdb_id):
    structure = get_object_or_404(Structure, pdb_code__index=pdb_id)
    protein = structure.protein
    canonical_sequence = protein.sequence

    protein_residues = Residue.objects.filter(protein_id=protein).order_by('sequence_number')
    protein_dict = {res.generic_number: res for res in protein_residues if res.generic_number is not None}

    structure_rotamers = Rotamer.objects.filter(structure_id=structure).order_by('sequence_number')
    rotamers_by_chain = {}
    for rotamer in structure_rotamers:
        chain_id = rotamer.chain
        if chain_id not in rotamers_by_chain:
            rotamers_by_chain[chain_id] = {}
        rotamers_by_chain[chain_id][rotamer.generic_number] = rotamer

    all_generic_numbers = sorted(set(protein_dict.keys()).union(*(chain_dict.keys() for chain_dict in rotamers_by_chain.values())))

    aligned_sequences = []
    for gn in all_generic_numbers:
        protein_residue = protein_dict.get(gn, None)
        structure_rotamers_for_gn = {chain_id: rotamers_by_chain[chain_id].get(gn, None) for chain_id in rotamers_by_chain}
        aligned_sequences.append((gn, protein_residue, structure_rotamers_for_gn))

    aln_content = 'CLUSTAL\n\n'

    uniprot_sequence = ''.join(protein_residue.amino_acid if protein_residue else '-' for gn, protein_residue, _ in aligned_sequences)
    aln_content += f'Unipr/ {uniprot_sequence}\n'

    for chain_id in rotamers_by_chain.keys():
        structure_sequence = ''.join(
            rotamers_by_chain[chain_id].get(gn, '-').amino_acid if gn in rotamers_by_chain[chain_id] and rotamers_by_chain[chain_id][gn] else '-'
            for gn, _, _ in aligned_sequences
        )
        aln_content += f'Chain_{chain_id}/ {structure_sequence}\n'

    response = HttpResponse(aln_content, content_type='text/plain')
    response['Content-Disposition'] = f'attachment; filename="{pdb_id}_alignment.aln"'
    return response


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


from common.models import ResiduePosition


class StructureDetails(TemplateView):
    template_name = 'structure/structure.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        structure_id = self.kwargs.get('structure_id')

        try:
            # 1) Retrieve this Structure and its Protein
            structure = Structure.objects.get(id=structure_id)
            protein = structure.protein
            canonical_sequence = protein.sequence

            # 2) Fetch all ResiduePosition objects (sorted by integer `position` field)
            all_positions = ResiduePosition.objects.all().order_by('position')
            pos_lookup = {pos_obj.ccn_number: pos_obj for pos_obj in all_positions}
            context['all_positions'] = all_positions

            # 3) Fetch all Residues for this protein
            protein_residues = Residue.objects.filter(
                protein=protein
            ).order_by('sequence_number')
            print(protein_residues)

            # Build a lookup: position (int) → Residue
            protein_dict = {}
            for res in protein_residues:
                if res.ccn_number and res.ccn_number in pos_lookup:
                    pos_obj = pos_lookup[res.ccn_number]
                    protein_dict[pos_obj.position] = res

            # 4) Fetch all Rotamers for this Structure
            structure_rotamers = Rotamer.objects.filter(
                structure=structure
            ).order_by('sequence_number')
            print(structure_rotamers)

            # Group rotamers by chain_id: rotamers_by_chain = { 'A': { pos_id: rotamer, … }, … }
            rotamers_by_chain = {}
            for rot in structure_rotamers:
                if hasattr(rot, 'ccn_number') and rot.ccn_number and rot.ccn_number in pos_lookup:
                    pos_obj = pos_lookup[rot.ccn_number]
                    pos_id = pos_obj.position
                    chain_id = rot.chain
                    if chain_id not in rotamers_by_chain:
                        rotamers_by_chain[chain_id] = {}
                    rotamers_by_chain[chain_id][pos_id] = rot
                    
            # 5) Build aligned_sequences as a list of triples:
            #      ( position_obj, protein_residue_or_None, { chain_id: Rotamer or None } )
            aligned_sequences = []
            for pos_obj in all_positions:
                pos_id = pos_obj.position
                protein_res = protein_dict.get(pos_id, None)
                rot_for_pos = {
                    chain_id: rotamers_by_chain.get(chain_id, {}).get(pos_id)
                    for chain_id in rotamers_by_chain
                }
                aligned_sequences.append((pos_obj, protein_res, rot_for_pos))

            context['aligned_sequences'] = aligned_sequences

            # 6) Also pass the list of chain IDs that appear in this structure
            context['chains'] = list(rotamers_by_chain.keys())

            # 7) Prepare “residues_by_segment” if you need by segment
            residues_by_segment = {}
            for residue in protein_residues:
                seg = getattr(residue, 'segment', None)
                if seg not in residues_by_segment:
                    residues_by_segment[seg] = []
                residues_by_segment[seg].append(residue)
            context['residues_by_segment'] = residues_by_segment

            # 8) Pass Structure / Protein objects into context
            context['structure'] = structure
            context['protein'] = protein
            context['canonical_sequence'] = canonical_sequence

            # 9) Entities (same as before)
            entity_types = EntityType.objects.prefetch_related(
                Prefetch(
                    'entity_set',
                    queryset=Entity.objects.filter(structure=structure).order_by('name'),
                    to_attr='filtered_entities'
                )
            )
            context['entity_types'] = entity_types

            # 10) Binding partners & interaction counts (same as before)
            binding_partners = ChemokineBindingPartner.objects.filter(
                structure=structure
            ).annotate(interaction_count=Count('chemokinepartnerinteraction'))

            interaction_urls = {}
            interactions_count = {}
            for partner in binding_partners:
                if partner.interaction_count > 0:
                    interaction_urls[partner] = f"/interaction/{structure.id}/{partner.id}/"
                    interactions_count[partner] = partner.interaction_count

            if not interactions_count:
                interactions_count = None

            context['binding_partners'] = binding_partners
            context['interaction_urls'] = interaction_urls
            context['interactions_count'] = interactions_count

        except Structure.DoesNotExist:
            return render(self.request, 'error.html')

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
            

            #try:
            #    chain = Chain.objects.get(structure=structure, chain=chain_id)
            #except Chain.DoesNotExist:
            #    residues = []   


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
    

#==============================================================================
class PDBDownload(View):
    """
    Serve the PDB file for a specific structure.
    """

    def get(self, request, structure_id, *args, **kwargs):
        # Fetch the structure by ID
        try:
            structure = Structure.objects.get(id=structure_id)
        except Structure.DoesNotExist:
            raise Http404("Structure not found.")

        # Check if the structure has associated PDB data
        if not structure.pdb_data or not structure.pdb_data.pdb:
            return HttpResponse("PDB data not available for this structure.", status=404)

        # Prepare PDB file content
        pdb_content = structure.pdb_data.pdb
        pdb_filename = f"{structure.pdb_code.index}.pdb"

        # Create HTTP response with PDB content
        response = HttpResponse(pdb_content, content_type="chemical/x-pdb")
        response['Content-Disposition'] = f'attachment; filename="{pdb_filename}"'

        return response