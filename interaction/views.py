import os
import re
import json
import glob
import zipfile
import subprocess
import base64
import io
from io import StringIO, BytesIO
from random import randint
from collections import defaultdict
import csv

import requests
from umap import UMAP
import numpy as np
import pandas as pd
import networkx as nx
import plotly.express as px
import plotly.graph_objects as go
import xlsxwriter
import seaborn as sns
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend for Django
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from Bio.PDB import PDBIO, PDBParser
from rdkit.DataStructs import ExplicitBitVect, TanimotoSimilarity

# Django imports
from django.shortcuts import render, get_object_or_404
from django.shortcuts import render, redirect
from django.http import HttpResponse, FileResponse, JsonResponse
from django.template.loader import render_to_string
from django.views.generic import View, TemplateView
from django.conf import settings
from django.utils.html import format_html
from django.views.decorators.http import require_GET
from django.views.decorators.csrf import csrf_exempt
from django.utils.decorators import method_decorator

# Local app imports
from .forms import PDBSelectionForm
from common.diagrams_chemokine import *
from structure.utils import *
from structure.models import Structure, PdbData, Rotamer, Entity, ChemokineBindingPartner
from protein.models import Protein
from residue.models import Residue
from interaction.models import (
    ChemokinePartnerInteraction,
    ChemokinePartnerIFP,
    LigandNetworkHTML,
)


from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

def bitstring_to_bitvect(bitstring):
    """
    Converts a bitstring to an ExplicitBitVect.
    """
    bv = ExplicitBitVect(len(bitstring))
    for i, bit in enumerate(bitstring):
        if bit == '1':
            bv.SetBit(i)
    return bv

def calculate_tanimoto_similarity(ifp1, ifp2):
    """
    Calculates the Tanimoto similarity between two bitstrings.
    """
    if len(ifp1) != len(ifp2):
        raise ValueError("Bit strings must be of the same length")

    bv1 = bitstring_to_bitvect(ifp1)
    bv2 = bitstring_to_bitvect(ifp2)
    return TanimotoSimilarity(bv1, bv2)

def pad_bitstrings(ifp_strings):
    """
    Pads all bitstrings to the length of the longest string with trailing zeros.
    """
    max_length = max(len(ifp) for ifp in ifp_strings)
    return [ifp.ljust(max_length, '0') for ifp in ifp_strings]



#==============================================================================
class PDBDownload(View):
    """
    Serve the PDB file for a specific structure.
    """

    def get(self, request, binding_partner_id, *args, **kwargs):
        # Fetch the structure by ID
        try:
            binding_partner = ChemokineBindingPartner.objects.get(id=binding_partner_id)
        except ChemokineBindingPartner.DoesNotExist:
            raise Http404("ChemokineBindingPartner not found.")

        # Check if the structure has associated PDB data
        if not binding_partner.pdb_data or not binding_partner.pdb_data.pdb:
            return HttpResponse("PDB data not available for this ChemokineBindingPartner.", status=404)

        # Prepare PDB file content
        pdb_content = binding_partner.pdb_data.pdb
        pdb_filename = f"{binding_partner.structure.pdb_code.index}.pdb"

        # Create HTTP response with PDB content
        response = HttpResponse(pdb_content, content_type="chemical/x-pdb")
        response['Content-Disposition'] = f'attachment; filename="{pdb_filename}"'

        return response




# =============================================================================
# CSV export
# =============================================================================

def csv_export(request, pdb_id, chain_id, binding_partner):
    try:
        structure = Structure.objects.get(pdb_code__index=pdb_id)
    except Structure.DoesNotExist:
        return render(
            request,
            'interaction/error.html',
            {'message': f"Structure with PDB ID {pdb_id} does not exist."}
        )

    # Prepare rotamer data (not needed for CSV but kept for compatibility)
    rotamers = Rotamer.objects.filter(structure=structure, chain=chain_id).select_related('residue')
    
    # Fetch interactions for these rotamers
    interactions = ChemokinePartnerInteraction.objects.filter(
        structure=structure, 
        chemokine_binding_partner=binding_partner
    ).select_related('chemokine_residue__residue')

    data = []
    for interaction in interactions:
        residue = interaction.chemokine_residue.residue
        row = {
            'Chemokine Residue': f"{residue.amino_acid_three_letter} {residue.sequence_number}",
            'CCN Number': getattr(residue, 'ccn_number', 'N/A'),
            'Segment': residue.segment,
            'Partner Residue': interaction.partner_residue,
            'Partner': interaction.partner_chain,
            'Interaction Type': interaction.interaction_type,
        }
        data.append(row)

    headers = ['Chemokine Residue', 'CCN Number', 'Segment', 'Partner Residue', 'Partner', 'Interaction Type']

    # Generate CSV file
    output = StringIO()
    writer = csv.DictWriter(output, fieldnames=headers)
    writer.writeheader()
    for row in data:
        writer.writerow(row)

    response = HttpResponse(
        output.getvalue(),
        content_type='text/csv'
    )
    response['Content-Disposition'] = f'attachment; filename=Interaction_data_{pdb_id}_{chain_id}.csv'
    return response

# =============================================================================
# PDB Structure Selection and Alignment Views
# =============================================================================

def SelectStructure(request):
    if request.method == 'POST':
        form = PDBSelectionForm(request.POST)
        if form.is_valid():
            selected_pdb_ids = form.cleaned_data['pdb_ids']
            # Redirect to the alignment results view
            return display_alignment_results(request, selected_pdb_ids)
    else:
        form = PDBSelectionForm()
    return render(request, 'interaction/select_structure.html', {'form': form})


def display_alignment_results(request, pdb_ids):
    alignment_result = align_sequences(pdb_ids)
    context = {
        'alignment_result': alignment_result,
        'pdb_ids': pdb_ids,
    }
    return render(request, 'pdbapp/display_alignment_results.html', context)


class ViewAlignment(TemplateView):
    template_name = 'interaction/view_alignment.html'

    def post(self, request, *args, **kwargs):
        pdb_indices = request.POST.get('pdb_indices', '')
        pdb_ids = [p.strip() for p in pdb_indices.split(',') if p.strip()]
        if not pdb_ids:
            return render(request, 'error.html', {'message': 'No PDB indices provided.'})

        alignment_result = align_sequences(pdb_ids)
        context = {
            'alignment_result': alignment_result,
            'pdb_ids': pdb_ids,
        }
        return self.render_to_response(context)

    def get(self, request, *args, **kwargs):
        ids = kwargs.get('ids', '')
        pdb_ids = [p.strip() for p in ids.split(',') if p.strip()]
        if not pdb_ids:
            return render(request, 'error.html', {'message': 'No PDB IDs provided.'})

        alignment_result = align_sequences(pdb_ids)
        context = {
            'alignment_result': alignment_result,
            'pdb_ids': pdb_ids,
        }
        return self.render_to_response(context)


def align_sequences(pdb_ids):
    # Placeholder implementation
    return "Aligned sequences for PDB IDs: " + ", ".join(pdb_ids)

# =============================================================================
# Interaction Details View
# =============================================================================

class InteractionDetails(TemplateView):
    template_name = 'interaction/interaction_details_staticresidues.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        structure_id = self.kwargs.get('structure_id', '')
        chain_id = self.kwargs.get('chain_id', '')
        context['chain'] = chain_id
        
        try:
            structure = Structure.objects.get(id=structure_id)
            context['structure'] = structure
            context['protein'] = structure.protein

            # Ligand network HTML entries
            context['ligand_network_html_entries'] = LigandNetworkHTML.objects.filter(structure=structure)

            # Rotamers and associated residues
            rotamers = Rotamer.objects.filter(structure=structure, chain=chain_id).select_related('residue')
            context['residues_data'] = rotamers
            rotamer_ids = rotamers.values_list('id', flat=True)

            interactions = ChemokinePartnerInteraction.objects.filter(chemokine_residue_id__in=rotamer_ids)
            context['interactions'] = interactions

            # Unique sequence numbers for display
            sequence_numbers = [
                f"{interaction.chemokine_residue.pdbseq_number}:{chain_id}"
                for interaction in interactions
                if interaction.chemokine_residue and interaction.chemokine_residue.residue
            ]
            context['display_res'] = list(set(sequence_numbers))

            # Organize interacting residues by chain
            interacting_residues_by_chain = defaultdict(list)
            residues_lookup = {}
            for interaction in interactions:
                residue_number = interaction.partner_residue.split()[1]
                key = f"{residue_number}:{interaction.partner_chain}"
                interacting_residues_by_chain[interaction.partner_chain].append(key)
                residues_lookup[int(residue_number)] = interaction.partner_residue

            context['interacting_residues_by_chain'] = dict(interacting_residues_by_chain)
            context['residues_lookup'] = json.dumps(residues_lookup)

            # IFP entries and their interactions
            ifp_entries = ChemokinePartnerIFP.objects.filter(structure=structure)
            ifp_data_by_structure = defaultdict(lambda: {'entity_names': '', 'chain_id': '', 'ifp_data': [], 'interactions': []})

            interactions = ChemokinePartnerInteraction.objects.filter(structure=structure).select_related('chemokine_residue__residue')
            for ifp_entry in ifp_entries:
                sid = ifp_entry.structure.id  
                if not ifp_data_by_structure[sid]['entity_names']:
                    ifp_data_by_structure[sid]['entity_names'] = ifp_entry.structure.pdb_code.index if ifp_entry.structure.pdb_code else "Unknown"
                    ifp_data_by_structure[sid]['chain_id'] = ifp_entry.structure.chain_id
                ifp_data_by_structure[sid]['ifp_data'].append({'ifp_string': ifp_entry.ifp_string})

            for interaction in interactions:
                sid = interaction.structure.id  
                rotamer = interaction.chemokine_residue
                if not ifp_data_by_structure[sid]['entity_names']:
                    ifp_data_by_structure[sid]['entity_names'] = interaction.structure.pdb_code.index if interaction.structure.pdb_code else "Unknown"
                    ifp_data_by_structure[sid]['chain_id'] = interaction.structure.chain_id
                generic_number = rotamer.generic_number if rotamer else None
                if generic_number:
                    ifp_data_by_structure[sid]['interactions'].append({
                        'generic_number': generic_number,
                        'chemokine_residue': f"{rotamer.residue.amino_acid_three_letter} {generic_number}" if rotamer else "Unknown",
                        'partner_residue': interaction.partner_residue,
                        'partner_chain': interaction.partner_chain,
                        'interaction_type': interaction.interaction_type,
                    })

            context['ifp_data_by_structure'] = dict(ifp_data_by_structure)

            # Process IFP data list
            possible_interactions = [
                "Hydrophobic", "Pistacking", "HBDonor", "HBAcceptor",
                "Anionic", "Cationic", "CationPi", "PiCation", "VdWContact"
            ]
            ifp_data_list = []
            for ifp_entry in ifp_entries:
                ifp_string = ifp_entry.ifp_string
                ifp_data = []
                for i, generic_number in enumerate(range(30, 151)):
                    interaction_bits = ifp_string[i*9:(i+1)*9]
                    interactions_list = [
                        possible_interactions[j] for j, bit in enumerate(interaction_bits) if bit == '1'
                    ]
                    ifp_data.append({
                        'generic_number': generic_number,
                        'ifp_string': interaction_bits,
                        'interactions': ", ".join(interactions_list) if interactions_list else "None"
                    })
                ifp_data_list.append(ifp_data)

            context['ifp_data_list'] = ifp_data_list

        except Structure.DoesNotExist:
            return render(self.request, 'home/error.html')

        return context


class ChemokineBindingPartnerDetailView(TemplateView):
    template_name = 'interaction/chemokine_binding_partner_detail.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        binding_partner_id = self.kwargs.get('pk')

        try:
            binding_partner = ChemokineBindingPartner.objects.get(id=binding_partner_id)
        except ChemokineBindingPartner.DoesNotExist:
            return render(self.request, 'home/error.html')

        context['binding_partner'] = binding_partner
        structure = binding_partner.structure
        context['structure'] = structure
        context['protein'] = structure.protein

        # ───────────────────────
        # Interactions for this partner
        # ───────────────────────
        interactions = ChemokinePartnerInteraction.objects.filter(
            chemokine_binding_partner=binding_partner
        ).select_related('chemokine_residue__residue')
        context['interactions'] = interactions
        print(interactions)

        # Build display_res (used later for coloring in the snake plot)
        display_residues = [
            f"{interaction.chemokine_residue.pdbseq_number}:{binding_partner.chemokine_chain}"
            for interaction in interactions
            if interaction.chemokine_residue and interaction.chemokine_residue.residue
        ]
        context['display_res'] = list(set(display_residues))

        # Build interacting_residues_by_chain & residues_lookup (for NGL)
        interacting_residues_by_chain = defaultdict(list)
        residues_lookup = {}
        for interaction in interactions:
            try:
                # partner_residue is like "ALA 123"
                residue_number = interaction.partner_residue.split()[1]
                key = f"{residue_number}:{interaction.partner_chain}"
                interacting_residues_by_chain[interaction.partner_chain].append(key)
                residues_lookup[int(residue_number)] = interaction.partner_residue
            except (IndexError, ValueError):
                continue
        context['interacting_residues_by_chain'] = dict(interacting_residues_by_chain)
        context['residues_lookup'] = json.dumps(residues_lookup)

        # ───────────────────────
        # Fetch IFP entries for this partner
        # ───────────────────────
        ifp_entries = ChemokinePartnerIFP.objects.filter(binding_pair=binding_partner)

        # — 1) Add the "first" IFP record into context as ifp_record —
        try:
            ifp_record = ifp_entries.first()  # or .get(...) if 1:1 guaranteed
        except ChemokinePartnerIFP.DoesNotExist:
            ifp_record = None
        context['ifp_record'] = ifp_record

        # — 2) The rest of your existing grouping logic can stay the same —
        ifp_data_by_structure = defaultdict(lambda: {
            'entity_names': '',
            'chain_id': '',
            'ifp_data': [],
            'interactions': []
        })
        interactions_all = ChemokinePartnerInteraction.objects.filter(
            chemokine_binding_partner=binding_partner
        ).select_related('chemokine_residue__residue')

        for ifp_entry in ifp_entries:
            sid = ifp_entry.structure.id
            if not ifp_data_by_structure[sid]['entity_names']:
                if ifp_entry.structure.pdb_code:
                    ifp_data_by_structure[sid]['entity_names'] = ifp_entry.structure.pdb_code.index
                else:
                    ifp_data_by_structure[sid]['entity_names'] = "Unknown"
                ifp_data_by_structure[sid]['chain_id'] = ifp_entry.structure.chain_id
            ifp_data_by_structure[sid]['ifp_data'].append({
                'ifp_string': ifp_entry.ifp_string
            })

        for interaction in interactions_all:
            sid = interaction.structure.id
            rotamer = interaction.chemokine_residue
            if not ifp_data_by_structure[sid]['entity_names']:
                if interaction.structure.pdb_code:
                    ifp_data_by_structure[sid]['entity_names'] = interaction.structure.pdb_code.index
                else:
                    ifp_data_by_structure[sid]['entity_names'] = "Unknown"
                ifp_data_by_structure[sid]['chain_id'] = interaction.structure.chain_id

            generic_number = rotamer.generic_number if rotamer else None
            if generic_number:
                ifp_data_by_structure[sid]['interactions'].append({
                    'generic_number': generic_number,
                    'chemokine_residue': f"{rotamer.residue.amino_acid_three_letter} {generic_number}"
                                          if rotamer and rotamer.residue else "Unknown",
                    'partner_residue': interaction.partner_residue,
                    'partner_chain': interaction.partner_chain,
                    'interaction_type': interaction.interaction_type,
                })

        context['ifp_data_by_structure'] = dict(ifp_data_by_structure)

        # Build ifp_data_list (your bit‐decomposition for 30–150)
        possible_interactions = [
            "Hydrophobic", "Pistacking", "HBDonor", "HBAcceptor",
            "Anionic", "Cationic", "CationPi", "PiCation", "VdWContact"
        ]
        ifp_data_list = []
        for ifp_entry in ifp_entries:
            ifp_string = ifp_entry.ifp_string or ""
            row_data = []
            for i, generic_number in enumerate(range(30, 151)):
                interaction_bits = ifp_string[i*9:(i+1)*9]
                interactions_list = [
                    possible_interactions[j]
                    for j, bit in enumerate(interaction_bits) if bit == '1'
                ]
                row_data.append({
                    'generic_number': generic_number,
                    'ifp_string': interaction_bits,
                    'interactions': ", ".join(interactions_list) if interactions_list else "None"
                })
            ifp_data_list.append(row_data)
        context['ifp_data_list'] = ifp_data_list

        # ───────────────────────
        # Interactions for this partner (raw, flat)
        # ───────────────────────
        interactions = ChemokinePartnerInteraction.objects.filter(
            chemokine_binding_partner=binding_partner
        ).select_related('chemokine_residue__residue')
        context['interactions'] = interactions

        # Group by ccn_number for the snake plot and tooltips
        grouped = defaultdict(lambda: {
            'ccn_number': None,
            'chemokine_residue': None,
            'interaction_types': set(),
        })

        for interaction in interactions:
            residue = interaction.chemokine_residue.residue if interaction.chemokine_residue else None
            if not residue:
                continue
            ccn = residue.ccn_number
            if not ccn:
                continue
            grouped[ccn]['ccn_number'] = ccn
            grouped[ccn]['chemokine_residue'] = residue.amino_acid_three_letter
            grouped[ccn]['interaction_types'].add(interaction.interaction_type)

        # Convert sets to sorted lists for JSON serialization and consistency
        grouped_interactions = []
        for g in grouped.values():
            grouped_interactions.append({
                'ccn_number': g['ccn_number'],
                'chemokine_residue': g['chemokine_residue'],
                'interaction_types': sorted(list(g['interaction_types'])),
            })

        context['interactions_for_partner'] = grouped_interactions

        print(context['interactions_for_partner'])

        return context


# =============================================================================
# AJAX View for Interaction Data
# =============================================================================

@require_GET
def ajax(request):
    structure_id = request.GET.get('structure')
    chemokine_partner_pair_id = request.GET.get('chemokine_partner_pair_id')

    try:
        structure = Structure.objects.get(pdb_code__index=structure_id)
    except Structure.DoesNotExist:
        return JsonResponse({'error': 'Structure not found'}, status=404)

    try:
        chemokine_partner_pair = ChemokinePartnerPair.objects.get(
            structure=structure,
            id=chemokine_partner_pair_id
        )
    except ChemokinePartnerPair.DoesNotExist:
        return JsonResponse({'error': 'ChemokinePartnerPair not found'}, status=404)

    interactions = ChemokinePartnerInteraction.objects.filter(
        chemokine_partner_pair=chemokine_partner_pair
    ).select_related('chemokine_residue', 'chemokine_residue__residue')

    bundled_interactions = defaultdict(lambda: {
        'sequence_number': '',
        'generic_number': '',
        'interaction_types': set(),
        'amino_acid': ''
    })

    for interaction in interactions:
        residue = interaction.chemokine_residue.residue
        key = residue.sequence_number
        if key:
            bundled_interactions[key]['sequence_number'] = key
            bundled_interactions[key]['generic_number'] = residue.generic_number if residue.generic_number else ''
            bundled_interactions[key]['interaction_types'].add(interaction.interaction_type)
            bundled_interactions[key]['amino_acid'] = residue.amino_acid

    bundled_interactions = [
        {
            'sequence_number': key,
            'generic_number': value['generic_number'],
            'interaction_types': list(value['interaction_types']),
            'amino_acid': value['amino_acid']
        }
        for key, value in bundled_interactions.items()
    ]

    return JsonResponse({'interactions': bundled_interactions})

# =============================================================================
# IFP Search Results View
# =============================================================================

@method_decorator(csrf_exempt, name='dispatch')
class IFPSearchResults(TemplateView):
    template_name = 'interaction/ifp_search_results.html'

    def post(self, request, *args, **kwargs):
        tanimoto_threshold = float(request.POST.get('tanimoto', 0.10))
        full_ifp = request.POST.get('full_ifp', '')

        print(f"Received full IFP: {full_ifp}") 
        if not full_ifp:
            return render(request, self.template_name, {'error': 'No IFP string provided.'})

        similar_ifps = self.find_similar_ifps(full_ifp, tanimoto_threshold)
        return render(request, self.template_name, {'similar_ifps': similar_ifps})

    def find_similar_ifps(self, full_ifp, tanimoto_threshold):
        all_ifps = ChemokinePartnerIFP.objects.all()
        similar_ifps = []

        for ifp_entry in all_ifps:
            ifp_string = ifp_entry.ifp_string.ljust(len(full_ifp), '0')
            similarity = calculate_tanimoto_similarity(full_ifp, ifp_string)
            if similarity >= tanimoto_threshold:
                similar_ifps.append({
                    'ifp_entry': ifp_entry,
                    'similarity': similarity,
                    'chemokine_entities': [],
                    'partner_entities': []
                })
        return similar_ifps

    #def calculate_tanimoto_similarity(self, ifp1, ifp2):
    #    if len(ifp1) != len(ifp2):
    #        raise ValueError("Bit strings must be of the same length")
#
    #    def bitstring_to_bitvect(bitstring):
    #        bv = ExplicitBitVect(len(bitstring))
    #        for i, bit in enumerate(bitstring):
    #            if bit == '1':
    #                bv.SetBit(i)
    #        return bv
#
    #    bv1 = bitstring_to_bitvect(ifp1)
    #    bv2 = bitstring_to_bitvect(ifp2)
    #    return TanimotoSimilarity(bv1, bv2)


# =============================================================================
# UMAP and IFP Visualization View
# =============================================================================
from django.views import View
from django.shortcuts import render
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import plotly.graph_objs as go
from umap import UMAP

class UMAPIFPPlotView(View):
    template_name = 'interaction/umap_plot.html'
    
    def compute_tanimoto_matrix(self, ifp_strings: list) -> np.ndarray:
        n = len(ifp_strings)
        sim_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                if i == j:
                    sim_matrix[i, j] = 1.0
                else:
                    sim = calculate_tanimoto_similarity(ifp_strings[i], ifp_strings[j])
                    sim_matrix[i, j] = sim
                    sim_matrix[j, i] = sim
        return sim_matrix    

    def get(self, request, *args, **kwargs):
        # 1. Fetch and filter IFP entries
        ifp_entries = list(ChemokinePartnerIFP.objects.all())
        if not ifp_entries:
            return render(request, 'interaction/error.html', {'message': 'No non-zero IFPs found in the database.'})

        filtered_entries = [
            entry for entry in ifp_entries
            if sum(int(bit) for bit in entry.ifp_string) >= 1
        ]

        ifp_strings = [entry.ifp_string for entry in filtered_entries]
        structures = [entry.structure for entry in filtered_entries]

        # 2. Simplify IFPs (1 bit per residue if any interaction occurred)
        simplified_ifps = []
        residue_count = 121  # from 30 to 150 inclusive
        interaction_types_count = 9
        for full_ifp in ifp_strings:
            simplified_bits = []
            for i in range(residue_count):
                residue_slice = full_ifp[i * interaction_types_count : (i + 1) * interaction_types_count]
                simplified_bits.append("1" if "1" in residue_slice else "0")
            simplified_ifps.append("".join(simplified_bits))

        # 3. UMAP projections
        ifp_array = np.array([[int(bit) for bit in fp] for fp in ifp_strings])
        umap_default = UMAP(metric="hamming", n_components=2, min_dist=0.25, random_state=42, n_neighbors=15).fit_transform(ifp_array)

        simplified_ifp_array = np.array([[int(bit) for bit in fp] for fp in simplified_ifps])
        umap_simplified = UMAP(metric="hamming", n_components=2, min_dist=0.25, random_state=42, n_neighbors=15).fit_transform(simplified_ifp_array)

        # UMAP excluding chemokine partners
        non_chemokine_entries = [
            (entry, simp_fp, structure)
            for entry, simp_fp, structure in zip(filtered_entries, simplified_ifps, structures)
            if entry.binding_pair and entry.binding_pair.partner_type != 'Chemokine'
        ]

        if non_chemokine_entries:
            filtered_no_chem_entries, simplified_no_chem_ifps, structures_no_chem = zip(*non_chemokine_entries)
            simplified_array_no_chem = np.array([[int(bit) for bit in fp] for fp in simplified_no_chem_ifps])
            umap_no_chem = UMAP(metric="hamming", n_components=2, min_dist=0.25, random_state=42, n_neighbors=15).fit_transform(simplified_array_no_chem)
        else:
            filtered_no_chem_entries, simplified_no_chem_ifps, structures_no_chem, umap_no_chem = [], [], [], None

        # 4. Metadata DataFrame for default and simplified
        df = pd.DataFrame({
            'UMAP Default 1': umap_default[:, 0],
            'UMAP Default 2': umap_default[:, 1],
            'UMAP Simplified 1': umap_simplified[:, 0],
            'UMAP Simplified 2': umap_simplified[:, 1],
            'State': [s.state for s in structures],
            'PDB Code': [s.pdb_code.index for s in structures],
            'Partner Type': [
                entry.binding_pair.partner_type if entry.binding_pair and entry.binding_pair.partner_type else 'Unknown'
                for entry in filtered_entries
            ],
            'Chemokine Chain': [
                entry.binding_pair.chemokine_chain if entry.binding_pair else 'N/A'
                for entry in filtered_entries
            ],
            'Partner Chain': [
                entry.binding_pair.partner_chain if entry.binding_pair else 'N/A'
                for entry in filtered_entries
            ],
            'Chemokine': [
                entry.binding_pair.structure.protein.gene_name if entry.binding_pair and entry.binding_pair.structure.protein else 'Unknown'
                for entry in filtered_entries
            ],
        })

        unique_partner_types = df['Partner Type'].unique()
        default_colors = go.Figure().layout.template.layout.colorway or [
            '#636EFA', '#EF553B', '#00CC96', '#AB63FA',
            '#FFA15A', '#19D3F3', '#FF6692', '#B6E880',
            '#FF97FF', '#FECB52'
        ]
        partner_type_colors = {pt: default_colors[i % len(default_colors)] for i, pt in enumerate(unique_partner_types)}
        df['PartnerTypeColor'] = df['Partner Type'].map(partner_type_colors)

        def make_trace(x_col, y_col, df_src, title):
            trace = go.Scatter(
                x=df_src[x_col], y=df_src[y_col], mode='markers',
                marker=dict(color=df_src['PartnerTypeColor'], size=10),
                text=df_src.apply(lambda row: (
                    f"Chemokine: {row['Chemokine']}<br>"
                    f"State: {row['State']}<br>"
                    f"Partner Type: {row['Partner Type']}<br>"
                    f"PDB Code: {row['PDB Code']}<br>"
                    f"Chemokine Chain: {row['Chemokine Chain']}<br>"
                    f"Partner Chain: {row['Partner Chain']}"
                ), axis=1),
                hoverinfo='text', showlegend=False
            )
            fig = go.Figure(data=[trace])
            for pt, color in partner_type_colors.items():
                if pt in df_src['Partner Type'].values:
                    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(color=color, size=10), name=pt))
            fig.update_layout(
                title=title,
                width=1000,
                height=800,
                legend=dict(title="Partner Type", x=1.05, y=1),
                xaxis_title=x_col,
                yaxis_title=y_col
            )
            return fig.to_html(full_html=False)

        umap_plot_html_default = make_trace('UMAP Default 1', 'UMAP Default 2', df, 'UMAP (2D) on Default IFPs')
        umap_plot_html_simplified = make_trace('UMAP Simplified 1', 'UMAP Simplified 2', df, 'UMAP (2D) on Simplified IFPs')

        # 5. UMAP plot for non-chemokine partners
        if umap_no_chem is not None:
            df_no_chem = pd.DataFrame({
                'UMAP No Chem 1': umap_no_chem[:, 0],
                'UMAP No Chem 2': umap_no_chem[:, 1],
                'State': [s.state for s in structures_no_chem],
                'PDB Code': [s.pdb_code.index for s in structures_no_chem],
                'Partner Type': [
                    entry.binding_pair.partner_type for entry in filtered_no_chem_entries
                ],
                'Chemokine Chain': [
                    entry.binding_pair.chemokine_chain for entry in filtered_no_chem_entries
                ],
                'Partner Chain': [
                    entry.binding_pair.partner_chain for entry in filtered_no_chem_entries
                ],
                'Chemokine': [
                    entry.binding_pair.structure.protein.gene_name for entry in filtered_no_chem_entries
                ]
            })
            df_no_chem['PartnerTypeColor'] = df_no_chem['Partner Type'].map(partner_type_colors)

            umap_plot_html_no_chem = make_trace('UMAP No Chem 1', 'UMAP No Chem 2', df_no_chem,
                                                'UMAP (2D) on Simplified IFPs – Excluding Chemokine Partners')
        else:
            umap_plot_html_no_chem = '<div>No non-chemokine data available for this projection.</div>'

        return render(
            request,
            self.template_name,
            {
                'umap_plot_html_default': umap_plot_html_default,
                'umap_plot_html_simplified': umap_plot_html_simplified,
                'umap_plot_html_no_chem': umap_plot_html_no_chem,
            }
        )



from django.shortcuts import render
from django.views.generic import TemplateView
from interaction.models import ChemokinePartnerIFP
import numpy as np
import plotly.graph_objects as go
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity  # Ensure TanimotoSimilarity is imported properly


class IFPSimilarityMatrixView(TemplateView):
    template_name = 'interaction/ifp_similarity_matrix.html'

    # GPCR Codes to filter the IFP entries
    GPCR_codes = {
        ('4RWS', 'C', 'A'), ('4XT1', 'B', 'A'), ('4XT3', 'B', 'A'), ('5UIW', 'B', 'A'), 
        ('5WB2', 'B', 'A'), ('6LFM', 'D', 'R'), ('6LFO', 'D', 'R'), ('6WWZ', 'C', 'R'), 
        ('7F1R', 'R', 'R'), ('7O7F', 'I', 'C'), ('7RKF', 'L', 'R'), ('7RKM', 'L', 'R'),
        ('7RKN', 'L', 'R'), ('7SK3', 'B', 'A'), ('7SK4', 'B', 'A'), ('7SK5', 'B', 'A'),
        ('7SK6', 'B', 'A'), ('7SK7', 'B', 'A'), ('7SK8', 'B', 'A'), ('7VL9', 'L', 'R'),
        ('7VLA', 'L', 'R'), ('7XA3', 'L', 'R'), ('7XBX', 'R', 'R'), ('8HNK', 'L', 'R'),
        ('8IC0', 'F', 'A'), ('8JPS', 'C', 'A'), ('8JPS', 'D', 'B'), ('8K3Z', 'D', 'A'),
        ('8K4O', 'F', 'E'), ('8U4O', 'J', 'R'), ('8XVU', 'D', 'R'), ('8XWA', 'D', 'R'),
        ('8XWF', 'D', 'R'), ('8XWM', 'D', 'R'), ('8XWN', 'D', 'R'), ('8XWS', 'D', 'R'),
        ('8XWS', 'B', 'C'), ('8XWV', 'D', 'R'), ('8XX3', 'D', 'R'), ('8XX6', 'D', 'R'),
        ('8XX7', 'D', 'R'), ('8XX7', 'B', 'C'), ('8XXH', 'D', 'R'), ('8XXX', 'D', 'R'),
        ('9AST', 'L', 'R')
    }

    GPCR_NAMES = {
        ('4RWS', 'C', 'A'): '(4RWS_C) vMIP-II - CXCR4',
        ('4XT1', 'B', 'A'): '(4XT1_B) CX3CL1 - US28',
        ('4XT3', 'B', 'A'): '(4XT3_B) CX3CL1 - US28',
        ('5UIW', 'B', 'A'): '(5UIQ_B) [5P7]-CCL5 - CCR5',
        ('5WB2', 'B', 'A'): '(5WB2_B) CX3CL1.35 - US28',
        ('6LFM', 'D', 'R'): '(6LFM_D) CXCL8 - CXCR2',
        ('6LFO', 'D', 'R'): '(6LFO_D) CXCL8 - CXCR2',
        ('6WWZ', 'C', 'R'): '(6WWZ_C) CCL20 - CCR6',
        ('7F1R', 'R', 'R'): '(7F1R_R) CCL5 - CCR5',
        ('7O7F', 'I', 'C'): '(7O7F_I) [6P4]-CCL5 - CCR5',
        ('7RKF', 'L', 'R'): '(7RKF_L) CX3CL1 - US28',
        ('7RKM', 'L', 'R'): '(7RKM_L) CX3CL1 - US28',
        ('7RKN', 'L', 'R'): '(7RKN_L) CX3CL1 - US28',
        ('7SK3', 'B', 'A'): '(7SK3_B) CXCL12 - ACKR3',
        ('7SK4', 'B', 'A'): '(7SK4_B) CXCL12 - ACKR3',
        ('7SK5', 'B', 'A'): '(7SK5_B) CXCL12 - ACKR3',
        ('7SK6', 'B', 'A'): '(7SK6_B) CXCL12 - ACKR3',
        ('7SK7', 'B', 'A'): '(7SK7_B) CXCL12 - ACKR3',
        ('7SK8', 'B', 'A'): '(7SK8_B) CXCL12 - ACKR3',
        ('7VL9', 'L', 'R'): '(7VL9_L) CCL15-(26-92) - CCR1',
        ('7VLA', 'L', 'R'): '(7VLA_L) CCL15-(27-92) - CCR1',
        ('7XA3', 'L', 'R'): '(7XA3_L) CCL2 - CCR2',
        ('7XBX', 'R', 'R'): '(7XBX_R) CX3CL1 - CX3CR1',
        ('8HNK', 'L', 'R'): '(8HNK_L) CXCL11 - CXCR3',
        ('8IC0', 'F', 'A'): '(8IC0_F) CXCL8 - CXCR1',
        ('8JPS', 'C', 'A'): '(8JPS_C) CCL7 - ACKR1',
        ('8JPS', 'D', 'B'): '(8JPS_D) CCL7 - ACKR1',
        ('8K3Z', 'D', 'A'): '(8K3Z_D) CXCL12 - CXCR4',
        ('8K4O', 'F', 'E'): '(8K4O_F) CXCL1 - KSHV-GPCR',
        ('8U4O', 'J', 'R'): '(8U4O_J) CXCL1 - KSHV-GPCR',
        ('8XVU', 'D', 'R'): '(8XVU_D) CXCL2 - CXCR2',
        ('8XWA', 'D', 'R'): '(8XWA_D) CXCL1 - CXCR2',
        ('8XWF', 'D', 'R'): '(8XWF_D) CXCL3 - CXCR2',
        ('8XWM', 'D', 'R'): '(8XWM_D) CXCL6 - CXCR2',
        ('8XWN', 'D', 'R'): '(8XWN_D) CXCL8 - CXCR2',
        ('8XWS', 'D', 'R'): '(8XWS_D) CXCL5 - CXCR2',
        ('8XWS', 'B', 'C'): '(8XWS_B) CXCL5 - CXCR2',
        ('8XWV', 'D', 'R'): '(8XWV_D) CXCL1 - CXCR2',
        ('8XX3', 'D', 'R'): '(8XX3_D) CXCL3 - CXCR2',
        ('8XX6', 'D', 'R'): '(8XX6_D) CXCL8 - CXCR2',
        ('8XX7', 'D', 'R'): '(8XX7_D) CXCL5 - CXCR2',
        ('8XX7', 'B', 'C'): '(8XX7_B) CXCL5 - CXCR2',
        ('8XXH', 'D', 'R'): '(8XXH_D) CXCL2 - CXCR2',
        ('8XXX', 'D', 'R'): '(8XXX_D) CXCL6 - CXCR2',
        ('9AST', 'L', 'R'): '(9AST_L) XCL1 - XCR1'
    }


    def get(self, request, *args, **kwargs):
        # Fetch and filter IFP entries based on GPCR_codes
        ifp_entries = [
            entry for entry in ChemokinePartnerIFP.objects.all()
            if entry.binding_pair and (
                entry.structure.pdb_code.index,
                entry.binding_pair.chemokine_chain,
                entry.binding_pair.partner_chain
            ) in self.GPCR_codes
        ]
        if not ifp_entries:
            return render(request, 'interaction/error.html', {'message': 'No matching IFPs found in the database.'})

        # Process IFP bitstrings
        ifp_strings = pad_bitstrings([entry.ifp_string for entry in ifp_entries])
        bit_vectors = [bitstring_to_bitvect(ifp) for ifp in ifp_strings]
        sim_matrix = self.compute_tanimoto_matrix(bit_vectors)
        reordered_matrix, labels, _ = self.perform_clustering(sim_matrix, ifp_entries)

        # Create heatmap
        heatmap_fig = go.Figure(
            data=go.Heatmap(
                z=reordered_matrix,
                x=labels,
                y=labels,
                colorscale='Viridis',
                colorbar=dict(title='Tanimoto Similarity'),
                hovertemplate="Similarity: %{z}<br>Click for details<extra></extra>",
            )
        )
        heatmap_fig.update_layout(
            title='Clustered Tanimoto Similarity Matrix (GPCR)',
            xaxis=dict(title='IFP Entries', tickangle=45),
            yaxis=dict(title='IFP Entries'),
            width=1000, height=1000
        )

        return render(request, self.template_name, {
            'heatmap_html': heatmap_fig.to_html(full_html=False, include_plotlyjs=False)
        })

    def perform_clustering(self, sim_matrix, ifp_entries):
        """ Perform hierarchical clustering and reorder similarity matrix """
        distance_matrix = 1 - sim_matrix
        linkage_matrix = linkage(squareform(distance_matrix, checks=False), method='average')
        dendro = dendrogram(linkage_matrix, no_plot=True)
        order = dendro['leaves']
        reordered_matrix = sim_matrix[np.ix_(order, order)]
        labels = [
            self.GPCR_NAMES.get(
                (ifp_entries[i].structure.pdb_code.index,
                 ifp_entries[i].binding_pair.chemokine_chain,
                 ifp_entries[i].binding_pair.partner_chain),
                f"{ifp_entries[i].structure.pdb_code.index} ({ifp_entries[i].binding_pair.chemokine_chain}-{ifp_entries[i].binding_pair.partner_chain})"
            ) for i in order
        ]
        return reordered_matrix, labels, dendro


    def compute_tanimoto_matrix(self, bit_vectors):
        """ Compute Tanimoto similarity matrix """
        n = len(bit_vectors)
        sim_matrix = np.eye(n)
        for i in range(n):
            for j in range(i + 1, n):
                sim = TanimotoSimilarity(bit_vectors[i], bit_vectors[j])
                sim_matrix[i, j] = sim_matrix[j, i] = sim
        return sim_matrix



from django.views.generic import TemplateView
from interaction.models import ChemokinePartnerIFP
import numpy as np
import plotly.graph_objects as go
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

class IFPSimilarityMatrixAntibodyView(TemplateView):
    template_name = 'interaction/ifp_similarity_matrix_antibody.html'

    def get(self, request, *args, **kwargs):
        # Filter IFPs where the partner type is "antibody"
        ifp_entries = [
            entry for entry in ChemokinePartnerIFP.objects.all()
            if entry.binding_pair and entry.binding_pair.partner_type == 'Antibody'
        ]

        print(len(ifp_entries))


        if not ifp_entries:
            return render(request, 'interaction/error.html', {'message': 'No antibody-related IFPs found.'})

        # Convert IFP strings to bit vectors
        ifp_strings = pad_bitstrings([entry.ifp_string for entry in ifp_entries])
        bit_vectors = [bitstring_to_bitvect(ifp) for ifp in ifp_strings]
        sim_matrix = self.compute_tanimoto_matrix(bit_vectors)

        reordered_matrix, labels, _ = self.perform_clustering(sim_matrix, ifp_entries)

        heatmap_fig = go.Figure(
            data=go.Heatmap(
                z=reordered_matrix,
                x=labels,
                y=labels,
                colorscale='Viridis',
                colorbar=dict(title='Tanimoto Similarity'),
                hovertemplate="Similarity: %{z}<br><extra></extra>",
            )
        )
        heatmap_fig.update_layout(
            title='Tanimoto Similarity Matrix (Antibody Complexes)',
            xaxis=dict(title='Complex'),
            yaxis=dict(title='Complex'),
            width=1000,
            height=1000,
        )

        return render(request, self.template_name, {
            'heatmap_html': heatmap_fig.to_html(full_html=False, include_plotlyjs=False)
        })

    def perform_clustering(self, sim_matrix, ifp_entries):
        distance_matrix = 1 - sim_matrix
        linkage_matrix = linkage(squareform(distance_matrix, checks=False), method='average')
        dendro = dendrogram(linkage_matrix, no_plot=True)
        order = dendro['leaves']
        reordered_matrix = sim_matrix[np.ix_(order, order)]
        labels = [
            f"{ifp.structure.pdb_code.index} ({ifp.binding_pair.chemokine_chain}-{ifp.binding_pair.partner_chain})"
            for i, ifp in enumerate(ifp_entries)
        ]
        labels = [labels[i] for i in order]
        return reordered_matrix, labels, dendro

    def compute_tanimoto_matrix(self, bit_vectors):
        n = len(bit_vectors)
        sim_matrix = np.eye(n)
        for i in range(n):
            for j in range(i + 1, n):
                sim = TanimotoSimilarity(bit_vectors[i], bit_vectors[j])
                sim_matrix[i, j] = sim_matrix[j, i] = sim
        return sim_matrix
