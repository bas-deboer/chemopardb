from django.shortcuts import render
from .forms import PDBSelectionForm
from django.shortcuts import render
from django.http import HttpResponse, FileResponse
from django.template.loader import render_to_string
from django.views.generic import View, TemplateView
from django.db.models import Count, Q, Prefetch, TextField
from django.conf import settings
from django.utils.html import format_html
from django.shortcuts import render, redirect
from django.views.generic import TemplateView

from common.diagrams_chemokine import *
from structure.utils import *

from structure.models import Structure, PdbData, Chain, Rotamer
from protein.models import Protein
from residue.models import Residue
from interaction.models import ChemokinePartnerInteraction, ChemokinePartnerIFP

import requests
import json
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
import xlsxwriter


def SelectStructure(request):
    if request.method == 'POST':
        form = PDBSelectionForm(request.POST)
        if form.is_valid():
            selected_pdb_ids = form.cleaned_data['pdb_ids']
            # Redirect to the results page with the selected PDB IDs
            return display_alignment_results(request, selected_pdb_ids)
    else:
        form = PDBSelectionForm()
    return render(request, 'interaction/select_structure.html', {'form': form})

def display_alignment_results(request, pdb_ids):
    # Assuming 'align_sequences' is a function that takes PDB IDs and returns an alignment result
    alignment_result = align_sequences(pdb_ids)
    context = {
        'alignment_result': alignment_result,
        'pdb_ids': pdb_ids,
    }
    return render(request, 'pdbapp/display_alignment_results.html', context)


class InteractionDetails(TemplateView):
    template_name = 'interaction/interaction_details_staticresidues.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        pdb_id = self.kwargs.get('pdb_id', '')
        chain_id = self.kwargs.get('chain_id', '')

        try:
            structure = Structure.objects.get(pdb_code__index=pdb_id)

            context['structure'] = structure
            context['protein'] = structure.protein
            context['partner'] = structure.partner.all()

            # Prepare residues data for the context
            rotamers = Rotamer.objects.filter(structure=structure, chain=chain_id).select_related('residue')
            context['residues_data'] = rotamers

            # Get the IDs of these rotamers
            rotamer_ids = rotamers.values_list('id', flat=True)
            # Fetch interactions that are linked to the fetched rotamers and filter out self-interactions
            interactions = ChemokinePartnerInteraction.objects.filter(
                chemokine_residue_id__in=rotamer_ids
            ).exclude(
                partner_chain=chain_id
            )
            context['interactions'] = interactions

            # Extract sequence numbers using a list comprehension
            sequence_numbers = [interaction.chemokine_residue.pdbseq_number for interaction in interactions if interaction.chemokine_residue and interaction.chemokine_residue.residue]
            unique_sequence_numbers = list(set(sequence_numbers))
            context['display_res'] = unique_sequence_numbers

            pdb_data = None

            interaction_data = None
            csv_file_path = None
            chart = None
            dist_matrix = None
            
            # Parse IFP data from structure and use as input for table
            base_path = os.path.join(settings.DATA_DIR, 'prepared_pdbs', f"{pdb_id}_prepared")
            if os.path.exists(base_path):
                csv_file = glob.glob(os.path.join(base_path, '*.csv'))
                if csv_file:
                    csv_file_path = csv_file[0]
                    data = pd.read_csv(csv_file_path, header=None)
                    interaction_data = []
                    for index, row in data.iterrows():
                        interaction_row = {
                            'ChemokineResidue': row[0],
                            'ChemokineSegment': '-',
                            'PartnerResidue': row[1],
                            'InteractionType': row[2]
                        }
                        interaction_data.append(interaction_row)

                    ### Make plotly contact map ###
                    binary_matrix = csv_to_matrix(csv_file)

                    fig = px.imshow(binary_matrix, 
                                    labels=dict(x="Residue Number", color="Contact"),
                                    x=np.arange(binary_matrix.shape[1]),
                                    y=np.arange(binary_matrix.shape[0]),
                                    color_continuous_scale="Blues",
                                    )

                    fig.update(layout_coloraxis_showscale=False)
                    fig.update_traces({'xgap': 5, 'ygap':5})
                    fig.update_layout(title="Contact Map", xaxis_nticks=36, yaxis_nticks=36)
                    fig.update_xaxes(side="top")
                    fig.update_yaxes(showticklabels=False)

                    # Disable interactive features like panning and zooming
                    fig.update_layout(dragmode=False, xaxis_fixedrange=True, yaxis_fixedrange=True)

                    chart = fig.to_html()   

                    #### Create distance matrix ###
                    #chemokine_file_path = os.path.join(base_path, f"{pdb_id}_chemokine.pdb")
                    #print(chemokine_file_path)
                    #distance_matrix = calc_distance_matrix(chemokine_file_path)
                    #print(distance_matrix)
#
                    #fig = px.imshow(distance_matrix,
                    #                labels=dict(x="Residue", y="Residue", color="Distance"),
                    #                x=np.arange(distance_matrix.shape[1]),
                    #                y=np.arange(distance_matrix.shape[0]),
                    #                color_continuous_scale="Blues_r")
#
                    #fig.update_layout(title="Cα Distance Matrix", xaxis_nticks=20, yaxis_nticks=20)
                    #fig.update_xaxes(side="top")  # Move the x-axis to the top if preferred
#
                    #dist_matrix = fig.to_html() 

            # Path to the network HTML file
            html_file_path = os.path.join(settings.DATA_DIR, 'interactions_output', pdb_id, f"{pdb_id}_network.html")
            context["network_html"] = None
            if os.path.exists(html_file_path):
                with open(html_file_path, 'r') as file:
                    context["network_html"] = file.read()     
        
                                     
            # Construct the path for the PNG file
            png_file_path = os.path.join(settings.DATA_DIR, 'interactions_output', pdb_id, f"{pdb_id}_interaction_matrix.png")

            # Check if the PNG file exists
            png_exists = os.path.exists(png_file_path)
            
            ## Process sequence for displaying in IFP table
            #input_sequence = "KPVSLSYRC--PCRFFES-HVA---RANVKHLKILNTP-NC-A-LQIVARLK----NNNRQVCIDPKLKWIQEYLEKALNKRFKM-----"
            #sequence_data = parse_sequence(input_sequence)
            #context['sequence_data'] = sequence_data
            
            
            
            # Fetch IFPs and organize them by chemokine_partner_pair
            ifp_entries = ChemokinePartnerIFP.objects.filter(
                chemokine_partner_pair__structure=structure,
                chemokine_partner_pair__chemokine_chain=chain_id
            ).select_related('chemokine_partner_pair')

            # Create a dictionary to store unique pairs
            unique_pairs = {}
            possible_interactions = [
                "Hydrophobic", "Pistacking", "HBDonor", "HBAcceptor",
                "Anionic", "Cationic", "CationPi", "PiCation", "VdWContact"
            ]
            for ifp_entry in ifp_entries:
                pair = ifp_entry.chemokine_partner_pair
                # Omit pairs where both chains have the same ID
                if pair.chemokine_chain != pair.partner_chain:
                    pair_key = f"{pair.chemokine_chain}-{pair.partner_chain}"
                    if pair_key not in unique_pairs:
                        unique_pairs[pair_key] = pair

            ifp_data_by_pair = {}
            for pair_key, pair in unique_pairs.items():
                ifp_entries = ChemokinePartnerIFP.objects.filter(chemokine_partner_pair=pair)
                ifp_data = []
                for ifp_entry in ifp_entries:
                    ifp_string = ifp_entry.ifp_string
                    for i, generic_number in enumerate(range(30, 151)):
                        interaction_bits = ifp_string[i*9:(i+1)*9]
                        interactions = [possible_interactions[j] for j, bit in enumerate(interaction_bits) if bit == '1']
                        ifp_data.append({
                            'generic_number': generic_number,
                            'ifp_string': interaction_bits,
                            'interactions': ", ".join(interactions) if interactions else "None"
                        })
                ifp_data_by_pair[pair_key] = ifp_data

            context['ifp_data_by_pair'] = ifp_data_by_pair


        except Structure.DoesNotExist:
            return render(self.request, 'error.html')
        
        return context
    
def excel(request, slug, **response_kwargs):
    structure = Structure.objects.get(pdb_code__index=slug)

    # Prepare residues data for the context
    rotamers = Rotamer.objects.filter(structure=structure).select_related('residue')
    
    # Get the IDs of these rotamers
    rotamer_ids = rotamers.values_list('id', flat=True)
    
    # Fetch interactions that are linked to the fetched rotamers
    interactions = ChemokinePartnerInteraction.objects.filter(chemokine_residue_id__in=rotamer_ids) #.order_by('rotamer__residue__sequence_number')
    print(interactions)
    data = []
    
    for interaction in interactions:
        row = {}
        row['Sequence Number'] = interaction.chemokine_residue.residue.sequence_number
        row['Amino Acid'] = interaction.chemokine_residue.residue.amino_acid
        row['Segment'] = interaction.chemokine_residue.residue.segment
        #if interaction.chemokine_residue.residue.display_generic_number:
        #    row['Generic Number'] = interaction.chemokine_residue.residue.display_generic_number.label
        #    row['Segment'] = interaction.chemokine_residue.residue.protein_segment.slug
        #else:
        #    row['Generic Number'] = 'N/A'
        #    row['Segment'] = '-'

        row['Interaction'] = interaction.interaction_type
        #row['Interaction Slug'] = interaction.interaction_type
        #row['Ligand'] = interaction.structure_ligand_pair.ligand.name
        row['Partner'] = "Partner"

        data.append(row)


    headers = ['Partner','Amino Acid','Sequence Number','Segment','Interaction']

    #EXCEL SOLUTION
    output = BytesIO()
    workbook = xlsxwriter.Workbook(output)
    worksheet = workbook.add_worksheet()

    col = 0
    for h in headers:
        worksheet.write(0, col, h)
        col += 1
    row = 1
    for d in data:
        col = 0
        for h in headers:
            worksheet.write(row, col, d[h])
            col += 1
        row += 1
    workbook.close()
    output.seek(0)
    xlsx_data = output.read()

    response = HttpResponse(xlsx_data,content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    response['Content-Disposition'] = 'attachment; filename=Interaction_data_%s.xlsx' % slug
    return response

def ajax(request, slug, **response_kwargs):

    residuelist = Residue.objects.filter(protein__name=slug)
    lookup = {}

    for r in residuelist:
        if r.generic_number:
            lookup[r.generic_number] = r.sequence_number

    interactions = ChemokinePartnerInteraction.objects.filter(chemokine_residue__structure__pdb_code__index="2N55")
    # return HttpResponse("Hello, world. You're at the polls index. "+slug)
    jsondata = {}
    for interaction in interactions:
        if interaction.chemokine_residue.residue.generic_number:
            sequence_number = interaction.chemokine_residue.residue.sequence_number
            if interaction.chemokine_residue.residue.generic_number in lookup:
                sequence_number = lookup[interaction.chemokine_residue.residue.generic_number]

            #label = interaction.rotamer.residue.generic_number.label
            aa = interaction.chemokine_residue.residue.amino_acid
            interactiontype = interaction.interaction_type
            if sequence_number not in jsondata:
                jsondata[sequence_number] = []
            jsondata[sequence_number].append([aa, interactiontype])

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)


from django.shortcuts import render
from django.views.generic import TemplateView
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import ExplicitBitVect, TanimotoSimilarity
import numpy as np

class IFPSearchResults(TemplateView):
    template_name = 'interaction/ifp_search_results.html'

    def post(self, request, *args, **kwargs):
        tanimoto_threshold = float(request.POST.get('tanimoto'))
        full_ifp = request.POST.get('full_ifp')

        # Find similar IFPs based on the Tanimoto threshold
        similar_ifps = self.find_similar_ifps(full_ifp, tanimoto_threshold)

        return render(request, self.template_name, {'similar_ifps': similar_ifps})

    def find_similar_ifps(self, full_ifp, tanimoto_threshold):
        all_ifps = ChemokinePartnerIFP.objects.all()
        similar_ifps = []

        for ifp_entry in all_ifps:
            chemokine_chain = ifp_entry.chemokine_partner_pair.chemokine_chain
            partner_chain = ifp_entry.chemokine_partner_pair.partner_chain

            if chemokine_chain == partner_chain:
                continue

            ifp_string = ifp_entry.ifp_string.ljust(len(full_ifp), '0')
            similarity = self.calculate_tanimoto_similarity(full_ifp, ifp_string)
            if similarity >= tanimoto_threshold:
                similar_ifps.append((ifp_entry, similarity))

        return similar_ifps

    def calculate_tanimoto_similarity(self, ifp1, ifp2):
        if len(ifp1) != len(ifp2):
            raise ValueError("Bit strings must be of the same length")

        def bitstring_to_bitvect(bitstring):
            bitvect = ExplicitBitVect(len(bitstring))
            for i, bit in enumerate(bitstring):
                if bit == '1':
                    bitvect.SetBit(i)
            return bitvect

        bv1 = bitstring_to_bitvect(ifp1)
        bv2 = bitstring_to_bitvect(ifp2)

        return TanimotoSimilarity(bv1, bv2)


from django.shortcuts import render
import pandas as pd


class ViewAlignment(TemplateView):
    template_name = 'interaction/view_alignment.html'

    def post(self, request, *args, **kwargs):
        pdb_indices = request.POST.get('pdb_indices', '')
        pdb_ids = pdb_indices.split(',')
        if not pdb_ids:
            return render(request, 'error.html', {'message': 'No PDB indices provided.'})

        # Assuming 'align_sequences' is a function that takes PDB IDs and returns an alignment result
        alignment_result = align_sequences(pdb_ids)

        context = {
            'alignment_result': alignment_result,
            'pdb_ids': pdb_ids,
        }
        return self.render_to_response(context)




    def get(self, request, *args, **kwargs):
        ids = kwargs.get('ids', '')
        pdb_ids = ids.split(',')
        if not pdb_ids:
            return render(request, {'message': 'No PDB IDs provided.'})

        # Assuming 'align_sequences' is a function that takes PDB IDs and returns an alignment result
        alignment_result = align_sequences(pdb_ids)

        context = {
            'alignment_result': alignment_result,
            'pdb_ids': pdb_ids,
        }
        return self.render_to_response(context)


def align_sequences(pdb_ids):
    # Implement your alignment logic here
    # For example, you might use an external tool or library to perform the alignment
    # This is a placeholder implementation
    alignment_result = "Aligned sequences for PDB IDs: " + ", ".join(pdb_ids)
    return alignment_result