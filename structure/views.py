import csv
from django.shortcuts import render, get_object_or_404
from django.http import HttpResponse, FileResponse, Http404
from django.template.loader import render_to_string
from django.views.generic import View, TemplateView, DetailView
from django.db.models import Count, Q, Prefetch, TextField
from django.conf import settings
from django.utils.html import format_html
from django.utils.text import slugify
from django.urls import reverse

from common.diagrams_chemokine import *
from structure.utils import *

from common.models import ResiduePosition
from structure.models import Structure, PdbData, Rotamer, Entity, EntityType, ChemokineBindingPartner, PredictedStructure, PredictedComplex
from protein.models import Protein
from residue.models import Residue
from interaction.models import ChemokinePartnerInteraction, PredictedChemokinePartnerInteraction
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

from collections import OrderedDict

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

        selected_method = self.request.GET.get('method')
        if selected_method is not None:
            selected_method = selected_method.strip() or None

        selected_state = self.request.GET.get('state')
        if selected_state is not None:
            selected_state = selected_state.strip() or None

        structures = Structure.objects.all().select_related(
            'protein',
            'structure_type',
            'pdb_code',
            'publication__web_link',
        )

        context.update({
            'structures': structures,
            'selected_method': selected_method,
            'selected_state': selected_state,
        })
        return context


class StructureDetails(TemplateView):
    template_name = 'structure/structure.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        structure_id = self.kwargs.get('structure_id')

        structure = get_object_or_404(
            Structure.objects.select_related(
                'protein',
                'pdb_code',
                'structure_type',
                'publication__web_link',
            ),
            id=structure_id,
        )

        protein = structure.protein
        protein_residues = list(
            Residue.objects.filter(protein=protein).order_by('sequence_number')
        ) if protein else []
        canonical_sequence = ''.join(res.amino_acid for res in protein_residues) if protein_residues else ''

        all_positions = list(ResiduePosition.objects.all().order_by('position'))
        pos_lookup = {pos.ccn_number: pos for pos in all_positions if pos.ccn_number}

        protein_dict = {}
        for res in protein_residues:
            ccn = getattr(res, 'ccn_number', None)
            if ccn and ccn in pos_lookup:
                protein_dict[pos_lookup[ccn].position] = res

        structure_rotamers = list(
            Rotamer.objects.filter(structure=structure).order_by('sequence_number')
        )
        rotamers_by_chain = {}
        for rot in structure_rotamers:
            ccn = getattr(rot, 'ccn_number', None)
            chain_id = getattr(rot, 'chain', None)
            if not ccn or ccn not in pos_lookup or not chain_id:
                continue
            position_id = pos_lookup[ccn].position
            rotamers_by_chain.setdefault(chain_id, {})[position_id] = rot

        aligned_sequences = []
        for pos in all_positions:
            position_id = pos.position
            protein_residue = protein_dict.get(position_id)
            structure_rotamers_map = {
                chain_id: chain_map.get(position_id)
                for chain_id, chain_map in rotamers_by_chain.items()
            }
            aligned_sequences.append((pos, protein_residue, structure_rotamers_map))

        residues_by_segment = OrderedDict()
        for res in protein_residues:
            segment = getattr(res, 'segment', None) or '-'
            residues_by_segment.setdefault(segment, []).append(res)

        chains = sorted(rotamers_by_chain.keys())

        entity_types = EntityType.objects.prefetch_related(
            Prefetch(
                'entity_set',
                queryset=Entity.objects.filter(structure=structure).order_by('name'),
                to_attr='filtered_entities',
            )
        )

        binding_partner_entries = []
        for partner in structure.chemokine_binding_partners.annotate(
            interaction_count=Count('chemokinepartnerinteraction')
        ):
            binding_partner_entries.append({
                'id': partner.id,
                'chemokine_name': partner.chemokine_name,
                'chemokine_chain': partner.chemokine_chain,
                'partner_name': partner.partner_name,
                'partner_chain': partner.partner_chain,
                'interaction_count': getattr(partner, 'interaction_count', 0),
                'url': reverse('interaction:chemokine_binding_partner_detail', args=[partner.id]),
            })

        context.update({
            'structure': structure,
            'protein': protein,
            'canonical_sequence': canonical_sequence,
            'aligned_sequences': aligned_sequences,
            'residues_by_segment': residues_by_segment,
            'chains': chains,
            'entity_types': entity_types,
            'binding_partners': binding_partner_entries,
            'binding_partner_count': len(binding_partner_entries),
            'chain_count': len(chains),
            'residue_count': len(protein_residues),
            'rotamer_count': len(structure_rotamers),
            'pdb_available': bool(structure.pdb_data and structure.pdb_data.pdb),
            'structure_method': structure.structure_type.name if structure.structure_type else None,
            'structure_resolution': structure.resolution,
            'structure_state': structure.state,
            'structure_publication_date': structure.publication_date,
            'all_positions': all_positions,
        })

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
            file_path = os.path.join(
                settings.DATA_DIR,
                'prepared_pdbs',
                f"{pdb_id}_prepared",
                f"{pdb_id}_all_interactions.csv",
            )

            if os.path.exists(file_path):
                df = pd.read_csv(file_path, header=None)
                first_filter = df.columns[df.iloc[0].str.endswith(selected_chain)]
                first_filtered_df = df[first_filter]
                second_filter = first_filtered_df.columns[~first_filtered_df.iloc[1].str.endswith(selected_chain)]
                final_filtered_df = first_filtered_df[second_filter]
                interactions = final_filtered_df.to_html(
                    classes=["table", "table-bordered", "table-striped"],
                    index=False,
                )
            else:
                interactions = "File not found."

            context['interactions'] = interactions

        context['form'] = form
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
    

class PredictedStructureBrowser(TemplateView):
    """
    Browse all predicted (e.g., AlphaFold) chemokine structures.
    """
    template_name = "structure/predicted_structure_browser.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['predicted_structures'] = PredictedStructure.objects.select_related('protein').all()
        return context



class PredictedPDBDownload(View):
    """Serve the PDB file for a predicted structure."""

    def get(self, request, pk, *args, **kwargs):
        predicted_structure = get_object_or_404(PredictedStructure, pk=pk)

        if not predicted_structure.pdb_data or not predicted_structure.pdb_data.pdb:
            return HttpResponse('PDB data not available for this predicted structure.', status=404)

        pdb_content = predicted_structure.pdb_data.pdb
        protein = predicted_structure.protein
        base_name = (getattr(protein, 'gene_name', None) or getattr(protein, 'name', None) or 'predicted_structure')
        filename = f"{slugify(base_name) or 'predicted_structure'}_{pk}.pdb"

        response = HttpResponse(pdb_content, content_type='chemical/x-pdb')
        response['Content-Disposition'] = f'attachment; filename="{filename}"'
        return response


class PredictedStructureDetails(DetailView):
    """
    Show details for a single predicted structure.
    """
    model = PredictedStructure
    template_name = "structure/predicted_structure_details.html"
    context_object_name = "predicted_structure"

class PredictedComplexBrowser(TemplateView):
    """
    Browse all predicted chemokine-GPCR complex models.
    """
    template_name = "structure/predicted_complex_browser.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['predicted_complexes'] = PredictedComplex.objects.select_related('chemokine_protein', 'gpcr_partner').all()
        return context

class PredictedComplexDetails(DetailView):
    """
    Show details for a single predicted complex.
    """
    model = PredictedComplex
    template_name = "structure/predicted_complex_details.html"
    context_object_name = "predicted_complex"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['interactions_url'] = reverse('predicted_complex_interactions', args=[self.object.pk])
        context['interactions_csv_url'] = reverse('predicted_complex_interactions_csv', args=[self.object.pk])

        metric_specs = [
            ('ipTM', self.object.iptm, 2),
            ('pTM', self.object.ptm, 2),
            ('Mean pLDDT', self.object.mean_plddt, 2),
            ('Mean pLDDT (chain A)', getattr(self.object, 'mean_plddt_A', None), 2),
            ('Mean pLDDT (chain B)', getattr(self.object, 'mean_plddt_B', None), 2),
            ('Mean PAE', self.object.mean_pae, 2),
            ('Mean iPAE', self.object.mean_ipae, 2),
        ]
        context['metric_values'] = [
            {'label': label, 'value': value, 'precision': precision}
            for label, value, precision in metric_specs
            if value is not None
        ]
        return context

class PredictedComplexInteractions(DetailView):
    """Display residue-level interaction fingerprints for a predicted complex."""

    model = PredictedComplex
    template_name = "structure/predicted_complex_interactions.html"
    context_object_name = "predicted_complex"

    def get_queryset(self):
        return (
            super()
            .get_queryset()
            .select_related("chemokine_protein", "gpcr_partner")
            .prefetch_related("predicted_interactions")
        )

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        interactions_qs = self.object.predicted_interactions.all()

        context["interactions"] = interactions_qs.order_by(
            "chemokine_chain",
            "chemokine_residue_number",
            "partner_chain",
            "partner_residue_number",
            "partner_generic_number",
            "interaction_type",
        )

        context["interaction_counts"] = (
            interactions_qs.values("interaction_type")
            .order_by("interaction_type")
            .annotate(total=Count("id"))
        )

        context["has_interactions"] = interactions_qs.exists()
        context["interaction_total"] = interactions_qs.count()
        context["generic_annotated"] = interactions_qs.filter(
            partner_generic_number__isnull=False
        ).exclude(partner_generic_number="").count()

        context["details_url"] = reverse('predicted_complex_details', args=[self.object.pk])

        chem_surface = set()
        partner_surface = set()

        for interaction in interactions_qs:
            chem_chain = (interaction.chemokine_chain or "").strip()
            partner_chain = (interaction.partner_chain or "").strip()

            chem_number = interaction.chemokine_residue_number
            partner_number = interaction.partner_residue_number

            try:
                chem_number = int(chem_number)
            except (TypeError, ValueError):
                chem_number = None

            try:
                partner_number = int(partner_number)
            except (TypeError, ValueError):
                partner_number = None

            if chem_chain and chem_number:
                chem_surface.add((chem_chain, chem_number))

            if partner_chain and partner_number:
                partner_surface.add((partner_chain, partner_number))

        context["chemokine_surface"] = [
            {"chain": chain, "number": number}
            for chain, number in sorted(chem_surface)
        ]
        context["partner_surface"] = [
            {"chain": chain, "number": number}
            for chain, number in sorted(partner_surface)
        ]

        context["show_viewer"] = bool(
            self.object.pdb_data and self.object.pdb_data.pdb
        )
        return context


def predicted_complex_interactions_csv(request, pk):
    predicted_complex = get_object_or_404(PredictedComplex, pk=pk)
    interactions = (
        PredictedChemokinePartnerInteraction.objects
        .filter(predicted_complex=predicted_complex)
        .select_related('chemokine_residue')
        .order_by(
            'chemokine_chain',
            'chemokine_residue_number',
            'partner_chain',
            'partner_residue_number',
            'partner_generic_number',
            'interaction_type',
        )
    )

    chemokine_protein = predicted_complex.chemokine_protein
    gpcr_partner = predicted_complex.gpcr_partner
    chemokine_label = slugify(getattr(chemokine_protein, 'gene_name', None) or getattr(chemokine_protein, 'name', None) or 'chemokine') or 'chemokine'
    partner_label = slugify(getattr(gpcr_partner, 'name', None) or 'gpcr') or 'gpcr'
    filename = f"{chemokine_label}_{partner_label}_predicted_interactions.csv"

    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = f'attachment; filename="{filename}"'

    writer = csv.writer(response)
    writer.writerow([
        'Chemokine residue',
        'Chemokine chain',
        'Chemokine CCN',
        'Partner residue',
        'Partner chain',
        'GPCRdb generic',
        'Interaction type',
    ])

    for interaction in interactions:
        chem_name = interaction.chemokine_residue_name or ''
        chem_number = interaction.chemokine_residue_number
        chem_chain = interaction.chemokine_chain or ''
        chem_ccn = ''
        if interaction.chemokine_residue and interaction.chemokine_residue.ccn_number:
            chem_ccn = interaction.chemokine_residue.ccn_number
        partner_name = interaction.partner_residue_name or ''
        partner_number = interaction.partner_residue_number
        partner_chain = interaction.partner_chain or ''
        partner_generic = interaction.partner_generic_number or ''
        interaction_type = interaction.interaction_type or ''

        chem_residue = str(chem_name).strip()
        if chem_number is not None:
            chem_residue = f"{chem_residue} {chem_number}".strip()
        partner_residue = str(partner_name).strip()
        if partner_number is not None:
            partner_residue = f"{partner_residue} {partner_number}".strip()

        writer.writerow([
            chem_residue or '-',
            chem_chain or '-',
            chem_ccn or '-',
            partner_residue or '-',
            partner_chain or '-',
            partner_generic or '-',
            interaction_type or '-',
        ])

    return response


