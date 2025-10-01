from django.shortcuts import get_object_or_404, render, redirect
from django.views import generic
from django.http import JsonResponse, HttpResponse
from django.db.models import Q, F, Func, Value, Prefetch
from django.core.cache import cache
from django.views.decorators.cache import cache_page
from django.urls import reverse
from django.views.generic import View, TemplateView
from django.urls import reverse

import requests
import xml.etree.ElementTree as ET
from datetime import date

from common.models import ResiduePosition
from protein.models import Protein, ProteinSegment, ProteinAlias
from interaction.models import ChemokineGpcrGtPBinding
from structure.models import Structure, Rotamer, PredictedStructure, PredictedComplex
from residue.models import Residue, ResidueGenericNumber


def index(request):
    return HttpResponse("Protein page")

class ProteinBrowser(TemplateView):
    """
    Fetching Protein data for browser
    """
    template_name = "protein_browser.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        proteins = Protein.objects.all()
        context['proteins'] = proteins
        return context



def protein(request, name):
    try:
        protein = Protein.objects.get(gene_name=name)
    except Protein.DoesNotExist:
        return render(request, 'error.html')

    # Retrieve related data
    structures = Structure.objects.filter(protein=protein)
    residues = Residue.objects.filter(protein=protein).order_by('sequence_number')
    # Build raw sequence strings for copy/download
    sequence_text = "".join(r.amino_acid for r in residues)
    sequence_fasta = f">{protein.gene_name}\n{sequence_text}"
    predicted_structure_models = list(PredictedStructure.objects.filter(protein=protein).order_by('-date_generated', '-id'))
    predicted_complexes = list(PredictedComplex.objects.filter(chemokine_protein=protein).select_related('gpcr_partner').order_by('-date_generated', '-id'))
    predicted_structures = sorted(
        predicted_structure_models + predicted_complexes,
        key=lambda item: ((item.date_generated or date.min), item.pk),
        reverse=True,
    )
    alternate_names = ProteinAlias.objects.filter(protein=protein)
    # GtP bindings for this chemokine
    gtp_bindings = ChemokineGpcrGtPBinding.objects.filter(chemokine=protein).select_related('gpcr_partner')

    # ResiduePosition map: ccn_number <-> position
    position_map = {rp.ccn_number: rp.position for rp in ResiduePosition.objects.all()}
    ccn_map = {rp.position: rp.ccn_number for rp in ResiduePosition.objects.all()}

    # Prepare aligned sequences: protein_dict and rotamers_by_structure use ccn_number as key
    protein_dict = {res.ccn_number: res for res in residues if res.ccn_number}
    rotamers = Rotamer.objects.filter(structure__in=structures).order_by("sequence_number")
    rotamers_by_structure = {}
    for rotamer in rotamers:
        structure = rotamer.structure
        if structure not in rotamers_by_structure:
            rotamers_by_structure[structure] = {}
        rotamers_by_structure[structure][rotamer.ccn_number] = rotamer

    # Get all unique alignment positions from both protein and rotamers
    all_positions = set()
    for ccn in protein_dict.keys():
        pos = position_map.get(ccn)
        if pos is not None:
            all_positions.add(pos)
    for structure in rotamers_by_structure:
        for ccn in rotamers_by_structure[structure].keys():
            pos = position_map.get(ccn)
            if pos is not None:
                all_positions.add(pos)
    all_positions = sorted(all_positions)

    # Build the aligned_sequences list: (position, ccn_number, protein_residue, structure_rotamers)
    aligned_sequences = []
    for pos in all_positions:
        ccn = ccn_map.get(pos)
        protein_residue = protein_dict.get(ccn)
        structure_rotamers = {structure: rotamers_by_structure[structure].get(ccn) for structure in structures}
        aligned_sequences.append((pos, ccn, protein_residue, structure_rotamers))

    alternate_names_str = ', '.join([alias.name for alias in alternate_names])

    # Group residues by segment and sort by position
    residues_by_segment = {}
    for residue in residues:
        segment = residue.segment
        if segment not in residues_by_segment:
            residues_by_segment[segment] = []
        residues_by_segment[segment].append(residue)
    # Sort each segment's residues by position
    for segment, seg_residues in residues_by_segment.items():
        residues_by_segment[segment] = sorted(
            seg_residues, key=lambda res: position_map.get(res.ccn_number, 0)
        )

    # Build compact sequence viewer chunks without horizontal scrolling
    chunk_size = 10
    residues_list = list(residues)
    seq_chunks = []
    for i in range(0, len(residues_list), chunk_size):
        chunk_res = residues_list[i:i+chunk_size]
        # Build run-length encoded segments for header row
        segments = []
        def seg_category(seg_name: str):
            if seg_name == 'N-term':
                return 'nterm'
            if seg_name == 'C-term':
                return 'cterm'
            if seg_name == 'CX':
                return 'cx'
            if seg_name == 'Helix':
                return 'helix'
            if seg_name and seg_name.startswith('B'):
                return 'beta'
            return 'loop'

        if chunk_res:
            current_seg = chunk_res[0].segment or '-'
            # Determine if the first segment in this chunk is a continuation from previous chunk
            prev_seg = (residues_list[i-1].segment if i > 0 else None) or '-'
            continuing_from_prev = (prev_seg == current_seg)

            run_len = 0
            first_span = True
            for r in chunk_res:
                seg = r.segment or '-'
                if seg == current_seg:
                    run_len += 1
                else:
                    segments.append({
                        'name': current_seg,
                        'colspan': run_len,
                        'category': seg_category(current_seg),
                        'show_label': (False if first_span and continuing_from_prev else True),
                    })
                    current_seg = seg
                    run_len = 1
                    first_span = False
            # append last span
            segments.append({
                'name': current_seg,
                'colspan': run_len,
                'category': seg_category(current_seg),
                'show_label': (False if first_span and continuing_from_prev else True),
            })

        seq_chunks.append({
            'residues': chunk_res,
            'segments': segments,
            'end_number': chunk_res[-1].sequence_number if chunk_res else None,
        })

    context = {
        'protein': protein,
        'structures': structures,
        'predicted_structures': predicted_structures,
        'predicted_structure_models': predicted_structure_models,
        'predicted_complexes': predicted_complexes,
        'aligned_sequences': aligned_sequences,
        'residues_by_segment': residues_by_segment,
        'alternate_names': alternate_names_str,
        'gtp_bindings': gtp_bindings,
        'seq_chunks': seq_chunks,
        'seq_chunk_size': chunk_size,
        'sequence_text': sequence_text,
        'sequence_fasta': sequence_fasta,
    }
    return render(request, 'proteindetail.html', context)





def Autocomplete(request):
    if 'term' in request.GET:
        query = request.GET.get('term')
        proteins = Protein.objects.filter(name__icontains=query)  # Adjust the query
        results = [{
            'label': protein.name,
            'url': reverse('protein:protein', args=[protein.name])  # Adjust as per your URL pattern
        } for protein in proteins]
        return JsonResponse(results, safe=False)
