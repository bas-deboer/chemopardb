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

from common.models import ResiduePosition
from protein.models import Protein, ProteinSegment, ProteinAlias
from structure.models import Structure, Rotamer
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
    alternate_names = ProteinAlias.objects.filter(protein=protein)

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

    context = {
        'protein': protein,
        'structures': structures,
        'aligned_sequences': aligned_sequences,
        'residues_by_segment': residues_by_segment,
        'alternate_names': alternate_names_str,
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
