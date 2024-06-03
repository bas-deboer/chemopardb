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

from protein.models import Protein, ProteinSegment
from structure.models import Structure
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


@cache_page(60 * 60 * 24 * 7)
def protein(request, name):
    try:
        protein = Protein.objects.get(name=name)
    except Protein.DoesNotExist:
        return render(request, 'error.html')

    fullname = "CXCL12_HUMAN"
    structures = Structure.objects.filter(protein=protein)
    residues = Residue.objects.filter(protein=protein).order_by('sequence_number')
    generic_numbers = ResidueGenericNumber.objects.all()

    # Group residues by segment
    residues_by_segment = {}
    for residue in residues:
        segment = residue.segment
        if segment not in residues_by_segment:
            residues_by_segment[segment] = []
        residues_by_segment[segment].append(residue)

    context = {
        'protein': protein,
        'fullname': fullname,
        'structures': structures,
        'residues_by_segment': residues_by_segment,
        'generic_numbers': generic_numbers,
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