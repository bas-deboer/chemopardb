from django.shortcuts import render
from django.http import HttpResponse, FileResponse
from django.template.loader import render_to_string
from django.views.generic import View, TemplateView
from django.db.models import Count, Q, Prefetch, TextField

from structure.models import Structure
from protein.models import Protein
from partner.models import Partner

from io import StringIO, BytesIO
from Bio.PDB import PDBIO, PDBParser
import zipfile

from datetime import datetime
import requests

import subprocess
from random import randint
import os


class PartnerBrowser(TemplateView):
    """
    Fetching Partner data for browser
    """
    template_name = "partner/partner_browser.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        partners = Partner.objects.all()
        context['partners'] = partners
        return context

def partner_page(request):
    return HttpResponse("Hello, world. This is ChemoPar partner page.")

def partner(request, partner):
    try:
        partner = Partner.objects.get(name=partner)

        #chains = partner.chains
        #protein = partner.protein
        
        # Assuming structure_pdb contains the ID linking to the PDBData model
        #structure_pdb = Structure_PDB.objects.get(name=pdb_id)

        # Retrieve the PDB text content from the related PDBData instance
        #pdb_data = structure_pdb.pdbdata.pdb

        
    except Partner.DoesNotExist:
        return render(request, 'error.html')

    return render(request, 'partner/partner.html', {'partner' : partner,})
