from django.shortcuts import render
from django.http import HttpResponse, FileResponse
from django.template.loader import render_to_string
from django.views.generic import View

from structure.models import Structure
from protein.models import Protein
from model.models import Model

from io import StringIO, BytesIO
from Bio.PDB import PDBIO, PDBParser
import zipfile

from datetime import datetime

import subprocess
from random import randint
import os

def index(request):
    return HttpResponse("Hello, this is the model page.")


def model(request, name):
    try:
        model = Model.objects.get(name=name)

        pdb_id = model.pdb_id
        
        # Assuming structure_pdb contains the ID linking to the PDBData model
        #structure_pdb = Structure_PDB.objects.get(name=name)

        # Retrieve the PDB text content from the related PDBData instance
        #pdb_data = structure_pdb.pdbdata.pdb

        
    except Structure.DoesNotExist:
        return render(request, 'error.html')

    return render(request, 'model.html', {'model' : model,
                                            'pdb' : pdb_id,
                                            'name': name,
                                            'checkbox': True,})