from django.shortcuts import render
from django.shortcuts import render, get_object_or_404
from django.http import HttpResponse, FileResponse
from django.template.loader import render_to_string
from django.views.generic import View, TemplateView
from django.db.models import Count, Q, Prefetch, TextField

from structure.models import Structure
from protein.models import Protein
from partner.models import Partner, PartnerEntity


def partner_list(request):
    partners = Partner.objects.prefetch_related('structures__entity_set__entity_type')
    partner_data = []

    for partner in partners:
        entity_types = set()
        for structure in partner.structures.all():
            for entity in structure.entity_set.all():
                if entity.entity_type:
                    entity_types.add(entity.entity_type.name)
        partner_data.append({
            'id': partner.id,
            'name': partner.name,
            'type': partner.type,
        })

    context = {'partners': partner_data}
    return render(request, 'partner/partner_browser.html', context)

def partner_detail(request, partner_id):
    partner = get_object_or_404(Partner, pk=partner_id)
    structures = partner.structures.all().select_related('protein')
    
    context = {
        'partner': partner,
        'structures': structures,
    }
    return render(request, 'partner/partner_detail.html', context)