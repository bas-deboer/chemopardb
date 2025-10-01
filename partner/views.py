from django.shortcuts import render, get_object_or_404
from django.http import HttpResponse, FileResponse
from django.template.loader import render_to_string
from django.views.generic import View, TemplateView
from django.db.models import Count, Q, Prefetch, TextField

from structure.models import Structure, PredictedComplex
from protein.models import Protein
from partner.models import Partner, PartnerEntity


def partner_list(request):
    selected_type = request.GET.get('type')
    if selected_type is not None:
        selected_type = selected_type.strip() or None

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

    context = {
        'partners': partner_data,
        'selected_type': selected_type,
    }
    return render(request, 'partner/partner_browser.html', context)

def partner_detail(request, partner_id):
    partner = get_object_or_404(Partner, pk=partner_id)

    structures_qs = (
        partner.structures.all()
        .select_related('protein', 'pdb_code', 'publication')
        .order_by('-publication_date', '-id')
    )
    structures = list(structures_qs)
    structure_count = len(structures)

    predicted_complexes_qs = (
        PredictedComplex.objects.filter(gpcr_partner=partner)
        .select_related('chemokine_protein')
        .order_by('-date_generated', '-id')
    )
    predicted_complexes = list(predicted_complexes_qs)
    predicted_complex_count = len(predicted_complexes)

    context = {
        'partner': partner,
        'structures': structures,
        'structure_count': structure_count,
        'predicted_complexes': predicted_complexes,
        'predicted_complex_count': predicted_complex_count,
    }
    return render(request, 'partner/partner_detail.html', context)
