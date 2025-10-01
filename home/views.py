from django.shortcuts import render
from datetime import timedelta
from collections import defaultdict
import plotly.graph_objects as go
from django.utils.html import escape
from django.utils import timezone
from django.db.models import Prefetch
import csv
from functools import lru_cache
from pathlib import Path
from django.conf import settings
from .forms import DateForm, ChemokineDiagramForm
from protein.models import Protein
from partner.models import Partner
from structure.models import Structure, Rotamer, PredictedComplex, ChemokineBindingPartner
from residue.models import Residue
from common.models import ResiduePosition
from common.diagrams_chemokine import DrawArrestinPlot
from interaction.models import (ChemokinePartnerInteraction,
                                PredictedChemokinePartnerInteraction,
                                )
from collections import Counter


def _build_homepage_stats():
    """Collect key ChemoPar-db metrics for the landing page."""
    now = timezone.now()
    quarter_ago = (now - timedelta(days=90)).date()

    structures_total = Structure.objects.count()
    structures_recent = Structure.objects.filter(
        publication_date__isnull=False,
        publication_date__gte=quarter_ago,
    ).count()

    partners_total = Partner.objects.count()
    partner_types = (
        Partner.objects.exclude(type__isnull=True)
        .exclude(type__exact="")
        .values_list("type", flat=True)
        .distinct()
        .count()
    )

    experimental_contacts = ChemokinePartnerInteraction.objects.count()
    predicted_contacts = PredictedChemokinePartnerInteraction.objects.count()
    interactions_total = experimental_contacts + predicted_contacts

    predicted_complex_total = PredictedComplex.objects.count()
    predicted_complex_recent = PredictedComplex.objects.filter(
        date_generated__isnull=False,
        date_generated__gte=quarter_ago,
    ).count()

    def format_delta(count, positive_label, fallback):
        if count > 0:
            return f"+{count:,} {positive_label}"
        return fallback

    return [
        {
            "label": "Total structures",
            "value": structures_total,
            "meta": format_delta(structures_recent, "in last 90 days", "Curated experimental entries"),
        },
        {
            "label": "Chemokine partners",
            "value": partners_total,
            "meta": (
                f"{partner_types:,} partner types"
                if partner_types
                else "Unique binding partners"
            ),
        },
        {
            "label": "Predicted complexes",
            "value": predicted_complex_total,
            "meta": format_delta(
                predicted_complex_recent,
                "updated last 90 days",
                "AlphaFold & docking models",
            ),
        },
    ]
    
def index(request):    
    stats = _build_homepage_stats()
    return render(request, 'home/home.html', {'homepage_stats': stats})

def about(request):
    chart = generate_structure_chart()
    return render(request, 'home/about.html', {'chart': chart})


def browse(request):
    headers = ["Family", "Gene name", "Species", "Uniprot"]
    data = Protein.objects.all()
    return render(request, 'home/browse.html', {'data': data, 'headers': headers})


def example(request):
    return render(request, 'home/example.html')


def charts(request):
    structures = list(
        Structure.objects.select_related("protein", "structure_type").all()
    )

    structure_chart = generate_structure_chart(structures)
    species_pie_chart = generate_species_pie_chart(structures)
    type_pie_chart = generate_structure_type_pie_chart(structures)
    chemokine_type_pie = generate_chemokine_type_pie_chart(structures)
    subfamily_pie = generate_chemokine_subfamily_pie_chart(structures)
    structure_state_pie = generate_structure_state_pie_chart(structures)
    partner_type_bar_chart = generate_partner_type_bar_chart(structures)
    ccn_position_barplot = generate_ccn_position_barplot()

    structure_count = len(structures)
    unique_species = len({
        s.protein.species
        for s in structures
        if getattr(s.protein, 'species', None)
    })
    unique_methods = len({
        s.structure_type.name
        for s in structures
        if getattr(s, 'structure_type', None) and s.structure_type.name
    })
    unique_states = len({
        s.state
        for s in structures
        if getattr(s, 'state', None)
    })
    latest_publication_date = None
    publication_dates = [
        s.publication_date
        for s in structures
        if getattr(s, 'publication_date', None)
    ]
    if publication_dates:
        latest_publication_date = max(publication_dates)

    return render(
        request,
        "home/charts.html",
        {
            "structure_chart": structure_chart,
            "species_pie_chart": species_pie_chart,
            "type_pie_chart": type_pie_chart,
            "chemokine_type_pie": chemokine_type_pie,
            "subfamily_pie": subfamily_pie,
            "structure_state_pie": structure_state_pie,
            "partner_type_bar_chart": partner_type_bar_chart,
            "ccn_position_barplot": ccn_position_barplot,
            "structure_count": structure_count,
            "unique_species": unique_species,
            "unique_methods": unique_methods,
            "unique_states": unique_states,
            "latest_publication_date": latest_publication_date,
        },
    )


SEGMENT_ORDER = [
    "N-term",
    "CX",
    "N-loop",
    "B1",
    "30s-loop",
    "B2",
    "40s-loop",
    "B3",
    "50s-loop",
    "Helix",
    "C-term",
]

SEGMENT_ALIAS_MAP = {
    "N-term": "N-term",
    "CX": "CX",
    "cxb1": "N-loop",
    "N-loop": "N-loop",
    "B1": "B1",
    "b1b2": "30s-loop",
    "30s-loop": "30s-loop",
    "B2": "B2",
    "b2b3": "40s-loop",
    "40s-loop": "40s-loop",
    "B3": "B3",
    "b3h": "50s-loop",
    "50s-loop": "50s-loop",
    "Helix": "Helix",
    "C-term": "C-term",
}

ALLOWED_SEQUENCE_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZ")


class SimpleResidue:
    __slots__ = ("sequence_number", "amino_acid", "segment", "ccn_number")

    def __init__(self, sequence_number, amino_acid, segment):
        self.sequence_number = sequence_number
        self.amino_acid = amino_acid
        self.segment = segment
        self.ccn_number = None


@lru_cache(maxsize=1)
def _load_protein_catalog():
    from django.db.models import Prefetch

    residues_qs = Residue.objects.order_by('sequence_number')
    proteins = (
        Protein.objects.exclude(sequence__isnull=True)
        .exclude(sequence__exact="")
        .prefetch_related(Prefetch('residue_set', queryset=residues_qs))
    )

    catalog = []
    for protein in proteins:
        raw_sequence = protein.sequence or ""
        sequence = "".join(ch for ch in raw_sequence.upper() if ch.isalpha())
        if not sequence:
            continue

        segment_map = {}
        for residue in protein.residue_set.all():
            if residue.sequence_number is None:
                continue
            seg = residue.segment
            if seg:
                segment_map[int(residue.sequence_number)] = SEGMENT_ALIAS_MAP.get(seg, seg)

        label = protein.uniprot_id or protein.gene_name or protein.name or f"Protein {protein.id}"
        catalog.append({
            "id": protein.id,
            "label": label,
            "sequence": sequence,
            "segments": segment_map,
            "length": len(sequence),
        })

    closest_pair = None
    best_identity = -1.0
    for i in range(len(catalog)):
        seq_i = catalog[i]["sequence"]
        for j in range(i + 1, len(catalog)):
            seq_j = catalog[j]["sequence"]
            aligned_i, aligned_j = _global_align(seq_i, seq_j)
            identity = _calculate_identity(aligned_i, aligned_j)
            if identity > best_identity:
                best_identity = identity
                closest_pair = {
                    "id1": catalog[i]["label"],
                    "id2": catalog[j]["label"],
                    "identity": round(identity * 100, 2),
                }

    return catalog, closest_pair


def _resolve_segment(segment_map, position, max_position):
    if not segment_map:
        return None

    pos = int(position)
    seg = segment_map.get(pos)
    if seg:
        return seg

    left = pos - 1
    right = pos + 1
    while left >= 1 or right <= max_position:
        if left >= 1:
            seg = segment_map.get(left)
            if seg:
                return seg
            left -= 1
        if right <= max_position:
            seg = segment_map.get(right)
            if seg:
                return seg
            right += 1
    return None


def _global_align(reference, query):
    if not reference or not query:
        return reference, query

    match = 2
    mismatch = -1
    gap = -2

    m = len(reference)
    n = len(query)

    dp = [[0] * (n + 1) for _ in range(m + 1)]
    trace = [[None] * (n + 1) for _ in range(m + 1)]

    for i in range(1, m + 1):
        dp[i][0] = i * gap
        trace[i][0] = "up"
    for j in range(1, n + 1):
        dp[0][j] = j * gap
        trace[0][j] = "left"

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag_score = dp[i - 1][j - 1] + (match if reference[i - 1] == query[j - 1] else mismatch)
            up_score = dp[i - 1][j] + gap
            left_score = dp[i][j - 1] + gap

            choices = [
                (diag_score, "diag", 0),
                (left_score, "left", 1),
                (up_score, "up", 2),
            ]
            best_score, best_move, _ = max(choices, key=lambda item: (item[0], -item[2]))
            dp[i][j] = best_score
            trace[i][j] = best_move

    aligned_reference = []
    aligned_query = []
    i, j = m, n
    while i > 0 or j > 0:
        direction = trace[i][j]
        if direction == "diag":
            aligned_reference.append(reference[i - 1])
            aligned_query.append(query[j - 1])
            i -= 1
            j -= 1
        elif direction == "up":
            aligned_reference.append(reference[i - 1])
            aligned_query.append("-")
            i -= 1
        else:  # "left" or None
            aligned_reference.append("-")
            aligned_query.append(query[j - 1])
            j -= 1

    aligned_reference.reverse()
    aligned_query.reverse()
    return "".join(aligned_reference), "".join(aligned_query)



def _calculate_identity(aligned_a, aligned_b):
    matches = 0
    aligned = 0
    for a, b in zip(aligned_a, aligned_b):
        if a == "-" and b == "-":
            continue
        if a != "-" and b != "-":
            aligned += 1
            if a == b:
                matches += 1
    return (matches / aligned) if aligned else 0.0


def _build_position_mapping(aligned_a, aligned_b):
    mapping = {}
    pos_a = 0
    pos_b = 0
    for char_a, char_b in zip(aligned_a, aligned_b):
        if char_a != "-":
            pos_a += 1
        if char_b != "-":
            pos_b += 1
        if char_a != "-" and char_b != "-":
            mapping[pos_b] = pos_a
    return mapping


def _combine_mappings(map_b_to_a, map_c_to_b):
    combined = {}
    for c_pos, b_pos in map_c_to_b.items():
        a_pos = map_b_to_a.get(b_pos)
        if a_pos:
            combined[c_pos] = a_pos
    return combined
def _assign_segments(sequence):
    catalog, closest_pair = _load_protein_catalog()

    if not sequence or not catalog:
        return {
            "ranges": [],
            "segment_string": "",
            "mapping": {},
            "segment_counts": Counter(),
            "alignment_rows": [],
            "reference_id": None,
            "used_reference": None,
            "mapping_reference": None,
            "canonical_reference": None,
            "pair_info": closest_pair,
            "pair_identity": None,
            "coverage": {"mapped": 0, "total": len(sequence or ""), "percent": 0.0},
        }

    best_entry = None
    best_alignment = None
    best_identity = -1.0

    for entry in catalog:
        aligned_candidate, aligned_user = _global_align(entry["sequence"], sequence)
        identity = _calculate_identity(aligned_candidate, aligned_user)
        if identity > best_identity:
            best_identity = identity
            best_entry = entry
            best_alignment = (aligned_candidate, aligned_user)

    if best_entry is None:
        return {
            "ranges": [],
            "segment_string": "",
            "mapping": {},
            "segment_counts": Counter(),
            "alignment_rows": [],
            "reference_id": None,
            "used_reference": None,
            "mapping_reference": None,
            "canonical_reference": None,
            "pair_info": closest_pair,
            "pair_identity": None,
            "coverage": {"mapped": 0, "total": len(sequence), "percent": 0.0},
        }

    aligned_candidate, aligned_user = best_alignment
    user_to_candidate = _build_position_mapping(aligned_candidate, aligned_user)

    segments_per_residue = []
    if user_to_candidate:
        first_mapped = min(user_to_candidate)
        last_mapped = max(user_to_candidate)
    else:
        first_mapped = 0
        last_mapped = 0

    max_candidate_pos = best_entry["length"]
    candidate_segments = best_entry["segments"]

    for pos in range(1, len(sequence) + 1):
        segment = None
        candidate_pos = user_to_candidate.get(pos)
        if candidate_pos:
            segment = _resolve_segment(candidate_segments, candidate_pos, max_candidate_pos)

        if not segment:
            if user_to_candidate:
                if pos < first_mapped:
                    anchor = user_to_candidate.get(first_mapped)
                    segment = _resolve_segment(candidate_segments, anchor, max_candidate_pos) if anchor else None
                    if not segment:
                        segment = "N-term"
                elif pos > last_mapped:
                    anchor = user_to_candidate.get(last_mapped)
                    segment = _resolve_segment(candidate_segments, anchor, max_candidate_pos) if anchor else None
                    if not segment:
                        segment = "C-term"
                else:
                    segment = segments_per_residue[-1] if segments_per_residue else "N-loop"
            else:
                segment = segments_per_residue[-1] if segments_per_residue else "N-term"
        segments_per_residue.append(segment or "N-term")

    if segments_per_residue:
        ranges = []
        current = segments_per_residue[0]
        start = 1
        for pos in range(2, len(segments_per_residue) + 1):
            seg = segments_per_residue[pos - 1]
            prev = segments_per_residue[pos - 2]
            if seg != prev:
                ranges.append((start, pos - 1, prev))
                start = pos
        ranges.append((start, len(segments_per_residue), segments_per_residue[-1]))
    else:
        ranges = []

    segment_counts = Counter()
    for start, end, seg in ranges:
        segment_counts[seg] += end - start + 1

    segment_string = ",".join(f"{start}-{end}:{seg}" for start, end, seg in ranges)

    mapped_count = len(user_to_candidate)
    coverage_percent = round((mapped_count / len(sequence)) * 100, 1) if sequence else 0.0

    display_reference = best_entry["label"]

    return {
        "ranges": ranges,
        "segment_string": segment_string,
        "mapping": user_to_candidate,
        "segment_counts": segment_counts,
        "alignment_rows": _format_alignment_rows(aligned_candidate, aligned_user),
        "reference_id": display_reference,
        "used_reference": display_reference,
        "mapping_reference": display_reference,
        "canonical_reference": display_reference,
        "pair_info": closest_pair,
        "pair_identity": round(best_identity * 100, 2) if best_identity >= 0 else None,
        "coverage": {"mapped": mapped_count, "total": len(sequence), "percent": coverage_percent},
    }

def _build_residue_objects(sequence, ranges):
    residues = []
    for start, end, segment in ranges:
        for position in range(start, end + 1):
            amino_acid = sequence[position - 1] if 0 <= position - 1 < len(sequence) else "X"
            residues.append(SimpleResidue(position, amino_acid, segment))
    return residues


def _format_alignment_rows(reference_aligned, query_aligned, width=60):
    rows = []
    for start in range(0, len(reference_aligned), width):
        ref_slice = reference_aligned[start:start + width]
        query_slice = query_aligned[start:start + width]
        match_line = "".join("|" if rc == qc and rc != "-" else " " for rc, qc in zip(ref_slice, query_slice))
        rows.append({
            "index": start + 1,
            "reference": ref_slice,
            "match": match_line,
            "query": query_slice,
        })
    return rows


def _summarize_segments(ranges, mapping):
    summary = []
    for start, end, segment in ranges:
        mapped_positions = [mapping.get(pos) for pos in range(start, end + 1) if mapping.get(pos)]
        if mapped_positions:
            reference_start = min(mapped_positions)
            reference_end = max(mapped_positions)
            coverage = round(len(mapped_positions) / (end - start + 1) * 100, 1)
        else:
            reference_start = None
            reference_end = None
            coverage = 0.0
        summary.append({
            "segment": segment,
            "start": start,
            "end": end,
            "length": end - start + 1,
            "reference_start": reference_start,
            "reference_end": reference_end,
            "coverage": coverage,
        })
    return summary


def _parse_segment_definition(segment_def, sequence_length):
    entries = []
    for part in segment_def.split(','):
        item = part.strip()
        if not item:
            continue
        if ':' not in item or '-' not in item:
            raise ValueError(f"Could not parse segment definition '{item}'. Use the format start-end:segment.")
        range_part, seg_name = item.split(':', 1)
        start_str, end_str = range_part.split('-', 1)
        try:
            start = int(start_str)
            end = int(end_str)
        except ValueError as exc:
            raise ValueError(f"Invalid residue numbers in '{item}'.") from exc
        if start < 1 or end < start or end > sequence_length:
            raise ValueError(f"Segment '{item}' is out of bounds for a sequence of length {sequence_length}.")
        normalized = SEGMENT_ALIAS_MAP.get(seg_name.strip(), seg_name.strip())
        entries.append((start, end, normalized))

    if not entries:
        return []

    entries.sort(key=lambda value: value[0])
    expected_start = 1
    for start, end, _ in entries:
        if start != expected_start:
            raise ValueError(f"Segments must cover the sequence without gaps. Expected residue {expected_start} but found {start}.")
        expected_start = end + 1
    if expected_start - 1 != sequence_length:
        raise ValueError(f"Segments must end at residue {sequence_length}.")

    return entries


def _sanitize_sequence(raw_sequence):
    if not raw_sequence:
        return ""
    return "".join(ch for ch in raw_sequence.upper() if ch.isalpha())


def chemokine_diagram(request):
    """Render an interactive chemokine diagram builder with automatic segmentation."""

    context = {}
    errors = []
    raw_sequence = request.POST.get("sequence", "") if request.method == "POST" else ""
    cleaned_sequence = _sanitize_sequence(raw_sequence)
    segment_def = request.POST.get("segments", "").strip() if request.method == "POST" else ""

    if request.method == "POST":
        if not cleaned_sequence:
            errors.append("Please provide a protein sequence containing alphabetic characters.")
        else:
            invalid_chars = sorted(set(cleaned_sequence) - ALLOWED_SEQUENCE_CHARS)
            if invalid_chars:
                errors.append("Unsupported characters in sequence: " + ", ".join(invalid_chars))

        assignment_result = None
        segment_ranges = []
        segment_string = segment_def

        if not errors and segment_def:
            try:
                segment_ranges = _parse_segment_definition(segment_def, len(cleaned_sequence))
            except ValueError as exc:
                errors.append(str(exc))
        elif not errors:
            assignment_result = _assign_segments(cleaned_sequence)
            if assignment_result['ranges']:
                segment_ranges = assignment_result['ranges']
                segment_string = assignment_result['segment_string']
            else:
                errors.append("Unable to align the submitted sequence to the reference chemokine alignment.")

        if not errors and segment_ranges:
            residues = _build_residue_objects(cleaned_sequence, segment_ranges)
            svg = str(DrawArrestinPlot(residues, "Custom Sequence", nobuttons='chemokine'))

            if assignment_result is None:
                assignment_result = _assign_segments(cleaned_sequence)

            segment_summary = _summarize_segments(segment_ranges, assignment_result['mapping'])

            context.update({
                'svg': svg,
                'sequence': cleaned_sequence,
                'segment_summary': segment_summary,
                'segment_string': segment_string,
                'segment_counts': assignment_result['segment_counts'],
                'alignment_rows': assignment_result['alignment_rows'],
                'reference_id': assignment_result['reference_id'],
                'used_reference': assignment_result['used_reference'],
                'mapping_reference': assignment_result['mapping_reference'],
                'canonical_reference': assignment_result['canonical_reference'],
                'closest_pair': assignment_result['pair_info'],
                'pair_identity': assignment_result['pair_identity'],
                'coverage': assignment_result['coverage'],
                'auto_generated': not bool(segment_def),
            })

    if errors:
        context['errors'] = errors

    if 'sequence' not in context:
        context['sequence'] = cleaned_sequence
        context['segment_string'] = segment_def

    if 'mapping_reference' not in context:
        context['mapping_reference'] = None
    if 'used_reference' not in context:
        context['used_reference'] = None
    if 'canonical_reference' not in context:
        context['canonical_reference'] = None
    if 'closest_pair' not in context:
        context['closest_pair'] = None

    context['form'] = ChemokineDiagramForm(initial={'sequence': context.get('sequence', '')})

    return render(request, "home/chemokine_diagram.html", context)


def generate_structure_chart(structures=None):
    """Return yearly and cumulative structure counts grouped by chemokine subfamily."""
    if structures is None:
        structures = Structure.objects.select_related("protein").all()

    year_subfamily_counts = defaultdict(lambda: defaultdict(int))
    subfamilies_present = set()

    for structure in structures:
        publication_date = getattr(structure, "publication_date", None)
        if not publication_date:
            continue

        protein = getattr(structure, "protein", None)
        raw_subfamily = getattr(protein, "subfamily", None) if protein else None
        if isinstance(raw_subfamily, str):
            cleaned = raw_subfamily.strip()
            subfamily = cleaned.upper() if cleaned else "Unassigned"
        else:
            subfamily = "Unassigned"

        year = publication_date.year
        year_subfamily_counts[year][subfamily] += 1
        subfamilies_present.add(subfamily)

    sorted_years = sorted(year_subfamily_counts.keys())
    preferred_order = [
        "CXCL",
        "CCL",
        "CX3CL",
        "XCL",
        "CXC",
        "CC",
        "CX3C",
        "XC",
        "Unassigned",
    ]
    remaining_subfamilies = sorted(subfamilies_present - set(preferred_order))
    sorted_subfamilies = [s for s in preferred_order if s in subfamilies_present] + [s for s in remaining_subfamilies if s not in preferred_order]

    annual_counts_by_subfamily = {subfamily: [] for subfamily in sorted_subfamilies}
    cumulative_counts_by_subfamily = {subfamily: [] for subfamily in sorted_subfamilies}
    running_totals = {subfamily: 0 for subfamily in sorted_subfamilies}

    for year in sorted_years:
        subfamily_counts = year_subfamily_counts[year]
        for subfamily in sorted_subfamilies:
            annual_count = subfamily_counts.get(subfamily, 0)
            annual_counts_by_subfamily[subfamily].append(annual_count)
            running_totals[subfamily] += annual_count
            cumulative_counts_by_subfamily[subfamily].append(running_totals[subfamily])

    return {
        "years": sorted_years,
        "subfamilies": sorted_subfamilies,
        "annual_counts_by_subfamily": annual_counts_by_subfamily,
        "cumulative_counts_by_subfamily": cumulative_counts_by_subfamily,
    }

def generate_species_pie_chart(structures=None):
    """Generates a pie chart showing the species distribution of proteins that have structures."""
    if structures is None:
        structures = Structure.objects.select_related('protein').all()

    species_counts = defaultdict(int)
    for structure in structures:
        if structure.protein and structure.protein.species:
            species_counts[structure.protein.species] += 1

    labels = list(species_counts.keys())
    values = list(species_counts.values())

    fig = go.Figure(
        data=[
            go.Pie(
                labels=labels,
                values=values,
                hole=0.3,
                textinfo="label+percent",
                textposition="inside",
            )
        ]
    )
    fig.update_layout(title=None, showlegend=True)

    html = fig.to_html(full_html=False, include_plotlyjs=False, div_id="species_pie")
    html += (
        "<script>"  # add click handler to filter browse page
        "document.getElementById('species_pie').on('plotly_click', function(data){"
        "var label=data.points[0].label;"
        "window.location='/protein?species=' + encodeURIComponent(label);"
        "});" "</script>"
    )
    return html


def count_publications_by_year(structures):
    """Helper function to count publications by year from a queryset of structures."""
    year_counts = defaultdict(int)
    for structure in structures:
        if structure.publication_date:
            year_counts[structure.publication_date.year] += 1
    return year_counts

def generate_structure_type_pie_chart(structures=None):
    """Return structure-method counts for Chart.js doughnut rendering."""
    if structures is None:
        structures = Structure.objects.select_related('structure_type').all()

    type_counts = defaultdict(int)
    for structure in structures:
        structure_type = getattr(structure, 'structure_type', None)
        name = getattr(structure_type, 'name', None)
        if name:
            type_counts[name] += 1

    labels = list(type_counts.keys())
    counts = [type_counts[label] for label in labels]

    return {
        "labels": labels,
        "counts": counts,
    }

def generate_chemokine_type_pie_chart(structures=None):
    """Generates a pie chart showing distribution of chemokine types for proteins with structures."""
    if structures is None:
        structures = Structure.objects.select_related('protein').all()

    type_counts = defaultdict(int)
    for structure in structures:
        if structure.protein and structure.protein.type:
            type_counts[structure.protein.type] += 1

    labels = list(type_counts.keys())
    values = list(type_counts.values())

    fig = go.Figure(
        data=[
            go.Pie(
                labels=labels,
                values=values,
                hole=0.3,
                textinfo="label+percent",
                textposition="inside",
            )
        ]
    )
    fig.update_layout(title=None, showlegend=True)

    html = fig.to_html(full_html=False, include_plotlyjs=False, div_id="chemokine_type_pie")
    html += (
        "<script>"
        "document.getElementById('chemokine_type_pie').on('plotly_click', function(data){"
        "var label=data.points[0].label;"
        "window.location='/protein?type=' + encodeURIComponent(label);"
        "});" "</script>"
    )
    return html


def generate_chemokine_subfamily_pie_chart(structures=None):
    """Generates a pie chart showing distribution of chemokine subfamilies for proteins with structures."""
    if structures is None:
        structures = Structure.objects.select_related('protein').all()

    subfamily_counts = defaultdict(int)
    for structure in structures:
        if structure.protein and structure.protein.subfamily:
            subfamily_counts[structure.protein.subfamily] += 1

    labels = list(subfamily_counts.keys())
    values = list(subfamily_counts.values())

    fig = go.Figure(
        data=[
            go.Pie(
                labels=labels,
                values=values,
                hole=0.3,
                textinfo="label+percent",
                textposition="inside",
            )
        ]
    )
    fig.update_layout(title=None, showlegend=True)

    html = fig.to_html(full_html=False, include_plotlyjs=False, div_id="subfamily_pie")
    html += (
        "<script>"
        "document.getElementById('subfamily_pie').on('plotly_click', function(data){"
        "var label=data.points[0].label;"
        "window.location='/protein?subfamily=' + encodeURIComponent(label);"
        "});" "</script>"
    )
    return html

def generate_structure_state_pie_chart(structures=None):
    """Return structure-state counts for Chart.js doughnut rendering."""
    if structures is None:
        structures = Structure.objects.all()

    state_counts = defaultdict(int)
    for structure in structures:
        state = getattr(structure, 'state', None)
        if state:
            state_counts[state] += 1

    labels = list(state_counts.keys())
    counts = [state_counts[label] for label in labels]

    return {
        "labels": labels,
        "counts": counts,
    }

def generate_partner_type_bar_chart(structures=None):
    """Return counts of structures grouped by binding partner type for Chart.js."""
    if structures is None:
        structures = Structure.objects.all()

    structure_ids = [structure.id for structure in structures if structure.id]
    if not structure_ids:
        return {"labels": [], "counts": []}

    type_to_structures = defaultdict(set)
    partner_entries = (
        ChemokineBindingPartner.objects.filter(structure_id__in=structure_ids)
        .exclude(partner_type__isnull=True)
        .exclude(partner_type__exact="")
        .values_list("partner_type", "structure_id")
    )

    for partner_type, structure_id in partner_entries:
        if not partner_type:
            continue
        normalized = partner_type.strip()
        if not normalized:
            continue
        type_to_structures[normalized].add(structure_id)

    if not type_to_structures:
        return {"labels": [], "counts": []}

    sorted_types = sorted(
        type_to_structures.keys(),
        key=lambda t: (-len(type_to_structures[t]), t.lower()),
    )
    counts = [len(type_to_structures[partner_type]) for partner_type in sorted_types]

    return {
        "labels": sorted_types,
        "counts": counts,
    }


def generate_ccn_position_barplot():
    """Generates a barplot showing how much each CCN position occurs across all chemokine chains (rotamers), ordered by canonical position."""

    # Get canonical CCN order from ResiduePosition
    canonical_positions = list(
        ResiduePosition.objects.order_by("position").values_list("ccn_number", flat=True)
    )

    # Count occurrences of each CCN number in Rotamer table
    ccn_qs = Rotamer.objects.exclude(ccn_number__isnull=True).exclude(ccn_number="").values_list("ccn_number", flat=True)
    ccn_counts = Counter(ccn_qs)

    # Order counts based on canonical order, fill missing with zero
    counts = [ccn_counts.get(ccn, 0) for ccn in canonical_positions]

    fig = go.Figure(
        data=[go.Bar(x=canonical_positions, y=counts, marker_color="indigo")]
    )
    fig.update_layout(
        title="Occurrence of Each CCN Number Position Across All Chemokine Chains",
        xaxis_title="CCN Number",
        yaxis_title="Occurrence Count",
        bargap=0.2,
        margin=dict(l=40, r=20, t=60, b=90),
        xaxis_tickangle=45,
    )
    return fig.to_html(full_html=False, include_plotlyjs=False, div_id="ccn_position_barplot")















