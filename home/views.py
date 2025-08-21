from django.shortcuts import render
from collections import defaultdict
import plotly.graph_objects as go
from django.utils.html import escape

from .forms import DateForm, ChemokineDiagramForm
from protein.models import Protein
from structure.models import Structure, Rotamer
from common.models import ResiduePosition
from common.diagrams_chemokine import DrawArrestinPlot
from collections import Counter

def index(request):
    return render(request, 'home/home.html')


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
    ccn_position_barplot = generate_ccn_position_barplot()

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
            "ccn_position_barplot": ccn_position_barplot,
        },
    )


def chemokine_diagram(request):
    """Allow users to submit a sequence and then interactively define segments."""

    if request.method == "POST":
        seq = request.POST.get("sequence", "").replace("\n", "").strip().upper()
        segment_def = request.POST.get("segments", "").strip()

        if segment_def:
            class ResidueObj:
                def __init__(self, seq_num, aa, seg):
                    self.sequence_number = seq_num
                    self.amino_acid = aa
                    self.segment = seg
                    self.ccn_number = None

            residues = []
            try:
                for item in segment_def.split(","):
                    if not item.strip():
                        continue
                    range_part, seg_name = item.split(":")
                    start, end = [int(x) for x in range_part.split("-")]
                    seg_name = seg_name.strip()
                    for i in range(start, end + 1):
                        aa = seq[i - 1] if 0 <= i - 1 < len(seq) else "X"
                        residues.append(ResidueObj(i, aa, seg_name))
                svg = str(DrawArrestinPlot(residues, "Custom Sequence", nobuttons=True))
            except Exception as exc:
                svg = f"<p>Error generating diagram: {escape(exc)}</p>"
            return render(request, "home/chemokine_diagram.html", {"svg": svg})
        else:
            return render(request, "home/chemokine_diagram.html", {"seq": seq})

    form = ChemokineDiagramForm()
    return render(request, "home/chemokine_diagram.html", {"form": form})



def generate_structure_chart(structures=None):
    """Helper function to generate chart data for structures published per year and cumulative total."""
    if structures is None:
        structures = Structure.objects.all()
    year_counts = count_publications_by_year(structures)
    
    sorted_years = sorted(year_counts.keys())
    publication_counts = [year_counts[year] for year in sorted_years]
    cumulative_counts = [sum(publication_counts[:i + 1]) for i in range(len(publication_counts))]

    fig = go.Figure(
        data=[
            go.Bar(
                x=sorted_years,
                y=publication_counts,
                name="Structures Published",
                marker_color="blue",
            ),
            go.Scatter(
                x=sorted_years,
                y=cumulative_counts,
                mode="lines+markers",
                name="Cumulative Structures",
                line=dict(color="red", width=2),
                marker=dict(color="red", size=6),
            ),
        ]
    )

    fig.update_layout(
        xaxis_title="Year",
        yaxis_title="Number of Structures",
        legend=dict(x=0.01, y=0.99, bordercolor="Black", borderwidth=1),
    )

    return fig.to_html(full_html=False, include_plotlyjs=False, div_id="structure_chart")

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
    """Generates a pie chart showing the distribution of structure types (X-ray, NMR, Cryo-EM)."""
    if structures is None:
        structures = Structure.objects.select_related('structure_type').all()

    type_counts = defaultdict(int)
    for structure in structures:
        if structure.structure_type and structure.structure_type.name:
            type_counts[structure.structure_type.name] += 1

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

    html = fig.to_html(full_html=False, include_plotlyjs=False, div_id="structure_type_pie")
    html += (
        "<script>"
        "document.getElementById('structure_type_pie').on('plotly_click', function(data){"
        "var label=data.points[0].label;"
        "window.location='/structure?method=' + encodeURIComponent(label);"
        "});" "</script>"
    )
    return html

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
    """Generates a pie chart showing the distribution of structure states."""
    if structures is None:
        structures = Structure.objects.all()

    state_counts = defaultdict(int)
    for structure in structures:
        if structure.state:
            state_counts[structure.state] += 1

    labels = list(state_counts.keys())
    values = list(state_counts.values())

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

    html = fig.to_html(full_html=False, include_plotlyjs=False, div_id="structure_state_pie")
    html += (
        "<script>"
        "document.getElementById('structure_state_pie').on('plotly_click', function(data){"
        "var label=data.points[0].label;"
        "window.location='/structure?state=' + encodeURIComponent(label);"
        "});" "</script>"
    )
    return html

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

