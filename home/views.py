from django.shortcuts import render
from collections import defaultdict
import plotly.graph_objects as go
from .forms import DateForm
from protein.models import Protein
from structure.models import Structure, Rotamer
from common.models import ResiduePosition
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
    structure_chart = generate_structure_chart()
    species_pie_chart = generate_species_pie_chart()
    type_pie_chart = generate_structure_type_pie_chart()
    chemokine_type_pie = generate_chemokine_type_pie_chart()
    subfamily_pie = generate_chemokine_subfamily_pie_chart()
    structure_state_pie = generate_structure_state_pie_chart()
    ccn_position_barplot = generate_ccn_position_barplot()

    return render(request, 'home/charts.html', {
        'structure_chart': structure_chart,
        'species_pie_chart': species_pie_chart,
        'type_pie_chart': type_pie_chart,
        'chemokine_type_pie': chemokine_type_pie,
        'subfamily_pie': subfamily_pie,
        'structure_state_pie': structure_state_pie,
        'ccn_position_barplot': ccn_position_barplot,
    })



def generate_structure_chart():
    """Helper function to generate chart data for structures published per year and cumulative total."""
    structures = Structure.objects.all()
    year_counts = count_publications_by_year(structures)
    
    sorted_years = sorted(year_counts.keys())
    publication_counts = [year_counts[year] for year in sorted_years]
    cumulative_counts = [sum(publication_counts[:i + 1]) for i in range(len(publication_counts))]

    fig = go.Figure(data=[
        go.Bar(x=sorted_years, y=publication_counts, name='Structures Published', marker_color='blue'),
        go.Scatter(x=sorted_years, y=cumulative_counts, mode='lines+markers', name='Cumulative Structures',
                   line=dict(color='red', width=2), marker=dict(color='red', size=6))
    ])

    fig.update_layout(
        #title={'text': "Number of Chemokine Structures Published per Year and Cumulative Total", 'font_size': 24, 'xanchor': 'center', 'x': 0.5},
        xaxis_title='Year', yaxis_title='Number of Structures',
        legend=dict(x=0.01, y=0.99, bordercolor='Black', borderwidth=1)
    )

    return fig.to_html()

def generate_species_pie_chart():
    """Generates a pie chart showing the species distribution of proteins that have structures."""
    structures = Structure.objects.select_related('protein').all()

    species_counts = defaultdict(int)
    for structure in structures:
        if structure.protein and structure.protein.species:
            species_counts[structure.protein.species] += 1

    labels = list(species_counts.keys())
    values = list(species_counts.values())

    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=0.3, textinfo='none', showlegend=False)])
    fig.update_layout(title=None)

    return fig.to_html()


def count_publications_by_year(structures):
    """Helper function to count publications by year from a queryset of structures."""
    year_counts = defaultdict(int)
    for structure in structures:
        if structure.publication_date:
            year_counts[structure.publication_date.year] += 1
    return year_counts

def generate_structure_type_pie_chart():
    """Generates a pie chart showing the distribution of structure types (X-ray, NMR, Cryo-EM)."""
    structures = Structure.objects.select_related('structure_type').all()

    type_counts = defaultdict(int)
    for structure in structures:
        if structure.structure_type and structure.structure_type.name:
            type_counts[structure.structure_type.name] += 1

    labels = list(type_counts.keys())
    values = list(type_counts.values())

    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=0.3, textinfo='none', showlegend=False)])
    fig.update_layout(
        title=None
    )

    return fig.to_html()

def generate_chemokine_type_pie_chart():
    """Generates a pie chart showing distribution of chemokine types for proteins with structures."""
    structures = Structure.objects.select_related('protein').all()

    type_counts = defaultdict(int)
    for structure in structures:
        if structure.protein and structure.protein.type:
            type_counts[structure.protein.type] += 1

    labels = list(type_counts.keys())
    values = list(type_counts.values())

    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=0.3, textinfo='none', showlegend=False)])
    fig.update_layout(title=None)

    return fig.to_html()


def generate_chemokine_subfamily_pie_chart():
    """Generates a pie chart showing distribution of chemokine subfamilies for proteins with structures."""
    structures = Structure.objects.select_related('protein').all()

    subfamily_counts = defaultdict(int)
    for structure in structures:
        if structure.protein and structure.protein.subfamily:
            subfamily_counts[structure.protein.subfamily] += 1

    labels = list(subfamily_counts.keys())
    values = list(subfamily_counts.values())

    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=0.3, textinfo='none', showlegend=False)])
    fig.update_layout(title=None)

    return fig.to_html()

def generate_structure_state_pie_chart():
    """Generates a pie chart showing the distribution of structure states."""
    structures = Structure.objects.all()

    state_counts = defaultdict(int)
    for structure in structures:
        if structure.state:
            state_counts[structure.state] += 1

    labels = list(state_counts.keys())
    values = list(state_counts.values())

    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=0.3, textinfo='none', showlegend=False)])
    fig.update_layout(title=None)

    return fig.to_html()

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
    return fig.to_html()

