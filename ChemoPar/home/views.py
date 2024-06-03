from django.shortcuts import render
from django.http import HttpResponse, FileResponse
from django.db.models import Count
from django.db.models.functions import TruncDate

import datetime
from collections import defaultdict
import plotly.graph_objects as go
import plotly.express as px
from .forms import DateForm
from collections import defaultdict
from datetime import datetime

from protein.models import Protein
from structure.models import Structure

def index(request):
    return render(request, 'home/home.html')

def about(request):
    structures = Structure.objects.all()

    # Manually process dates to count publications per year and cumulative sum
    year_counts = defaultdict(int)
    for structure in structures:
        if structure.publication_date:
            year = structure.publication_date.year
            year_counts[year] += 1

    # Sort the years
    sorted_years = sorted(year_counts.keys())
    publication_counts = [year_counts[year] for year in sorted_years]
    cumulative_counts = [sum(publication_counts[:i+1]) for i in range(len(publication_counts))]

    # Create the bar plot
    bar_trace = go.Bar(
        x=sorted_years,
        y=publication_counts,
        name='Structures Published',
        marker_color='blue'
    )

    # Create the line plot for cumulative counts
    line_trace = go.Scatter(
        x=sorted_years,
        y=cumulative_counts,
        mode='lines+markers',
        name='Cumulative Structures',
        line=dict(color='red', width=2),
        marker=dict(color='red', size=6)
    )

    # Combine the plots
    fig = go.Figure(data=[bar_trace, line_trace])

    # Update layout
    fig.update_layout(
        title={
            'text': "Number of Chemokine Structures Published per Year and Cumulative Total",
            'font_size': 24,
            'xanchor': 'center',
            'x': 0.5
        },
        xaxis_title='Year',
        yaxis_title='Number of Structures',
        legend=dict(
            x=0.01,
            y=0.99,
            bordercolor='Black',
            borderwidth=1
        )
    )

    chart = fig.to_html()
    context = {'chart': chart}

    return render(request, 'home/about.html', context)

def browse(request):

    headers = ["Family", "Gene name", "Species", "Uniprot"]
    
    data = Protein.objects.all()
    
    return render(request, 'home/browse.html', {'data' : data, 'headers' : headers})

def search(request):
    return render(request, 'search.html')

def partners(request):
    return render(request, 'home/partners.html')

def documentation(request):
    return render(request, 'home/documentation.html')

def contact(request):
    return render(request, 'home/contact.html')

def example(request):
    return render(request, 'home/example.html')