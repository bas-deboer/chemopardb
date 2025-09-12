from django.urls import path
from structure.views import (
    StructureBrowser, StructureDetails, StructureInteractions, PDBDownload,
    PredictedStructureBrowser, PredictedStructureDetails,
    PredictedComplexBrowser, PredictedComplexDetails,
)
from structure import views
from django.conf import settings
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page

urlpatterns = [
    # Experimental structure pages
    path('', cache_page(60*60*24)(StructureBrowser.as_view()), name='structure_browser'),
    path('<int:structure_id>/', StructureDetails.as_view(), name='structure_details'),
    path('<str:pdb_id>/download_alignment/', views.download_alignment, name='download_alignment'),
    path('<int:structure_id>/pdb_download/', PDBDownload.as_view(), name='pdb_download'),

    # Predicted structure pages
    path('predicted/', PredictedStructureBrowser.as_view(), name='predicted_structure_browser'),
    path('predicted/<int:pk>/', PredictedStructureDetails.as_view(), name='predicted_structure_details'),

    # Predicted complex pages
    path('predicted_complexes/', PredictedComplexBrowser.as_view(), name='predicted_complex_browser'),
    path('predicted_complexes/<int:pk>/', PredictedComplexDetails.as_view(), name='predicted_complex_details'),
]