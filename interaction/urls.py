from django.urls import path
from django.conf import settings
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page

from interaction import views

app_name = 'interaction' 

urlpatterns = [
    path('binding-partner/<int:pk>/', views.ChemokineBindingPartnerDetailView.as_view(), name='chemokine_binding_partner_detail'),
    path('csv/<str:pdb_id>/<str:chain_id>/<int:binding_partner>/', views.csv_export, name='csv'),
    path('', views.SelectStructure, name='selectstructure'),
    path('<str:structure_id>/<str:chain_id>/', cache_page(60*60*24)(views.InteractionDetails.as_view()), name='interaction_details'),
    path('ajax/', views.ajax, name='ajax'),
    path('ifp_search/', views.IFPSearchResults.as_view(), name='ifp_search'),
    path('view_alignment/', views.ViewAlignment.as_view(), name='view_alignment'),
    path('umap-ifp/', views.UMAPIFPPlotView.as_view(), name='umap_ifp'),
    path('similarity-matrix/', views.IFPSimilarityMatrixView.as_view(), name='ifp_similarity_matrix'),
    path('binding-partner/<int:binding_partner_id>/pdb_download/', views.PDBDownload.as_view(), name='pdb_download'),
    path('ifp/similarity/antibody/', views.IFPSimilarityMatrixAntibodyView.as_view(), name='ifp_similarity_antibody'),
]