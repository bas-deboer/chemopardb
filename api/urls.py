from django.urls import path, re_path
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page
from rest_framework.schemas import get_schema_view

from drf_spectacular.views import SpectacularAPIView, SpectacularSwaggerView
from rest_framework.routers import DefaultRouter

from api import views


urlpatterns = [
    # Documentation
    path('', SpectacularSwaggerView.as_view(url_name='schema'), name='swagger-ui'),
    path('schema/', SpectacularAPIView.as_view(), name='schema'),
    # Protein
    path('protein_details/', cache_page(3600*24*7)(views.ProteinDetailList.as_view()), name='protein-details'),
    # Structure
    path('structure/', cache_page(3600*24*7)(views.StructureDetailList.as_view()), name='structure-detail-list'),
    path('structure/search/', cache_page(3600*24*7)(views.StructureSearch.as_view()), name='structure-search'),
    path('structure/<int:pk>/', cache_page(3600*24*7)(views.StructureDetail.as_view()), name='structure-detail'),
    # Binding Partners
    path('binding-pairs/', views.ChemokineBindingPartnerList.as_view(), name='binding-pair-list'),
    # Partners
    path('partners_list', cache_page(3600*24*7)(views.PartnersList.as_view()), name='partners-list'),
    #path('partner__pair_list/', views.ChemokinePartnerPairList.as_view(), name='chemokine-partner-pair-list'),
    # Interactions
    path('interactions/', views.InteractionListView.as_view(), name='interaction-list'),
    path('interactions_IFPs/', views.ChemokinePartnerIFPList.as_view(), name='chemokine-partner-ifp-list'),
    # Autocomplete search bar
    path('protein-autocomplete/', views.protein_autocomplete, name='protein-autocomplete'),
    path('structure-autocomplete/', views.structure_autocomplete, name='structure-autocomplete'),
]

