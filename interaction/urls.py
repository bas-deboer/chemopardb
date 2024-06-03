from django.urls import path
from django.conf import settings
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page

from interaction import views


app_name = 'interaction' 

urlpatterns = [
    path('', views.SelectStructure, name='selectstructure'),
    path('<str:pdb_id>/<str:chain_id>/', cache_page(60*60*24)(views.InteractionDetails.as_view()), name='interaction_details'),
    path('excel/<slug>/', views.excel, name='download_excel'),
    path('ajax/', views.ajax, name='ajax'),
    path('ifp_search/', views.IFPSearchResults.as_view(), name='ifp_search'),
    path('view_alignment/', views.ViewAlignment.as_view(), name='view_alignment'),
]