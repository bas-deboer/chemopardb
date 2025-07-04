from django.urls import path
from structure.views import StructureBrowser, StructureDetails, StructureInteractions, PDBDownload
from structure import views
from django.conf import settings
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page

#urlpatterns = [
#    url(r'^$', cache_page(60*60*24)(StructureBrowser.as_view()), name='structure_browser'),]

urlpatterns = [
    path('', cache_page(60*60*24)(StructureBrowser.as_view()), name='structure_browser'),
    #path('<str:pdb_id>/', cache_page(60*60*24)(StructureDetails.as_view()), name='structure'),
    path('<int:structure_id>/', StructureDetails.as_view(), name='structure_details'),
    path('<str:pdb_id>/download_alignment/', views.download_alignment, name='download_alignment'),
    #path('<int:structure_id>/<str:chain>', cache_page(60*60*24)(StructureInteractions.as_view()), name='structure_interactions'),
    path('<int:structure_id>/pdb_download/', PDBDownload.as_view(), name='pdb_download'),
]