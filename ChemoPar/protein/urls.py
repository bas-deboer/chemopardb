from django.urls import path
from protein.views import ProteinBrowser
from protein import views
from django.conf import settings
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page

app_name = 'protein' 

urlpatterns = [
    path('', cache_page(60*60*24)(ProteinBrowser.as_view()), name='browse'),
    path('<str:name>/', views.protein, name='protein'),
    path('autocomplete', views.Autocomplete, name='autocomplete'),
]