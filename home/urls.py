from django.urls import path
from django.views.decorators.cache import cache_page

from home import views

urlpatterns = [
    path('', views.index, name='index'),
    path('index/', views.index, name='index'),
    path('about/', views.about, name='about'),
    path('search/', views.search, name='search'),
    path('browse/', views.browse, name='browse'),
    path('partners/', views.partners, name='partners'),
    path('documentation/', views.documentation, name='documentation'),
    path('contact/', views.contact, name='contact'),
    path('example/', views.example, name='example'),
]