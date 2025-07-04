from django.urls import path
from django.views.decorators.cache import cache_page

from home import views

urlpatterns = [
    path('', views.index, name='index'),
    path('index/', views.index, name='index'),
    path('about/', views.about, name='about'),
    path('example/', views.example, name='example'),
    path('charts/', views.charts, name='charts'),
]