from django.urls import path
from partner import views
from partner.views import *
from django.conf import settings
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page


urlpatterns = [
    path('', views.partner_list, name='partner_list'),
    path('<int:partner_id>/info/', views.partner_detail, name='partner_detail'),
]