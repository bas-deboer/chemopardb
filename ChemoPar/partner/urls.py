from django.urls import path
from partner import views
from partner.views import *
from django.conf import settings
from django.views.generic import TemplateView
from django.views.decorators.cache import cache_page


urlpatterns = [
    path('', cache_page(60*60*24)(PartnerBrowser.as_view()), name='partner_browser'),
    path('<str:partner>/', views.partner, name='partner'),
]