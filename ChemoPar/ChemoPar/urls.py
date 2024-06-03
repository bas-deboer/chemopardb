"""
URL configuration for ChemoPar project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
from django.conf.urls.static import static
from django.conf import settings

urlpatterns = [
    path('admin/', admin.site.urls),
    path('protein/', include('protein.urls')),
    path('structure/', include('structure.urls')),
    path('partner/', include('partner.urls')),
    path('model/', include('model.urls')),
    path('', include('home.urls')),
    path('index/', include('home.urls')),
    path('about/', include('home.urls')),
    path('search/', include('structure.urls')),
    path('browse/', include('home.urls')),
    path('partners/', include('home.urls')),
    path('documentation/', include('home.urls')),
    path('interaction/', include('interaction.urls')),
]

if not settings.PRODUCTION:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
    