from django.conf.urls import url

from . import views


urlpatterns = [
    url(r'^$', views.basic, name='basic'),
    url(r'^submit', views.submit)
]
