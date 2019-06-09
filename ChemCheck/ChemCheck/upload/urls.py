from django.urls import path
from upload import views

urlpatterns = [
    path('home/', views.Home.as_view(), name='home'),
    path('upload/', views.upload, name='upload'),
    path('ace/', views.ace, name='ace-editor'),
    path('list/', views.upload_list, name='list'),
]