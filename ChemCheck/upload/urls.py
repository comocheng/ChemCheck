from django.urls import path
from upload import views

urlpatterns = [
    path('home/', views.Home.as_view(), name='home'),
    path('upload/', views.upload, name='upload'),
    path('mechanism/<int:pk>/', views.MechanismDetailView.as_view(), name='mechanism-detail'),
    path('ck2yaml/<int:pk>/', views.ck2yaml, name='ck2yaml'),
    path('ace/', views.ace, name='ace-editor'),
    path('list/', views.mechanisms_list, name='list'),
]