from django.urls import path
from upload import views

urlpatterns = [
    path('home/', views.Home.as_view(), name='home'),
    path('upload/', views.upload, name='upload'),
    path('mechanism/<int:pk>/', views.MechanismDetailView.as_view(), name='mechanism-detail'),
    path('ck2yaml/<int:pk>/', views.ck2yaml, name='ck2yaml'),
    path('mechanism/<int:pk>/ace/', views.ace, name='ace-editor'),
    path('list/', views.mechanisms_list, name='list'),
    path('mechanism/<int:pk>/delete_mechanism', views.MechanismDeleteView.as_view(), name='delete_file'),
    path('mechanism/<int:pk>/delete_thermo', views.MechanismthermoDeleteView.as_view(), name='delete_file_thermo'),
    path('mechanism/<int:pk>/delete_transport', views.MechanismtransportDeleteView.as_view(), name='delete_file_transport'),
    path('mechanism/<int:pk>/delete_surface', views.MechanismsurfaceDeleteView.as_view(), name='delete_file_surface'),
    path('mechanism/<int:pk>/update', views.MechanismUpdateView.as_view(), name='update_file_transport'),
]