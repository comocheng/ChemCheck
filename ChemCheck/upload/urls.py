from django.urls import path, register_converter, re_path
from upload import views

app_name='upload'

urlpatterns = [
    path('home/', views.Home.as_view(), name='home'),
    path('upload/', views.upload, name='upload'),
    path('mechanism/<int:pk>/', views.mechanism_detail, name='mechanism_detail'),
    path('mechanism/ck2yaml/<int:pk>/', views.ck2yaml, name='ck2yaml'),
    path('list/', views.mechanisms_list, name='list'),
    path('mechanism/<int:pk>/ace/<str:filetype>', views.ace, name='ace_editor'),
    path('mechanism/<int:pk>/delete_thermo', views.MechanismthermoDeleteView.as_view(), name='delete_file_thermo'),
    path('mechanism/<int:pk>/delete_mechanism', views.MechanismDeleteView.as_view(), name='delete_file'),
    path('mechanism/<int:pk>/delete_transport', views.MechanismtransportDeleteView.as_view(), name='delete_file_transport'),
    path('mechanism/<int:pk>/delete_surface', views.MechanismsurfaceDeleteView.as_view(), name='delete_file_surface'),
    path('mechanism/<int:pk>/update', views.MechanismUpdateView.as_view(), name='update_file'),
    path('mechanism/<int:pk>/chemcheck', views.chemcheck, name='chemcheck'),
    path('mechanism/<int:pk>/dup_negative_A', views.check_negative_dup_rxns_negative_A, name='check_negative_dup_rxns_negative_A'),
    path('mechanism/<int:pk>/pdep_negativeA', views.check_pdep_negative_A, name='check_pdep_negative_A'),
    path('mechanism/<int:pk>/reaction_condition', views.reaction_condition, name='reaction_condition'),
    path('mechanism/<int:pk>/collision_violation', views.collision_violation_check, name='collision_violation'),
    path('mechanism/<int:pk>/collision_violation/bokeh_chart', views.bokeh_chart, name='bokeh_chart'),
]
