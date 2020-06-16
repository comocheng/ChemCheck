from django.urls import path, register_converter, re_path
from upload import views

app_name='upload'

urlpatterns = [
    path('home/', views.Home.as_view(), name='home'),
    path('<username>/upload/', views.upload, name='upload'),
    path('<username>/mechanism/<int:pk>/', views.mechanism_detail, name='mechanism_detail'),
    path('<username>/ck2yaml/<int:pk>/', views.ck2yaml, name='ck2yaml'),
    path('<username>/mechanism/<int:pk>/ace/<str:filetype>', views.ace, name='ace-editor'),
    path('<str:username>/list/', views.mechanisms_list, name='list'),
    path('<username>/mechanism/<int:pk>/delete_mechanism', views.MechanismDeleteView.as_view(), name='delete_file'),
    path('<username>/mechanism/<int:pk>/delete_thermo', views.MechanismthermoDeleteView.as_view(), name='delete_file_thermo'),
    path('<username>/mechanism/<int:pk>/delete_transport', views.MechanismtransportDeleteView.as_view(), name='delete_file_transport'),
    path('<username>/mechanism/<int:pk>/delete_surface', views.MechanismsurfaceDeleteView.as_view(), name='delete_file_surface'),
    path('<username>/mechanism/<int:pk>/update', views.MechanismUpdateView.as_view(), name='update_file_transport'),
    path('<username>/mechanism/<int:pk>/chemcheck', views.chemcheck, name='chemcheck'),
    path('<username>/mechanism/<int:pk>/pdep_negativeA', views.check_pdep_negative_A, name='check_pdep_negative_A'),
    path('<username>/mechanism/<int:pk>/dup_negative_A', views.check_negative_dup_rxns_negative_A, name='check_negative_dup_rxns_negative_A'),
    path('<username>/mechanism/<int:pk>/reaction_condition', views.reaction_condition, name='reaction_condition'),
    path('<username>/mechanism/<int:pk>/collision_violation', views.collision_violation_check, name='collision_violation'),
    path('<username>/mechanism/<int:pk>/collision_violation/bokeh_chart', views.bokeh_chart, name='bokeh_chart'),
]
