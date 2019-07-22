from django.conf.urls import url
from accounts import views
from django.urls import path


app_name = 'accounts'

urlpatterns = [
    path('signup/', views.signup_view, name='signup'),
    path('logout/', views.logout_view, name='logout'),
]
