"""canteradebugger URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.2/topics/http/urls/
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
from django.urls import path, include, reverse_lazy
from django.conf.urls.static import static
from django.conf import settings
from django.contrib.auth import views as auth_views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('upload.urls')),
    path('accounts/', include('accounts.urls')),
    
    path('accounts/login/', 
         auth_views.LoginView.as_view(template_name='accounts/login.html'), 
         name='login'),
    
    path('accounts/password_reset/', 
         auth_views.PasswordResetView.as_view(template_name='accounts/password_reset_form.html'), 
         name='password_reset'),

    path('accounts/password_reset/done', 
          auth_views.PasswordResetDoneView.as_view(template_name='accounts/password_reset_done.html'), 
          name='password_reset_done'),

    path('accounts/reset/<uid64>/<token>/', 
          auth_views.PasswordResetConfirmView.as_view(template_name='accounts/password_reset_form.html'), 
          name='password_reset_confirm'),
    
    path('accounts/reset_password_complete', 
         auth_views.PasswordResetCompleteView.as_view(template_name='accounts/password_reset_complete.html'),
         name='password_reset_complete'),
    
    path('accounts/password_change', 
          auth_views.PasswordChangeView.as_view(template_name='accounts/password_change_form.html'), 
          name='password_change'),
    
    path('accounts/password_change/done', 
          auth_views.PasswordChangeDoneView.as_view(template_name='accounts/password_change_done.html'), 
          name='password_change_done'),
    ]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
