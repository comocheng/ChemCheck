from django.shortcuts import render, HttpResponseRedirect
from django.contrib.auth.forms import UserCreationForm, AuthenticationForm, PasswordChangeForm
from django.contrib.auth import login, authenticate, logout
from django.contrib.auth import update_session_auth_hash

# Create your views here.
def signup_view(request):
    if request.method=='POST':
        form = UserCreationForm(request.POST)
        if form.is_valid():
            form.save()
            user = form.save()
            login(request, user, backend='django.contrib.auth.backends.ModelBackend')
            return HttpResponseRedirect('/home/')
    else:
        form = UserCreationForm()
    return render(request, 'accounts/signup.html', {'form':form})

def logout_view(request):
    logout(request)
    return HttpResponseRedirect('/accounts/login/')
