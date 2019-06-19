from django.shortcuts import render
from django.views.generic import TemplateView, DetailView
from .forms import ChemkinUpload
from .models import Mechanism
from django.http import HttpResponseRedirect
from django.core.files.storage import FileSystemStorage

# Create your views here.

class Home(TemplateView):
    template_name = 'home.html'
   

def upload(request):
    if request.method == 'POST':
        form = ChemkinUpload(request.POST, request.FILES)
        if form.is_valid:
            form.save()
            return HttpResponseRedirect('/list/')
    else:
        form = ChemkinUpload()
    return render(request, 'upload.html', {
        'form': form
    })


def mechanisms_list(request):
    mechanisms = Mechanism.objects.all()
    return render(request, 'list.html', {
       'mechanisms': mechanisms
    })
    
class MechanismDetailView(DetailView):
    model = Mechanism

def ace(request):
    return render(request, 'ace.html', {})
