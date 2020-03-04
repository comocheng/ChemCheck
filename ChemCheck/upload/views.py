from django.shortcuts import render, get_object_or_404
from django.views.generic import TemplateView, DetailView, View
from .forms import ChemkinUpload
from .models import Mechanism
from .models import ChemError, cp_calculate, s_calculate, h_calculate
from django.http import HttpResponseRedirect, Http404
from django.core.files.storage import FileSystemStorage
import os
import re
from django.core.files.base import File
from django.urls import reverse_lazy
from .ck2yaml import strip_nonascii
import linecache
from canteradebugger.settings import MEDIA_ROOT

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

def ck2yaml(request, pk):
    mechanism = get_object_or_404(Mechanism, pk=pk)

    conversion_log = "Going to try this...\n"

    from .ck2yaml import Parser
    import traceback
    import tempfile
    from canteradebugger.settings import MEDIA_ROOT

    input_file = mechanism.ck_mechanism_file.path
    thermo_file =  mechanism.ck_thermo_file.path if mechanism.ck_thermo_file else None
    transport_file = mechanism.ck_transport_file.path if mechanism.ck_transport_file else None
    surface_file = mechanism.ck_surface_file.path if mechanism.ck_surface_file else None
    phase_name = 'gas' # will default to 'gas'
    out_name = os.path.join(os.path.split(input_file)[0], 'cantera.yaml')
    error_filename = os.path.join(os.path.split(input_file)[0], 'error.txt')
    with open(error_filename, 'w') as err_content:
        err_content.write('This is the error generated from Mechanism{0}\n'.format(mechanism.id))
    parser = Parser()
    suggestion = ''

    try:
        parser.convert_mech(input_file, 
                        thermo_file, 
                        transport_file,
                        surface_file,
                        phase_name,
                        out_name = out_name,
                        quiet = False,
                        permissive = True,
                        )
    except Exception as e:
        with open(error_filename, "r") as err:
            content = err.read()
        conversion_log += str(content)
        conversion_log += str(e)                      
        error_message = traceback.format_exc()
        conversion_log += error_message
        syn_error = re.search('Section starts with unrecognized keyword', str(e))
        ec_error = re.search("Error parsing elemental composition for "
                             "species (?P<name>...)", str(e))
        missing_end_number = re.search('Error while reading thermo entry starting on line (\d+)', conversion_log)
        value_error = re.search("ValueError: could not convert string to float: (?P<name>...)", conversion_log)
        transport_error = re.search("No transport data for species ", str(e))
        duplicate_reaction_type = re.search('Reaction entry contains parameters for more than one reaction type', str(e))
        
        match = re.search('Unable to parse .* near line (\d+):', content)
        if match:
            conversion_log += '\n\n'
            line_number = int(match.group(1))
            error_path = match.group(0).split("'")[1]
            with open(error_path, 'r', errors='ignore') as ck:
                lines = ck.readlines()
            context = 4
            excerpt = lines[ max(line_number-context,0):min(line_number+context, len(lines)) ]
            conversion_log += '\n'.join(excerpt)
            error_file_name = os.path.split(error_path)[1]
            if syn_error:
                suggestion += 'Suggestion: Please replace or delete any special character or redundant word in {0} line {1}'.format(error_file_name, line_number)
            
            elif missing_end_number:
                with open(error_path, 'r', errors='ignore'):
                    line_num = int(missing_end_number.group(1))
                    err_line = linecache.getline(error_path, line_num)
                    position = int(err_line.rfind('1', 74, 84))
                    if position == 79:
                        line_num += 1
                        err_line = linecache.getline(error_path, line_num)
                        position = int(err_line.rfind('2', 74, 84))
                        if position == 79:
                            line_num += 1
                            err_line = linecache.getline(error_path, line_num)
                            position = int(err_line.rfind('3', 74, 84))
                            if position == 79:
                                line_num += 1
                                err_line = linecache.getline(error_path, line_num)
                                position = int(err_line.rfind('4', 74, 84))
                                if position == 79:
                                    if ec_error:
                                        species = str(e).split()[-1]
                                        suggestion += 'Suggestion: Please make sure there is no indent error and typo in the error species {0} data in \n{1}(You can do this by comparing the error species with other species in the file).\nYou can also delete the data of species {0} and manually add them into converted file'.format(species, error_file_name)
                                    #suggestion += 'Suggestion: Please make sure your NASA data are neatly aligned in the 5 columns \n and the format of your first line is correct'
                                elif position == -1:
                                    suggestion += 'Suggestion: You are missing the index number 4 at the end of the line {}!'.format(line_num)
                                elif position != 79 and position != -1:
                                    suggestion += 'Suggestion: The index number 4 at the end of line {} is not in the same column with other lines'.format(line_num)
                                # else:
                                #     suggestion += 'Suggestion: Please make sure your NASA data are neatly aligned in the 5 columns \n and the format of your first line is correct'
                            elif position == -1:
                                suggestion += 'Suggestion: You are missing the index number 3 at the end of the line {}!'.format(line_num)
                            elif position != 79 and position != -1:
                                suggestion += 'Suggestion: The index number 3 at the end of line {} is not in the same column with other lines'.format(line_num)
                            # else:
                            #     suggestion += 'Suggestion: Please make sure your NASA data are neatly aligned in the 5 columns \n and the format of your first line is correct'
                        elif position == -1:
                            suggestion += 'Suggestion: You are missing the index number 2 at the end of the line {}!'.format(line_num)
                        elif position != 79 and position != -1:
                            suggestion += 'Suggestion: The index number 2 at the end of line {} is not in the same column with other lines'.format(line_num)
                        # else:
                        #     suggestion += 'Suggestion: Please make sure your NASA data are neatly aligned in the 5 columns \n and the format of your first line is correct'
                    elif position == -1:
                        suggestion += 'Suggestion: You are missing the index number 1 at the end of the line {}!'.format(line_num)
                    elif position != 79 and position != -1:
                        suggestion += 'Suggestion: The index number 1 at the end of line {} is not in the same column with other lines'.format(line_num)
                    # else:
                    #     suggestion += 'Suggestion: Please make sure your NASA data are neatly aligned in the 5 columns \n and the format of your first line is correct'
            elif ec_error:
                species = str(e).split()[-1]
                suggestion += 'Suggestion: Please make sure there is no indent error and typo in the error species {0} data in \n{1}(You can do this by comparing the error species with other species in the file).\nYou can also delete the data of species {0} and manually add them into converted file'.format(species, error_file_name)
            elif value_error:
                suggestion += 'Suggestion: Here is expecting a number instead a string, \nYou can check the source to make sure the data is correct.\nThere could be an indentation error or missing E or unexpected character in that string which confused the system. \nPlease make sure you have got the indents and data format correctly in line {}.'.format(int(match.group(1)))            
            elif duplicate_reaction_type:
                suggestion += 'Suggestion: You may have two set of parameters for one reaction \n Try to delete the duplicate parameters and convert again .'
        elif transport_error:
            suggestion += 'Suggestion: You can manually add the transport data for that species \nor you can delete the transport file and do the conversion again.'
        mechanism.ct_conversion_errors = error_message
        mechanism.ct_mechanism_file = None
        mechanism.save()
        
    else:
        mechanism.ct_mechanism_file = out_name
        mechanism.save()
        conversion_log += f"\nConversion successful!\nCantera yaml file saved to {out_name}"
    return render(request, 'ck2yaml.html', {
       'mech': mechanism,
       'conversion_log': conversion_log,
       'suggestion': suggestion
    })

    
class MechanismDetailView(DetailView):
    model = Mechanism
    
def ace(request, pk, filetype):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    
    try:
        f = {
            'mechanism': mechanism.ck_mechanism_file,
            'thermo': mechanism.ck_thermo_file,
            'transport': mechanism.ck_transport_file,
            'surface': mechanism.ck_surface_file,
        }[filetype]
    except KeyError:
        raise Http404("Unknown file type {}.".format(filetype))

    if not f:
        raise Http404("No {} file in this model".format(filetype))
    f = f.path

    filename = os.path.split(f)[-1]
    with open(f, 'r', errors='ignore') as file_content:
        content = strip_nonascii(file_content.read())
    return render(request, 'ace.html', {
        'content': content,
        'mechanism': mechanism,
        'filename': filename,
        })

def mechanisms_list(request):
    mechanisms = Mechanism.objects.all()
    return render(request, 'list.html', {
       'mechanisms': mechanisms
    })
    
class MechanismObjectMixin(object):
    model = Mechanism
    def get_object(self):
        pk = self.kwargs.get('pk')
        obj = None
        if pk is not None:
            obj = get_object_or_404(self.model, pk=pk)
        return obj

class MechanismDeleteView(MechanismObjectMixin, View):
    template_name="file_delete_mechanism.html"
    def get(self, request, id=id, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj is not None:
            context['object'] = obj
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        context = {}
        obj = self.get_object().ck_mechanism_file
        if obj is not None:
            obj.delete()
            context['object'] = None
            return HttpResponseRedirect('/list/')
        return render(request, self.template_name, context)

class MechanismthermoDeleteView(MechanismObjectMixin, View):
    template_name="file_delete_thermo.html"
    def get(self, request, id=id, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj is not None:
            context['object'] = obj
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        context = {}
        obj = self.get_object().ck_thermo_file
        if obj is not None:
            obj.delete()
            context['object'] = None
            return HttpResponseRedirect('/list/')
        return render(request, self.template_name, context)


class MechanismtransportDeleteView(MechanismObjectMixin, View):
    template_name="file_delete_transport.html"
    def get(self, request, id=id, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj is not None:
            context['object'] = obj
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        context = {}
        obj = self.get_object().ck_transport_file
        if obj is not None:
            obj.delete()
            context['object'] = None
            return HttpResponseRedirect('/list/')
        return render(request, self.template_name, context)

class MechanismsurfaceDeleteView(MechanismObjectMixin, View):
    template_name="file_delete_surface.html"
    def get(self, request, id=id, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj is not None:
            context['object'] = obj
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        context = {}
        obj = self.get_object().ck_surface_file
        if obj is not None:
            obj.delete()
            context['object'] = None
            return HttpResponseRedirect('/list/')
        return render(request, self.template_name, context)


class MechanismUpdateView(MechanismObjectMixin, View):
    template_name="file_update.html"
    def get(self, request, id=id, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj is not None:
            form = ChemkinUpload(instance=obj)
            context['object'] = obj
            context['form'] = form
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        obj = self.get_object()
        form = ChemkinUpload(request.POST, request.FILES, instance=obj)
        if form.is_valid():
            form.save()
            url = reverse_lazy('mechanism-detail', args=[obj.pk])
            return HttpResponseRedirect(url)
def chemcheck(request, pk):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    path = os.path.join(MEDIA_ROOT,'uploads/',str(mechanism.pk),'cantera.yaml')
    name = 'cantera'
    species_name = ChemError(path, name).check_continuity()
    return render(request, 'chemcheck.html', {
        'species_name':species_name
    })
