from django.shortcuts import render, get_object_or_404
from django.views.generic import TemplateView, DetailView, View
from .forms import ChemkinUpload, ReactionCondition
from .models import Mechanism
from .chemcheck import ChemError, CheckNegativeA, err_line_without_comment, check_collision_violation
from django.http import HttpResponseRedirect, Http404, request
from django.core.files.storage import FileSystemStorage
import os, re, sys, linecache
from django.core.files.base import File
from django.urls import reverse_lazy
from .ck2yaml import strip_nonascii
from canteradebugger.settings import MEDIA_ROOT, BASE_DIR
import cantera as ct
import numpy as np
from bokeh.plotting import figure, ColumnDataSource
from bokeh.embed import components
from bokeh.models import CustomJS, Slider
from bokeh.layouts import column, row, layout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.utils.decorators import method_decorator

# Create your views here.

class Home(TemplateView):
    template_name = 'home.html'
   
@login_required
def upload(request):
    # if username != request.user.username:
    #     raise Http404('You are not logged in with this account!')
    if request.method == 'POST':
        form = ChemkinUpload(request.POST, request.FILES, request.user)
        if form.is_valid:
            files = form.save(commit=False)
            files.user = request.user
            form.save()
            return HttpResponseRedirect(reverse_lazy('upload:list'))
    else:
        form = ChemkinUpload()
    return render(request, 'upload.html', {
        'form': form,
    })

@login_required
def ck2yaml(request, pk):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    if mechanism.user != request.user or request.user.is_authenticated == False:
        raise Http404('You do not have the access to this file!')
    conversion_log = "Going to try this...\n"

    from .ck2yaml import Parser
    import traceback
    import tempfile

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
        ct.Solution(out_name)
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
                    err_line, line_num = err_line_without_comment(error_path, line_num)
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
                                    elif value_error:
                                        suggestion += 'Suggestion: Here is expecting a number instead a string, \nYou can check the source to make sure the data is correct.\nThere could be an indentation error or missing E or unexpected character in that string which confused the system. \nPlease make sure you have got the indents and data format correctly in line {}.'.format(int(match.group(1)))
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
                suggestion += ('Suggestion: You may have two series of parameters for one reaction,'
                              '\nor the given parameters does not match to your reaction type'
                              '\nso the reaction type cannot be determined.'
                              '\nPlease check your data source to correct your data .')
        elif transport_error:
            suggestion += 'Suggestion: You can manually add the transport data for that species\nor delete the species from mechanism file \nor you can delete the transport file and do the conversion again.'
        mechanism.ct_conversion_errors = error_message
        mechanism.ct_mechanism_file = None
        mechanism.save()
        
    else:
        mechanism.ct_mechanism_file = out_name.split('media/')[-1]
        mechanism.save()
        conversion_log += f"\nConversion successful!\nCantera yaml file saved to {out_name}"
    return render(request, 'ck2yaml.html', {
    'mech': mechanism,
    'conversion_log': conversion_log,
    'suggestion': suggestion
    })

@login_required
def mechanism_detail(request, pk):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    if mechanism.user != request.user or request.user.is_authenticated == False:
        raise Http404('You do not have the access to this file!')
    mech_name = mechanism.ck_mechanism_file.name.split('/')[-1]
    therm_name = mechanism.ck_thermo_file.name.split('/')[-1]
    trans_name = mechanism.ck_transport_file.name.split('/')[-1]
    surf_name = mechanism.ck_surface_file.name.split('/')[-1]
    return render(request, 'upload/mechanism_detail.html', {
                  'object':mechanism,
                  'mech_name':mech_name,
                  'therm_name':therm_name,
                  'trans_name':trans_name,
                  'surf_name':surf_name})



@login_required   
def ace(request, pk, filetype):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    if mechanism.user != request.user or request.user.is_authenticated == False:
        raise Http404('You do not have the access to this file!')
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
        'filetype':filetype,
        })

@login_required
def mechanisms_list(request):
    user_id = request.user.id
    mechanisms = Mechanism.objects.filter(user=user_id)
    return render(request, 'list.html', {
       'mechanisms': mechanisms,
    })
    
class MechanismObjectMixin(object):
    model = Mechanism
    def get_object(self):
        pk = self.kwargs.get('pk')
        obj = None
        if pk is not None:
            obj = get_object_or_404(self.model, pk=pk)
        return obj

@method_decorator(login_required, name='dispatch')
class MechanismDeleteView(MechanismObjectMixin, View):
    template_name="file_delete_mechanism.html"
    def get(self, request, id=id, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj is not None:
            if obj.user != request.user or request.user.is_authenticated == False:
                raise Http404('You do not have the access to this file!')
            context['object'] = obj
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj.user != request.user or request.user.is_authenticated == False:
            raise Http404('You do not have the access to this file!')
        mech = obj.ck_mechanism_file
        if mech is not None:
            mech.delete()
            context['object'] = None
            username = request.user.username
            return HttpResponseRedirect(reverse_lazy('upload:list'))
        return render(request, self.template_name, context)

@method_decorator(login_required, name='dispatch')
class MechanismthermoDeleteView(MechanismObjectMixin, View):
    template_name="file_delete_thermo.html"
    def get(self, request, id=id, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj is not None:
            if obj.user != request.user or request.user.is_authenticated == False:
                raise Http404('You do not have the access to this file!')
            context['object'] = obj
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj.user != request.user or request.user.is_authenticated == False:
            raise Http404('You do not have the access to this file!')
        thermo = obj.ck_thermo_file
        if thermo is not None:
            thermo.delete()
            context['object'] = None
            username = request.user.username
            return HttpResponseRedirect(reverse_lazy('upload:list'))
        return render(request, self.template_name, context)

@method_decorator(login_required, name='dispatch')
class MechanismtransportDeleteView(MechanismObjectMixin, View):
    template_name="file_delete_transport.html"
    def get(self, request, id=id, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj is not None:
            if obj.user != request.user or request.user.is_authenticated == False:
                raise Http404('You do not have the access to this file!')
            context['object'] = obj
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj.user != request.user or request.user.is_authenticated == False:
            raise Http404('You do not have the access to this file!')
        transport = obj.ck_transport_file
        if transport is not None:
            transport.delete()
            context['object'] = None
            username = request.user.username
            return HttpResponseRedirect(reverse_lazy('upload:list'))
        return render(request, self.template_name, context)

@method_decorator(login_required, name='dispatch')
class MechanismsurfaceDeleteView(MechanismObjectMixin, View):
    template_name="file_delete_surface.html"
    def get(self, request, id=id, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj is not None:
            if obj.user != request.user or request.user.is_authenticated == False:
                raise Http404('You do not have the access to this file!')
            context['object'] = obj
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj.user != request.user or request.user.is_authenticated == False:
            raise Http404('You do not have the access to this file!')
        surface = obj.ck_surface_file
        if surface is not None:
            surface.delete()
            context['object'] = None
            uername = request.user.username
            return HttpResponseRedirect(reverse_lazy('upload:list'))
        return render(request, self.template_name, context)

@method_decorator(login_required, name='dispatch')
class MechanismUpdateView(MechanismObjectMixin, View):
    template_name="file_update.html"
    def get(self, request, id=id, *args, **kwargs):
        context = {}
        obj = self.get_object()
        if obj is not None:
            if obj.user != request.user or request.user.is_authenticated == False:
                raise Http404('You do not have the access to this file!')
            form = ChemkinUpload(instance=obj)
            context['object'] = obj
            context['form'] = form
        return render(request, self.template_name, context)

    def post(self, request, *args, **kwargs):
        obj = self.get_object()
        if obj.user != request.user or request.user.is_authenticated == False:
                raise Http404('You do not have the access to this file!')
        mech = obj.ck_mechanism_file
        thermo = obj.ck_thermo_file
        trans = obj.ck_transport_file
        surface = obj.ck_surface_file
        form = ChemkinUpload(request.POST, request.FILES, instance=obj)
        if form.is_valid():
            mech_update = form.cleaned_data['ck_mechanism_file']
            thermo_update = form.cleaned_data['ck_thermo_file']
            trans_update = form.cleaned_data['ck_transport_file']
            surface_update = form.cleaned_data['ck_surface_file']
            file_update_list = [mech_update, thermo_update, trans_update, surface_update]
            if mech.name:
                mech_path = mech.path
                if mech_update == False or mech_update != mech.name:
                    os.remove(mech_path)
            if thermo.name:
                thermo_path = thermo.path
                if thermo_update == False or thermo_update != thermo.name:
                    os.remove(thermo_path)
            if trans.name:
                trans_path = trans.path
                if trans_update == False or trans_update != trans.name:
                    os.remove(trans_path)
            if surface.name:
                surface_path = surface.path
                if surface_update == False or surface_update != surface.name:
                    os.remove(surface_path)
            form.save()
            return HttpResponseRedirect(reverse_lazy('upload:mechanism_detail', args=[obj.id]))

@login_required
def chemcheck(request, pk):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    if mechanism.user != request.user or request.user.is_authenticated == False:
        raise Http404('You do not have the access to this file!')
    path = mechanism.ct_mechanism_file.path
    name = 'cantera'
    species_plots = ChemError(path, name).check_continuity()
    if len(species_plots) == 0:
        return render(request, 'chemcheck.html')
    else: 
        l = layout(species_plots)
        script, div = components(l)
        return render(request, 'chemcheck1.html', {'script':script, 'div':div})

@login_required
def check_pdep_negative_A(request, pk):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    if mechanism.user != request.user or request.user.is_authenticated == False:
        raise Http404('You do not have the access to this file!')
    path = mechanism.ct_mechanism_file.path
    new_list_pdep = CheckNegativeA(path).new_arrhenius_dict()
    negative_A_reactions = CheckNegativeA(path).check_negative_A_factor(new_list_pdep)
    arr_sum_error_dict = {}
    T = [200, 500, 1000, 2000, 5000, 10000]
    for t in T:
        arr_sum_error = CheckNegativeA(path).check_sum_of_k(new_list_pdep, t)
        arr_sum_error_dict['{} K'.format(t)] = arr_sum_error
    return render(request, 'pdep_negative_A.html',{
        'negative_A_reactions':negative_A_reactions,
        'arr_sum_error_dict':arr_sum_error_dict,
    })

@login_required
def check_negative_dup_rxns_negative_A(request, pk):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    if mechanism.user != request.user or request.user.is_authenticated == False:
        raise Http404('You do not have the access to this file!')
    path = mechanism.ct_mechanism_file.path
    duplicate_reactions = CheckNegativeA(path).duplicate_reactions()
    pdep_duplicate_reactions = CheckNegativeA(path).duplicate_reactions_multi_P()
    dup_rxns_err_dict = {}
    pdep_dup_err_dict = {}
    T = [200, 500, 1000, 2000, 5000, 10000]
    for t in T:
        duplicate_reactions_err = CheckNegativeA(path).check_sum_of_k(duplicate_reactions, t)
        pdep_duplicate_reactions_err = CheckNegativeA(path).check_sum_of_k(pdep_duplicate_reactions, t)
        dup_rxns_err_dict['{} K'.format(t)] = duplicate_reactions_err
        pdep_dup_err_dict['{} K'.format(t)] = pdep_duplicate_reactions_err
    return render(request, 'dup_negative_A.html', {
        'dup_rxns_err_dict':dup_rxns_err_dict,
        'pdep_dup_err_dict':pdep_dup_err_dict
    })

@login_required
def reaction_condition(request, pk):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    if mechanism.user != request.user or request.user.is_authenticated == False:
        raise Http404('You do not have the access to this file!')
    if request.method == 'POST':
        form = ReactionCondition(request.POST, instance=mechanism)
        if form.is_valid:
            form.save()
            return HttpResponseRedirect('collision_violation')
    else:
        form = ReactionCondition()
    return render(request, 'reaction_condition.html', {
        'form': form
    })

@login_required
def collision_violation_check(request, pk):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    if mechanism.user != request.user or request.user.is_authenticated == False:
        raise Http404('You do not have the access to this file!')
    log = ''
    violators = 0
    T = mechanism.temperature
    P = mechanism.pressure
    path = mechanism.ct_mechanism_file.path
    try:
        violators = check_collision_violation(path,T, P).check_collision_violation()
    except Exception as e:
        error = "Unexpected error:" + str(e)
        print(e)
        log += error
    if violators == 0:
        return render(request, 'collision_violation.html', {
            'log': log
        })        
    else:
        return render(request, 'collision_violation.html', {
            'Temperature': T,
            'Pressure': P,
            'violators': violators,
            'mechanism': mechanism
        })

@login_required
def bokeh_chart(request, pk):
    mechanism = get_object_or_404(Mechanism, pk=pk)
    if mechanism.user != request.user or request.user.is_authenticated == False:
        raise Http404('You do not have the access to this file!')
    try:
        kf_ = []
        kr_ = []
        path = mechanism.ct_mechanism_file.path
        gas = ct.Solution(path)
        T = mechanism.temperature
        P = mechanism.pressure
        gas.TP = T, P
        kf = gas.forward_rate_constants
        kr = gas.reverse_rate_constants
        for i in range(len(kf)):
            if kf[i] > 0:
                kf_.append(kf[i])
        for i in range(len(kr)):
            if kr[i] > 0:
                kr_.append(kr[i])
        lgkf_min = np.log10(min(kf_))
        lgkr_min = np.log10(min(kr_))
        lgkf_max = np.log10(max(kf_))
        lgkr_max = np.log10(max(kr_))
        lgk_min = min([lgkf_min, lgkr_min, lgkf_max, lgkr_max])
        lgk_max = max([lgkf_min, lgkr_min, lgkf_max, lgkr_max])
        rxn_index = np.arange(len(kf)) + 1
        reactions = gas.reactions()
        for i in range(len(reactions)):
            reactions[i] = str(reactions[i])
        kf = np.log10(kf)
        kr = np.log10(kr)
        source = ColumnDataSource(data=dict(
            x = rxn_index,
            y1 = kf,
            y2 = kr,
            reaction = reactions,   
        ))
        TOOLTIPS = [
            ("index", "@x"),
            ("(x,y)", "($x, $y)"),
            ("reaction", "@reaction")
        ]
        plot = figure(plot_width=600, plot_height=600, tooltips=TOOLTIPS,
                    title="Reaction Rate Constants Check", x_axis_label="reaction index", 
                    y_axis_label="logrithm of reaction rate constant")
        plot.scatter("x", "y1", color= "blue", source=source, legend_label="forward reaction rate constant")
        plot.scatter("x", "y2", color = "red", source=source, legend_label="reverse reaction rate constant")
        plot.legend.location = "top_left"
        plot.legend.click_policy = "hide"

        k_slider = Slider(start=lgk_min, end=lgk_max, value=lgk_min, step=0.1, title="reaction rate constant")
        callback = CustomJS(args=dict(y_range=plot.y_range, source=source), code="""
            const data = source['data'];
            var start = cb_obj.value;
            y_range.setv({"start": start});
            const reactions = data['reaction'];
            const kf = data['y1'];
            const kr = data['y2'];
            const index = data['x'];
            var rxns_data = [];
            for (var i = 0; i < index.length; i++){
                if ( kf[i] >= start){
                    var rxn = {};
                    rxn["index"] = i + 1;
                    rxn["reaction"] = reactions[i];
                    if(rxns_data.includes(rxn) === false) {
                        rxns_data.push(rxn)
                    }                    
                } else if (kr[i] >= start){
                    var rxn = {};
                    rxn["index"] = i + 1;
                    rxn["reaction"] = reactions[i];
                    if(rxns_data.includes(rxn) === false) {
                        rxns_data.push(rxn)
                    }   
                };
            };
       //     document.getElementById("demo").innerHTML = JSON.stringify(rxns_data);
            document.getElementById("table").innerHTML = "";
            let table = document.querySelector("table");
            let thead = table.createTHead();
            let row = thead.insertRow();
            for (var key in rxns_data[0]) {
                let th = document.createElement("th");
                let text = document.createTextNode(key);
                th.appendChild(text);
                row.appendChild(th);
            }
            for (var element of rxns_data){
                let row_data = table.insertRow();
                for (key in element){
                    let cell = row_data.insertCell();
                    let text_data = document.createTextNode(element[key]);
                    cell.appendChild(text_data);
                }
            }
        """)
        k_slider.js_on_change('value', callback)
        layout = column(plot, column(k_slider))
        script, div = components(layout)
    except:
        error = "Unexpected error:" + str(sys.exc_info()[1])
        raise Http404("error when loading the yaml file by Cantera: {}".format(error))
    return render(request, 'bokeh_chart.html', {'script':script, 'div':div, 'mechanism':mechanism})
