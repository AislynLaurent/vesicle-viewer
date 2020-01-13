## IMPORTS
# From django
from django.shortcuts import get_object_or_404, render
from django.shortcuts import redirect

# Other libraries
import json as js
from matplotlib import pyplot as plt
import mpld3
from datetime import datetime
import numpy as np
import lmfit as lsq
from copy import deepcopy

# Models
from .models import ExtendedUser
from .models import Lipid
from .models import Atom
from .models import Project
from .models import Project_Lipid
from .models import Symmetrical_Parameters
from .models import Data_Set
from .models import Data_Lipid
from .models import Data_Lipid_Atom

# Forms
from .forms import Project_Form
from .forms import Project_Lipid_Form
from .forms import Parameter_Form
from .forms import Parameter_Fit_Form
from .forms import Data_Upload_Form
from .forms import Data_Form
from .forms import Data_Lipid_Form
from .forms import Data_Lipid_Atom_Form

# Other imports
from .symfit import symmetrical_fit
from .symfit import symmetrical_graph

## STATIC PAGES
# Home
def index(request):
    context = {}
    return render(request, 'viewer/index.html', context)
# About
def about(request):
    context = {}
    return render(request, 'viewer/about.html', context)
# Help
def get_help(request):
    context = {}
    return render(request, 'viewer/help.html', context)
# Privacy Policy
def privacy(request):
    context = {}
    return render(request, 'viewer/privacy.html', context)

## MODEL PAGES

## DYNAMIC PAGES
## Project stuff
# Start a new project
def project_new(request):
    if request.method == 'POST':
        form = Project_Form(request.POST)
        if form.is_valid():
            project = form.save(commit=False)
            project.owner = request.user
            project.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Project_Form()
    return render(request, 'viewer/form.html', {'form': form})

# List existing projects
def project_list(request):
    symmetrical_projects = Project.objects.filter(owner=request.user, model_type='SM').order_by('project_title')
    return render(request, 'viewer/project_list.html', {'symmetrical_projects':symmetrical_projects})

# Specific project details
def project_detail(request, project_id):
    project = get_object_or_404(Project, id=project_id)
    project_lipids = Project_Lipid.objects.filter(project_title_id=project_id)
    parameters = Symmetrical_Parameters.objects.filter(project_title_id=project_id)
    datas = Data_Set.objects.filter(project_title_id=project_id)
    data_lipids = Data_Lipid.objects.filter(data_lipid_name__project_title__id=project_id).select_related('data_lipid_name')

    percentage = 0
    for lipid in project_lipids:
        percentage = percentage + lipid.lipid_mol_fraction

    return render(
        request,
        'viewer/project_detail.html', {
            'project':project,
            'project_lipids':project_lipids,
            'datas':datas,
            'data_lipids':data_lipids,
            'parameters':parameters,
            'total_lipid_percentage':percentage,
        }
    )

# Edit an existing project
def project_edit(request, project_id):
    project = get_object_or_404(Project, id=project_id)

    if request.method == 'POST':
        form = Project_Form(request.POST, instance=project)
        if form.is_valid():
            post = form.save(commit=False)
            post.owner = request.user
            post.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Project_Form(instance=project)

    return render(request, 'viewer/form.html', {'form': form})

## Project Lipids
# Add lipids to a project
def project_lipid_new(request, project_id):
    project = get_object_or_404(Project, id=project_id)

    if request.method == 'POST':
        form = Project_Lipid_Form(request.POST)
        if form.is_valid():
            lipids = form.save(commit=False)
            lipids.project_title = project
            lipids.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Project_Lipid_Form()

    return render(request, 'viewer/form.html', {'form': form})

# Edit lipid
def project_lipid_edit(request, project_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    lipid = get_object_or_404(Project_Lipid, id=lipid_id)

    if request.method == 'POST':
        form = Project_Lipid_Form(request.POST, instance=lipid)
        if form.is_valid():
            lipids = form.save(commit=False)
            lipids.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Project_Lipid_Form(instance=lipid)

    return render(request, 'viewer/form.html', {'form': form})

# Delete project lipid
def project_lipid_delete_warning(request, project_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    lipid = get_object_or_404(Project_Lipid, id=lipid_id)
    
    if request.method == 'POST':
        lipid.delete()
        return redirect('viewer:project_detail', project_id=project.id)

    return render(request, 'viewer/delete_warning.html', {'project':project})

## Parameters
# Add parameters to a project
def parameters_new(request, project_id):
    project = get_object_or_404(Project, id=project_id)
    project_lipids = Project_Lipid.objects.filter(project_title_id=project_id).select_related('project_lipid_name')

    combined_head_volume = 0
    combined_tail_volume = 0

    x = project.system_tempurature

    for lipid in project_lipids:
        combined_head_volume = combined_head_volume + (
                eval(lipid.project_lipid_name.hg_volume_equation)
                * lipid.lipid_mol_fraction
            )
        combined_tail_volume = combined_tail_volume + (
                eval(lipid.project_lipid_name.tg_volume_equation)
                * lipid.lipid_mol_fraction
            )

    if request.method == 'POST':
        form = Parameter_Form(request.POST)
        if form.is_valid():
            parameters = form.save(commit=False)

            # overall
            parameters.project_title = project

            # Math'ed values
            parameters.headgroup_volume = combined_head_volume
            parameters.headgroup_volume_upperbound = combined_head_volume*1.5
            parameters.headgroup_volume_lowerbound = -(combined_head_volume*1.5)

            parameters.chain_volume = combined_tail_volume
            parameters.chain_volume_upperbound = combined_tail_volume*1.5
            parameters.chain_volume_lowerbound = -(combined_tail_volume*1.5)

            # Upper & lowerbounds
            bilayer_thickness = form.cleaned_data['bilayer_thickness']
            parameters.bilayer_thickness_upperbound = bilayer_thickness*1.5
            parameters.bilayer_thickness_lowerbound = -(bilayer_thickness*1.5)

            terminal_methyl_volume = form.cleaned_data['terminal_methyl_volume']
            parameters.terminal_methyl_volume_upperbound = terminal_methyl_volume*1.5
            parameters.terminal_methyl_volume_lowerbound = -(terminal_methyl_volume*1.5)

            lipid_area = form.cleaned_data['lipid_area']
            parameters.lipid_area_upperbound = lipid_area*1.5
            parameters.lipid_area_lowerbound = -(lipid_area*1.5)

            headgroup_thickness = form.cleaned_data['headgroup_thickness']
            parameters.headgroup_thickness_upperbound = headgroup_thickness*1.5
            parameters.headgroup_thickness_lowerbound = -(headgroup_thickness*1.5)

            sigma = form.cleaned_data['sigma']
            parameters.sigma_upperbound = sigma*1.5
            parameters.sigma_lowerbound = -(sigma*1.5)

            parameters.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Parameter_Form()

    return render(request, 'viewer/form.html', {'form': form})

# Edit parameters
def parameters_edit(request, project_id, parameter_id):
    project = get_object_or_404(Project, id=project_id)
    parameters = get_object_or_404(Symmetrical_Parameters, id=parameter_id)

    if request.method == 'POST':
        form = Parameter_Form(request.POST, instance=parameters)
        if form.is_valid():
            parameters = form.save(commit=False)
            parameters.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Parameter_Form(instance=parameters)

    return render(request, 'viewer/form.html', {'form': form})

# Delete parameters
def parameter_delete_warning(request, project_id, parameter_id):
    project = get_object_or_404(Project, id=project_id)
    parameter = get_object_or_404(Symmetrical_Parameters, id=parameter_id)
    
    if request.method == 'POST':
        parameter.delete()
        return redirect('viewer:project_detail', project_id=project.id)

    return render(request, 'viewer/delete_warning.html', {'project':project})

## Data
# Input data
def data_upload(request, project_id):
    project = get_object_or_404(Project, id=project_id)

    q = []
    i = []
    e = []

    # Upload file
    if "data_upload" in request.POST:
        data_upload_form = Data_Upload_Form(request.POST, request.FILES)
        if data_upload_form.is_valid():
            data_info = data_upload_form.save(commit=False)
            data_info.project_title = project

            data_file = request.FILES["data_file"]

            # Incorrect file type?
            if not ( data_file.name.endswith('.txt') or data_file.name.endswith('.dat') ):
                print("incorrect type")
                # messages.error(request,'File is not CSV type')
                # return HttpResponseRedirect(reverse("myapp:upload_csv"))

            # File too big?
            if data_file.multiple_chunks():
                print("too big")
                # messages.error(request,"Uploaded file is too big (%.2f MB)." % (data_file.size/(1000*1000),))
                # return redirect('viewer:project_detail', project_id=project.id)

            file_data = data_file.read().decode("utf-8")

            lines = file_data.split("\n")

            #loop over the lines and save them in db. If error , store as string and then display
            for line in lines:

                line = line.strip()

                if not line[:1].isdigit():
                    pass

                else:
                    fields = line.split()

                    q.append(float(fields[0]))
                    i.append(float(fields[1]))
                    try:
                        e.append(float(fields[2]))
                    except IndexError:
                        e.append(1)

            data_info.q_value = q
            data_info.intensity_value = i
            data_info.error_value = e

            data_upload_form.save()

        return redirect('viewer:project_detail', project_id=project.id)
    else:
        data_upload_form = Data_Upload_Form()

    return render(
        request, 'viewer/data_upload.html',
        {
            'project':project,
            'data_upload_form':data_upload_form
        }
    )

# Edit data files
def data_edit(request, project_id, data_id):
    project = get_object_or_404(Project, id=project_id)
    data = get_object_or_404(Data_Set, id=data_id)

    if request.method == 'POST':
        form = Data_Form(request.POST, instance=data)
        if form.is_valid():
            data = form.save(commit=False)
            data.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Data_Form(instance=data)

    return render(request, 'viewer/form.html', {'form': form})

# Delete data
def data_delete_warning(request, project_id, data_id):

    project = get_object_or_404(Project, id=project_id)
    data = get_object_or_404(Data_Set, id=data_id)
    
    if request.method == 'POST':
        data.delete()
        return redirect('viewer:project_detail', project_id=project.id)

    return render(request, 'viewer/delete_warning.html', {'project':project})

## Data Lipids
# Add lipids adjustment
def data_lipid_new(request, project_id, data_id):
    project = get_object_or_404(Project, id=project_id)
    data = get_object_or_404(Data_Set, id=data_id)

    number_project_lipids = Project_Lipid.objects.filter(project_title_id=project_id).count()
    number_data_lipids = Data_Lipid.objects.filter(data_lipid_name__project_title__id=project_id).count()

    enough_adjustments = False

    if number_data_lipids == number_project_lipids:
        enough_adjustments = True

    if "lipid_info" in request.POST:
        lipid_info_form = Data_Lipid_Form(request.POST)

        if lipid_info_form.is_valid():
            # lipid_info = lipid_info_form.save(commit=False)
            # lipid_info.data_set_title = data
            # lipid_info.save()

            in_lipid_name = lipid_info_form.cleaned_data['data_lipid_name']
            in_lipid_suffix = lipid_info_form.cleaned_data['data_lipid_suffix']

            existing_lipid, created_lipid = Data_Lipid.objects.get_or_create(
                data_lipid_name=in_lipid_name,
                defaults={
                    'data_set_title':data,
                    'data_lipid_name':in_lipid_name,
                    'data_lipid_suffix':in_lipid_suffix
                }
            )

            if existing_lipid:
                return redirect('viewer:data_lipid_edit', project_id=project.id, data_id=data.id, lipid_id=existing_lipid.id)
            elif created_lipid:
                return redirect('viewer:data_lipid_edit', project_id=project.id, data_id=data.id, lipid_id=created_lipid.id)
    else:
        lipid_info_form = Data_Lipid_Form()

    return render(
        request,
        'viewer/data_lipid_adjustment.html', {
            'project':project,
            'data':data,
            'lipid_info_form':lipid_info_form,
            'enough_adjustments':enough_adjustments
        }
    )

# Edit lipid adjustments
def data_lipid_edit(request, project_id, data_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    data = get_object_or_404(Data_Set, id=data_id)
    data_lipid = get_object_or_404(Data_Lipid, id=lipid_id)
    data_lipid_atoms = Data_Lipid_Atom.objects.filter(data_lipid_name=data_lipid)

    if "lipid_info" in request.POST:
        lipid_info_form = Data_Lipid_Form(request.POST, instance=data_lipid)
    else:
        lipid_info_form = Data_Lipid_Form(instance=data_lipid)

    if "atom_ammount" in request.POST:
        atom_ammount_form = Data_Lipid_Atom_Form(request.POST)
        if atom_ammount_form.is_valid():

            in_atom_name = atom_ammount_form.cleaned_data['data_lipid_atom_name']
            in_atom_location = atom_ammount_form.cleaned_data['atom_location']
            in_atom_ammount = atom_ammount_form.cleaned_data['data_lipid_atom_ammount']

            existing_atom, created_atom = Data_Lipid_Atom.objects.update_or_create(
                data_lipid_atom_name=in_atom_name,
                atom_location=in_atom_location,
                defaults={
                    'data_lipid_name':data_lipid,
                    'data_lipid_atom_name':in_atom_name,
                    'data_lipid_atom_ammount':in_atom_ammount,
                    'atom_location':in_atom_location,
                }
            )
    else:
        atom_ammount_form = Data_Lipid_Atom_Form()

    if "done" in request.POST:
        return redirect('viewer:project_detail', project_id=project.id)

    return render(
        request, 
        'viewer/data_lipid_adjustment.html', {
            'data_lipid_atoms':data_lipid_atoms,
            'lipid_info_form':lipid_info_form,
            'atom_ammount_form':atom_ammount_form,
            'data':data
        }
    )

# Delete lipid adjustment
def data_lipid_delete_warning(request, project_id, data_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    data_lipid = get_object_or_404(Data_Lipid, id=lipid_id)
    
    if request.method == 'POST':
        data_lipid.delete()
        return redirect('viewer:project_detail', project_id=project.id)

    return render(request, 'viewer/delete_warning.html', {'project':project})

## Fit
# Main fit page
def fit_main(request, project_id, parameter_id):
    ## Import
    x_user = get_object_or_404(ExtendedUser, user=request.user)

    # Overall
    project = get_object_or_404(Project, id=project_id)
    parameter = get_object_or_404(Symmetrical_Parameters, id=parameter_id)
    project_lipids = Project_Lipid.objects.filter(project_title_id=project_id).select_related('project_lipid_name')

    # Data
    datas = Data_Set.objects.filter(project_title_id=project_id)
    xray_datas = datas.filter(data_type='XR')
    neutron_datas = datas.filter(data_type='NU')

    # Decalre
    now = datetime.now()
    fit_result = None
    show_statistics = False

    # Forms
    # Dismiss the tutorial
    if "dismiss" in request.POST:
        x_user.display_tutorial = False
        x_user.save()

    # Update parameters
    if "parameter_update" in request.POST:
        parameter_update_form = Parameter_Fit_Form(request.POST, instance=parameter)
        if parameter_update_form.is_valid():
            parameter = parameter_update_form.save(commit=False)
            parameter.save()
    else:
        parameter_update_form = Parameter_Fit_Form(instance=parameter)

    # Save parameters
    if "parameter_save" in request.POST:
        # Copy current instance
        new_parameter = deepcopy(parameter)

        # Set new values
        new_parameter.description = now.strftime("%m/%d/%H:%M")
        new_parameter.id = None

        # Save
        new_parameter.save()

        return redirect('viewer:fit_main', project_id=project.id, parameter_id=new_parameter.id)

    # Do the fit
    if "fit" in request.POST:
        # Do fit
        fit_result = symmetrical_fit(parameter, project_lipids, datas)
        fit_parameters = fit_result.params

        # Copy current instance
        new_parameter = deepcopy(parameter)

        # Set title
        new_parameter.description = now.strftime("Fit %m/%d/%H:%M")

        # Set params
        new_parameter.terminal_methyl_volume = round(fit_parameters['terminal_methyl_volume'].value, 3)
        new_parameter.lipid_area = round(fit_parameters['area_per_lipid'].value, 3)
        new_parameter.headgroup_thickness = round(fit_parameters['headgroup_thickness'].value, 3)

        # Set report
        fit_report = lsq.fit_report(fit_result)
        new_parameter.fit_report = fit_report.split('\n')

        new_parameter.id = None
        new_parameter.save()

        return redirect('viewer:fit_main', project_id=project.id, parameter_id=new_parameter.id)

    # Show stats / graph
    if "statistics" in request.POST:
        show_statistics = True

    if "graphs" in request.POST:
        show_statistics = False

    # X-Ray graphs
    xray_fig, ax = plt.subplots(xray_datas.count(), squeeze=False)
    plt.subplots_adjust(hspace=0.5)

    i = 0
    for xray_data in xray_datas:
        # Data scatter plot
        ax[i, 0].errorbar(
            xray_data.q_value, 
            xray_data.intensity_value, 
            yerr=xray_data.error_value, 
            fmt='.k',
            color='c',
            ecolor='gray', 
            elinewidth=1, 
            capsize=2,
            zorder=0
        )
        ax[i, 0].set_xscale('log')
        ax[i, 0].set_yscale('log')
        ax[i, 0].set_title(xray_data.data_set_title)
        
        # Fit line
        if parameter.fit_report:
            ax[i, 0].plot(
                xray_data.q_value,
                symmetrical_graph(parameter, project_lipids, xray_data),
                color='r',
                label='Best Fit',
                zorder=1
            )
            ax[i, 0].legend()

        # Move to next dataset
        i = i+1

    xray_graph = mpld3.fig_to_html(xray_fig)

    # Neutron graphs
    neutron_fig, ax = plt.subplots(neutron_datas.count(), squeeze=False)
    plt.subplots_adjust(hspace=0.5)

    j = 0
    for neutron_data in neutron_datas:
        # Data scatter plot
        ax[j, 0].errorbar(
            neutron_data.q_value,
            neutron_data.intensity_value,
            yerr=neutron_data.error_value, 
            fmt='.k',
            color='c',
            ecolor='gray', 
            elinewidth=1, 
            capsize=2,
            zorder=0
        )
        ax[j, 0].set_xscale('log')
        ax[j, 0].set_yscale('log')
        ax[j, 0].set_title(neutron_data.data_set_title)

        # Fit line
        if parameter.fit_report:
            ax[j, 0].plot(
                neutron_data.q_value,
                symmetrical_graph(parameter, project_lipids, neutron_data),
                color='r',
                label='Best Fit',
                zorder=1
            )
            ax[j, 0].legend()

        # Move to next dataset
        j = j+1

    neutron_graph = mpld3.fig_to_html(neutron_fig)

    return render(request, 'viewer/fit_main.html', {
        'x_user':x_user,
        'project':project, 
        'parameter':parameter, 
        'xray_datas':xray_datas, 
        'neutron_datas':neutron_datas, 
        'xray_graph': [xray_graph],
        'neutron_graph': [neutron_graph],
        'parameter_update_form':parameter_update_form,
        'fit_result':fit_result,
        'show_stats':show_statistics
    })
