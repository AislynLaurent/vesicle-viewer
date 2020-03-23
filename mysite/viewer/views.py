## IMPORTS
# From django
from django.shortcuts import get_object_or_404, render
from django.shortcuts import redirect
from django.http import HttpResponse

# Other libraries
import csv
from matplotlib import pyplot as plt
import mpld3
from datetime import datetime
import numpy as np
import lmfit as lsq
from copy import deepcopy
import re

# Models
from .models import *

# Forms
from .forms import *

# Other imports
from .symfit import symmetrical_fit
from .symfit import symmetrical_graph
from .asymfit import asymmetrical_fit
from .asymfit import asymmetrical_graph

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
# Lipid
def lipid_detail(request, lipid_name):
    lipid = get_object_or_404(Lipid, lipid_name=lipid_name)
    return render(request, 'viewer/lipid_detail.html', {'lipid':lipid})

def molecule_detail(request, molecule_name):
    molecule = get_object_or_404(Molecule, compound_name=molecule_name)
    return render(request, 'viewer/lipid_detail.html', {'molecule':molecule})

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
    asymmetrical_projects = Project.objects.filter(owner=request.user, model_type='AS').order_by('project_title')
    return render(request, 'viewer/project_list.html', {'symmetrical_projects':symmetrical_projects, 'asymmetrical_projects':asymmetrical_projects})

# Specific project details
def project_detail(request, project_id):
    project = get_object_or_404(Project, id=project_id)
    samples = Sample.objects.filter(project_title_id=project_id)
    
    if project.model_type == "SM":
        project_lipids = Project_Lipid.objects.filter(project_title_id=project_id)
        in_project_lipids = 0
        out_project_lipids = 0
    elif project.model_type == "AS":
        project_lipids = 0
        in_project_lipids = Project_Lipid.objects.filter(project_title_id=project_id, lipid_location='IN')
        out_project_lipids = Project_Lipid.objects.filter(project_title_id=project_id, lipid_location='OUT')

    # mol_fraction = 0
    # for lipid in project_lipids:
    #     mol_fraction = mol_fraction + lipid.lipid_mol_fraction

    return render(
        request,
        'viewer/project_detail.html', {
            'project':project,
            'project_lipids':project_lipids,
            'in_project_lipids':in_project_lipids,
            'out_project_lipids':out_project_lipids,
            # 'total_lipid_mol_fraction':mol_fraction,
            'samples':samples,
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

# Delete project
def project_delete_warning(request, project_id):
    project = get_object_or_404(Project, id=project_id)
    
    if request.method == 'POST':
        project.delete()
        return redirect('viewer:project_list')

    return render(request, 'viewer/project_delete_warning.html', {'project':project})

## Project Lipids
# Add lipids to a project
def project_lipid_new(request, project_id):
    project = get_object_or_404(Project, id=project_id)

    if project.model_type == "SM":
        if request.method == 'POST':
            form = Sym_Project_Lipid_Form(request.POST)
            if form.is_valid():
                lipids = form.save(commit=False)
                lipids.project_title = project
                lipids.lipid_location = "BOTH"
                lipids.save()
                return redirect('viewer:project_detail', project_id=project.id)
        else:
            form = Sym_Project_Lipid_Form()

    elif project.model_type == "AS":
        if request.method == 'POST':
            form = Asym_Project_Lipid_Form(request.POST)
            if form.is_valid():
                lipids = form.save(commit=False)
                
                in_lipid_location = form.cleaned_data['lipid_location']

                lipids.lipid_location = in_lipid_location
                lipids.project_title = project
                lipids.save()
                return redirect('viewer:project_detail', project_id=project.id)
        else:
            form = Asym_Project_Lipid_Form()

    return render(request, 'viewer/form.html', {'form': form})

# Edit lipid
def project_lipid_edit(request, project_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    lipid = get_object_or_404(Project_Lipid, id=lipid_id)

    if project.model_type == "SM":
        if request.method == 'POST':
            form = Sym_Project_Lipid_Form(request.POST, instance=lipid)
            if form.is_valid():
                lipids = form.save(commit=False)
                lipids.save()
                return redirect('viewer:project_detail', project_id=project.id)
        else:
            form = Sym_Project_Lipid_Form(instance=lipid)

    elif project.model_type == "AS":
        if request.method == 'POST':
            form = Asym_Project_Lipid_Form(request.POST, instance=lipid)
            if form.is_valid():
                lipids = form.save(commit=False)

                in_lipid_location = form.cleaned_data['lipid_location']

                lipids.lipid_location = in_lipid_location
                lipids.save()
                return redirect('viewer:project_detail', project_id=project.id)
        else:
            form = Asym_Project_Lipid_Form(instance=lipid)

    return render(request, 'viewer/form.html', {'form': form})

# Delete project lipid
def project_lipid_delete_warning(request, project_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    lipid = get_object_or_404(Project_Lipid, id=lipid_id)
    
    if request.method == 'POST':
        lipid.delete()
        return redirect('viewer:project_detail', project_id=project.id)

    return render(request, 'viewer/delete_warning.html', {'project':project})

## Samples
# Add a new sample
def sample_new(request, project_id):
    project = get_object_or_404(Project, id=project_id)

    if request.method == 'POST':
        form = Sample_Form(request.POST)
        if form.is_valid():
            sample = form.save(commit=False)
            sample.project_title = project
            sample.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Sample_Form()
    return render(request, 'viewer/form.html', {'form': form})

# Specific sample details
def sample_detail(request, project_id, sample_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    datas = Data_Set.objects.filter(sample_title_id=sample_id)
    data_lipids = Data_Lipid.objects.filter(data_set_title__sample_title__id=sample_id).select_related('data_lipid_name')

    if project.model_type == "SM":
        parameters = Symmetrical_Parameters.objects.filter(sample_title_id=sample_id)
    elif project.model_type == "AS":
        parameters = Asymmetrical_Parameters.objects.filter(sample_title_id=sample_id)

    return render(
        request,
        'viewer/sample_detail.html', {
            'project':project,
            'sample':sample,
            'datas':datas,
            'data_lipids':data_lipids,
            'parameters':parameters,
        }
    )

# Edit an existing sample
def sample_edit(request, project_id, sample_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)

    if request.method == 'POST':
        form = Sample_Form(request.POST, instance=sample)
        if form.is_valid():
            post = form.save(commit=False)
            post.owner = request.user
            post.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Sample_Form(instance=sample)

    return render(request, 'viewer/form.html', {'form': form})

# Delete sample
def sample_delete_warning(request, project_id, sample_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    
    if request.method == 'POST':
        sample.delete()
        return redirect('viewer:project_detail', project_id=project.id)

    return render(request, 'viewer/delete_warning.html', {'project':project})

## Sym Parameters
# Add parameters to a project
def symmetrical_parameters_new(request, project_id, sample_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    project_lipids = Project_Lipid.objects.filter(project_title_id=project_id).select_related('project_lipid_name')

    combined_head_volume = 0
    combined_tail_volume = 0

    x = project.system_tempurature

    # Calculate combined volume
    for lipid in project_lipids:
        combined_head_volume = combined_head_volume + (
                lipid.project_lipid_name.hg_volume * lipid.lipid_mol_fraction
            )
        combined_tail_volume = combined_tail_volume + (
                (
                    eval(
                        lipid.project_lipid_name.total_volume_equation
                    ) 
                    - lipid.project_lipid_name.hg_volume
                )
                * lipid.lipid_mol_fraction
            )

    if request.method == 'POST':
        form = Symmetrical_Parameter_Form(request.POST)
        if form.is_valid():
            parameters = form.save(commit=False)

            # overall
            parameters.sample_title = sample

            # Math'ed values
            parameters.headgroup_volume = combined_head_volume
            parameters.chain_volume = combined_tail_volume

            parameters.save()
            return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)
    else:
        form = Symmetrical_Parameter_Form()

    return render(request, 'viewer/form.html', {'form': form})

# Edit parameters
def symmetrical_parameters_edit(request, project_id, sample_id, parameter_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    parameters = get_object_or_404(Symmetrical_Parameters, id=parameter_id)

    if request.method == 'POST':
        form = Symmetrical_Parameter_Form(request.POST, instance=parameters)
        if form.is_valid():
            parameters = form.save(commit=False)
            parameters.save()
            return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)
    else:
        form = Symmetrical_Parameter_Form(instance=parameters)

    return render(request, 'viewer/form.html', {'form': form})

# Delete parameters
def symmetrical_parameter_delete_warning(request, project_id, sample_id, parameter_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    parameter = get_object_or_404(Symmetrical_Parameters, id=parameter_id)
    
    if request.method == 'POST':
        parameter.delete()
        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)

    return render(request, 'viewer/delete_warning.html', {'project':project, 'sample':sample})

## Asym Parameters
# Add parameters to a project
def asymmetrical_parameters_new(request, project_id, sample_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)

    in_project_lipids = Project_Lipid.objects.filter(project_title_id=project_id, lipid_location='IN').select_related('project_lipid_name')
    out_project_lipids = Project_Lipid.objects.filter(project_title_id=project_id, lipid_location='OUT').select_related('project_lipid_name')

    in_combined_head_volume = 0
    in_combined_tail_volume = 0

    out_combined_head_volume = 0
    out_combined_tail_volume = 0

    x = project.system_tempurature

    # Calculate combined volume
    for lipid in in_project_lipids:
        in_combined_head_volume = in_combined_head_volume + (
                lipid.project_lipid_name.hg_volume * lipid.lipid_mol_fraction
            )
        in_combined_tail_volume = in_combined_tail_volume + (
                (
                    eval(lipid.project_lipid_name.total_volume_equation) 
                    - lipid.project_lipid_name.hg_volume
                )
                * lipid.lipid_mol_fraction
            )

    for lipid in out_project_lipids:
        out_combined_head_volume = out_combined_head_volume + (
                lipid.project_lipid_name.hg_volume * lipid.lipid_mol_fraction
            )
        out_combined_tail_volume = out_combined_tail_volume + (
                (
                    eval(lipid.project_lipid_name.total_volume_equation) 
                    - lipid.project_lipid_name.hg_volume
                )
                * lipid.lipid_mol_fraction
            )

    if request.method == 'POST':
        form = Asymmetrical_Parameter_Form(request.POST)
        if form.is_valid():
            parameters = form.save(commit=False)

            # overall
            parameters.sample_title = sample

            # Math'ed values
            parameters.in_headgroup_volume = in_combined_head_volume
            parameters.in_chain_volume = in_combined_tail_volume

            parameters.out_headgroup_volume = out_combined_head_volume
            parameters.out_chain_volume = out_combined_tail_volume

            parameters.save()
            return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)
    else:
        form = Asymmetrical_Parameter_Form()

    return render(request, 'viewer/form.html', {'form': form})

# Edit parameters
def asymmetrical_parameters_edit(request, project_id, sample_id, parameter_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    parameters = get_object_or_404(Asymmetrical_Parameters, id=parameter_id)

    if request.method == 'POST':
        form = Asymmetrical_Parameter_Form(request.POST, instance=parameters)
        if form.is_valid():
            parameters = form.save(commit=False)
            parameters.save()
            return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)
    else:
        form = Asymmetrical_Parameter_Form(instance=parameters)

    return render(request, 'viewer/form.html', {'form': form})

# Delete parameters
def asymmetrical_parameter_delete_warning(request, project_id, sample_id, parameter_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    parameter = get_object_or_404(Asymmetrical_Parameters, id=parameter_id)
    
    if request.method == 'POST':
        parameter.delete()
        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)

    return render(request, 'viewer/delete_warning.html', {'project':project, 'sample':sample})

## Data
# Input data
def data_upload(request, project_id, sample_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)

    q = []
    i = []
    e = []

    # Upload file
    if "data_upload" in request.POST:
        data_upload_form = Data_Upload_Form(request.POST, request.FILES)
        if data_upload_form.is_valid():
            data_info = data_upload_form.save(commit=False)
            data_info.sample_title = sample

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
                
                # check for letters / words / headers / footers
                if not line[:1].isdigit() or re.search('[a-df-zA-DF-Z]', line):
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

            data_info.min_index = 0
            data_info.max_index = len(q)-1

            data_upload_form.save()

        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)
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
def data_edit(request, project_id, sample_id, data_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    data = get_object_or_404(Data_Set, id=data_id)

    if request.method == 'POST':
        form = Data_Form(request.POST, instance=data)
        if form.is_valid():
            data = form.save(commit=False)
            data.save()
            return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)
    else:
        form = Data_Form(instance=data)

    return render(request, 'viewer/form.html', {'form': form})

# Delete data
def data_delete_warning(request, project_id, sample_id, data_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    data = get_object_or_404(Data_Set, id=data_id)
    
    if request.method == 'POST':
        data.delete()
        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)

    return render(request, 'viewer/delete_warning.html', {'project':project, 'sample':sample})

## Data Lipids
# Add lipids adjustment
def data_lipid_new(request, project_id, sample_id, data_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    data = get_object_or_404(Data_Set, id=data_id)

    number_project_lipids = Project_Lipid.objects.filter(project_title_id=project_id).count()
    number_data_lipids = Data_Lipid.objects.filter(data_lipid_name__project_title__id=project_id).count()

    enough_adjustments = False

    if number_data_lipids > number_project_lipids:
        enough_adjustments = True

    if "lipid_info" in request.POST:
        lipid_info_form = Data_Lipid_Form(project_id, request.POST)

        if lipid_info_form.is_valid():
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
                return redirect('viewer:data_lipid_edit', project_id=project.id, sample_id=sample.id, data_id=data.id, lipid_id=existing_lipid.id)
            elif created_lipid:
                return redirect('viewer:data_lipid_edit', project_id=project.id, sample_id=sample.id, data_id=data.id, lipid_id=created_lipid.id)
    else:
        lipid_info_form = Data_Lipid_Form(project_id)

    return render(
        request,
        'viewer/data_lipid_adjustment.html', {
            'project':project,
            'sample':sample,
            'data':data,
            'lipid_info_form':lipid_info_form,
            'enough_adjustments':enough_adjustments
        }
    )

# Edit lipid adjustments
def data_lipid_edit(request, project_id, sample_id, data_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    data = get_object_or_404(Data_Set, id=data_id)
    data_lipid = get_object_or_404(Data_Lipid, id=lipid_id)
    data_lipid_atoms = Data_Lipid_Atom.objects.filter(data_lipid_name=data_lipid)

    if "lipid_info" in request.POST:
        lipid_info_form = Data_Lipid_Form(project_id, request.POST, instance=data_lipid)
    else:
        lipid_info_form = Data_Lipid_Form(project_id, instance=data_lipid)

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
        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)

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
def data_lipid_delete_warning(request, project_id, sample_id, data_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    data_lipid = get_object_or_404(Data_Lipid, id=lipid_id)
    
    if request.method == 'POST':
        data_lipid.delete()
        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)

    return render(request, 'viewer/delete_warning.html', {'project':project, 'sample':sample})

## Fit
# Main fit page
def fit_main(request, project_id, sample_id, parameter_id):
    ## Import
    x_user = get_object_or_404(ExtendedUser, user=request.user)

    # Overall
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)

    if project.model_type == "SM":
        project_lipids = Project_Lipid.objects.filter(project_title_id=project_id).select_related('project_lipid_name')
        parameter = get_object_or_404(Symmetrical_Parameters, id=parameter_id)
    elif project.model_type == "AS":
        in_project_lipids = Project_Lipid.objects.filter(project_title_id=project_id, lipid_location="IN").select_related('project_lipid_name')
        out_project_lipids = Project_Lipid.objects.filter(project_title_id=project_id, lipid_location="OUT").select_related('project_lipid_name')
        parameter = get_object_or_404(Asymmetrical_Parameters, id=parameter_id)
        

    # Data
    datas = Data_Set.objects.filter(sample_title_id=sample_id)
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
    if project.model_type == "SM":
        if "parameter_update" in request.POST:
            parameter_update_form = Symmetrical_Parameter_Fit_Form(request.POST, instance=parameter)
            if parameter_update_form.is_valid():
                parameter = parameter_update_form.save(commit=False)
                parameter.save()
        else:
            parameter_update_form = Symmetrical_Parameter_Fit_Form(instance=parameter)
    elif project.model_type == "AS":
        if "parameter_update" in request.POST:
            parameter_update_form = Asymmetrical_Parameter_Fit_Form(request.POST, instance=parameter)
            if parameter_update_form.is_valid():
                parameter = parameter_update_form.save(commit=False)
                parameter.save()
        else:
            parameter_update_form = Asymmetrical_Parameter_Fit_Form(instance=parameter)

    # Ranges
    xray_ranges = []
    xray_scales = []

    # Update q range for all x-ray datasets
    for xray_data in xray_datas:
        if xray_data.data_set_title in request.POST:
            xray_range_form = Data_Range_Form(request.POST)
            xray_scale_form = Data_Scale_Form(request.POST, instance=xray_data)
            if xray_range_form.is_valid() and xray_scale_form.is_valid():
                # Get range values
                max_value = xray_range_form.cleaned_data['max_value']
                min_value = xray_range_form.cleaned_data['min_value']

                try:
                    # Find the indexes for the closes value in q_value
                    max_index = min(enumerate(xray_data.q_value), key=lambda x: abs(max_value - x[1]))
                    min_index = min(enumerate(xray_data.q_value), key=lambda x: abs(min_value - x[1]))

                    # Set the indexes in the db - +1 to take the closest on the inside for the low end
                    xray_data.max_index = max_index[0]
                    xray_data.min_index = min_index[0] + 1

                    xray_data.save()
                except TypeError:
                    xray_data.save()
        else:
            xray_range_form = Data_Range_Form()
            xray_scale_form = Data_Scale_Form(instance=xray_data)

        xray_ranges.append(xray_range_form)
        xray_scales.append(xray_scale_form)

    # Ranges
    neutron_ranges = []
    neutron_scales = []

    # Update q range for all x-ray datasets
    for neutron_data in neutron_datas:
        if neutron_data.data_set_title in request.POST:
            neutron_range_form = Data_Range_Form(request.POST)
            neutron_scale_form = Data_Scale_Form(request.POST, instance=neutron_data)
            if neutron_range_form.is_valid() and neutron_scale_form.is_valid():
                # Get range values
                max_value = neutron_range_form.cleaned_data['max_value']
                min_value = neutron_range_form.cleaned_data['min_value']

                try:
                    # Find the indexes for the closes value in q_value
                    max_index = min(enumerate(neutron_data.q_value), key=lambda x: abs(max_value - x[1]))
                    min_index = min(enumerate(neutron_data.q_value), key=lambda x: abs(min_value - x[1]))

                    # Set the indexes in the db - +1 to take the closest on the inside for the low end
                    neutron_data.max_index = max_index[0]
                    neutron_data.min_index = min_index[0] + 1

                    neutron_data.save()
                except TypeError:
                    neutron_data.save()
        else:
            neutron_range_form = Data_Range_Form()
            neutron_scale_form = Data_Scale_Form(instance=neutron_data)

        neutron_ranges.append(neutron_range_form)
        neutron_scales.append(neutron_scale_form)

    # Do the fit
    if project.model_type == "SM":
        if "fit" in request.POST:
            # Do fit
            fit_result = symmetrical_fit(parameter, project_lipids, datas, project.system_tempurature)
            fit_parameters = fit_result.params

            # Copy current instance
            new_parameter = deepcopy(parameter)

            # Set title
            new_parameter.name = now.strftime("%m/%d/%H:%M")

            # Set params
            new_parameter.terminal_methyl_volume = round(fit_parameters['terminal_methyl_volume'].value, 6)
            new_parameter.lipid_area = round(fit_parameters['area_per_lipid'].value, 6)
            new_parameter.headgroup_thickness = round(fit_parameters['headgroup_thickness'].value, 6)

            # Set report
            fit_report = lsq.fit_report(fit_result)
            new_parameter.fit_report = fit_report.split('\n')

            new_parameter.id = None
            new_parameter.save()

            for data in datas:
                data.scale = fit_parameters['scale_%i' % data.id].value
                data.background = fit_parameters['background_%i' % data.id].value

                data.save()

            # print(lsq.fit_report(fit_result))

            return redirect('viewer:fit_main', project_id=project.id, sample_id=sample.id, parameter_id=new_parameter.id)

    elif project.model_type == "AS":
        if "fit" in request.POST:
            # Do fit
            fit_result = asymmetrical_fit(parameter, in_project_lipids, out_project_lipids, datas, project.system_tempurature)
            fit_parameters = fit_result.params

            # Copy current instance
            new_parameter = deepcopy(parameter)

            # Set title
            new_parameter.name = now.strftime("%m/%d/%H:%M")

            # Set params
            new_parameter.in_terminal_methyl_volume = round(fit_parameters['in_terminal_methyl_volume'].value, 6)
            new_parameter.in_lipid_area = round(fit_parameters['in_area_per_lipid'].value, 6)
            new_parameter.in_headgroup_thickness = round(fit_parameters['in_headgroup_thickness'].value, 6)

            new_parameter.out_terminal_methyl_volume = round(fit_parameters['out_terminal_methyl_volume'].value, 6)
            new_parameter.out_lipid_area = round(fit_parameters['out_area_per_lipid'].value, 6)
            new_parameter.out_headgroup_thickness = round(fit_parameters['out_headgroup_thickness'].value, 6)

            # Set report
            fit_report = lsq.fit_report(fit_result)
            new_parameter.fit_report = fit_report.split('\n')

            new_parameter.id = None
            new_parameter.save()

            for data in datas:
                data.scale = fit_parameters['scale_%i' % data.id].value
                data.background = fit_parameters['background_%i' % data.id].value

                data.save()

            # print(lsq.fit_report(fit_result))

            return redirect('viewer:fit_main', project_id=project.id, sample_id=sample.id, parameter_id=new_parameter.id)

    # caluclated values
    calculated_i_values = []

    # Download data
    if "download" in request.POST:
        # Create the HttpResponse object with the appropriate CSV header.
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="vv_fitdata_download.csv"'

        writer = csv.writer(response)
        writer.writerow(['VesicleViewer Data output', now])
        writer.writerow(['Project Name', 'Sample Name', 'Parameter Set'])
        writer.writerow([project.project_title, sample.sample_title, parameter.name])
        writer.writerow([])

        writer.writerow({'Fit Statistics'})
        for line in parameter.fit_report:
            writer.writerow([line])

        writer.writerow([])

        writer.writerow(['Q', 'Experimental i', 'Calculated i'])

        for xray_data in xray_datas:
            writer.writerow([xray_data.data_set_title])
            if project.model_type == "SM":
                calculated_i_values = symmetrical_graph(parameter, project_lipids, xray_data, project.system_tempurature)
            elif project.model_type == "AS":
                calculated_i_values = asymmetrical_graph(parameter, in_project_lipids, out_project_lipids, xray_data, project.system_tempurature)

            j = 0
            for i in range(xray_data.min_index, xray_data.max_index):
                writer.writerow([xray_data.q_value[i], xray_data.intensity_value[i], calculated_i_values[j]])
                j = j+1

        for neutron_data in neutron_datas:
            writer.writerow([neutron_data.data_set_title])
            if project.model_type == "SM":
                calculated_i_values = symmetrical_graph(parameter, project_lipids, neutron_data, project.system_tempurature)
            elif project.model_type == "AS":
                calculated_i_values = asymmetrical_graph(parameter, in_project_lipids, out_project_lipids, neutron_data, project.system_tempurature)

            j = 0
            for i in range(neutron_data.min_index, neutron_data.max_index):
                writer.writerow([neutron_data.q_value[i], neutron_data.intensity_value[i], calculated_i_values[j]])
                j = j+1

        return response

    # Show stats / graph
    if "statistics" in request.POST:
        show_statistics = True

    if "graphs" in request.POST:
        show_statistics = False

    # Graphs
    xray_figures = []
    neutron_figures = []

    # X-Ray graphs
    for xray_data in xray_datas:
        xray_fig = plt.figure(figsize=(6,4.5))

        # Data scatter plot
        plt.errorbar(
            xray_data.q_value[xray_data.min_index:xray_data.max_index], 
            xray_data.intensity_value[xray_data.min_index:xray_data.max_index],
            xray_data.error_value[xray_data.min_index:xray_data.max_index],
            fmt='o',
            color='c',
            mfc='w',
            ecolor='gray', 
            elinewidth=0.5, 
            capsize=1,
            zorder=0
        )
        plt.xscale('log')
        plt.xlabel('q(A-1)')

        plt.yscale('log')
        plt.ylabel('Intensity (A.U.)')

        plt.title(xray_data.data_set_title)

        # Fit line
        if parameter.fit_report:
            if project.model_type == "SM":
                plt.plot(
                    xray_data.q_value[xray_data.min_index:xray_data.max_index],
                    symmetrical_graph(parameter, project_lipids, xray_data, project.system_tempurature),
                    color='r',
                    label='Best Fit',
                    zorder=1
                )
            elif project.model_type == "AS":
                plt.plot(
                    xray_data.q_value[xray_data.min_index:xray_data.max_index],
                    asymmetrical_graph(parameter, in_project_lipids, out_project_lipids, xray_data, project.system_tempurature),
                    color='r',
                    label='Best Fit',
                    zorder=1
                )

        xray_figures.append(mpld3.fig_to_html(xray_fig))

    # Neutron graphs
    for neutron_data in neutron_datas:
        neutron_fig = plt.figure(figsize=(6,4.5))
        # Data scatter plot
        plt.errorbar(
            neutron_data.q_value[neutron_data.min_index:neutron_data.max_index],
            neutron_data.intensity_value[neutron_data.min_index:neutron_data.max_index],
            neutron_data.error_value[neutron_data.min_index:neutron_data.max_index],
            fmt='o',
            color='c',
            mfc='w',
            ecolor='gray', 
            elinewidth=1, 
            capsize=2,
            zorder=0
        )
        plt.xscale('log')
        plt.xlabel('q(A-1)')

        plt.yscale('log')
        plt.ylabel('Intensity (A.U.)')

        plt.title(neutron_data.data_set_title)

        # Fit line
        if parameter.fit_report:
            if project.model_type == "SM":
                plt.plot(
                    neutron_data.q_value[neutron_data.min_index:neutron_data.max_index],
                    symmetrical_graph(parameter, project_lipids, neutron_data, project.system_tempurature),
                    color='r',
                    label='Best Fit',
                    zorder=1
                )
            elif project.model_type == "AS":
                plt.plot(
                    neutron_data.q_value[neutron_data.min_index:neutron_data.max_index],
                    asymmetrical_graph(parameter, in_project_lipids, out_project_lipids, neutron_data, project.system_tempurature),
                    color='r',
                    label='Best Fit',
                    zorder=1
                )

        neutron_figures.append(mpld3.fig_to_html(neutron_fig))

    xray_graphs_and_forms = zip(xray_figures, xray_ranges, xray_scales, xray_datas)
    neutron_graphs_and_forms = zip(neutron_figures, neutron_ranges, neutron_scales, neutron_datas)

    return render(request, 'viewer/fit_main.html', {
        'x_user':x_user,
        'project':project,
        'sample':sample,
        'parameter':parameter,
        'parameter_update_form':parameter_update_form,
        'fit_result':fit_result,
        'show_stats':show_statistics,
        'xray_graphs_and_forms':xray_graphs_and_forms,
        'neutron_graphs_and_forms':neutron_graphs_and_forms
    })
