## IMPORTS
# From django
from django.shortcuts import get_object_or_404, render
from django.shortcuts import redirect
from django.http import HttpResponse
from django.utils import timezone
from django.utils import text
from django.contrib import messages 

# Other libraries
import csv

import matplotlib
matplotlib.use('Agg')
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
from .symfit import *
from .asymfit import *
from .probabilities import *

from viewer.fit import generate_fit_main

## STATIC PAGES
# Home
def index(request):
    ## Import
    if request.user.is_anonymous:
        xuser_tutorial = False
    else:
        xuser = ExtendedUser.objects.get(user=request.user)
        xuser_tutorial = xuser.display_tutorial
    
    tutorial = True

    # Forms
    # Dismiss the tutorial
    if "dismiss_this" in request.POST:
        tutorial = False

    # Dismiss all tutorials
    if "dismiss_all" in request.POST:
        xuser.display_tutorial = False
        xuser.save()

    return render(
        request,
        'viewer/index.html', {
            'tutorial':tutorial,
            'xuser_tutorial':xuser_tutorial,
        }
    )

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

# Tutorial Form
def enable_tutorials(request):
    ## Import
    x_user = get_object_or_404(ExtendedUser, user=request.user)

    if request.method == 'POST':
        form = Tutorial_Form(request.POST, instance=x_user)
        if form.is_valid():
            post = form.save(commit=False)
            post.owner = request.user
            post.save()
            return redirect('viewer:index')
    else:
        form = Tutorial_Form(instance=x_user)

    return render(request, 'viewer/form.html', {'form': form})

## MODEL PAGES
# Lipid
def lipid_detail(request, lipid_name):
    lipid = get_object_or_404(Lipid, slug=lipid_name)
    return render(request, 'viewer/lipid_detail.html', {'lipid':lipid})

# Custom lipid
def user_lipid_detail(request, owner, lipid_name):
    lipid = get_object_or_404(User_Lipid, slug=lipid_name, owner=owner)
    return render(request, 'viewer/lipid_detail.html', {'lipid':lipid})

# New custom lipid
def user_lipid_new(request, owner):
    if request.method == 'POST':
        form = User_Lipid_Form(request.POST)
        if form.is_valid():
            lipids = form.save(commit=False)
            lipids.owner = request.user
            lipids.save()
            return redirect('viewer:project_list')
    else:
        form = User_Lipid_Form()

    return render(request, 'viewer/form.html', {'form': form})

# Edit user lipid
def user_lipid_edit(request, owner, lipid_name):
    lipid = get_object_or_404(User_Lipid, slug=lipid_name, owner=owner)

    if request.method == 'POST':
        form = User_Lipid_Form(request.POST, instance=lipid)
        if form.is_valid():
            lipids = form.save(commit=False)
            lipids.owner = request.user
            lipids.save()
            return redirect('viewer:project_list')
    else:
        form = User_Lipid_Form(instance=lipid)

    return render(request, 'viewer/form.html', {'form': form})

# Molecule
def molecule_detail(request, molecule_name):
    molecule = get_object_or_404(Molecule, slug=molecule_name)
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
    ## Tutorials
    if request.user.is_anonymous:
        xuser_tutorial = False
    else:
        xuser = ExtendedUser.objects.get(user=request.user)
        xuser_tutorial = xuser.display_tutorial
    tutorial = True
    # Dismiss the tutorial
    if "dismiss_this" in request.POST:
        tutorial = False

    # Dismiss all tutorials
    if "dismiss_all" in request.POST:
        xuser.display_tutorial = False
        xuser.save()
    
    ## List
    symmetrical_projects = Project.objects.filter(owner=request.user, model_type='SM').order_by('project_title')
    asymmetrical_projects = Project.objects.filter(owner=request.user, model_type='AS').order_by('project_title')
    user_lipids = User_Lipid.objects.filter(owner=request.user).order_by('user_lipid_name')

    return render(
        request,
        'viewer/project_list.html', {
            'tutorial':tutorial,
            'xuser_tutorial':xuser_tutorial,
            'symmetrical_projects':symmetrical_projects,
            'asymmetrical_projects':asymmetrical_projects,
            'user_lipids':user_lipids,
        }
    )

# Specific project details
def project_detail(request, project_id):
    ## Tutorials
    if request.user.is_anonymous:
        xuser_tutorial = False
    else:
        xuser = ExtendedUser.objects.get(user=request.user)
        xuser_tutorial = xuser.display_tutorial

    tutorial = True
    # Dismiss the tutorial
    if "dismiss_this" in request.POST:
        tutorial = False

    # Dismiss all tutorials
    if "dismiss_all" in request.POST:
        xuser.display_tutorial = False
        xuser.save()

    project = get_object_or_404(Project, id=project_id)
    samples = Sample.objects.filter(project_title_id=project_id)
    project_lipids = Project_Lipid.objects.filter(project_title_id=project_id)
    project_user_lipid_volumes = Project_User_Lipid_Volume.objects.filter(project_title_id=project_id)

    x = project.system_tempurature

    # Calculate or retrieve the volume and pair it with the lipid for display
        # Workaround for PostgreSQL type issue
    lipids_and_volumes = []
    for lipid in project_lipids:
        if lipid.project_lipid_name:
            calc_volume = round(eval(lipid.project_lipid_name.total_volume_equation) - lipid.project_lipid_name.hg_volume, 2)
            lipids_and_volumes.append([lipid, calc_volume])
        elif project_user_lipid_volumes:
            check = False
            for volume in project_user_lipid_volumes:
                if volume.project_user_lipid_name == lipid.project_user_lipid_name:
                    lipids_and_volumes.append([lipid, volume])
                    check = True
                    break
                
            if not check:
                lipids_and_volumes.append([lipid, 0])
        else:
            lipids_and_volumes.append([lipid, 0])

    return render(
        request,
        'viewer/project_detail.html', {
            'tutorial':tutorial,
            'xuser_tutorial':xuser_tutorial,
            'project':project,
            'lipids_and_volumes':lipids_and_volumes,
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

    return render(request, 'viewer/form.html', {'project':project, 'form': form})

# Allow users to use advanced options (vary volume parameters mostly)
def project_advanced_options(request, project_id):
    project = get_object_or_404(Project, id=project_id)

    if request.method == 'POST':
        form = Advanced_Options(request.POST, instance=project)
        if form.is_valid():
            post = form.save(commit=False)
            post.owner = request.user
            post.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Advanced_Options(instance=project)

    return render(request, 'viewer/form.html', {'project':project, 'form': form})

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
    owner = request.user

    if request.method == 'POST':
        form = Project_Lipid_Form(owner, request.POST)
        if form.is_valid():
            lipids = form.save(commit=False)
            lipids.project_title = project
            lipids.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Project_Lipid_Form(owner)

    return render(request, 'viewer/form.html', {'project':project, 'form': form})

# Delete project lipid
def project_lipid_delete_warning(request, project_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    lipid = get_object_or_404(Project_Lipid, id=lipid_id)
    
    if request.method == 'POST':
        lipid.delete()
        return redirect('viewer:project_detail', project_id=project.id)

    return render(request, 'viewer/delete_warning.html', {'project':project, 'project_lipid':lipid})

# Project custom lipid volume
def project_user_lipid_volume_new(request, project_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    lipid = get_object_or_404(User_Lipid, id=lipid_id)

    if request.method == 'POST':
        form = Project_User_Lipid_Volume_Form(request.POST)
        if form.is_valid():
            user_lipid = form.save(commit=False)
            user_lipid.project_title = project
            user_lipid.project_user_lipid_name = lipid
            user_lipid.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Project_User_Lipid_Volume_Form()
    return render(request, 'viewer/form.html', {'form': form})

# Edit project custom lipid volume
def project_user_lipid_volume_edit(request, project_id, volume_id):
    project = get_object_or_404(Project, id=project_id)
    lipid_volume = get_object_or_404(Project_User_Lipid_Volume, id=volume_id)

    if request.method == 'POST':
        form = Project_User_Lipid_Volume_Form(request.POST, instance=lipid_volume)
        if form.is_valid():
            post = form.save(commit=False)
            post.owner = request.user
            post.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Project_User_Lipid_Volume_Form(instance=lipid_volume)
    return render(request, 'viewer/form.html', {'form': form})

## Samples
# Add a new sample
def sample_new(request, project_id):
    project = get_object_or_404(Project, id=project_id)

    if request.method == 'POST':
        form = Sample_Form(project.id, request.POST)
        if form.is_valid():
            sample = form.save(commit=False)
            sample.project_title = project
            sample.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Sample_Form(project.id)
    return render(request, 'viewer/form.html', {'project':project, 'form': form})

# Specific sample details
def sample_detail(request, project_id, sample_id):
    ## Tutorials
    if request.user.is_anonymous:
        xuser_tutorial = False
    else:
        xuser = ExtendedUser.objects.get(user=request.user)
        xuser_tutorial = xuser.display_tutorial

    tutorial = True
    # Dismiss the tutorial
    if "dismiss_this" in request.POST:
        tutorial = False

    # Dismiss all tutorials
    if "dismiss_all" in request.POST:
        xuser.display_tutorial = False
        xuser.save()

    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)

    data_xr = Data_Set.objects.filter(sample_title_id=sample_id, data_type='XR')
    data_nu = Data_Set.objects.filter(sample_title_id=sample_id, data_type='NU')
    
    lipids_augments_both = []
    lipids_augments_in = []
    lipids_augments_out = []

    sample_lipids_both = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='BOTH')
    sample_lipids_in = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='IN')
    sample_lipids_out = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='OUT')

    for lipid in sample_lipids_both:
        augments = Data_Sample_Lipid_Augment.objects.filter(sample_lipid_name=lipid)
        combo = [lipid, augments]
        lipids_augments_both.append(combo)
    for lipid in sample_lipids_in:
        augments = Data_Sample_Lipid_Augment.objects.filter(sample_lipid_name=lipid)
        combo = [lipid, augments]
        lipids_augments_in.append(combo)
    for lipid in sample_lipids_out:
        augments = Data_Sample_Lipid_Augment.objects.filter(sample_lipid_name=lipid)
        combo = [lipid, augments]
        lipids_augments_out.append(combo)

    total_mols = 0
    total_mols_in = 0
    total_mols_out = 0

    if project.model_type == "SM":
        parameters = Symmetrical_Parameters.objects.filter(sample_title_id=sample_id)

        total_mols_in = 1
        total_mols_out = 1

        for lipid in sample_lipids_both:
            total_mols = total_mols + lipid.lipid_mol_fraction

    elif project.model_type == "AS":
        parameters = Asymmetrical_Parameters.objects.filter(sample_title_id=sample_id)

        total_mols = 1

        for lipid in sample_lipids_in:
            total_mols_in = total_mols_in + lipid.lipid_mol_fraction

        for lipid in sample_lipids_out:
            total_mols_out = total_mols_out + lipid.lipid_mol_fraction

    return render(
        request,
        'viewer/sample_detail.html', {
            'tutorial':tutorial,
            'xuser_tutorial':xuser_tutorial,
            'project':project,
            'sample':sample,
            'sample_lipids_both':sample_lipids_both,
            'sample_lipids_in':sample_lipids_in,
            'sample_lipids_out':sample_lipids_out,
            'lipids_augments_both': lipids_augments_both,
            'lipids_augments_in': lipids_augments_in,
            'lipids_augments_out': lipids_augments_out,
            'data_xr':data_xr,
            'data_nu':data_nu,
            'parameters':parameters,
            'total_mols':total_mols,
            'total_mols_in':total_mols_in,
            'total_mols_out':total_mols_out,
        }
    )

# Edit an existing sample
def sample_edit(request, project_id, sample_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)

    if request.method == 'POST':
        form = Sample_Form(project.id, request.POST, instance=sample)
        if form.is_valid():
            post = form.save(commit=False)
            post.owner = request.user
            post.save()
            return redirect('viewer:project_detail', project_id=project.id)
    else:
        form = Sample_Form(project.id, instance=sample)

    return render(request, 'viewer/form.html', {'project':project, 'sample':sample, 'form': form})

# Delete sample
def sample_delete_warning(request, project_id, sample_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    
    if request.method == 'POST':
        sample.delete()
        return redirect('viewer:project_detail', project_id=project.id)

    return render(request, 'viewer/delete_warning.html', {'project':project, 'sample':sample})

# Add lipids to a sample
def sample_lipid_new(request, project_id, sample_id):
    ## Tutorials
    if request.user.is_anonymous:
        xuser_tutorial = False
    else:
        xuser = ExtendedUser.objects.get(user=request.user)
        xuser_tutorial = xuser.display_tutorial

    tutorial = True
    # Dismiss the tutorial
    if "dismiss_this" in request.POST:
        tutorial = False

    # Dismiss all tutorials
    if "dismiss_all" in request.POST:
        xuser.display_tutorial = False
        xuser.save()

    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)

    if project.model_type == "SM":
        if "lipid_info" in request.POST:
            lipid_form = Sym_Sample_Lipid_Form(project.id, request.POST)
            if lipid_form.is_valid():
                in_sample_lipid_name = lipid_form.cleaned_data['sample_lipid_name']
                in_lipid_mol_fraction = lipid_form.cleaned_data['lipid_mol_fraction']
                in_lipid_location = "BOTH"

                existing_lipid, created_lipid = Sample_Lipid.objects.update_or_create(
                    sample_lipid_name=in_sample_lipid_name,
                    defaults={
                        'sample_title':sample,
                        'lipid_mol_fraction':in_lipid_mol_fraction,
                        'lipid_location':in_lipid_location,
                })

                if existing_lipid:
                    return redirect('viewer:sample_lipid_edit', project_id=project.id, sample_id=sample.id, lipid_id=existing_lipid.id)
                elif created_lipid:
                    return redirect('viewer:sample_lipid_edit', project_id=project.id, sample_id=sample.id, lipid_id=created_lipid.id)
        else:
            lipid_form = Sym_Sample_Lipid_Form(project.id)

    elif project.model_type == "AS":
        if "lipid_info" in request.POST:
            lipid_form = Asym_Sample_Lipid_Form(project.id, request.POST)
            if lipid_form.is_valid():
                in_sample_lipid_name = lipid_form.cleaned_data['sample_lipid_name']
                in_lipid_mol_fraction = lipid_form.cleaned_data['lipid_mol_fraction']
                in_lipid_location = lipid_form.cleaned_data['location']

                existing_lipid, created_lipid = Sample_Lipid.objects.update_or_create(
                    sample_lipid_name=in_sample_lipid_name,
                    lipid_location = in_lipid_location,
                    defaults={
                        'sample_title':sample,
                        'lipid_mol_fraction':in_lipid_mol_fraction,
                })

                if existing_lipid:
                    return redirect('viewer:sample_lipid_edit', project_id=project.id, sample_id=sample.id, lipid_id=existing_lipid.id)
                elif created_lipid:
                    return redirect('viewer:sample_lipid_edit', project_id=project.id, sample_id=sample.id, lipid_id=created_lipid.id)
        else:
            lipid_form = Asym_Sample_Lipid_Form(project.id)

    return render(
        request,
        'viewer/sample_lipid_form.html', {
            'tutorial':tutorial,
            'xuser_tutorial':xuser_tutorial,
            'project':project,
            'sample':sample,
            'lipid_form': lipid_form,
        })

# Edit sample lipid
def sample_lipid_edit(request, project_id, sample_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    lipid = get_object_or_404(Sample_Lipid, id=lipid_id)
    custom_augment = Sample_Lipid_Augmentation.objects.filter(sample_lipid_name=lipid)
    data_augments = Data_Sample_Lipid_Augment.objects.filter(sample_lipid_name=lipid)

    data_augment_forms = []
    
    if project.model_type == "SM":
        if "lipid_info" in request.POST:
            lipid_form = Sym_Sample_Lipid_Form(project.id, request.POST, instance=lipid)
            if lipid_form.is_valid():
                in_sample_lipid_name = lipid_form.cleaned_data['sample_lipid_name']
                in_lipid_mol_fraction = lipid_form.cleaned_data['lipid_mol_fraction']
                in_lipid_location = "BOTH"

                existing_lipid, created_lipid = Sample_Lipid.objects.update_or_create(
                    sample_lipid_name=in_sample_lipid_name,
                    defaults={
                        'sample_title':sample,
                        'lipid_mol_fraction':in_lipid_mol_fraction,
                        'lipid_location':in_lipid_location,
                })
        else:
            lipid_form = Sym_Sample_Lipid_Form(project.id, instance=lipid)

    elif project.model_type == "AS":
        if "lipid_info" in request.POST:
            lipid_form = Asym_Sample_Lipid_Form(project.id, request.POST, instance=lipid)
            if lipid_form.is_valid():
                in_sample_lipid_name = lipid_form.cleaned_data['sample_lipid_name']
                in_lipid_mol_fraction = lipid_form.cleaned_data['lipid_mol_fraction']
                in_lipid_location = lipid_form.cleaned_data['lipid_location']

                existing_lipid, created_lipid = Sample_Lipid.objects.update_or_create(
                    sample_lipid_name=in_sample_lipid_name,
                    lipid_location = in_lipid_location,
                    defaults={
                        'sample_title':sample,
                        'lipid_mol_fraction':in_lipid_mol_fraction,
                })
        else:
            lipid_form = Asym_Sample_Lipid_Form(project.id, instance=lipid)

    if "augment" in request.POST:
        augment_form = Lipid_Augmentation_Form(lipid.sample_lipid_name, sample.id, request.POST)
        if augment_form.is_valid():
            sample_lipid_augment_in = augment_form.cleaned_data['sample_lipid_augment']
            sample_lipid_custom_augment_in = augment_form.cleaned_data['sample_lipid_custom_augment']
            data_set_title_in = augment_form.cleaned_data['data_set_title']

            existing_augment, created_augment = Data_Sample_Lipid_Augment.objects.update_or_create(
                sample_lipid_name=lipid,
                data_set_title=data_set_title_in,
                defaults={
                    'sample_lipid_augment':sample_lipid_augment_in,
                    'sample_lipid_custom_augment':sample_lipid_custom_augment_in,
            })

            return redirect('viewer:sample_lipid_edit', project_id=project.id, sample_id=sample.id, lipid_id=lipid.id)
    else:
        augment_form = Lipid_Augmentation_Form(lipid.sample_lipid_name, sample.id)
    
    if data_augments:
        for data_augment in data_augments:
            if 'augment'+str(data_augment.id) in request.POST:
                update_augment_form = Lipid_Augmentation_Form(lipid.sample_lipid_name, sample.id, request.POST, instance=data_augment)
                if update_augment_form.is_valid():
                    sample_lipid_augment_in = update_augment_form.cleaned_data['sample_lipid_augment']
                    sample_lipid_custom_augment_in = update_augment_form.cleaned_data['sample_lipid_custom_augment']
                    data_set_title_in = update_augment_form.cleaned_data['data_set_title']

                    existing_augment, created_augment = Data_Sample_Lipid_Augment.objects.update_or_create(
                        sample_lipid_name=lipid,
                        data_set_title=data_set_title_in,
                        defaults={
                            'sample_lipid_augment':sample_lipid_augment_in,
                            'sample_lipid_custom_augment':sample_lipid_custom_augment_in,
                    })

                    return redirect('viewer:sample_lipid_edit', project_id=project.id, sample_id=sample.id, lipid_id=lipid.id)
            else:
                update_augment_form = Lipid_Augmentation_Form(lipid.sample_lipid_name, sample.id, instance=data_augment)
                
            data_augment_forms.append(update_augment_form)

            if 'delete'+str(data_augment.id) in request.POST:
                data_augment.delete()
                return redirect('viewer:sample_lipid_edit', project_id=project.id, sample_id=sample.id, lipid_id=lipid.id)

    augment_form_data_set = zip(data_augment_forms, data_augments)

    if "done" in request.POST:
        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)

    return render(
        request,
        'viewer/sample_lipid_form.html', {
            'project':project,
            'sample':sample,
            'lipid': lipid,
            'custom_augment': custom_augment,
            'lipid_form': lipid_form,
            'augment_form': augment_form,
            'augment_form_data_set': augment_form_data_set,
        })

# Add a custom augment
def sample_custom_lipid_edit(request, project_id, sample_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    lipid = get_object_or_404(Sample_Lipid, id=lipid_id)

    try:
        custom_augment = Sample_Lipid_Augmentation.objects.get(sample_lipid_name=lipid)
    except Sample_Lipid_Augmentation.DoesNotExist:
        custom_augment = False

    if request.method == 'POST':
        if custom_augment:
            form = Custom_Lipid_Augmentation_Form(request.POST, instance=custom_augment)
        else:
            form = Custom_Lipid_Augmentation_Form(request.POST)
        if form.is_valid():

            aug_suffix = form.cleaned_data['augmentation_suffix']
            hg_change = form.cleaned_data['hg_scattering_net_change']
            tg_change = form.cleaned_data['tg_scattering_net_change']
            tm_change = form.cleaned_data['tmg_scattering_net_change']

            existing_augment, created_augment = Sample_Lipid_Augmentation.objects.update_or_create(
                sample_lipid_name=lipid,
                defaults={
                    'augmentation_suffix':aug_suffix,
                    'hg_scattering_net_change':hg_change,
                    'tg_scattering_net_change':tg_change,
                    'tmg_scattering_net_change':tm_change,
            })

            return redirect('viewer:sample_lipid_edit', project_id=project.id, sample_id=sample.id, lipid_id=lipid.id)
    else:
        if custom_augment:
            form = Custom_Lipid_Augmentation_Form(instance=custom_augment)
        else:
            form = Custom_Lipid_Augmentation_Form()

    return render(request, 'viewer/form.html', {'project':project, 'sample':sample, 'form': form})

# Delete sample lipid
def sample_lipid_delete_warning(request, project_id, sample_id, lipid_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    lipid = get_object_or_404(Sample_Lipid, id=lipid_id)
    
    if request.method == 'POST':
        lipid.delete()
        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)

    return render(request, 'viewer/delete_warning.html', {'project':project, 'sample':sample, 'sample_lipid':lipid})

## Sym Parameters
# Add parameters to a project
def symmetrical_parameters_new(request, project_id, sample_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    sample_lipids = Sample_Lipid.objects.filter(sample_title_id=sample_id)
    volumes = Project_User_Lipid_Volume.objects.filter(project_title_id=project_id)

    combined_head_volume = 0
    combined_tail_volume = 0

    x = project.system_tempurature

    # Calculate combined volume
    for lipid in sample_lipids:
        if lipid.sample_lipid_name.project_lipid_name:
            combined_head_volume = combined_head_volume + (
                    lipid.sample_lipid_name.project_lipid_name.hg_volume * lipid.lipid_mol_fraction
                )
            combined_tail_volume = combined_tail_volume + (
                    (
                        eval(lipid.sample_lipid_name.project_lipid_name.total_volume_equation) - lipid.sample_lipid_name.project_lipid_name.hg_volume
                    )* lipid.lipid_mol_fraction
                )
        else:
            for volume in volumes:
                if volume.project_user_lipid_name == lipid.sample_lipid_name.project_user_lipid_name:
                    combined_head_volume = combined_head_volume + (
                            lipid.sample_lipid_name.project_user_lipid_name.hg_volume * lipid.lipid_mol_fraction
                        )
                    combined_tail_volume = combined_tail_volume + (
                            (
                                volume.user_lipid_volume - lipid.sample_lipid_name.project_user_lipid_name.hg_volume
                            ) * lipid.lipid_mol_fraction
                        )
                    break
    
    combined_head_volume = round(combined_head_volume, 2)
    combined_tail_volume = round(combined_tail_volume, 2)

    if request.method == 'POST':
        form = Symmetrical_Parameter_Form(request.POST) 
        if form.is_valid():
            parameters = form.save(commit=False)

            # overall
            parameters.sample_title = sample

            # Math'ed values
            # Head
            parameters.headgroup_volume = combined_head_volume
            if combined_head_volume == 0:
                parameters.headgroup_volume_upperbound = 1
                parameters.headgroup_volume_lowerbound = -1
            else:
                parameters.headgroup_volume_upperbound = (abs(combined_head_volume)*1.20)
                parameters.headgroup_volume_lowerbound = -(abs(combined_head_volume)*1.20)

            # Chain
            parameters.chain_volume = combined_tail_volume
            if combined_tail_volume == 0:
                parameters.chain_volume_upperbound = 1
                parameters.chain_volume_lowerbound = -1
            else:
                parameters.chain_volume_upperbound = (abs(combined_tail_volume)*1.10)
                parameters.chain_volume_lowerbound = -(abs(combined_tail_volume)*1.10)

            parameters.save()
            return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)
    else:
        form = Symmetrical_Parameter_Form()

    return render(request, 'viewer/form.html', {'project':project, 'sample':sample, 'form': form})

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

    return render(request, 'viewer/form.html', {'project':project, 'sample':sample, 'parameters':parameters, 'form': form})

# Delete parameters
def symmetrical_parameter_delete_warning(request, project_id, sample_id, parameter_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    parameter = get_object_or_404(Symmetrical_Parameters, id=parameter_id)
    
    if request.method == 'POST':
        parameter.delete()
        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)

    return render(request, 'viewer/delete_warning.html', {'project':project, 'sample':sample, 'parameters':parameter})

## Asym Parameters
# Add parameters to a project
def asymmetrical_parameters_new(request, project_id, sample_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    volumes = Project_User_Lipid_Volume.objects.filter(project_title_id=project_id)

    sample_lipids_in = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='IN')
    sample_lipids_out = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='OUT')

    in_combined_head_volume = 0
    in_combined_tail_volume = 0

    out_combined_head_volume = 0
    out_combined_tail_volume = 0

    x = project.system_tempurature

    # Calculate combined volume
    for lipid in sample_lipids_in:
        if lipid.sample_lipid_name.project_lipid_name:
            in_combined_head_volume = in_combined_head_volume + (
                    lipid.sample_lipid_name.project_lipid_name.hg_volume * lipid.lipid_mol_fraction
                )
            in_combined_tail_volume = in_combined_tail_volume + (
                    (
                        eval(lipid.sample_lipid_name.project_lipid_name.total_volume_equation) - lipid.sample_lipid_name.project_lipid_name.hg_volume
                    )
                    * lipid.lipid_mol_fraction
                )
        else:
            for volume in volumes:
                if volume.project_user_lipid_name == lipid.sample_lipid_name.project_user_lipid_name:
                    in_combined_head_volume = in_combined_head_volume + (
                            lipid.sample_lipid_name.project_user_lipid_name.hg_volume * lipid.lipid_mol_fraction
                        )
                    in_combined_tail_volume = in_combined_tail_volume + (
                            (
                                volume.user_lipid_volume - lipid.sample_lipid_name.project_user_lipid_name.hg_volume
                            ) * lipid.lipid_mol_fraction
                        )
                    break

    for lipid in sample_lipids_out:
        if lipid.sample_lipid_name.project_lipid_name:
            out_combined_head_volume = out_combined_head_volume + (
                    lipid.sample_lipid_name.project_lipid_name.hg_volume * lipid.lipid_mol_fraction
                )
            out_combined_tail_volume = out_combined_tail_volume + (
                    (
                        eval(lipid.sample_lipid_name.project_lipid_name.total_volume_equation) - lipid.sample_lipid_name.project_lipid_name.hg_volume
                    )
                    * lipid.lipid_mol_fraction
                )
        else:
            for volume in volumes:
                if volume.project_user_lipid_name == lipid.sample_lipid_name.project_user_lipid_name:
                    out_combined_head_volume = out_combined_head_volume + (
                            lipid.sample_lipid_name.project_user_lipid_name.hg_volume * lipid.lipid_mol_fraction
                        )
                    out_combined_tail_volume = out_combined_tail_volume + (
                            (
                                volume.user_lipid_volume - lipid.sample_lipid_name.project_user_lipid_name.hg_volume
                            ) * lipid.lipid_mol_fraction
                        )
                    break

    in_combined_head_volume = round(in_combined_head_volume, 2)
    in_combined_tail_volume = round(in_combined_tail_volume, 2)

    out_combined_head_volume = round(out_combined_head_volume, 2)
    out_combined_tail_volume = round(out_combined_tail_volume, 2)

    if request.method == 'POST':
        form = Asymmetrical_Parameter_Form(request.POST)
        if form.is_valid():
            parameters = form.save(commit=False)

            # overall
            parameters.sample_title = sample

            # Math'ed values
            # Head
            parameters.in_headgroup_volume = in_combined_head_volume
            if in_combined_head_volume == 0:
                parameters.in_headgroup_volume_upperbound = 1
                parameters.in_headgroup_volume_lowerbound = -1
            else:
                parameters.in_headgroup_volume_upperbound = (abs(in_combined_head_volume)*1.20)
                parameters.in_headgroup_volume_lowerbound = -(abs(in_combined_head_volume)*1.20)

            parameters.out_headgroup_volume = out_combined_head_volume
            if out_combined_head_volume == 0:
                parameters.out_headgroup_volume_upperbound = 1
                parameters.out_headgroup_volume_lowerbound = -1
            else:
                parameters.out_headgroup_volume_upperbound = (abs(out_combined_head_volume)*1.20)
                parameters.out_headgroup_volume_lowerbound = -(abs(out_combined_head_volume)*1.20)

            # Chain
            parameters.in_chain_volume = in_combined_tail_volume
            if in_combined_tail_volume == 0:
                parameters.in_chain_volume_upperbound = 1
                parameters.in_chain_volume_lowerbound = -1
            else:
                parameters.in_chain_volume_upperbound = (abs(in_combined_tail_volume)*1.10)
                parameters.in_chain_volume_lowerbound = -(abs(in_combined_tail_volume)*1.10)

            parameters.out_chain_volume = out_combined_tail_volume
            if out_combined_tail_volume == 0:
                parameters.out_chain_volume_upperbound = 1
                parameters.out_chain_volume_lowerbound = -1
            else:
                parameters.out_chain_volume_upperbound = (abs(out_combined_tail_volume)*1.10)
                parameters.out_chain_volume_lowerbound = -(abs(out_combined_tail_volume)*1.10)

            parameters.save()
            return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)
    else:
        form = Asymmetrical_Parameter_Form()

    return render(request, 'viewer/form.html', {'project':project, 'sample':sample, 'form': form})

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

    return render(request, 'viewer/form.html', {'project':project, 'sample':sample, 'parameters':parameters, 'form': form})

# Delete parameters
def asymmetrical_parameter_delete_warning(request, project_id, sample_id, parameter_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    parameter = get_object_or_404(Asymmetrical_Parameters, id=parameter_id)
    
    if request.method == 'POST':
        parameter.delete()
        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)

    return render(request, 'viewer/delete_warning.html', {'project':project, 'sample':sample, 'parameters':parameter})

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
        data_upload_form = Data_Upload_Form(sample.id, request.POST, request.FILES)
        if data_upload_form.is_valid():
            data_info = data_upload_form.save(commit=False)
            data_info.sample_title = sample

            data_file = request.FILES["data_file"]

            file_data = data_file.read().decode("utf-8")

            lines = file_data.split("\n")

            #loop over the lines and save them in db. If error, store as string and then display
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
        data_upload_form = Data_Upload_Form(sample.id)

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
        form = Data_Edit_Form(request.POST, instance=data)
        if form.is_valid():
            data = form.save(commit=False)
            data.save()
            return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)
    else:
        form = Data_Edit_Form(instance=data)

    return render(request, 'viewer/form.html', {'project':project, 'sample':sample, 'data':data, 'form': form})

# Delete data
def data_delete_warning(request, project_id, sample_id, data_id):
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)
    data = get_object_or_404(Data_Set, id=data_id)
    
    if request.method == 'POST':
        data.delete()
        return redirect('viewer:sample_detail', project_id=project.id, sample_id=sample.id)

    return render(request, 'viewer/delete_warning.html', {'project':project, 'sample':sample, 'data':data})

## Fit
# Main fit page
def fit_main(request, project_id, sample_id, parameter_id):
## Tutorials
    if request.user.is_anonymous:
        xuser_tutorial = False
    else:
        xuser = ExtendedUser.objects.get(user=request.user)
        xuser_tutorial = xuser.display_tutorial

    tutorial = True
    # Dismiss the tutorial
    if "dismiss_this" in request.POST:
        tutorial = False

    # Dismiss all tutorials
    if "dismiss_all" in request.POST:
        xuser.display_tutorial = False
        xuser.save()

# Overall
    project = get_object_or_404(Project, id=project_id)
    sample = get_object_or_404(Sample, id=sample_id)

    # Data
    datas = Data_Set.objects.filter(sample_title_id=sample_id)
    data_exists = True if datas else False

    # with the given project and sample, generate a fit main and return it. 
    #generate_fit_main(project, sample, datas)

    if project.model_type == "SM":
        sample_lipids = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='BOTH')
        parameter = get_object_or_404(Symmetrical_Parameters, id=parameter_id)
    elif project.model_type == "AS":
        sample_lipids_in = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='IN')
        sample_lipids_out = Sample_Lipid.objects.filter(sample_title_id=sample_id, lipid_location='OUT')
        parameter = get_object_or_404(Asymmetrical_Parameters, id=parameter_id)

# Check for 0 value parameters
    zero_parameter = False
    if project.model_type == "SM":
        if parameter.chain_volume == 0  or parameter.headgroup_volume == 0 or parameter.terminal_methyl_volume == 0 or parameter.lipid_area == 0 or parameter.sigma == 0:
            zero_parameter = True
    if project.model_type == "AS":
        if parameter.in_chain_volume == 0  or parameter.in_headgroup_volume == 0 or parameter.in_terminal_methyl_volume == 0 or parameter.in_lipid_area == 0 or parameter.out_chain_volume == 0  or parameter.out_headgroup_volume == 0 or parameter.out_terminal_methyl_volume == 0 or parameter.out_lipid_area == 0 or parameter.sigma == 0:
            zero_parameter = True

    xray_datas = datas.filter(data_type='XR')
    neutron_datas = datas.filter(data_type='NU')

# Decalre
    now = timezone.now()
    fit_result = None
    show_statistics = False
    show_probabilities = False

## Forms
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
                # Get scale values
                scale_value = xray_scale_form.cleaned_data['scale']

                if xray_data.scale_upperbound <= scale_value:
                    xray_data.scale_upperbound = (abs(scale_value)*1.5)
                elif xray_data.scale_lowerbound >= scale_value:
                    xray_data.scale_lowerbound = -(abs(scale_value)*1.5)
                
                # Get range values
                max_value = xray_range_form.cleaned_data['max_value']
                min_value = xray_range_form.cleaned_data['min_value']

                try:
                    # Find the indexes for the closest value in q_value
                    max_index = min(enumerate(xray_data.q_value), key=lambda x: abs(max_value - x[1]))
                    min_index = min(enumerate(xray_data.q_value), key=lambda x: abs(min_value - x[1]))

                    # Set the indexes in the db - +1 to take the closest on the inside for the low end
                    xray_data.max_index = max_index[0]
                    xray_data.min_index = min_index[0] + 1

                    xray_data.save()
                except TypeError:
                    xray_data.save()

                return redirect('viewer:fit_main', project_id=project.id, sample_id=sample.id, parameter_id=parameter.id)

        else:
            xray_range_form = Data_Range_Form()
            xray_scale_form = Data_Scale_Form(instance=xray_data)

        xray_ranges.append(xray_range_form)
        xray_scales.append(xray_scale_form)

    neutron_ranges = []
    neutron_scales = []

    # Update q range for all neutron datasets
    for neutron_data in neutron_datas:
        if neutron_data.data_set_title in request.POST:
            neutron_range_form = Data_Range_Form(request.POST)
            neutron_scale_form = Data_Scale_Form(request.POST, instance=neutron_data)
            if neutron_range_form.is_valid() and neutron_scale_form.is_valid():
                # Get scale values
                scale_value = neutron_scale_form.cleaned_data['scale']

                if neutron_data.scale_upperbound <= scale_value:
                    neutron_data.scale_upperbound = (abs(scale_value)*1.5)
                elif neutron_data.scale_lowerbound >= scale_value:
                    neutron_data.scale_lowerbound = -(abs(scale_value)*1.5)

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

                return redirect('viewer:fit_main', project_id=project.id, sample_id=sample.id, parameter_id=parameter.id)

        else:
            neutron_range_form = Data_Range_Form()
            neutron_scale_form = Data_Scale_Form(instance=neutron_data)

        neutron_ranges.append(neutron_range_form)
        neutron_scales.append(neutron_scale_form)

# Do the fit
    if project.model_type == "SM":
        if "fit" in request.POST:
            # Do fit
            fit_result = symmetrical_fit(parameter, sample_lipids, datas, project.system_tempurature, project.advanced_options)
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
            fit_result = asymmetrical_fit(parameter, sample_lipids_in, sample_lipids_out, datas, project.system_tempurature, project.advanced_options)
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

    # Show stats / probabilities / graph
    if "statistics" in request.POST:
        show_statistics = True
        show_probabilities = False

    if "probabilities" in request.POST:
        show_probabilities = True
        show_statistics = False

    if "graphs" in request.POST:
        show_statistics = False
        show_probabilities = False

    ## GRAPHS
    xray_figures = []
    neutron_figures = []

# X-Ray fit graphs
    for xray_data in xray_datas:
        xray_fig = plt.figure(figsize=(5.5,4.3))

        x = np.asarray(xray_data.q_value[xray_data.min_index:xray_data.max_index])
        y = np.asarray(xray_data.intensity_value[xray_data.min_index:xray_data.max_index])
        error = np.asarray(xray_data.error_value[xray_data.min_index:xray_data.max_index])

        # Data scatter plot
        plt.errorbar(
            x, 
            y,
            fmt='o',
            color='c',
            mfc='w',
            zorder=1
        )
        plt.errorbar(x, y-error, fmt='_', color='grey', zorder=0)
        plt.errorbar(x, y+error, fmt='_', color='grey', zorder=0)

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
                    symmetrical_graph(parameter, sample_lipids, xray_data, project.system_tempurature, project.advanced_options),
                    color='r',
                    label='Best Fit',
                    zorder=2
                )
            elif project.model_type == "AS":
                plt.plot(
                    xray_data.q_value[xray_data.min_index:xray_data.max_index],
                    asymmetrical_graph(parameter, sample_lipids_in, sample_lipids_out, xray_data, project.system_tempurature, project.advanced_options),
                    color='r',
                    label='Best Fit',
                    zorder=2
                )

        xray_figures.append(mpld3.fig_to_html(xray_fig))
        plt.cla()

# Neutron fit graphs
    for neutron_data in neutron_datas:
        neutron_fig = plt.figure(figsize=(5.5,4.3))

        x = np.asarray(neutron_data.q_value[neutron_data.min_index:neutron_data.max_index])
        y = np.asarray(neutron_data.intensity_value[neutron_data.min_index:neutron_data.max_index])
        error = np.asarray(neutron_data.error_value[neutron_data.min_index:neutron_data.max_index])

        # Data scatter plot
        plt.errorbar(
            x,
            y,
            fmt='o',
            color='c',
            mfc='w',
            zorder=1
        )
        plt.errorbar(x, y-error, fmt='_', color='grey', zorder=0)
        plt.errorbar(x, y+error, fmt='_', color='grey', zorder=0)

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
                    symmetrical_graph(parameter, sample_lipids, neutron_data, project.system_tempurature, project.advanced_options),
                    color='r',
                    label='Best Fit',
                    zorder=2
                )
            elif project.model_type == "AS":
                plt.plot(
                    neutron_data.q_value[neutron_data.min_index:neutron_data.max_index],
                    asymmetrical_graph(parameter, sample_lipids_in, sample_lipids_out, neutron_data, project.system_tempurature, project.advanced_options),
                    color='r',
                    label='Best Fit',
                    zorder=2
                )

        neutron_figures.append(mpld3.fig_to_html(neutron_fig))
        plt.cla()

# Probability graphs
    prob_fig = plt.figure(figsize=(6,5))
    xray_sdp_graphs = []
    xray_sdp_data = {}
    neutron_sdp_graphs = []
    neutron_sdp_data = {}

    # Symmetrical - combined halfs
    if project.model_type == "SM":
        x_values = np.arange(-40, 40, 0.2)
        
        # Calculate probabilities
        head_prob = head(
                parameter.chain_volume,
                parameter.headgroup_volume,
                parameter.lipid_area,
                parameter.headgroup_thickness,
                parameter.sigma,
                x_values
            )
        chain_prob = chain(
                parameter.chain_volume,
                parameter.lipid_area,
                parameter.sigma,
                x_values
            )
        tm_prob = terminal(
                parameter.terminal_methyl_volume,
                parameter.lipid_area,
                parameter.sigma,
                x_values
            )
        methylene_prob = methylene(
                parameter.chain_volume,
                parameter.terminal_methyl_volume,
                parameter.lipid_area,
                parameter.sigma,
                x_values
            )
        water_prob = water(
                parameter.chain_volume,
                parameter.headgroup_volume,
                parameter.lipid_area,
                parameter.headgroup_thickness,
                parameter.sigma,
                x_values
            )

        if zero_parameter:
            plt.plot()
            plt.title('! DIVIDE BY ZERO ERROR !')
        else:
            # headgroup
            plt.plot(
                x_values,
                head_prob,
                color='c',
                marker='.',
                markersize='5',
                label='Headgroup',
                zorder=0
            )
            # chain
            plt.plot(
                x_values,
                chain_prob,
                color='g',
                marker='v',
                markersize='5',
                label='Chains',
                zorder=1
            )
            # terminal methyl
            plt.plot(
                x_values,
                tm_prob,
                color='m',
                marker='s',
                markersize='5',
                label='Terminal Methyl',
                zorder=2
            )
            # methylene
            plt.plot(
                x_values,
                methylene_prob,
                color='k',
                marker='p',
                markersize='5',
                label='Methylene',
                zorder=3
            )
            # water
            plt.plot(
                x_values,
                water_prob,
                color='b',
                marker='x',
                markersize='5',
                label='Water',
                zorder=4
            )

            plt.legend(loc=1)
            plt.xlabel('Distance from bilayer center [Å]')
            plt.ylabel('Volume probability')

        for xray_data in xray_datas:
            xray_sdp_fig = plt.figure(figsize=(5.5,4.3))

            sdp_results = symmetrical_sdp(
                    parameter,
                    head_prob,
                    methylene_prob,
                    tm_prob,
                    water_prob,
                    sample_lipids,
                    xray_data,
                    project.system_tempurature,
                    project.advanced_options
                )

            additional_parameters = sym_additional_parameters(
                    parameter,
                    sample_lipids,
                    xray_data,
                    project.system_tempurature,
                    np.asarray(x_values),
                    np.asarray(head_prob),
                    project.advanced_options
                )

            xray_sdp_data[xray_data] = sdp_results

            if zero_parameter:
                plt.plot()
                plt.title(str(xray_data.data_set_title)+' ! DIVIDE BY ZERO ERROR !')
            else:
                plt.plot(
                    x_values,
                    sdp_results[0],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    label='Combined SDP',
                    zorder=5
                )
                # headgroup
                plt.plot(
                    x_values,
                    sdp_results[1],
                    color='c',
                    marker='.',
                    markersize='5',
                    label='Headgroup',
                    zorder=0
                )
                # methyl
                plt.plot(
                    x_values,
                    sdp_results[2],
                    color='g',
                    marker='v',
                    markersize='5',
                    label='Methylene',
                    zorder=1
                )
                # terminal methyl
                plt.plot(
                    x_values,
                    sdp_results[3],
                    color='m',
                    marker='s',
                    markersize='5',
                    label='Terminal Methyl',
                    zorder=2
                )
                # water
                plt.plot(
                    x_values,
                    sdp_results[4],
                    color='b',
                    marker='x',
                    markersize='5',
                    label='Water',
                    zorder=4
                )

                plt.legend(loc=1)
                plt.xlabel('Distance from bilayer center [Å]')
                plt.ylabel('ED (e Å-3 )')

                plt.title(xray_data.data_set_title)

            xray_sdp_graphs.append(mpld3.fig_to_html(xray_sdp_fig))
            plt.cla()

        for neutron_data in neutron_datas:
            neutron_sdp_fig = plt.figure(figsize=(5.5,4.3))

            sdp_results = symmetrical_sdp(
                    parameter,
                    head_prob,
                    methylene_prob,
                    tm_prob,
                    water_prob,
                    sample_lipids,
                    neutron_data,
                    project.system_tempurature,
                    project.advanced_options
                )

            additional_parameters = sym_additional_parameters(
                    parameter,
                    sample_lipids,
                    neutron_data,
                    project.system_tempurature,
                    np.asarray(x_values),
                    np.asarray(head_prob),
                    project.advanced_options
                )

            neutron_sdp_data[neutron_data] = sdp_results

            if zero_parameter:
                plt.plot()
                plt.title(str(neutron_data.data_set_title)+' ! DIVIDE BY ZERO ERROR !')
            else:
                plt.plot(
                    x_values,
                    sdp_results[0],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    label='Combined SDP',
                    zorder=5
                )
                # headgroup
                plt.plot(
                    x_values,
                    sdp_results[1],
                    color='c',
                    marker='.',
                    markersize='5',
                    label='Headgroup',
                    zorder=0
                )
                # methyl
                plt.plot(
                    x_values,
                    sdp_results[2],
                    color='g',
                    marker='v',
                    markersize='5',
                    label='Methylene',
                    zorder=1
                )
                # terminal methyl
                plt.plot(
                    x_values,
                    sdp_results[3],
                    color='m',
                    marker='s',
                    markersize='5',
                    label='Terminal Methyl',
                    zorder=2
                )
                # water
                plt.plot(
                    x_values,
                    sdp_results[4],
                    color='b',
                    marker='x',
                    markersize='5',
                    label='Water',
                    zorder=4
                )

                plt.legend(loc=1)
                plt.xlabel('Distance from bilayer center [Å]')
                plt.ylabel('NSLD (Å-3 x 10-5)')

                plt.title(neutron_data.data_set_title)

            neutron_sdp_graphs.append(mpld3.fig_to_html(neutron_sdp_fig))
            plt.cla()

    # Asymmetrical - separate halfs
    if project.model_type == "AS":
        in_x_values = np.arange(-40, 0.2, 0.2)
        out_x_values = np.arange(-0.2, 40, 0.2)

        # Calculate probabilities
        in_head_prob = head(
                parameter.in_chain_volume,
                parameter.in_headgroup_volume,
                parameter.in_lipid_area,
                parameter.in_headgroup_thickness,
                parameter.sigma,
                in_x_values
            )
        out_head_prob = head(
                parameter.out_chain_volume,
                parameter.out_headgroup_volume,
                parameter.out_lipid_area,
                parameter.out_headgroup_thickness,
                parameter.sigma,
                out_x_values
            )
        in_chain_prob = chain(
                parameter.in_chain_volume,
                parameter.in_lipid_area,
                parameter.sigma,
                in_x_values
            )
        out_chain_prob = chain(
                parameter.out_chain_volume,
                parameter.out_lipid_area,
                parameter.sigma,
                out_x_values
            )
        in_tm_prob = terminal(
                parameter.in_terminal_methyl_volume,
                parameter.in_lipid_area,
                parameter.sigma,
                in_x_values
            )
        out_tm_prob = terminal(
                parameter.out_terminal_methyl_volume,
                parameter.out_lipid_area,
                parameter.sigma,
                out_x_values
            )
        in_methylene_prob = methylene(
                parameter.in_chain_volume,
                parameter.in_terminal_methyl_volume,
                parameter.in_lipid_area,
                parameter.sigma,
                in_x_values
            )
        out_methylene_prob = methylene(
                parameter.out_chain_volume,
                parameter.out_terminal_methyl_volume,
                parameter.out_lipid_area,
                parameter.sigma,
                out_x_values
            )
        in_water_prob = water(
                parameter.in_chain_volume,
                parameter.in_headgroup_volume,
                parameter.in_lipid_area,
                parameter.in_headgroup_thickness,
                parameter.sigma,
                in_x_values
            )
        out_water_prob = water(
                parameter.out_chain_volume,
                parameter.out_headgroup_volume,
                parameter.out_lipid_area,
                parameter.out_headgroup_thickness,
                parameter.sigma,
                out_x_values
            )

        if zero_parameter:
            plt.plot()
            plt.title('! DIVIDE BY ZERO ERROR !')
        else:
            # in headgroup
            plt.plot(
                in_x_values,
                in_head_prob,
                color='c',
                marker='.',
                markersize='5',
                label='Headgroup',
                zorder=0
            )
            # out headgroup
            plt.plot(
                out_x_values,
                out_head_prob,
                color='c',
                marker='.',
                markersize='5',
                zorder=0
            )
            # in chain
            plt.plot(
                in_x_values,
                in_chain_prob,
                color='g',
                marker='v',
                markersize='5',
                label='Chains',
                zorder=1
            )
            # out chain
            plt.plot(
                out_x_values,
                out_chain_prob,
                color='g',
                marker='v',
                markersize='5',
                zorder=1
            )
            # in terminal methyl
            plt.plot(
                in_x_values,
                in_tm_prob,
                color='m',
                marker='s',
                markersize='5',
                label='Terminal Methyl',
                zorder=2
            )
            # out terminal methyl
            plt.plot(
                out_x_values,
                out_tm_prob,
                color='m',
                marker='s',
                markersize='5',
                zorder=2
            )
            # in methylene
            plt.plot(
                in_x_values,
                in_methylene_prob,
                color='k',
                marker='p',
                markersize='5',
                label='Methylene',
                zorder=3
            )
            # out methylene
            plt.plot(
                out_x_values,
                out_methylene_prob,
                color='k',
                marker='p',
                markersize='5',
                zorder=3
            )
            # in water
            plt.plot(
                in_x_values,
                in_water_prob,
                color='b',
                marker='x',
                markersize='5',
                label='Water',
                zorder=4
            )
            # out water
            plt.plot(
                out_x_values,
                out_water_prob,
                color='b',
                marker='x',
                markersize='5',
                zorder=4
            )

            plt.legend(loc=1)
            plt.xlabel('Distance from bilayer center [Å]')
            plt.ylabel('Volume probability')

        for xray_data in xray_datas:
            xray_sdp_fig = plt.figure(figsize=(5.5,4.3))

            sdp_results = asymmetrical_sdp(
                    parameter,
                    in_head_prob,
                    in_methylene_prob,
                    in_tm_prob,
                    in_water_prob,
                    out_head_prob,
                    out_methylene_prob,
                    out_tm_prob,
                    out_water_prob,
                    sample_lipids_in,
                    sample_lipids_out, 
                    xray_data,
                    project.system_tempurature,
                    project.advanced_options
                )
            
            additional_parameters = asym_additional_parameters(
                    parameter,
                    sample_lipids_in,
                    sample_lipids_out, 
                    xray_data,
                    project.system_tempurature,
                    np.asarray(in_head_prob),
                    np.asarray(out_head_prob),
                    in_x_values,
                    out_x_values,
                    project.advanced_options
                )

            xray_sdp_data[xray_data] = sdp_results

            if zero_parameter:
                plt.plot()
                plt.title(str(xray_data.data_set_title)+' ! DIVIDE BY ZERO ERROR !')
            else:
                plt.plot(
                    in_x_values,
                    sdp_results[0],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    label='Combined SDP',
                    zorder=5
                )
                plt.plot(
                    out_x_values,
                    sdp_results[1],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    zorder=5
                )
                # headgroup
                plt.plot(
                    in_x_values,
                    sdp_results[2],
                    color='c',
                    marker='.',
                    markersize='5',
                    label='Headgroup',
                    zorder=0
                )
                plt.plot(
                    out_x_values,
                    sdp_results[3],
                    color='c',
                    marker='.',
                    markersize='5',
                    zorder=0
                )
                # methyl
                plt.plot(
                    in_x_values,
                    sdp_results[4],
                    color='g',
                    marker='v',
                    markersize='5',
                    label='Methylene',
                    zorder=1
                )
                plt.plot(
                    out_x_values,
                    sdp_results[5],
                    color='g',
                    marker='v',
                    markersize='5',
                    zorder=1
                )
                # terminal methyl
                plt.plot(
                    in_x_values,
                    sdp_results[6],
                    color='m',
                    marker='s',
                    markersize='5',
                    label='Terminal Methyl',
                    zorder=2
                )
                plt.plot(
                    out_x_values,
                    sdp_results[7],
                    color='m',
                    marker='s',
                    markersize='5',
                    zorder=2
                )
                # water
                plt.plot(
                    in_x_values,
                    sdp_results[8],
                    color='b',
                    marker='x',
                    markersize='5',
                    label='Water',
                    zorder=4
                )
                plt.plot(
                    out_x_values,
                    sdp_results[9],
                    color='b',
                    marker='x',
                    markersize='5',
                    zorder=4
                )

                plt.legend(loc=1)
                plt.xlabel('Distance from bilayer center [Å]')
                plt.ylabel('ED (e Å-3 )')

                plt.title(xray_data.data_set_title)

            xray_sdp_graphs.append(mpld3.fig_to_html(xray_sdp_fig))
            plt.cla()

        for neutron_data in neutron_datas:
            neutron_sdp_fig = plt.figure(figsize=(5.5,4.3))

            sdp_results = asymmetrical_sdp(
                    parameter,
                    in_head_prob,
                    in_methylene_prob,
                    in_tm_prob,
                    in_water_prob,
                    out_head_prob,
                    out_methylene_prob,
                    out_tm_prob,
                    out_water_prob,
                    sample_lipids_in,
                    sample_lipids_out, 
                    neutron_data,
                    project.system_tempurature,
                    project.advanced_options
                )

            additional_parameters = asym_additional_parameters(
                    parameter,
                    sample_lipids_in,
                    sample_lipids_out, 
                    neutron_data,
                    project.system_tempurature,
                    np.asarray(in_head_prob),
                    np.asarray(out_head_prob),
                    in_x_values,
                    out_x_values,
                    project.advanced_options
                )

            neutron_sdp_data[neutron_data] = sdp_results

            if zero_parameter:
                plt.plot()
                plt.title(str(neutron_data.data_set_title)+' ! DIVIDE BY ZERO ERROR !')
            else:
                plt.plot(
                    in_x_values,
                    sdp_results[0],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    label='Combined SDP',
                    zorder=5
                )
                plt.plot(
                    out_x_values,
                    sdp_results[1],
                    linewidth=4,
                    color='k',
                    markersize='5',
                    zorder=5
                )
                # headgroup
                plt.plot(
                    in_x_values,
                    sdp_results[2],
                    color='c',
                    marker='.',
                    markersize='5',
                    label='Headgroup',
                    zorder=0
                )
                plt.plot(
                    out_x_values,
                    sdp_results[3],
                    color='c',
                    marker='.',
                    markersize='5',
                    zorder=0
                )
                # methyl
                plt.plot(
                    in_x_values,
                    sdp_results[4],
                    color='g',
                    marker='v',
                    markersize='5',
                    label='Methylene',
                    zorder=1
                )
                plt.plot(
                    out_x_values,
                    sdp_results[5],
                    color='g',
                    marker='v',
                    markersize='5',
                    zorder=1
                )
                # terminal methyl
                plt.plot(
                    in_x_values,
                    sdp_results[6],
                    color='m',
                    marker='s',
                    markersize='5',
                    label='Terminal Methyl',
                    zorder=2
                )
                plt.plot(
                    out_x_values,
                    sdp_results[7],
                    color='m',
                    marker='s',
                    markersize='5',
                    zorder=2
                )
                # water
                plt.plot(
                    in_x_values,
                    sdp_results[8],
                    color='b',
                    marker='x',
                    markersize='5',
                    label='Water',
                    zorder=4
                )
                plt.plot(
                    out_x_values,
                    sdp_results[9],
                    color='b',
                    marker='x',
                    markersize='5',
                    zorder=4
                )

                plt.legend(loc=1)
                plt.xlabel('Distance from bilayer center [Å]')
                plt.ylabel('NSLD (Å-3 x 10-5)')

                plt.title(neutron_data.data_set_title)

            neutron_sdp_graphs.append(mpld3.fig_to_html(neutron_sdp_fig))
            plt.cla()

    prob_graph = mpld3.fig_to_html(prob_fig)

# Fit data download
    if "fit_download" in request.POST:
        # Filename
        file_name = str(project.project_title).replace(' ','-').replace(':','.')+'_FIT_download_'+'_'+str(sample.sample_title).replace(' ','-').replace(':','.')+'_'+str(parameter.name).replace(' ','-').replace(':','.')+now.strftime("%m-%d-%H.%M")+'.csv'

        # Create the HttpResponse object with the appropriate CSV header.
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename={0}'.format(file_name)

        writer = csv.writer(response)
        writer.writerow(['VesicleViewer Fit output', now])
        writer.writerow(['Project Name', 'Sample Name', 'Parameter Set'])
        writer.writerow([project.project_title, sample.sample_title, parameter.name])
        writer.writerow([])

        if project.model_type == "SM":
            writer.writerow(['Calculated Parameters'])
            writer.writerow(['Db', additional_parameters[0]])
            writer.writerow(['2Dc', additional_parameters[1]])
            writer.writerow(['Dhh', additional_parameters[2]])
            writer.writerow(['Dh', parameter.headgroup_thickness])
            writer.writerow(['Al', parameter.lipid_area])
            writer.writerow([])
        elif project.model_type == "AS":
            writer.writerow(['Calculated Parameters'])
            writer.writerow([])

            writer.writerow(['Inner'])
            writer.writerow(['Db', additional_parameters[0]])
            writer.writerow(['2Dc', additional_parameters[1]])
            writer.writerow(['Dhh', additional_parameters[4]])
            writer.writerow(['Dh', parameter.in_headgroup_thickness])
            writer.writerow(['Al', parameter.in_lipid_area])
            writer.writerow([])

            writer.writerow(['Outter'])
            writer.writerow(['Db', additional_parameters[2]])
            writer.writerow(['2Dc', additional_parameters[3]])
            writer.writerow(['Dhh', additional_parameters[5]])
            writer.writerow(['Dh', parameter.out_headgroup_thickness])
            writer.writerow(['Al', parameter.out_lipid_area])
            writer.writerow([])

        writer.writerow({'Fit Statistics'})
        for line in parameter.fit_report:
            writer.writerow([line])

        writer.writerow([])

        writer.writerow(['Q', 'Experimental i', 'Experimental Error', 'Calculated i'])

        for xray_data in xray_datas:
            writer.writerow([xray_data.data_set_title])
            if project.model_type == "SM":
                calculated_i_values = symmetrical_graph(parameter, sample_lipids, xray_data, project.system_tempurature, project.advanced_options)
            elif project.model_type == "AS":
                calculated_i_values = asymmetrical_graph(parameter, sample_lipids_in, sample_lipids_out, xray_data, project.system_tempurature, project.advanced_options)

            j = 0
            for i in range(xray_data.min_index, xray_data.max_index):
                writer.writerow([xray_data.q_value[i], xray_data.intensity_value[i], xray_data.error_value[i], calculated_i_values[j]])
                j = j+1

        for neutron_data in neutron_datas:
            writer.writerow([neutron_data.data_set_title])
            if project.model_type == "SM":
                calculated_i_values = symmetrical_graph(parameter, sample_lipids, neutron_data, project.system_tempurature, project.advanced_options)
            elif project.model_type == "AS":
                calculated_i_values = asymmetrical_graph(parameter, sample_lipids_in, sample_lipids_out, neutron_data, project.system_tempurature, project.advanced_options)

            j = 0
            for i in range(neutron_data.min_index, neutron_data.max_index):
                writer.writerow([neutron_data.q_value[i], neutron_data.intensity_value[i], neutron_data.error_value[i], calculated_i_values[j]])
                j = j+1

        return response

# SDP data download
    if "sdp_download" in request.POST:
        # Filename
        file_name = str(project.project_title).replace(' ','-').replace(':','.')+'_SDP_download_'+'_'+str(sample.sample_title).replace(' ','-').replace(':','.')+'_'+str(parameter.name).replace(' ','-').replace(':','.')+now.strftime("%m-%d-%H.%M")+'.csv'

        # Create the HttpResponse object with the appropriate CSV header.
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename={0}'.format(file_name)

        writer = csv.writer(response)
        writer.writerow(['VesicleViewer SDP output', now])
        writer.writerow(['Project Name', 'Sample Name', 'Parameter Set'])
        writer.writerow([project.project_title, sample.sample_title, parameter.name])
        writer.writerow([])

        # Water probablilities
        writer.writerow(['Volume Probabilities'])

        if project.model_type == "SM":
            writer.writerow([])
            writer.writerow(['Headgroup'])
            writer.writerow(['z', 'Ph(z)'])
            for z, ph in zip (x_values, head_prob):
                writer.writerow([z, ph])

            writer.writerow([])
            writer.writerow(['Chains'])
            writer.writerow(['z', 'Phc(z)'])
            for z, pc in zip (x_values, chain_prob):
                writer.writerow([z, pc])

            writer.writerow([])
            writer.writerow(['Terminal Methyl'])
            writer.writerow(['z', 'Ptm(z)'])
            for z, ptm in zip (x_values, tm_prob):
                writer.writerow([z, ptm])

            writer.writerow([])
            writer.writerow(['Methylene'])
            writer.writerow(['z', 'Pch(z)'])
            for z, pm in zip (x_values, methylene_prob):
                writer.writerow([z, pm])

            writer.writerow([])
            writer.writerow(['Water'])
            writer.writerow(['z', 'Pw(z)'])
            for z, pw in zip (x_values, water_prob):
                writer.writerow([z, pw])

        if project.model_type == "AS":
            writer.writerow([])
            writer.writerow(['Headgroup'])
            writer.writerow(['z', 'Ph(z)'])
            writer.writerow([])
            writer.writerow(['INNER'])
            for z, ph in zip (in_x_values, in_head_prob):
                writer.writerow([z, ph])
            writer.writerow([])
            writer.writerow(['OUTTER'])
            for z, ph in zip (out_x_values, out_head_prob):
                writer.writerow([z, ph])

            writer.writerow([])
            writer.writerow(['Chains'])
            writer.writerow(['z', 'Phc(z)'])
            writer.writerow([])
            writer.writerow(['INNER'])
            for z, pc in zip (in_x_values, in_chain_prob):
                writer.writerow([z, pc])
            writer.writerow([])
            writer.writerow(['OUTTER'])
            for z, pc in zip (out_x_values, out_chain_prob):
                writer.writerow([z, pc])

            writer.writerow([])
            writer.writerow(['Terminal Methyl'])
            writer.writerow(['z', 'Ptm(z)'])
            writer.writerow([])
            writer.writerow(['INNER'])
            for z, ptm in zip (in_x_values, in_tm_prob):
                writer.writerow([z, ptm])
            writer.writerow([])
            writer.writerow(['OUTTER'])
            for z, ptm in zip (out_x_values, out_tm_prob):
                writer.writerow([z, ptm])

            writer.writerow([])
            writer.writerow(['Methylene'])
            writer.writerow(['z', 'Pch(z)'])
            writer.writerow([])
            writer.writerow(['INNER'])
            for z, pm in zip (in_x_values, in_methylene_prob):
                writer.writerow([z, pm])
            writer.writerow([])
            writer.writerow(['OUTTER'])
            for z, pm in zip (out_x_values, out_methylene_prob):
                writer.writerow([z, pm])

            writer.writerow([])
            writer.writerow(['Water'])
            writer.writerow(['z', 'Pw(z)'])
            writer.writerow([])
            writer.writerow(['INNER'])
            for z, pw in zip (in_x_values, in_water_prob):
                writer.writerow([z, pw])
            writer.writerow([])
            writer.writerow(['OUTTER'])
            for z, pw in zip (out_x_values, out_water_prob):
                writer.writerow([z, pw])
            
        # SDP and Scaled Probabilities
        writer.writerow([])
        writer.writerow(['Scattering Density Profile'])

        for xray_data in xray_datas:
            sdp_results = xray_sdp_data[xray_data]

            writer.writerow([])
            writer.writerow([xray_data.data_set_title])
            if project.model_type == "SM":
                writer.writerow([])
                writer.writerow(['Combined SDP'])
                writer.writerow(['z', ''])
                for z, sdp in zip (x_values, sdp_results[0]):
                    writer.writerow([z, sdp])

                writer.writerow([])
                writer.writerow(['Headgroup'])
                writer.writerow(['z', ''])
                for z, sdph in zip (x_values, sdp_results[1]):
                    writer.writerow([z, sdph])

                writer.writerow([])
                writer.writerow(['Methylene'])
                writer.writerow(['z', ''])
                for z, sdpc in zip (x_values, sdp_results[2]):
                    writer.writerow([z, sdpc])

                writer.writerow([])
                writer.writerow(['Terminal Methyl'])
                writer.writerow(['z', ''])
                for z, sdptm in zip (x_values, sdp_results[3]):
                    writer.writerow([z, sdptm])

                writer.writerow([])
                writer.writerow(['Water'])
                writer.writerow(['z', ''])
                for z, sdpw in zip (x_values, sdp_results[4]):
                    writer.writerow([z, sdpw])

            if project.model_type == "AS":
                writer.writerow([])
                writer.writerow(['Combined SDP'])
                writer.writerow(['z', 'Pch(z)'])
                writer.writerow([])
                writer.writerow(['INNER'])
                for z, sdp in zip (in_x_values, sdp_results[0]):
                    writer.writerow([z, sdp])
                writer.writerow([])
                writer.writerow(['OUTTER'])
                for z, sdp in zip (out_x_values, sdp_results[1]):
                    writer.writerow([z, sdp])

                writer.writerow([])
                writer.writerow(['Headgroup'])
                writer.writerow(['z', 'Ph(z)'])
                writer.writerow([])
                writer.writerow(['INNER'])
                for z, sdph in zip (in_x_values, sdp_results[2]):
                    writer.writerow([z, sdph])
                writer.writerow([])
                writer.writerow(['OUTTER'])
                for z, sdph in zip (out_x_values, sdp_results[3]):
                    writer.writerow([z, sdph])

                writer.writerow([])
                writer.writerow(['Methylene'])
                writer.writerow(['z', 'Phc(z)'])
                writer.writerow([])
                writer.writerow(['INNER'])
                for z, sdpc in zip (in_x_values, sdp_results[4]):
                    writer.writerow([z, sdpc])
                writer.writerow([])
                writer.writerow(['OUTTER'])
                for z, sdpc in zip (out_x_values, sdp_results[5]):
                    writer.writerow([z, sdpc])

                writer.writerow([])
                writer.writerow(['Terminal Methyl'])
                writer.writerow(['z', 'Ptm(z)'])
                writer.writerow([])
                writer.writerow(['INNER'])
                for z, sdptm in zip (in_x_values, sdp_results[6]):
                    writer.writerow([z, sdptm])
                writer.writerow([])
                writer.writerow(['OUTTER'])
                for z, sdptm in zip (out_x_values, sdp_results[7]):
                    writer.writerow([z, sdptm])

                writer.writerow([])
                writer.writerow(['Water'])
                writer.writerow(['z', 'Pw(z)'])
                writer.writerow([])
                writer.writerow(['INNER'])
                for z, sdpw in zip (in_x_values, sdp_results[8]):
                    writer.writerow([z, sdpw])
                writer.writerow([])
                writer.writerow(['OUTTER'])
                for z, sdpw in zip (out_x_values, sdp_results[9]):
                    writer.writerow([z, sdpw])

        for neutron_data in neutron_datas:
            sdp_results = neutron_sdp_data[neutron_data]

            writer.writerow([])
            writer.writerow([neutron_data.data_set_title])
            if project.model_type == "SM":
                writer.writerow([])
                writer.writerow(['Combined SDP'])
                writer.writerow(['z', ''])
                for z, sdp in zip (x_values, sdp_results[0]):
                    writer.writerow([z, sdp])

                writer.writerow([])
                writer.writerow(['Headgroup'])
                writer.writerow(['z', ''])
                for z, sdph in zip (x_values, sdp_results[1]):
                    writer.writerow([z, sdph])

                writer.writerow([])
                writer.writerow(['Methylene'])
                writer.writerow(['z', ''])
                for z, sdpc in zip (x_values, sdp_results[2]):
                    writer.writerow([z, sdpc])

                writer.writerow([])
                writer.writerow(['Terminal Methyl'])
                writer.writerow(['z', ''])
                for z, sdptm in zip (x_values, sdp_results[3]):
                    writer.writerow([z, sdptm])

                writer.writerow([])
                writer.writerow(['Water'])
                writer.writerow(['z', ''])
                for z, sdpw in zip (x_values, sdp_results[4]):
                    writer.writerow([z, sdpw])

            if project.model_type == "AS":
                writer.writerow([])
                writer.writerow(['Combined SDP'])
                writer.writerow(['z', 'Pch(z)'])
                writer.writerow([])
                writer.writerow(['INNER'])
                for z, sdp in zip (in_x_values, sdp_results[0]):
                    writer.writerow([z, sdp])
                writer.writerow([])
                writer.writerow(['OUTTER'])
                for z, sdp in zip (out_x_values, sdp_results[1]):
                    writer.writerow([z, sdp])

                writer.writerow([])
                writer.writerow(['Headgroup'])
                writer.writerow(['z', 'Ph(z)'])
                writer.writerow([])
                writer.writerow(['INNER'])
                for z, sdph in zip (in_x_values, sdp_results[2]):
                    writer.writerow([z, sdph])
                writer.writerow([])
                writer.writerow(['OUTTER'])
                for z, sdph in zip (out_x_values, sdp_results[3]):
                    writer.writerow([z, sdph])

                writer.writerow([])
                writer.writerow(['Methylene'])
                writer.writerow(['z', 'Phc(z)'])
                writer.writerow([])
                writer.writerow(['INNER'])
                for z, sdpc in zip (in_x_values, sdp_results[4]):
                    writer.writerow([z, sdpc])
                writer.writerow([])
                writer.writerow(['OUTTER'])
                for z, sdpc in zip (out_x_values, sdp_results[5]):
                    writer.writerow([z, sdpc])

                writer.writerow([])
                writer.writerow(['Terminal Methyl'])
                writer.writerow(['z', 'Ptm(z)'])
                writer.writerow([])
                writer.writerow(['INNER'])
                for z, sdptm in zip (in_x_values, sdp_results[6]):
                    writer.writerow([z, sdptm])
                writer.writerow([])
                writer.writerow(['OUTTER'])
                for z, sdptm in zip (out_x_values, sdp_results[7]):
                    writer.writerow([z, sdptm])

                writer.writerow([])
                writer.writerow(['Water'])
                writer.writerow(['z', 'Pw(z)'])
                writer.writerow([])
                writer.writerow(['INNER'])
                for z, sdpw in zip (in_x_values, sdp_results[8]):
                    writer.writerow([z, sdpw])
                writer.writerow([])
                writer.writerow(['OUTTER'])
                for z, sdpw in zip (out_x_values, sdp_results[9]):
                    writer.writerow([z, sdpw])

        return response

    xray_graphs_and_forms = zip(xray_figures, xray_ranges, xray_scales, xray_datas)
    neutron_graphs_and_forms = zip(neutron_figures, neutron_ranges, neutron_scales, neutron_datas)

    plt.close('all')

## Done
    return render(request, 'viewer/fit_main.html', {
        'tutorial':tutorial,
        'xuser_tutorial':xuser_tutorial,
        'project':project,
        'sample':sample,
        'data_exists':data_exists,
        'parameter':parameter,
        'zero_parameter':zero_parameter,
        'parameter_update_form':parameter_update_form,
        'fit_result':fit_result,
        'show_stats':show_statistics,
        'show_probs':show_probabilities,
        'xray_graphs_and_forms':xray_graphs_and_forms,
        'neutron_graphs_and_forms':neutron_graphs_and_forms,
        'prob_graph':prob_graph,
        'xray_sdp_graphs':xray_sdp_graphs,
        'neutron_sdp_graphs':neutron_sdp_graphs,
        'additional_parameters':additional_parameters
    })
