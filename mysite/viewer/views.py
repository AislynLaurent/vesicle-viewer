## IMPORTS
# From django
from django.shortcuts import get_object_or_404, render
from django.shortcuts import redirect

# Other libraries
import re

# Models
from .models import *

# Forms
from .forms import *

# Other imports
from .symfit import *
from .asymfit import *
from .probabilities import *

# note: test if you need the above imports, as most of these imports are needed in fit.py as well
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
                    sample_title=sample,
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
                    sample_title=sample,
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
                    sample_title=sample,
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
                print(lipid_form)
                in_sample_lipid_name = lipid_form.cleaned_data['sample_lipid_name']
                in_lipid_mol_fraction = lipid_form.cleaned_data['lipid_mol_fraction']
                in_lipid_location = lipid_form.cleaned_data['location']

                existing_lipid, created_lipid = Sample_Lipid.objects.update_or_create(
                    sample_title=sample,
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
    return generate_fit_main(request, project_id, sample_id, parameter_id)