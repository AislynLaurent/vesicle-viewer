from django.urls import path
from . import views

app_name = 'viewer'

urlpatterns = [
    path('', views.index, name='index'),
    path('about', views.about, name='about'),
    path('help', views.get_help, name='help'),
    path('privacy', views.privacy, name='privacy'),

    ## Tutorials
    path('tutorials', views.enable_tutorials, name='enable_tutorials'),
    
    ## Models
    # Lipid
    path('lipid/<slug:lipid_name>', views.lipid_detail, name='lipid_detail'),
    # Molecule
    path('molecule/<slug:molecule_name>', views.molecule_detail, name='molecule_detail'),

    ## Pojects
    # New project
    path('project/new', views.project_new, name='project_new'),
    # List projects
    path('project/list', views.project_list, name='project_list'),
    # Project details
    path('project/<int:project_id>', views.project_detail, name='project_detail'),
    # Turn on advanced options for this project
    path('project/<int:project_id>/edit/', views.project_advanced_options, name='project_advanced_options'),
    # Edit a project
    path('project/<int:project_id>/edit/', views.project_edit, name='project_edit'),
    # Warning & delete project
    path('project/<int:project_id>/delete_warning/', views.project_delete_warning, name='project_delete_warning'),

    ## Project Lipids
    # Add project lipid
    path('project/<int:project_id>/lipid/add/', views.project_lipid_new, name='project_lipid_new'),
    # # Edit a project lipid
    # path('project/<int:project_id>/lipid/<int:lipid_id>/edit/', views.project_lipid_edit, name='project_lipid_edit'),
    # Warning & delete lipid
    path('project/<int:project_id>/lipid/<int:lipid_id>/delete_warning/', views.project_lipid_delete_warning, name='project_lipid_delete_warning'),

    ## Samples
    # New project
    path('project/<int:project_id>/sample/new', views.sample_new, name='sample_new'),
    # Sample details
    path('project/<int:project_id>/sample/<int:sample_id>', views.sample_detail, name='sample_detail'),
    # Edit a sample
    path('project/<int:project_id>/sample/<int:sample_id>/edit/', views.sample_edit, name='sample_edit'),
    # Warning & delete sample
    path('project/<int:project_id>/sample/<int:sample_id>/delete_warning/', views.sample_delete_warning, name='sample_delete_warning'),

    ## Sample Lipids
    # Add sample lipid
    path('project/<int:project_id>/sample/<int:sample_id>/lipid/add/', views.sample_lipid_new, name='sample_lipid_new'),
    # Edit a sample lipid
    path('project/<int:project_id>/sample/<int:sample_id>/lipid/<int:lipid_id>/edit/', views.sample_lipid_edit, name='sample_lipid_edit'),
    # Edit custom augments
    path('project/<int:project_id>/sample/<int:sample_id>/lipid/<int:lipid_id>/custom/', views.sample_custom_lipid_edit, name='sample_custom_lipid_edit'),
    # Warning & delete lipid
    path('project/<int:project_id>/sample/<int:sample_id>/lipid/<int:lipid_id>/delete_warning/', views.sample_lipid_delete_warning, name='sample_lipid_delete_warning'),

    ## Sym Parameters
    # Add project parameters
    path('project/<int:project_id>/sample/<int:sample_id>/symparameters/new/', views.symmetrical_parameters_new, name='symmetrical_parameters_new'),
    # Edit a projects parameters
    path('project/<int:project_id>/sample/<int:sample_id>/symparameters/<int:parameter_id>/edit/', views.symmetrical_parameters_edit, name='symmetrical_parameters_edit'),
    # Warning & delete parameters
    path('project/<int:project_id>/sample/<int:sample_id>/symparameters/<int:parameter_id>/delete_warning/', views.symmetrical_parameter_delete_warning, name='symmetrical_parameter_delete_warning'),

    ## Asym Parameters
    # Add project parameters
    path('project/<int:project_id>/sample/<int:sample_id>/asymparameters/new/', views.asymmetrical_parameters_new, name='asymmetrical_parameters_new'),
    # Edit a projects parameters
    path('project/<int:project_id>/sample/<int:sample_id>/asymparameters/<int:parameter_id>/edit/', views.asymmetrical_parameters_edit, name='asymmetrical_parameters_edit'),
    # Warning & delete parameters
    path('project/<int:project_id>/sample/<int:sample_id>/asymparameters/<int:parameter_id>/delete_warning/', views.asymmetrical_parameter_delete_warning, name='asymmetrical_parameter_delete_warning'),

    ## Data
    # Upload data
    path('project/<int:project_id>/sample/<int:sample_id>/data/upload/', views.data_upload, name='data_upload'),
    # Edit data files
    path('project/<int:project_id>/sample/<int:sample_id>/data/<int:data_id>/edit/', views.data_edit, name='data_edit'),
    # Warning & delete data
    path('project/<int:project_id>/sample/<int:sample_id>/data/<int:data_id>/delete_warning/', views.data_delete_warning, name='data_delete_warning'),

    ## Fitting
    # Sym fit page
    path('project/<int:project_id>/sample/<int:sample_id>/fit/<int:parameter_id>/', views.fit_main, name='fit_main'),
]