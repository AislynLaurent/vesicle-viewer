from django.urls import path
from . import views

app_name = 'viewer'

urlpatterns = [
    path('', views.index, name='index'),
    path('about', views.about, name='about'),
    path('help', views.get_help, name='help'),
    path('privacy', views.privacy, name='privacy'),
    
    ## Pojects
    # New project
    path('project/new', views.project_new, name='project_new'),
    # List projects
    path('project/list', views.project_list, name='project_list'),
    # Project details
    path('project/<int:project_id>', views.project_detail, name='project_detail'),
    # Edit a project
    path('project/<int:project_id>/edit/', views.project_edit, name='project_edit'),

    ## Project Lipids
    # Add project lipid
    path('project/<int:project_id>/lipid/add/', views.project_lipid_new, name='project_lipid_new'),
    # Edit a project lipid
    path('project/<int:project_id>/lipid/<int:lipid_id>/edit/', views.project_lipid_edit, name='project_lipid_edit'),
    # Warning & delete lipid
    path('project/<int:project_id>/lipid/<int:lipid_id>/delete_warning/', views.project_lipid_delete_warning, name='project_lipid_delete_warning'),

    ## Parameters
    # Add project parameters
    path('project/<int:project_id>/parameters/new/', views.parameters_new, name='parameters_new'),
    # Edit a projects parameters
    path('project/<int:project_id>/parameters/<int:parameter_id>/edit/', views.parameters_edit, name='parameters_edit'),
    # Warning & delete parameters
    path('project/<int:project_id>/parameters/<int:parameter_id>/delete_warning/', views.parameter_delete_warning, name='parameter_delete_warning'),

    ## Data
    # Upload data
    path('project/<int:project_id>/data/upload/', views.data_upload, name='data_upload'),
    # Edit data files
    path('project/<int:project_id>/data/<int:data_id>/edit/', views.data_edit, name='data_edit'),
    # Warning & delete data
    path('project/<int:project_id>/data/<int:data_id>/delete_warning/', views.data_delete_warning, name='data_delete_warning'),

    ## Data Lipids
    # Add lipid adjustments
    path('project/<int:project_id>/data/<int:data_id>/lipid/adjust/', views.data_lipid_new, name='data_lipid_new'),
    # Edit lipid adjustments
    path('project/<int:project_id>/data/<int:data_id>/lipid/<int:lipid_id>/edit/', views.data_lipid_edit, name='data_lipid_edit'),
    # Warning & delete lipid adjustments
    path('project/<int:project_id>/data/<int:data_id>/lipid/<int:lipid_id>/delete_warning/', views.data_lipid_delete_warning, name='data_lipid_delete_warning'),

    ## Fitting
    # Main fit page
    path('project/<int:project_id>/fit/<int:parameter_id>/', views.fit_main, name='fit_main'),
]