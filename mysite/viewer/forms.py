from django import forms
from django.forms import formset_factory

from .models import Atom
from .models import Project
from .models import Project_Lipid
from .models import Symmetrical_Parameters
from .models import Data_Set
from .models import Data_Lipid
from .models import Data_Lipid_Atom

class Project_Form(forms.ModelForm):
    class Meta:
        model = Project
        fields = [
            'project_title',
            'model_type',
            'system_tempurature',
        ]

class Project_Lipid_Form(forms.ModelForm):
    class Meta:
        model = Project_Lipid
        fields = [
            'project_lipid_name',
            'lipid_mol_fraction'
        ]

class Parameter_Form(forms.ModelForm):
    class Meta:
        model = Symmetrical_Parameters
        fields = [
            'description',
            'bilayer_thickness',
            'lipid_area',
            'headgroup_thickness',
            'terminal_methyl_volume',
            'sigma',
        ]

class Parameter_Fit_Form(forms.ModelForm):
    class Meta:
        model = Symmetrical_Parameters

        fields = [
            # LA
            'lipid_area',
            'lipid_area_upperbound',
            'lipid_area_lowerbound',
            'lipid_area_lock',

            # HT
            'headgroup_thickness',
            'headgroup_thickness_upperbound',
            'headgroup_thickness_lowerbound',
            'headgroup_thickness_lock',

            # TMV
            'terminal_methyl_volume',
            'terminal_methyl_volume_upperbound',
            'terminal_methyl_volume_lowerbound',
            'terminal_methyl_volume_lock',

            # Sigma
            'sigma',
            'sigma_upperbound',
            'sigma_lowerbound',
            'sigma_lock',

        #     # Scale
        #     'scale',
        #     'scale_upperbound',
        #     'scale_lowerbound',
        #     'scale_lock',

        #     # BG
        #     'background',
        #     'background_upperbound',
        #     'background_lowerbound',
        #     'background_lock',
        ]

        widgets = {
            # Values
            'lipid_area' : forms.NumberInput(attrs={'class' : 'value'}),
            'headgroup_thickness' : forms.NumberInput(attrs={'class' : 'value'}),
            'terminal_methyl_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'sigma' : forms.NumberInput(attrs={'class' : 'value'}),

            # Upperbound
            'lipid_area_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'headgroup_thickness_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'terminal_methyl_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'sigma_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),

            # Lowerbound
            'lipid_area_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'headgroup_thickness_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'terminal_methyl_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'sigma_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),

            # Lock
            'lipid_area_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'headgroup_thickness_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'terminal_methyl_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'sigma_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
        }

class Data_Form(forms.ModelForm):
    class Meta:
        model = Data_Set
        fields = [
            'data_set_title',
            'd2o_mol_fraction',
            'data_type',
        ]

class Data_Upload_Form(Data_Form):

    data_file = forms.FileField()

    class Meta(Data_Form.Meta):
        fields = Data_Form.Meta.fields + ['data_file',]

class Data_Range_Form(forms.Form):
    max_value = forms.FloatField(
        widget=forms.NumberInput(attrs={'class' : 'value'})
    )
    min_value = forms.FloatField(
        widget=forms.NumberInput(attrs={'class' : 'value'})
    )

class Data_Scale_Form(forms.ModelForm):
    class Meta:
        model = Data_Set
        fields = [
            'scale',
            'scale_upperbound',
            'scale_lowerbound',
            'scale_lock',
            'background',
            'background_upperbound',
            'background_lowerbound',
            'background_lock',
        ]
        widgets = {
            # Values
            'scale' : forms.NumberInput(attrs={'class' : 'value'}),
            'background' : forms.NumberInput(attrs={'class' : 'value'}),

            # Upperbound
            'scale_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'background_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),

            # Lowerbound
            'scale_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'background_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),

            # Lock
            'scale_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'background_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
        }

class Data_Lipid_Form(forms.ModelForm):
    class Meta:
        model = Data_Lipid
        fields = [
            'data_lipid_name',
            'data_lipid_suffix',
        ]

class Data_Lipid_Atom_Form(forms.Form):
    # Choices
    LOCATION_CHOICES = (
        ('HG', 'Headgroup'),
        ('TG', 'Tailgroup')
    )

    data_lipid_atom_name = forms.ModelChoiceField(queryset=Atom.objects.all())
    data_lipid_atom_ammount = forms.IntegerField()
    atom_location = forms.ChoiceField(choices=LOCATION_CHOICES)
