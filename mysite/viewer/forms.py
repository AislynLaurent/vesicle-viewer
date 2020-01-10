from django import forms

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
            'lipid_area_upperbound',
            'lipid_area_lowerbound',
            'lipid_area_lock',
            'headgroup_thickness_upperbound',
            'headgroup_thickness_lowerbound',
            'headgroup_thickness_lock',
            # 'headgroup_volume_upperbound',
            # 'headgroup_volume_lowerbound',
            # 'headgroup_volume_lock',
            'terminal_methyl_volume_upperbound',
            'terminal_methyl_volume_lowerbound',
            'terminal_methyl_volume_lock',
            # 'chain_volume_upperbound',
            # 'chain_volume_lowerbound',
            # 'chain_volume_lock',
            'sigma_upperbound',
            'sigma_lowerbound',
            'sigma_lock',
        ]

        widgets = { 

            'lipid_area_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'headgroup_thickness_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            # 'headgroup_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'terminal_methyl_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            # 'chain_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'sigma_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            

            'lipid_area_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'headgroup_thickness_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            # 'headgroup_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'terminal_methyl_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            # 'chain_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'sigma_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),

            'lipid_area_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'headgroup_thickness_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            # 'headgroup_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'terminal_methyl_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            # 'chain_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
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
