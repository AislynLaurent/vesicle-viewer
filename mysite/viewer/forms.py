from django import forms
from django.forms import formset_factory
from django.core.exceptions import ValidationError

from .models import *

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
        ]

class Sample_Form(forms.ModelForm):
    class Meta:
        model = Sample
        fields = [
            'sample_title',
        ]

class Sym_Sample_Lipid_Form(forms.ModelForm):
    class Meta:
        model = Sample_Lipid
        fields = [
            'sample_lipid_name',
            'lipid_mol_fraction',
        ]

    def __init__(self, project_id, *args, **kwargs):
        super(Sym_Sample_Lipid_Form, self).__init__(*args, **kwargs)
        self.fields['sample_lipid_name'].queryset = Project_Lipid.objects.filter(project_title__id=project_id)

class Asym_Sample_Lipid_Form(forms.ModelForm):

    class Meta:
        model = Sample_Lipid
        fields = [
            'sample_lipid_name',
            'lipid_mol_fraction',
            'lipid_location',
        ]

    def __init__(self, project_id, *args, **kwargs):
        super(Asym_Sample_Lipid_Form, self).__init__(*args, **kwargs)
        self.fields['sample_lipid_name'].queryset = Project_Lipid.objects.filter(project_title__id=project_id)

class Lipid_Augmentation_Form(forms.ModelForm):
    class Meta:
        model = Sample_Lipid
        fields = [
            'sample_lipid_augment',
        ]

    def __init__(self, lipid_name, *args, **kwargs):
        super(Lipid_Augmentation_Form, self).__init__(*args, **kwargs)
        self.fields['sample_lipid_augment'].queryset = Lipid_Augmentation.objects.filter(original_lipid_name__lipid_name=lipid_name)

class Custom_Lipid_Augmentation_Form(forms.ModelForm):
    class Meta:
        model = Sample_Lipid_Augmentation
        fields = [
            'augmentation_suffix',
            'hg_scattering_net_change',
            'tg_scattering_net_change',
            'tmg_scattering_net_change',
        ]

class Symmetrical_Parameter_Form(forms.ModelForm):
    class Meta:
        model = Symmetrical_Parameters
        fields = [
            'name',
            'lipid_area',
            'headgroup_thickness',
            'terminal_methyl_volume',
            'sigma',
        ]

class Symmetrical_Parameter_Fit_Form(forms.ModelForm):
    class Meta:
        model = Symmetrical_Parameters

        fields = [
            # SFF
            'separated',
            
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

            # Average vesicle radius
            'average_vesicle_radius',
            'average_vesicle_radius_upperbound',
            'average_vesicle_radius_lowerbound',
            'average_vesicle_radius_lock',

            # Relative size polydispersity
            'relative_size',
            'relative_size_upperbound',
            'relative_size_lowerbound',
            'relative_size_lock',

        ]

        widgets = {
            # Values
            'lipid_area' : forms.NumberInput(attrs={'class' : 'value'}),
            'headgroup_thickness' : forms.NumberInput(attrs={'class' : 'value'}),
            'terminal_methyl_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'sigma' : forms.NumberInput(attrs={'class' : 'value'}),
            'average_vesicle_radius' : forms.NumberInput(attrs={'class' : 'value'}),
            'relative_size' : forms.NumberInput(attrs={'class' : 'value'}),

            # Upperbound
            'lipid_area_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'headgroup_thickness_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'terminal_methyl_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'sigma_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'average_vesicle_radius_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'relative_size_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),

            # Lowerbound
            'lipid_area_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'headgroup_thickness_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'terminal_methyl_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'sigma_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'average_vesicle_radius_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'relative_size_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),

            # Lock
            'lipid_area_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'headgroup_thickness_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'terminal_methyl_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'sigma_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'average_vesicle_radius_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'relative_size_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
        }

class Asymmetrical_Parameter_Form(forms.ModelForm):
    class Meta:
        model = Asymmetrical_Parameters
        fields = [
            'name',
            'in_lipid_area',
            'in_headgroup_thickness',
            'in_terminal_methyl_volume',
            'out_lipid_area',
            'out_headgroup_thickness',
            'out_terminal_methyl_volume',
            'sigma',
        ]

class Asymmetrical_Parameter_Fit_Form(forms.ModelForm):
    class Meta:
        model = Asymmetrical_Parameters

        fields = [
            # SFF
            'separated',
            
            # LA
            'in_lipid_area',
            'in_lipid_area_upperbound',
            'in_lipid_area_lowerbound',
            'in_lipid_area_lock',

            # HT
            'in_headgroup_thickness',
            'in_headgroup_thickness_upperbound',
            'in_headgroup_thickness_lowerbound',
            'in_headgroup_thickness_lock',

            # TMV
            'in_terminal_methyl_volume',
            'in_terminal_methyl_volume_upperbound',
            'in_terminal_methyl_volume_lowerbound',
            'in_terminal_methyl_volume_lock',

            # LA
            'out_lipid_area',
            'out_lipid_area_upperbound',
            'out_lipid_area_lowerbound',
            'out_lipid_area_lock',

            # HT
            'out_headgroup_thickness',
            'out_headgroup_thickness_upperbound',
            'out_headgroup_thickness_lowerbound',
            'out_headgroup_thickness_lock',

            # TMV
            'out_terminal_methyl_volume',
            'out_terminal_methyl_volume_upperbound',
            'out_terminal_methyl_volume_lowerbound',
            'out_terminal_methyl_volume_lock',

            # Sigma
            'sigma',
            'sigma_upperbound',
            'sigma_lowerbound',
            'sigma_lock',

            # Average vesicle radius
            'average_vesicle_radius',
            'average_vesicle_radius_upperbound',
            'average_vesicle_radius_lowerbound',
            'average_vesicle_radius_lock',

            # Relative size polydispersity
            'relative_size',
            'relative_size_upperbound',
            'relative_size_lowerbound',
            'relative_size_lock',

        ]

        widgets = {
            # Values
            'in_lipid_area' : forms.NumberInput(attrs={'class' : 'value'}),
            'in_headgroup_thickness' : forms.NumberInput(attrs={'class' : 'value'}),
            'in_terminal_methyl_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'out_lipid_area' : forms.NumberInput(attrs={'class' : 'value'}),
            'out_headgroup_thickness' : forms.NumberInput(attrs={'class' : 'value'}),
            'out_terminal_methyl_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'sigma' : forms.NumberInput(attrs={'class' : 'value'}),
            'average_vesicle_radius' : forms.NumberInput(attrs={'class' : 'value'}),
            'relative_size' : forms.NumberInput(attrs={'class' : 'value'}),

            # Upperbound
            'in_lipid_area_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'in_headgroup_thickness_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'in_terminal_methyl_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'out_lipid_area_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'out_headgroup_thickness_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'out_terminal_methyl_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'sigma_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'average_vesicle_radius_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'relative_size_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),

            # Lowerbound
            'in_lipid_area_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'in_headgroup_thickness_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'in_terminal_methyl_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'out_lipid_area_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'out_headgroup_thickness_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'out_terminal_methyl_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'sigma_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'average_vesicle_radius_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'relative_size_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),

            # Lock
            'lipid_area_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'headgroup_thickness_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'terminal_methyl_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'sigma_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'average_vesicle_radius_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'relative_size_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
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
        required=False,
        widget=forms.NumberInput(attrs={'class' : 'upper_bound'})
    )
    min_value = forms.FloatField(
        required=False,
        widget=forms.NumberInput(attrs={'class' : 'lower_bound'})
    )

    def clean(self):
        # q range
        if self.data['max_value'] == self.data['min_value']:
            raise ValidationError('Maximum and minimum values for a dataset cannot be equal')

        if self.data['max_value'] < self.data['min_value']:
            raise ValidationError('Minimum value for a dataset cannot be greater than it\'s maximum')

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