## IMPORT
# Django
from django import forms
from django.forms import formset_factory
from django.core.exceptions import ValidationError
from django.core.exceptions import NON_FIELD_ERRORS
# # Other
# from upload_validator import FileTypeValidator
# Models
from .models import *
from .validators import *

class User_Lipid_Form(forms.ModelForm):
    class Meta:
        model = User_Lipid
        fields = [
            'user_lipid_name',
            # Head group
            'hg_scattering',
            'hg_electrons',
            'hg_volume',
            # Tail group
            'tg_scattering',
            'tg_electrons',
            # Terminal methyl
            'tm_scattering',
            'tm_electrons',
        ]

class Tutorial_Form(forms.ModelForm):
    class Meta:
        model = ExtendedUser
        fields = [
            'display_tutorial',
        ]
        labels = {
            'display_tutorial': 'Show all tutorials',
        }

class Project_Form(forms.ModelForm):
    class Meta:
        model = Project
        fields = [
            'project_title',
            'model_type',
            'system_tempurature',
        ]
        labels = {
            'system_tempurature': 'System tempurature (in \N{DEGREE SIGN}C)',
        }

class Advanced_Options(forms.ModelForm):
    class Meta:
        model = Project
        fields = [
            'advanced_options',
        ]
        labels = {
            'advanced_options': 'Allow me to vary the headgroup and chain volume',
        }

class Project_Lipid_Form(forms.ModelForm):
    class Meta:
        model = Project_Lipid
        fields = [
            'project_lipid_name',
            'project_user_lipid_name',
        ]
        labels = {
            'project_lipid_name': 'Add a standard lipid',
            'project_user_lipid_name': 'Or, add one of your custom lipids'
        }

    def __init__(self, owner, *args, **kwargs):
        super(Project_Lipid_Form, self).__init__(*args, **kwargs)
        self.fields['project_lipid_name'].queryset = self.fields['project_lipid_name'].queryset.order_by('lipid_name')
        self.fields['project_user_lipid_name'].queryset = User_Lipid.objects.filter(owner=owner).order_by('user_lipid_name')

class Project_User_Lipid_Volume_Form(forms.ModelForm):
    class Meta:
        model = Project_User_Lipid_Volume
        fields = [
            'user_lipid_volume',
        ]

class Sample_Form(forms.ModelForm):
    project_id = 0

    def __init__(self, project_id, *args, **kwargs):
        super(Sample_Form, self).__init__(*args, **kwargs)
        self.project_id = project_id
        self.instance = kwargs["instance"]

    def clean(self):
        samples = Sample.objects.filter(project_title_id=self.project_id, sample_title=self.data['sample_title'])
        # Unique name per sample (except if instance is passed, which indicates an edit)
        if samples and self.data['sample_title'] != self.instance.sample_title:
            raise ValidationError('A sample with that name already exists for this project - please choose another')

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
        self.fields['sample_lipid_name'].queryset = Project_Lipid.objects.filter(project_title__id=project_id).order_by('project_lipid_name', 'project_user_lipid_name')

class Asym_Sample_Lipid(forms.ModelForm):
    class Meta:
        model = Sample_Lipid
        fields = [
            'sample_lipid_name',
            'lipid_mol_fraction',
        ]

# Workaround - limit choices to hide 'both' option
class Asym_Sample_Lipid_Form(Asym_Sample_Lipid):
    CHOICES = (
        ('IN','Inner'),
        ('OUT','Outer'),
    )

    location = forms.ChoiceField(choices=CHOICES)

    class Meta(Asym_Sample_Lipid.Meta):
        fields = Asym_Sample_Lipid.Meta.fields + ['location',]

    def __init__(self, project_id, *args, **kwargs):
        super(Asym_Sample_Lipid_Form, self).__init__(*args, **kwargs)

        self.fields['sample_lipid_name'].queryset = Project_Lipid.objects.filter(project_title__id=project_id)

class Lipid_Augmentation_Form(forms.ModelForm):
    class Meta:
        model = Data_Sample_Lipid_Augment
        fields = [
            'sample_lipid_augment',
            'sample_lipid_custom_augment',
            'data_set_title',
        ]

    def __init__(self, lipid_name, sample_id, *args, **kwargs):
        super(Lipid_Augmentation_Form, self).__init__(*args, **kwargs)
        self.sample_id = sample_id
        self.lipid_name = lipid_name
        self.fields['sample_lipid_augment'].queryset = Lipid_Augmentation.objects.filter(original_lipid_name__lipid_name=lipid_name)
        self.fields['sample_lipid_custom_augment'].queryset = Sample_Lipid_Augmentation.objects.filter(sample_lipid_name__sample_lipid_name=lipid_name)
        self.fields['data_set_title'].queryset = Data_Set.objects.filter(sample_title_id=sample_id)

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
            'use_structure_factor',
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

            # HV
            'headgroup_volume',
            'headgroup_volume_upperbound',
            'headgroup_volume_lowerbound',
            'headgroup_volume_lock',

            # TMV
            'terminal_methyl_volume',
            'terminal_methyl_volume_upperbound',
            'terminal_methyl_volume_lowerbound',
            'terminal_methyl_volume_lock',

            # CV
            'chain_volume',
            'chain_volume_upperbound',
            'chain_volume_lowerbound',
            'chain_volume_lock',

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
            'headgroup_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'terminal_methyl_volume' : forms.NumberInput(attrs={'class' : 'value'}),            'chain_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'sigma' : forms.NumberInput(attrs={'class' : 'value'}),
            'average_vesicle_radius' : forms.NumberInput(attrs={'class' : 'value'}),
            'relative_size' : forms.NumberInput(attrs={'class' : 'value'}),

            # Upperbound
            'lipid_area_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'headgroup_thickness_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'headgroup_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'terminal_methyl_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'chain_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'sigma_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'average_vesicle_radius_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'relative_size_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),

            # Lowerbound
            'lipid_area_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'headgroup_thickness_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'headgroup_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'terminal_methyl_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'chain_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'sigma_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'average_vesicle_radius_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'relative_size_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),

            # Lock
            'lipid_area_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'headgroup_thickness_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'headgroup_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'terminal_methyl_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'chain_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
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

            # HV
            'in_headgroup_volume',
            'in_headgroup_volume_upperbound',
            'in_headgroup_volume_lowerbound',
            'in_headgroup_volume_lock',

            # TMV
            'in_terminal_methyl_volume',
            'in_terminal_methyl_volume_upperbound',
            'in_terminal_methyl_volume_lowerbound',
            'in_terminal_methyl_volume_lock',

            # CV
            'in_chain_volume',
            'in_chain_volume_upperbound',
            'in_chain_volume_lowerbound',
            'in_chain_volume_lock',

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

            # HV
            'out_headgroup_volume',
            'out_headgroup_volume_upperbound',
            'out_headgroup_volume_lowerbound',
            'out_headgroup_volume_lock',

            # TMV
            'out_terminal_methyl_volume',
            'out_terminal_methyl_volume_upperbound',
            'out_terminal_methyl_volume_lowerbound',
            'out_terminal_methyl_volume_lock',

            # CV
            'out_chain_volume',
            'out_chain_volume_upperbound',
            'out_chain_volume_lowerbound',
            'out_chain_volume_lock',

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
            'in_headgroup_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'in_terminal_methyl_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'in_chain_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'out_lipid_area' : forms.NumberInput(attrs={'class' : 'value'}),
            'out_headgroup_thickness' : forms.NumberInput(attrs={'class' : 'value'}),
            'out_headgroup_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'out_terminal_methyl_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'out_chain_volume' : forms.NumberInput(attrs={'class' : 'value'}),
            'sigma' : forms.NumberInput(attrs={'class' : 'value'}),
            'average_vesicle_radius' : forms.NumberInput(attrs={'class' : 'value'}),
            'relative_size' : forms.NumberInput(attrs={'class' : 'value'}),

            # Upperbound
            'in_lipid_area_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'in_headgroup_thickness_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'in_headgroup_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'in_terminal_methyl_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'in_chain_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'out_lipid_area_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'out_headgroup_thickness_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'out_headgroup_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'out_terminal_methyl_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'out_chain_volume_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'sigma_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'average_vesicle_radius_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),
            'relative_size_upperbound' : forms.NumberInput(attrs={'class' : 'upper_bound'}),

            # Lowerbound
            'in_lipid_area_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'in_headgroup_thickness_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'in_headgroup_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'in_terminal_methyl_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'in_chain_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'out_lipid_area_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'out_headgroup_thickness_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'out_headgroup_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'out_terminal_methyl_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'out_chain_volume_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'sigma_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'average_vesicle_radius_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),
            'relative_size_lowerbound' : forms.NumberInput(attrs={'class' : 'lower_bound'}),

            # Lock
            'lipid_area_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'headgroup_thickness_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'headgroup_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'terminal_methyl_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'chain_volume_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'sigma_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'average_vesicle_radius_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
            'relative_size_lock' : forms.CheckboxInput(attrs={'class' : 'lock'}),
        }

class Data_Form(forms.ModelForm):
    sample_id = 0

    def __init__(self, sample_id, *args, **kwargs):
        super(Data_Form, self).__init__(*args, **kwargs)
        self.sample_id = sample_id

    def clean(self):
        datas = Data_Set.objects.filter(sample_title_id=self.sample_id, data_set_title=self.data['data_set_title'])
        # Unique name per sample
        if datas:
            raise ValidationError('A dataset with that name already exists for this sample - please choose another')

    class Meta:
        model = Data_Set
        fields = [
            'data_set_title',
            'd2o_mol_fraction',
            'data_type',
        ]

class Data_Edit_Form(forms.ModelForm):
    class Meta:
        model = Data_Set
        fields = [
            'data_set_title',
            'd2o_mol_fraction',
            'data_type',
        ]

class Data_Upload_Form(Data_Form):

    validate_file = FileValidator(max_size=1024 * 100, content_types=('text/plain',))
    data_file = forms.FileField(validators=[validate_file])

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
        if self.data['max_value'] !='' and self.data['min_value'] == '':
            raise ValidationError('Cannot set a maximum value without a minimum')

        if self.data['max_value'] =='' and self.data['min_value'] != '':
            raise ValidationError('Cannot set a minimum value without a maximum')

        if self.data['max_value'] and self.data['min_value'] != '' and self.data['max_value'] == self.data['min_value']:
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