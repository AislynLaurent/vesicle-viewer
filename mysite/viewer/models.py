# Basics
from __future__ import unicode_literals
from django.db import models
from django.conf import settings
from django.db.models.signals import *
from django.core.exceptions import ValidationError
from django.contrib.contenttypes.fields import GenericForeignKey
from django.contrib.contenttypes.models import ContentType

# addons
from autoslug import AutoSlugField

# Add-ons
from django.contrib.auth.models import AbstractUser
from django.contrib.postgres.fields import ArrayField

## User
class User(AbstractUser):
    pass

class ExtendedUser(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE)

    # Institution
    institution = models.CharField(verbose_name='institution', max_length=500, default='None')

    # Tutorials?
    display_tutorial = models.BooleanField(verbose_name='display_tutorial', default=True)

def create_user_profile(sender, instance, created, **kwargs):
    if created:
        ExtendedUser.objects.create(user=instance)

post_save.connect(create_user_profile, sender=User)

## Constants
# Lipids
class Lipid(models.Model):
    # Names
    lipid_name = models.CharField(verbose_name='lipid name', max_length=100, primary_key=True, unique=True)

    # Head group
    hg_scattering = models.FloatField(verbose_name='head group scattering length', default=0)
    hg_electrons = models.FloatField(verbose_name='head group electrons', default=0)
    hg_volume = models.FloatField(verbose_name='head group volume', default=0)

    # Tail group
    tg_scattering = models.FloatField(verbose_name='tail group scattering length', default=0)
    tg_electrons = models.FloatField(verbose_name='tail group electrons', default=0)

    # Terminal methyl
    tm_scattering = models.FloatField(verbose_name='terminal methyl scattering length', default=0)
    tm_electrons = models.FloatField(verbose_name='terminal methyl electrons', default=0)

    # Total volume
    total_volume_equation = models.CharField(verbose_name='total volume equation', max_length=800, default='x')
    
    # Slug
    slug = AutoSlugField(populate_from='lipid_name', always_update=True)

    # Meta
    class Meta:
        verbose_name = 'lipid'
        verbose_name_plural = 'lipids'

    # Methods
    def __str__(self):
        return '%s' % (self.lipid_name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.lipid_name)])

# Lipid Augmentations
class Lipid_Augmentation(models.Model):
    # Lipid
    original_lipid_name = models.ForeignKey(Lipid, related_name='original_lipid_name', on_delete=models.CASCADE)
    augmentation_suffix = models.CharField(max_length=100)

    # Headgroup
    hg_scattering_net_change = models.FloatField(verbose_name='head group scattering length net change', default=0)

    # Tailgroup
    tg_scattering_net_change = models.FloatField(verbose_name='tail group scattering length net change', default=0)

    # Terminal methyl
    tmg_scattering_net_change = models.FloatField(verbose_name='terminal methyl scattering length net change', default=0)

    # Meta
    class Meta:
        verbose_name = 'lipid augmentation'
        verbose_name_plural = 'lipid augmentations'

    # Methods
    def __str__(self):
        return '-%s' % (self.augmentation_suffix)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Misc
class Molecule(models.Model):
    # Names
    compound_name = models.CharField(verbose_name='compound name', max_length=100, primary_key=True, unique=True)
    # Total volume
    total_volume_equation = models.CharField(verbose_name='total volume equation', max_length=200, default='x')
    # Scattering length
    scattering_length = models.FloatField(verbose_name='scattering length', default=0)
    # Electrons
    electrons = models.FloatField(verbose_name='electrons', default=0)

    # Slug
    slug = AutoSlugField(populate_from='compound_name', always_update=True)

    # Meta
    class Meta:
        verbose_name = 'molecule'
        verbose_name_plural = 'molecules'

    # Methods
    def __str__(self):
        return '%s' % (self.compound_name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.compound_name)])

## User Inputs
# User Lipids
class User_Lipid(models.Model):
    # User
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='user_lipid_owner', default='admin', on_delete=models.CASCADE)

    # Names
    user_lipid_name = models.CharField(verbose_name='user lipid name', max_length=100)

    # Head group
    hg_scattering = models.FloatField(verbose_name='user head group scattering length', default=0)
    hg_electrons = models.FloatField(verbose_name='user head group electrons', default=0)
    hg_volume = models.FloatField(verbose_name='user head group volume', default=0)

    # Tail group
    tg_scattering = models.FloatField(verbose_name='user tail group scattering length', default=0)
    tg_electrons = models.FloatField(verbose_name='user tail group electrons', default=0)

    # Terminal methyl
    tm_scattering = models.FloatField(verbose_name='user terminal methyl scattering length', default=0)
    tm_electrons = models.FloatField(verbose_name='user terminal methyl electrons', default=0)

    # Slug
    slug = AutoSlugField(populate_from='user_lipid_name', always_update=True)

    # Meta
    class Meta:
        verbose_name = 'user lipid'
        verbose_name_plural = 'user lipids'

    # Methods
    def __str__(self):
        return '%s' % (self.user_lipid_name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.user_lipid_name)])

# Project
class Project(models.Model):
    # User
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='project_owner', default='admin', on_delete=models.CASCADE)

    # Choices
    MODEL_CHOICES = (
        ('SM', 'Symmetrical'),
        ('AS', 'Asymmetrical'),
    )
    # Title
    project_title = models.CharField(verbose_name='project title', max_length=200)
    # Symmetry
    model_type = models.CharField(verbose_name='model type', choices=MODEL_CHOICES, max_length=3)
    # Temp
    system_tempurature = models.FloatField(verbose_name='system tempurature', default=0)
    # Advanced
    advanced_options = models.BooleanField(verbose_name="allow advanced options", default=False)

    # Meta
    class Meta:
        verbose_name = 'project'
        verbose_name_plural = 'projects'

    # # Methods
    def __str__(self):
        return '%s' % (self.project_title)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.project_title)])

# Project Lipids
class Project_Lipid(models.Model):
    # Project
    project_title = models.ForeignKey(Project, related_name='project_lipid', on_delete=models.CASCADE)

    # Lipid
    project_lipid_name = models.ForeignKey(Lipid, related_name='project_lipid_name', blank=True, null=True, on_delete=models.CASCADE)
    project_user_lipid_name = models.ForeignKey(User_Lipid, related_name='project_user_lipid_name', blank=True, null=True, on_delete=models.CASCADE)

    # Meta
    class Meta:
        verbose_name = 'project lipid'
        verbose_name_plural = 'project lipids'

    # Methods
    def clean(self):
        # At least one lipid
        if self.project_lipid_name is None and self.project_user_lipid_name is None:
            raise ValidationError('Please select a lipid.')

        if self.project_lipid_name is not None and self.project_user_lipid_name is not None:
            raise ValidationError('Please select only one lipid.')
    
    def __str__(self):
        if self.project_lipid_name is not None:
            return '%s' % (self.project_lipid_name)
        else:
            return '%s' % (self.project_user_lipid_name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Project User Lipid Volume
class Project_User_Lipid_Volume(models.Model):
    # Project
    project_title = models.ForeignKey(Project, related_name='project_user_lipid', on_delete=models.CASCADE)

    # Lipid
    project_user_lipid_name = models.ForeignKey(User_Lipid, related_name='project_user_lipid_volume_name', on_delete=models.CASCADE)

    # Volume
    user_lipid_volume = models.FloatField(verbose_name='user lipid volume', default=0)

    # Meta
    class Meta:
        verbose_name = 'user lipid volume'
        verbose_name_plural = 'user lipid volumes'

    # Methods
    def __str__(self):
        return '%s-%s' % (self.project_title, self.project_user_lipid_name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Sample
class Sample(models.Model):
    # Project
    project_title = models.ForeignKey(Project, related_name='sample', on_delete=models.CASCADE)
    #Title
    sample_title = models.CharField(verbose_name='sample title', max_length=200)

    # Meta
    class Meta:
        verbose_name = 'sample'
        verbose_name_plural = 'samples'

    # # Methods
    def __str__(self):
        return '%s-%s' % (self.project_title, self.sample_title)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.project_title)])

# Sample Lipids
class Sample_Lipid(models.Model):
    # Choices
    LOCATION_CHOICES = (
        ('IN', 'Inner'),
        ('OUT', 'Outer'),
        ('BOTH', 'Both'),
    )

    # Sample
    sample_title = models.ForeignKey(Sample, related_name='sample_lipid', on_delete=models.CASCADE)
    # Lipid
    sample_lipid_name = models.ForeignKey(Project_Lipid, related_name='sample_lipid_name', on_delete=models.CASCADE)
    
    # Percentage
    lipid_mol_fraction = models.FloatField(verbose_name='sample_lipid_mol_fraction', default=0)
    # Leaflet
    lipid_location = models.CharField(verbose_name='model', choices=LOCATION_CHOICES, max_length=4)

    # Meta
    def clean(self):
        # Mol fraction
        if self.lipid_mol_fraction < 0:
            raise ValidationError('Mol fraction cannot be less than 0')

        if self.lipid_mol_fraction > 1:
            raise ValidationError('Mol fraction cannot be greater than 1')

    class Meta:
        verbose_name = 'sample lipid'
        verbose_name_plural = 'sample lipids'

    # Methods
    def __str__(self):
        return '%s' % (self.sample_lipid_name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Custom lipid Augmentations
class Sample_Lipid_Augmentation(models.Model):
    # Lipid
    sample_lipid_name = models.ForeignKey(Sample_Lipid, related_name='augment_sample_lipid_name', on_delete=models.CASCADE)

    # suffix
    augmentation_suffix = models.CharField(max_length=100)

    # Headgroup
    hg_scattering_net_change = models.FloatField(verbose_name='head group scattering length net change', default=0)

    # Tailgroup
    tg_scattering_net_change = models.FloatField(verbose_name='tail group scattering length net change', default=0)

    # Terminal methyl
    tmg_scattering_net_change = models.FloatField(verbose_name='terminal methyl scattering length net change', default=0)

    # Meta
    class Meta:
        verbose_name = 'sample lipid augmentation'
        verbose_name_plural = 'sample lipid augmentations'

    # Methods
    def __str__(self):
        return '-%s' % (self.augmentation_suffix)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Symmetrical Parameters
class Symmetrical_Parameters(models.Model):
    ## Project
    sample_title = models.ForeignKey(Sample, related_name='sym_parameters', on_delete=models.CASCADE)

    ## Description
    name = models.CharField(verbose_name="name", max_length=255)

    ## Use the separated form factor?
    separated = models.BooleanField(verbose_name="separated form factor", default=False)

    ## Fixed
    # Chain volume
    chain_volume = models.FloatField(verbose_name='chain volume', default=0)
    chain_volume_upperbound = models.FloatField(verbose_name='cv upper bound', default=100)
    chain_volume_lowerbound = models.FloatField(verbose_name='cv lower bound', default=0)
    chain_volume_lock = models.BooleanField(verbose_name='cv lock', default=True)

    # Headgroup volume
    headgroup_volume = models.FloatField(verbose_name='headgroup volume', default=0)
    headgroup_volume_upperbound = models.FloatField(verbose_name='hv upper bound', default=100)
    headgroup_volume_lowerbound = models.FloatField(verbose_name='hv lower bound', default=0)
    headgroup_volume_lock = models.BooleanField(verbose_name='hv lock', default=True)

    ## Varied
    # Terminal methyl volume
    terminal_methyl_volume = models.FloatField(verbose_name='terminal methyl volume', default=55)
    terminal_methyl_volume_upperbound = models.FloatField(verbose_name='tmv upper bound', default=60)
    terminal_methyl_volume_lowerbound = models.FloatField(verbose_name='tmv lower bound', default=45)
    terminal_methyl_volume_lock = models.BooleanField(verbose_name='tmv lock', default=False)

    # Lipid area
    lipid_area = models.FloatField(verbose_name='lipid area', default=60)
    lipid_area_upperbound = models.FloatField(verbose_name='la upper bound', default=80)
    lipid_area_lowerbound = models.FloatField(verbose_name='la lower bound', default=40)
    lipid_area_lock = models.BooleanField(verbose_name='la lock', default=False)

    # Headgroup thickness
    headgroup_thickness = models.FloatField(verbose_name='headgroup thickness',default=8)
    headgroup_thickness_upperbound = models.FloatField(verbose_name='ht upper bound', default=15)
    headgroup_thickness_lowerbound = models.FloatField(verbose_name='ht lower bound', default=5)
    headgroup_thickness_lock = models.BooleanField(verbose_name='ht lock', default=False)

    ## Other
    # SIG
    sigma = models.FloatField(verbose_name='sigma', default=2.5)
    sigma_upperbound = models.FloatField(verbose_name='sig upper bound', default=4)
    sigma_lowerbound = models.FloatField(verbose_name='sig lower bound', default=2)
    sigma_lock = models.BooleanField(verbose_name='sig lock', default=True)

    ## Separated form factor
    # Average vesicle radius
    average_vesicle_radius = models.FloatField(verbose_name='average vesicle radius', default=500)
    average_vesicle_radius_upperbound = models.FloatField(verbose_name='avr upper bound', default=1000)
    average_vesicle_radius_lowerbound = models.FloatField(verbose_name='avr lower bound', default=250)
    average_vesicle_radius_lock = models.BooleanField(verbose_name='avr lock', default=True)

    # Relative size polydispersity
    relative_size = models.FloatField(verbose_name='relative size', default=0.25)
    relative_size_upperbound = models.FloatField(verbose_name='rs upper bound', default=0.7)
    relative_size_lowerbound = models.FloatField(verbose_name='rs lower bound', default=0.1)
    relative_size_lock = models.BooleanField(verbose_name='rs lock', default=True)

    ## Report
    fit_report = ArrayField(models.CharField(max_length=500, blank=True), verbose_name='fitted parameter report', null=True, blank=True)

    # Meta
    class Meta:
        verbose_name = 'sym parameter'
        verbose_name_plural = 'sym parameters'

    # Methods
    def clean(self):
        ## Fixed
        # Chain volume
        if self.chain_volume_upperbound == self.chain_volume_lowerbound:
            raise ValidationError('Chain volume upper and lower bound cannot be equal')

        if self.chain_volume_upperbound < self.chain_volume_lowerbound:
            raise ValidationError('Chain volume lower bound cannot be greater than it\'s upper bound')

        # Headgroup volume
        if self.headgroup_volume_upperbound == self.headgroup_volume_lowerbound:
            raise ValidationError('Headgroup volume upper and lower bound cannot be equal')

        if self.headgroup_volume_upperbound < self.headgroup_volume_lowerbound:
            raise ValidationError('Headgroup volume lower bound cannot be greater than it\'s upper bound')

        ## Varied
        # Terminal methyl volume
        if self.terminal_methyl_volume_upperbound == self.terminal_methyl_volume_lowerbound:
            raise ValidationError('Terminal methyl volume upper and lower bound cannot be equal')

        if self.terminal_methyl_volume_upperbound < self.terminal_methyl_volume_lowerbound:
            raise ValidationError('Terminal methyl volume lower bound cannot be greater than it\'s upper bound')

        # Lipid area
        if self.lipid_area <= 0:
            raise ValidationError('Lipid area cannot be less than or equal to zero')

        if self.lipid_area_upperbound == self.lipid_area_lowerbound:
            raise ValidationError('Lipid area upper and lower bound cannot be equal')

        if self.lipid_area_upperbound < self.lipid_area_lowerbound:
            raise ValidationError('Lipid area lower bound cannot be greater than it\'s upper bound')

        # Headgroup thickness
        if self.headgroup_thickness_upperbound == self.headgroup_thickness_lowerbound:
            raise ValidationError('Headgroup thickness upper and lower bound cannot be equal')

        if self.headgroup_thickness_upperbound < self.headgroup_thickness_lowerbound:
            raise ValidationError('Headgroup thickness lower bound cannot be greater than it\'s upper bound')

        ## Other
        # SIG
        if self.sigma_upperbound == self.sigma_lowerbound:
            raise ValidationError('Sigma upper and lower bound cannot be equal')

        if self.sigma_upperbound < self.sigma_lowerbound:
            raise ValidationError('Sigma lower bound cannot be greater than it\'s upper bound')

        ## Separated form factor
        # Average vesicle radius
        if self.average_vesicle_radius_upperbound == self.average_vesicle_radius_lowerbound:
            raise ValidationError('Average radius upper and lower bound cannot be equal')

        if self.average_vesicle_radius_upperbound < self.average_vesicle_radius_lowerbound:
            raise ValidationError('Average radius lower bound cannot be greater than it\'s upper bound')

        # Relative size polydispersity
        if self.relative_size_upperbound == self.relative_size_lowerbound:
            raise ValidationError('Relative size upper and lower bound cannot be equal')

        if self.relative_size_upperbound < self.relative_size_lowerbound:
            raise ValidationError('Relative size lower bound cannot be greater than it\'s upper bound')

    def __str__(self):
        return '%s-%s' % (self.sample_title, self.name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Asymmetrical Parameters
class Asymmetrical_Parameters(models.Model):
    ## Project
    sample_title = models.ForeignKey(Sample, related_name='asym_parameters', on_delete=models.CASCADE)

    ## Description
    name = models.CharField(verbose_name="name", max_length=255)

    ## Use the separated form factor?
    separated = models.BooleanField(verbose_name="separated form factor", default=False)

    ### Inside
    ## Fixed
    # Chain volume
    in_chain_volume = models.FloatField(verbose_name='in chain volume', default=0)
    in_chain_volume_upperbound = models.FloatField(verbose_name='in cv upper bound', default=100)
    in_chain_volume_lowerbound = models.FloatField(verbose_name='in cv lower bound', default=0)
    in_chain_volume_lock = models.BooleanField(verbose_name='in cv lock', default=True)

    # Headgroup volume
    in_headgroup_volume = models.FloatField(verbose_name='in headgroup volume', default=0)
    in_headgroup_volume_upperbound = models.FloatField(verbose_name='in hv upper bound', default=100)
    in_headgroup_volume_lowerbound = models.FloatField(verbose_name='in hv lower bound', default=0)
    in_headgroup_volume_lock = models.BooleanField(verbose_name='in hv lock', default=True)

    ## Varied
    # Terminal methyl volume
    in_terminal_methyl_volume = models.FloatField(verbose_name='in terminal methyl volume', default=55)
    in_terminal_methyl_volume_upperbound = models.FloatField(verbose_name='in tmv upper bound', default=60)
    in_terminal_methyl_volume_lowerbound = models.FloatField(verbose_name='in tmv lower bound', default=45)
    in_terminal_methyl_volume_lock = models.BooleanField(verbose_name='in tmv lock', default=False)

    # Lipid area
    in_lipid_area = models.FloatField(verbose_name='in lipid area', default=60)
    in_lipid_area_upperbound = models.FloatField(verbose_name='in la upper bound', default=80)
    in_lipid_area_lowerbound = models.FloatField(verbose_name='in la lower bound', default=40)
    in_lipid_area_lock = models.BooleanField(verbose_name='in la lock', default=False)

    # Headgroup thickness
    in_headgroup_thickness = models.FloatField(verbose_name='in headgroup thickness',default=8)
    in_headgroup_thickness_upperbound = models.FloatField(verbose_name='in ht upper bound', default=15)
    in_headgroup_thickness_lowerbound = models.FloatField(verbose_name='in ht lower bound', default=5)
    in_headgroup_thickness_lock = models.BooleanField(verbose_name='in ht lock', default=False)

    ### Outside
    ## Fixed
    # Chain volume
    out_chain_volume = models.FloatField(verbose_name='out chain volume', default=0)
    out_chain_volume_upperbound = models.FloatField(verbose_name='out cv upper bound', default=100)
    out_chain_volume_lowerbound = models.FloatField(verbose_name='out cv lower bound', default=0)
    out_chain_volume_lock = models.BooleanField(verbose_name='out cv lock', default=True)

    # Headgroup volume
    out_headgroup_volume = models.FloatField(verbose_name='out headgroup volume', default=0)
    out_headgroup_volume_upperbound = models.FloatField(verbose_name='out hv upper bound', default=100)
    out_headgroup_volume_lowerbound = models.FloatField(verbose_name='out hv lower bound', default=0)
    out_headgroup_volume_lock = models.BooleanField(verbose_name='out hv lock', default=True)

    ## Varied
    # Terminal methyl volume
    out_terminal_methyl_volume = models.FloatField(verbose_name='out terminal methyl volume', default=55)
    out_terminal_methyl_volume_upperbound = models.FloatField(verbose_name='out tmv upper bound', default=60)
    out_terminal_methyl_volume_lowerbound = models.FloatField(verbose_name='out tmv lower bound', default=45)
    out_terminal_methyl_volume_lock = models.BooleanField(verbose_name='out tmv lock', default=False)

    # Lipid area
    out_lipid_area = models.FloatField(verbose_name='out lipid area', default=60)
    out_lipid_area_upperbound = models.FloatField(verbose_name='out la upper bound', default=80)
    out_lipid_area_lowerbound = models.FloatField(verbose_name='out la lower bound', default=40)
    out_lipid_area_lock = models.BooleanField(verbose_name='out la lock', default=False)

    # Headgroup thickness
    out_headgroup_thickness = models.FloatField(verbose_name='out headgroup thickness',default=8)
    out_headgroup_thickness_upperbound = models.FloatField(verbose_name='out ht upper bound', default=15)
    out_headgroup_thickness_lowerbound = models.FloatField(verbose_name='out ht lower bound', default=5)
    out_headgroup_thickness_lock = models.BooleanField(verbose_name='out ht lock', default=False)

    ## Other
    # SIG
    sigma = models.FloatField(verbose_name='sigma', default=2.5)
    sigma_upperbound = models.FloatField(verbose_name='sig upper bound', default=4)
    sigma_lowerbound = models.FloatField(verbose_name='sig lower bound', default=2)
    sigma_lock = models.BooleanField(verbose_name='sig lock', default=True)

    ## Separated form factor
    # Average vesicle radius
    average_vesicle_radius = models.FloatField(verbose_name='average vesicle radius', default=500)
    average_vesicle_radius_upperbound = models.FloatField(verbose_name='avr upper bound', default=1000)
    average_vesicle_radius_lowerbound = models.FloatField(verbose_name='avr lower bound', default=250)
    average_vesicle_radius_lock = models.BooleanField(verbose_name='avr lock', default=True)

    # Relative size polydispersity
    relative_size = models.FloatField(verbose_name='relative size', default=0.25)
    relative_size_upperbound = models.FloatField(verbose_name='rs upper bound', default=0.7)
    relative_size_lowerbound = models.FloatField(verbose_name='rs lower bound', default=0.1)
    relative_size_lock = models.BooleanField(verbose_name='rs lock', default=True)

    ## Report
    fit_report = ArrayField(models.CharField(max_length=500, blank=True), verbose_name='fitted parameter report', null=True, blank=True)

    # Meta
    class Meta:
        verbose_name = 'asym parameter'
        verbose_name_plural = 'asym parameters'

    # Methods
    def clean(self):
        ### Inside
        ## Fixed
        # Chain volume
        if self.in_chain_volume_upperbound == self.in_chain_volume_lowerbound:
            raise ValidationError('Inner chain volume upper and lower bound cannot be equal')

        if self.in_chain_volume_upperbound < self.in_chain_volume_lowerbound:
            raise ValidationError('Inner chain volume lower bound cannot be greater than it\'s upper bound')

        # Headgroup volume
        if self.in_headgroup_volume_upperbound == self.in_headgroup_volume_lowerbound:
            raise ValidationError('Inner head group volume upper and lower bound cannot be equal')

        if self.in_headgroup_volume_upperbound < self.in_headgroup_volume_lowerbound:
            raise ValidationError('Inner head group volume lower bound cannot be greater than it\'s upper bound')

        ## Varied
        # Terminal methyl volume
        if self.in_terminal_methyl_volume_upperbound == self.in_terminal_methyl_volume_lowerbound:
            raise ValidationError('Inner head group volume upper and lower bound cannot be equal')

        if self.in_terminal_methyl_volume_upperbound < self.in_terminal_methyl_volume_lowerbound:
            raise ValidationError('Inner head group volume lower bound cannot be greater than it\'s upper bound')

        # Lipid area
        if self.in_lipid_area <= 0:
            raise ValidationError('Lipid area cannot be less than or equal to zero')

        if self.in_lipid_area_upperbound == self.in_lipid_area_lowerbound:
            raise ValidationError('Inner lipid area upper and lower bound cannot be equal')

        if self.in_lipid_area_upperbound < self.in_lipid_area_lowerbound:
            raise ValidationError('Inner lipid area lower bound cannot be greater than it\'s upper bound')

        # Headgroup thickness
        if self.in_headgroup_thickness_upperbound == self.in_headgroup_thickness_lowerbound:
            raise ValidationError('Inner head group thickness upper and lower bound cannot be equal')

        if self.in_headgroup_thickness_upperbound < self.in_headgroup_thickness_lowerbound:
            raise ValidationError('Inner head group thickness lower bound cannot be greater than it\'s upper bound')

        ### Outside
        ## Fixed
        # Chain volume
        if self.out_chain_volume_upperbound == self.out_chain_volume_lowerbound:
            raise ValidationError('Outer chain volume upper and lower bound cannot be equal')

        if self.out_chain_volume_upperbound < self.out_chain_volume_lowerbound:
            raise ValidationError('Outer chain volume lower bound cannot be greater than it\'s upper bound')

        # Headgroup volume
        if self.out_headgroup_volume_upperbound == self.out_headgroup_volume_lowerbound:
            raise ValidationError('Outer head group volume upper and lower bound cannot be equal')

        if self.out_headgroup_volume_upperbound < self.out_headgroup_volume_lowerbound:
            raise ValidationError('Outer head group volume lower bound cannot be greater than it\'s upper bound')

        ## Varied
        # Terminal methyl volume
        if self.out_terminal_methyl_volume_upperbound == self.out_terminal_methyl_volume_lowerbound:
            raise ValidationError('Outer terminal methyl volume upper and lower bound cannot be equal')

        if self.out_terminal_methyl_volume_upperbound < self.out_terminal_methyl_volume_lowerbound:
            raise ValidationError('Outer terminal methyl volume lower bound cannot be greater than it\'s upper bound')

        # Lipid area
        if self.out_lipid_area <= 0:
            raise ValidationError('Lipid area cannot be less than or equal to zero')

        if self.out_lipid_area_upperbound == self.out_lipid_area_lowerbound:
            raise ValidationError('Outer lipid area upper and lower bound cannot be equal')

        if self.out_lipid_area_upperbound < self.out_lipid_area_lowerbound:
            raise ValidationError('Outer lipid area lower bound cannot be greater than it\'s upper bound')

        # Headgroup thickness
        if self.out_headgroup_thickness_upperbound == self.out_headgroup_thickness_lowerbound:
            raise ValidationError('Outer head group thickness upper and lower bound cannot be equal')

        if self.out_headgroup_thickness_upperbound < self.out_headgroup_thickness_lowerbound:
            raise ValidationError('Outer head group thickness lower bound cannot be greater than it\'s upper bound')

        ## Other
        # SIG
        if self.sigma_upperbound == self.sigma_lowerbound:
            raise ValidationError('Sigma upper and lower bound cannot be equal')

        if self.sigma_upperbound < self.sigma_lowerbound:
            raise ValidationError('Sigma lower bound cannot be greater than it\'s upper bound')

        ## Separated form factor
        # Average vesicle radius
        if self.average_vesicle_radius_upperbound == self.average_vesicle_radius_lowerbound:
            raise ValidationError('Average radius upper and lower bound cannot be equal')

        if self.average_vesicle_radius_upperbound < self.average_vesicle_radius_lowerbound:
            raise ValidationError('Average radius lower bound cannot be greater than it\'s upper bound')

        # Relative size polydispersity
        if self.relative_size_upperbound == self.relative_size_lowerbound:
            raise ValidationError('Relative size upper and lower bound cannot be equal')

        if self.relative_size_upperbound < self.relative_size_lowerbound:
            raise ValidationError('Relative size lower bound cannot be greater than it\'s upper bound')


    def __str__(self):
        return '%s-%s' % (self.sample_title, self.name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Data Set
class Data_Set(models.Model):
    # Project
    sample_title = models.ForeignKey(Sample, related_name='datas', on_delete=models.CASCADE)

    # Choices
    DATA_TYPE_CHOICES = (
        ('XR', 'X-Ray'),
        ('NU', 'Neutron'),
    )

    # File information
    data_set_title = models.CharField(max_length=255, blank=True)
    upload_time = models.DateTimeField(auto_now_add=True)

    # Data information
    d2o_mol_fraction = models.FloatField(verbose_name='d2o mol fraction', default=0)
    data_type = models.CharField(verbose_name='data type', choices=DATA_TYPE_CHOICES, max_length=2, default='XR')

    # Data values
    # Values
    q_value = ArrayField(models.FloatField(), verbose_name='q_values', blank=True)
    intensity_value = ArrayField(models.FloatField(), verbose_name='intensity_values', blank=True)
    error_value = ArrayField(models.FloatField(), verbose_name='error_values', blank=True)

    # q range
    max_index = models.IntegerField(verbose_name='max index', blank=True, null=True)
    min_index = models.IntegerField(verbose_name='min index', blank=True, null=True)

    # Tweaks
    # SC
    scale = models.FloatField(verbose_name='scale', default=1)
    scale_upperbound = models.FloatField(verbose_name='scale upper bound', default=10)
    scale_lowerbound = models.FloatField(verbose_name='scale lower bound', default=1e-6)
    scale_lock = models.BooleanField(verbose_name='scale lock', default=False)

    # BG
    background = models.FloatField(verbose_name='background', default=0)
    background_upperbound = models.FloatField(verbose_name='bg upper bound', default=1)
    background_lowerbound = models.FloatField(verbose_name='bg lower bound', default=-1)
    background_lock = models.BooleanField(verbose_name='bg lock', default=False)

    # Meta
    class Meta:
        verbose_name = 'data_set'
        verbose_name_plural = 'data_sets'

    # Methods
    def clean(self):
        # D2O
        if self.d2o_mol_fraction < 0:
            raise ValidationError('Mol fraction cannot be less than 0')

        if self.d2o_mol_fraction > 1:
            raise ValidationError('Mol fraction cannot be greater than 1')

        ## Tweaks
        # SC
        if self.scale_upperbound == self.scale_lowerbound:
            raise ValidationError('Scale upper and lower bound for \"' + self.data_set_title + '\" cannot be equal')

        if self.scale_upperbound < self.scale_lowerbound:
            raise ValidationError('Scale lower bound for \"' + self.data_set_title + '\" cannot be greater than it\'s upper bound')

        # BG
        if self.background_upperbound == self.background_lowerbound:
            raise ValidationError('Background upper and lower bound for \"' + self.data_set_title + '\" cannot be equal')

        if self.background_upperbound < self.background_lowerbound:
            raise ValidationError('Background lower bound for \"' + self.data_set_title + '\" cannot be greater than it\'s upper bound')

    def __str__(self):
        return '%s' % (self.data_set_title)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Data augmentation - internal contrast
class Data_Sample_Lipid_Augment(models.Model):
    # Lipid
    sample_lipid_name = models.ForeignKey(Sample_Lipid, related_name='data_sample_lipid_name', on_delete=models.CASCADE)
    #Dataset
    data_set_title = models.ForeignKey(Data_Set, related_name='data_set_augment_title', on_delete=models.CASCADE)

    # Lipid augmentation
    sample_lipid_augment = models.ForeignKey(Lipid_Augmentation, related_name='sample_lipid_augment', blank=True, null=True, on_delete=models.CASCADE)
    sample_lipid_custom_augment = models.ForeignKey(Sample_Lipid_Augmentation, related_name='data_sample_lipid_custom_augment', blank=True, null=True, on_delete=models.CASCADE)

    # Meta
    class Meta:
        verbose_name = 'data lipid augment'
        verbose_name_plural = 'data lipid augments'

    # Methods
    def clean(self):
        # At least one augment
        if self.sample_lipid_augment is None and self.sample_lipid_custom_augment is None:
            raise ValidationError('Please select an augment.')
        # Only one augment
        if self.sample_lipid_augment is not None and self.sample_lipid_custom_augment is not None:
            raise ValidationError('Please select only one augmentation.')
    
    def __str__(self):
        return '%s-%s' % (self.sample_lipid_name, self.data_set_title)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])