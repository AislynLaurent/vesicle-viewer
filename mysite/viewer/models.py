# Basics
from django.db import models
from django.conf import settings

# Add-ons
from django.contrib.auth.models import AbstractUser
from django.contrib.postgres.fields import ArrayField

## User
class User(AbstractUser):
    pass

class ExtendedUser(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE)

    # Tutorials?
    display_tutorial = models.BooleanField(verbose_name='display_tutorial', default=True)

## Constants
# Lipids
class Lipid(models.Model):
    # Names
    lipid_name = models.CharField(verbose_name='lipid name', max_length=100, primary_key=True, unique=True)

    # Volume
    total_volume = models.FloatField(verbose_name='total volume', default=0)

    # Head group
    hg_scattering = models.FloatField(verbose_name='head group scattering length', default=0)
    hg_electrons = models.FloatField(verbose_name='head group electrons', default=0)
    hg_volume_equation = models.CharField(verbose_name='head group volume equation', max_length=200, default='x')

    # Tail group
    tg_scattering = models.FloatField(verbose_name='tail group scattering length', default=0)
    tg_electrons = models.FloatField(verbose_name='tail group electrons', default=0)
    tg_volume_equation = models.CharField(verbose_name='tail group volume equation', max_length=200, default='x')

    # Terminal methyl
    tm_scattering = models.FloatField(verbose_name='terminal methyl scattering length', default=0)
    tm_electrons = models.FloatField(verbose_name='terminal methyl electrons', default=0)
    tg_volume_equation = models.CharField(verbose_name='tail group volume equation', max_length=200, default='x')

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

# Misc
class Molecule(models.Model):
    # Names
    compound_name = models.CharField(max_length=100, primary_key=True, unique=True)
    # Volume
    total_volume = models.FloatField(verbose_name='total volume', default=0)
    # Scattering length
    scattering_length = models.FloatField(verbose_name='scattering length', default=0)
    # Electrons
    electrons = models.FloatField(verbose_name='electrons', default=0)

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

# Atoms
class Atom(models.Model):
    # Names
    atom_name = models.CharField(max_length=100, primary_key=True, unique=True)
    # Scattering length
    scattering_length_adj = models.FloatField(verbose_name='scattering length adjustment', default=0)

    # Meta
    class Meta:
        verbose_name = 'atom'
        verbose_name_plural = 'atoms'   

    # Methods
    def __str__(self):
        return '%s' % (self.atom_name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.atom_name)])    

## User Inputs
# Project
class Project(models.Model):
    # User
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, related_name='project_owner', default='admin', on_delete=models.CASCADE)

    # Choices
    MODEL_CHOICES = (
        ('SM', 'Symmetrical'),
    )
    # Title
    project_title = models.CharField(verbose_name='project title', max_length=200)
    # Symmetry
    model_type = models.CharField(verbose_name='model', choices=MODEL_CHOICES, max_length=3)
    # Temp
    system_tempurature = models.FloatField(verbose_name='system tempurature', default=0)

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

# Parameters
class Symmetrical_Parameters(models.Model):
    ## Project
    project_title = models.ForeignKey(Project, related_name='parameters', on_delete=models.CASCADE)

    ## Description
    description = models.CharField(max_length=255, blank=True)

    ## Overall
    # Thickness
    bilayer_thickness = models.FloatField(verbose_name='bilayer thickness', default=0)
    bilayer_thickness_upperbound = models.FloatField(verbose_name='bilayer upper bound', default=0)
    bilayer_thickness_lowerbound = models.FloatField(verbose_name='bilayer lower bound', default=0)
    bilayer_thickness_lock = models.BooleanField(verbose_name='bilayer lock', default=False)

    #Radius
    vesicle_radius = models.FloatField(verbose_name='vesicle radius', default=0)
        # Calculated from the lipid volume along with the thickness

    ## Fixed
    # Chain volume
    chain_volume = models.FloatField(verbose_name='chain volume', default=0)
    chain_volume_upperbound = models.FloatField(verbose_name='cv upper bound', default=1)
    chain_volume_lowerbound = models.FloatField(verbose_name='cv lower bound', default=-1)
    chain_volume_lock = models.BooleanField(verbose_name='cv lock', default=True)

    # Headgroup volume
    headgroup_volume = models.FloatField(verbose_name='headgroup volume', default=0)
    headgroup_volume_upperbound = models.FloatField(verbose_name='hv upper bound', default=1)
    headgroup_volume_lowerbound = models.FloatField(verbose_name='hv lower bound', default=-1)
    headgroup_volume_lock = models.BooleanField(verbose_name='hv lock', default=True)

    ## Varied
    # Terminal methyl volume
    terminal_methyl_volume = models.FloatField(verbose_name='terminal methyl volume', default=0)
    terminal_methyl_volume_upperbound = models.FloatField(verbose_name='tmv upper bound', default=1)
    terminal_methyl_volume_lowerbound = models.FloatField(verbose_name='tmv lower bound', default=-1)
    terminal_methyl_volume_lock = models.BooleanField(verbose_name='tmv lock', default=False)

    # Lipid area
    lipid_area = models.FloatField(verbose_name='lipid area', default=60)
    lipid_area_upperbound = models.FloatField(verbose_name='la upper bound', default=65)
    lipid_area_lowerbound = models.FloatField(verbose_name='la lower bound', default=-65)
    lipid_area_lock = models.BooleanField(verbose_name='la lock', default=False)

    # Headgroup thickness
    headgroup_thickness = models.FloatField(verbose_name='headgroup thickness', default=6)
    headgroup_thickness_upperbound = models.FloatField(verbose_name='ht upper bound', default=7)
    headgroup_thickness_lowerbound = models.FloatField(verbose_name='ht lower bound', default=-7)
    headgroup_thickness_lock = models.BooleanField(verbose_name='ht lock', default=False)

    ## Other
    # SIG
    sigma = models.FloatField(verbose_name='sigma', default=2.5)
    sigma_upperbound = models.FloatField(verbose_name='sig upper bound', default=3)
    sigma_lowerbound = models.FloatField(verbose_name='sig lower bound', default=-3)
    sigma_lock = models.BooleanField(verbose_name='sig lock', default=True)

    ## Report
    fit_report = ArrayField(models.CharField(max_length=500, blank=True), verbose_name='fitted parameter report', null=True, blank=True)

    # Meta
    class Meta:
        verbose_name = 'parameter'
        verbose_name_plural = 'parameters'

    # Methods
    def __str__(self):
        return '%s-%s' % (self.project_title, self.description)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Project Lipids
class Project_Lipid(models.Model):
    # Project
    project_title = models.ForeignKey(Project, related_name='project_lipid', on_delete=models.CASCADE)

    # Lipid
    project_lipid_name = models.ForeignKey(Lipid, related_name='project_lipid_name', on_delete=models.CASCADE)

    # Percentage
    lipid_mol_fraction = models.FloatField(verbose_name='project_lipid_mol_fraction', default=0)

    # Meta
    class Meta:
        verbose_name = 'project lipid'
        verbose_name_plural = 'project lipids'

    # Methods
    def __str__(self):
        return '%s' % (self.project_lipid_name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Data Set
class Data_Set(models.Model):
    # Project
    project_title = models.ForeignKey(Project, related_name='datas', on_delete=models.CASCADE)

    # Choices
    DATA_TYPE_CHOICES = (
        ('XR', 'X-Ray'),
        ('NU', 'Neutron'),
    )
    
    # File information
    data_set_title = models.CharField(max_length=255, blank=True)
    upload_time = models.DateTimeField(auto_now_add=True)

    # Data information
    d2o_mol_fraction = models.FloatField(verbose_name='d2o percentage', default=0)
    data_type = models.CharField(verbose_name='data type', choices=DATA_TYPE_CHOICES, max_length=2, default='XR')

    # Data values
    # Values
    q_value = ArrayField(models.FloatField(), verbose_name='q_values', blank=True)
    intensity_value = ArrayField(models.FloatField(), verbose_name='intensity_values', blank=True)
    error_value = ArrayField(models.FloatField(), verbose_name='error_values', blank=True)

    # Tweaks
    # SC
    scale = models.FloatField(verbose_name='scale', default=1)
    scale_upperbound = models.FloatField(verbose_name='sig upper bound', default=0)
    scale_lowerbound = models.FloatField(verbose_name='sig lower bound', default=0)
    scale_lock = models.BooleanField(verbose_name='sig lock', default=True)

    # BG
    background = models.FloatField(verbose_name='background', default=0)
    background_upperbound = models.FloatField(verbose_name='sig upper bound', default=0)
    background_lowerbound = models.FloatField(verbose_name='sig lower bound', default=0)
    background_lock = models.BooleanField(verbose_name='sig lock', default=True)

    # Meta
    class Meta:
        verbose_name = 'data_set'
        verbose_name_plural = 'data_sets'

    # Methods
    def __str__(self):
        return '%s' % (self.data_set_title)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Data Lipids
class Data_Lipid(models.Model):
    # Project
    data_set_title = models.ForeignKey(Data_Set, related_name='data_lipid', on_delete=models.CASCADE)

    # Lipid
    data_lipid_name = models.ForeignKey(Project_Lipid, related_name='data_lipid_name', on_delete=models.CASCADE)
    data_lipid_suffix = models.CharField(max_length=100, blank=True)

    # Meta
    class Meta:
        verbose_name = 'data lipid'
        verbose_name_plural = 'data lipids'

    # Methods
    def __str__(self):
        return '%s-%s' % (self.data_lipid_name, self.data_lipid_suffix)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])

# Data Lipids
class Data_Lipid_Atom(models.Model):
    # Choices
    LOCATION_CHOICES = (
        ('HG', 'Headgroup'),
        ('TG', 'Tailgroup')
    )

    # Lipid
    data_lipid_name = models.ForeignKey(Data_Lipid, related_name='data_atom_lipid_name', on_delete=models.CASCADE)

    # Atom info
    data_lipid_atom_name = models.ForeignKey(Atom, related_name='data_lipid_atom_name', on_delete=models.CASCADE)
    data_lipid_atom_ammount = models.IntegerField(verbose_name='data_lipid_atom_ammount', default=0)

    # Location
    atom_location = models.CharField(verbose_name='atom_location', choices=LOCATION_CHOICES, max_length=3)

    # Meta
    class Meta:
        verbose_name = 'data lipid atom'
        verbose_name_plural = 'data lipid atoms'

    # Methods
    def __str__(self):
        return '%s-%s' % (self.data_lipid_name, self.data_lipid_atom_name)

    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('viewer.views.details', args=[str(self.id)])