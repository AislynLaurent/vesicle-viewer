from django.contrib import admin
from django.contrib.auth.admin import UserAdmin

## Models
from .models import *

## Overall
admin.site.register(Lipid)
admin.site.register(Lipid_Augmentation)
admin.site.register(Molecule)
# User
admin.site.register(User, UserAdmin)
admin.site.register(ExtendedUser)

## Project
admin.site.register(Project)
admin.site.register(Sample)
# Parameters
admin.site.register(Symmetrical_Parameters)
admin.site.register(Asymmetrical_Parameters)
admin.site.register(Project_Lipid)
admin.site.register(Sample_Lipid)
admin.site.register(Sample_Lipid_Augmentation)
# Data
admin.site.register(Data_Set)
