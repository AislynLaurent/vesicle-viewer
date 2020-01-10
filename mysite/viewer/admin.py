from django.contrib import admin
from django.contrib.auth.admin import UserAdmin

## Overall
from .models import Lipid
from .models import Molecule
from .models import Atom
# User
from .models import User
from .models import ExtendedUser

## Project
from .models import Project
# Parameters
from .models import Symmetrical_Parameters
from .models import Project_Lipid
# Data
from .models import Data_Set
from .models import Data_Lipid
from .models import Data_Lipid_Atom

## Overall
admin.site.register(Lipid)
admin.site.register(Molecule)
admin.site.register(Atom)
# User
admin.site.register(User, UserAdmin)
admin.site.register(ExtendedUser)

## Project
admin.site.register(Project)
# Parameters
admin.site.register(Symmetrical_Parameters)
admin.site.register(Project_Lipid)
# Data
admin.site.register(Data_Set)
admin.site.register(Data_Lipid)
admin.site.register(Data_Lipid_Atom)
