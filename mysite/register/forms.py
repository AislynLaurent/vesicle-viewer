## Import
# From Django
from django import forms
from django.contrib.auth import login, authenticate
from django.contrib.auth.forms import UserCreationForm

# Models
from viewer.models import User

class RegisterForm(UserCreationForm):
    institution = forms.CharField(label='Institution', max_length=250, required=True)

    class Meta:
        model = User
        fields = [
            'username',
            'first_name',
            'last_name',
            'institution',
            'email',
            'password1',
            'password2'
        ]