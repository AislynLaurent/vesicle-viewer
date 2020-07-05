## Imports
# From Django
from django.shortcuts import get_object_or_404, render
from django.shortcuts import redirect
from django.shortcuts import render

# Models
from viewer.models import ExtendedUser

# Forms
from .forms import *

# Register main
def register(request):
    # Form
    if request.method == 'POST':
        form = RegisterForm(request.POST)
        if form.is_valid():
            in_username = form.cleaned_data['username']
            in_institution = form.cleaned_data['institution']

            post = form.save(commit=False)
            post.save()

            xuser = ExtendedUser.objects.get(user__username=in_username)
            xuser.institution = in_institution
            xuser.save()

            return redirect('viewer:index')
    else:
        form = RegisterForm()

    return render(
        request,
        'register/register.html', {
            'form':form,
        }
    )