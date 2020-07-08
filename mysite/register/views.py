## Imports
# From Django
# Render stuff
from django.shortcuts import get_object_or_404, render
from django.shortcuts import redirect
from django.shortcuts import render
from django.http import HttpResponse
# Email stuff
from django.core.mail import send_mail
from django.core.mail import EmailMessage
# Authentication stuff
from django.utils.encoding import force_bytes, force_text
from django.utils.http import urlsafe_base64_encode, urlsafe_base64_decode
# Other stuff
from django.contrib.sites.shortcuts import get_current_site
from django.template.loader import render_to_string

# Models
from viewer.models import User
from viewer.models import ExtendedUser

# Forms
from .forms import *

# Other
from .tokens import account_activation_token

# Register main
def register(request):
    # Form
    if request.method == 'POST':
        form = RegisterForm(request.POST)
        if form.is_valid():
            # Get info from the form for later
            in_username = form.cleaned_data['username']
            in_institution = form.cleaned_data['institution']

            # Save the user
            post = form.save(commit=False)
            post.save()

            # Get user and extended user
            user = User.objects.get(username=in_username)
            xuser = ExtendedUser.objects.get(user=user)

            # Save etxended user info
            xuser.institution = in_institution
            xuser.save()

            # Deactivate the user until authenticated
            user.is_active = False
            user.save()

            # Verification email
            current_site = get_current_site(request)
            mail_subject = 'Activate your blog account.'
            message = render_to_string('register/verification_email.html', {
                'user': user,
                'domain': current_site.domain,
                'uid':urlsafe_base64_encode(force_bytes(user.id)),
                'token':account_activation_token.make_token(user),
            })
            to_email = form.cleaned_data.get('email')
            email = EmailMessage(
                        mail_subject, message, to=[to_email]
            )
            email.send()

            return redirect('please_check_email')
    else:
        form = RegisterForm()

    return render(
        request,
        'register/register.html', {
            'form':form,
        }
    )

def activate(request, uidb64, token):
    try:
        uid = force_text(urlsafe_base64_decode(uidb64))
        user = User.objects.get(id=uid)
    except(TypeError, ValueError, OverflowError, User.DoesNotExist):
        user = None
    if user is not None and account_activation_token.check_token(user, token):
        user.is_active = True
        user.save()
        login(request, user)
        # return redirect('home')
        return redirect('thank_you')
    else:
        return redirect('invalid_link')

# Click activation link message
def please_check_email(request):
    context = {}
    return render(request, 'register/acc_please_email.html', context)
# The link was invalid message
def invalid_link(request):
    context = {}
    return render(request, 'register/acc_invalid.html', context)
# Registration complete message
def thank_you(request):
    context = {}
    return render(request, 'register/acc_thank_you.html', context)