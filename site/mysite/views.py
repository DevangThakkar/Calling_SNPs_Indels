from django.contrib.auth.forms import UserCreationForm
from django.shortcuts import render

from .forms import ContactForm

def _form_view(request, template_name='basic.html', form_class=ContactForm):
    if request.method == 'POST':
        form = form_class(request.POST)
        if form.is_valid():
            pass  # does nothing, just trigger the validation
    else:
        form = form_class()
    return render(request, template_name, {'form': form})

def basic(request):
    return _form_view(request)

def submit(request):
    info=request.POST['info']
    print('Hello')
    print(info)