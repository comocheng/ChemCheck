from django import forms
from .models import Chemkin

class Chemkinupload(forms.ModelForm):
    class Meta:
        model = Chemkin
        fields = ('upload_file',)

