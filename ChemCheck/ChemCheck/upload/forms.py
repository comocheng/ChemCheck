from django import forms
from .models import Chemkin


class Chemkinupload(forms.ModelForm):
    class Meta:
        model = Chemkin
        fields = ('mechanism_file',
                  'thermo_file',
                  'transport_file',
                  'surface_file')


