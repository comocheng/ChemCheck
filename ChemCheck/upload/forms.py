from django import forms
from .models import Mechanism


class ChemkinUpload(forms.ModelForm):
    class Meta:
        model = Mechanism
        fields = ('ck_mechanism_file',
                  'ck_thermo_file',
                  'ck_transport_file',
                  'ck_surface_file',)

class ReactionCondition(forms.ModelForm):
    class Meta:
        model = Mechanism
        fields = ('temperature',
                  'pressure')
