from django.db import models
from django.utils import timezone
import datetime

# Create your models here.

def chemkin_upload_path(instance, filename):
    
    return 'uploads/%Y%m%d-%H%M%S/chem.inp'.format(instance)

class Chemkin(models.Model):

    mechanism_file = models.FileField(upload_to='uploads/%Y%m%d-%H%M%S/chem.inp', max_length=100, blank=True)
    thermo_file = models.FileField(upload_to='uploads/%Y%m%d-%H%M%S/therm.dat', max_length=100, blank=True)
    transport_file = models.FileField(upload_to='uploads/%Y%m%d-%H%M%S/tran.dat', max_length=100, blank=True)
    surface_file = models.FileField(upload_to='uploads/%Y%m%d-%H%M%S/surf.inp', max_length=100, blank=True)
    timestamps = models.DateTimeField(auto_now_add=True)

