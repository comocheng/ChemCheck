from django.db import models
from django.utils import timezone
import datetime
import yaml
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import linecache


# Create your models here.

def upload_to(instance, filename):
    return 'uploads/{id}/{fn}'.format(id=instance.pk,fn=filename)


class Mechanism(models.Model):
    """
    A chemical kinetic mechanism, from Chemkin or Cantera
    """
    ck_mechanism_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True,
                                       verbose_name="Chemkin mechanism file")
    ck_thermo_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True,
                                       verbose_name="Chemkin thermo file")
    ck_transport_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True,
                                       verbose_name="Chemkin transport file")
    ck_surface_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True,
                                       verbose_name="Chemkin surface file")
    ct_mechanism_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True,
                                        verbose_name='Cantera yaml file')
    ct_conversion_errors = models.TextField(verbose_name='Errors from the ck2yaml conversion')
    timestamps = models.DateTimeField(auto_now_add=True)

    def get_absolute_url(self):
        return 'uploads/{self.id}/'
    


    def save(self, *args, **kwargs):
        """ 
        The folder used to upload the files to, depends on the id (primary key)
        of the Mechanism object. If that is newly created and not yet saved, it 
        doesn't have an id yet. So to save the Mechanism object and all its files,
        you need to save the Mechanism object first with no files if it has no id,
        before then saving it again with the files re-attached.

        This solution is based on john.don83 answer here
        https://stackoverflow.com/questions/9968532/django-admin-file-upload-with-current-model-id

        """
        if self.id is None:
            files_saved_for_later = []
            for f in self.__class__._meta.get_fields():
                if isinstance(f, models.FileField):
                    files_saved_for_later.append((f.name, getattr(self, f.name)))
                    setattr(self, f.name, None)
            # Save the model once witout files to create an id
            super(self.__class__, self).save(*args, **kwargs)
            for name, val in files_saved_for_later:
                setattr(self, name, val)
        # Save the model, with all the files
        super(self.__class__, self).save(*args, **kwargs)
