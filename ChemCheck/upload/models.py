from django.db import models
from django.utils import timezone
import datetime

# Create your models here.

def upload_to(instance, filename):
    return 'uploads/{id}/{fn}'.format(id=instance.pk,fn=filename)

class Chemkin(models.Model):

    mechanism_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True)
    thermo_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True)
    transport_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True)
    surface_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True)
    timestamps = models.DateTimeField(auto_now_add=True)

    def get_absolute_url(self):
        return 'uploads/%Y%m%d-%H%M%S/chem.inp'


    def save(self, *args, **kwargs):
        """ 
        The folder used to upload the files to, depends on the id (primary key)
        of the Chemkin object. If that is newly created and not yet saved, it 
        doesn't have an id yet. So to save the Chemkin object and all its files,
        you need to save the Chemkin object first with no files if it has no id,
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

