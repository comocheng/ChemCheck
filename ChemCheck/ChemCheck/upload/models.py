from django.db import models
from django.utils import timezone
import datetime

# Create your models here.

class Chemkin(models.Model):
    upload_file = models.FileField(upload_to='uploads/', max_length=100, blank=True)
    timestamps = models.DateTimeField(auto_now_add=True)