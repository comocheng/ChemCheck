'''This file tests if the forms are valid.'''

from django.test import TestCase
from upload.forms import *
from unittest import mock
from django.core.files.uploadedfile import SimpleUploadedFile

class SetupForms(TestCase):
    def Setup(self):
        fake_file = SimpleUploadedFile(
            'test-file.txt',
            'this is what\'s in the file that isn\'t a file'.encode(),
            content_type='test/plain')
        self.mech = Mechanism.objects.create(
            id=3,
            ck_mechanism_file=fake_file,
            ck_thermo_file=fake_file,
            ck_transport_file=fake_file,
            ck_surface_file=fake_file
            )
class TestForms(TestCase):
    def test_forms(self):
        form = ChemkinUpload(data={'ck_mechanism_file':"fake_file",
                                    'ck_thermo_file':"fake_file",
                                    'ck_transport_file':"fake_file",
                                    'ck_surface_file':"fake_file"})
        self.assertTrue(form.is_valid())
