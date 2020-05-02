'''This file tests if the forms are valid.'''

from django.test import TestCase
from upload.forms import ChemkinUpload, ReactionCondition
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
    
    def test_reaction_condition_forms(self):
        form = ReactionCondition(data={'temperature':298,
                                        'pressure':1e5})
        self.assertTrue(form.is_valid())

        form = ReactionCondition(data={'temperature':99,
                                        'pressure':1e5})
        self.assertFalse(form.is_valid())
