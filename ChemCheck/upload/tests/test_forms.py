'''This file tests if the forms are valid.'''

from django.test import TestCase, Client
from upload.forms import ChemkinUpload, ReactionCondition
from unittest import mock
from django.core.files.uploadedfile import SimpleUploadedFile
from accounts.forms import SignupForm
class SetupForms(TestCase):
    def Setup(self):
        fake_file = SimpleUploadedFile(
            'test-file.txt',
            'this is a test'.encode(),
            content_type='test/plain')
class TestForms(TestCase):
    def test_signup(self):
        form =  SignupForm(data={'username':'testuser',
                                 'email':'testuser@gmail.com',
                                 'password1':'Ab89898980',
                                 'password2':'Ab89898980'})
        self.assertTrue(form.is_valid())
        form.save()
        Client().login(username='testuser', password='Ab89898980') 
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
