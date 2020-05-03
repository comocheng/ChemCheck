''' This is a test for fields in models and the upload directory'''

from django.test import TestCase
from upload.models import Mechanism, upload_to
from unittest import mock
from django.core.files import File
from django.core.files.uploadedfile import SimpleUploadedFile
from canteradebugger.settings import MEDIA_ROOT
import os

class MechanismTest(TestCase):
    def test_file_field(self):
        file_test = mock.Mock(spec=File)
        file_test.name = 'test.txt'
        file_get_tested = Mechanism(ck_mechanism_file=file_test,
                                    ck_thermo_file=file_test,
                                    ck_transport_file=file_test,
                                    ck_surface_file=file_test,
                                    ct_mechanism_file=file_test,
                                    temperature = 298,
                                    pressure = 1e5)
        self.assertEqual(file_test.name, file_get_tested.ck_mechanism_file)
        self.assertEqual(file_test.name, file_get_tested.ck_thermo_file)
        self.assertEqual(file_test.name, file_get_tested.ck_transport_file)
        self.assertEqual(file_test.name, file_get_tested.ck_surface_file)
        self.assertEqual(file_test.name, file_get_tested.ct_mechanism_file)
        self.assertEqual(298, file_get_tested.temperature)
        self.assertEqual(1e5, file_get_tested.pressure)

class MechDictTest(TestCase):
    def test_ffield(self):
        fake_file = SimpleUploadedFile(
            'test-file.txt',
            'this is what\'s in the file that isn\'t a file'.encode(),
            content_type='test/plain')
        mech = Mechanism.objects.create(
            id=3,
            ck_mechanism_file=fake_file
            )
        path_mech = mech.ck_mechanism_file.path
        name_mech = str(mech.ck_mechanism_file.name)
        self.assertEqual(os.path.join(MEDIA_ROOT, name_mech), path_mech)
        print(path_mech)
        os.remove(path_mech)
        if len(os.listdir(os.path.join(MEDIA_ROOT,'uploads/3/'))) == 0:
            os.rmdir(os.path.join(MEDIA_ROOT,'uploads/3/'))
        else:
            pass
