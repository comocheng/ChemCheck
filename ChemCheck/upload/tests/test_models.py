''' This is a test for fields in models and the upload directory'''

from django.test import TestCase, Client
from upload.models import Mechanism, upload_to
from unittest import mock
from django.core.files import File
from django.core.files.uploadedfile import SimpleUploadedFile
from canteradebugger.settings import MEDIA_ROOT
import os
from django.contrib.auth.models import User
# from django.contrib.auth import authenticate, login

# user = User.objects.get(username='testuser')
# user = authenticate(username=user.username, password=user.password)
# user = login(user=user.username, password=user.password)
class MechanismTest(TestCase):
    # def SetUp(self):
    #     self.user.save()
    #     Client().login(username='testuser', password='Ab89898980')
    #     print(User.objects.all())
        #self.user = User.objects.get(username='testuser')
    def test_file_field(self):
        file_test = mock.Mock(spec=File)
        file_test.name = 'test.txt'
        user = User.objects.create_user(username='testuser', 
                                        email='testuser@gmail.com',
                                        password='Ab89898980')
        file_get_tested = Mechanism(user= user,
                                    ck_mechanism_file=file_test,
                                    ck_thermo_file=file_test,
                                    ck_transport_file=file_test,
                                    ck_surface_file=file_test,
                                    ct_mechanism_file=file_test,
                                    temperature = 298,
                                    pressure = 1e5)
        self.assertEqual(user, file_get_tested.user)
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
            'this is a test'.encode(),
            content_type='test/plain')
        user = User(username='testuser', password='Ab89898980')
        user.save()
        user = User.objects.get(username='testuser')
        mech = Mechanism.objects.create(user=user,
                                        id=3,
                                        ck_mechanism_file=fake_file)
        path_mech = mech.ck_mechanism_file.path
        name_mech = str(mech.ck_mechanism_file.name)
        self.assertEqual(os.path.join(MEDIA_ROOT, name_mech), path_mech)
        print(path_mech)
        os.remove(path_mech)
        if len(os.listdir(os.path.join(MEDIA_ROOT,'uploads/testuser/3/'))) == 0:
            os.rmdir(os.path.join(MEDIA_ROOT,'uploads/testuser/3/'))
