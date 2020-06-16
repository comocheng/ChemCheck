from django.test import TestCase, Client
from django.urls import reverse_lazy
from django.core.files.uploadedfile import SimpleUploadedFile, UploadedFile
from upload.models import Mechanism
import os
from canteradebugger.settings import MEDIA_ROOT
from django.contrib.auth.models import User

        

class TestViews(TestCase):
    def SetUp(self):
        self.client = Client()
    def test_upload_view(self):
        user = User.objects.create_user(username='testuser', 
                                        email='testuser@gmail.com',
                                        password='Ab89898980')
        self.client.login(username='testuser', password='Ab89898980')
        fake_file = SimpleUploadedFile(
            'test-file.txt',
            'this is what\'s in the file that isn\'t a file'.encode(),
            content_type='text/plain')
        post_response = self.client.post(reverse_lazy('upload:upload', args=['testuser']),{
                        'ck_mechanism_file':fake_file
        })
        post_file_path = str(Mechanism.objects.get(id=1).ck_mechanism_file)
        self.assertEquals(os.path.split(post_file_path)[1], 'test-file.txt')
        self.assertEquals(post_response.status_code, 302)
    
    #start of detail view test
        get_response = self.client.get(reverse_lazy('upload:mechanism_detail', args=['testuser', 1]))
        self.assertEquals(get_response.status_code, 200)
        self.assertTemplateUsed(get_response, 'upload/mechanism_detail.html', 'base.html')
    
    #start of ck2yaml view test
        ck2yaml_response = self.client.get(reverse_lazy('upload:ck2yaml', args=['testuser', 1]))
        self.assertEqual(ck2yaml_response.status_code, 200)
        self.assertTemplateUsed(ck2yaml_response, 'ck2yaml.html', 'base.html')
        self.assertTrue(os.path.isfile(str(MEDIA_ROOT)+'/uploads/testuser/1/error.txt'))
        os.remove(str(MEDIA_ROOT)+'/uploads/testuser/1/error.txt')
    #start of update view test
        update_file = SimpleUploadedFile(
            'update-file.txt',
            'Only test'.encode(),
            content_type='test/plain'
        )
        update_get_response = self.client.get(reverse_lazy('upload:update_file_transport', args=['testuser', 1]))
        self.assertEquals(update_get_response.status_code, 200)
        self.assertTemplateUsed(update_get_response, 'base.html', 'file_update.html')

        update_post_response = self.client.post(reverse_lazy('upload:update_file_transport', args=['testuser', 1]), {
            'ck_mechanism_file':update_file
        })
        update_path = str(Mechanism.objects.get(id=1).ck_mechanism_file.path)
        self.assertEquals(update_post_response.status_code, 302)
        self.assertEquals(os.path.split(update_path)[1], 'update-file.txt')
    
    #start of delete view test
        delete_get_response = self.client.get(reverse_lazy('upload:delete_file', args=['testuser', 1]))
        self.assertEquals(delete_get_response.status_code, 200)
        self.assertTemplateUsed(delete_get_response, 'base.html', 'file_delete_mechanism.html')
        
        delete_post_response = self.client.post(reverse_lazy('upload:delete_file', args=['testuser', 1]))
        self.assertEquals(delete_post_response.status_code, 302)
        self.assertFalse(Mechanism.objects.get(id=1).ck_mechanism_file)
        self.assertFalse(os.path.isfile(str(MEDIA_ROOT)+'/uploads/testuser/1/update-file.txt'))
    
    #start of discontinuity check test
        file_dir = os.path.join(os.getcwd(), 'upload/tests/cantera.yaml')
        with open(file_dir, 'rb') as f:
            test_file = SimpleUploadedFile('cantera.yaml',
                                           f.read(),
                                           content_type='text/yaml')
        mech = Mechanism.objects.create(user=user,
                                        id=2,
                                        ct_mechanism_file=test_file
                                        )
        chemcheck_get_response = self.client.get(reverse_lazy('upload:chemcheck', args=['testuser', 2]))
        self.assertEqual(chemcheck_get_response.status_code, 200)
    
    # start of pdep negative A test
        get_pdep_neg_A_response = self.client.get(reverse_lazy('upload:check_pdep_negative_A',  args=['testuser', 2]))
        self.assertEqual(get_pdep_neg_A_response.status_code, 200)
       
    
    # start of duplicate reaction negative A test
        get_dup_neg_A_response = self.client.get(reverse_lazy('upload:check_negative_dup_rxns_negative_A',  args=['testuser', 2]))
        self.assertEqual(get_dup_neg_A_response.status_code, 200)
    
    # start of reaction condition test
        get_rxn_condition_response = self.client.get(reverse_lazy('upload:reaction_condition',  args=['testuser', 2]))
        self.assertEqual(get_rxn_condition_response.status_code, 200)
    
    # start of collision violation check test
        get_collision_check_response = self.client.get(reverse_lazy('upload:collision_violation',  args=['testuser', 2]))
        self.assertEqual(get_collision_check_response.status_code, 200)
        path_mech = mech.ct_mechanism_file.path
        os.remove(path_mech)
        folder_path = os.path.join(MEDIA_ROOT,'uploads/testuser/2/')

        if len(os.listdir(os.path.join(MEDIA_ROOT,'uploads/testuser/2/'))) == 0:
            os.rmdir(os.path.join(MEDIA_ROOT,'uploads/testuser/2/'))
        if len(os.listdir(os.path.join(MEDIA_ROOT,'uploads/testuser/'))) == 0:
            os.rmdir(os.path.join(MEDIA_ROOT,'uploads/testuser/1/'))
            User.objects.get(username='testuser').delete()
