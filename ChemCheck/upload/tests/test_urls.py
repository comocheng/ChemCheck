"""
This test is aimed to test the connection of url and view function
"""
from django.test import SimpleTestCase
from django.urls import reverse, resolve
from upload.views import upload, mechanisms_list, MechanismDetailView, ck2yaml, MechanismUpdateView, MechanismDeleteView

class TestUrls(SimpleTestCase):
    def test_upload_url(self):
        upload_url = reverse('upload')
        print(resolve(upload_url))
        self.assertEqual(resolve(upload_url).func, upload)
    def test_list_url(self):
        list_url = reverse('list')
        print(resolve(list_url))
        self.assertEqual(resolve(list_url).func, mechanisms_list)
    def test_mechanism_detail_url(self):
        detail_url = reverse('mechanism-detail', args=['1'])
        print(resolve(detail_url))
        self.assertEqual(resolve(detail_url).func.view_class, MechanismDetailView)
    def test_convert_to_yaml_url(self):
        convert_to_yaml_url = reverse('ck2yaml', args=['1'])
        print(resolve(convert_to_yaml_url))
        self.assertEqual(resolve(convert_to_yaml_url).func, ck2yaml)
    def test_update_url(self):
        update_url = reverse('update_file_transport', args=['1'])
        print(resolve(update_url))
        self.assertEqual(resolve(update_url).func.view_class, MechanismUpdateView)
    def test_delete_url(self):
        delete_url = reverse('delete_file', args=['1'])
        print(resolve(delete_url))
        self.assertEqual(resolve(delete_url).func.view_class, MechanismDeleteView)