"""
This test is aimed to test the connection of url and view function
"""
from django.test import SimpleTestCase
from django.urls import reverse_lazy, resolve
from upload.views import (upload, mechanisms_list, mechanism_detail, ck2yaml, MechanismUpdateView, MechanismDeleteView, 
                          chemcheck, check_pdep_negative_A, check_negative_dup_rxns_negative_A, reaction_condition, collision_violation_check)                      

class TestUrls(SimpleTestCase):
    def test_upload_url(self):
        upload_url = reverse_lazy('upload:upload')
        print(resolve(upload_url))
        self.assertEqual(resolve(upload_url).func, upload)
    
    def test_list_url(self):
        list_url = reverse_lazy('upload:list')
        print(resolve(list_url))
        self.assertEqual(resolve(list_url).func, mechanisms_list)
    
    def test_mechanism_detail_url(self):
        detail_url = reverse_lazy('upload:mechanism_detail', args=[1])
        print(resolve(detail_url))
        self.assertEqual(resolve(detail_url).func, mechanism_detail)
    
    def test_convert_to_yaml_url(self):
        convert_to_yaml_url = reverse_lazy('upload:ck2yaml', args=[1])
        print(resolve(convert_to_yaml_url))
        self.assertEqual(resolve(convert_to_yaml_url).func, ck2yaml)
    
    def test_update_url(self):
        update_url = reverse_lazy('upload:update_file', args=[1])
        print(resolve(update_url))
        self.assertEqual(resolve(update_url).func.view_class, MechanismUpdateView)
    
    def test_delete_url(self):
        delete_url = reverse_lazy('upload:delete_file', args=[1])
        print(resolve(delete_url))
        self.assertEqual(resolve(delete_url).func.view_class, MechanismDeleteView)
    
    def test_chemcheck_url(self):
        chemcheck_url = reverse_lazy('upload:chemcheck', args=[1])
        print(resolve(chemcheck_url))
        self.assertEqual(resolve(chemcheck_url).func, chemcheck)
    
    def test_check_pdep_negative_A_url(self):
        check_pdep_negative_A_url = reverse_lazy('upload:check_pdep_negative_A', args=[1])
        print(resolve(check_pdep_negative_A_url))
        self.assertEqual(resolve(check_pdep_negative_A_url).func, check_pdep_negative_A)

    def test_check_negative_dup_rxns_negative_A_url(self):
        check_negative_dup_rxns_negative_A_url = reverse_lazy('upload:check_negative_dup_rxns_negative_A', args=[1])
        print(resolve(check_negative_dup_rxns_negative_A_url))
        self.assertEqual(resolve(check_negative_dup_rxns_negative_A_url).func, check_negative_dup_rxns_negative_A)

    def test_reaction_condition_url(self):
        reaction_condition_url = reverse_lazy('upload:reaction_condition', args=[1])
        print(resolve(reaction_condition_url))
        self.assertEqual(resolve(reaction_condition_url).func, reaction_condition)

    def test_collision_violation_check_url(self):
        collision_violation_check_url = reverse_lazy('upload:collision_violation', args=[1])
        print(resolve(collision_violation_check_url))
        self.assertEqual(resolve(collision_violation_check_url).func, collision_violation_check)
