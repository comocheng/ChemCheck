from django.db import models
from django.utils import timezone
import datetime
import yaml
import math
import numpy as np
import matplotlib.pyplot as plt
import os


# Create your models here.

def upload_to(instance, filename):
    return 'uploads/{id}/{fn}'.format(id=instance.pk,fn=filename)


class Mechanism(models.Model):
    """
    A chemical kinetic mechanism, from Chemkin or Cantera
    """
    ck_mechanism_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True,
                                       verbose_name="Chemkin mechanism file")
    ck_thermo_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True,
                                       verbose_name="Chemkin thermo file")
    ck_transport_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True,
                                       verbose_name="Chemkin transport file")
    ck_surface_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True,
                                       verbose_name="Chemkin surface file")
    ct_mechanism_file = models.FileField(upload_to=upload_to, max_length=100, blank=True, null=True,
                                        verbose_name='Cantera yaml file')
    ct_conversion_errors = models.TextField(verbose_name='Errors from the ck2yaml conversion')
    timestamps = models.DateTimeField(auto_now_add=True)

    def get_absolute_url(self):
        return 'uploads/{self.id}/'
    


    def save(self, *args, **kwargs):
        """ 
        The folder used to upload the files to, depends on the id (primary key)
        of the Mechanism object. If that is newly created and not yet saved, it 
        doesn't have an id yet. So to save the Mechanism object and all its files,
        you need to save the Mechanism object first with no files if it has no id,
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

# class Images(models.Model):
#     """
#     Images with polynomial discontinuities
#     """
#     mechanism = models.ForeignKey(Mechanism, on_delete=models.CASCADE)
#     discontinuous_species = models.ImageField(upload_to=upload_to, verbose_name='discontinuous_species_images')

def cp_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high):
    if T < T_mid:
        cp_low = Nasa_poly_low[0] + Nasa_poly_low[1] * T + Nasa_poly_low[2] * T ** 2 + Nasa_poly_low[3] * T ** 3 + Nasa_poly_low[4] * T ** 4
        return cp_low
    elif T > T_mid:
        cp_high = Nasa_poly_high[0] + Nasa_poly_high[1] * T + Nasa_poly_high[2] * T ** 2 + Nasa_poly_high[3] * T ** 3 + Nasa_poly_high[4] * T ** 4
        return cp_high
    else:
        cp_low = Nasa_poly_low[0] + Nasa_poly_low[1] * T + Nasa_poly_low[2] * T ** 2 + Nasa_poly_low[3] * T ** 3 + Nasa_poly_low[4] * T ** 4
        cp_high = Nasa_poly_high[0] + Nasa_poly_high[1] * T + Nasa_poly_high[2] * T ** 2 + Nasa_poly_high[3] * T ** 3 + Nasa_poly_high[4] * T ** 4
        if cp_low - cp_high == 0:
            return cp_low
        else:
            return cp_low, cp_high

def h_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high):
    if T < T_mid:
       h_low = Nasa_poly_low[0] + Nasa_poly_low[1] / 2 * T + Nasa_poly_low[2] / 3 * T ** 2 + Nasa_poly_low[3] / 4 * T ** 3 + Nasa_poly_low[4] / 5 * T ** 4 + Nasa_poly_low[5] / T
       return h_low
    elif T > T_mid:
        h_high = Nasa_poly_high[0] + Nasa_poly_high[1] / 2 * T + Nasa_poly_high[2] / 3 * T ** 2 + Nasa_poly_high[3] / 4 * T ** 3 + Nasa_poly_high[4] / 5 * T ** 4 + Nasa_poly_high[5] / T
        return h_high
    else:
        h_low = Nasa_poly_low[0] + Nasa_poly_low[1] / 2 * T + Nasa_poly_low[2] / 3 * T ** 2 + Nasa_poly_low[3] / 4 * T ** 3 + Nasa_poly_low[4] / 5 * T ** 4 + Nasa_poly_low[5] / T
        h_high = Nasa_poly_high[0] + Nasa_poly_high[1] / 2 * T + Nasa_poly_high[2] / 3 * T ** 2 + Nasa_poly_high[3] / 4 * T ** 3 + Nasa_poly_high[4] / 5 * T ** 4 + Nasa_poly_high[5] / T
        if h_low - h_high == 0:
            return h_low
        else:
            return h_low, h_high

def s_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high):
    if T < T_mid:
       s_low = Nasa_poly_low[0] * math.log(T) + Nasa_poly_low[1] * T + Nasa_poly_low[2] / 2 * T ** 2 + Nasa_poly_low[3] / 3 * T ** 3 + Nasa_poly_low[4] / 4 * T ** 4 + Nasa_poly_low[6]
       return s_low
    elif T > T_mid:
        s_high = Nasa_poly_high[0] * math.log(T) + Nasa_poly_high[1] * T + Nasa_poly_high[2] / 2 * T ** 2 + Nasa_poly_high[3] / 3 * T ** 3 + Nasa_poly_high[4] / 4 * T ** 4 + Nasa_poly_high[6]
        return s_high
    else:
        s_low = Nasa_poly_low[0] * math.log(T) + Nasa_poly_low[1] * T + Nasa_poly_low[2] / 2 * T ** 2 + Nasa_poly_low[3] / 3 * T ** 3 + Nasa_poly_low[4] / 4 * T ** 4 + Nasa_poly_low[6]
        s_high = Nasa_poly_high[0] * math.log(T) + Nasa_poly_high[1] * T + Nasa_poly_high[2] / 2 * T ** 2 + Nasa_poly_high[3] / 3 * T ** 3 + Nasa_poly_high[4] / 4 * T ** 4 + Nasa_poly_high[6]
        if s_low - s_high == 0:
            return s_low
        else:
            return s_low, s_high

class ChemError:
    def __init__(self, path, name):
        self.path = path
        self.name = name
        
    
    def check_continuity(self):
        discontinuous_species = []
        with open(self.path, 'r') as f:
            chem_data = yaml.load(f, Loader=yaml.FullLoader)
        for species in chem_data['species']:
            #T_low = species['thermo']['temperature-ranges'][0]
            T_mid = species['thermo']['temperature-ranges'][1]
            #T_high = species['thermo']['temperature-ranges'][2]
            Nasa_poly_high = species['thermo']['data'][1]
            Nasa_poly_low = species['thermo']['data'][0]
            cp_low = Nasa_poly_low[0] + Nasa_poly_low[1] * T_mid + Nasa_poly_low[2] * T_mid ** 2 + Nasa_poly_low[3] * T_mid ** 3 + Nasa_poly_low[4] * T_mid ** 4
            cp_high = Nasa_poly_high[0] + Nasa_poly_high[1] * T_mid + Nasa_poly_high[2] * T_mid ** 2 + Nasa_poly_high[3] * T_mid ** 3 + Nasa_poly_high[4] * T_mid ** 4
            h_low = Nasa_poly_low[0] + Nasa_poly_low[1] / 2 * T_mid + Nasa_poly_low[2] / 3 * T_mid ** 2 + Nasa_poly_low[3] / 4 * T_mid ** 3 + Nasa_poly_low[4] / 5 * T_mid ** 4 + Nasa_poly_low[5] / T_mid
            h_high = Nasa_poly_high[0] + Nasa_poly_high[1] / 2 * T_mid + Nasa_poly_high[2] / 3 * T_mid ** 2 + Nasa_poly_high[3] / 4 * T_mid ** 3 + Nasa_poly_high[4] / 5 * T_mid ** 4 + Nasa_poly_high[5] / T_mid
            s_low = Nasa_poly_low[0] * math.log(T_mid) + Nasa_poly_low[1] * T_mid + Nasa_poly_low[2] / 2 * T_mid ** 2 + Nasa_poly_low[3] / 3 * T_mid ** 3 + Nasa_poly_low[4] / 4 * T_mid ** 4 + Nasa_poly_low[6]
            s_high = Nasa_poly_high[0] * math.log(T_mid) + Nasa_poly_high[1] * T_mid + Nasa_poly_high[2] / 2 * T_mid ** 2 + Nasa_poly_high[3] / 3 * T_mid ** 3 + Nasa_poly_high[4] / 4 * T_mid ** 4 + Nasa_poly_high[6]
            if abs(cp_high - cp_low) <= 0.1 and abs(h_low - h_high) <= 0.1 and abs(s_low - s_high) <= 0.1:
                pass
            else:
                image = '{}.png'.format(species['name'])
                img_path = os.path.join(os.path.split(self.path)[0], image)
                img_url = img_path.split('ChemCheck')[2]
                discontinuous_species.append(img_url)
                #gas = ct.Solution(self.path)
                #sp = gas.species(species['name'])
                T = np.linspace(500, 3000, 200)
                fig,ax = plt.subplots(1,3,figsize=(8,3.5))
                cp = [cp_calculate(tt, T_mid, Nasa_poly_low, Nasa_poly_high) for tt in T]
                h = [h_calculate(tt, T_mid, Nasa_poly_low, Nasa_poly_high) for tt in T]
                s = [s_calculate(tt, T_mid, Nasa_poly_low, Nasa_poly_high) for tt in T]
                ax[0].plot(T,cp)
                ax[0].set_title('$c_p/R$')
                ax[1].plot(T,h)
                ax[1].set_title('$h/RT$')
                ax[2].plot(T,s)
                ax[2].set_title('$s/R$')
                fig.tight_layout()
                fig.suptitle('{} thermo property'.format(species['name']), y=1)
                fig.savefig(img_path)
        return discontinuous_species


                # T_mid = symbols('x')
                # cp_low = Nasa_poly_low[0] + Nasa_poly_low[1] * T_mid + Nasa_poly_low[2] * T_mid ** 2 + Nasa_poly_low[3] * T_mid ** 3 + Nasa_poly_low[4] * T_mid ** 4
                # cp_high = Nasa_poly_high[0] + Nasa_poly_high[1] * T_mid + Nasa_poly_high[2] * T_mid ** 2 + Nasa_poly_high[3] * T_mid ** 3 + Nasa_poly_high[4] * T_mid ** 4
                # h_low = Nasa_poly_low[0] + Nasa_poly_low[1] / 2 * T_mid + Nasa_poly_low[2] / 3 * T_mid ** 2 + Nasa_poly_low[3] / 4 * T_mid ** 3 + Nasa_poly_low[4] / 5 * T_mid ** 4 + Nasa_poly_low[5] / T_mid
                # h_high = Nasa_poly_high[0] + Nasa_poly_high[1] / 2 * T_mid + Nasa_poly_high[2] / 3 * T_mid ** 2 + Nasa_poly_high[3] / 4 * T_mid ** 3 + Nasa_poly_high[4] / 5 * T_mid ** 4 + Nasa_poly_high[5] / T_mid
                # s_low = Nasa_poly_low[0] * ln(T_mid) + Nasa_poly_low[1] * T_mid + Nasa_poly_low[2] / 2 * T_mid ** 2 + Nasa_poly_low[3] / 3 * T_mid ** 3 + Nasa_poly_low[4] / 4 * T_mid ** 4 + Nasa_poly_low[6]
                # s_high = Nasa_poly_high[0] * ln(T_mid) + Nasa_poly_high[1] * T_mid + Nasa_poly_high[2] / 2 * T_mid ** 2 + Nasa_poly_high[3] / 3 * T_mid ** 3 + Nasa_poly_high[4] / 4 * T_mid ** 4 + Nasa_poly_high[6]
                # T_cp = solveset(Abs(cp_high - cp_low) < 1, T_mid, domain=S.Reals)
                # T_h = solveset(Abs(h_high - h_low) < 1, T_mid, domain=S.Reals)
                # T_s = solveset(Abs(s_high - s_low) < 1, T_mid, domain=S.Reals)
                # T_range = Interval(T_low, T_high)
                # T_continuous = T_cp.intersect(T_h).intersect(T_s).intersect(T_range)
                # T_continuouse = str(T_continuous)[str(T_continuous).find('Interval'):len(str(T_continuous))]
                # if not T_continuouse:
                #     print(species['name'])
                # else:
                #     print(species['name'], T_continuouse)
