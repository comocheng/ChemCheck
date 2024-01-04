# ChemCheck
## Please use branch `fix_auth_bugs` for now, the `main` branch is not working at this point
This project creates a web-based application to visualize and diagnose errors in the conversion of CHEMKIN files to Cantera input files.
ChemCheck is developed in Django 2.2 and python 3.6, and it allows users edit and convert files on the website. The current version can
report error messages and codes around the error line when conversion fails. In this moment, we have not covered enough conversion errors.

Next Step: we will test different error files reported in Cantera users group to provide more error diagnosis.

### Initital setup

Initial setup, installing dependencies using conda:
```
    $ cd ChemCheck
    $ conda env create -f environment.yml
    $ source activate chemcheck_env
    $ cd ChemCheck
    $ python manage.py makemigrations
    $ python manage.py migrate
    $ python manage.py test
```
### Updates

First run or every time someone changes the models:

    $ python manage.py makemigrations
    $ python manage.py migrate

### Running

To launch the server:

    $ python manage.py runserver

then point a browser at http://127.0.0.1:8000/upload/

### Citation

If you want to cite this project, please use https://doi.org/10.17760/D20399888
```
@phdthesis{Xu_2020, place={Boston}, title={ChemCheck : a cantera debugging tool to detect chemical and syntax errors in kinetic models}, author={Xu, Chao}, year={2020}, pages={1–53}} 
```
