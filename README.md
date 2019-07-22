# ChemCheck

A tool to check chemistry, in the form of CHEMKIN and Cantera files.

Initial setup, installing dependencies using conda:

    $ conda create -n chemcheck -c conda-forge django\>=2.2.2 ruamel.yaml numpy django-crispy-forms
    $ conda activate chemcheck
    $ cd ChemCheck

Initial setup, installing dependencies using `pipenv`:

    $ pipenv install
    $ pipenv shell
    $ cd ChemCheck

First run or every time someone changes the models:

    $ python3 manage.py makemigrations
    $ python3 manage.py migrate

To launch the server:

    $ python3 manage.py runserver

then point a browser at http://127.0.0.1:8000/upload/
