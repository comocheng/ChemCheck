# ChemCheck
A tool to check chemistry, in the form of chemkin and cantera files.


Initial setup, installing dependencies (or use the Pipfile):

    $ conda install django sqlparse ruamel_yaml

First run or every time someone changes your Models:

    $ python3 manage.py makemigrations
    $ python3 manage.py migrate

To launch the server:

    $ python3 manag.py runserver

then point a browser at http://127.0.0.1:8000/upload/
