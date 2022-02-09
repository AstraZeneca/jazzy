.. image:: https://raw.githubusercontent.com/f3rmion/jazzy/main/assets/jazzy.png
  :width: 400
  :alt: Jazzy

|PyPI| |Status| |Python Version| |License|

|Read the Docs| |Tests| |Codecov|

|pre-commit| |Black|

.. |PyPI| image:: https://img.shields.io/pypi/v/jazzy.svg
   :target: https://pypi.org/project/jazzy/
   :alt: PyPI
.. |Status| image:: https://img.shields.io/pypi/status/jazzy.svg
   :target: https://pypi.org/project/jazzy/
   :alt: Status
.. |Python Version| image:: https://img.shields.io/pypi/pyversions/jazzy
   :target: https://pypi.org/project/jazzy
   :alt: Python Version
.. |License| image:: https://img.shields.io/pypi/l/jazzy
   :target: https://opensource.org/licenses/Apache-2.0
   :alt: License
.. |Read the Docs| image:: https://img.shields.io/readthedocs/jazzy/latest.svg?label=Read%20the%20Docs
   :target: https://jazzy.readthedocs.io/
   :alt: Read the documentation at https://jazzy.readthedocs.io/
.. |Tests| image:: https://github.com/f3rmion/jazzy/workflows/Tests/badge.svg
   :target: https://github.com/f3rmion/jazzy/actions?workflow=Tests
   :alt: Tests
.. |Codecov| image:: https://codecov.io/gh/f3rmion/jazzy/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/f3rmion/jazzy
   :alt: Codecov
.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit
.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Black

Full Author List
----------------
* Eike Caldeweyher
* Gian Marco Ghiandoni


Introduction
------------
*Jazzy* is an efficient computational tool for the calculation of hydration free energies and hydrogen bond acceptor and donor strengths.


Features
--------

* Hydration free energies
* Hydrogen bond strengths
* Visualisation functions (SVG)
* Application programming interface
* Command-line interface


Requirements
------------

.. code:: console

   click==7.1.2 Composable command line interface toolkit
   kallisto==1.0.7 Atomic and molecular featurizer
   numpy==1.21.1  NumPy is the fundamental package for array computing with Python.


Installation
------------

You can install *Jazzy* via pip_ from PyPI_:

.. code:: console

   $ pip install jazzy


Installation from Source
------------------------

Requirements to install *Jazzy* from sources:

- `poetry`_
- `pyenv`_ or `conda`_
- python >=3.6

First check that ``poetry`` is running correctly (v1.0.10 at the time of writing)

.. code:: console

    $ poetry --version
    Poetry version 1.0.10

Create a virtual environment (via ``pyenv`` or ``conda``) and activate it. Afterwards, clone the *Jazzy* project from GitHub and install it using ``poetry``

.. code:: console

    $ git clone git@github.com:f3rmion/jazzy.git
    $ cd jazzy
    $ poetry install


Usage
-----

Please see the `Command-line Reference <Usage_>`_ for details.


Contributing
------------

Contributions are very welcome.
To learn more, see the `Contributor Guide`_.


License
-------

Distributed under the terms of the `Apache 2.0 license`_,
*Jazzy* is free and open source software.


Issues
------

If you encounter any problems,
please `file an issue`_ along with a detailed description.


Credits
-------

This project was generated from `@cjolowicz`_'s `Hypermodern Python Cookiecutter`_ template.

.. _@cjolowicz: https://github.com/cjolowicz
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _Apache 2.0 license: https://opensource.org/licenses/Apache-2.0
.. _PyPI: https://pypi.org/
.. _Hypermodern Python Cookiecutter: https://github.com/cjolowicz/cookiecutter-hypermodern-python
.. _file an issue: https://github.com/f3rmion/jazzy/issues
.. _pip: https://pip.pypa.io/
.. github-only
.. _Contributor Guide: CONTRIBUTING.rst
.. _Usage: https://jazzy.readthedocs.io/en/latest/usage.html
.. _poetry: https://python-poetry.org/docs/#installation
.. _pyenv: https://github.com/pyenv/pyenv#installation
.. _conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
