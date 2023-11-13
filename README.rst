.. image:: https://raw.githubusercontent.com/AstraZeneca/jazzy/master/docs/_static/jazzy.png
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
.. |Tests| image:: https://github.com/AstraZeneca/jazzy/workflows/Tests/badge.svg
   :target: https://github.com/AstraZeneca/jazzy/actions?workflow=Tests
   :alt: Tests
.. |Codecov| image:: https://codecov.io/gh/AstraZeneca/jazzy/branch/master/graph/badge.svg?token=4HCWYH61S5
   :target: https://codecov.io/gh/AstraZeneca/jazzy
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
*Jazzy* is an efficient computational tool for the calculation of hydration free energies and hydrogen-bond acceptor and donor strengths.
A publication describing the implementation, fitting, and validation of *Jazzy* can be found at `doi.org/10.1038/s41598-023-30089-x`_.

| If you are using *Jazzy* in your research, please remember to cite our publication as:
| *Ghiandoni, G.M., Caldeweyher, E. Fast calculation of hydrogen-bond strengths and free energy of hydration of small molecules. Sci Rep 13, 4143 (2023)*


Features
--------

* Hydration free energies
* Hydrogen-bond strengths
* Visualisation functions (SVG)
* Application programming interface
* Command-line interface


Installation via PyPI
---------------------

You can install *Jazzy* via pip_ from PyPI_:

.. code:: console

   $ pip install jazzy


Installation via Poetry (for development)
-----------------------------------------

Requirements to install *Jazzy* from the code repository:

- `poetry`_
- `pyenv`_ or `conda`_
- python>=3.6

First check that ``poetry`` is running (v1.4.1 at March 2023)

.. code:: console

   $ poetry --version
   Poetry version v1.4.1

Create a virtual environment (via ``pyenv`` or ``conda``) and activate it. Afterwards, clone the *Jazzy* project from GitHub and install it using ``poetry``.

.. code:: console

   $ git clone git@github.com:AstraZeneca/jazzy.git
   $ cd jazzy
   $ poetry install

If you wish to replicate our parameter fitting (see ``data/optuna_fitting`` and ``optimisation``), you need to install the specific version of ``optuna``.
At the time of the study, we used ``optuna==2.3.0``. You can install that manually using ``pip`` or ``poetry``.

.. code:: console

   $ poetry install --with optuna
   $ pip freeze | grep optuna
   optuna==2.3.0

Usage and Cookbook
------------------

Please see the `Usage <Usage_>`_ and `Cookbook <Cookbook_>`_ sections for details.


Contributing
------------

Jazzy is an open project in every shape and form, thus feedback on how to improve its documentation or functionalities is always welcome.
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
.. _poetry: https://python-poetry.org/docs/#installation
.. _pyenv: https://github.com/pyenv/pyenv#installation
.. _conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
.. _PyPI: https://pypi.org/
.. _Hypermodern Python Cookiecutter: https://github.com/cjolowicz/cookiecutter-hypermodern-python
.. _file an issue: https://github.com/AstraZeneca/jazzy/issues
.. _pip: https://pip.pypa.io/
.. _doi.org/10.1038/s41598-023-30089-x: https://doi.org/10.1038/s41598-023-30089-x
.. github-only
.. _Contributor Guide: contributing.rst
.. _Cookbook: https://jazzy.readthedocs.io/en/latest/cookbook.html
.. _Usage: https://jazzy.readthedocs.io/en/latest/usage.html
