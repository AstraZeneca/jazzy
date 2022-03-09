Cookbook
========

Introduction
------------
What is Jazzy?
^^^^^^^^^^^^^^
Jazzy is a Python library that allows you to calculate a set of atomic/molecular descriptors which include the Gibbs free energy of hydration (kJ/mol), its polar/apolar components, and the hydrogen bond strength of donor and acceptor atoms using either SMILES or MOL/SDF inputs. Jazzy is easy to use, does not require expensive hardware, and produces accurate estimations within milliseconds to seconds for drug-like molecules. The library also exposes functionalities to depict molecules with atomistic hydrogen bond strengths in two or three dimensions.

What can I use Jazzy for?
^^^^^^^^^^^^^^^^^^^^^^^^^
The library was originally designed to support the processes of drug discovery and development but its applicability goes beyond life sciences. These are just some examples of the application of Jazzy:

* To score compounds based on their molecular or atomic hydrogen bond strengths or free energies of hydrations.
* To describe molecules for machine learning purposes (e.g. physicochemical or ADME modelling).
* To discriminate compounds on their C-H and X-H hydrogen bond donor strengths.
* To determine polar and apolar contributions of free energies of hydration.

What is this document for?
^^^^^^^^^^^^^^^^^^^^^^^^^^
This document contains a comprehensive list of examples on how to use Jazzy. The library provides three levels of programmatic access which require increasing programming skills: (1) command-line interface (CLI), (2) API functions, and (3) core functions. This architecture aims to maximise accessibility, ease of integration, and compatibility across different versions of Jazzy - but it also encourages the development of new functionalities and the reutilisation of existing components. As a note for developers, Jazzy relies on `RDKit <https://www.rdkit.org/>`_ and `kallisto <https://github.com/AstraZeneca/kallisto>`_ for the handling of molecule objects and the calculation of their atomic features.

Which functions should I use?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Depending on your programming expertise and what you are aiming to do with Jazzy, you should select an appropriate method to interact with it. Here, we are summarising what each level provides:

1. :ref:`cli-label`: Terminal commands to predict properties against an individual SMILES string or a file containing a set of molecules. The CLI is the easiest way to use Jazzy and ensures the highest compatibility across versions. This is the level you might want to use if you are, for example, interested in just running Jazzy against a data set of molecules.

2. :ref:`api-label`: Python functions that are simple to use directly or to integrate within scripts or other software. These methods ensure high compatibility across versions and high performance without the need to implement the calculation logic from scratch. This is how you might want to access Jazzy, for example, if you are implementing a script where molecules are described before feeding a machine learning regressor.

3. :ref:`core-label`: Python functions that allow fine-grained configuration of the parameters used to calculate the descriptors and direct control on RDKit and kallisto. These methods are not necessarily cross-compatible in different versions of Jazzy, and they need to be integrated using appropriate exception handling. This might how you want to use Jazzy if you are planning to use only some of its functionalities or if you have already your chemoinformatics methods in place.

Feedback and Contributions
^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to include new examples or review the existing ones, please refer to the `Contributor Guide`_ or submit a request to eike.caldeweyher@astrazeneca.com or gianmarco.ghiandoni@astrazeneca.com.

Examples
--------

.. _cli-label:

Command-line Interface
^^^^^^^^^^^^^^^^^^^^^^
.. include:: cli.rst

.. _api-label:

APIs
^^^^
.. include:: api.rst

.. _core-label:

Core Functions
^^^^^^^^^^^^^^
.. include:: core.rst

.. _Contributor Guide: https://jazzy.readthedocs.io/en/latest/contributing.html
