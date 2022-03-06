Molecular Descriptors from SMILES
"""""""""""""""""""""""""""""""""
Example of calculation of Jazzy descriptors for an individual SMILES string from the command line. Features include *C-H donor strength* (sdc), *X-H donor strength* (sdx) where X includes any non-carbon atoms, *acceptor strength* (sda), *apolar contribution to delta g of hydration* (dga), *polar contribution to delta g of hydration* (dgp), *total delta G of hydration* (dgtot) which also accounts for an interaction term that is not included in the results.

.. code:: console

   $ jazzy vec 'NC1=CC=C(C=C1)O' --opt MMFF94

.. code:: console

   {'sdc': 2.2437, 'sdx': 2.111, 'sa': 1.999, 'dga': -3.432102037745161, 'dgp': -39.642437224760336, 'tot': -43.074539262505496, 'status': 'success', 'smiles': 'NC1=CC=C(C=C1)O'}

**TODO: CLI vis needs modifications**

**TODO: Need to make examples of how to read SDF of SMILES and get results (batches and individual). No programming language. Very high level**
