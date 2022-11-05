Molecular Descriptors from SMILES
"""""""""""""""""""""""""""""""""
**Please note that the CLI functionalities are beta and will be subjected to changes in the future.**

Example of calculation of Jazzy descriptors for an individual SMILES string from the command line. Features include *C-H donor strength* (sdc), *X-H donor strength* (sdx) where X includes any non-carbon atoms, *acceptor strength* (sda), *apolar contribution to delta g of hydration* (dga), *polar contribution to delta g of hydration* (dgp), *total delta G of hydration* (dgtot) which also accounts for an interaction term that is not included in the results.

.. code:: console

   $ jazzy vec --opt MMFF94 'NC1=CC=C(C=C1)O'

.. code:: console

   {'sdc': 2.2437, 'sdx': 2.111, 'sa': 1.999, 'dga': -3.4321, 'dgp': -39.6424, 'tot': -43.0745, 'status': 'success', 'smiles': 'NC1=CC=C(C=C1)O'}


Atomistic Strength Visualisation from SMILES
""""""""""""""""""""""""""""""""""""""""""""
**Please note that the CLI functionalities are beta and will be subjected to changes in the future.**

Example of calculation of atomic hydrogen-bond strengths for an individual SMILES string from the command line. Creates an SVG string of the molecule with its atomistic hydrogen-bond donor and acceptor strengths within a dictionary from an input SMILES string.

.. code:: console

   $ jazzy vis --opt MMFF94 'NC1=CC=C(C=C1)O'

.. code:: console

   {'svg': "<?xml version='1.0' encoding='iso-8859-1'?>\n<svg version='1.1' baseProfile='full'\nxmlns='http://www.w3.org/2000/svg'\nxmlns:rdkit='http://www.rdkit.org/xml'\nxmlns:xlink='http://www.w3.org/1999/xlink'\nxml:space='preserve'\nwidth='500px' height='500px' viewBox='0 0 500 500'>\n<!-- END OF HEADER -->\n<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='500.0' height='500.0' x='0.0' y='0.0'> </rect>\n<path class='bond-0 atom-0 atom-1' d='M 406.4,222.5 L 369.6,230.9' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n<path class='bond-0 atom-0 atom-1' d='M 369.6,230.9 L 332.8,239.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n...</svg>\n", 'smiles': 'NC1=CC=C(C=C1)O', 'status': 'success'}
