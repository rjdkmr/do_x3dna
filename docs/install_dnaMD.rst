.. |numpy| raw:: html

  <a href="http://www.numpy.org/" target="_blank"> numpy </a>

.. |scipy| raw:: html

  <a href="http://www.scipy.org/" target="_blank"> scipy </a>

.. |h5py| raw:: html

  <a href="http://www.h5py.org/" target="_blank"> h5py </a>

.. |PyPI| raw:: html

  <a href="https://pypi.python.org/pypi" target="_blank"> PyPI </a>

Installation of dnaMD
=====================

dnaMD can be installed from source or directly from |PyPI| repository (pip).
We recommend to use pip.

**Package required after installation:** These packages are installed automatically
during dnaMD installation.

* |numpy|
* |scipy|
* |h5py|


Install using pip
-----------------

Execute following command to install with **Python-2.x** version:

.. code-block:: bash

  sudo pip install dnaMD


Execute following command to install with **Python-3.x** version:

.. code-block:: bash

  sudo pip3 install dnaMD


Updating dnaMD using pip
~~~~~~~~~~~~~~~~~~~~~~~~
To update the dnaMD use following command:

.. code-block:: bash

    pip install --upgrade --no-deps dnaMD


**OR**

.. code-block:: bash

    pip3 install --upgrade --no-deps dnaMD


``--upgrade`` flag is used to update the package and ``--no-deps`` prevents
update of dependent packages like numpy, scipy, matplotlib etc.


Install from source
-------------------

Download do_x3dna package as described `here <install_do_x3dna.html#download>`_.

Follow these steps to install with **Python-2.x** version:

.. code-block:: bash

    cd do_x3dna
    cd dnaMD
    sudo python setup.py install


Follow these steps to install with **Python-3.x** version:

.. code-block:: bash

    cd do_x3dna
    cd dnaMD
    sudo python3 setup.py install
