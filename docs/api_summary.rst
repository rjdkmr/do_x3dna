dnaMD Module
============

**dnaMD** is written in Python and it is also available as Python module. It
can be used in Python scripting for fast and flexible analysis
of `do_x3dna <do_x3dna_usage.html>`_ data.

**dnaMD** module contain a :class:`dnaMD.DNA` class and other functions.

.. currentmodule:: dnaMD

dnaMD.DNA class
~~~~~~~~~~~~~~~

.. autosummary::
    DNA.get_parameters
    DNA.set_base_pair_parameters
    DNA.set_helical_radius
    DNA.set_base_step_parameters
    DNA.get_mean_error
    DNA.set_helical_axis
    DNA.time_vs_parameter
    DNA.parameter_distribution
    DNA.set_major_minor_groove
    DNA.set_backbone_dihedrals
    DNA.generate_smooth_axis
    DNA.calculate_curvature_tangent
    DNA.calculate_angle_bw_tangents
    DNA.write_haxis_pdb


.. toctree::
   :maxdepth: 1

   DNA class <dna_class_api>


dnaMD module
~~~~~~~~~~~~
.. currentmodule:: dnaMD
.. autosummary::
    dnaMD.setParametersFromFile
    dnaMD.read_param_file
    dnaMD.get_error
    dnaMD.dev_parameters_vs_axis
    dnaMD.dev_bps_vs_parameter
    dnaMD.vector_angle


.. toctree::
   :maxdepth: 1

   Other functions <dnaMD_api>


dnaMD.dnaEY class
~~~~~~~~~~~~~~~~~

.. currentmodule:: dnaMD.dnaEY

.. autosummary::
    dnaEY.getStretchTwistBendModulus
    dnaEY.getStretchTwistModulus
    dnaEY.getModulusByTime
    dnaEY.calcDeformationEnergy

.. toctree::
   :maxdepth: 1

   dnaMD.dnaEY class for elasticity calculation <dnaEY_api>
