How to use dnaMD tool?
======================

`do_x3dna <do_x3dna_usage.html>`_ generates `large amount of data <do_x3dna_usage.html#output-files-table>`_
and analyzing these data is difficult. Therefore, ``dnaMD`` tools are
developed to analyze data obtained from `do_x3dna output files <do_x3dna_usage.html#output-files-table>`_.
It contains set of tools to analyze the data.

These tools are particularly developed for the users with no programming experiences.
These can be used as commands and the resultant output files can be used
in external plotting program like xmgrace, gnuplot, excel etc to plot the results.

**Other tools are presently in development.**

.. list-table:: List of sub-commands available in dnaMD
    :widths: 1, 4
    :header-rows: 1
    :name: commands-table

    * - Command
      - Function

    * - `vsTime <commands/vsTime.html>`_
      - Extract a parameter as a function of time

    * - `saveAsH5 <commands/saveAsH5.html>`_
      - Save parameters to a HDF5 file

    * - `axisCurv <commands/axisCurv.html>`_
      - Calculate global helical-axis, curvatures and tangents from local helical axis.


.. toctree::
    vsTime : extract a parameter as a function of time <commands/vsTime>
    saveAsH5 : Save parameters to a HDF5 file <commands/saveAsH5>
    axisCurv : Calculate global helical-axis, curvatures and tangents <commands/axisCurv>
