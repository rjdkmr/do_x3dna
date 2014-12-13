do_x3dna
========

do\_x3dna has been developed for analysis of the DNA/RNA dynamics during the molecular dynamics simulations. It uses [3DNA package](http://x3dna.org) to calculate several structural descriptors of DNA/RNA from the GROMACS MD trajectory. It executes 3DNA tools to calculate these descriptors and subsequently, extracts these output and saves in to external output files as a function of time.

From the MD trajectory, several [output files](http://rjdkmr.github.io/do_x3dna/usage.html#output-files) are generated. For easy analysis of the obtained data, do\_x3dna package includes a Python code [dnaMD.py](https://github.com/rjdkmr/do_x3dna/blob/master/Python_API/dnaMD.py), which contains several [methods](http://rjdkmr.github.io/do_x3dna/apidoc.html) for the analysis of structural descriptors.


<strong> Note: </strong> do\_x3dna can be used with the trajectory files that are obtained from other MD packages such as NAMD and AMBER. Input trajectory files should be converted in to Gromacs format trajectory files. A PDB file could be used in place of a GROMACS _tpr_ file.

To execute do_x3dna, 3DNA package should be installed and `$X3DNA` environment variable (Detail is given in 3DNA manual) should be defined.

**Last Update: 12 Dec. 2014**

<strong> For detailed documentation about the do\_x3dna, please visit  [do\_x3dna home-page](http://rjdkmr.github.io/do_x3dna). </strong>


<strong> Please cite the follwoing publication and the current github page:</strong>               
Xiang-Jun Lu & Wilma K. Olson (2003)                                    
3DNA: a software package for the analysis, rebuilding and visualization 
of three-dimensional nucleic acid structures.                           
_Nucleic Acids Res._ 31(17), 5108-21.                                     
