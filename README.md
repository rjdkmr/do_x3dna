do_x3dna
========

do\_x3dna uses 3DNA package to calculate several structural descriptors or parameters of DNA/RNA using the GROMACS MD trajectory. It extracts output of the 3DNA package,and saves these parameters to external output files with function of time. 

<strong> Note: </strong> do\_x3dna can be used for analyzing DNA and RNA from the trajectory files that are obtained from other MD packages such as NAMD and AMBER. A PDB file could be used in place of a GROMACS _tpr_ file.

To execute do_x3dna, 3DNA package should be installed and $X3DNA environment variable (Detail is given in 3DNA manual) should be defined.

<strong> Please cite the follwoing publication and the current github page:</strong>               
Xiang-Jun Lu & Wilma K. Olson (2003)                                    
3DNA: a software package for the analysis, rebuilding and visualization 
of three-dimensional nucleic acid structures.                           
_Nucleic Acids Res._ 31(17), 5108-21.                                     

###Requirements

##### do_x3dna
* [3DNA package](http://x3dna.org)

* GROMACS libraries <code> libgmx, libmd, libgmxana </code> are required.

##### Python APIs
* numpy

***

###Download
<pre><code>git clone https://github.com/rjdkmr/do_x3dna
</code></pre>
***

###Installation
<pre><code>cd do_x3dna
mkdir build
cd build
cmake ..  -DGMX_PATH=/opt/gromacs -DCMAKE_INSTALL_PREFIX=/opt/do_x3dna
make
make install
</code></pre>

Directory <code>/opt/gromacs</code> should contains <code>include</code> and <code> lib </code> directories. If these directories are in seprate locations, use followings:
<pre><code>cmake ..  -DGMX_LIB=/path/to/lib -GMX_INCLUDE=/path/to/include -DCMAKE_INSTALL_PREFIX=/opt/do_x3dna
</code></pre>

If fftw library <code> libfftw3f.so or libfftw3f.a </code> are not present in standard locations:
<pre><code>-DFFTW_LIB=/path/to/fftw3/lib</code></pre>
***

###Usage
#### do_x3dna
<pre><code>do_x3dna -h
</code></pre>

<pre><code>Option     Filename  Type         Description
------------------------------------------------------------
  -f       traj.xtc  Input        Trajectory: xtc trr trj gro g96 pdb cpt
  -s      topol.tpr  Input        Structure+mass(db): tpr tpb tpa gro g96 pdb
  -n      index.ndx  Input, Opt.  Index file
  -o   BP_count.xvg  Output       xvgr/xmgr file
-map     BP_map.xpm  Output, Opt. X PixMap compatible matrix file
  -g        map.log  Output, Opt. Log file

Option       Type   Value   Description
------------------------------------------------------
-[no]h       bool   yes     Print help info and quit
-[no]version bool   no      Print version info and quit
-nice        int    19      Set the nicelevel
-b           time   0       First frame (ps) to read from trajectory
-e           time   0       Last frame (ps) to read from trajectory
-dt          time   0       Only use frame when t MOD dt = first time (ps)
-tu          enum   ps      Time unit: fs, ps, ns, us, ms or s
-xvg         enum   xmgrace  xvg plot formatting: xmgrace, xmgr or none
-[no]noisy   bool   no      Generate information from the X3DNA package
-[no]hbond   bool   no      Hydrogen bond map for base pairs
-[no]ref     bool   no      Base pair parameters will be calculated from base
                            pair of the reference frame
-name        string g       Output files will be prefixed by this string 
-[no]fit     bool   yes     Fitting frames on reference structure
-[no]mwa     bool   yes     Mass weighted fitting
-[no]lbpm    bool   no      To calculate local base pair parameters
-[no]lbpsm   bool   no      To calculate local base-pair step parameters
-[no]lbphm   bool   no      To calculate local base-pair helical parameters
-[no]avg     bool   yes     Average and Standard Deviation over all the frames
</code></pre>


do_x3dna uses 3DNA package to calculate several structural descriptors or parameters of DNA/RNA using the GROMACS MD trajectory. It extracts output of the 3DNA package,and saves these parameters to external output files with function of time.

To execute do_x3dna, 3DNA package should be installed and $X3DNA environment variable (Detail is given in 3DNA manual) should be defined.

<strong>NOTE 1:</strong> Before running do\_x3dna, make sure 3DNA package is correctly working for the input DNA/RNA structure. To check it, run do_x3dna with a input gro/pdb file containing only one frame instead of xtc/trr trajectory file as follows : 
<pre><code>do_x3dna -s topol.tpr -f check.gro -n index.ndx -noisy
</code></pre>
<code>-noisy</code> option switch on the output message from the 3DNA package on the display, that would be necessary to troubleshoot problems related to the 3DNA package. Most common problem is the residue names mismatch in input DNA/RNA structure and the 3DNA package dictionary.

<strong>NOTE 2:</strong> Only PBC corrected trajectory and tpr files should be used as inputs. PBC corrected PDB/GRO file can be used in place of tpr file. Index file <strong>SHOULD BE</strong> used to select ONLY DNA/RNA in the input trajectory file.

<strong>NOTE 3:</strong> If <code> -ref </code> option is used, all base-pairs/steps parameters will be calculated on the basis of the structure of the DNA/RNA present in the input tpr/pdb file with <code> -s </code> option. <code> -ref </code> option <strong> SHOULD BE </strong> used if output files are further used as input files in Python APIs or Scripts that are provided with do_x3dna package. To analyze the formation or breaking of base-pairs during MD simulations, either <code> -noref </code> could be used or do not include <code> -ref </code> option because this option is switched off by default.

<strong>NOTE 4:</strong> If <code> -fit </code> is enabled, during fitting procedure the DNA/RNA is translated to origin such that its center of mass is located at the origin. Most of the parameters are unaffected by this fitting, however coordinates of the local helical axis could mismatch with the input coordinates with the DNA/RNA.

<code> -hbond </code> option extracts hydrogen bonds for each base pair. A map.log ( <code>-g</code> )
file is generated containing the base pair information as per index of the hydrogen bond map ( <code>-map</code> ).

<code> -lbpm </code> option calculates Local Base Pair Parameters (Shear, Stretch, Stagger, Buckle, Propeller and Opening) with function of time, and average (with <code>-avg</code> ) of these parameters with function of the base-pairs.

<code> -lbpsm </code> option calculates Local Base Pair-Step Parameters (Shift, Slide, Rise, Tilt, Roll and Twist) with function of time, and average (with <code>-avg</code> ) of these parameters with function of the base-steps.

<code> -lbphm </code> option calculates Local Base Pair-Helical Parameters (X-displacement, Y-displacement, H-rise, Inclination, Tip and H-twist) with function of time, and average (with <code>-avg</code> ) of these parameters with function of the base-steps.

Also, all the above parameters including local helical axis, major and minor grooves, local helical radius, backbone dihedral angles (alpha, beta, gamma, delta, epsilon, zeta and chi) and sugar dihedral angles (v0, v1, v2, v3 and v4) of both strands are calculated using 3DNA package for each frame and written in separate files as a function of time.

<strong>OUTPUT FILE LIST:</strong>

| File name                 |   output contents                                                                   |
|:--------------------------|:------------------------------------------------------------------------------------|
|base_pairs_g.dat           | Base-pairs                                                                          |
|h-bond_g.dat               | hydrogen bonds between base-pairs                                                   |
|L-BP_g.dat                 | Base-pairs parameters                                                               |
|L-BPS_g.dat                | Base-steps parameters                                                               |
|L-BPH_g.dat                | Helical Base-steps parameters                                                       |
|HelAxis_g.dat              | Local helical axis coordinates                                                      |
|MGroove_g.dat              | Major and Minor grooves                                                             |
|HelixRad_g.dat             | Local helical radius                                                                |
|BackBoneCHiDihedrals_g.dat | Backbone dihederal angles including Chi-dihedral                                    |
|SugarDihedrals_g.dat       | Sugar dihederal angles including puckring type                                      |
|Stretch_g.xvg              | Stretch of base-pairs with function of time                                         |
|Shear_g.xvg                | Shear of base-pairs with function of time                                           |
|Stagger_g.xvg              | Stagger of base-pairs with function of time                                         |
|Buckle_g.xvg               | Buckle of base-pairs with function of time                                          |
|Propeller_g.xvg            | Propeller of base-pairs with function of time                                       |
|Opening_g.xvg              | Opening of base-pairs with function of time                                         |
|Shift_g.xvg                | Shift of base-steps with function of time                                           |
|Slide_g.xvg                | Slide of base-steps with function of time                                           |
|Rise_g.xvg                 | Rise of base-steps with function of time                                            |
|Tilt_g.xvg                 | Tilt of base-steps with function of time                                            |
|Roll_g.xvg                 | Roll of base-steps with function of time                                            |
|Twist_g.xvg                | Twist of base-steps with function of time                                           |
|X-displacement_g.xvg       | Helical X-displacment of helical base-steps with function of time                   |
|Y-displacement_g.xvg       | Helical Y-displacment of base-steps with function of time                           |
|Tip_g.xvg                  | Tip of base-steps with function of time                                             |
|Inclination_g.xvg          | Helical inclination of base-steps with function of time                             |
|H-twist_g.xvg              | Helical twist of helical base-steps with function of time                           |
|H-rise_g.xvg               | Helical rise of base-steps with function of time                                    |
|Avg_Local_BP_param_g.xvg   | Average and standard deviations of Base-pairs parameters with respect to base-pairs |
|Avg_bp_step_param_g.xvg    | Average and standard deviations of Base-steps parameters with respect to base-steps |
|Avg_bp_helical_param_g.xvg | Average and standard deviations of helical Base-steps parameters with respect to base-steps|

Name of these files could be change by setting different prefix instead of "g" using <code> -name </code> option. These files could be used with the Python APIs or scripts for further analysis.
***

#### Python APIs
These scripts are still in the development phase. Soon tutorial with the complete documentaion will be available. In the mean time, one could write own script to analyze the data obtained from the do\_x3dna.
***

