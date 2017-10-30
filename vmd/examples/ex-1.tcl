# Required when do_x3dna is not installed as VMD plugin
# lappend auto_path [pwd]

package require do_x3dna

set pfile /home/rajendra/icelab/workspace/DNA_elasticity/fis_DNA/F1/analysis/complex_ion.pdb
set trajectories {/home/rajendra/icelab/workspace/DNA_elasticity/fis_DNA/F1/analysis/S1.xtc}

do_x3dna -tp $pfile -trajs $trajectories -fitAtoms "name P"
quit
