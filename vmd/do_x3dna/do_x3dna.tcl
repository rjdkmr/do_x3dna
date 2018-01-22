#
#
# This file is part of do_x3dna
#
# Author: Rajendra Kumar
# Copyright (C) 2014-2018  Rajendra Kumar
#
# do_x3dna uses 3DNA package (http://x3dna.org).
# Please cite the original publication of the 3DNA package:
# Xiang-Jun Lu & Wilma K. Olson (2003)
# 3DNA: a software package for the analysis, rebuilding and visualization of
# three-dimensional nucleic acid structures
# Nucleic Acids Res. 31(17), 5108-21.
#
# do_x3dna is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# do_x3dna is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with do_x3dna.  If not, see <http://www.gnu.org/licenses/>.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ============================================================================

package provide do_x3dna 1.0

source bigdcd.tcl

global eBasePairs eHbond  eLBP   eLBPS   eLBPH   eHelAxis eMgroove eHelixRad eBBnDihedral eSugarConf bS1 bS2
global True False
global bProp

set True 1
set False 0

# Boolean flags to read file
lassign {         0          1       2      3       4       5        6        7         8            9          10      11     } \
                  eBasePairs eHbond  eLBP   eLBPS   eLBPH   eHelAxis eMgroove eHelixRad eBBnDihedral eSugarConf bS1     bS2


proc parseArguments {args} {
  global True False

  global tp ftp dnaAtoms fitAtoms refFrame refFrameType trajs trajType dt suffix loadAll noisy
  set opts_value [list -tp -ftp -dnaAtoms -fitAtoms -refFrame -refFrameType -trajs -trajType -dt -suffix -help]
  set opts_bool [list  -loadAll -noisy -help]
  set opts [concat $opts_value $opts_bool]

  set i 0

  set loadAll $False
  set noisy $False

  set optLength [llength $args]

  while {$i < $optLength} {
    set key [lindex $args $i]
    incr i

    if {$key ni $opts_bool}  {
      set value [lindex $args $i]
      incr i
    }

    if {$key in $opts} {
      switch $key {
        -tp {set tp $value}
        -ftp {set ftp $value}
        -dnaAtoms {set dnaAtoms $value}
        -fitAtoms {set fitAtoms $value}
        -refFrame {set refFrame $value}
        -refFrameType {set refFrameType $value}
        -trajs {set trajs $value}
        -trajType {set dna $value}
        -dt {set dt $value}
        -suffix {set suffix $value}
        -noisy {set noisy $True}
        -loadAll {set loadAll $True}
        -help {do_x3dna_help
              exit 0}
      }
    } else {
      puts "========= ERROR ================"
      puts " $key is not a valid option.    "
      puts " See below for usage.           "
      puts "================================"
      do_x3dna_help
      exit 2
    }
  }

  # Set some default variables when not present
  if { ! [info exists dnaAtoms ] } {
    set dnaAtoms "nucleic"
  }

  if { ! [info exists dt ] } {
    set dt 1
  }

  if { ! [info exists suffix ] } {
    set suffix g
  }

  if { ! [info exists ftp ] } {
    set ftp "auto"
  }

  if { ! [info exists refFrameType ] } {
    set fframeType "auto"
  }

  if { ! [info exists trajType ] } {
    set trajType "auto"
  }

}

proc do_x3dna_help {} {
  puts {
    --------------------------------------------------------------------------------------------------
    Description
    ===========
    do_x3dna uses 3DNA package to calculate several structural descriptors or parameters
    of DNA/RNA using the MD trajectory. It extracts output of the 3DNA package
    and saves these parameters to external output files with function of time.

    To execute do_x3dna, 3DNA package should be installed and \$X3DNA environment
    variable (Detail is given in 3DNA manual) should be defined.

    NOTE 1: Before running do_x3dna, make sure 3DNA package is correctly working for
    the input DNA/RNA structure. To check it, run do_x3dna with only reference frame
    instead of trajectory file.
    "-noisy" option switch-on the output message from the 3DNA package on the display
    that would be necessary to troubleshoot problems related to the 3DNA package
    Most common problem is the residue names mismatch in input DNA/RNA structure
    and the 3DNA package dictionary.

    NOTE 2: Only PBC corrected trajectory  should be used as inputs.

    NOTE 3: If "-fitAtoms" is provided, this atoms group is used to superimpose current
    frame on reference frame. Most of the parameters are unaffected by this fitting, however
    coordinates of the local helical axis could mismatch with the input coordinates with the
    DNA/RNA.

    NOTE 4: By default, bigdcd (http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/bigdcd/)
    routine is used because it enables the trajectory reading frame by frame and allows to use
    very large trajectory. It can be switched-off by using "-loadAll" option.

    Also, all the above parameters including local helical axis, major and minor grooves,
    local helical radius, backbone dihedral angles (alpha, beta, gamma, delta, epsilon, zeta
    and chi) of both strands are calculated using 3DNA package for each frame and written in
    separate files as a function of time.

    OUTPUT FILE LIST:
      base_pairs_g.dat             => Base-pairs
      h-bond_g.dat                 => Hydrogen bonds between base-pairs
      L-BP_g.dat                   => Base-pairs parameters
      L-BPS_g.dat                  => Base-steps parameters
      L-BPH_g.dat                  => Helical Base-steps parameters
      HelAxis_g.dat                => Local helical axis coordinates
      MGroove_g.dat                => Major and Minor grooves
      HelixRad_g.dat               => Local helical radius
      BackBoneCHiDihedrals_g.dat   => Backbone dihedral angles including Chi-dihedral
      SugarDihedrals_g.dat         => Sugar dihedral angles including puckering type

    Name of these files could be change by setting different suffix instead of "g" using
    \"-suffix\" option. These files could be used with dnaMD for further analysis.

    USAGE
    =====
    do_x3dna -help -tp <file> -ftp <file type> -dnaAtoms <DNA atoms> -fitAtoms <Atoms to fit>
             -refFrame <Reference file> -refFrameType <reference file type> -suffix <output suffix>
             -trajs <List of trajectory files> -trajType <trajectory file type> -dt <time-step ps>
             -noisy -loadAll

    OPTIONS DESCRIPTION
    ===================

    -help          [Optional]   Show help and exit.

    -tp            [Required]   Topology-Parameter file with or without coordinates (psf/prmtop/pdb/gro).
                                If it is pdb/gro containing coordinates of atoms, this will be used as
                                reference frame for fitting and base-pair calculation. In case of coordinate
                                files such as pdb/gro, -refFrame will be only considered as an additional
                                frame as a part of trajectory.

    -ftp           [Optional]   File format of Topology-Parameter file. If not provided, file format
                                is automatically determined by VMD using extension.

    -dnaAtoms      [Optional]   Atoms belonging to DNA/RNA. If not provided, by default "nucleic" is
                                used to select atoms.

    -fitAtoms      [Optional]   Atoms for fitting to reference frame. If not provided, fitting will be
                                not done.

    -refFrame      [Required/   Reference coordinate file for fitting and base-pair calculation. In case
                    Optional]   of psf and prmtop file, this is required because these files do not
                                contain coordinates.

    -refFrameType  [Optional]   File format of Reference frame file. If not provided, file format
                                is automatically determined by VMD using extension.

    -trajs         [Optional]   List of trajectory files. If not provided, only first frame will be used
                                for calculation. It is useful to check whether 3DNA is working correctly
                                with -noisy option.

    -trajType      [Optional]   File format of trajectory file. If not provided, file format
                                is automatically determined by VMD using extension.

    -suffix        [Optional]   Suffix for all output file names. If not provided, "g" is used by default.

    -dt            [Optional]   Time-step between frames in ps. It is dumped in output files and make
                                easier for plotting. If not provided, by default, its value is 1 ps.

    -noisy         [Optional]   Display messages and errors from 3DNA programs find_pair and analyze.

    -loadAll       [Optional]   With this option, at first all trajectory will be loaded and then 3DNA
                                will be executed for each frame. Therefore, large memory will be consumed
                                for large trajectory.

    --------------------------------------------------------------------------------------------------
    }

}

proc do_x3dna {args} {

  # Variables for arguments
  global tp ftp dnaAtoms fitAtoms refFrame refFrameType trajs trajType dt suffix loadAll noisy

  # Variables
  global tmpFileBase tmpPDB inp3DNA find_pair_cmd analyze_cmd x3dnaOut outFileNames outFilePointers frameNum

  # vmd Atomselect variable
  global mol allAtomSelect dnaAtomSelect fitAtomSelect

  # check $X3DNA environment variable is available
  if { ! [info exists ::env(X3DNA) ] }  {
    puts "=================== ERROR ======================"
    puts " \$X3DNA environment variable is not available. "
    puts " See below for usage.                           "
    puts "================================================"
    do_x3dna_help
    exit 2
  }

  # Parse arguments and store in variables
  parseArguments {*}$args

  # Loading molecule
  if {$ftp == "auto"} {
    set mol [mol new $tp waitfor all]
  } else {
    set mol [mol new $tp type $ftp waitfor all]
  }

  # Check if topology-parameter file contains coordinate
  if { [molinfo top get numframes] == 0  && ! [info exists refFrame ] } {
    puts "=================== ERROR ======================"
    puts " Reference Frame is not provided. Use -refFrame "
    puts " and load reference/first coordinate file.      "
    puts " See below for usage.                           "
    puts "================================================"
    do_x3dna_help
    exit 2
  }

  # Load reference frame
  if { [info exists refFrame ] } {
    if {$refFrameType == "auto"} {
      mol addfile $refFrame waitfor all
    } else {
      mol addfile $refFrame type $refFrameType waitfor all
    }
  }

  # Set temporary filenames
  set tmpFileBase [randomString 8]
  set tmpPDB $tmpFileBase.pdb
  set inp3DNA $tmpFileBase.inp
  set x3dnaOut $tmpFileBase.out

  # set 3DNA command
  set find_pair_cmd [list [file join $::env(X3DNA) bin find_pair] $tmpPDB $inp3DNA ]
  set analyze_cmd [list [file join $::env(X3DNA) bin analyze] $inp3DNA ]
  if { ! $noisy } {
    lappend find_pair_cmd >/dev/null 2>/dev/null
    lappend analyze_cmd >/dev/null 2>/dev/null
  }

  # Initialize other variables
  set frameNum 0
  set outFileNames {}
  set outFilePointers {}

  # Generate output file names
  getOutFileNames

  # Initialize output files
  initOutFiles

  # run 3DNA for first frame
  run3DNA 0

  # Without trajectories end here
  if {! [info exists trajs ] } {
    remove3DnaTemporaryFiles
    exec rm $tmpPDB
    return 0
  }

  if { $loadAll } {
    # At first load all trajectories and run 3DNA
    foreach trj $trajs {
      if { $trajType == "auto" } {
          mol addfile $trj waitfor all
      } else {
          mol addfile $trj type $trajType waitfor all
      }
    }

    # run 3DNA for each frame
    set num_steps [molinfo $mol get numframes]
    set all [atomselect $mol all]
    for {set frame 1} {$frame < $num_steps} {incr frame} {
      animate goto $frame
      run3DNA $frame
    }
  } else {

    # If reference file contains one or more frames include here
    set num_steps [molinfo $mol get numframes]
    if { $num_steps > 1 } {
      for {set frame 1} {$frame < $num_steps} {incr frame} {
        animate goto $frame
        run3DNA $frame
      }
    }

    # Using bigdcd for trajectory
    bigdcd run3DNA $trajType {*}$trajs
    bigdcd_wait
  }

  # Remove temporary 3DNA files
  remove3DnaTemporaryFiles
  exec rm $tmpPDB
}

proc run3DNA { currentFrame }  {
  # Variables for arguments
  global tp ftp dnaAtoms fitAtoms refFrame refFrameType trajs trajType dt

  global tmpPDB inp3DNA find_pair_cmd analyze_cmd x3dnaOut outFileNames suffix outFilePointers frameNum

  # vmd Atomselect variable, set in main proc
  global mol allAtomSelect dnaAtomSelect fitAtomSelect

  # Information current frame
  set frameNum $currentFrame
  puts [format "\n CURRENT FRAME: %d\n" $currentFrame]

  # Fitting if fitting is enabled
  if { [info exists fitAtoms] } {
    set allAtomSelect [atomselect $mol all frame now]
    set ref [atomselect $mol $fitAtoms frame 0]
    set sel [atomselect $mol $fitAtoms frame now]
    $allAtomSelect move [measure fit $sel $ref]
  }

  # Select DNA atoms
  set dnaAtomSelect [atomselect $mol $dnaAtoms frame now]
  $dnaAtomSelect writepdb $tmpPDB

  # Run find_pair
  if {$currentFrame == 0} {
    if { [catch { exec {*}$find_pair_cmd } msg] } {
      puts $::errorInfo

      # Remove temporary 3DNA files
      remove3DnaTemporaryFiles

      exit 2
    }
  }

  # Run analyze
  if { [catch { exec {*}$analyze_cmd } msg] } {
    puts $::errorInfo

    # Remove temporary 3DNA files
    remove3DnaTemporaryFiles

    exit 2
  }

  # Write time to all output files
  writeTimeToFiles

  # Parse 3DNA output file and write extracted data to main output files
  parse3DnaOutput

}

proc getOutFileNames {}  {
  global eBasePairs eHbond eLBP   eLBPS   eLBPH  eHelAxis
  global eMgroove eHelixRad eBBnDihedral eSugarConf bS1  bS2

  global outFileNames suffix

  set i 0
  while {$i < 10} {
    lappend outFileNames dummy
    incr i
  }
  lset outFileNames $eBasePairs base_pairs_$suffix.dat
  lset outFileNames $eHbond h-bond_$suffix.dat
  lset outFileNames $eLBP L-BP_$suffix.dat
  lset outFileNames $eLBPS L-BPS_$suffix.dat
  lset outFileNames $eLBPH L-BPH_$suffix.dat
  lset outFileNames $eHelAxis HelAxis_$suffix.dat
  lset outFileNames $eMgroove MGroove_$suffix.dat
  lset outFileNames $eHelixRad HelixRad_$suffix.dat
  lset outFileNames $eBBnDihedral BackBoneCHiDihedrals_$suffix.dat
  lset outFileNames $eSugarConf SugarDihedrals_$suffix.dat

}

proc remove3DnaTemporaryFiles {}  {
  global inp3DNA
  # Remove all intermediate files
  exec rm $inp3DNA bestpairs.pdb hel_regions.pdb col_helices.scr col_chains.scr bp_step.par
  exec rm bp_helical.par stacking.pdb hstacking.pdb cf_7methods.par auxiliary.par bp_order.dat ref_frames.dat
}

proc initOutFiles {} {
  global eBasePairs eHbond  eLBP   eLBPS   eLBPH  eHelAxis
  global eMgroove eHelixRad eBBnDihedral eSugarConf bS1  bS2

  global outFileNames outFilePointers

  set i 0
  while {$i < 10} {
    lappend outFilePointers [open [lindex $outFileNames $i] w]
    incr i
  }

  puts [lindex $outFilePointers $eBasePairs] "#resid    chain   resid    chain\n"
  puts [lindex $outFilePointers $eHbond] "#Hydrogen bond counts\n"
  puts [lindex $outFilePointers $eLBP] "#Shear    Stretch   Stagger    Buckle  Propeller  Opening\n"
  puts [lindex $outFilePointers $eLBPS] "#Shift     Slide      Rise      Tilt      Roll     Twist\n"
  puts [lindex $outFilePointers $eLBPH] "#X-disp    Y-disp   h-Rise     Incl.       Tip   h-Twist\n"
  puts [lindex $outFilePointers $eHelAxis] "#Position (Px, Py, Pz) and local helical axis vector (Hx, Hy, Hz)for each dinucleotide step\n"
  puts [lindex $outFilePointers $eMgroove] "#Minor Groove        Major Groove\n"
  puts [lindex $outFilePointers $eMgroove] "#P-P     Refined     P-P     Refined\n"
  puts [lindex $outFilePointers $eHelixRad] "#   Strand I Atoms                Strand II Atoms\n"
  puts [lindex $outFilePointers $eHelixRad] "#P        O4'       C1'        P        O4'        C1'\n"
  puts [lindex $outFilePointers $eBBnDihedral] "#Strand I                                                    Strand II\n"
  puts [lindex $outFilePointers $eBBnDihedral] "#alpha    beta   gamma   delta  epsilon   zeta    chi   |||  alpha    beta   gamma   delta  epsilon   zeta    chi\n"
  puts [lindex $outFilePointers $eSugarConf] "#Strand I                                                              Strand II "
  puts [lindex $outFilePointers $eSugarConf] "#v0      v1      v2      v3      v4      tm       P    Puckering  |||  v0      v1      v2      v3      v4      tm       P    Puckering\n"
}

proc writeTimeToFiles {} {
  global eBasePairs eHbond  eLBP   eLBPS   eLBPH  eHelAxis
  global eMgroove eHelixRad eBBnDihedral eSugarConf bS1  bS2

  global dt frameNum outFilePointers outFileNames

  set i 0
  while {$i < 10} {
    set time_string [format "\n# Time = %15.5f\n" [expr $frameNum*$dt]]
    puts [lindex $outFilePointers $i] $time_string
    incr i
  }

}


proc parse3DnaOutput {}  {

  global True False
  global eBasePairs eHbond  eLBP   eLBPS   eLBPH  eHelAxis
  global eMgroove eHelixRad eBBnDihedral eSugarConf bS1  bS2

  global x3dnaOut outFilePointers

  set S1_Params {}
  set S2_Params {}

  set bProp {}
  set i 0
  while {$i < 12} {
    lappend bProp $False
    incr i
  }

  #  Slurp up the data file
  set finp [open $x3dnaOut r]
  #set fp [open "ddymEKpr.out" r]
  set file_data [read $finp]
  close $finp

  #  Process data file
  set data [split $file_data "\n"]
  foreach line $data {

    set line [string trim $line]

    # Read number of base-pairs from file
    if [regexp {Number of base-pairs} $line] {
      set num_bp [ expr [lindex [split $line ":"] 1] ]
    }


    if [regexp {\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*} $line] {
      lset bProp $eBasePairs $False
      lset bProp $eHbond $False
      lset bProp $eMgroove $False
      lset bProp $eHelixRad $False

      if [lindex $bProp $eBBnDihedral]	{
        foreach a $S1_Params b $S2_Params {
          set result [join [concat $a $b] "    "]
          puts [lindex $outFilePointers $eBBnDihedral] $result
        }
        set S1_Params {}
        set S2_Params {}
        lset bProp $eBBnDihedral $False
      }

      if [lindex $bProp $eSugarConf]	{
        foreach a $S1_Params b $S2_Params {
          set result [join [concat $a $b] "    "]
          puts [lindex $outFilePointers $eSugarConf] $result
        }
        set S1_Params {}
        set S2_Params {}
        lset bProp $eSugarConf $False
      }

    }

    if [regexp {~~~~~~~~~~~~~~~~~} $line] {
      lset bProp $eLBP $False
      lset bProp $eLBPS $False
      lset bProp $eLBPH $False
    }


    if [regexp {RMSD of the bases} $line] {
      lset bProp $eBasePairs $True
      continue
    }
    if [regexp {Detailed H-bond information} $line] {
      lset bProp $eHbond $True
      continue
    }
    if [regexp {Local base-pair parameters} $line] {
      lset bProp $eLBP $True
      continue
    }
    if [regexp {Local base-pair step parameters} $line] {
      lset bProp $eLBPS $True
      continue
    }
    if [regexp {Local base-pair helical parameters} $line] {
      lset bProp $eLBPH $True
      continue
    }
    if [regexp {Position \(Px, Py, Pz\) and local helical axis} $line] {
      lset bProp $eHelAxis $True
      continue
    }
    if [regexp {Minor and major groove widths} $line] {
      lset bProp $eMgroove $True
      continue
    }
    if [regexp {Helix radius} $line] {
      lset bProp $eHelixRad $True
      continue
    }
    if [regexp {Main chain and chi torsion angles} $line] {
      lset bProp $eBBnDihedral $True
      continue
    }
    if [regexp {Sugar conformational parameters} $line] {
      lset bProp $eSugarConf $True
      continue
    }

    if [lindex $bProp $eBBnDihedral]	{
      if { [regexp {Strand I} $line] && [lindex $bProp $bS1]==$False } {
        lset bProp $bS1 $True
        lset bProp $bS2 $False
        set bp 0
        continue
      }

      if { [regexp {Strand II} $line] &&  [lindex $bProp $bS2]==$False } {
        lset bProp $bS1 $False
        lset bProp $bS2 $True
        set bp 0
        continue
      }
    }

    if [lindex $bProp $eSugarConf]	{
      if { [regexp {Strand I} $line] &&  [lindex $bProp $bS1]==$False } {
        lset bProp $bS1 $True
        lset bProp $bS2 $False
        set bp 0
        continue
      }

      if { [regexp {Strand II} $line] &&  [lindex $bProp $bS2]==$False } {
        lset bProp $bS1 $False
        lset bProp $bS2 $True
        set bp 0
        continue
      }
    }

    if ![regexp {^[0-9]} $line]  {
      continue
    }

    if [lindex $bProp $eBasePairs]	{
      set temp [split $line  ":"]
      set temp1 [regexp -inline -- {([0-9]+)} [lindex $temp 1] ]
      set bp1 [expr [lindex $temp1 1] ]
      set temp1 [regexp -inline -- {([0-9]+)} [lindex $temp 3] ]
      set bp2 [expr [lindex $temp1 1] ]

      set temp1 [regexp -inline -- {([A-Za-z|\-]+$)} [lindex $temp 0] ]
      set chain1 [lindex $temp1 1]
      set temp1 [regexp -inline -- {(^[A-Za-z|\-])} [lindex $temp 4] ]
      set chain2 [lindex $temp1 1]

      puts [lindex $outFilePointers $eBasePairs] "$bp1   $chain1   $bp2   $chain2"
    }

    if [lindex $bProp $eHbond]	{
      set temp [regexp -inline -all -- {\S+} $line ]
      set temp1 [regexp -inline -- {(\d)} [lindex $temp 2] ]
      set hbond [expr [lindex $temp1 1] ]
      puts [lindex $outFilePointers $eHbond] $hbond
    }

    if [lindex $bProp $eLBP]	{
      set temp [regexp -inline -all -- {\S+} $line ]
      set result [join [lrange $temp 2 7] "    "]
      puts [lindex $outFilePointers $eLBP] $result
    }

    if [lindex $bProp $eLBPS]	{
      set temp [regexp -inline -all -- {\S+} $line ]
      set result [join [lrange $temp 2 7] "    "]
      puts [lindex $outFilePointers $eLBPS] $result
    }

    if [lindex $bProp $eLBPH]	{
      set temp [regexp -inline -all -- {\S+} $line ]
      set result [join [lrange $temp 2 7] "    "]
      puts [lindex $outFilePointers $eLBPH] $result
    }

    if [lindex $bProp $eMgroove]	{
      set temp [regexp -inline -all -- {\S+} $line ]
      set result [join [lrange $temp 2 5] "    "]
      puts [lindex $outFilePointers $eMgroove] $result
    }

    if [lindex $bProp $eHelAxis]	{
      set temp [regexp -inline -all -- {\S+} $line ]
      set result [join [lrange $temp 2 7] "    "]
      puts [lindex $outFilePointers $eHelAxis] $result
    }

    if [lindex $bProp $eHelixRad]	{
      set temp [regexp -inline -all -- {\S+} $line ]
      set result [join [lrange $temp 2 7] "    "]
      puts [lindex $outFilePointers $eHelixRad] $result
    }

    if [lindex $bProp $eBBnDihedral]	{
      if [lindex $bProp $bS1] {
        set temp [regexp -inline -all -- {\S+} $line ]
        set tmp_str [lrange $temp 2 8]
        lappend S1_Params $tmp_str
      }
      if [lindex $bProp $bS2] {
        set temp [regexp -inline -all -- {\S+} $line ]
        set tmp_str [lrange $temp 2 8]
        lappend S2_Params $tmp_str
      }
    }

    if [lindex $bProp $eSugarConf]	{
      if [lindex $bProp $bS1] {
        set temp [regexp -inline -all -- {\S+} $line ]
        set tmp_str [lrange $temp 2 9]
        lappend S1_Params $tmp_str
      }
      if [lindex $bProp $bS2] {
        set temp [regexp -inline -all -- {\S+} $line ]
        set tmp_str [lrange $temp 2 9]
        lappend S2_Params $tmp_str
      }
    }

  }
  # For loop end here

}

#
# Proc to generate a string of (given) characters
# Range defaults to "ABCDEF...wxyz'
#
proc randomString {length {chars "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"}} {
    set range [expr {[string length $chars]-1}]

    set txt ""
    for {set i 0} {$i < $length} {incr i} {
       set pos [expr {int(rand()*$range)}]
       append txt [string range $chars $pos $pos]
    }
    return $txt
}
