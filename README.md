Copyright 2010, Torbjorn Bjorkman


## Python2.7


## LAMMPS Support
I have added support for LAMMPS .data files with orthogonal and triclinic boxes. To convert a CIF-file to LAMMPS data file, simply run (in this directory)
```
./cif2cell /path/to/file.cif --no-reduce -p lammps -o cif_file.data
```

## CIF2CELL

A tool to generate the geometrical setup for various electronic 
structure codes from a CIF (Crystallographic Information 
Framework) file. The code will generate the crystal structure for 
the primitive cell or the conventional cell.

## CURRENTLY SUPPORTS

| code             | alloy support | output files (--setup-all) |
| -----------------|:-------------:|--------------|
| ASE              |  no           | positions.py |
| ATAT             |  no           | [compoundname].in |
| CPMD             |  no           | [compoundname].inp |
| CP2k             |  no           | [compoundname].inp |
| CASTEP           | VCA           | [compoundname].cell |
| Crystal09        |  no           | [compoundname].d12 |
| VASP             | VCA           | POSCAR (INCAR, POTCAR and KPOINTS) |
| ABINIT           |  no           | [compoundname].in ([compoundname].files) |
| quantum espresso |  no           | [compoundname].in ([compoundname].in) |
| Siesta           |  no           | [compoundname].fdf |
| OpenMX           |  no           | [compoundname].dat ([compoundname].dat)|
| FDMNES           |  no           | [compoundname].txt (fdmfile.txt) |
| FHI-aims         |  no           | geometry.in |
| RSPt             |  no           | symt.inp |
| Fleur            |  no           | inp_[compoundname] |
| elk              |  no           | GEOMETRY.OUT (GEOMETRY.OUT)|
| exciting         |  no           | input.xml |
| hutsepot         |  no           | [compoundname].sys |
| cellgen          |  no           | cellgen.inp |
| spacegroup       |  no           | spacegroup.in |
| ncol             |  no           | [spacegroupname/compoundname].dat  for bstr. |
| emto             |  yes          | [spacegroupname/compoundname].dat for kstr, bmdl, shape, kgrn and kfcd in separate directories. |
| Akai-KKR         |  yes          | [compoundname].in ([compoundname].in) |
| SPR-KKR          |  yes          | [compoundname].sys |
| MOPAC            |  no           | [compoundname].mop |
| DFTB+            |  no           | [compoundname].gen |
| Lammps           |  no           | [compoundname].data |
| xyz              |  no           | [compoundname].xyz |
| cfg              |  no           | [compoundname].cfg |
| coo              |  no           | [compoundname].coo |
| spc              |  no           | [compoundname].dat |

## CONTENTS
The distribution includes:
* This README file.
* The file LICENSE with the GPLv3 license.
* The python files cif2cell, uctools.py and spacegroupdata.py
* Installation files, setup.py and MANIFEST.
* A manual.
* The directory cifs/ containing a set of example CIF files 
  as well as the crystal structures of the full periodic table 
  from COD, the Crystallography Open Database <http://www.crystallography.net>
  and also a few from ICSD (with permission).
* The file PyCifRW-3.3.tar.gz, containing the PyCifRW package needed for
  parsing CIF files.


## INSTALLATION INSTRUCTIONS

Prerequisites: The program requires Python 2.4 or higher and the
               PyCIFRW python package (which will be installed 
               automatically if not present, see below for manual 
   	       installation instructions). Note however that the output
               may be slightly different (but formally equivalent) 
	       with Python 2.4 than with later versions.

To install the program in your systems standard location, simply type:
python setup.py install 
To choose a different location, add 
--prefix=where/you/want/it 
to the above line. For help and more options type
python setup.py --help

The installation will also create a directory $PREFIX/lib/cif2cell
that contains the manual and sample cif files.


## DOCUMENTATION

The setup will install the manual, cif2cell.pdf, into the 
$PREFIX/lib/cif2cell/docs directory. 


## RUNNING

Run 'cif2cell -h' to get a listing of the different options.
Example:
cif2cell Ni20Mn3P6.cif -p vasp --vasp-cartesian-positions
will generate a POSCAR file for VASP with the positions in cartesian format.


## LICENSE INFORMATION

cif2cell is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

cif2cell is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with cif2cell.  If not, see <http://www.gnu.org/licenses/>.



Happy computing!

Torbjorn Bjorkman
COMP, Aaalto University School ofScience and Technology, 
Department of Applied Physics, 
Espoo, Finland
torbjorn@cc.hut.fi


## References
	[1] https://github.com/andeplane/cif2cell-lammps
	[2] https://github.com/kmu/cif2cell


cif2cell-informal
======

	cif2cell support cif file (from VESTA version 3.1.0 - 3.1.7 in 2013-2014. you can get VESTA old version from this URL: http://jp-minerals.org/vesta/archives/)
	Or, expert user recommend to change "_symmetry_" part on cif file. new version cif remove "_symmetry_" and add "alt."


## Install (on ubuntu 18.04 LTS or on Debian 10.0)
1. sudo apt install -y git python python-setuptools python-dev gcc
2. git clone https://github.com/by-student-2017/cif2cell-informal.git
3. cd cif2cell-informal
4. tar zxvf PyCifRW-3.3.tar.gz
5. cd PyCifRW-3.3
6. sudo python setup.py install
7. cd ..
8. sudo python setup.py install


## Install (on google colaboratory)
	!apt update
	!apt install -y git python python-setuptools python-dev gcc
	%cd /content
	!git clone https://github.com/by-student-2017/cif2cell-informal.git
	%cd cif2cell-informal
	!tar zxvf PyCifRW-3.3.tar.gz
	%cd PyCifRW-3.3
	!python2 setup.py install
	%cd /content/cif2cell-informal
	!python2 setup.py install


## Usage (help)


*	cif2cell -h


## Usage (examples)


*	VASP


	cif2cell -p vasp --vasp-format=5 --vasp-encutfac=1.0 --vasp-pseudo-libdr="/home/username/vasp.5.4.1/potpaw_PBE" --vasp-cartesian-lattice-vectors --setup-all  -f *.cif


*	Abinit


	cif2cell -p abinit --abinit-pseudo-JTH-libdr='/home/username/JTH' --setup-all -f *.cif


*	PWscf


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary' --setup-all -f *.cif


*	Elk


	cif2cell -p elk --setup-all -f *.cif


*	Akai-KKR


	cif2cell -p akaikkr --setup-all -f *.cif


*	OpenMX


	cif2cell -p openmx --openmx-pseudo-libdr=/home/username/openmx3.8/DFT_DATA13/PAO --setup-all -f *.cif


*	Lammps


	cif2cell -p lammps --no-reduce -f *.cif


*	DFTB+


	cif2cell -p dftb --no-reduce -f *.cif



*	FDMNES


	cif2cell -p fdmnes --no-reduce -f *.cif


## Usage (Expert mode)


*	PWscf


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary' --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.48 --pwscf-brav --pwscf-spin=no --run-calc --pwscf-run-type=scf -f *.cif


	export PWscf_PAWLIB = $HOME/PSLibrary


	Number of k point * lattice constant = 10 - 13 Angstrom for SCF calculation.  This range is from --k-resolution=2*3.1415/13=0.48 to --k-resolution=2*3.1415/10=0.62.  The default setting is --k-resolution=0.2. This setting is enough to calculate DOS. Hence, nk(DOS) = nk(SCF)*2. nk = number of k point.  This code set automatically [--k-resolution value / 2] for DOS calculation.


*	PWscf (TDDFT)


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary' --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.48 --pwscf-run-type=tddft --pwscf-spin=no --run-calc -f *.cif


	If you would calculate it with gamma point only for B3LYP or HSE, set large value to --k-resolution, e.g. --k-resolution=10 for molecule.


*	PWscf (EELS)


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary' --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.48 --pwscf-run-type=eels --pwscf-spin=no --run-calc -f *.cif


	If you would calculate it with gamma point only for B3LYP or HSE, set large value to --k-resolution, e.g. --k-resolution=10 for molecule.


*	PWscf (PWcond)


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary' --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.48 --pwscf-run-type=lead --run-calc -f *.cif


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary' --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.48 --pwscf-run-type=scat --run-calc -f *.cif


*	PWscf (NEB)


	cif2cell -p pwscf  --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary' --setup-all --k-resolution=0.48 --pwscf-bin-dir=$HOME/q-e-qe-*/bin --pwscf-run-type=neb-start --pwscf-neb-ems=bc2 --pwscf-neb-mu=0.0 --pwscf-neb-bias-voltage=0.5 -f start-structure.cif


	cif2cell -p pwscf  --pwscf-pseudo-PSLibrary-libdr='/home/' --setup-al --pwscf-run-type=neb-end --pwscf-neb-ems=bc2 --pwscf-neb-mu=0.0 --pwscf-neb-bias-voltage=0.5 --run-calc -f end-structure.cif


*	PWscf (relax for interface)


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/student/psl' --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.40 --pwscf-run-type=opt --pwscf-fix-all-pos -f *.cif


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/student/psl' --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.40 --pwscf-run-type=relax --pwscf-fix-atomic-species="Pd" -f *.cif


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/student/psl' --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.40 --pwscf-run-type=opt --pwscf-fix-all-pos --pwscf-move-atomic-species="Pd" -f *.cif


    cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/student/psl' --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.40 --pwscf-run-type=opt --pwscf-fix-all-pos --pwscf-move-atomic-species="Pd" --pwscf-space-group -f *.cif


*	Abinit (Phonon)


	cif2cell -p abinit --abinit-pseudo-JTH-libdr='/home/student/jth' --setup-all --abinit-run-type=phonon --k-resolution=0.4 --abinit-q-resolution=0.8 -f *.cif


*	Akai-KKR (brvtyp mode)


	cif2cell -p akaikkr --setup-all --akaikkr-brvtyp --akaikkr-collect-atoms --akaikkr-run-level=2 -f *.cif


*	FDMNES (TDDFT, L23 edge, Photon energy, Green function, Eimag 0.2, LDA, DOS, Molecular)


	cif2cell -p fdmnes --no-reduce --setup-all --fdmnes-dir=$HOME/fdmnes  --fdmnes-run-type=tddft --fdmnes-edge=L23 --fdmnes-photon-energy --fdmnes-green --fdmnes-eimag=0.2 --fdmnes-dos  --fdmnes-dft-type=LDA  --fdmnes-molecular --run-calc -f *.cif


*	FDMNES (XES, Photon energy, 1.0 core hole, Green function, Eimg 0.2, Crystal structure, auto run)


	cif2cell -p fdmnes --no-reduce --setup-all --fdmnes-dir=$HOME/fdmnes  --fdmnes-run-type=xes --fdmnes-photon-energy --fdmnes-core-hole=1.0  --fdmnes-green --fdmnes-eimag=0.2 --run-calc  -f *.cif


*	FDMNES (XANES, FDM, Crystal structure, 0.5 core hole, DOS for all atom, Crystal structure, auto run)


	cif2cell -p fdmnes --no-reduce --setup-all --fdmnes-dir=$HOME/fdmnes  --fdmnes-photon-energy --fdmnes-core-hole=0.5 --fdmnes-dos-all --run-calc  -f *.cif


*	FDMNES (Photon energy, 1.0 core hole, Green function, Eimag 0.2, DOS for all atom, LDA,  auto run, DAFS, Crystal structure)


	cif2cell -p fdmnes --no-reduce --setup-all --fdmnes-dir=$HOME/fdmnes  --fdmnes-photon-energy --fdmnes-core-hole=1.0 --fdmnes-green --fdmnes-eimag=0.2 --fdmnes-dos-all --fdmnes-dft-type=LDA --run-calc --fdmnes-run-type=dafs -f *.cif


*	FDMNES (Photon energy, 1.0 core hole, Green function, Eimag 0.2, DOS for all atom, LDA,  auto run, RXS, Crystal structure)


	cif2cell -p fdmnes --no-reduce --setup-all --fdmnes-dir=$HOME/fdmnes  --fdmnes-photon-energy  --fdmnes-core-hole=1.0  --fdmnes-green --fdmnes-eimag=0.2 --fdmnes-dos-all --fdmnes-dft-type=LDA --run-calc --fdmnes-run-type=rxs -f *.cif


*	FDMNES (Photon energy, 1.0 core hole, Green function, Eimag 0.2, DOS_comp, LDA,  auto run, Optic, Crystal structure)


	cif2cell -p fdmnes --no-reduce --setup-all --fdmnes-dir=$HOME/fdmnes  --fdmnes-photon-energy  --fdmnes-core-hole=1.0  --fdmnes-green --fdmnes-eimag=0.2 --fdmnes-dos-comp --fdmnes-dft-type=LDA --run-calc --fdmnes-run-type=opitc  -f *.cif


*	FDMNES (Photon energy, 1.0 core hole, Green function, Eimag 0.2, DOS_comp, LDA,  auto run, Circular, Crystal structure)


	cif2cell -p fdmnes --no-reduce --setup-all --fdmnes-dir=$HOME/fdmnes  --fdmnes-photon-energy  --fdmnes-core-hole=1.0  --fdmnes-green --fdmnes-eimag=0.2 --fdmnes-dos-comp --fdmnes-dft-type=LDA --run-calc --fdmnes-run-type=circular  -f *.cif


*	FDMNES (Photon energy, 1.0 core hole, Green function, Eimag 0.2, DOS for all atom, LDA,  auto run, Dichroism, Crystal structure)


	cif2cell -p fdmnes --no-reduce --setup-all --fdmnes-dir=$HOME/fdmnes  --fdmnes-photon-energy  --fdmnes-core-hole=1.0  --fdmnes-green --fdmnes-eimag=0.2 --fdmnes-dos-all --fdmnes-dft-type=LDA --run-calc  --fdmnes-run-type=dichroism -f *.cif


## Usage (examples, Expert mode)


*	VASP


	cif2cell -p vasp --vasp-format=5 --vasp-encutfac=1.0 --vasp-cartesian-lattice-vectors --setup-all  -f *.cif


	export VASP_PAWLIB = $HOME/vasp.5.4.1/potpaw_PBE




*	Abinit


	cif2cell -p abinit --setup-all -f *.cif



	export ABINIT_PAWLIB = $HOME/JTH




*	PWscf


	cif2cell -p pwscf --setup-all -f *.cif


	export PWscf_PAWLIB = $HOME/PSLibrary


## Version: comments


* 1.2.39+PyCifRW-4.4: add FDMNES option.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.38+PyCifRW-4.4: modified k-point settings for PWscf.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.37+PyCifRW-4.4: modified phonon option and k mesh (dos) for Abinit.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.36+PyCifRW-4.4: add --abinit-run-type option for Phonon (trial version).
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.35+PyCifRW-4.4: add --abinit-k-point-even=yes and --elk-k-point-even=yes option (default: yes) for Abinit and Elk.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.34+PyCifRW-4.4: modified DFT+U option.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.33+PyCifRW-4.4: add --pwscf-k-point-even=yes option (default: yes) for PWscf.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.32+PyCifRW-4.4: add NEB option for PWscf.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.31+PyCifRW-4.4: add tddft and eels option for PWscf.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.30+PyCifRW-4.4: add pwcond-lead and pwcond-scat option for PWscf (PWcond calulation).
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.29+PyCifRW-4.4: add cif2cell-lammps code
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.28+PyCifRW-4.4: modified input_dft part for PWscf.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.27+PyCifRW-4.4: add LDA+U option for PWscf, Abinit, Elk, OpenMX and VASP.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.26+PyCifRW-4.4: refine k-point setting for band dispersion on PWscf. (under construction)
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.25+PyCifRW-4.4: add --run-calc option. 
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.24+PyCifRW-4.4: add ibrav option for PWscf. get number of cup (cif2cell.py)
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.23+PyCifRW-4.4: show calculation command. make automatically dos.in file for PWscf.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.22+PyCifRW-4.4: show calculation command. make automatically bands.in and bands.plot file for PWscf.
	Please, input ef value from picene.scf.out file (fermi energy). [the Fermi energy is     X.XXXXX ev] on *.scf.out.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.21+PyCifRW-4.4: add run-type option for pwscf.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.20+PyCifRW-4.4: bug fix. reading pslibrary (H atom) for pwscf.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.19+PyCifRW-4.4: replace PycifRW-3.3 to PyCifRW-4.4.
	could use new verion cif (e.g. new version vesta in 2018-2019)


* 1.2.19: automatically set energy cutoff from pseudo-potential to Abinit and PWscf input files.


* 1.2.18: bug fix making k-point setting.  Add vasp option (run type setting) (test version).
	Set total energy convergence to 1 meV/atom for Elk, Abinit, PWscf, OpenMX and VASP.
	(change part: ESPInterfaces.py)


* 1.2.17: support OpenMX (test version), bug fix Abinit part. (test version) 
	(change part: cif2cell, ESPInterfaces.py, utils.py, elementdata.py)


* 1.2.16: potential file option for Elk, Abinit(ONCVPSP, GBRV, JTH), PWscf(GBRV, PSLibrary(0.3.1 and 1.0.0)) (test version) 
	(change part: cif2cell, ESPInterfaces.py)


* 1.2.15: potential file option for Elk, Abinit(ONCVPSP, GBRV, JTH), PWscf(GBRV, PSLibrary-0.3.1) (test version). 
	(change part: cif2cell, ESPInterfaces.py)


* 1.2.14: modified input file for abinit and elk.  Add run option (scf, dos, band, etc.) (test version).  
	(change part: cif2cell, ESPInterfaces.py)


* 1.2.13: support Akai-KKR option (test version). 
	(change part: cif2cell, ESPInterfaces.py)


* 1.2.12: support Akai-KKR (test version). 
	(change part: cif2cell, ESPInterfaces.py)


* 1.2.11: mofified incarfile on cif2cell and modified VASP(INCAR) and PWscf inputf file on ESPInterfaces.py. 
	(change part: cif2cell, ESPInterfaces.py)


