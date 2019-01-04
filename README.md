cif2cell-informal
======

	cif2cell support cif file (from VESTA version 3.1.0 - 3.1.7 in 2013-2014. you can get VESTA old version from this URL: http://jp-minerals.org/vesta/archives/)
	Or, expert user recommend to change "_symmetry_" part on cif file. new version cif remove "_symmetry_" and add "alt."


## Compiling


*	tar zxvf cif2cell-code-1.2.35+PyCifRW-4.4.tar.gz


*	cd cif2cell-code-1.2.35+PyCifRW-4.4


*	sudo python setup.py install




## Usage (help)


*	cif2cell -h


## Usage (examples)


*	VASP


	cif2cell -p vasp --vasp-format=5 --vasp-encutfac=1.0 --vasp-pseudo-libdr="/home/username/vasp.5.4.1/potpaw_PBE" --vasp-cartesian-lattice-vectors --setup-all  -f *.cif


	export VASP_PAWLIB = $HOME/vasp.5.4.1/potpaw_PBE




*	Abinit


	cif2cell -p abinit --abinit-pseudo-JTH-libdr='/home/username/JTH' --setup-all -f *.cif



	export ABINIT_PAWLIB = $HOME/JTH




*	PWscf


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary" --setup-all -f *.cif


	export PWscf_PAWLIB = $HOME/PSLibrary




*	Elk


	cif2cell -p elk --setup-all -f *.cif




*	Akai-KKR


	cif2cell -p akaikkr --setup-all -f *.cif




*	OpenMX


	cif2cell -p openmx --openmx-seudo-libdr=/home/username/openmx3.8/DFT_DATA13/PAO --setup-all -f *.cif




*	Lammps


	cif2cell -p lammps --no-reduce -f *.cif


## Usage (Expert mode)


*	PWscf


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary" --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.48 --pwscf-brav --pwscf-spin=no --run-calc --pwscf-run-type=scf -f *.cif


	export PWscf_PAWLIB = $HOME/PSLibrary


	Number of k point * lattice constant = 10 - 13 Angstrom for SCF calculation.  This range is from --k-resolution=2*3.1415/13=0.48 to --k-resolution=2*3.1415/10=0.62.  The default setting is --k-resolution=0.2. This setting is enough to calculate DOS. Hence, nk(DOS) = nk(SCF)*2. nk = number of k point.  This code set automatically [--k-resolution value / 2] for DOS calculation.


*	PWscf (TDDFT)


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary" --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.48 --pwscf-run-type=tddft --pwscf-spin=no --run-calc -f *.cif


	If you would calculate it with gamma point only for B3LYP or HSE, set large value to --k-resolution, e.g. --k-resolution=10 for molecule.


*	PWscf (EELS)


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary" --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.48 --pwscf-run-type=eels --pwscf-spin=no --run-calc -f *.cif


	If you would calculate it with gamma point only for B3LYP or HSE, set large value to --k-resolution, e.g. --k-resolution=10 for molecule.


*	PWscf (PWcond)


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary" --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.48 --pwscf-run-type=lead --run-calc -f *.cif


	cif2cell -p pwscf --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary" --pwscf-bin-dir=$HOME/q-e-qe-6.3/bin --setup-all --k-resolution=0.48 --pwscf-run-type=scat --run-calc -f *.cif


*	PWscf (NEB)


	cif2cell -p pwscf  --pwscf-pseudo-PSLibrary-libdr='/home/username/PSLibrary' --setup-all --k-resolution=0.48 --pwscf-bin-dir=$HOME/q-e-qe-*/bin --pwscf-run-type=neb-start --pwscf-neb-ems=bc2 --pwscf-neb-mu=0.0 --pwscf-neb-bias-voltage=0.5 -f start-structure.cif


	cif2cell -p pwscf  --pwscf-pseudo-PSLibrary-libdr='/home/' --setup-al --pwscf-run-type=neb-end --pwscf-neb-ems=bc2 --pwscf-neb-mu=0.0 --pwscf-neb-bias-voltage=0.5 --run-calc -f end-structure.cif


## Version: comments


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


