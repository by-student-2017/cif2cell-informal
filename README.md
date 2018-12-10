### cif2cell-informal


	Now, this cif2cell support cif file (from VESTA version 3.1.0 - 3.1.7 in 2013-2014. you can get VESTA old version from this URL: http://jp-minerals.org/vesta/archives/)
	Or, expert user recommend to change "_symmetry_" part on cif file. new version cif remove "_symmetry_" and add "alt."


### Compiling


	tar zxvf cif2cell-code-1.2.1X


	sudo python setup.py install




### Usage (help)


	cif2cell -h




### Version: comments


1.2.19+PyCifRW-4.4: replace PycifRW-3.3 to PyCifRW-4.4.
	could use new verion cif (e.g. new version vesta in 2018-2019)


1.2.19: automatically set energy cutoff from pseudo-potential to Abinit and PWscf input files.


1.2.18: bug fix making k-point setting.  Add vasp option (run type setting) (test version).
	Set total energy convergence to 1 meV/atom for Elk, Abinit, PWscf, OpenMX and VASP.
	(change part: ESPInterfaces.py)


1.2.17: support OpenMX (test version), bug fix Abinit part. (test version) 
	(change part: cif2cell, ESPInterfaces.py, utils.py, elementdata.py)


1.2.16: potential file option for Elk, Abinit(ONCVPSP, GBRV, JTH), PWscf(GBRV, PSLibrary(0.3.1 and 1.0.0)) (test version) 
	(change part: cif2cell, ESPInterfaces.py)


1.2.15: potential file option for Elk, Abinit(ONCVPSP, GBRV, JTH), PWscf(GBRV, PSLibrary-0.3.1) (test version). 
	(change part: cif2cell, ESPInterfaces.py)


1.2.14: modified input file for abinit and elk.  Add run option (scf, dos, band, etc.) (test version).  
	(change part: cif2cell, ESPInterfaces.py)


1.2.13: support Akai-KKR option (test version). 
	(change part: cif2cell, ESPInterfaces.py)


1.2.12: support Akai-KKR (test version). 
	(change part: cif2cell, ESPInterfaces.py)


1.2.11: mofified incarfile on cif2cell and modified VASP(INCAR) and PWscf inputf file on ESPInterfaces.py. 
	(change part: cif2cell, ESPInterfaces.py)


