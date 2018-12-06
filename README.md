### cif2cell-informal


	+This cif2cell support cif file (from VESTA version 3.1.0 - 3.1.7 in 2013-2014 year. you can get VESTA old version from this URL: http://jp-minerals.org/vesta/archives/)


### compile


	+tar zxvf cif2cell-code-1.2.1X


	+sudo python setup.py install




### Usage (help)


	+cif2cell -h




### version: comments


1.2.17: support OpenMX (test version), bug fix Abinit part. (test version) 
	+(change part: cif2cell, ESPInterfaces.py, utils.py, elementdata.py)


1.2.16: potential file option for Elk, Abinit(ONCVPSP, GBRV, JTH), PWscf(GBRV, PSLibrary(0.3.1 and 1.0.0)) (test version) 
	+(change part: cif2cell, ESPInterfaces.py)


1.2.15: potential file option for Elk, Abinit(ONCVPSP, GBRV, JTH), PWscf(GBRV, PSLibrary-0.3.1) (test version). 
	+(change part: cif2cell, ESPInterfaces.py)


1.2.14: modified input file for abinit and elk.  Add run option (scf, dos, band, etc.) (test version).  
	+(change part: cif2cell, ESPInterfaces.py)


1.2.13: support Akai-KKR option (test version). 
	+(change part: cif2cell, ESPInterfaces.py)


1.2.12: support Akai-KKR (test version). 
	+(change part: cif2cell, ESPInterfaces.py)


1.2.11: mofified incarfile on cif2cell and modified VASP(INCAR) and PWscf inputf file on ESPInterfaces.py. 
	+(change part: cif2cell, ESPInterfaces.py)


