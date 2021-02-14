# Copyright 2010 Torbjorn Bjorkman
# This file is part of cif2cell
#
# cif2cell is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cif2cell is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cif2cell.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************************************************************
#  Description: Interfaces for a number of electronic structure programs. Currently only
#               reads CIF and outputs to the ESP's. Supported programs are: ABINIT, ATAT,
#               CASTEP, CPMD, Crystal09, DFTB+, Elk, EMTO, Exciting, Fleur, Hutsepot, NCOL, 
#               Quantum Espresso, RSPt, Siesta, VASP, xyz, AkaiKKR, OpenMX, Lammps, fdmnes
#               
#  Author:      Torbjorn Bjorkman, torbjorn.bjorkman(at)abo.fi
#  Affiliation: Abo Akademi University
#               Physics/Department of Natural Sciences
#               Porthansgatan 3
#               20500 Turku, Finland
#  ORCID:       0000-0002-1154-9846
#******************************************************************************************
import copy
import os
import math
import string
from utils import *
from elementdata import *
from uctools import *
from re import search


################################################################################################
ed = ElementData()
suspiciouslist = set(["Cr", "Mn", "Fe", "Co", "Ni",
        "Ce","Pr","Nd","Pm","Sm","Eu",
        "Gd","Tb","Dy","Ho","Er","Tm",
        "Th","Pa","U","Np","Pu"])
initialmoments = {"Cr" : 3, "Mn" : 3, "Fe" : 3, "Co" : 3, "Ni" : 1,
        "Ce" : 1, "Pr" : 2, "Nd" : 3, "Pm" : 4, "Sm" : 5, "Eu" : 6,
        "Gd" : 7, "Tb" : 8, "Dy" : 9, "Ho" : 10, "Er" : 11, "Tm" : 12,
        "Th" : 1, "Pa" : 2, "U" : 3, "Np" : 4, "Pu" : 5 }

################################################################################################
class GeometryOutputFile:
    """
    Parent class for electronic struture code files generated from geometrical information.
    A CrystalStructure object and a documentation string are required input.
    """
    def __init__(self, crystalstructure, string):
        self.cell = crystalstructure
        self.docstring = string

################################################################################################
# ATAT FILE
class ATATFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting an ATAT input file
    and the method __str__ that outputs the contents of the file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Document string on first line after '//'
        self.programdoc = string.rstrip("\n")
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                tmp.add(b.spcstring())
        self.species = list(tmp)
        # make sure the docstring goes on one line
        self.cell.newunit("bohr")
    def __str__(self):
        filestring = str(self.cell.a)+" "+str(self.cell.b)+" "+str(self.cell.c)+" "+str(self.cell.alpha)+" "+str(self.cell.beta)+" "+str(self.cell.gamma)+"\n"
        filestring += str(self.cell.lattrans)
        for a in self.cell.atomdata:
            for b in a:
                filestring += str(b.position)
                if b.alloy():
                    for k,v in b.species.iteritems():
                        filestring += k+"="+str(v)+","
                    filestring = filestring.rstrip(",")+"\n"
                else:
                    filestring += " "+b.spcstring(separator=',')+"\n"
        return filestring
        
################################################################################################
# HUTSEPOT FILE
class HUTSEPOTFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting a hutsepot input file
    and the method __str__ that outputs the contents of the file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Document string on first line after '//'
        self.programdoc = string.rstrip("\n")
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                tmp.add(b.spcstring())
        self.species = list(tmp)
        # make sure the docstring goes on one line
        self.cell.newunit("bohr")
    def __str__(self):
        streck="-------------------------------------------------------------------------------"
        filestring = streck+"\n"
        filestring += "------------------------- Generated by cif2cell -------------------------------\n"
        filestring = streck+"\n"
        t = self.cell.lengthscale
        filestring += "3D unit cell alat=%18.12f blat=%18.12f clat=%18.12f\n"%(t,t,t)
        filestring += "   rb="+str(self.cell.latticevectors[0].scalmult(t))+"ascale=1.0\n"
        filestring += "      "+str(self.cell.latticevectors[1].scalmult(t))+"bscale=1.0\n"
        filestring += "      "+str(self.cell.latticevectors[2].scalmult(t))+"cscale=1.0\n"
        # positions
        filestring += "-------------------------------------------------------------------------------\n"
        filestring += "------------------------------- atomic positions ------------------------------\n"
        filestring += "-------------------------------------------------------------------------------\n"
        positionstring = ""
        species = 0
        atom = 0
        for sp in self.species:
            nr = 0
            species += 1
            for a in self.cell.atomdata:
                for b in a:
                    if b.spcstring() == sp:
                        atom += 1
                        nr += 1
                        p = Vector(mvmult3(self.cell.latticevectors,b.position.scalmult(self.cell.lengthscale)))
                        positionstring += str(species)+"."+b.spcstring()+"_"+str(nr)+" type=%i"%(atom)+" nat=60 tau="+str(p)
                        positionstring += "\n"
        filestring += positionstring
        filestring += "-------------------------------------------------------------------------------\n"
        filestring += "---------------------------- atomic configurations ----------------------------\n"
        filestring += "-------------------------------------------------------------------------------\n"
        species = 0
        for sp in self.species:
            species += 1
            filestring += str(species)+". "+ed.hutsepotelements[sp]+"\n"
        filestring += "-------------------------------------------------------------------------------\n"
        filestring += "------------------------------- atomic options --------------------------------\n"
        filestring += "-------------------------------------------------------------------------------\n"
        species = 0
        atom = 0
        for sp in self.species:
            species += 1
            for a in self.cell.atomdata:
                for b in a:
                    if b.spcstring() == sp:
                        atom += 1
                        filestring += str(species)+". atom="+sp+" type="+str(atom)
                        filestring += " fix=F lmax=3 lmaxv=0 conc=1.0 mtz=T sort="+str(atom)+"\n"
        filestring += "-------------------------------------------------------------------------------\n"
        filestring += "--------------------------------- potentials ----------------------------------\n"
        filestring += "-------------------------------------------------------------------------------\n"
        species = 0
        atom = 0
        for sp in self.species:
            species += 1
            nr = 0
            for a in self.cell.atomdata:
                for b in a:
                    if b.spcstring() == sp:
                        atom += 1
                        nr += 1
                        filestring += str(species)+". type="+str(atom)
                        filestring += " np=1001 r1=1.0E-05 rnp=-2 pfile="+sp+str(nr)+".pot\n"
        filestring += "-------------------------------------------------------------------------------\n"
        return filestring

################################################################################################
# ASE FILE
class ASEFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting data for ASE
    and the method __str__ that outputs the contents of the file as a string of
    python code.
    """
    def __init__(self,crystalstructure,docstring):
        GeometryOutputFile.__init__(self,crystalstructure,docstring)
        # Variables
        self.cartesian = True  # Cartesian coordinates?
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                tmp.add(b.spcstring())
        self.species = list(tmp)
    def __str__(self):
        filestring = "from ase import *\n\n"
        # Cartesian or lattice coordinates?
        if self.cartesian:
            transmtx = []
            for i in range(3):
                transmtx.append([])
                for j in range(3):
                    transmtx[i].append(self.cell.latticevectors[i][j] * self.cell.lengthscale)
        else:
            transmtx = [[1,0,0],[0,1,0],[0,0,1]]
        # positions and number of species
        nspcs = []
        positionstring = ""
        for sp in self.species:
            nsp = 0
            for a in self.cell.atomdata:
                for b in a:
                    if b.spcstring() == sp:
                        nsp += 1
                        p = Vector(mvmult3(transmtx,b.position))
                        positionstring += "%11f, %11f, %11f),\n               ("%(p[0],p[1],p[2])
            nspcs.append(nsp)
        positionstring = positionstring.rstrip("\n (,")+"],\n"

        # Atoms object
        filestring += "atoms = Atoms("
        # Species
        for i in range(len(self.species)):
            filestring += "['"+self.species[i]+"' for i in range("+str(nspcs[i])+")]+"
        filestring = filestring.rstrip("+")+",\n"
        # Positions
        filestring += "              [("
        filestring += positionstring
        # Boundary conditions
        filestring += "              pbc = (True,True,True))\n"
        # Set lattice vectors
        filestring += "atoms.set_cell([["
        for i in range(3):
            for j in range(3):
                filestring += "%12f, "%(self.cell.latticevectors[i][j]*self.cell.lengthscale)
            filestring = filestring.rstrip(", ")+"],\n                ["
        filestring = filestring.rstrip(",[ ]\n")+"]],\n"
        filestring += "                scale_atoms = %s)\n"%str(not self.cartesian)
        
        return filestring
    
################################################################################################
# CFG FILE
class CFGFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting a .cfg file
    and the method __str__ that outputs the contents of the .coo file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
    def __str__(self):
        # Set up atom list for printing.
        tmplist = list(self.cell.atomset)
        atomlist = []
        for a in tmplist:
            if a.alloy():
                for sp,occ in a.species.iteritems():
                    atomlist.append(AtomSite(position=a.position,species={sp : occ},charges={sp : a.charges[sp]}))
            else:
                atomlist.append(a)
        atomlist.sort(key = lambda x: max([ed.elementnr[sp] for sp in x.species]),reverse=True)
        prevsp = ""
        natoms = len(atomlist)
        # Make string
        filestring = self.docstring
        filestring += "Number of particles = %i \n"%(natoms)
        filestring += "A = 1.0 Angstrom\n"
        filestring += "H0(1,1) = %f A\n"%(self.cell.lengthscale*self.cell.latticevectors[0][0])
        filestring += "H0(1,2) = %f A\n"%(self.cell.lengthscale*self.cell.latticevectors[0][1])
        filestring += "H0(1,3) = %f A\n"%(self.cell.lengthscale*self.cell.latticevectors[0][2])
        filestring += "H0(2,1) = %f A\n"%(self.cell.lengthscale*self.cell.latticevectors[1][0])
        filestring += "H0(2,2) = %f A\n"%(self.cell.lengthscale*self.cell.latticevectors[1][1])
        filestring += "H0(2,3) = %f A\n"%(self.cell.lengthscale*self.cell.latticevectors[1][2])
        filestring += "H0(3,1) = %f A\n"%(self.cell.lengthscale*self.cell.latticevectors[2][0])
        filestring += "H0(3,2) = %f A\n"%(self.cell.lengthscale*self.cell.latticevectors[2][1])
        filestring += "H0(3,3) = %f A\n"%(self.cell.lengthscale*self.cell.latticevectors[2][2])
        filestring += ".NO_VELOCITY.\n"
        ## # Cut the fancy stuff for now, stick with just the positions
        ## filestring += "entry_count = 3\n"
        filestring += "entry_count = 6\n"
        for a in atomlist:
                for sp,occ in a.species.iteritems():
                        if prevsp != sp:
                            filestring += "%i\n"%(int(round(ed.elementweight[sp])))
                            filestring += sp+"\n"
                        prevsp = sp
                        DW = 0.45*ed.elementnr['Si']/ed.elementnr[sp] # Debye-Waller factor, QSTEM prescription
                        filestring += str(a.position)+" %f "%(DW)+" %f "%(occ)+" %f\n"%(a.charges[sp])
                        ## filestring += str(a.position)+"\n"
        return filestring
    
################################################################################################
# COO FILE
class COOFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting a .coo file
    and the method __str__ that outputs the contents of the .coo file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Document string on first line after '//'
        self.programdoc = string.rstrip("\n")
    def __str__(self):
        filestring = "//"+self.programdoc+"\n"
        a = self.cell.latticevectors[0].length()*self.cell.lengthscale
        b = self.cell.latticevectors[1].length()*self.cell.lengthscale
        c = self.cell.latticevectors[2].length()*self.cell.lengthscale
        alpha = abs(self.cell.latticevectors[1].angle(self.cell.latticevectors[2]))*180/pi
        beta = abs(self.cell.latticevectors[2].angle(self.cell.latticevectors[0]))*180/pi
        gamma = abs(self.cell.latticevectors[0].angle(self.cell.latticevectors[1]))*180/pi
        filestring += " %10.7f %10.7f %10.7f"%(a,b,c)
        filestring += " %10.7f %10.7f %10.7f %i\n"%(alpha,beta,gamma,len(self.cell.atomset))
        for a in self.cell.atomdata:
            for b in a:
                filestring += str(b.position)+" %3i  0.500 0.000 1.000\n"%(ed.elementnr[b.spcstring()])
        return filestring
    
################################################################################################
# XYZ FILE
class XYZFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting an .xyz file
    and the method __str__ that outputs the contents of the .xyz file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # To be put on the second line
        self.programdoc = ""
    def __str__(self):
        filestring = ""
        filestring += "%i \n"%sum([len(v) for v in self.cell.atomdata])
        filestring += self.docstring+"\n"
        lv = []
        for i in range(3):
            lv.append([])
            for j in range(3):
                lv[i].append(self.cell.lengthscale*self.cell.latticevectors[i][j])
        for a in self.cell.atomdata:
            for b in a:
                t = Vector(mvmult3(lv,b.position))
                filestring += str(b).split()[0]+"  "+str(t)+"\n"
        return filestring
    
################################################################################################
# NCOL FILES
class OldNCOLFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the ncol program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.jobnam = "default"
        self.bstrjobnam = "default"
        # To be put on the first line
        self.programdoc = ""
        # Set atomic units for length scale
        self.cell.newunit("bohr")
    def __str__(self):
        # Element data
        ed = ElementData()
        # l quantum number setup (same as from bstr)
        l = { "s" : 2, "p" : 2, "d" : 3, "f" : 4 }
        filestring = ""
        tmpstring = "BULK      IDSYST=  7 SCRATCH=R"
        tmpstring = tmpstring.ljust(25)+"    "+deletenewline(self.programdoc,replace=" ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(10)+" MSGL.=  1 BSAVE..=N COLD...=Y DOS...=N SPO...=N ISM...=G RCLCR...=Y\n"
        filestring += tmpstring
        filestring += "FOR001=./"+self.bstrjobnam+".tfm\n"
        filestring += "FOR002=\n"
        filestring += "FOR003=\n"
        filestring += "FOR004=\n"
        filestring += "FOR006=\n"
        filestring += "FOR010=\n"
        filestring += "Band: 4 lines, "+deletenewline(self.docstring,replace=" ")+"\n"
        filestring += "NITER.=200 NOB..=  2 NPRN.=  0 NFIX.=  0 MIXKEY=  2 NCOL.=Y  PMODE=K\n"
        filestring += "REP.....=B FIXD...=Y CRT....=S NB...= 16 CLSIZE= 32 NPROW= 0 NPCOL= 0\n"
        filestring += "NKX...=  1 NKY..=  1 NKZ..=  1 TFERMI..= 2000.0(K)\n"
        filestring += "AMIX.....=     0.100 TOLE....= 0.0000100 TOLEL...= 0.0000010\n"
        # average wigner-seitz radius
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        volume = abs(det3(self.cell.latticevectors))
        wsr = self.cell.lengthscale * (3*volume/(nosites * 4 * pi))**third
        filestring += "SWS......= %9f NP...=  1 SMIX.= 0.500 TMIX.= 0.0000\n" % wsr
        filestring += "Setup: 3 + NQ*NS*NP lines\n"
        filestring += "EFGS.....=    0.0000 EFGS....=   0.00000 FTMAG...=  0.000000\n"
        filestring += "DEO(l)...=     0.020     0.010     0.005     0.001      0.02\n"
        filestring += "Symb IQ IT NL IP NSP   SWP  QTRO  SPLT NFIX NDWF     Eny(spdf)\n"
        # set first species
        if self.cell.atomdata[0][0].alloy():
            prevspecies = "??"
        else:
            for v in self.cell.atomdata[0][0].species:
                prevspecies = v
        # type loop
        iq = 1
        it = 1
        nsp = 1
        for a in self.cell.atomdata:
            for b in a:
                if b.alloy():
                    species = "??"
                else:
                    species = b.spcstring()
                if species != prevspecies:
                    prevspecies = species
                    nsp += 1
                tmpstring = species.ljust(2)+"    "+str(iq).ljust(3)+str(it).ljust(3)
                try:
                    tmpstring += str(l[ed.elementblock[species]]).ljust(3)+str(1).ljust(3)
                except KeyError:
                    tmpstring += "  ?  1"
                tmpstring += str(nsp).ljust(3)
                tmpstring += "  1.000 .000 0.00 0000 1111   .0   .0   .0   .0"
                if b.alloy():
                    # print alloy components at the end of the line
                    tmpstring += "       "+b.spcstring()
                filestring += tmpstring+"\n"
                iq += 1
            it += 1
        for a in self.cell.atomdata:
            filestring += "Theta....=     90.00 Phia....=      0.00 FIXMOM..=         N moment..=      0.0\n"
        filestring += "PQX......=      0.00 PQY.....=      0.00 PQZ.....=   0.00000 COORD...=L\n"
        filestring += "Atom: 4 lines + NT*6 lines\n"
        filestring += "IEX...=  4  NP..=500 NES..= 15 NITER=250 IWAT.=  0\n"
        filestring += "VMIX.....=  0.300000 RWAT....=  3.500000 RMAX....= 20.000000\n"
        filestring += "DPAS.....=  0.049000 DR1.....=  1.00E-08 TEST....=  1.00E-08\n"
        filestring += "TESTE....=  1.00E-07 TESTY...=  1.00E-08 TESTV...=  1.00E-07\n"
        for a in self.cell.atomdata:
            for comp in a[0].species:
                filestring += comp+"\n"
                try:
                    filestring += ed.emtoelements[comp]
                except KeyError:
                    filestring += "\n\n\n\n\n"
        return filestring

class BSTRFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the bstr program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.jobnam = "default"
        self.a = 1
        self.b = 1
        self.c = 1
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        ed = ElementData()
        filestring = ""
        tmpstring = "BSTR      IDSYST=  7"
        tmpstring = tmpstring.ljust(40)+deletenewline(self.programdoc,replace=" ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(9)+" MSGL.=   1 \n"
        filestring += tmpstring
        filestring += "MODE....=B STORE..=Y SCREEN.=B CMBC...=Y\n"
        filestring += "FOR001=\n"
        filestring += "FOR006=\n"
        filestring += deletenewline(self.docstring,replace=" ")+"\n"
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        # Setting the real space summation cutoff to 4.5*(wigner-seitz radius)
        volume = abs(det3(self.cell.latticevectors))
        wsr = (3*volume/(nosites * 4 * pi))**third
        tmpstring = "IALF...= 0 NPRN..= 1 DMAX....=%10.5f \n" % (wsr*4.5)
        filestring += tmpstring
        filestring += "ALF(spdf)= 0.3205350 0.0413320 0.0084290 0.0015370\nDKAPPA...= 0.00010\n"
        tmpstring = "NQ3....=%3i LAT...= 0 IPRIM.= 0" % nosites
        filestring += tmpstring
        # Set up basis functions. Just setting lmax = 2 for s-/p-, 3 for d- and 4 for f- blocks
        tmpstring = "\nNLX(IQ)..="
        for a in self.cell.atomdata:
            for b in a:
                for k in b.species:
                    l = 1
                    if ed.elementblock[k] == "s" or ed.elementblock[k] == "p":
                        l = max(l,2)
                    elif ed.elementblock[k] == "d":
                        l = max(l,3)
                    elif ed.elementblock[k] == "f":
                        l = max(l,4)
                tmpstring += " %1i" % l
                if len(tmpstring) % 69 == 0:
                    tmpstring += "\n          "
        # Need to strip newline character if the last line was 69 characters long...
        tmpstring = tmpstring.rstrip(string.whitespace)
        tmpstring = tmpstring+"\n"
        filestring += tmpstring
        # Print lattice vectors
        coa = self.c / self.a
        boa = self.b / self.a
        filestring += "A........=  1.00000000 B.......=  1.00000000 C.......=  1.00000000\n"
        tmpstring = ""
        lv = self.cell.latticevectors
        for i in range(3):
            tmpstring += "BSX......=%12.7f BSY.....=%12.7f BSZ.....=%12.7f\n" % (lv[i][0],lv[i][1],lv[i][2])
        filestring += tmpstring
        # All positions
        it = 1
        for a in self.cell.atomdata:
            for b in a:
                pos = mvmult3(lv,b.position)
                tmpstring = "QX.......=%12.7f QY......=%12.7f QZ......=%12.7f" % (pos[0],pos[1],pos[2])
                tmpstring += "      "+b.spcstring()+"\n"
                filestring += tmpstring
            it += 1
        filestring += "LAMDA....=    2.5000 AMAX....=    5.5000 BMAX....=    5.5000\n"
        return filestring
    
################################################################################################
# RSPT FILES
class CellgenFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a cellgen.inp file and the method
    __str__ that outputs the contents of an cellgen.inp file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.supercellmap = [[1,0,0],[0,1,0],[0,0,1]]
        self.referencevector = [0,0,0]
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
    def __str__(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        # Add lattice constant
        filestring += "# Lattice constant in a.u.: "+str(self.cell.lengthscale)+"\n"
        # RSPt reads the lattice vectors as columns...
        filestring +="# Lattice vectors (columns)\n"
        tmpstring = ""
        for i in range(3):
            for j in range(3):
                tmpstring += "%19.15f "%self.cell.latticevectors[j][i]
            tmpstring += "\n"
        filestring += tmpstring
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        filestring += "# Sites\n"
        filestring += str(nosites)+"\n"
        it = 1
        for a in self.cell.atomdata:
            for b in a:
                tmpstring = ""
                tmpstring += str(b.position)+" "
                if b.alloy():
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    tmpstring += "%3i"%ed.elementnr[b.spcstring()]
                tmpstring += " l "+chr(it+96)+"   # "+b.spcstring()+"\n"
                filestring += tmpstring
            it += 1
        filestring += "# Supercell map\n"
        tmpstring = ""
        for i in self.supercellmap:
            for j in i:
                tmpstring += str(j).rjust(4)
            tmpstring += "\n"
        filestring += tmpstring
        filestring += "# Reference vector\n"
        tmpstring = ""
        for i in self.referencevector:
            tmpstring += "%19.15f "%i
        tmpstring += "\n"
        filestring += tmpstring
        return filestring

################################################################################################
class SymtFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an old format symt.inp file and the method
    __str__ that outputs the contents of an symt.inp file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        # Default spin axis is [0,0,0]
        self.spinaxis = [0.0, 0.0, 0.0]
        self.rsptcartlatvects = False
        self.passwyckoff = False
        self.printlabels = False
    def __str__(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        # Add lattice constant
        filestring += "# Lattice constant in a.u.: "+str(self.cell.lengthscale)+"\n"
        # RSPt reads the lattice vectors as columns...
        filestring +="# Lattice vectors (columns)\n"
        if self.rsptcartlatvects:
            fac = self.cell.lengthscale
        else:
            fac = 1.0
        tmpstring = ""
        for i in range(3):
            for j in range(3):
                tmpstring += "%19.15f "%(self.cell.latticevectors[j][i]*fac)
            tmpstring += "\n"
        filestring += tmpstring
        filestring += "# Spin axis\n"
        filestring += "%19.15f %19.15f %19.15f  l\n"%(self.spinaxis[0],self.spinaxis[1],self.spinaxis[2])
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        filestring += "# Sites\n"
        filestring += str(nosites)+"\n"
        it = 1
        label = "a"
        for a in self.cell.atomdata:
            for b in a:
                tmpstring = ""
                tmpstring += str(b.position)+" "
                if b.alloy():
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    tmpstring +=  "%3i"%ed.elementnr[b.spcstring()]
                if self.passwyckoff:
                    label = chr(it+96)
                if self.printlabels:
                    tmpstring += " l "+label+"   # "+b.label+"\n"
                else:
                    tmpstring += " l "+label+"   # "+b.spcstring()+"\n"
                filestring += tmpstring
            it += 1
        return filestring

################################################################################################
class SymtFile2(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a new format symt.inp file and the method
    __str__ that outputs the contents of an symt.inp file as a string.
    """
    def __init__(self,crystalstructure,string,kresolution=0.1):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        # Default spin axis is [0,0,0]
        self.spinaxis = Vector([0.0, 0.0, 0.0])
        # parameters for spin polarization
        self.spinpol = False
        self.relativistic = False
        self.forcenospin = False
        self.rsptcartlatvects = False
        self.mtradii = 0
        self.passwyckoff = False
        # k-mesh generation etc.
        self.setupall = False
        self.kresolution = kresolution
        self.nokshifts = False
        self.kshifts = [[0,0,0],[1,1,1]]
        self.printlabels = False
    def __str__(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        # Add lattice constant
        filestring += "# Lattice constant in a.u.\n"
        filestring += "lengthscale\n"
        if self.rsptcartlatvects:
            filestring += "1.000 \n"
        else:
            filestring += str(self.cell.lengthscale)+"\n"
        if self.spinpol and not self.forcenospin:
            filestring += "# Spin polarized calculation\nspinpol\n"
            filestring += "# Spin polarize atomic densities\nspinpol_atomdens\n"
        if self.relativistic:
            if self.setupall:
                filestring += "# Relativistic symmetries\nspinorbit\n"
            else:
                filestring += "# Relativistic symmetries\nfullrel\n"
            # Default to z-direction for relativistic calculations...
            t = self.spinaxis - Vector([0.,0.,0.])
            if t.length() < 1e-7:
                self.spinaxis = mvmult3(minv3(self.cell.latticevectors),Vector([0.,0.,1.]))
                # ... unless these space group settings, when we pick a more likely high-symmetry axis
                if self.cell.spacegroupsetting == "A":
                    self.spinaxis = mvmult3(minv3(self.cell.latticevectors),Vector([1.,0.,0.]))
                elif self.cell.spacegroupsetting == "B":
                    self.spinaxis = mvmult3(minv3(self.cell.latticevectors),Vector([0.,1.,0.]))
        if self.mtradii != 0:
            filestring += "# Choice of MT radii\n"
            filestring += "mtradii\n"+str(self.mtradii)+"\n"
        # RSPt reads the lattice vectors as columns...
        filestring += "# Lattice vectors (columns)\n"
        filestring += "latticevectors\n"
        tmpstring = ""
        if self.rsptcartlatvects:
            fac = self.cell.lengthscale
        else:
            fac = 1.0
        for i in range(3):
            for j in range(3):
                tmpstring += "%19.15f "%(self.cell.latticevectors[j][i]*fac)
            tmpstring += "\n"
        filestring += tmpstring            
        filestring += "# Spin axis\n"
        filestring += "spinaxis\n"
        filestring += "%19.15f %19.15f %19.15f  l\n"%(self.spinaxis[0],self.spinaxis[1],self.spinaxis[2])
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        filestring += "# Sites\n"
        filestring += "atoms\n"
        filestring += str(nosites)+"\n"
        it = 1
        label = "a"
        for a in self.cell.atomdata:
            for b in a:
                label = "a"
                tmpstring = ""
                tmpstring += str(b.position)+" "
                if b.alloy():
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    tmpstring += "%3i"%ed.elementnr[b.spcstring()]
                if self.passwyckoff:
                    label = chr(it+96)
                if self.setupall and b.spcstring() in suspiciouslist and not self.forcenospin:
                    label = "up"
                if self.printlabels:
                    tmpstring += " l "+label+"   # "+b.label+"\n"
                else:
                    tmpstring += " l "+label+"   # "+b.spcstring()+"\n"
                filestring += tmpstring
            it += 1
        # k-mesh setup for new input
        if self.setupall:
            # Using k-resolution together with Froyen map needs supervised choice of mesh,
            # or they easily become unnecessarily dense, so don't use this feature.
            ## filestring += "\n"
            ## filestring += "# k space resolution\n"
            ## filestring += "kresolution\n"
            ## filestring += "  %f\n"%(self.kresolution)
            # Guess a suitable Froyen map !!! Column vectors for RSPt !!!
            mapmatrix = LatticeMatrix([[1,0,0],[0,1,0],[0,0,1]])
            if self.cell.primcell:
                if self.cell.spacegroupsetting == 'F':
                    mapmatrix = LatticeMatrix([[1,1,0],[1,0,1],[0,1,1]])
                elif self.cell.spacegroupsetting == 'I':
                    if self.cell.crystal_system() == 'cubic':
                        mapmatrix = LatticeMatrix([[-1,1,1],[1,-1,1],[1,1,-1]])
                    else:
                        mapmatrix = LatticeMatrix([[2,0,1],[0,2,1],[0,0,2]])
                elif self.cell.spacegroupsetting == 'A':
                    mapmatrix = LatticeMatrix([[1,0,0],[0,1,-1],[0,1,1]])
                elif self.cell.spacegroupsetting == 'B':
                    mapmatrix = LatticeMatrix([[1,0,-1],[0,1,0],[1,0,1]])
                elif self.cell.spacegroupsetting == 'C':
                    mapmatrix = LatticeMatrix([[1,-1,0],[1,1,0],[0,0,1]])
                elif self.cell.spacegroupsetting == 'R' and abs(self.cell.latticevectors[0].angle(self.cell.latticevectors[1])*180/pi) > 10:
                    # Generate in hexagonal supercell unless the rhombohedral angle is close to 90 degrees.
                    mapmatrix = LatticeMatrix([[1,0,1],[-1,1,1],[0,-1,1]])
            # Determine mesh
            reclatvect = LatticeMatrix(mmmult3(self.cell.reciprocal_latticevectors().transpose(),mapmatrix)).transpose()
            for j in range(3):
                for i in range(3):
                    reclatvect[j][i] = reclatvect[j][i] / self.cell.lengthscale
            # Lengths of reciprocal lattice vectors
            reclatvectlen = [elem.length() for elem in reclatvect]
            kgrid = [max(1,int(round(elem/self.kresolution))) for elem in reclatvectlen]
            # Manual adjustments to make the choice work well with the Froyen mesh.
            # Some centerings should have even meshes, for rhombohedral it should be dividable by 3
            # along c.
            if self.cell.primcell:
                if self.cell.spacegroupsetting == 'F' or self.cell.spacegroupsetting == "I":
                    for i in range(3):
                        kgrid[i] += kgrid[i]%2
                elif self.cell.spacegroupsetting == 'A':
                    kgrid[1] += kgrid[1]%2
                    kgrid[2] += kgrid[2]%2
                elif self.cell.spacegroupsetting == 'B':
                    kgrid[0] += kgrid[0]%2
                    kgrid[2] += kgrid[2]%2
                elif self.cell.spacegroupsetting == 'C':
                    kgrid[0] += kgrid[0]%2
                    kgrid[1] += kgrid[1]%2
                elif self.cell.spacegroupsetting == 'R' and abs(self.cell.latticevectors[0].angle(self.cell.latticevectors[1])*180/pi) > 10:
                    for i in range(3):
                        # This rounds to nearest multiple of 3
                        if kgrid[i]%3 == 1:
                            kgrid[i] -= 1
                        elif kgrid[i]%3 == 2:
                            kgrid[i] += 1
            filestring += "\n# k-points\n"
            filestring += "kpoints\n"
            filestring += " %i %i %i\n\n"%(kgrid[0], kgrid[1], kgrid[2])
            # Write Froyen map.
            filestring += "# Froyen map\n"
            filestring += "kmapmatrix\n"
            for v in mapmatrix:
                filestring += " %4i %4i %4i\n"%(v[0],v[1],v[2])
        # Return
        return filestring

################################################################################################
class Crystal09File(GeometryOutputFile):
    """
    Class for storing the geometrical data needed by Crystal09 and the method
    __str__ that outputs the contents of an Crystal09 input file as a string.
    Presently only handles standard settings (space group numbers, not H-M symbols),
    and the special case of rhombohedral settings for relevant trigonal space groups.
    """
    def __init__(self,crystalstructure,string,rhombohedral=False):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("angstrom")
        # Rhombohedral cell setting
        self.rhombohedral = rhombohedral
        # Make sure the docstring has the form of a f90 comment
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("!")
            string = "!"+string+"\n"
            self.docstring += string
    def __str__(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        filestring += "CRYSTAL\n"
        # Space group setting and crystal parameters
        if self.cell.is_spacegroup("triclinic"):
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f %13.8f %13.8f %13.8f %13.8f\n"%(self.cell.a, self.cell.b, self.cell.c, self.cell.alpha, self.cell.beta, self.cell.gamma)
        elif self.cell.is_spacegroup("monoclinic"):
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f %13.8f %13.8f\n"%(self.cell.a, self.cell.b, self.cell.c, self.cell.beta)
        elif self.cell.is_spacegroup("orthorhombic"):
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f %13.8f\n"%(self.cell.a, self.cell.b, self.cell.c)
        elif self.cell.is_spacegroup("tetragonal"):
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f\n"%(self.cell.a, self.cell.c)
        elif self.cell.is_spacegroup("trigonal") and not (self.cell.is_spacegroup("rhombohedral") and self.rhombohedral):
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f\n"%(self.cell.a, self.cell.c)
        elif self.cell.is_spacegroup("rhombohedral") and self.rhombohedral:
                filestring += "0 1 0\n"
                filestring += str(self.spacegroupnr)+"\n"
                filestring += "%13.8f %13.8f\n"%(self.cell.latticevectors[0].length()*self.cell.lengthscale, self.cell.latticevectors[1].angle(self.cell.latticevectors[2])*180/pi)
        elif self.cell.is_spacegroup("hexagonal"):
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f\n"%(self.cell.a, self.cell.c)
        elif self.cell.is_spacegroup("cubic"):
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f\n"%(self.cell.a)
        else:
            if self.force:
                sys.stderr.write("***Warning: Could not determine crystal system corresponding to space group "+str(self.spacegroupnr)+".")
                filestring += "0 0 0\n"
                filestring += str(self.spacegroupnr)+"\n"
                filestring += "%13.8f %13.8f %13.8f %13.8f %13.8f %13.8f\n"%(self.cell.a, self.cell.b, self.cell.c, self.cell.alpha, self.cell.beta, self.cell.gamma)
            else:
                return "***Error: Could not determine crystal system corresponding to space group "+str(self.spacegroupnr)+"."
        # Number of atoms
        filestring += str(len(self.cell.ineqsites))+"\n"
        # Atomic numbers and representative positions
        for a in self.cell.atomdata:
            if len(a[0].species) > 1:
                # don't know what to put for an alloy
                filestring += "??"
            else:
                for k in a[0].species:
                    filestring += str(ed.elementnr[k]).rjust(2)
            filestring += "  "+str(a[0].position)+"      ! "+a[0].spcstring()+"\n"
        filestring += "END\n"
        return filestring

################################################################################################
class SpacegroupFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a spacegroup.in file and the method
    __str__ that outputs the contents of an spacegroup.in file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.HermannMauguin = ""
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.supercelldims = [1, 1, 1]
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has the form of a f90 comment
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("!")
            string = "!"+string+"\n"
            self.docstring += string
    def __str__(self):
        filestring = ""
        if (self.HermannMauguin[-1] == 'R' or self.HermannMauguin[-1] == 'H') and self.HermannMauguin[-2] != ':':
            self.HermannMauguin = self.HermannMauguin[0:len(self.HermannMauguin)-1]+':'+self.HermannMauguin[-1]
        tmpstring=" '"+self.HermannMauguin+"'"
        tmpstring = tmpstring.ljust(50)+": hrmg\n"
        filestring += tmpstring
        tmpstring = ""
        tmpstring += " %15.11f" % (self.a)
        tmpstring += " %15.11f" % (self.b)
        tmpstring += " %15.11f" % (self.c)
        tmpstring = tmpstring.ljust(50)+": a, b, c\n"
        filestring += tmpstring
        tmpstring = " %15.9f %15.9f %15.9f"%(self.gamma,self.beta,self.alpha)
        tmpstring = tmpstring.ljust(50)+": ab, ac, bc\n"
        filestring += tmpstring
        tmpstring = ""
        for i in self.supercelldims:
            tmpstring += str(i)+"  "
        tmpstring = tmpstring.ljust(50)
        tmpstring += ": ncell\n"
        filestring += tmpstring
        filestring += ".true.".ljust(50)+": primcell\n"
        # Get species info
        species = set([])
        for occ in self.cell.occupations:
            spcstring = ""
            for k in occ:
                spcstring += k+"/"
            spcstring = spcstring.rstrip("/")
            species.add(spcstring)
        tmpstring = str(len(species)).ljust(50)+": nspecies\n"
        filestring += tmpstring
        for spcs in species:
            # find number of representative sites for this species
            spcsites = 0
            positionstring = ""
            i = 0
            for occ in self.cell.occupations:
                spcstring = ""
                for k in occ:
                    spcstring += k+"/"
                spcstring = spcstring.rstrip("/")
                if spcstring == spcs:
                    spcsites += 1
                    positionstring += str(self.cell.ineqsites[i])+"\n"
                i += 1
            # output species info
            if len(spcs) > 2:
                # alloy
                spcsheader = "'??'".ljust(50)+"! "+spcs+"\n"+str(spcsites)+"\n"
            else:
                spcsheader = "'"+spcs+"'\n"+str(spcsites)+"\n"
            filestring += spcsheader
            filestring += positionstring
        filestring += "\n"+self.docstring
        return filestring

################################################################################################
class ElkFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an elk.in file and the method
    __str__ that outputs the contents of an elk.in file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has the form of a f90 comment
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("!")
            string = "!"+string+"\n"
            self.docstring += string
    def __str__(self):
        filestring = self.docstring
        # Lattice vectors
        filestring += "avec\n"
        tmpstring = ""
        for pos in self.cell.latticevectors:
            for coord in pos:
                tmpstring += "  %13.10f"%coord
            tmpstring += "\n"
        tmpstring += "\n"
        filestring += tmpstring
        # Scale factor
        filestring += "scale\n"
        filestring += "  %13.10f\n\n"%self.cell.lengthscale
        # Atoms
        filestring += "atoms\n"
        # Get number of species
        species = set([])
        for a in self.cell.atomdata:
            for b in a:
                species.add(b.spcstring())
        tmpstring = ("  "+str(len(species))).ljust(37)+": nspecies\n"
        filestring += tmpstring
        # local B-field string
        bfcmtstring = "   0.00000000  0.00000000  0.00000000"
        # initialize some stuff
        natoms = 0
        spcstring = self.cell.atomdata[0][0].spcstring()
        positionstring = ""
        for a in self.cell.atomdata:
            for b in a:
                spcs = b.spcstring()
                # Accumulate for this species
                if spcs == spcstring:
                    ## natoms += len(a)
                    natoms += 1
                    positionstring += str(b.position)+bfcmtstring+"\n"
                else:
                    # Print species
                    if len(spcstring) > 2:
                        # alloy
                        filestring += "'??.in'".ljust(37)+": spfname = "+spcstring+"\n"
                    else:
                        filestring += ("'"+spcstring+".in'").ljust(37)+": spfname \n"
                    filestring += "  "+str(natoms)+"\n"
                    filestring += positionstring
                    # Initialize next species
                    spcstring = spcs
                    natoms = 1
                    positionstring = str(b.position)+bfcmtstring+"\n"
        # Print last species
        if len(spcstring) > 2:
            # alloy
            filestring += "'??.in'".ljust(37)+": spfname = "+spcstring+"\n"
        else:
            filestring += ("'"+spcstring+".in'").ljust(37)+": spfname\n"
        filestring += "  "+str(natoms)+"\n"
        filestring += positionstring
        return filestring


class Elk_input_File(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an elk.in file and the method
    __str__ that outputs the contents of an elk.in file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has the form of a f90 comment
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("!")
            string = "!"+string+"\n"
            self.docstring += string
    def __str__(self):
        #filestring = self.docstring
        tasks = self.tasks
        times = 1
        filestring  = "tasks\n"
        if tasks == "scf":
            filestring += "  0 : SCF calculation\n"
        elif tasks == "opt":
            filestring += "  2 : structure optimization\n"
            filestring += "\n"
            filestring += "mixtype\n"
            filestring += "  3 : Broyden mixing\n"
            filestring += "\n"
            filestring += "latvopt\n"
            filestring += "  1\n"
        elif tasks == "dos":
            filestring += "   0 : SCF calculation\n"
            filestring += "  10 : DOS calculation\n"
            times = 2
        elif tasks == "band":
            filestring += "   0 : SCF calculation\n"
            filestring += "  21 : Band calculation for every atom\n"
        elif tasks == "rescf+dos":
            filestring += "   1 : SCF recalculation (use STATE.out)\n"
            filestring += "  10 : DOS calculation\n"
            times = 2
        elif tasks == "rescf+band":
            filestring += "   1 : SCF recalculation (use STATE.out)\n"
            filestring += "  21 : Band calculation for every atom\n"
        elif tasks == "3dplot":
            filestring += "   0 : SCF calculation\n"
            filestring += "  33 : 3D plot\n"
            filestring += "\n"
            filestring += "plot3d\n"
            filestring += "  0.0 0.0 0.0 : vclp3d"
            filestring += "  1.0 0.0 0.0"
            filestring += "  0.0 1.0 0.0"
            filestring += "  0.0 0.0 1.0"
            filestring += "  56  56  56"
            filestring += "\n"
            filestring += "scale\n"
            filestring += "  1.0\n"
            filestring += "\n"
            filestring += "scal1\n"
            filestring += "  1.0\n"
            filestring += "\n"
            filestring += "scale2\n"
            filestring += "  1.0\n"
            filestring += "\n"
            filestring += "scale3\n"
            filestring += "  1.0\n"
        elif tasks == "rescf+3dplot":
            filestring += "   1 : SCF recalculation (use STATE.out)\n"
            filestring += "  33 : 3D plot\n"
            filestring += "\n"
            filestring += "plot3d\n"
            filestring += "  0.0 0.0 0.0 : vclp3d"
            filestring += "  1.0 0.0 0.0"
            filestring += "  0.0 1.0 0.0"
            filestring += "  0.0 0.0 1.0"
            filestring += "  56  56  56"
            filestring += "\n"
            filestring += "scale\n"
            filestring += "  1.0\n"
            filestring += "\n"
            filestring += "scal1\n"
            filestring += "  1.0\n"
            filestring += "\n"
            filestring += "scale2\n"
            filestring += "  1.0\n"
            filestring += "\n"
            filestring += "scale3\n"
            filestring += "  1.0\n"
        #
        if self.elk_dft_u or self.elk_dft_u_j:
            filestring += "\n"
            filestring += "dft+u\n"
            filestring += "  2 1 : dftu, inpdftu\n"     
            natoms = 0
            spcstring = ""
            for a in self.cell.atomdata:
                for b in a:
                    spcs = b.spcstring()
                    if spcs != spcstring:
                        natoms += 1
                        atomic_number = int(ed.elementnr[b.spcstring()])
                        if 21 <= atomic_number and atomic_number <= 29:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.81 + 0.08*(atomic_number-21)
                            U = F0
                            Ueff = U - J
                        elif 39 <= atomic_number and atomic_number <= 47:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.59 + 0.056*(atomic_number-39)
                            U = F0
                            Ueff = U - J
                        elif 71 <= atomic_number and atomic_number <= 79:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.860 + 0.053*(atomic_number-71)
                            U = F0
                            Ueff = U - J
                        elif 57 <= atomic_number and atomic_number <= 70:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.90 + 0.036*(atomic_number-57)
                            U = F0
                            Ueff = U - J
                        elif 89 <= atomic_number and atomic_number <= 103:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.66 + 0.035*(atomic_number-57)
                            U = F0
                            Ueff = U - J
                        else:
                            l = 0
                            U = 0.0
                            J = 0.0
                            Ueff = 0.0
                        #
                        if self.elk_dft_u:
                            Ueff = Ueff/27.12*float(self.elk_dftu_times)
                            J = 0.0
                            filestring += "%4i %i %5.3f %5.3f : is, l, U, J\n"%(natoms,l,Ueff,J)
                        else:
                            U = U/27.12*float(self.elk_dftu_times)
                            J = J/27.12*float(self.elk_dftu_times)
                            filestring += "%4i %i %5.3f %5.3f : is, l, U, J\n"%(natoms,l,U,J)
                    spcstring = spcs 
            filestring += "\n"
            filestring += "nempty\n"
            filestring += "  8 : number of empty states\n"
        elif tasks == "tddft-alda":
            filestring += "  0 : SCF calculation\n"
            filestring += "120 : Momentum matrix elements\n"
            filestring += "320 : TDDFT calculation\n"
            filestring += "\n"
            filestring += "vecql\n"
            filestring += "  0.5 0.5 0.5\n"
            filestring += "\n"
            filestring += "xctype\n"
            filestring += "  3 : 3=LDA, 20=PBE\n"
            filestring += "\n"
            filestring += "fxctype\n"
            filestring += "  3\n"
            filestring += "\n"
            filestring += "gmaxrf\n"
            filestring += "  2.0\n"
            filestring += "\n"
            filestring += "wplot\n"
            filestring += "  500 100 1\n"
            filestring += "  0.0 1.5\n"
        elif tasks == "tddft-lrc":
            filestring += "  0 : SCF calculation\n"
            filestring += "120 : Momentum matrix elements\n"
            filestring += "121 : RPA calculation\n"
            filestring += "320 : TDDFT calculation\n"
            filestring += "\n"
            filestring += "fxctype\n"
            filestring += "200\n"
            filestring += "\n"
            filestring += "fxclrc\n"
            filestring += "  5.5\n"
            filestring += "\n"
            filestring += "gmaxrf\n"
            filestring += "  1.0\n"
            filestring += "\n"
            filestring += "swidth\n"
            filestring += "  0.01\n"
            filestring += "\n"
            filestring += "nempty\n"
            filestring += "  8\n"
            filestring += "\n"
            filestring += "lradstp\n"
            filestring += "  2\n"
            filestring += "\n"
            filestring += "wplot\n"
            filestring += "  800 100 1\n"
            filestring += "  0.0 1.5\n"
        elif tasks == "tddft-bootstrap":
            filestring += "  0 : SCF calculation\n"
            filestring += "120 : Momentum matrix elements\n"
            filestring += "121 : RPA calculation\n"
            filestring += "320 : TDDFT calculation\n"
            filestring += "\n"
            filestring += "scissor\n"
            filestring += "0.192\n"
            filestring += "\n"
            filestring += "xctype\n"
            filestring += "  3\n"
            filestring += "\n"
            filestring += "fxctype\n"
            filestring += "210\n"
            filestring += "\n"
            filestring += "fxclrc\n"
            filestring += "  5.5\n"
            filestring += "\n"
            filestring += "gmaxrf\n"
            filestring += "  0.0\n"
            filestring += "\n"
            filestring += "swidth\n"
            filestring += "  0.01\n"
            filestring += "\n"
            filestring += "nempty\n"
            filestring += "  8\n"
            filestring += "\n"
            filestring += "lradstp\n"
            filestring += "  2\n"
            filestring += "\n"
            filestring += "wplot\n"
            filestring += "  800 100 1\n"
            filestring += "  0.0 1.5\n"
        elif tasks == "xas":
            filestring += "  0 : SCF calculation\n"
            filestring += "120 : Momentum matrix elements\n"
            filestring += "121 : RPA calculation\n"
            filestring += "\n"
            filestring += "wplot\n"
            filestring += " 400 200 0\n"
            filestring += " 25.0 26.5\n"
            filestring += "\n"
            filestring += "optcomp\n"
            filestring += " 1 1 : x x\n"
            filestring += " 1 2 : x y\n"
            filestring += " 3 3 : z z\n"
            filestring += "\n"
            filestring += "spinpol\n"
            filestring += "  .ture.\n"
            filestring += "\n"
            filestring += "spinorb\n"
            filestring += "  .ture.\n"
            filestring += "\n"
            filestring += "xctype\n"
            filestring += "  20 : 3=LDA, 20=PBE\n"
            filestring += "\n"
            filestring += "bfieldc\n"
            filestring += "  0.0 0.0 2.0\n"
            filestring += "\n"
            filestring += "reducebf\n"
            filestring += "  0.8\n"
            filestring += "\n"
            filestring += "nempty\n"
            filestring += "  20\n"
            filestring += "\n"
            filestring += "swidth\n"
            filestring += "  0.01\n"
        elif tasks == "bse":
            filestring += "  0 : SCF calculation\n"
            filestring += "120 : Momentum matrix elements\n"
            filestring += "121 : RPA calculation\n"
            filestring += "180 : RPA calculation\n"
            filestring += "185 : BSE Hamiltonian matrix\n"
            filestring += "186 : Diagonalize BSE matrix\n"
            filestring += "187 : Generate BSE dielectric function\n"
            filestring += "\n"
            filestring += "scissor\n"
            filestring += "  0.21\n"
            filestring += "\n"
            filestring += "lmaxvr\n"
            filestring += "  5\n"
            filestring += "\n"
            filestring += "gmaxvr\n"
            filestring += "  0.0\n"
            filestring += "\n"
            filestring += "mvbse\n"
            filestring += "  3\n"
            filestring += "\n"
            filestring += "gmaxrf\n"
            filestring += "  3.0\n"
            filestring += "\n"
            filestring += "swidth\n"
            filestring += "  0.005\n"
            filestring += "\n"
            filestring += "nempty\n"
            filestring += "  20\n"
            filestring += "\n"
            filestring += "wplot\n"
            filestring += "  800 100 0 : nwplot, ngrkf, nswplot\n"
            filestring += "  0.0 1.5 : wplot\n"
        elif tasks == "bse-core":
            filestring += "  0 : SCF calculation\n"
            filestring += "120 : Momentum matrix elements\n"
            filestring += "121 : RPA calculation\n"
            filestring += "180 : RPA calculation\n"
            filestring += "185 : BSE Hamiltonian matrix\n"
            filestring += "186 : Diagonalize BSE matrix\n"
            filestring += "187 : Generate BSE dielectric function\n"
            filestring += "\n"
            filestring += "nvbser\n"
            filestring += "  0\n"
            filestring += "\n"
            filestring += "istxbse\n"
            filestring += "  1\n"
            filestring += "  2\n"
            filestring += "  3\n"
            filestring += "  4\n"
            filestring += "  5\n"
            filestring += "  6\n"
            filestring += "\n"
            filestring += "ncbse\n"
            filestring += "  20\n"
            filestring += "\n"
            filestring += "nempty\n"
            filestring += "  20\n"
            filestring += "\n"
            filestring += "gmaxrf\n"
            filestring += "  5.0\n"
            filestring += "\n"
            filestring += "gmaxvr\n"
            filestring += "  0.0\n"
            filestring += "\n"
            filestring += "vkloff\n"
            filestring += "  0.05 0.15 0.25\n"
            filestring += "\n"
            filestring += "mixtype\n"
            filestring += "  3 : Broyden mixing\n"
            filestring += "\n"
            filestring += "wplot\n"
            filestring += "  5000 100 0 : nwplot, ngrkf, nswplot\n"
            filestring += "  0.0 14.0 : wplot\n"
            filestring += "\n"
            filestring += "swidth\n"
            filestring += "  0.01\n"
            filestring += "spinorb\n"
            filestring += "  .ture.\n"
        elif tasks == "elnes":
            filestring += "  0 : SCF calculation\n"
            filestring += "140\n"
            filestring += "\n"
            filestring += "spinorb\n"
            filestring += "  .ture.\n"
            filestring += "\n"
            filestring += "emaxelnes\n"
            filestring += "  -30.0\n"
            filestring += "\n"
            filestring += "wplot\n"
            filestring += "  500 200 0 : nwplot, ngrkf, nswplot\n"
            filestring += "  0.0 30.0 : wplot\n"
            filestring += "\n"
            filestring += "vecql\n"
            filestring += "  0.0 0.0 0.1\n"
            filestring += "\n"
            filestring += "rgkmax\n"
            filestring += "  8.0\n"
            filestring += "\n"
            filestring += "gkmaxvr\n"
            filestring += "  14.0\n"
            filestring += "\n"
            filestring += "nempty\n"
            filestring += "  20\n"
        elif tasks == "hybrid":
            filestring += "  0 : SCF calculation\n"
            filestring += "  5\n"
            filestring += "\n"
            filestring += "xctype\n"
            filestring += " 20\n"
            filestring += "\n"
            filestring += "hybrid\n"
            filestring += "  .ture.\n"
            filestring += "\n"
            filestring += "! hybridc\n"
            filestring += "\n"
            filestring += "hybmix\n"
            filestring += "  0.25\n"
            filestring += "nempty\n"
            filestring += " 8\n"
        elif tasks == "oep+band":
            filestring += "  0 : SCF calculation\n"
            filestring += " 20 : Band calculation\n"
            filestring += "\n"
            filestring += "xctype\n"
            filestring += " -1\n"
            filestring += "\n"
            filestring += "maxitoep\n"
            filestring += "  200\n"
            filestring += "\n"
            filestring += "nempty\n"
            filestring += "  15\n"
        else:
            pass
        filestring += "\n"
        #
        filestring += "spinpol\n"
        filestring += "  .true. : Spin Polarized calculation\n"
        filestring += "\n"
        #
        atomic_number = 0
        spinorb = ".false."
        for a in self.cell.atomdata:
            for b in a:
                atomic_number = int(ed.elementnr[b.spcstring()])
                if atomic_number >= 50:
                    spinorb = ".true."
        if spinorb == ".true.":
            filestring += "spinorb\n"
            filestring += "  .true. : atomic number >= 50\n"
            filestring += "\n"
        #
        # POT library
        Elk_pot_dir = self.elkpotlib
        if self.elkpotlib == "":
            try:
                Elk_pot_dir = os.environ['Elk_POTLIB']
            except:
                #filestring += "! same as Example file adress\n"
                Elk_pot_dir = "../../species/"
        filestring += "sppath\n"
        filestring += "   '"+str(Elk_pot_dir)+"'\n"
        #
        filestring += "\n"
        # Lattice vectors
        filestring += "avec\n"
        tmpstring = ""
        for pos in self.cell.latticevectors:
            for coord in pos:
                tmpstring += "  %13.10f"%coord
            tmpstring += "\n"
        tmpstring += "\n"
        filestring += tmpstring
        # Scale factor
        filestring += "scale\n"
        filestring += "  %13.10f\n\n"%self.cell.lengthscale
        # Atoms
        filestring += "atoms\n"
        # Get number of species
        species = set([])
        for a in self.cell.atomdata:
            for b in a:
                species.add(b.spcstring())
        tmpstring = ("  "+str(len(species))).ljust(37)+": nspecies\n"
        filestring += tmpstring
        # local B-field string
        bfcmtstring = "   0.00000000  0.00000000  0.00000000"
        # initialize some stuff
        natoms = 0
        spcstring = self.cell.atomdata[0][0].spcstring()
        positionstring = ""
        for a in self.cell.atomdata:
            for b in a:
                spcs = b.spcstring()
                # Accumulate for this species
                if spcs == spcstring:
                    ## natoms += len(a)
                    natoms += 1
                    positionstring += str(b.position)+bfcmtstring+"\n"
                else:
                    # Print species
                    if len(spcstring) > 2:
                        # alloy
                        filestring += "'??.in'".ljust(37)+": spfname = "+spcstring+"\n"
                    else:
                        filestring += ("'"+spcstring+".in'").ljust(37)+": spfname \n"
                    filestring += "  "+str(natoms)+"\n"
                    filestring += positionstring
                    # Initialize next species
                    spcstring = spcs
                    natoms = 1
                    positionstring = str(b.position)+bfcmtstring+"\n"
        # Print last species
        if len(spcstring) > 2:
            # alloy
            filestring += "'??.in'".ljust(37)+": spfname = "+spcstring+"\n"
        else:
            filestring += ("'"+spcstring+".in'").ljust(37)+": spfname\n"
        filestring += "  "+str(natoms)+"\n"
        filestring += positionstring
        #
        filestring += "\n"
        filestring += "epsengy\n"
        filestring += "  %-4.1e : Ry unit (ca. 1 meV/atom)\n"%(7.35e-5*self.cell.natoms())
        #
        filestring += "\n"
        filestring += "ngridk\n"
        reclatvect = self.cell.reciprocal_latticevectors()
        for j in range(3):
            for i in range(3):
                reclatvect[j][i] = reclatvect[j][i] / self.cell.lengthscale
        # Lengths of reciprocal lattice vectors
        reclatvectlen = [elem.length() for elem in reclatvect]
        self.kgrid = [max(1,int(round(elem/0.52917/self.kresolution))) for elem in reclatvectlen]
        filestring += "  "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"\n"
        #
        if self.kpeven == "yes":
            if self.cell.crystal_system() == "hexagonal":
                kshift_list = [2]
            elif self.cell.crystal_system() == "tetragonal":
                kshift_list = []
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                kshift_list = [0,1,2]
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                kshift_list = [0,1,2]
            elif self.cell.spacegroupsetting == "P":
                kshift_list = [0,1,2]
            else:
                kshift_list = []
            #
            for k_no in kshift_list:
                if self.kgrid[k_no]%2 != 0 and self.kgrid[k_no] > 1:
                    self.kgrid[k_no] = int(self.kgrid[k_no]-0.1)
                #if self.qgrid[k_no] >= self.kgrid[k_no]:
                #    self.qgrid[k_no] = self.kgrid[k_no]
                #elif self.qgrid[k_no]%2 != 0:
                #    self.qgrid[k_no] = self.kgrid[k_no]/2
        #
        filestring += "\n"
        filestring += "vkloff\n"
        if self.cell.crystal_system() == "hexagonal":
            filestring += "  0.0 0.0 0.5\n"
        elif self.cell.crystal_system() == "tetragonal":
            filestring += "  0.0 0.0 0.0\n"
        elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
            filestring += "  0.5 0.5 0.5\n"
        elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
            filestring += "  0.5 0.5 0.5\n"
        elif self.cell.spacegroupsetting == "P":
            filestring += "  0.5 0.5 0.5\n"
        else:
            filestring += "  0.0 0.0 0.0\n"
        filestring += "\n"        
        #
        # symmetry
        if tasks == "band" or tasks == "oep+band":
            filestring += "plot1d\n"
            if self.cell.crystal_system() == "hexagonal":
                #brvtyp = "hcp", primitive basis
                filestring += "5 250 : nvp1d, npp1d\n"
                filestring += "0.0   0.0   0.0  : vlvp1d, G-point\n"
                filestring += "0.5   0.0   0.0  : M-point\n"
                filestring += "0.333 0.333 0.0  : K-point\n"
                filestring += "0.0   0.0   0.0  : G-point\n"
                filestring += "0.0   0.0   0.5  : A-point\n"
            elif self.cell.spacegroupsetting == "F":
                #brvtyp = "fcc", primitive basis
                filestring += "6 300 : nvp1d, npp1d\n"
                filestring += "0.5  0.25  0.75 : vlvp1d, W-point\n"
                filestring += "0.5  0.0   0.0  : L-point\n"
                filestring += "0.0  0.0   0.0  : G-point\n"
                filestring += "0.5  0.5   0.0  : X-point\n"
                filestring += "0.75 0.5   0.25 : W-point\n"
                filestring += "0.75 0.375 0.375 : K-point\n"
            elif self.cell.spacegroupsetting == "I":
                #brvtyp = "bcc", primitive basis
                filestring += "5 250 : nvp1d, npp1d\n"
                filestring += "0.0  0.0  0.0  : vlvp1d, G-point\n"
                filestring += "0.5 -0.5  0.5  : H-point\n"
                filestring += "0.0  0.0  0.5  : N-point\n"
                filestring += "0.0  0.0  0.0  : G-point\n"
                filestring += "0.25 0.25 0.25 : P-point\n"
            elif self.cell.spacegroupsetting == "P":
                #brvtyp = "sc", primitive basis
                filestring += "5 250 : nvp1d, npp1d\n"
                filestring += "0.5  0.5  0.5  : vlvp1d, R-point\n"
                filestring += "0.0  0.0  0.0  : G-point\n"
                filestring += "0.5  0.0  0.0  : X-point\n"
                filestring += "0.5  0.5  0.0  : M-point\n"
                filestring += "0.0  0.0  0.0  : G-point\n"
            else:
                #brvtyp = "sc", primitive basis
                filestring += "5 250 : nvp1d, npp1d\n"
                filestring += "0.5  0.5  0.5  : vlvp1d, R-point\n"
                filestring += "0.0  0.0  0.0  : G-point\n"
                filestring += "0.5  0.0  0.0  : X-point\n"
                filestring += "0.5  0.5  0.0  : M-point\n"
                filestring += "0.0  0.0  0.0  : G-point\n"
            #
            filestring += "\n"
        #
        return filestring

################################################################################################
class ExcitingFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an input.xml file and the method
    __str__ that outputs the contents of an input.xml file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.title = ""
        self.docstring = self.docstring.rstrip("\n")+"\n"
    def __str__(self):
        filestring = "<input>\n"
        filestring += "  <title>\n"
        filestring += self.docstring
        filestring += "  </title>\n"
        # Add title if there is one
        if self.title != "":
            filestring += "  <title>"+self.title+"</title>\n"
        filestring += "  <structure>\n"
        # scale factor
        filestring += "    <crystal scale="+str(self.cell.lengthscale)+">\n"
        # Lattice vectors
        tmpstring = ""
        for pos in self.cell.latticevectors:
            tmpstring += "      <basevect>"
            for coord in pos:
                tmpstring += " %13.10f"%coord
            tmpstring += "</basevect>\n"
        filestring += tmpstring
        filestring += "    </crystal>\n"
        # Atoms
        # local B-field string
        bfcmtstring = "   0.00000000  0.00000000  0.00000000"
        # initialize some stuff
        spcstring = self.cell.atomdata[0][0].spcstring()
        positionstring = ""
        for a in self.cell.atomdata:
            for b in a:
                spcs = b.spcstring()
                # Accumulate for this species
                if spcs == spcstring:
                    positionstring += "      <atom coord=\""
                    positionstring += str(b.position)+"\"/>\n"
                else:
                    # Print species
                    if len(spcstring) > 2:
                        # alloy
                        filestring += "    <species speciesfile=\"??.xml\">"
                        filestring += "       <!-- "+spcstring+" -->\n"
                    else:
                        filestring += "    <species speciesfile=\""+spcstring+".xml\">\n"
                    filestring += positionstring+"    </species>\n"
                    # Initialize next species
                    spcstring = spcs
                    positionstring = "      <atom coord=\""+str(b.position)+"\"/>\n"
        # Print last species
        if len(spcstring) > 2:
            # alloy
            filestring += "    <species speciesfile=\"??.xml\">"
            filestring += "       <!-- "+spcstring+" -->\n"
        else:
            filestring += "    <species speciesfile=\""+spcstring+".xml\">\n"
        filestring += positionstring
        filestring += "    </species>\n"
        filestring += "  </structure>\n"
        filestring += "</input>\n"
        return filestring

################################################################################################
class FleurFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a Fleur input generator input file (how about
    that, we generate input for the generator of the input...) and the method
    __str__ that outputs the contents as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # make sure the docstring goes on one line
        self.docstring = self.docstring.replace("\n"," ")
        if len(self.docstring) > 80:
            self.docstring = self.docstring[0:78]+"...\n"
    def __str__(self):
        ed = ElementData()
        filestring = self.docstring+"\n"
        filestring += "&input cartesian=f oldfleur=f\n\n"
        # Lattice vectors
        tmpstring = ""
        n = 1
        for pos in self.cell.latticevectors:
            tmpstring += str(pos)
            tmpstring += "    !  a%1i\n"%n
            n += 1
        filestring += tmpstring
        # Scale factor
        filestring += "%13.9f    ! aa\n"%self.cell.lengthscale
        filestring += "1.0  1.0  1.0 ! scale(1), scale(2), scale(3)\n"
        # Atoms
        natom = 0
        for a in self.cell.atomdata:
            natom += len(a)
        filestring += str(natom)+"\n"
        nspcs = 0
        spcs = ""
        coordstring = ""
        for a in self.cell.atomdata:
            for b in a:
                # Check for alloy
                if b.alloy():
                    prestring = "??"
                    poststring = "  ! "
                    for k in b.species:
                        poststring += str(ed.elementnr[k])+"/"
                    poststring = poststring.rstrip("/")+" "+str(b.spcstring())+"\n"
                else:
                    prestring = str(ed.elementnr[b.spcstring()]).ljust(2)
                    poststring = "  ! "+b.spcstring()+"\n"
                coordstring += prestring+str(b.position)+poststring
        # To filestring
        filestring += coordstring
        return filestring

################################################################################################
class CASTEPFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a CASTEP run and the method
    __str__ that outputs to a .cell file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Cartesian units?
        self.cartesian = False
        # What units?
        self.unit = "angstrom"
        self.cell.newunit("angstrom")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        # VCA calculation?
        self.vca = False
        # Print labels?
        self.printlabels = False
    def __str__(self):
        # Set units
        self.cell.newunit(self.unit)
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        filestring = self.docstring+"\n"
        filestring += "%BLOCK LATTICE_CART\n"
        # units
        if self.cell.unit == "angstrom":
            filestring += "ang    # angstrom units\n"
        elif self.cell.unit == "bohr":
            filestring += "bohr   # atomic units\n"
        # lattice
        for vec in lattice:
            for coord in vec:
                filestring += " %19.15f"%(coord*a)
            filestring += "\n"
        # Cutoff
        filestring += "%ENDBLOCK LATTICE_CART\n\n"
        # The atom position info
        if self.cartesian:
            # Correct block name and units
            filestring += "%BLOCK POSITIONS_ABS\n"
            if self.cell.unit == "angstrom":
                filestring += "ang    # angstrom units\n"
            elif self.cell.unit == "bohr":
                filestring += "bohr   # atomic units\n"
            # Set transformation matrix
            transmat = LatticeMatrix(self.cell.latticevectors)
            scalfac = self.cell.a
        else:
            filestring += "%BLOCK POSITIONS_FRAC\n"
            transmat = LatticeMatrix([[1,0,0],[0,1,0],[0,0,1]])
            scalfac = 1.0
        i = 0
        for a in self.cell.atomdata:
            for b in a:
                pos = Vector(mvmult3(transmat,b.position)).scalmult(scalfac)
                # Check for VCA calculation
                if self.cell.alloy and self.vca:
                    if len(b.species) > 1:
                        i = i + 1
                        for sp,conc in b.species.iteritems():
                            filestring += sp.ljust(2)+" "+str(pos)+"  MIXTURE:( %i %6.5f )"%(i,conc)
                    else:
                        filestring += b.spcstring().ljust(2)+" "+str(pos)
                else:
                    filestring += b.spcstring().ljust(2)+" "+str(pos)
                if self.printlabels and b.label != "":
                    filestring += " ID="+b.label
                filestring += "\n"
        if self.cartesian:
            filestring += "%ENDBLOCK POSITIONS_ABS\n"
        else:
            filestring += "%ENDBLOCK POSITIONS_FRAC\n"
        # pseudo-potential block
        species = set([])
        for a in self.cell.atomdata:
            if self.vca:
                for sp,conc in a[0].species.iteritems():
                    species.add(sp)
            else:
                species.add(a[0].spcstring())
        filestring += "\n"
        filestring += "# Commented out pseudopotential block for easy editing\n"
        filestring += "#%BLOCK SPECIES_POT\n"
        for sp in species:
            filestring += "# "+sp.ljust(2)+"  "+sp+"_00.usp\n"
        filestring += "#%ENDBLOCK SPECIES_POT\n"
        # Put in the symmetry operations
        filestring += "\n%BLOCK SYMMETRY_OPS\n"
        latvect = self.cell.conventional_latticevectors()
        # make list and make sure that identity comes first
        symoplist = sorted(list(self.cell.symops))
        k = 1
        for op in symoplist:
            filestring += "# Symm. op. %i\n"%k
            filestring += str(op)
            k += 1
        filestring += "%ENDBLOCK SYMMETRY_OPS\n"        
        return filestring

################################################################################################
# PWSCF (Quantum Espresso)
class PWSCFFile(GeometryOutputFile):
    """
    Class for storing the geometrical data for a PWSCF run and the method
    __str__ that outputs to a .in file as a string.
    """
    def __init__(self, crystalstructure, string, kresolution=0.2):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        #
        self.setupall = False
        # Cartesian units?
        self.cartesian = False
        self.cartesianpositions = False
        self.cartesianlatvects = False
        self.scaledcartesianpositions = False
        # What units?
        self.unit = "angstrom"
        self.cell.newunit("angstrom")
        # Pseudopotential string
        self.pseudostring = "_PSEUDO"
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                tmp.add(b.spcstring())
        self.species = list(tmp)
        # k-space information
        reclatvect = self.cell.reciprocal_latticevectors()
        for j in range(3):
            for i in range(3):
                reclatvect[j][i] = reclatvect[j][i] / self.cell.lengthscale
        # Lengths of reciprocal lattice vectors
        reclatvectlen = [elem.length() for elem in reclatvect]
        self.kgrid = [max(1,int(round(elem/kresolution))) for elem in reclatvectlen]
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        self.docstring += "\n"
    def __str__(self):
        filestring = self.docstring
        # Set current units and stuff
        if self.cartesian:
            self.cartesianpositions = True
            self.cartesianlatvects = True
        self.cell.newunit(self.unit)
        # Determine max width of spcstring
        width = 0
        for a in self.cell.atomdata:
            for b in a:
                width = max(width, len(b.spcstring()))
        #
        filestring += "&SYSTEM\n"
        filestring += "  ibrav = %i\n"%(0)
        if self.unit == "bohr":
            filestring += "  celldm(1) = %10.5f\n"%(self.cell.lengthscale)
        elif self.unit == "angstrom":
            filestring += "  A = %10.5f\n"%(self.cell.lengthscale)
        filestring += "  nat = %i\n"%(self.cell.natoms())
        filestring += "  ntyp = %i\n"%(len(self.species))
        filestring += "/\n"
        if self.cartesianlatvects:
            if self.unit == "bohr":
                filestring += "CELL_PARAMETERS {bohr}\n"
            elif self.unit == "angstrom":
                filestring += "CELL_PARAMETERS {angstrom}\n"
            t = LatticeMatrix(self.cell.latticevectors)
            for i in range(3):
                for j in range(3):
                    t[i][j] = self.cell.latticevectors[i][j]*self.cell.lengthscale
            filestring += str(t)
        else:
            filestring += "CELL_PARAMETERS {alat}\n"
            filestring += str(self.cell.latticevectors)
        filestring += "ATOMIC_SPECIES\n"
        for sp in self.species:
                filestring += "  %2s"%(sp.rjust(width))
                try:
                    filestring += ("  %8.5f"%(ed.elementweight[sp])).rjust(11)
                except:
                    filestring += "   ???".rjust(11)
                filestring += "  %2s%s\n"%(sp.rjust(width),self.pseudostring)
        if self.cartesianpositions:
            if self.scaledcartesianpositions:
                filestring += "ATOMIC_POSITIONS {alat}\n"
            else:
                if self.unit == "bohr":
                    filestring += "ATOMIC_POSITIONS {bohr}\n"
                elif self.unit == "angstrom":
                    filestring += "ATOMIC_POSITIONS {angstrom}\n"
        else:
            if self.scaledcartesianpositions:
                filestring += "ATOMIC_POSITIONS {alat}\n"
            else:
                filestring += "ATOMIC_POSITIONS {crystal}\n"
        for a in self.cell.atomdata:
            for b in a:
                if self.cartesianpositions:
                    t = Vector(mvmult3(self.cell.latticevectors,b.position))
                    if self.scaledcartesianpositions:
                        filestring += b.spcstring().rjust(width)+" "+str(t)+"\n"
                    else:
                        for i in range(3):
                            t[i] = self.cell.lengthscale*t[i]
                        filestring += b.spcstring().rjust(width)+" "+str(t)+"\n"
                else:
                    if self.scaledcartesianpositions:
                        t = Vector(mvmult3(self.cell.latticevectors,b.position))
                        filestring += b.spcstring().rjust(width)+" "+str(t)+"\n"
                    else:
                        filestring += b.spcstring().rjust(width)+" "+str(b.position)+"\n"
        # Add k-space mesh
        if self.setupall:
            filestring += "\n# k-space resolution ~"+str(self.kresolution)+"/A.\n"
            # Opt for gamma-point run if possible
            if self.kgrid[0]*self.kgrid[1]*self.kgrid[2] == 1:                
                filestring += "K_POINTS gamma\n"
            else:
                filestring += "K_POINTS automatic\n"
                filestring += str(self.kgrid[0])+" "+str(self.kgrid[1])+" "+str(self.kgrid[2])+"  0 0 0\n"
        return filestring
    # Return the PWscf internal bravais lattice number
    def ibrav(self):
        system = self.cell.crystal_system()
        setting = self.cell.spacegroupsetting
        if self.supercell:
            return 14
        if system == 'cubic':
            if self.primcell:
                if setting == 'P':
                    return 1
                elif setting == 'F':
                    return 2
                elif setting == 'I':
                    return 3
            else:
                return 1
        if system == 'hexagonal':
            if self.primcell:
                if setting == 'P':
                    return 4
                elif setting == 'R':
                    return 5


class PWSCF_Input_File(GeometryOutputFile):
    """
    Class for storing the geometrical data for a PWSCF run and the method
    __str__ that outputs to a .in file as a string.
    """
    def __init__(self, crystalstructure, string, kresolution=0.2, qresolution=0.4, kpeven="yes"):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        #
        self.setupall = False
        # Cartesian units?
        self.cartesian = False
        self.cartesianpositions = False
        self.cartesianlatvects = False
        self.scaledcartesianpositions = False
        # What units?
        self.unit = "angstrom"
        self.cell.newunit("angstrom")
        # Pseudopotential string
        #self.pseudostring = "_PSEUDO"
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                tmp.add(b.spcstring())
        self.species = list(tmp)
        # k-space information
        reclatvect = self.cell.reciprocal_latticevectors()
        for j in range(3):
            for i in range(3):
                reclatvect[j][i] = reclatvect[j][i] / self.cell.lengthscale
        # Lengths of reciprocal lattice vectors
        reclatvectlen = [elem.length() for elem in reclatvect]
        self.kgrid = [max(1,int(round(elem/kresolution))) for elem in reclatvectlen]
        self.qgrid = [max(1,int(round(elem/qresolution))) for elem in reclatvectlen]
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        self.docstring += "\n"
        self.ibrav
    def __str__(self):
        #filestring = self.docstring
        filestring = ""
        if str(self.run_type) == "neb-start":
            filestring += "BEGIN\n"
            filestring += "BEGIN_PATH_INPUT\n"
            filestring += "&PATH\n"
            filestring += "  restart_mode      = 'from_scratch'\n"
            filestring += "  string_method     = 'neb',\n"
            filestring += "  nstep_path        = 20,\n"
            filestring += "  ds                = 2.D0,\n"
            filestring += '  opt_scheme        = "broyden",\n'
            filestring += "  num_of_images     = 7,\n"
            filestring += "  k_max             = 0.3D0,\n"
            filestring += "  k_min             = 0.2D0,\n"
            filestring += '  CI_scheme         = "auto",\n'
            filestring += "  path_thr          = 0.1D0,\n"
            if self.neb_mu:
                filestring += "  lfcpopt           = .TRUE.\n"
                filestring += "  fcp_mu            = %-8.6f,\n"%((float(self.neb_mu)+float(self.neb_bias_v))/13.60568)
                filestring += "  fcp_tot_charge_first = %-8.6f,\n"%(-0.021744/0.5*float(self.neb_bias_v))
                filestring += "  fcp_tot_charge_last  = %-8.6f,\n"%(-0.021744/0.5*float(self.neb_bias_v))
            filestring += "/\n"
            filestring += "END_PATH_INPUT\n"
            filestring += "BEGIN_ENGINE_INPUT\n"
        # Set current units and stuff
        if self.cartesian:
            self.cartesianpositions = True
            self.cartesianlatvects = True
        self.cell.newunit(self.unit)
        # Determine max width of spcstring
        width = 0
        for a in self.cell.atomdata:
            for b in a:
                width = max(width, len(b.spcstring()))
        #
        filestring += "&CONTROL\n"
        prefix_name = self.filename
        if str(self.run_type) == "dos":
            filestring += "  calculation  = 'nscf' ,\n"
        elif str(self.run_type) == "tddft" or str(self.run_type) == "eels":
            filestring += "  calculation  = 'scf' ,\n"            
        elif str(self.run_type) == "phonon":
            filestring += "  calculation  = 'scf' ,\n"
            filestring += "  tstress = .true. ,\n"
            filestring += "  tprnfor = .true. ,\n"
        elif str(self.run_type) == "lead" or str(self.run_type) == "scat" or str(self.run_type) == "lead_left" or str(self.run_type) == "lead-right":
            filestring += "  calculation  = 'scf' ,\n"
            prefix_name = str(self.run_type)
        elif str(self.run_type) == "neb-start" or str(self.run_type) == "neb-end":
            prefix_name = "neb"
        elif str(self.run_type) == "cp" or str(self.run_type) == "vc-cp":
            filestring += "  calculation  = '"+self.run_type+"' ,\n"    
            filestring += "  restart_mode='from_scratch',\n"
            filestring += "  nstep=1000, iprint=100, isave=1000,\n"
            if self.eDFT == "yes":
                filestring += "  dt=125.0,\n"
                filestring += "  ndr=50, ndw=51,\n"
            else:
                filestring += "  dt=3.0,\n"
                filestring += "  ndr=90, ndw=91,\n"
            filestring += "  verbosity='medium' ,\n"
        elif str(self.run_type) == "relax" or str(self.run_type) == "vc-relax":
            filestring += "  calculation  = '"+self.run_type+"' ,\n"
            filestring += "  nstep=100,\n"
            filestring += "  forc_conv_thr = 1.0e-03,\n"
        elif str(self.run_type) == "md":
            filestring += "  calculation  = '"+self.run_type+"' ,\n"
            filestring += "  dt=20,\n"
            filestring += "  nstep=50,\n"
            filestring += "  forc_conv_thr = 1.0e-03,\n"
        elif str(self.run_type) == "vc-md":
            filestring += "  calculation  = '"+self.run_type+"' ,\n"
            vcmd_dt = 50.0 / 70.0 # Wentzcovitch
            #vcmd_dt = 100.0 / 150.0 # damping algorithm
            filestring += "  dt="+str(vcmd_dt)+",\n"
            filestring += "  nstep=50,\n"
            filestring += "  forc_conv_thr = 1.0e-03,\n"
        else:
            filestring += "  calculation  = '"+self.run_type+"' ,\n"
        if self.stress == "yes" and str(self.run_type) != "phonon":
            filestring += "  tstress = .true. ,\n"
        if self.force == "yes" and str(self.run_type) != "phonon":
            filestring += "  tprnfor = .true. ,\n"
        filestring += "  prefix  = '"+prefix_name+"' ,\n"
        filestring += "  outdir  = './work/' ,\n"
        #
        # POT library
        pwscf_pot_dir = self.pseudolib
        if self.pseudolib == "":
            try:
                pwscf_pot_dir = os.environ['PWscf_PSEUDOLIB']
            except:
                try:
                    pwscf_pot_dir = os.environ['PWscf_PAWLIB']
                except:
                    #filestring += "! same as Example file adress\n"
                    pwscf_pot_dir = "../pseudo/"
        if not(self.setupall):
           filestring += "  pseudo_dir = './' ,\n"
        else:
           filestring += "  pseudo_dir = '"+str(pwscf_pot_dir)+"' ,\n"
        filestring += "/\n"
        #
        # POTCAR libraryself.pseudostring
        directory = str(pwscf_pot_dir)
        if directory != "":
            self.dir = directory
        else:
            try:
                self.dir = os.environ['PWscf_PSEUDOLIB']
            except:
                try:
                    self.dir = os.environ['PWscf_PAWLIB']
                except:
                    self.dir = ""
        # check directory
        if self.dir == "":
            print "No path to the PWscf pseudopotential library specified.\n"
        if not os.path.exists(self.dir):
            print "The specified path to the PWscf pseudopotential library does not exist.\n"+self.dir
        # get all files (GBRV, PSLibrary)
        potname = ""
        potfile = ""
        potlist = []
        potnamelist = []
        ATOMIC_SPECIES_string = ""
        if self.dfttype.lower() == "lda":
            if self.pseudotype == "PSLibrary":
                dft = "pz"
            if self.pseudotype == "GBRV":
                dft = "lda"
        elif self.dfttype.lower() == "pbe":
            dft = "pbe"
        elif self.dfttype.lower() == "pbesol":
            dft = "pbesol"
        elif self.pseudotype == "PSLibrary" and self.dfttype.lower() == "bp" or self.dfttype.lower() == "pw91" or self.dfttype.lower() == "revpbe" or self.dfttype.lower() == "pz" or self.dfttype.lower() == "wc":
            dft = self.dfttype.lower()
        else:
            dft = "pbe"
        if str(self.pslrel) == "yes":
            spslrel_list = ["rel-",""]
        else:
            spslrel_list = [""]
        for sp in self.species:
            for version in self.prioritylist:
                if self.pseudotype == "GBRV":
                    potname = sp.lower()+"_"+dft+"_"+version+str(self.pseudostring)
                    potfile = self.dir+"/"+potname
                    if os.path.exists(potfile):
                        potlist.append(potfile)
                        potnamelist.append(potfile)
                        break
                elif self.pseudotype == "PSLibrary":
                    for srel in spslrel_list:
                        for Next_version in self.prioritylist:
                            potname = sp+"."+srel+dft+version+Next_version+"_psl.1.0.0"+str(self.pseudostring)
                            potfile = self.dir+"/"+potname
                            if os.path.exists(potfile):
                                potlist.append(potfile)
                                potnamelist.append(potfile)
                                break
                            for n2 in range(3,-1,-1):
                                potname = sp+"."+srel+dft+version+Next_version+"_psl.0."+str(n2)+str(self.pseudostring)
                                potfile = self.dir+"/"+potname
                                if os.path.exists(potfile):
                                    potlist.append(potfile)
                                    potnamelist.append(potfile)
                                    break
                                for n3 in range(3,-1,-1):
                                    potname = sp+"."+srel+dft+version+Next_version+"_psl.0."+str(n2)+"."+str(n3)+str(self.pseudostring)
                                    potfile = self.dir+"/"+potname
                                    if os.path.exists(potfile):
                                        potlist.append(potfile)
                                        potnamelist.append(potfile)
                                        break
                                else:
                                    continue
                                break
                            else:
                                continue
                            break
                        else:
                            continue
                        break
                    else:
                        continue
                    break
            #
            ATOMIC_SPECIES_string += "  %2s"%(sp.rjust(width))
            try:
                ATOMIC_SPECIES_string += ("  %8.5f"%(ed.elementweight[sp])).rjust(11)
            except:
                ATOMIC_SPECIES_string += "   ???".rjust(11)
            if self.dir == "" or not os.path.exists(self.dir) or potname == "":
                ATOMIC_SPECIES_string += "  %2s%s\n"%(sp.rjust(width),self.pseudostring)
            else:
                ATOMIC_SPECIES_string += "  %s\n"%(potname)
        #
        # read potcar files and put in outstring
        outstring = ""
        i = 0
        max_ecutwfc = 0
        max_ecutrho = 0
        print ""
        print "setting pseudo-potential:"
        for f in potlist:
            print str(f)
            pot = open(f,"r")
            outstring = pot.read()
            pot.close()
            #
            pwscf_pot = open(potnamelist[i], "w")
            pwscf_pot.write(str(outstring))
            pwscf_pot.close()
            #
            fe = open(potnamelist[i],"r")
            lines = fe.readlines()
            fe.close()
            for line in lines:
                if search("wavefunctions:",line):
                    ecutwfc = float(line.split("wavefunctions:")[1].lstrip(" ").split()[0].strip(string.punctuation))
                    if ecutwfc >= max_ecutwfc:
                        max_ecutwfc = ecutwfc
                if search("density:",line):
                    ecutrho = float(line.split("density:")[1].lstrip(" ").split()[0].strip(string.punctuation))
                    if ecutrho >= max_ecutrho:
                        max_ecutrho = ecutrho
            i += 1
        #
        #
        filestring += "&SYSTEM\n"
        if self.space_group:
            filestring += "  space_group = "+str(self.cell.spacegroupnr)+",\n"
            filestring += "  origin_choice = 2,\n"
        if self.brav:
            filestring += "  ibrav = %i\n"%(self.ibrav())
            # http://web.mit.edu/espresso_v6.1/i386_linux26/qe-6.1/Doc/INPUT_PW.html
            if self.unit == "bohr":
                filestring += "  celldm(1) = %10.5f\n"%(self.cell.a)
                filestring += "  celldm(2) = %10.5f\n"%(self.cell.b/self.cell.a)
                filestring += "  celldm(3) = %10.5f\n"%(self.cell.c/self.cell.a)
                filestring += "  celldm(4) = %10.5f\n"%(cos(math.radians(self.cell.alpha)))
                filestring += "  celldm(5) = %10.5f\n"%(cos(math.radians(self.cell.beta)))
                filestring += "  celldm(6) = %10.5f\n"%(cos(math.radians(self.cell.gamma)))
            elif self.unit == "angstrom":
                filestring += "  A = %10.5f\n"%(self.cell.a)
                filestring += "  B = %10.5f\n"%(self.cell.b)
                filestring += "  C = %10.5f\n"%(self.cell.c)
                filestring += "  cosAB = %10.5f\n"%(cos(math.radians(self.cell.gamma)))
                filestring += "  cosAC = %10.5f\n"%(cos(math.radians(self.cell.beta)))
                filestring += "  cosBC = %10.5f\n"%(cos(math.radians(self.cell.alpha)))
        else:
            filestring += "  ibrav = %i\n"%(0)
            if self.unit == "bohr":
                filestring += "  celldm(1) = %10.5f\n"%(self.cell.lengthscale)
            elif self.unit == "angstrom":
                filestring += "  A = %10.5f\n"%(self.cell.lengthscale)
        filestring += "  nat = %i\n"%(self.cell.natoms())
        filestring += "  ntyp = %i\n"%(len(self.species))
        if max_ecutrho == 0:
            filestring += "  ecutwfc  = 40 ,\n"
            filestring += "  ecutrho  = 160 ,\n"
        else:
            filestring += "  ecutwfc  =  "+str(max_ecutwfc)+" ,\n"
            filestring += "  ecutrho  = "+str(max_ecutrho)+" ,\n"
        if self.dfttype.lower() == "b3lyp" or self.dfttype.lower() == "b3lyp-v1r" or self.dfttype.lower() == "pbe0" or self.dfttype.lower() == "hse" or self.dfttype.lower() == "gaup":
            filestring += "  input_dft = '"+self.dfttype+"',\n"
            filestring += "  nqx1=1 ,\n"
            filestring += "  nqx2=1 ,\n"
            filestring += "  nqx3=1 ,\n"
        elif self.dfttype.lower() != "pbe":
            filestring += "  input_dft = '"+self.dfttype+"',\n"
        if self.run_type == "dos":
            filestring += "  occupations = 'tetrahedra', \n"    
        elif self.run_type == "cp" or self.run_type == "vc-cp":  
            if self.eDFT == "yes":
                filestring += "  occupations = 'ensemble' ,\n"
                filestring += "  smearing = 'cs',\n"
                filestring += "  degauss  = 0.018 ,\n"
            filestring += "  nr1b=24, nr2b=24, nr3b=24,\n"
            for sp1 in self.species:
                if sp1=="Si":
                    for sp2 in self.species:
                        elem_flag = 1    
                        if sp2=="O":
                            filestring += "  qcutz=150., q2sigma=2.0, ecfixed=16.0,\n"
                            elem_flag = 0  
                    if elem_flag==1:
                        filestring += "  qcutz=12., q2sigma=4.0, ecfixed=12.0,\n"
        else:
            filestring += "  occupations = 'smearing' ,\n"
            filestring += "  degauss  = 0.02 ,\n"
            filestring += "  smearing = 'mp',\n"
        #
        if str(self.pslrel) == "yes":
            filestring += "  noncolin = .ture. ,\n"
            filestring += "  lspinorb = .ture. ,\n"    
            for i in range(len(self.species)):
                 filestring += "  starting_magnetization("+str(i+1)+") = 1 ,\n"
            filestring += "  tot_magnetization = -1 ,\n" 
        elif self.spin == "yes":
            filestring += "  nspin = 2 ,\n"
            for i in range(len(self.species)):
                 filestring += "  starting_magnetization("+str(i+1)+") = 1 ,\n"
            filestring += "  tot_magnetization = -1 ,\n"
        #
        if self.pwscf_dft_u or self.pwscf_dft_u_j:
            filestring += "  lda_plus_u = .true.,\n"
            spcstring = ""
            natoms = 0
            for a in self.cell.atomdata:
                for b in a:
                    spcs = b.spcstring()
                    if spcs != spcstring:
                        natoms += 1
                        atomic_number = int(ed.elementnr[b.spcstring()])
                        if 21 <= atomic_number and atomic_number <= 29:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.81 + 0.08*(atomic_number-21)
                            U = F0
                            Ueff = F0 - J
                        elif 39 <= atomic_number and atomic_number <= 47:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.59 + 0.056*(atomic_number-39)
                            U = F0
                            Ueff = F0 - J
                        elif 71 <= atomic_number and atomic_number <= 79:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.860 + 0.053*(atomic_number-71)
                            U = F0
                            Ueff = F0 - J
                        elif 57 <= atomic_number and atomic_number <= 70:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.90 + 0.036*(atomic_number-57)
                            U = F0
                            Ueff = F0 - J
                        elif 89 <= atomic_number and atomic_number <= 103:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.66 + 0.035*(atomic_number-57)
                            U = F0
                            Ueff = F0 - J
                        else:
                            l = 0
                            U = 0.0
                            J = 0.0
                            Ueff = 0.0
                        #
                        if self.pwscf_dft_u:
                            Ueff = Ueff*float(self.pwscf_dftu_times)
                            J = 0.0
                            filestring += "  Hubbard_U("+str(natoms)+") = %5.3f , \n"%Ueff
                        else:
                            U = U*float(self.pwscf_dftu_times)
                            J = J*float(self.pwscf_dftu_times)
                            filestring += "  Hubbard_U("+str(natoms)+") = %5.3f , \n"%U
                            filestring += "  Hubbard_J("+str(natoms)+") = %5.3f , \n"%J
                    spcstring = spcs
        #
        if self.run_type == "cp" or self.run_type == "vc-cp": 
            pass
        elif self.run_type == "relax" or self.run_type == "vc-relax":
            pass
        elif self.run_type == "md" or self.run_type == "vc-md":
            filestring += "  nosym=.true.,\n" # use MD or isolated atoms
        else:
            filestring += "  exxdiv_treatment = 'gygi-baldereschi' ,\n"
            filestring += "  x_gamma_extrapolation = .true. ,\n"
        #
        if self.neb_ems:
            filestring += "  assume_isolated='esm',\n"
            filestring += "  esm_bc='"+str(self.neb_ems)+"',\n"
        #
        filestring += "/\n"
        #
        #
        filestring += "&ELECTRONS\n"
        if self.run_type == "cp" or self.run_type == "vc-cp": 
            if self.eDFT == "yes":
                filestring += "  orthogonalization = 'Gram-Schmidt',\n"
                filestring += "  startingwfc='random', ampre=0.02, \n"
                filestring += "  tcg=.true., passop=0.3, maxiter=250,\n"
                filestring += "  emass=400., emass_cutoff=3.,\n"
                filestring += "  conv_thr = %-4.1e ,\n"%(7.35e-5*self.cell.natoms()) # Ry unit (ca. 1 meV/atom)
                filestring += "  n_inner = 2,\n"
                filestring += "  niter_cold_restart = 2,\n"
                filestring += "  lambda_cold = 0.03,\n"
            else:
                filestring += "  electron_dynamics='damp', electron_damping=0.15,\n"
                filestring += "  startingwfc='random', ampre=0.01,\n"
                filestring += "  emass=400., emass_cutoff=3.,\n"
        else:
            filestring += "  conv_thr = %-4.1e ,\n"%(7.35e-5*self.cell.natoms()) # Ry unit (ca. 1 meV/atom)
            filestring += "  mixing_beta = 0.7 ,\n"
        filestring += "/\n"
        #
        #
        if self.run_type == "cp" or self.run_type == "vc-cp":
            natoms = 0
            sum_atomic_mass_bc = 0.0
            for a in self.cell.atomdata:
                for b in a:
                    sp_b = b.spcstring()
                    for c in a:
                        sp_c = c.spcstring()
                        natoms += 1
                        atomic_mass_b   = float(ed.elementweight[sp_b])
                        atomic_mass_c   = float(ed.elementweight[sp_c])
                        sum_atomic_mass_bc  = sum_atomic_mass_bc + (atomic_mass_b*atomic_mass_c)/(atomic_mass_b+atomic_mass_c)
            avg_atomic_mass_bc = sum_atomic_mass_bc / float(natoms)
            fnosep =  0.01024+36.189*math.exp(-(avg_atomic_mass_bc-0.9482)/5.0085) + 1.0 
            fnoseh = (0.01024+36.189*math.exp(-(avg_atomic_mass_bc-0.9482)/5.0085))*(1.1**3.0) + 1.0
        if self.run_type == "cp":
            filestring += "&IONS\n"
            filestring += "  ion_dynamics='verlet', ion_temperature='nose',\n"
            filestring += "  fnosep = %-7.3f\n"%float(fnosep)
            filestring += "  tempw = "+str(self.temperature)+"\n"
            for i in range(len(self.species)):
                 filestring += "  ion_radius("+str(i+1)+") = 1.0 \n"
            filestring += "/\n"
        elif self.run_type == "vc-cp":
            filestring += "&IONS\n"
            filestring += "  ion_dynamics='verlet', ion_temperature='nose',\n"
            filestring += "  tempw = "+str(self.temperature)+", fnosep = %-7.3f\n"%float(fnosep)
            for i in range(len(self.species)):
                 filestring += "  ion_radius("+str(i+1)+") = 1.0,\n"
            filestring += "/\n"
            #
            filestring += "&CELL\n"
            filestring += "  cell_dynamics='pr', cell_temperature='nose',\n"
            filestring += "  temph = "+str(self.temperature)+", fnoseh = %-7.3f\n"%float(fnoseh)
            filestring += "  press=0.0,\n"
            filestring += "/\n"
        elif self.run_type == "relax":
            filestring += "&IONS\n"
            filestring += "  ion_dynamics = 'bfgs' ,\n"
            #filestring += "  ion_dynamics='damp',\n"
            filestring += "/\n"
        elif self.run_type == "vc-relax":
            filestring += "&IONS\n"
            filestring += "  ion_dynamics = 'bfgs' ,\n"
            #filestring += "  ion_dynamics='damp',\n"
            filestring += "/\n"
            #
            filestring += "&CELL\n"
            filestring += "  cell_dynamics='bfgs',\n"
            filestring += "  press_conv_thr = 0.5,\n" # Kbar
            filestring += "/\n"
        elif self.run_type == "md":
            filestring += "&IONS\n"
            filestring += "  tempw="+str(self.temperature)+", \n"
            filestring += "  ion_temperature = 'rescale-v' ,\n"
            #filestring += "  nraise = 1 ,\n"
            filestring += "  ion_dynamics = 'verlet' ,\n"
            #filestring += "  pot_extrapolation = 'second-order' ,\n"
            #filestring += "  wfc_extrapolation = 'second-order' ,\n"
            filestring += "/\n"
        elif self.run_type == "vc-md":
            filestring += "&IONS\n"
            filestring += "  tempw="+str(self.temperature)+", \n"
            filestring += "  ion_temperature = 'rescale-v' ,\n"
            #filestring += "  tolp = 100 ,\n"
            #filestring += "  nraise = 1 ,\n"
            #filestring += "  ion_dynamics = 'beeman' ,\n"
            #filestring += "  pot_extrapolation = 'second-order' ,\n"
            #filestring += "  wfc_extrapolation = 'second-order' ,\n"
            filestring += "/\n"
            #
            filestring += "&CELL\n"
            #filestring += "  cell_dynamics='pr',\n"
            filestring += "  cell_dynamics='w',\n"
            filestring += "  press_conv_thr = 0.5,\n" # Kbar
            filestring += "/\n"
        elif self.pressure:
            filestring += "&IONS\n"
            filestring += "  ion_dynamics = 'bfgs' ,\n"
            filestring += "/\n"
            #
            filestring += "&CELL\n"
            filestring += "  cell_dynamics = 'bfgs' ,\n"
            filestring += "  press = %7.4f ,\n"%float(self.pressure)
            filestring += "  press_conv_thr = 0.5 ,\n" # Kbar
            filestring += "/\n"
        #
        if self.brav != True:
            if self.cartesianlatvects:
                if self.unit == "bohr":
                    filestring += "CELL_PARAMETERS {bohr}\n"
                elif self.unit == "angstrom":
                    filestring += "CELL_PARAMETERS {angstrom}\n"
                t = LatticeMatrix(self.cell.latticevectors)
                for i in range(3):
                    for j in range(3):
                        t[i][j] = self.cell.latticevectors[i][j]*self.cell.lengthscale
                filestring += str(t)
            else:
                filestring += "CELL_PARAMETERS {alat}\n"
                filestring += str(self.cell.latticevectors)
        #
        filestring += "ATOMIC_SPECIES\n"
        filestring += ATOMIC_SPECIES_string
        filestring = filestring[:-1]+" \n"
        #
        if str(self.run_type) == "neb-start":
            filestring += "BEGIN_POSITIONS\n"
            filestring += "FIRST_IMAGE\n"
        if str(self.run_type) == "neb-end":
            filestring  = "LAST_IMAGE\n"
        if self.cartesianpositions:
            if self.scaledcartesianpositions:
                filestring += "ATOMIC_POSITIONS {alat}\n"
            else:
                if self.unit == "bohr":
                    filestring += "ATOMIC_POSITIONS {bohr}\n"
                elif self.unit == "angstrom":
                    filestring += "ATOMIC_POSITIONS {angstrom}\n"
        else:
            if self.scaledcartesianpositions:
                filestring += "ATOMIC_POSITIONS {alat}\n"
            else:
                filestring += "ATOMIC_POSITIONS {crystal}\n"
        if self.fix_all_pos:
            if_pos = " 0 0 0 "
        else:
            if_pos = ""
        for a in self.cell.atomdata:
            for b in a:
                for mv_atom_sp in self.move_atom_sp_list:
                    if b.spcstring() == mv_atom_sp:
                        if_pos = " 1 1 1 "
                        break
                    elif self.fix_all_pos:
                        if_pos = " 0 0 0 "
                for mv_atom_sp in self.fix_atom_sp_list:
                    if b.spcstring() == mv_atom_sp:
                        if_pos = " 0 0 0 "
                        break
                if self.cartesianpositions:
                    t = Vector(mvmult3(self.cell.latticevectors,b.position))
                    if self.scaledcartesianpositions:
                        filestring += "  "+b.spcstring().rjust(width)+" "+str(t)+if_pos+"\n"
                    else:
                        for i in range(3):
                            t[i] = self.cell.lengthscale*t[i]
                        filestring += "  "+b.spcstring().rjust(width)+" "+str(t)+if_pos+"\n"
                else:
                    if self.scaledcartesianpositions:
                        t = Vector(mvmult3(self.cell.latticevectors,b.position))
                        filestring += "  "+b.spcstring().rjust(width)+" "+str(t)+if_pos+"\n"
                    else:
                        filestring += "  "+b.spcstring().rjust(width)+" "+str(b.position)+if_pos+"\n"
        #
        if str(self.run_type) == "neb-start":
            neb_start_filestring = filestring
        if str(self.run_type) == "neb-end":
            filestring  += "END_POSITIONS\n"
        # Opt for gamma-point run if possible
        if self.run_type != "bands":
            if self.kgrid[0]*self.kgrid[1]*self.kgrid[2] == 1:                
                filestring += "K_POINTS gamma\n"
            else:
                filestring += "K_POINTS automatic\n"
                #filestring += "  "+str(self.kgrid[0])+" "+str(self.kgrid[1])+" "+str(self.kgrid[2])+"  1 1 1\n"
                if self.run_type == "dos":
                    times = 2
                    #
                    f = open("dos.in", "w")
                    dos_filestring  = "&dos\n"
                    dos_filestring += "    outdir = './work/',\n"
                    dos_filestring += "    prefix  = '"+self.filename+"' ,\n"
                    dos_filestring += "    fildos = '"+self.filename+".dos.data',\n"
                    dos_filestring += "/\n"
                    f.write(str(dos_filestring))
                    f.close()
                    #
                    f = open("dos.plot", "w")
                    dos_plot_filestring  = "#!/usr/local/bin/gnuplot -persist\n"
                    dos_plot_filestring += "# Last modified: 2018/12/17\n"
                    dos_plot_filestring += "set terminal postscript eps enhanced 28 lw 2\n"
                    dos_plot_filestring += 'set output "'+self.filename+'.dos.eps"\n'
                    dos_plot_filestring += "unset key\n"
                    dos_plot_filestring += "set xlabel 'Energy (eV)'\n"
                    dos_plot_filestring += "set ylabel 'DOS (states/eV)'\n"
                    dos_plot_filestring += 'dos_filename = "'+self.filename+'.dos.data"\n'
                    dos_plot_filestring += 'ef = system("cat " . dos_filename . " | grep " . "EFermi" . "'+" | awk '{print $9}'"+'")\n'
                    dos_plot_filestring += "set yzeroaxis lt 2\n"
                    dos_plot_filestring += "plot dos_filename using ($1-ef):2 w l\n"
                    f.write(str(dos_plot_filestring))
                    f.close()
                    #
                else:
                    times = 1
                #
                kshift_list = []
                if self.kpeven == "yes":
                    if self.cell.crystal_system() == "hexagonal":
                        kshift_list = [2]
                    elif self.cell.crystal_system() == "tetragonal":
                        kshift_list = []
                    elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                        kshift_list = [0,1,2]
                    elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                        kshift_list = [0,1,2]
                    elif self.cell.spacegroupsetting == "P":
                        kshift_list = [0,1,2]
                    else:
                        kshift_list = []
                    #
                    for k_no in kshift_list:
                        if self.kgrid[k_no]%2 != 0 and self.kgrid[k_no] > 1:
                            self.kgrid[k_no] = int(self.kgrid[k_no]-0.1)
                        if self.qgrid[k_no] >= self.kgrid[k_no]:
                            self.qgrid[k_no] = self.kgrid[k_no]
                        elif self.qgrid[k_no]%2 != 0:
                            self.qgrid[k_no] = self.kgrid[k_no]/2
                #
                kshift = [0,0,0]
                for k_no in kshift_list:
                    if self.kgrid[k_no] == 1:
                        kshift[k_no] = 0
                    else:
                        kshift[k_no] = 1
                if self.cell.crystal_system() == "hexagonal":
                    filestring += "  "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"  0 0 "+str(kshift[2])+"\n"
                elif self.cell.crystal_system() == "tetragonal":
                    filestring += "  "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"  0 0 0\n"
                elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                    filestring += "  "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"  "+str(kshift[0])+" "+str(kshift[1])+" "+str(kshift[2])+"\n"
                elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                    filestring += "  "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"  "+str(kshift[0])+" "+str(kshift[1])+" "+str(kshift[2])+"\n"
                elif self.cell.spacegroupsetting == "P":
                    filestring += "  "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"  "+str(kshift[0])+" "+str(kshift[1])+" "+str(kshift[2])+"\n"
                else:
                    filestring += "  "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"  0 0 0\n"
        else:
            filestring += "K_POINTS crystal_b\n"
            #filestring += "  "+str(self.kgrid[0])+" "+str(self.kgrid[1])+" "+str(self.kgrid[2])+"  1 1 1\n"
            times = 1
            if self.ibrav()==1: #sc
                nk_symbol_ibrav = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_list         = ["20","40","40","20","1 "]
                #
                nk_symbol       = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_pos          = ["0.5  0.5  0.5","0.0  0.0  0.0","0.0  0.5  0.0","0.5  0.5  0.0","0.0  0.0  0.0"] 
            elif self.ibrav()==2: # fcc (cubic)
                nk_symbol_ibrav = ["W ","L ","gG","X ","W ","K "] # Modules/bz_form.f90 init_bz_2 ibz=2
                nk_list         = ["20","40","40","20","10","1 "]
                #
                nk_symbol       = ["W ","L ","gG","X ","W ","K "] # Modules/bz_form.f90 init_bz_2 ibz=2
                nk_pos          = ["0.5   0.25  0.75","0.0   0.5   0.0 ","0.0   0.0   0.0 ","0.0   0.5   0.5 ","0.75  0.5   0.25","0.75  0.375 0.375"]
            elif self.ibrav()==3: # bcc (cubic) 
                nk_symbol_ibrav = ["gG","H ","N ","gG","P "]  # Modules/bz_form.f90 init_bz_3 ibz=3
                nk_list         = ["20","40","40","20","1 "]
                #
                nk_symbol       = ["gG","H ","N ","gG","P "]  # Modules/bz_form.f90 init_bz_3 ibz=3
                nk_pos          = ["0.0  0.0  0.0 ","0.5 -0.5  0.5 ","0.0  0.0  0.5 ","0.0  0.0  0.0 ","0.25 0.25 0.25"]
            elif self.ibrav()==4: # hcp 
                nk_symbol_ibrav = ["gG","M ","K ","gG","A "]  # init_bz_13 hexagonal ibz=13
                nk_list         = ["20","40","40","20","1 "]
                #
                nk_symbol       = ["gG","M ","K ","gG","A "]  # init_bz_13 hexagonal ibz=13
                nk_pos          = ["0.0   0.0   0.0","0.5   0.0   0.0","0.333 0.333 0.0","0.0   0.0   0.0","0.0   0.0   0.5"]
            elif self.ibrav()==5 and self.cell.alpha < 90: # init_bz_14 trigonal alpha < 90 bz ibz=14
                nk_symbol_ibrav = ["gG","L ","X ","L1","P2","F ","P1","Z ","B ","Q ","B1"]
                nk_list         = ["10","10","10","10","10","10","10","10","10","10","1 "]
                # 
                nk_symbol       = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_pos          = ["0.5  0.5  0.5","0.0  0.0  0.0","0.0  0.5  0.0","0.5  0.5  0.0","0.0  0.0  0.0"] 
            elif self.ibrav()==5 and 90 < self.cell.alpha and self.cell.alpha < 120: # init_bz_15 rigonal alpha > 90 bz ibz=15
                nk_symbol_ibrav = ["gG","P ","F ","P1","Q1","Z ","L ","Q ","gG"]
                nk_list         = ["10","10","10","10","10","10","10","10","1 "]
                #
                nk_symbol       = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_pos          = ["0.5  0.5  0.5","0.0  0.0  0.0","0.0  0.5  0.0","0.5  0.5  0.0","0.0  0.0  0.0"] 
            elif self.ibrav()==6: # Modules/bz_form.f90 init_bz_4 simple tetragonal bz ibz=4
                nk_symbol_ibrav = ["gG","X ","R ","Z ","gG","M ","A ","Z "]
                nk_list         = ["10","10","10","10","10","10","10","1 "]
                #
                nk_symbol       = ["gG","X ","R ","Z ","gG","M ","A ","Z "]
                nk_pos          = ["0.0   0.0   0.0","0.0   0.5   0.0","0.0   0.5   0.5","0.0   0.0   0.5","0.0   0.0   0.0","0.5   0.5   0.0","0.5   0.5   0.5","0.0   0.0   0.5"]
            elif self.ibrav()==7 and (self.cell.c/self.cell.a) < 1: # init_bz_5 centered tetragonal (c<a) bz ibz=5
                nk_symbol_ibrav = ["gG","X ","P ","N ","Z ","gG","M ","Z1","N"]
                nk_list         = ["10","10","10","10","10","10","10","10","1 "]
                #
                nk_symbol       = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_pos          = ["0.5  0.5  0.5","0.0  0.0  0.0","0.0  0.5  0.0","0.5  0.5  0.0","0.0  0.0  0.0"] 
            elif self.ibrav()==7 and (self.cell.c/self.cell.a) > 1: # init_bz_6 centered tetragonal (c>a) bz ibz=6
                nk_symbol_ibrav = ["gG","gS","Y ","X ","P ","N ","gS1","Z ","Y1"]
                nk_list         = ["10","10","10","10","10","10","10","10","1 "]
                #
                nk_symbol       = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_pos          = ["0.5  0.5  0.5","0.0  0.0  0.0","0.0  0.5  0.0","0.5  0.5  0.0","0.0  0.0  0.0"] 
            elif self.ibrav()==8: # init_bz_7 simple orthorombic bz ibz=7
                nk_symbol_ibrav = ["gG","X ","U ","Z ","gG","Y ","S ","R ","T ","Y "]
                nk_list         = ["10","10","10","10","10","10","10","10","10","1 "]
                #
                nk_symbol       = ["gG","X ","U ","Z ","gG","Y ","S ","R ","T ","Y "]
                nk_pos          = ["0.0   0.0   0.0","0.5   0.0   0.0","0.5   0.0   0.5","0.0   0.0   0.5","0.0   0.0   0.0","0.0   0.5   0.0","0.5   0.5   0.0","0.5   0.5   0.5","0.0   0.5   0.5","0.0   0.5   0.0"]
            elif self.ibrav()==9 and (self.cell.b/self.cell.a) < 1: # init_bz_12 ibz=12
                nk_symbol_ibrav = ["gG","X ","Y1","S ","Y ","gG","Z ","A ","A1","R ","T ","Z "]
                nk_list         = ["10","10","10","10","10","10","10","10","10","10","10","1 "]
                #
                nk_symbol       = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_pos          = ["0.5  0.5  0.5","0.0  0.0  0.0","0.0  0.5  0.0","0.5  0.5  0.0","0.0  0.0  0.0"] 
            elif self.ibrav()==9 and (self.cell.b/self.cell.a) > 1: # init_bz_12 ibz=12
                nk_symbol_ibrav = ["gG","X ","S ","X1","Y ","gG","Z ","A ","R ","A1","T ","Z "]
                nk_list         = ["10","10","10","10","10","10","10","10","10","10","10","1 "]
                #
                nk_symbol       = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_pos          = ["0.5  0.5  0.5","0.0  0.0  0.0","0.0  0.5  0.0","0.5  0.5  0.0","0.0  0.0  0.0"] 
            elif self.ibrav()==10: #init_bz_8,9,10 ibz=8,9,10
                nk_symbol_ibrav = ["L ","gG","X ","Y ","gG","Z "]
                nk_list         = ["10","10","10","10","10","1 "]
                #
            elif self.ibrav()==11: # init_bz_11 ibz=11
                nk_symbol_ibrav = ["L ","gG","X ","Y ","gG","Z ","R ","W ","T ","gG","S "]
                nk_list         = ["10","10","10","10","10","10","10","10","10","10","1 "]
                #
                nk_symbol       = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_pos          = ["0.5  0.5  0.5","0.0  0.0  0.0","0.0  0.5  0.0","0.5  0.5  0.0","0.0  0.0  0.0"] 
            elif self.ibrav()==12: # init_bz_16 Simple monoclinic lattice c unique ibz=16
                nk_symbol_ibrav = ["gG","X ","R ","gG","Y ","D ","Z "]
                nk_list         = ["10","10","10","10","10","10","1 "]
                #
                nk_symbol       = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_pos          = ["0.5  0.5  0.5","0.0  0.0  0.0","0.0  0.5  0.0","0.5  0.5  0.0","0.0  0.0  0.0"] 
            elif self.ibrav()==-12: # Simple monoclinic lattice b unique
                nk_symbol_ibrav = ["gG","X ","A ","Y ","gG","Z ","D ","Y "]
                nk_list         = ["10","10","10","10","10","10","10","1 "]
                #
                nk_symbol       = ["R ","gG","X ","M ","gG"] # Modules/bz_form.f90 init_bz_1 ibz=1
                nk_pos          = ["0.5  0.5  0.5","0.0  0.0  0.0","0.0  0.5  0.0","0.5  0.5  0.0","0.0  0.0  0.0"] 
            elif self.ibrav()==13: # not support
                print 'find_bz_type','This ibrav is not supported'
            elif self.ibrav()==14: # not support
                print 'find_bz_type','This ibrav is not supported'
            else:
                print 'find_bz_type','Wrong ibrav'
            #
            filestring += "  "+str(len(nk_list))+"\n"
            bands_plot_filestring_symmbol = 'set xtics ('
            if self.brav: 
                nk = len(nk_symbol_ibrav)
            else:
                nk = len(nk_symbol)
            for nk_no in range(nk):
                if self.brav: 
                    filestring += "  "+nk_symbol[nk_no]+"  "+nk_list[nk_no]+"\n"
                    now_nk_symbol = nk_symbol_ibrav[nk_no]
                else:
                    filestring += "  "+nk_pos[nk_no]+"  "+nk_list[nk_no]+"\n"
                    now_nk_symbol = nk_symbol[nk_no]
                plot_symbol = now_nk_symbol.replace(" ","").replace("gG","{/Symbol G}").replace("gS","{/Symbol S}").replace("1","_1").replace("2","_2")
                bands_plot_filestring_symmbol += '"'+plot_symbol+'" x'+str(nk_no)+', '
            bands_plot_filestring_symmbol = bands_plot_filestring_symmbol[:-2]
            bands_plot_filestring_symmbol += ")\n"
            #
            f = open("bands.in", "w")
            bands_filestring  = "&bands\n"
            bands_filestring += "   outdir = './work/',\n"
            bands_filestring += "   prefix  = '"+self.filename+"' ,\n"
            bands_filestring += "   filband = '"+self.filename+".bands' ,\n"
            bands_filestring += "   lsym=.true.\n"
            bands_filestring += "/\n"
            f.write(str(bands_filestring))
            f.close()
            #
            f = open("bands.plot", "w")
            bands_plot_filestring  = "#!/usr/local/bin/gnuplot -persist\n"
            bands_plot_filestring += "# Last modified: 2018/12/17\n"
            bands_plot_filestring += "set terminal postscript eps enhanced 28 lw 2\n"
            bands_plot_filestring += 'set output "'+self.filename+'.bands.eps"\n'
            bands_plot_filestring += "\n"
            bands_plot_filestring += "#set size 0.45,1\n"
            bands_plot_filestring += "\n"
            bands_plot_filestring += "set ylabel 'Energy (eV)'\n"
            bands_plot_filestring += "\n"
            bands_plot_filestring += "unset key\n"
            bands_plot_filestring += "\n"
            bands_plot_filestring += "#ry = 13.60569193\n"
            bands_plot_filestring += "ry = 1\n"
            bands_plot_filestring += 'filename = "'+self.filename+'.bands.gnu"\n'
            bands_plot_filestring += "x0 = 0\n"
            xno = 0
            nkp = 1
            for nk in range(len(nk_list)-1):
                nkp += int(nk_list[xno])
                xno += 1
                bands_plot_filestring += 'x'+str(xno)+' = system("cat " . filename . " | awk \\'+'\'NR==" . '+str(nkp)+' ."{printf(\\"%8.4f\\", $" . 1 . ")}\\\'")\n'
            bands_plot_filestring += "ymin = -15.0\n"
            bands_plot_filestring += "ymax =   5.0\n"
            if os.path.isfile(self.filename+".dos.data"):
                bands_plot_filestring += 'dos_filename = "'+self.filename+'.dos.data"\n'
                bands_plot_filestring += 'ef = system("cat " . dos_filename . " | grep " . "EFermi" . "'+" | awk '{print $9}'"+'")\n'
            elif os.path.isfile(self.filename+".scf.out"):
                bands_plot_filestring += 'scf_filename = "'+self.filename+'.scf.out"\n'
                bands_plot_filestring += 'ef = system("cat " . scf_filename . " | grep " . "Fermi" . "'+" | awk '{print $5}'"+'")\n'
            else:
                bands_plot_filestring += "ef = 0.0\n"
            bands_plot_filestring += "\n"
            bands_plot_filestring += "set xrange [0:x"+str(len(nk_list)-1)+"]\n"
            bands_plot_filestring += "set yrange [ymin:ymax]\n"
            bands_plot_filestring += "\n"
            bands_plot_filestring += bands_plot_filestring_symmbol
            bands_plot_filestring += "set xzeroaxis lt 2\n"
            xno = 0
            for nk in range(len(nk_list)-1):
                xno += 1
                bands_plot_filestring += "set arrow "+str(xno)+" nohead from x"+str(xno)+",ymin to x"+str(xno)+",ymax lt 1\n"
            bands_plot_filestring += "\n"
            bands_plot_filestring += "plot filename using 1:($2*ry-ef) w l\n"
            f.write(str(bands_plot_filestring))
            f.close()
        #
        if str(self.run_type) == "neb-end":
            filestring += "END_ENGINE_INPUT\n"
            filestring += "END\n"
            neb_end_filestring = filestring
        #
        if self.run_type == "tddft" or self.run_type == "eels":
            if self.run_type == "tddft":
                f = open(self.filename+".tddfpt.in", "w")
            else:
                f = open(self.filename+".tddfpt-eels.in", "w")
            tddfpt_filestring  = "&lr_input\n"
            tddfpt_filestring += "  outdir = './work/',\n"
            tddfpt_filestring += "  prefix  = '"+self.filename+"' ,\n"
            tddfpt_filestring += "  restart_step=250,\n"
            tddfpt_filestring += "  restart=.false.\n"
            tddfpt_filestring += "/\n"
            tddfpt_filestring += "&lr_control\n"
            tddfpt_filestring += "  itermax=500,\n"
            if self.run_type == "tddft":
                tddfpt_filestring  += "  ipol=4\n"
            else:
                tddfpt_filestring  += "  q1 = 0.0d0,\n"
                tddfpt_filestring  += "  q2 = 0.0d0,\n"
                tddfpt_filestring  += "  q3 = 0.15d0,\n"
            tddfpt_filestring += "/\n"
            f.write(str(tddfpt_filestring))
            f.close()
            #
            f = open(self.filename+".tddfpt_pp.in", "w")
            tddfpt_pp_filestring  = "&lr_input\n"
            tddfpt_pp_filestring += "  outdir = './work/',\n"
            tddfpt_pp_filestring += "  prefix  = '"+self.filename+"' ,\n"
            tddfpt_pp_filestring += "  itermax0 = 500\n"
            tddfpt_pp_filestring += "  itermax  = 10000\n"
            tddfpt_pp_filestring += '  extrapolation="osc"\n'
            tddfpt_pp_filestring += "  units=1\n"
            tddfpt_pp_filestring += "  start=0.0d0\n"
            tddfpt_pp_filestring += "  increment=0.001d0\n"
            if self.run_type == "tddft":
                tddfpt_pp_filestring += "  end=3.50d0\n"
                tddfpt_pp_filestring += "  epsil=0.01\n"
                tddfpt_pp_filestring += "  ipol=4\n"
            else:
                tddfpt_pp_filestring += "  end=50.d0\n"
                tddfpt_pp_filestring += "  epsil=0.035\n"
            tddfpt_pp_filestring += "/\n"
            f.write(str(tddfpt_pp_filestring))
            f.close()
        #
        if self.run_type == "phonon":
            f = open("ph.in", "w")
            phonon_filestring  = "&inputph\n"
            phonon_filestring += "   alpha_mix(1) = 0.7,\n"
            phonon_filestring += "   niter_ph = 99,\n"
            phonon_filestring += "   tr2_ph = 1.0d-12,\n"
            phonon_filestring += "   prefix  = '"+self.filename+"' ,\n"
            phonon_filestring += "   outdir = './work/',\n"
            phonon_filestring += "   fildyn = '"+self.filename+".dyn',\n"
            phonon_filestring += "   lsym = .true.\n"
            phonon_filestring += "   nq1 = "+str(self.qgrid[0])+",\n"
            phonon_filestring += "   nq2 = "+str(self.qgrid[1])+",\n"
            phonon_filestring += "   nq3 = "+str(self.qgrid[2])+",\n"
            phonon_filestring += "/\n"
            f.write(str(phonon_filestring))
            f.close()
            #
            f = open("q2r.in", "w")
            q2r_filestring  = "&input\n"
            q2r_filestring += "   fildyn = 'ph.dyn',\n"
            q2r_filestring += "   zasr = 'simple',\n"
            q2r_filestring += "   flfrc = 'ph.fc',\n"
            q2r_filestring += "/\n"
            f.write(str(q2r_filestring))
            f.close()
            #
            f = open("Plot_input", "w")
            Plot_input_filestring  = "#!/bin/bash\n"
            Plot_input_filestring += ". ../environment_variables\n"
            Plot_input_filestring += "FC_name='"+self.filename+"'\n"
            Plot_input_filestring += "cat > Atomic_mass <<EOF'\n"
            ph_natom = 0
            for sp in self.species:
                ph_natom += 1
                Plot_input_filestring += ("    amass("+str(ph_natom)+") = %8.5f,\n"%(ed.elementweight[sp])).rjust(11)  
            Plot_input_filestring += "EOF\n"
            Plot_input_filestring += "freq=meV\n"
            f.write(str(Plot_input_filestring))
            f.close()
            #
            f = open("matdyn.dos.in", "w")
            times = 3
            matdyn_dos_filestring  = "&input\n"
            matdyn_dos_filestring += "   asr = 'simple',\n"
            matdyn_dos_filestring += "   flfrc = '"+self.filename+".fc'\n"
            matdyn_dos_filestring += "   dos = .ture.,\n"
            matdyn_dos_filestring += "   fldos='"+self.filename+".phdos',\n"
            matdyn_dos_filestring += "   nq1 = "+str(self.qgrid[0]*times)+",\n"
            matdyn_dos_filestring += "   nq2 = "+str(self.qgrid[1]*times)+",\n"
            matdyn_dos_filestring += "   nq3 = "+str(self.qgrid[2]*times)+",\n"
            matdyn_dos_filestring += "/\n"
            f.write(str(matdyn_dos_filestring))
            f.close()
            #
            f = open("matdyn.bands.in", "w")
            matdyn_band_filestring  = "&input\n"
            matdyn_band_filestring += "   asr = 'simple',\n"
            matdyn_band_filestring += "   flfrc = '"+self.filename+".fc'\n"
            matdyn_band_filestring += "   flfrq = '"+self.filename+".freq',\n"
            matdyn_band_filestring += "   q_in_band_form = .true.,',\n"
            matdyn_band_filestring += "/\n"
            matdyn_band_filestring += "2\n"
            matdyn_band_filestring += "0.0 0.0 0.0 5\n"
            matdyn_band_filestring += "0.5 0.0 0.0 1\n"
            f.write(str(matdyn_band_filestring))
            f.close()
        #
        if self.run_type == "lead" or self.run_type == "scat" or self.run_type == "lead_left" or self.run_type == "lead_right":
            if self.run_type == "lead_left":
                f = open(self.filename+".cond_left.in", "w")
            elif self.run_type == "lead_right":
                f = open(self.filename+".cond_right.in", "w")
            else:
                f = open(self.filename+".cond.in", "w")
            cond_filestring  = " &inputcond\n"
            cond_filestring += "  outdir  = './work/' ,\n"
            if self.run_type == "lead_left":
                cond_filestring += "    ikind = 2\n"
                cond_filestring += "    prefixl='lead_left',\n"
                cond_filestring += "    prefixs='scat',\n"
                cond_filestring += "    tran_file='trans."+self.filename+"_left.out',\n"
            elif self.run_type == "lead_right":
                cond_filestring += "    ikind = 2\n"
                cond_filestring += "    prefixl='lead_left',\n"
                cond_filestring += "    prefixs='scat',\n"
                cond_filestring += "    prefixr='lead_right',\n"
                cond_filestring += "    tran_file='trans."+self.filename+"_right.out',\n"
            elif self.run_type == "lead":
                cond_filestring += "    ikind = 0\n"
                cond_filestring += "    prefixl='lead',\n"
                if self.spin == "yes":
                    cond_filestring += "    band_file='bands_up',\n"
                else:
                    cond_filestring += "    band_file='bands',\n"
                    #cond_filestring += "    prefixs='lead',\n"
                cond_filestring += "    tran_file='trans."+self.filename+".out_up',\n"
            else:
                cond_filestring += "    ikind = 1\n"
                cond_filestring += "    prefixl='lead',\n"
                cond_filestring += "    prefixs='scat',\n"
                cond_filestring += "    tran_file='trans."+self.filename+".out',\n"
            if self.spin == "yes":
                cond_filestring += "    iofspin = 2,\n"
            else:
                cond_filestring += "    iofspin = 1,\n"
            cond_filestring += "    energy0 =  0.01\n"
            cond_filestring += "    denergy = -0.01\n"
            cond_filestring += "    ewind = 3.d0\n"
            cond_filestring += "    epsproj = 1.d-5\n"
            cond_filestring += "    nz1 = 11\n"
            cond_filestring += "    delgep = 1.d-10\n"
            cond_filestring += " /\n"
            cond_filestring += "  1\n"
            cond_filestring += "  0.0 0.0 0.1250000\n"
            cond_filestring += "  11\n"
            cond_filestring += "\n"
            f.write(str(cond_filestring))
            f.close()
        #
        filestring += "\n"
        if str(self.run_type) == "neb-start":
            filestring = neb_start_filestring
        if str(self.run_type) == "neb-end":
            filestring = neb_end_filestring
        #
        return filestring
    #
    # Return the PWscf internal bravais lattice number
    def ibrav(self):
        system = self.cell.crystal_system()
        setting = self.cell.spacegroupsetting
        #if self.cell.supercell:
        #    return 14
        if system == 'cubic':
            if self.cell.primcell:
                if setting == 'P':
                    return 1 # cubic P (sc)
                elif setting == 'F':
                    return 2 # cubic F (fcc)
                elif setting == 'I':
                    if self.cell.latticevectors[0][0] >= 0:
                        return 3 # cubic I (bcc)
                        # v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1)
                    else:
                        return -3  #cubic I (bcc), more symmetric axis:
                        # v1 = (a/2)(-1,1,1), v2 = (a/2)(1,-1,1),  v3 = (a/2)(1,1,-1)
            else:
                return 1
        if system == 'hexagonal' or system == 'trigonal':
            if self.primcell:
                if setting == 'P':
                    return 4 # Hexagonal and Trigonal P
                elif setting == 'R':
                    if self.cell.latticevectors[1][0] == 0:
                        return 5 # Trigonal R, 3fold axis c
                        # v1 = a(tx,-ty,tz),   v2 = a(0,2ty,tz),   v3 = a(-tx,-ty,tz)
                    else:
                        return -5 # Trigonal R, 3fold axis <111>
                        # v1 = a' (u,v,v),   v2 = a' (v,u,v),   v3 = a' (v,v,u)
        #
        if system == 'tetragonal':
            if self.primcell:
                if setting == 'P':
                    return 6 # Tetragonal P (st)
                elif setting == 'I':
                    return 7 # Tetragonal I (bct) 
        if system == 'orthorhombic':
            if self.primcell:
                if setting == 'P':
                    return 8 # Orthorhombic P 
                elif setting == 'C':
                    if self.cell.latticevectors[0][1] >= 0:
                        return 9 # Orthorhombic base-centered(bco)
                        # v1 = (a/2, b/2,0),  v2 = (-a/2,b/2,0),  v3 = (0,0,c)
                    else:
                        return -9 # as 9, alternate description
                        # v1 = (a/2,-b/2,0),  v2 = (a/2, b/2,0),  v3 = (0,0,c)
                elif setting == 'A': 
                    return 91 # Orthorhombic one-face base-centered A-type
                elif setting == 'F':
                    return 10 # Orthorhombic face-centered
                elif setting == 'I':
                    return 11 # Orthorhombic body-centered
        if system == 'monoclinic':
            if self.primcell:
                if self.cell.gamma != 90: # c unique
                    return 12 # Monoclinic P, unique axis c 
                elif self.cell.beta != 90: # b unique
                    return -12 # Monoclinic P, unique axis b
            else:
                if self.cell.gamma != 90: # c unique
                    return 13 # Monoclinic base-centered (unique axis c)
                elif self.cell.beta != 90: # b unique
                    return -13 # Monoclinic base-centered (unique axis b)
        if system == 'triclinic':
            return 14
    #

################################################################################################
# CP2K
class CP2KFile(GeometryOutputFile):
    """
    Class for storing the geometrical data for a CP2k run and the method
    __str__ that outputs to a .inp file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        self.docstring += "\n"
    def __str__(self):
        filestring = self.docstring
        filestring += "&CELL\n"
        filestring += "  PERIODIC XYZ\n"
        filestring += "  A "+str(self.cell.latticevectors[0].scalmult(self.cell.lengthscale))+"\n"
        filestring += "  B "+str(self.cell.latticevectors[1].scalmult(self.cell.lengthscale))+"\n"
        filestring += "  C "+str(self.cell.latticevectors[2].scalmult(self.cell.lengthscale))+"\n"
        filestring += "&END CELL\n\n"
        filestring += "&COORD\n"
        for a in self.cell.atomdata:
            for b in a:
                filestring += b.spcstring()+str(Vector(mvmult3(self.cell.latticevectors,b.position)).scalmult(self.cell.lengthscale))+"\n"
        filestring += "&END COORD\n"
        return filestring 

################################################################################################
class CPMDFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a CPMD run and the method
    __str__ that outputs to a .inp file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("bohr")
        self.cutoff = 100.0
    def __str__(self):
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # Transformation to cartesian coordinates
        transmtx = []
        for i in range(3):
            transmtx.append([])
            for j in range(3):
                transmtx[i].append(lattice[i][j] * a)
        # docstring
        filestring = self.docstring+"\n"
        filestring += "&SYSTEM\n"
        # lattice
        filestring += " CELL VECTORS\n"
        for vec in transmtx:
            for coord in vec:
                filestring += " %19.15f"%coord
            filestring += "\n"
        # Cutoff
        filestring += " CUTOFF\n"
        filestring += " "+str(self.cutoff)+"\n"
        filestring += "&END\n\n"
        # The atom position info
        filestring += "&ATOMS\n"
        # get all species
        species = set([])
        for a in self.cell.atomdata:
            for b in a:
                species.add(b.spcstring())
        for spc in species:
            filestring += "*[pseudopotential file for "+spc+" here]\n"
            # Find maximal angular momentum
            spcs = spc.split("/")
            l = "s"
            for s in spcs:
                 if ed.angularmomentum[ed.elementblock[s]] > ed.angularmomentum[l]:
                     l = ed.elementblock[s]
            natoms = 0
            posstring = ""
            for a in self.cell.atomdata:
                for b in a:
                    if b.spcstring() == spc:
                        natoms +=1
                        posstring += str(Vector(mvmult3(transmtx,b.position)))+"\n"
            # Print
            filestring += str(natoms)+"\n"
            filestring += posstring
        filestring += "&END\n"
        return filestring
################################################################################################
# DFTB
class DFTBFile(GeometryOutputFile):
    """
    Class for storing the geometrical data for a DFTB run and the method
    __str__ that outputs to a .gen file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
    def __str__(self):
        natoms = sum([len(a) for a in self.cell.atomdata])
        filestring = self.docstring
        # Assumes the user wants cartesian coordinates ('S'), not fractional ('F')
        filestring += str(natoms) + " S\n"
        species = set([])
        natom = 0
        for a in self.cell.atomdata:
            natom += len(a)
            for b in a:
                species.add(b.spcstring())
        filestring += ' '.join(species) + "\n"
        # Numbered dictionary from entry 1
        nameDict = dict(zip(species, range(1, len(species) + 1)))
        nAt = 0
        for a in self.cell.atomdata:
            for b in a:
                nAt += 1
                filestring += str(nAt)+" "+str(nameDict[b.spcstring()])+str(Vector(mvmult3(self.cell.latticevectors,b.position)).scalmult(self.cell.lengthscale))+"\n"
        filestring += "  0.0 0.0 0.0\n"
        filestring += str(self.cell.latticevectors[0].scalmult(self.cell.lengthscale))+"\n"
        filestring += str(self.cell.latticevectors[1].scalmult(self.cell.lengthscale))+"\n"
        filestring += str(self.cell.latticevectors[2].scalmult(self.cell.lengthscale))+"\n"
        return filestring

################################################################################################
class SiestaFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a Siesta run and the method
    __str__ that outputs to a .fdf file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
    def __str__(self):
        # Assign some local variables
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        filestring = self.docstring
        filestring += "AtomicCoordinatesFormat".ljust(28)+"Fractional\n"
        species = set([])
        natom = 0
        for a in self.cell.atomdata:
            natom += len(a)
            for b in a:
                species.add(b.spcstring())
        species = list(species)
        nspcs = len(species)
        filestring += "LatticeConstant".ljust(28)+str(self.cell.lengthscale)+" Ang\n"
        filestring += "NumberOfAtoms".ljust(28)+str(natom)+"\n"
        filestring += "NumberOfSpecies".ljust(28)+str(nspcs)+"\n"
        # lattice
        filestring += "%block LatticeVectors\n"
        for vec in lattice:
            filestring += str(vec)+"\n"
        filestring += "%endblock LatticeVectors\n"
        # Atomic coordinates
        filestring += "%block AtomicCoordinatesAndAtomicSpecies\n"
        i = 1
        for sp in species:
            for a in self.cell.atomdata:
                for b in a:
                    if b.spcstring() == sp:
                        filestring += str(b.position)
                        filestring += "   %i\n"%i
            i += 1
        filestring += "%endblock AtomicCoordinatesAndAtomicSpecies\n"
        # Chemical species
        filestring += "%block ChemicalSpeciesLabel\n"
        i = 1
        for sp in species:
            filestring += str(i).ljust(8)
            if len(sp) > 2:
                filestring += "??      ??    # "
                tsp = sp.split("/")
                for t in tsp:
                    filestring += str(ed.elementnr[t])+"/"
                filestring = filestring.rstrip("/")
                filestring += "      "+sp+"\n"
            else:
                filestring += str(ed.elementnr[sp]).ljust(8)+sp.ljust(8)+"\n"
            i += 1
        filestring += "%endblock ChemicalSpeciesLabel\n"
        return filestring

################################################################################################
class ABINITFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an abinit run and the method
    __str__ that outputs the contents of a abinit input file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("bohr")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        self.printbraces = False
    def __str__(self):
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        filestring = self.docstring
        # VASP needs lattice vector matrix to have positive triple product
        if det3(lattice) < 0:
            if lattice[0].length() == lattice[1].length() == lattice[2].length():
                # Shift the first and last for cubic lattices
                transmtx = [[0, 0, 1],
                            [0, 1, 0],
                            [1, 0, 0]]
            else:
                # Else shift the two shortest
                if lattice[0].length() > lattice[1].length() and lattice[0].length() > lattice[2].length():
                    transmtx = [[1, 0, 0],
                                [0, 0, 1],
                                [0, 1, 0]]
                elif lattice[1].length() > lattice[2].length() and lattice[1].length() > lattice[0].length():
                    transmtx = [[0, 0, 1],
                                [0, 1, 0],
                                [1, 0, 0]]
                else:
                    transmtx = [[0, 1, 0],
                                [1, 0, 0],
                                [0, 0, 1]]
        else:
            transmtx = [[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
        lattice = mmmult3(transmtx,lattice)
        # Print braces around values (or not)
        if self.printbraces:
            lbrace = "{"
            rbrace = "}"
        else:
            lbrace = " "
            rbrace = ""
        # length scale and lattice
        filestring += "# Structural parameters\n"
        filestring += "acell "+lbrace+"  3*"+str(a)+" "+rbrace+" \n\n"
        filestring += "rprim "+lbrace
        for vec in lattice:
            filestring += str(Vector(vec))+"\n       "
        filestring = filestring[:-1]
        filestring += rbrace+" \n"
        # The atom position info
        alloy = False
        spcs = ""
        typatstring = "typat  "+lbrace+" "
        natom = 0
        ntypat = 0
        znuclstring = "znucl  "+lbrace+" "
        alloystring = ""
        xredstring = "xred "+lbrace+" "
        for a in self.cell.atomdata:
            for b in a:
                natom += 1
                if spcs != b.spcstring():
                    ntypat += 1
                    if b.alloy():
                        znuclstring += "?? "
                        alloystring += b.spcstring()+" "
                        alloy = True
                    else:
                        znuclstring += str(ed.elementnr[b.spcstring()])+" "
                typatstring += str(ntypat)+" "
                xredstring += str(Vector(mvmult3(transmtx,b.position)))+"\n       "
                spcs = b.spcstring()
        filestring += "natom  "+lbrace+" "+str(natom)+" "+rbrace+" \n"
        filestring += "ntypat "+lbrace+" "+str(ntypat)+" "+rbrace+" \n"
        filestring += typatstring+rbrace+" \n"
        filestring += znuclstring+rbrace+" "
        if alloy:
            filestring += "    # "+alloystring
        filestring += "\n"
        filestring += xredstring
        filestring = filestring[:-2]+rbrace+" \n"
        return filestring

class Abinit_Files_File(GeometryOutputFile):
    """
    Class for storing the input setting for a ABINIT run and the method
    __str__ that outputs to a .files file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Pseudopotential string
        #self.pseudostring = ".GGA_PBE-JTH.xml"
        # setup-all option
        self.setupall = False
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                tmp.add(b.spcstring())
        self.species = list(tmp)
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        self.filename = ""
    def __str__(self):
        # Determine max width of spcstring
        width = 0
        for a in self.cell.atomdata:
            for b in a:
                width = max(width, len(b.spcstring()))
        #
        filename = self.filename.replace(".in","")
        filestring  = filename+".in\n"
        sufix = ""
        tasks = self.run_type
        if tasks == "nonscf+dos":
            sufix = "_dos"
            filestring += filename+sufix+".out\n"
            filestring += filename+"o\n"
            filestring += filename+sufix+"\n"
        elif tasks == "nonscf+band":
            sufix = "_band"
            filestring += filename+sufix+".out\n"
            filestring += filename+"o\n"
            filestring += filename+sufix+"\n"
        else:
            filestring += filename+".out\n"
            filestring += filename+"i\n"
            filestring += filename+"o\n"
        filestring += filename+sufix+"tmp\n"
        #
        # POTCAR library
        directory = self.pseudolib
        if directory != "":
            self.dir = directory
        else:
            try:
                self.dir = os.environ['ABINIT_PSEUDOLIB']
            except:
                try:
                    self.dir = os.environ['ABINIT_PAWLIB']
                except:
                    self.dir = "../pseudo/"
        # check directory
        if self.dir == "":
            print "No path to the ABINIT pseudopotential library specified.\n"
        if not os.path.exists(self.dir):
            print "The specified path to the ABINIT pseudopotential library does not exist.\n"+self.dir
        #
        spcs = ""
        # get all files (GBRV, JTH)
        potlist = []
        potnamelist = []
        pp_filestring = ""
        for a in self.cell.atomdata:
            for b in a:
                if spcs != b.spcstring():
                    if b.alloy():
                        znuclstring += "?? "
                        alloystring += b.spcstring()+self.pseudostring+"\n"
                        alloy = True
                    else:
                        # get all files (GBRV, JTH)
                        potname = ""
                        exists_potname = ""
                        for version in self.prioritylist:
                            if self.pseudotype == "ONCVPSP":
                                potfile = self.dir+"/"+b.spcstring()+"."+str(self.pseudostring)
                                potname = b.spcstring()+"."+str(self.pseudostring)
                            elif self.pseudotype == "GBRV":
                                potfile = self.dir+"/"+b.spcstring().lower()+"_"+str(self.dfttype).lower()+"_"+version+str(self.pseudostring)
                                potname = b.spcstring().lower()+"_"+str(self.dfttype).lower()+"_"+version+str(self.pseudostring)
                            elif self.pseudotype == "JTH":
                                potfile = self.dir+"/"+b.spcstring()+"."+str(self.dfttype).upper()+str(self.pseudostring)
                                potname = b.spcstring()+"."+str(self.dfttype).upper()+str(self.pseudostring)
                            else:
                                potfile = self.dir+"/"+b.spcstring()+"."+str(self.dfttype).upper()+str(self.pseudostring)
                                potname = b.spcstring()+"."+str(self.dfttype).upper()+str(self.pseudostring)     
                            if os.path.exists(potfile):
                                exists_potname = potname
                                potlist.append(potfile)
                                potnamelist.append(potname)
                                break 
                        if self.dir == "" or not os.path.exists(self.dir) or exists_potname == "":
                            pp_filestring += b.spcstring()+"."+str(self.dfttype).upper()+str(self.pseudostring)  +"\n"
                            filestring += pp_filestring
                        else:
                            pp_filestring += "%s\n"%(exists_potname)
                            filestring += pp_filestring
                spcs = b.spcstring()
        #
        # read potcar files and put in outstring
        outstring = ""
        i = 0
        print ""
        print "setting pseudo-potential:"
        for f in potlist:
            print str(f)
            pot = open(f,"r")
            outstring = pot.read()
            pot.close()
            #
            abinit_pot = open(potnamelist[i], "w")
            abinit_pot.write(str(outstring))
            abinit_pot.close()
            #
            i += 1
        #
        tasks = self.run_type
        if tasks == "phonon":
            filename = self.filename.replace(".in","")
            f = open(self.filename.replace(".in",".kpt.files"), "w")
            kpt_filestring  = filename+".kpt.in\n"
            kpt_filestring += filename+".kpt.out\n"
            kpt_filestring += filename+".kpti\n"
            kpt_filestring += filename+".kpto\n"
            kpt_filestring += filename+".kpttmp\n"
            kpt_filestring += pp_filestring
            kpt_filestring += "\n"
            f.write(str(kpt_filestring))
            f.close()
            #
            filename = self.filename.replace(".in","")
            f = open(self.filename.replace(".in",".anaddb.files"), "w")
            anaddb_filestring  = filename+".anaddb.in\n"
            anaddb_filestring += filename+".anaddb.out\n"
            anaddb_filestring += filename+".mrgddb.out\n"
            anaddb_filestring += filename+".anaddb_band2eps\n"
            anaddb_filestring += filename+".anaddb_dummy1\n"
            anaddb_filestring += filename+".anaddb_dummy2\n"
            anaddb_filestring += filename+".anaddb_dummy3\n"
            anaddb_filestring += "\n"
            f.write(str(anaddb_filestring))
            f.close()
            #
            filename = self.filename.replace(".in","")
            f = open(self.filename.replace(".in",".thermo.files"), "w")
            thermo_filestring  = filename+".thermo.in\n"
            thermo_filestring += filename+".thermo.out\n"
            thermo_filestring += filename+".mrgddb.out\n"
            thermo_filestring += filename+".thermo_dummy\n"
            thermo_filestring += filename+".thermo_dummy1\n"
            thermo_filestring += filename+".thermo_dummy2\n"
            thermo_filestring += filename+".thermo_dummy3\n"
            thermo_filestring += "\n"
            f.write(str(thermo_filestring))
            f.close()
        #
        tasks = self.run_type
        if tasks == "nonscf+dos":
            f = open("dos.plot", "w")
            dos_plot_filestring  = "#!/usr/local/bin/gnuplot -persist\n"
            dos_plot_filestring += "# Last modified: 2018/12/17\n"
            dos_plot_filestring += "set terminal postscript eps enhanced 28 lw 2\n"
            dos_plot_filestring += 'set output "'+self.filename.replace(".in","")+'.dos.eps"\n'
            dos_plot_filestring += "unset key\n"
            dos_plot_filestring += "set xrange [-12:5]\n"
            dos_plot_filestring += "set xlabel 'Energy / eV'\n"
            dos_plot_filestring += "set ylabel 'Density of states / eV'\n"
            dos_plot_filestring += 'dos_filename = "'+self.filename.replace(".in","")+'_dos_DOS"\n'
            dos_plot_filestring += 'ef = system("cat " . dos_filename . " | grep " . "Fermi" . "'+" | awk '{print $5}'"+'")\n'
            dos_plot_filestring += "set yzeroaxis lt 2\n"
            dos_plot_filestring += "plot dos_filename using ($1-ef)*13.6058*2:($2/13.6058/2) w l \n"
            f.write(str(dos_plot_filestring))
            f.close()
        #
        return filestring


class Abinit_Input_File(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an abinit run and the method
    __str__ that outputs the contents of a abinit input file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("bohr")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        self.printbraces = False
    def __str__(self):
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        #filestring = self.docstring
        # VASP needs lattice vector matrix to have positive triple product
        if det3(lattice) < 0:
            if lattice[0].length() == lattice[1].length() == lattice[2].length():
                # Shift the first and last for cubic lattices
                transmtx = [[0, 0, 1],
                            [0, 1, 0],
                            [1, 0, 0]]
            else:
                # Else shift the two shortest
                if lattice[0].length() > lattice[1].length() and lattice[0].length() > lattice[2].length():
                    transmtx = [[1, 0, 0],
                                [0, 0, 1],
                                [0, 1, 0]]
                elif lattice[1].length() > lattice[2].length() and lattice[1].length() > lattice[0].length():
                    transmtx = [[0, 0, 1],
                                [0, 1, 0],
                                [1, 0, 0]]
                else:
                    transmtx = [[0, 1, 0],
                                [1, 0, 0],
                                [0, 0, 1]]
        else:
            transmtx = [[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
        lattice = mmmult3(transmtx,lattice)
        # Print braces around values (or not)
        if self.printbraces:
            lbrace = "{"
            rbrace = "}"
        else:
            lbrace = " "
            rbrace = ""
        #
        # length scale and lattice
        filestring  = "# Structural parameters\n"
        filestring += "acell "+lbrace+"  3*"+str(a)+" "+rbrace+" \n\n"
        filestring += "rprim "+lbrace
        for vec in lattice:
            filestring += str(Vector(vec))+"\n       "
        filestring = filestring[:-1]
        filestring += rbrace+" \n"
        # The atom position info
        alloy = False
        spcs = ""
        typatstring = "typat  "+lbrace+" "
        natom = 0
        ntypat = 0
        znuclstring = "znucl  "+lbrace+" "
        alloystring = ""
        xredstring = "xred "+lbrace+" "
        for a in self.cell.atomdata:
            for b in a:
                natom += 1
                if spcs != b.spcstring():
                    ntypat += 1
                    if b.alloy():
                        znuclstring += "?? "
                        alloystring += b.spcstring()+" "
                        alloy = True
                    else:
                        znuclstring += str(ed.elementnr[b.spcstring()])+" "
                typatstring += str(ntypat)+" "
                xredstring += str(Vector(mvmult3(transmtx,b.position)))+"\n       "
                spcs = b.spcstring()
        filestring += "natom  "+lbrace+" "+str(natom)+" "+rbrace+" \n"
        filestring += "ntypat "+lbrace+" "+str(ntypat)+" "+rbrace+" \n"
        filestring += typatstring+rbrace+" \n"
        filestring += znuclstring+rbrace+" "
        if alloy:
            filestring += "    # "+alloystring
        filestring += "\n"
        filestring += xredstring
        filestring = filestring[:-2]+rbrace+" \n"
        #
        filestring += "# cut off enegy (Ha)\n"
        #
        # read potcar files and put in outstring
        files = open(self.filename.replace(".in",".files"),"r")
        lines = files.readlines()
        nline = 0
        potlist = []
        max_ecut = 0
        for line in lines: 
            if nline >= 5:
                potlist.append(line)
            nline += 1
        for f in potlist:
            fe = open(f.replace("\n",""),"r")
            lines = fe.readlines()
            fe.close()
            for line in lines:
                if search("pw_ecut",line):
                    ecut = float(line.split("high")[1].lstrip(" ").split()[0].strip(string.punctuation))
                    if ecut >= max_ecut:
                        max_ecut = ecut
        #
        if max_ecut == 0 and self.pseudotype == "ONCVPSP":
            filestring += "ecut 50\n"
        elif max_ecut == 0 and self.pseudotype == "GBRV":
            filestring += "ecut 18\n"
            filestring += "pawecutdg 18\n"
        elif max_ecut == 0 and self.pseudotype == "JTH":
            filestring += "ecut "+str(max_ecut)+"\n"
        else:
            filestring += "ecut 18\n"
            filestring += "pawecutdg 18\n"
        filestring += "\n"
        filestring += "# 1: non-spin, 2: spin\n"
        filestring += "nsppol 2\n"
        filestring += "\n"
        filestring += "# occuupation option\n"
        filestring += "occopt 7 # Gaussian smearing\n"
        filestring += "\n"
        #
        atomic_number = 0
        spinorb = ".false."
        for a in self.cell.atomdata:
            for b in a:
                atomic_number = int(ed.elementnr[b.spcstring()])
                if atomic_number >= 50:
                    spinorb = ".true."
        if spinorb == ".true.":
            filestring += "# atomic number >= 50\n"
            filestring += "nspinor 2\n"
            filestring += "\n"
        #
        tasks = self.run_type
        times = 1
        if tasks == "scf":
             filestring += "# energy convrgence (ca. 1 meV/atom)\n"
             filestring += "toldfe %-4.1e  # Ha unit\n"%(3.675e-5*natom)
             filestring += "\n"
        elif tasks == "opt":
             filestring += "# energy convrgence (ca. 1 meV/atom)\n"
             filestring += "toldfe %-4.1e  # Ha unit\n"%(3.675e-5*natom)
             filestring += "\n"
             #
             filestring += "#Number of Data Sets\n"
             filestring += "ndtset 2\n"
             filestring += "\n"
             #
             filestring += "#Cell Optimization\n"
             filestring += "#Data Set 1\n"
             filestring += "optcell1 0\n"
             filestring += "ionmov1 2\n"
             filestring += "getcell1 0\n"
             filestring += "getxred1 0\n"
             filestring += "ntime1 10\n"
             filestring += "dilatmx1 1.05\n"
             filestring += "ecutsm1 0.5\n"
             filestring += "tolmxf1 5.00e-005\n"
             filestring += "\n"
             #
             filestring += "#Data Set 2\n"
             filestring += "optcell2 1\n"
             filestring += "ionmov2 3\n"
             filestring += "getcell2 0\n"
             filestring += "getxred2 -1\n"
             filestring += "ntime2 10\n"
             filestring += "dilatmx2 1.05\n"
             filestring += "ecutsm2 0.5\n"
             filestring += "tolmxf2 5.00e-005\n"
             filestring += "\n"
             #
             filestring += "prtposcar2 1 \n"
             filestring += "\n"
        elif tasks == "nonscf+dos":
            filestring += "iscf -2 # non-scf calculation\n"
            filestring += "prtdos 2 # 2:TDOS, 3:PDOS\n"
            filestring += "pawprtdos 0 # 0:TDOS, 3: PODS\n"
            filestring += "dosdeltae 0.001 # Ha unit\n"
            filestring += "tolwfr 1.0e-10\n"
            filestring += "\n"
            times = 2
        elif tasks == "nonscf+band":
            filestring += "tolwfr 1.0e-10\n"
            filestring += "enunit 1 # eV\n"
            filestring += "\n"
        elif tasks == "hybrid":
            filestring += "# Hybrid functional calculation\n"
            filestring += "# in a self-consistent approach\n"
            filestring += "# Dataset 1: ground state calculation with WFK output\n"
            filestring += "# Dataset 2: calculation of the HSE06 first iteration\n"
            filestring += "\n"
            filestring += "ndtset 2\n"
            filestring += "gwpara 2\n"
            filestring += "enunit 1\n"
            filestring += "gw_qprange -14  # Compute correction for all the bands\n"
            filestring += "\n"
            filestring += "# Dataset1: usual self-consistent ground-state calculation\n"
            filestring += "# Definition of the k-point grid\n"
            filestring += "#occopt 1           # Semiconductor case\n"
            filestring += "occopt 3\n"     
            filestring += "tsmear 0.002\n"
            filestring += "\n"
            #filestring += "# Definition of the k-point grid\n"
            #filestring += "kptopt 1             # Option for the automatic generation of k points, taking\n"
            #filestring += "                     # into account the symmetry\n"
            #filestring += "ngkpt   4 4 4\n"
            #filestring += "nshiftk 1\n"
            #filestring += "shiftk  0.0 0.0 0.0  # The mesh contains the Gamma point\n"
            #filestring += "                     # so that we can evaluate the QP correction for this point.\n"
            #filestring += "istwfk      *1       # Option needed for Gamma\n"
            #filestring += "\n"
            filestring += "# Common to all hybrid calculations\n"
            filestring += "getkss       1        # Obtain KSS file from previous dataset\n"
            filestring += "ecutwfn      11.5     # Planewaves to be used to represent the wavefunctions\n"
            filestring += "ecutsigx     11.5     # Planewaves to be used to represent the exchange operator\n"
            filestring += "#nkptgw      1\n"
            filestring += "#bdgw        1  8\n"
            filestring += "#kptgw       0.0 0.0 0.0\n"
            filestring += "#symsigma    1\n"
            filestring += "getqps      -1\n"
            filestring += "\n"
            filestring += "# Dataset2: Calculation of the HSE06 iteration\n"
            filestring += "optdriver2   4\n"
            filestring += "gwcalctyp    125 # self-consistent\n"
            filestring += "icutcoul2    5     # short-range exchange only\n"
            filestring += "rcut2  9.090909  # corresponds to omega = 1/rc = 0.11 bohr^1\n"
            filestring += "\n"
            filestring += "## Dataset2: Calculation of the PBE0 or B3LYP band gap\n"
            filestring += "# optdriver2 4\n"
            filestring += "# gwcalctyp2 225 # 225:PBE0, 325:B3LYP\n"
            filestring += "# icutcoul2  6   # full range exchange\n"
            filestring += "# Definition of the planewave basis set\n"    
            filestring += "\n"
            filestring += "# Definition of the SCF procedure\n"
            filestring += "nstep   250      # Maximal number of SCF cycles\n"
            filestring += "diemac  12.0     # Although this is not mandatory, it is worth to\n"
            filestring += "                 # precondition the SCF cycle. The model dielectric\n"
            filestring += "                 # function used as the standard preconditioner\n"
            filestring += "                 # is described in the dielng input variable section.\n"
            filestring += "                 # Here, we follow the prescription for bulk silicon.\n"
            filestring += "\n"
            filestring += "tolvrs   1.0d-15\n"
            filestring += "\n"
        elif tasks == "tddft":
            filestring += "# model dielectric macroscopic constant\n"
            filestring += "diemac 2.0d0\n"
            filestring += "diemix 0.5d0\n"
            filestring += "\n"
            filestring += "# number of DATASET\n"
            filestring += "ndtset 2\n"
            filestring += "\n"
            filestring += "# DATASET 1 SCF\n"
            filestring += "tolwfr1 1.0d-9\n"
            filestring += "nband1 5\n"
            filestring += "prtden1 1\n"
            filestring += "getwfk1 0\n"
            filestring += "\n"
            filestring += "# DATASET 2 TDDFT\n"
            filestring += "iscf2 -1\n"
            filestring += "tolwfr2 1.0d-9\n"
            filestring += "getden2 1\n"
            filestring += "getwfk2 1\n"
            filestring += "\n"
        elif tasks == "optic":
            filestring += "# model dielectric macroscopic constant\n"
            filestring += "diemac 2.0d0\n"
            filestring += "\n"
            filestring += "# number of DATASET\n"
            filestring += "ndtset 6\n"
            filestring += "\n"
            filestring += "# DATASET 1 SCF\n"
            filestring += "kptopt1 1\n"
            filestring += "prtden1 1\n"
            filestring += "\n"
            filestring += "# DATASET 2 NSCF with IBZ\n"
            filestring += "iscf2 -1\n"
            filestring += "kptopt2 1\n"
            filestring += "getwfk2 1\n"
            filestring += "getden2 1\n"
            filestring += "\n"
            filestring += "# DATASET 3 NSCF with BZ\n"
            filestring += "iscf3 -1\n"
            filestring += "kptopt3 3\n"
            filestring += "getwfk3 2\n"
            filestring += "getden3 1\n"
            filestring += "\n"
            axis_list = ["x","y","z"]
            for nset in range(3):
                filestring += "# DATASET "+str(4+nset)+" DDK (axis "+axis_list[nset]+")\n"
                filestring += "iscf"+str(4+nset)+" -3\n"
                filestring += "nline"+str(4+nset)+" 0\n"
                filestring += "prtwf"+str(4+nset)+" 3\n"
                filestring += "kptopt"+str(4+nset)+" 3\n"
                filestring += "nqpt"+str(4+nset)+" 1\n"
                if axis_list[nset] == "x":
                    filestring += "rfdir"+str(4+nset)+" 1 0 0\n"
                if axis_list[nset] == "y":
                    filestring += "rfdir"+str(4+nset)+" 0 1 0\n"
                if axis_list[nset] == "z":
                    filestring += "rfdir"+str(4+nset)+" 0 0 1\n"
                filestring += "rfelfd"+str(4+nset)+" 2\n"
                filestring += "getwfk"+str(4+nset)+" 3\n"
                filestring += "\n"
        #
        if self.dftu or self.dftuj:
            filestring += "usepawu 2 #1:FULL, 2:AMF\n"
            tmp_lpawu  = "lpawu "
            tmp_upawu  = "upawu "
            tmp_jpawu  = "jpawu "
            natoms = 0
            spcstring = ""
            for a in self.cell.atomdata:
                for b in a:
                    spcs = b.spcstring()
                    if spcs != spcstring:
                        natoms += 1
                        atomic_number = int(ed.elementnr[b.spcstring()])
                        if 21 <= atomic_number and atomic_number <= 29:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.81 + 0.08*(atomic_number-21)
                            U = F0
                            Ueff = U - J
                        elif 39 <= atomic_number and atomic_number <= 47:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.59 + 0.056*(atomic_number-39)
                            U = F0
                            Ueff = U - J
                        elif 71 <= atomic_number and atomic_number <= 79:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.860 + 0.053*(atomic_number-71)
                            U = F0
                            Ueff = U - J
                        elif 57 <= atomic_number and atomic_number <= 70:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.90 + 0.036*(atomic_number-57)
                            U = F0
                            Ueff = U - J
                        elif 89 <= atomic_number and atomic_number <= 103:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.66 + 0.035*(atomic_number-57)
                            U = F0
                            Ueff = U - J
                        else:
                            l = -1
                            U = 0.00 
                            J = 0.00 
                            Ueff = 0.00
                        #
                        tmp_lpawu += str(l)+" "
                        if self.dftu:
                            Ueff = Ueff * float(self.dftutimes)
                            J = 0.0
                            tmp_upawu += "%5.3f "%Ueff
                            tmp_jpawu += "%5.3f "%J
                        else:
                            U = U * float(self.dftutimes)
                            J = J * float(self.dftutimes)
                            tmp_upawu += "%5.3f "%U
                            tmp_jpawu += "%5.3f "%J
                    spcstring = spcs 
            filestring += tmp_lpawu+"\n"
            filestring += tmp_upawu+"eV\n"
            filestring += tmp_jpawu+"eV\n"
            filestring += "\n"
        #
        if tasks != "nonscf+band": 
            #
            filestring += "# k-point settings\n" # I think of (SCF: 0.6/Ang, DOS: 0.2/Ang)
            #
            kresolution = self.kresolution
            reclatvect = self.cell.reciprocal_latticevectors()
            for j in range(3):
                for i in range(3):
                    reclatvect[j][i] = reclatvect[j][i] / self.cell.lengthscale
            # Lengths of reciprocal lattice vectors
            reclatvectlen = [elem.length() for elem in reclatvect]
            self.kgrid = [max(1,int(round(elem/0.52917/kresolution))) for elem in reclatvectlen]
            #
            qresolution = self.qresolution
            self.qgrid = [max(1,int(round(elem/0.52917/qresolution))) for elem in reclatvectlen]
            #
            if self.kpeven == "yes":
                if self.cell.crystal_system() == "hexagonal":
                    kshift_list = [2]
                elif self.cell.crystal_system() == "tetragonal":
                    kshift_list = []
                elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                    kshift_list = [0,1,2]
                elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                    kshift_list = [0,1,2]
                elif self.cell.spacegroupsetting == "P":
                    kshift_list = [0,1,2]
                else:
                    kshift_list = []
                #
                for k_no in kshift_list:
                    if self.kgrid[k_no]%2 != 0 and self.kgrid[k_no] > 1:
                        self.kgrid[k_no] = int(self.kgrid[k_no]-0.1)
                    if self.qgrid[k_no] >= self.kgrid[k_no]:
                        self.qgrid[k_no] = self.kgrid[k_no]
                    elif self.qgrid[k_no]%2 != 0:
                        self.qgrid[k_no] = self.kgrid[k_no]/2
            #
            if tasks != "phonon":
                filestring += "ngkpt "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"\n"
            else:
                filestring += "ngkpt %3i %3i %3i\n"%(self.kgrid[0], self.kgrid[1], self.kgrid[2])     
            filestring += "\n"
            #
            # Monkhorst-Pack method (Special point method)
            filestring += "# k-point shift settings\n"
            if tasks == "phonon":
                filestring += "nshiftk 1\n"
                filestring += "shiftk 0 0 0\n"
                filestring += "\n"
                filestring += "nstep 1\n"
                filestring += "nline 1\n"
                phonon_pre_filestring = filestring
                f = open(self.filename.replace(".in",".kpt")+".in", "w")
                f.write(str(phonon_pre_filestring))
                f.close()
            elif self.cell.crystal_system() == "hexagonal":
                #brvtyp = "hcp"
                filestring += "nshiftk 1\n"
                filestring += "shiftk 0.0  0.0  0.5\n"
            elif self.cell.crystal_system() == "tetragonal":
                #brvtyp = "tetra"
                filestring += "nshiftk 1\n"
                filestring += "shiftk 0.0  0.0  0.0\n"
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                #brvtyp = "fcc"
                filestring += "nshiftk 4\n"
                filestring += "shiftk 0.5 0.5 0.5\n"
                filestring += "       0.5 0.0 0.0\n"
                filestring += "       0.0 0.5 0.0\n"
                filestring += "       0.0 0.0 0.5\n"
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                #brvtyp = "bcc"
                #filestring += "nshiftk 2\n"
                #filestring += " 0.25  0.25  0.25\n"
                #filestring += "-0.25 -0.25 -0.25\n"
                filestring += "nshiftk 1\n"
                filestring += "shiftk 0.5  0.5  0.5\n"
            elif self.cell.spacegroupsetting == "P":
                #brvtyp = "sc"
                filestring += "nshiftk 1\n"
                filestring += "shiftk 0.5  0.5  0.5\n"
            else:
                filestring += "nshiftk 1\n"
                filestring += "shiftk 0.0  0.0  0.0\n"  
            filestring += "\n"
        else:
            if self.cell.crystal_system() == "hexagonal":
                # hcp, primitive cell
                filestring += "# band dispersion\n"
                filestring += "kptopt -4 #kptbounds and ndivk\n"
                filestring += "kptbounds #kptopt -4 -> 4+1 = 5 points\n"
                filestring += "0.0   0.0   0.0  : vlvp1d, G-point\n"
                filestring += "0.5   0.0   0.0  : M-point\n"
                filestring += "0.333 0.333 0.0  : K-point\n"
                filestring += "0.0   0.0   0.0  : G-point\n"
                filestring += "0.0   0.0   0.5  : A-point\n"
                filestring += "ndivk 20 40 40 20 # divisions of the 4 segments, delimited by 5 points.\n"
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                # fcc, primitive cell
                filestring += "# band dispersion\n"
                filestring += "kptopt -5 #kptbounds and ndivk\n"
                filestring += "kptbounds #kptopt -5 -> 5+1 = 6 points\n"
                filestring += "0.5   0.25  0.75  : W-point\n"
                filestring += "0.5   0.0   0.0   : L-point\n"
                filestring += "0.0   0.0   0.0   : G-point\n"
                filestring += "0.5   0.5   0.0   : X-point\n"
                filestring += "0.75  0.5   0.25  : W-point\n"
                filestring += "0.75  0.375 0.375 : K-point\n"
                filestring += "ndivk 20 40 40 20 10 # divisions of the 5 segments, delimited by 6 points.\n"
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                # bcc, primitive cell
                filestring += "# band dispersion\n"
                filestring += "kptopt -4 #kptbounds and ndivk\n"
                filestring += "kptbounds #kptopt -4 -> 4+1 = 5 points\n"
                filestring += "0.0  0.0  0.0  : vlvp1d, G-point\n"
                filestring += "0.5 -0.5  0.5  : H-point\n"
                filestring += "0.0  0.0  0.5  : N-point\n"
                filestring += "0.0  0.0  0.0  : G-point\n"
                filestring += "0.25 0.25 0.25 : P-point\n"
                filestring += "ndivk 20 40 40 20 # divisions of the 4 segments, delimited by 5 points.\n"
            elif self.cell.spacegroupsetting == "P":
                # sc, primitive cell
                filestring += "# band dispersion\n"
                filestring += "kptopt -4 #kptbounds and ndivk\n"
                filestring += "kptbounds #kptopt -4 -> 4+1 = 5 points\n"
                filestring += "0.5  0.5  0.5  : vlvp1d, R-point\n"
                filestring += "0.0  0.0  0.0  : G-point\n"
                filestring += "0.5  0.0  0.0  : X-point\n"
                filestring += "0.5  0.5  0.0  : M-point\n"
                filestring += "0.0  0.0  0.0  : G-point\n"
                filestring += "ndivk 20 40 40 20 # divisions of the 4 segments, delimited by 5 points.\n"
            else:
                # sc, primitive cell
                filestring += "# band dispersion\n"
                filestring += "kptopt -4 #kptbounds and ndivk\n"
                filestring += "kptbounds #kptopt -4 -> 4+1 = 5 points\n"
                filestring += "0.5  0.5  0.5  : vlvp1d, R-point\n"
                filestring += "0.0  0.0  0.0  : G-point\n"
                filestring += "0.5  0.0  0.0  : X-point\n"
                filestring += "0.5  0.5  0.0  : M-point\n"
                filestring += "0.0  0.0  0.0  : G-point\n"
                filestring += "ndivk 20 40 40 20 # divisions of the 4 segments, delimited by 5 points.\n"
            #
            filestring += "\n"
        #
        if tasks == "phonon":
            filestring += "# model dielectric macroscopic constant\n"
            filestring += "diemac 2.0\n"
            filestring += "\n"
            filestring += "# number of DATASET\n"
            noq = self.qgrid[0]*self.qgrid[1]*self.qgrid[2]
            filestring += "ndtset "+str(2+noq)+"\n"
            filestring += "\n"
            filestring += "# DATASET 1 SCF\n"
            filestring += "getwfk1 0 # Cancel default\n"
            filestring += "kptopt1 1 # Automatic generation of k points, taking\n"
            filestring += "          # into account the symmetry\n"
            filestring += "nqpt1   0 # Cancel default\n"
            filestring += "tolvrs1 1.0d-18 # SCF stopping criterion (modify default)\n"
            filestring += "rfphon1 0 # Cancel default\n"
            filestring += "\n"
            filestring += "# One qpt for each dataset\n"
            filestring += "nqpt 1\n"
            filestring += "\n"
            filestring += "# DATASET 2 DDK, Response function calculation of d/dk wave function \n"
            filestring += "iscf2  -3 # Need this non-self-consistent option for d/dk\n"
            filestring += "kptopt2 2 # Modify default to use time-reversal symmetry\n"
            filestring += "rfphon2 0 # Cancel default\n"
            filestring += "rfelfd2 2 # Calculate d/dk wave function only\n"
            filestring += "tolwfr2 1.0d-9 # Use wave function residual criterion instead\n"
            filestring += "qpt2   %7.5e  %7.5e  %7.5e\n"%(0.0, 0.0, 0.0)
            filestring += "\n"
            filestring += "# DATASET 3 Response function calculation of q=0 phonons and electric field pert.\n"
            filestring += "getddk3 2 # d/dk wave functions from last dataset\n"
            filestring += "kptopt3 2 # Modify default to use time-reversal symmetry\n"
            filestring += "rfelfd3 3 # Electric-field perturbation response only\n"
            filestring += "qpt3   %7.5e  %7.5e  %7.5e\n"%(0.0, 0.0, 0.0)
            filestring += "\n"
            filestring += "# DATASET 4-"+str(2+noq)+" Finite-wave-vector phonon calculations (defaults for all datasets)\n"
            filestring += "getwfk  1 # Use GS wave functions from dataset1\n"
            filestring += "kptopt  3 # Need full k-point set for finite-q response\n"
            filestring += "rfphon  1 # Do phonon response\n"
            filestring += "rfatpol 1 "+str(natom)+" # Treat displacements of all atoms\n"
            filestring += "rfdir   1 1 1 # Do all directions (symmetry will be used)\n"
            filestring += "tolvrs  1.0d-8 # This default is active for sets 3-"+str(2+noq)+"\n"
            filestring += "\n"
            filestring += "# q point for # DATASET 4-"+str(2+noq)+"\n"
            #noq = 0
            #for noa in range(self.qgrid[0]):
            #    for nob in range(self.qgrid[1]):
            #        for noc in range(self.qgrid[2]):
            #            noq += 1
            #            if noq > 1:
            #                filestring += "qpt%-3i %7.5e  %7.5e  %7.5e\n"%(int(2+noq), float(noa)/float(self.qgrid[0]), float(nob)/float(self.qgrid[1]), float(noc)/float(self.qgrid[2]))
        #
        if tasks == "phonon":
            f = open(self.filename.replace(".in",".mrgddb")+".in", "w")
            mrgddb_filestring  = self.filename.replace(".in",".mrgddb")+".out\n"
            mrgddb_filestring += self.filename.replace(".in","")+" phonons on %3i %3i %3i  mesh\n"%(self.qgrid[0], self.qgrid[1], self.qgrid[2])
            mrgddb_filestring += str(noq)+"\n"
            for dsno in range(noq):
                mrgddb_filestring += "trf2_1o_DS"+str(3+dsno)+"_DDB\n"
            mrgddb_filestring += "\n"
            f.write(str(mrgddb_filestring))
            f.close()
            #
            f = open(self.filename.replace(".in",".anaddb")+".in", "w")
            anaddb_filestring  = "!Input file for the anaddb code.\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Flags\n"
            anaddb_filestring += "  ifcflag 1 ! Interatomic force constant flag\n"
            anaddb_filestring += "  ifcout 0\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Wavevector grid number 1 (coarse grid, from DDB)\n"
            if self.cell.crystal_system() == "hexagonal":
                kshift_list = [2]
                brav = 4
            elif self.cell.crystal_system() == "tetragonal":
                kshift_list = []
                brav = 1
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                kshift_list = [0,1,2]
                brav = 2
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                kshift_list = [0,1,2]
                brav = 3
            elif self.cell.spacegroupsetting == "P":
                kshift_list = [0,1,2]
                brav = 1
            anaddb_filestring += "  brav "+str(brav)+" ! Bravais Lattice : 1-S.C., 2-F.C., 3-B.C., 4-Hex.)\n"
            anaddb_filestring += "  ngqpt %3i %3i %3i ! Monkhorst-Pack indices\n"%(self.qgrid[0], self.qgrid[1], self.qgrid[2])
            anaddb_filestring += "  nqshft 1 ! number of q-points in repeated basic q-cell\n"
            anaddb_filestring += "  q1shft 3*0.0\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Effective charges\n"
            anaddb_filestring += "  chneut 1 ! Charge neutrality requirement for effective charges.\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Interatomic force constant info\n"
            anaddb_filestring += "  dipdip 1 ! Dipole-dipole interaction treatment\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Phonon band structure output for band2eps - See note near end for\n"
            anaddb_filestring += "! dealing with gamma LO-TO splitting issue.\n"
            anaddb_filestring += "  eivec  4\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Wavevector list number 1 (Reduced coordinates and normalization factor), FCC\n"  
            anaddb_filestring += "  nph1l 8 ! number of phonons in list 1\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "  qph1l   0.0000  0.0000  0.0000   1.0    !(gamma point)\n"
            anaddb_filestring += "          0.3750  0.3750  0.7500   1.0    !(K point)\n"
            anaddb_filestring += "          0.5000  0.5000  1.0000   1.0    !(X point)\n"
            anaddb_filestring += "          1.0000  1.0000  1.0000   1.0    !(gamma point)\n"
            anaddb_filestring += "          0.5000  0.5000  0.5000   1.0    !(L point)\n"
            anaddb_filestring += "          0.5000  0.0000  0.5000   1.0    !(X point)\n"
            anaddb_filestring += "          0.5000  0.2500  0.7500   1.0    !(W point)\n"
            anaddb_filestring += "          0.5000  0.5000  0.5000   1.0    !(L point)\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Wavevector list number 2 (Cartesian directions for non-analytic gamma phonons)\n"
            anaddb_filestring += "  nph2l 1 ! number of directions in list 2\n"
            anaddb_filestring += "  qph2l   1.0  0.0  0.0    0.0\n"
            anaddb_filestring += "\n"
            f.write(str(anaddb_filestring))
            f.close()
            #
            f = open(self.filename.replace(".in",".thermo")+".in", "w")
            anaddb_filestring  = "!Input file for the anaddb code.\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Flags\n"
            anaddb_filestring += "  ifcflag 1 ! Interatomic force constant flag\n"
            anaddb_filestring += "  ifcout 0\n"
            anaddb_filestring += "  thmflag 1 ! Thermodynamical properties flag\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Wavevector grid number 1 (coarse grid, from DDB)\n"
            if self.cell.crystal_system() == "hexagonal":
                kshift_list = [2]
                brav = 4
            elif self.cell.crystal_system() == "tetragonal":
                kshift_list = []
                brav = 1
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                kshift_list = [0,1,2]
                brav = 2
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                kshift_list = [0,1,2]
                brav = 3
            elif self.cell.spacegroupsetting == "P":
                kshift_list = [0,1,2]
                brav = 1
            anaddb_filestring += "  brav "+str(brav)+" ! Bravais Lattice : 1-S.C., 2-F.C., 3-B.C., 4-Hex.)\n"
            anaddb_filestring += "  ngqpt %3i %3i %3i ! Monkhorst-Pack indices\n"%(self.qgrid[0], self.qgrid[1], self.qgrid[2])
            anaddb_filestring += "  nqshft 1 ! number of q-points in repeated basic q-cell\n"
            anaddb_filestring += "  q1shft 3*0.0\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Effective charges\n"
            anaddb_filestring += "  chneut 1 ! Charge neutrality requirement for effective charges.\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Interatomic force constant info\n"
            anaddb_filestring += "  dipdip 1 ! Dipole-dipole interaction treatment\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Wavevector grid number 2 (series of fine grids, extrapolated from interat forces)\n"
            anaddb_filestring += "  ng2qpt 20 20 20  ! sample the BZ up to ngqpt2\n"
            anaddb_filestring += "  ngrids 5         ! number of grids of increasing size\n"
            anaddb_filestring += "  q2shft 3*0.0\n"
            anaddb_filestring += "\n"
            anaddb_filestring += "!Thermal information\n"
            anaddb_filestring += "  nchan   1250   ! # of channels for the DOS with channel width 1 cm-1\n"
            anaddb_filestring += "  nwchan  5      ! # of different channel widths from this integer down to 1 cm-1\n"
            anaddb_filestring += "  thmtol  0.120  ! Tolerance on thermodynamical function fluctuations\n"
            anaddb_filestring += "  ntemper 10     ! Number of temperatures\n"
            anaddb_filestring += "  temperinc 20.  ! Increment of temperature in K for temperature dependency\n"
            anaddb_filestring += "  tempermin 20.  ! Minimal temperature in Kelvin\n"
            anaddb_filestring += "\n"
            f.write(str(anaddb_filestring))
            f.close() 
        #
        return filestring

################################################################################################
class AIMSFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a FHI-AIMS run and the method
    __str__ that outputs the contents of a abinit input file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        self.cartesian = False
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
    def __str__(self):
        filestring = self.docstring+"\n"
        latvecs = self.cell.latticevectors
        for i in range(3):
            for j in range(3):
                latvecs[i][j] = latvecs[i][j] * self.cell.lengthscale
        for vec in latvecs:
            filestring += "lattice_vector  "+str(vec)+"\n"
        for a in self.cell.atomdata:
            for b in a:
                if self.cartesian:
                    filestring += "atom  "+str(Vector(mvmult3(latvecs,b.position)))+" "+b.spcstring()+"\n"
                else:
                    filestring += "atom_frac  "+str(b.position)+" "+b.spcstring()+"\n"
        return filestring
        
################################################################################################
#UNFINISHED!
class MCSQSFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed by the mcsqs SQS generator and the method
    __str__ that outputs the contents of a mcsqs input file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        self.cartesian = False
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
    def __str__(self):
        l = self.cell.lengthscale
        filestring = "%12.8f %12.8f %12.8f 90.0 90.0 90.0"
        filestring += str(self.cell.latticevectors)
        for a in self.cell.atomdata:
            for b in a:
                filestring += +str(b.position)+" "
                for k,v in b.species.iteritems():
                    filestring += k+"="+str(v)
        return filestring
        
################################################################################################
class POSCARFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a POSCAR file and the method
    __str__ that outputs the contents of a POSCAR file as a string.
    If you want POSCAR to be printed with the atomic positions in Cartesian form,
    then set
    POSCARFile.printcartpos = True
    If you want to put the overall length scale on the lattice vectors and print 1.0
    for the length scale, then set
    POSCARFile.printcartvecs = True
    """
    def __init__(self, crystalstructure, string, vca=False):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        self.printcartvecs = False
        self.printcartpos = False
        self.vasp5format = False
        self.selectivedyn = False
        self.vca = vca
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                if self.vca:
                    for k,v in b.species.iteritems():
                        tmp.add(k)
                else:
                    tmp.add(b.spcstring())
        self.species = list(tmp)
        # make sure the docstring goes on one line
        self.docstring = self.docstring.replace("\n"," ")
    def SpeciesOrder(self):
        """
        Return a string with the species in the order they appear in POSCAR.
        """
        returnstring = ""
        for sp in self.species:
            returnstring += sp+" "
        return returnstring
    def __str__(self):
        # Assign some local variables
        lattice = self.cell.latticevectors
        # VASP needs lattice vector matrix to have positive triple product
        if det3(lattice) < 0:
            if lattice[0].length() == lattice[1].length() == lattice[2].length():
                # Shift the first and last for cubic lattices
                transmtx = [[0, 0, 1],
                            [0, 1, 0],
                            [1, 0, 0]]
            else:
                # Else shift the two shortest
                if lattice[0].length() > lattice[1].length() and lattice[0].length() > lattice[2].length():
                    transmtx = [[1, 0, 0],
                                [0, 0, 1],
                                [0, 1, 0]]
                elif lattice[1].length() > lattice[2].length() and lattice[1].length() > lattice[0].length():
                    transmtx = [[0, 0, 1],
                                [0, 1, 0],
                                [1, 0, 0]]
                else:
                    transmtx = [[0, 1, 0],
                                [1, 0, 0],
                                [0, 0, 1]]
        else:
            transmtx = [[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
        lattice = mmmult3(transmtx,lattice)
                
        # For output of atomic positions
        a = self.cell.lengthscale
        positionunits = ""
        if self.selectivedyn:
            positionunits += "Selective dynamics\n"
        if self.printcartpos:
            positionunits += "Cartesian\n"
            coordmat = []
            for i in range(3):
                coordmat.append([])
                for j in range(3):
                    coordmat[i].append(lattice[i][j] * a)
        else:
            coordmat = [[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
            positionunits += "Direct\n"
        # The first line with info from input docstring
        filestring = self.docstring
        if not self.vasp5format:
            filestring += " Species order: "
            for sp in self.species:
                filestring += sp+" "
        filestring += "\n"
        # Lattice parameter and vectors
        if self.printcartvecs:
            latticestring = " 1.0\n"
            for i in range(3):
                latticestring += "%19.15f %19.15f %19.15f\n" % (a*lattice[i][0], a*lattice[i][1], a*lattice[i][2])
        else:
            latticestring = " %10f\n" % a
            for i in range(3):
                latticestring += "%19.15f %19.15f %19.15f\n" % (lattice[i][0], lattice[i][1], lattice[i][2])
        filestring += latticestring
        # print species here if vasp 5 format
        if self.vasp5format:
            for sp in self.species:
                filestring += (" "+sp).rjust(4)
            filestring += "\n"
        # positions and number of species
        nspstring = ""
        positionstring = ""
        for sp in self.species:
            nsp = 0
            for a in self.cell.atomdata:
                for b in a:
                    if self.vca:
                        for k,v in b.species.iteritems():
                            if k == sp:
                                nsp += 1
                                p = Vector(mvmult3(coordmat,mvmult3(transmtx,b.position)))
                                positionstring += str(p)
                                if self.selectivedyn:
                                    positionstring += "   T  T  T"
                                positionstring += "\n"
                    else:
                        if b.spcstring() == sp:
                            nsp += 1
                            p = Vector(mvmult3(coordmat,mvmult3(transmtx,b.position)))
                            positionstring += str(p)
                            if self.selectivedyn:
                                positionstring += "   T  T  T"
                            positionstring += "\n"
            nspstring += (" "+str(nsp)).rjust(4)
        filestring += nspstring+"\n"
        filestring += positionunits
        filestring += positionstring
        return filestring

class POTCARFile:
    """
    Class for representing and outputting a POTCAR file for VASP.
    """
    def __init__(self, crystalstructure, directory="",vca=False,
                 prioritylist=["_d","_pv","_sv","","_h","_s"]):
        self.cell = crystalstructure
        self.vca = vca
        self.prioritylist = prioritylist
        # POTCAR library
        if directory != "":
            self.dir = directory
        else:
            try:
                self.dir = os.environ['VASP_PSEUDOLIB']
            except:
                try:
                    self.dir = os.environ['VASP_PAWLIB']
                except:
                    self.dir = ""
        # check directory
        if self.dir == "":
            raise SetupError("No path to the VASP pseudopotential library specified.\n")
        if not os.path.exists(self.dir):
            raise SetupError("The specified path to the VASP pseudopotential library does not exist.\n"+self.dir)
        # set up species list
        poscarfile = POSCARFile(self.cell, "", vca=self.vca)
        self.species = poscarfile.species
    def __str__(self):
        # get all files
        potcarlist = []
        for a in self.species:
            for version in self.prioritylist:
                potcarfile = self.dir+"/"+a+version+"/POTCAR"
                if os.path.exists(potcarfile):
                    potcarlist.append(potcarfile)
                    break
        # read potcar files and put in outstring
        outstring = ""
        for f in potcarlist:
            potcar = open(f,"r")
            outstring += potcar.read()
            potcar.close()
        return outstring

class KPOINTSFile:
    """
    Class for representing and outputting a KPOINTS file for VASP.
    """
    def __init__(self, crystalstructure, docstring="",kresolution=0.2,runtype="",kpeven=""):
        self.cell = crystalstructure
        self.docstring = docstring
        self.kresolution = kresolution
        self.runtype = runtype
        self.kpeven = kpeven
        # set reciprocal lattice vectors in reciprocal angstroms
        reclatvect = crystalstructure.reciprocal_latticevectors()
        for j in range(3):
            for i in range(3):
                reclatvect[j][i] = reclatvect[j][i] / crystalstructure.lengthscale
        # Lengths of reciprocal lattice vectors
        reclatvectlen = [elem.length() for elem in reclatvect]
        self.kgrid = [max(1,int(round(elem/self.kresolution))) for elem in reclatvectlen]
    def __str__(self):
        tmp = self.docstring
        if self.runtype != "nonscf+band":
            tmp += " k-space resolution ~"+str(self.kresolution)+"/A\n"
            tmp += " 0\n"
            if self.kpeven == "yes":
                if self.cell.crystal_system() == "hexagonal":
                    kshift_list = [2]
                elif self.cell.crystal_system() == "tetragonal":
                    kshift_list = []
                elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                    kshift_list = [0,1,2]
                elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                    kshift_list = [0,1,2]
                elif self.cell.spacegroupsetting == "P":
                    kshift_list = [0,1,2]
                else:
                    kshift_list = []
                #
                for k_no in kshift_list:
                    if self.kgrid[k_no]%2 != 0 and self.kgrid[k_no] > 1:
                        self.kgrid[k_no] = int(self.kgrid[k_no]-0.1)
                #
                tmp += "Monk ! Gamma|Monk|Auto\n"
                tmp += str(self.kgrid[0])+" "+str(self.kgrid[1])+" "+str(self.kgrid[2])+"\n"
                if self.cell.crystal_system() == "hexagonal":
                    tmp += "0 0 1\n"
                elif self.cell.crystal_system() == "tetragonal":
                    tmp += "0 0 0\n"
                elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                    tmp += "1 1 1\n"
                elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                    tmp += "1 1 1\n"
                elif self.cell.spacegroupsetting == "P":
                    tmp += "1 1 1\n"
                else:
                    tmp += "0 0 0\n"
            else:
                tmp += "Gamma ! Gamma|Monk|Auto\n"
                tmp += str(self.kgrid[0])+" "+str(self.kgrid[1])+" "+str(self.kgrid[2])+"\n"
                tmp += "0 0 0\n"
        else:
            if self.cell.crystal_system() == "hexagonal":
                # hcp, primitive cell
                tmp += "# band dispersion (hcp type)\n"
                tmp += "k-points for bandstructure G-M-K-G-A\n"
                tmp += " 40               ! 40 intersections \n"
                tmp += "line              ! Line-mode\n"
                tmp += "reciprocal        ! reciprocal|cart\n"
                tmp += "0.0   0.0   0.0   ! G-point\n"
                tmp += "0.5   0.0   0.0   ! M-point\n"
                tmp += "\n"
                tmp += "0.5   0.0   0.0   ! M-point\n"
                tmp += "0.333 0.333 0.0   ! K-point\n"
                tmp += "\n"
                tmp += "0.333 0.333 0.0   ! K-point\n"
                tmp += "0.0   0.0   0.0   ! G-point\n"
                tmp += "\n"
                tmp += "0.0   0.0   0.0   ! G-point\n"
                tmp += "0.0   0.0   0.5   ! A-point\n"
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "F":
                # fcc, primitive cell
                tmp += "# band dispersion (fcc type)\n"
                tmp += "k-points for bandstructure W-L-G-X-W-K\n"
                tmp += " 40               ! 40 intersections \n"
                tmp += "line              ! Line-mode\n"
                tmp += "reciprocal        ! reciprocal|cart\n"
                tmp += "0.5   0.25  0.75  ! W-point\n"
                tmp += "0.5   0.0   0.0   ! L-point\n"
                tmp += "\n"
                tmp += "0.5   0.0   0.0   ! L-point\n"
                tmp += "0.0   0.0   0.0   ! G-point\n"
                tmp += "\n"
                tmp += "0.0   0.0   0.0   ! G-point\n"
                tmp += "0.5   0.5   0.0   ! X-point\n"
                tmp += "\n"
                tmp += "0.5   0.5   0.0   ! X-point\n"
                tmp += "0.75  0.5   0.25  ! W-point\n"
                tmp += "\n"
                tmp += "0.75  0.5   0.25  ! W-point\n"
                tmp += "0.75  0.375 0.375 ! K-point\n"
            elif self.cell.crystal_system() == "cubic" and self.cell.spacegroupsetting == "I":
                # bcc, primitive cell
                tmp += "# band dispersion (bcc type)\n"
                tmp += "k-points for bandstructure G-H-N-G-P\n"
                tmp += " 40               ! 40 intersections \n"
                tmp += "line              ! Line-mode\n"
                tmp += "reciprocal        ! reciprocal|cart\n"
                tmp += "0.0  0.0  0.0     ! G-point\n"
                tmp += "0.5 -0.5  0.5     ! H-point\n"
                tmp += "\n"
                tmp += "0.5 -0.5  0.5     ! H-point\n"
                tmp += "0.0  0.0  0.5     ! N-point\n"
                tmp += "\n"
                tmp += "0.0  0.0  0.5     ! N-point\n"
                tmp += "0.0  0.0  0.0     ! G-point\n"
                tmp += "\n"
                tmp += "0.0  0.0  0.0     ! G-point\n"
                tmp += "0.25 0.25 0.25    ! P-point\n"
            elif self.cell.spacegroupsetting == "P":
                # sc, primitive cell
                tmp += "# band dispersion (simple cubic type)\n"
                tmp += "k-points for bandstructure G-H-N-G-P\n"
                tmp += " 40               ! 40 intersections \n"
                tmp += "line              ! Line-mode\n"
                tmp += "reciprocal        ! reciprocal|cart\n"
                tmp += "0.5  0.5  0.5     ! R-point\n"
                tmp += "0.0  0.0  0.0     ! G-point\n"
                tmp += "\n"
                tmp += "0.0  0.0  0.0     ! G-point\n"
                tmp += "0.5  0.0  0.0     ! X-point\n"
                tmp += "\n"
                tmp += "0.5  0.0  0.0     ! X-point\n"
                tmp += "0.5  0.5  0.0     ! M-point\n"
                tmp += "\n"
                tmp += "0.5  0.5  0.0     ! M-point\n"
                tmp += "0.0  0.0  0.0     ! G-point\n"
            else:
                # sc, primitive cell
                tmp += "# band dispersion (simple cubic type)\n"
                tmp += "k-points for bandstructure G-H-N-G-P\n"
                tmp += " 40               ! 40 intersections \n"
                tmp += "line              ! Line-mode\n"
                tmp += "reciprocal        ! reciprocal|cart\n"
                tmp += "0.5  0.5  0.5     ! R-point\n"
                tmp += "0.0  0.0  0.0     ! G-point\n"
                tmp += "\n"
                tmp += "0.0  0.0  0.0     ! G-point\n"
                tmp += "0.5  0.0  0.0     ! X-point\n"
                tmp += "\n"
                tmp += "0.5  0.0  0.0     ! X-point\n"
                tmp += "0.5  0.5  0.0     ! M-point\n"
                tmp += "\n"
                tmp += "0.5  0.5  0.0     ! M-point\n"
                tmp += "0.0  0.0  0.0     ! G-point\n"
            #
            tmp += "\n"
        #
        return tmp
        
class INCARFile:
    """
    Class for representing and outputting a INCAR file for VASP.
    """
    def __init__(self, crystalstructure, docstring="",potcardir="",vca=False,
                 prioritylist=["_d","_pv","_sv","","_h","_s"], encutfac=1.5,
                 dfttype="", runtype="",ldautype=2,ldautimes=0.25,ncore=1,spin="yes"):
        self.cell = crystalstructure
        self.docstring = "# "+docstring.lstrip("#").rstrip("\n")+"\n"
        self.prioritylist = prioritylist
        self.vca = vca
        self.vcaspecies = None
        # ecut = max(encuts found in potcars)*encutfac
        self.encutfac = encutfac
        self.dfttype = dfttype
        self.runtype = runtype
        self.ldautype = ldautype
        self.ldautimes = ldautimes
        self.ncore = ncore
        self.spin = spin
        poscarfile = POSCARFile(self.cell, "", vca=self.vca)
        # we need the potcar directory
        if potcardir != "":
            self.potcardir = potcardir
        else:
            try:
                self.potcardir = os.environ['VASP_PSEUDOLIB']
            except:
                try:
                    self.potcardir = os.environ['VASP_PAWLIB']
                except:
                    self.potcardir = ""
        # check directory
        if self.potcardir == "":
            raise SetupError("No path to the VASP pseudopotential library specified.\n")
        if not os.path.exists(self.potcardir):
            raise SetupError("The specified path to the VASP pseudopotential library does not exist.\n"+self.dir)

        if self.vca:
            # set up species list
            tmp = set([])
            for a in self.cell.atomdata:
                for b in a:
                    for k,v in b.species.iteritems():
                        tmp.add((k,v))
            tmp = list(tmp)
            self.vcaspecies = []
            for s in poscarfile.species:
                for t in tmp:
                    if t[0] == s:
                        self.vcaspecies.append(t)
        # set up species dict
        speciesdict = dict([])
        for a in self.cell.atomdata:
            for b in a:
                if self.vca:
                    for k,v in b.species.iteritems():
                        spcstr = k
                        if spcstr in speciesdict:
                            t = speciesdict[spcstr] + 1
                            speciesdict[spcstr] = t
                        else:
                            speciesdict[spcstr] = 1
                else:
                    spcstr = b.spcstring()
                    if spcstr in speciesdict:
                        t = speciesdict[spcstr] + 1
                        speciesdict[spcstr] = t
                    else:
                        speciesdict[spcstr] = 1
        # species list in the same order as poscar
        self.species = []
        for s in poscarfile.species:
            for k,v in speciesdict.iteritems():
                if k == s:
                    self.species.append((k,v))
        # get potcar list
        potcars = dict([])
        specieslist = []
        for a in speciesdict:
            for version in self.prioritylist:
                potcarfile = self.potcardir+"/"+a+version+"/POTCAR"
                if os.path.exists(potcarfile):
                    potcars[a] = potcarfile
                    specieslist.append(a)
                    break        
        # get maximal encut and number of electrons from potcars
        enmaxs = dict([])
        zvals = dict([])
        for a,f in potcars.iteritems():
            potcar = open(f,"r")
            for line in potcar:
                if search("ZVAL",line):
                    zvals[a] = float(line.split("ZVAL")[1].lstrip("= ").split()[0].strip(string.punctuation))
                if search("ENMAX",line):
                    enmaxs[a] = float(line.split("ENMAX")[1].lstrip("= ").split()[0].strip(string.punctuation))
                if search("END of PSCTR",line):
                    break
            potcar.close()
        self.maxencut = max([k for v,k in enmaxs.iteritems()])
        #self.maxencut = 300.0 # For test setting
        # do we suspect that this might be magnetic?
        self.magnetic = False
        self.magmomlist = []
        for s in self.species:
            if s[0] in suspiciouslist:
                self.magnetic = True
                self.magmomlist.append(str(s[1])+"*"+str(initialmoments[s[0]]))
            else:
                self.magmomlist.append(str(s[1])+"*0")
        # Determine NBANDS
        nmag = sum([eval(i) for i in self.magmomlist])
        nelect = 0.0
        for sp,z in zvals.iteritems():
            for a in self.cell.atomdata:
                for b in a:
                    if sp == b.spcstring():
                        nelect += z
        if nmag > 0:
            nstates = int(math.ceil(nelect))
        else:
            nstates = int(math.ceil(nelect/2))
        # NBANDS is max of the default VASP definition and occupied bands+20
        natoms = sum([len(a) for a in self.cell.atomdata])
        self.nbands = max(max(max(int(math.ceil(nelect/2))+int(natoms/2),3), math.ceil(0.6*nelect))+nmag, nstates+20)
    def __str__(self):
        tmp = self.docstring
        tmp += "\n"
        tmp += "##### General settings #####\n"
        tmp += "#NPAR = %i \n"%int(4.0-float(self.ncore)**(1/2))
        tmp += "#KPAR = %i # parallelization of k-points\n"%int(self.ncore-(4.0-float(self.ncore)**(1/2)))
        if self.runtype == "" or self.runtype == "scf" or self.runtype == "opt":
            if self.dfttype != "HSE03" and self.dfttype != "HSE06" and self.dfttype != "PBE0":
                tmp += "ALGO = Fast\n"
        ## tmp += "NELMIN = 4\n"
        tmp += "PREC = Accurate\n"
        tmp += "LREAL = Auto\n"
        tmp += "#ENCUT = "+str(self.maxencut*self.encutfac)+"\n"
        tmp += "#NBANDS = "+str(self.nbands)+"\n"
        tmp += "\n"
        tmp += "## ca. 1 meV/atom for total (free) energy change\n"
        natoms = sum([len(a) for a in self.cell.atomdata])
        tmp += "#EDIFF = %-4.1e\n"%(1.0e-4*natoms) # Ry unit ?
        tmp += "#EDIFFG = %-5.1e\n"%(1.0e-4*natoms*10)
        tmp += "\n"
        tmp += "##ISMEAR: -5 = tetra, -1 = Fermi, 0 = Gauss, 1 = MP\n"
        tmp += "#ISMEAR = 1 # Methfessel-Paxton order 1.\n"
        tmp += "#SIGMA = 0.2 # width (eV)\n"
        tmp += "\n"
        tmp += "##ISPIN: 1 = non-spin, 2 = spin\n"
        if self.spin == "yes":
            tmp += "ISPIN = 2\n" 
            if self.magnetic:
                tmp += "MAGMOM = "
                for species in self.magmomlist:
                    tmp += species+" "
                tmp += "\n"
        else:
            tmp += "ISPIN = 1\n" 
        #
        tmp += "\n"
        if self.vca:
            tmp += "##### Virtual crystal approximation settings #####\n"
            tmp += "VCA = "
            for s in self.vcaspecies:
                tmp += str(s[1])+" "
            tmp += "\n"
            tmp += "LVCADER = .True.\n"
            tmp += "\n"
            tmp += "\n"
        #
        if self.runtype == "opt":
            tmp += "#### Structure Optimization settings #####\n"
            tmp += "#IBRION: 0 = MD, 1 = RMM-DIIS, 2 = cg\n"
            tmp += "IBRION = 1\n"
            ## tmp += "POTIM = 0.4\n"
            tmp += "ISIF = 3\n"
            tmp += "NSW = 40\n"
            tmp += "\n"
        else:
            tmp += "##### Structure Optimization settings #####\n"
            tmp += "##IBRION: 0 = MD, 1 = RMM-DIIS, 2 = cg\n"
            tmp += "#IBRION = 1\n"
            ## tmp += "POTIM = 0.4\n"
            tmp += "#ISIF = 3\n"
            tmp += "#NSW = 40\n"
            tmp += "\n"
        #
            tmp += "\n"
        if self.dfttype != "" and self.dfttype != "GGA+U":
            tmp += "#### DFT setting #####\n"
            tmp += "#GGA: 91 = PW91, PE = PBE, RP = rPBE, PS = PBEsol\n"
            tmp += "GGA = "+self.dfttype+"\n"
            tmp += "\n"
        else:
            tmp += "##### DFT setting #####\n"
            tmp += "##GGA: 91 = PW91, PE = PBE, RP = rPBE, PS = PBEsol\n"
            tmp += "#GGA = PE\n"
            tmp += "\n"  
        #
        if self.dfttype == "GGA+U":
            tmp += "##### GGA+U calculation #####\n"
            tmp += "LDAU = .TURE.\n"
            tmp += "LDAUTYPE = "+str(self.ldautype)+"\n"
            tmp_LDAUL  = "LDAUL = "
            tmp_LDAUU  = "LDAUU = "
            tmp_LDAUJ  = "LDAUJ = "
            natoms = 0
            spcstring = ""
            for a in self.cell.atomdata:
                for b in a:
                    spcs = b.spcstring()
                    if spcs != spcstring:
                        natoms += 1
                        atomic_number = int(ed.elementnr[b.spcstring()])
                        if 21 <= atomic_number and atomic_number <= 11:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.81 + 0.08*(atomic_number-21)
                            U = F0
                            Ueff = U - J
                        elif 39 <= atomic_number and atomic_number <= 47:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.59 + 0.056*(atomic_number-39)
                            U = F0
                            Ueff = U - J
                        elif 71 <= atomic_number and atomic_number <= 79:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.860 + 0.053*(atomic_number-71)
                            U = F0
                            Ueff = U - J
                        elif 57 <= atomic_number and atomic_number <= 70:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.90 + 0.036*(atomic_number-57)
                            U = F0
                            Ueff = U - J
                        elif 89 <= atomic_number and atomic_number <= 103:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.66 + 0.035*(atomic_number-57)
                            U = F0
                            Ueff = U - J
                        else:
                            l = -1
                            U = 0.0
                            J = 0.0
                            Ueff = 0.0
                        #
                        tmp_LDAUL += str(l)+" "
                        if str(self.ldautype) == "1":
                            U = U * float(self.ldautimes)
                            J = J * float(self.ldautimes)
                            tmp_LDAUU += "%5.3f "%U
                            tmp_LDAUJ += "%5.3f "%J
                        else:
                            Ueff = Ueff * float(self.ldautimes)
                            J = 0.0
                            tmp_LDAUU += "%5.3f "%Ueff
                            tmp_LDAUJ += "%5.3f "%J
                    spcstring = spcs 
            tmp += tmp_LDAUL+"\n"
            tmp += tmp_LDAUU+"\n"
            tmp += tmp_LDAUJ+"\n"
            tmp += "LDAUPRINT = 2\n"
            tmp += "LMAXMIX = 4\n"
            tmp += "\n"
        #
        if self.dfttype == "B3LYP":
            tmp += "#### Hybrid function settings #####\n"
            tmp += "# Selects the B3LYP hybrid function\n"
            tmp += "LHFCALC = .TRUE. ; GGA     = B3 ;\n"
            tmp += "AEXX    = 0.2  ; AGGAX   = 0.72 ;\n"
            tmp += "AGGAC   = 0.81 ; ALDAC   = 0.19 ;\n"
            tmp += "\n"
        else:
            tmp += "##### Hybrid function settings #####\n"
            tmp += "## Selects the B3LYP hybrid function\n"
            tmp += "#LHFCALC = .TRUE. ; GGA     = B3 ;\n"
            tmp += "#AEXX    = 0.2  ; AGGAX   = 0.72 ;\n"
            tmp += "#AGGAC   = 0.81 ; ALDAC   = 0.19 ;\n"
            tmp += "\n"
        #
        if self.dfttype == "HSE03":
            tmp += "# Selects the HSE03 hybrid function\n"
            tmp += "LHFCALC = .TRUE. ; HFSCREEN = 0.3 ;\n"
            tmp += "ALGO = D ; TIME = 0.4 ; LSUBROT = .TRUE.\n"
            tmp += "\n"
        else:
            tmp += "## Selects the HSE03 hybrid function\n"
            tmp += "#LHFCALC = .TRUE. ; HFSCREEN = 0.3 ;\n"
            tmp += "#ALGO = D ; TIME = 0.4 ; LSUBROT = .TRUE.\n"
            tmp += "\n"
        #
        if self.dfttype == "HSE06":
            tmp += "# Selects the HSE06 hybrid function\n"
            tmp += "LHFCALC = .TRUE. ; HFSCREEN = 0.2 ;\n"
            tmp += "ALGO = D ; TIME = 0.4 ; LSUBROT = .TRUE.\n"
            tmp += "\n"
        else:
            tmp += "## Selects the HSE06 hybrid function\n"
            tmp += "#LHFCALC = .TRUE. ; HFSCREEN = 0.2 ;\n"
            tmp += "#ALGO = D ; TIME = 0.4 ; LSUBROT = .TRUE.\n"
            tmp += "\n"
        #
        if self.dfttype == "PBE0":
            tmp += "# Selects the PBE0  hybrid function\n"
            tmp += "LHFCALC = .TRUE. ;\n"
            tmp += "ALGO = D ; TIME = 0.4 ; LSUBROT = .TRUE.\n"
            tmp += "\n"
        else:
            tmp += "## Selects the PBE0  hybrid function\n"
            tmp += "#LHFCALC = .TRUE. ;\n"
            tmp += "#ALGO = D ; TIME = 0.4 ; LSUBROT = .TRUE.\n"
            tmp += "\n"
        #
        tmp += "\n"
        if self.runtype == "gw":
            tmp += "#### GW settings #####\n"
            tmp += "ISMEAR =  0\n"
            tmp += "ALGO = GW0\n"
            tmp += "NELM = 1\n"
            tmp += "PRECFOCK = Fast \n"
            tmp += "NOMEGA = 200\n"
            tmp += "NOMEGA = 200\n"
            tmp += "\n"
        else:
            tmp += "##### GW settings #####\n"
            tmp += "#ISMEAR =  0\n"
            tmp += "#ALGO = GW0\n"
            tmp += "#NELM = 1\n"
            tmp += "#PRECFOCK = Fast \n"
            tmp += "#NOMEGA = 200\n"
            tmp += "#NOMEGA = 200\n"
            tmp += "\n"
        #
        if self.runtype == "gw0+wannier":
            tmp += "#### GW0+Wannier settings #####\n"
            tmp += "ALGO = GW0\n"
            tmp += "LSPECTRAL = .TRUE.\n"
            tmp += "NOMEGA = 50\n"
            tmp += "LWANNIER90_RUN = .TRUE.\n"
            tmp += "\n"
        else:
            tmp += "##### GW0+Wannier settings #####\n"
            tmp += "#ALGO = GW0\n"
            tmp += "#LSPECTRAL = .TRUE.\n"
            tmp += "#NOMEGA = 50\n"
            tmp += "#LWANNIER90_RUN = .TRUE.\n"
            tmp += "\n"
        #
        tmp += "\n"
        if self.runtype == "nonscf+dos":
            tmp += "#### DOS settings #####\n"
            tmp += "ALGO = NONE\n"
            tmp += "NELM = 1\n"
            tmp += "ISMEAR = -5\n"
            tmp += "NSW = 0\n"
            tmp += "IBRION = -1\n"
            tmp += "ICHARGE = 11\n"
            tmp += "LORBIT = 11 # l-m decomposed DOS\n"
            tmp += "LWAVE = .FALSE. # do not overwrite WAVECAR\n"
            tmp += "LCHARG = .FALSE. # do not overwrite CHGCAR\n"
            tmp += "\n"
            tmp += "DOS plot settings\n"
            tmp += "EMIN = -20.0\n"
            tmp += "EMAX =  20.0\n"
            tmp += "NEDOS = 1082\n"
            tmp += "\n"
        else:
            tmp += "##### DOS settings #####\n"
            tmp += "#ALGO = NONE\n"
            tmp += "#NELM = 1\n"
            tmp += "#ISMEAR = -5\n"
            tmp += "#NSW = 0\n"
            tmp += "#IBRION = -1\n"
            tmp += "#ICHARGE = 11\n"
            tmp += "#LORBIT = 11 # l-m decomposed DOS\n"
            tmp += "#LWAVE = .FALSE. # do not overwrite WAVECAR\n"
            tmp += "#LCHARG = .FALSE. # do not overwrite CHGCAR\n"
            tmp += "\n"
            tmp += "#DOS plot settings\n"
            tmp += "#EMIN = -20.0\n"
            tmp += "#EMAX =  20.0\n"
            tmp += "#NEDOS = 1082\n"
            tmp += "\n"
        #
        tmp += "\n"
        if self.runtype == "nonscf+band":
            tmp += "#### band dispersion settings #####\n"
            tmp += "ALGO = NONE\n"
            tmp += "NELM = 1\n"
            tmp += "ISMEAR = 0\n" 
            tmp += "LWAVE = .FALSE. # do not overwrite WAVECAR\n"
            tmp += "LCHARG = .FALSE. # do not overwrite CHGCAR\n"
            tmp += "LWANNIER90_RUN = .TRUE. # run wannier90 in library mode\n"
            tmp += "\n"
        else:
            tmp += "##### band dispersion settings #####\n"
            tmp += "#ALGO = NONE\n"
            tmp += "#NELM = 1\n"
            tmp += "#ISMEAR = 0\n" 
            tmp += "#LWAVE = .FALSE. # do not overwrite WAVECAR\n"
            tmp += "#LCHARG = .FALSE. # do not overwrite CHGCAR\n"
            tmp += "#LWANNIER90_RUN = .TRUE. # run wannier90 in library mode\n"
            tmp += "\n"
        #
        tmp += "\n"
        if self.runtype == "phonopy":
            tmp  = "#### phonopy or phono3py settings #####\n"
            tmp += "PREC = Accurate\n"
            tmp += "ENCUT = 500\n"
            tmp += "IBRION = 8\n"
            tmp += "EDIFF = 1.0e-08\n"
            tmp += "IALGO = 38\n"
            tmp += "ISMEAR = 0; SIGMA = 0.1\n"
            tmp += "LREAL = .FALSE.\n"
            tmp += "ADDGRID = .TRUE.\n"
            tmp += "LWAVE = .FALSE.\n"
            tmp += "LCHARG = .FALSE.\n"
            tmp += "\n"
        else:
            tmp += "##### phonopy or phono3py settings #####\n"
            tmp += "#PREC = Accurate\n"
            tmp += "#ENCUT = 500\n"
            tmp += "#IBRION = 8\n"
            tmp += "#EDIFF = 1.0e-08\n"
            tmp += "#IALGO = 38\n"
            tmp += "#ISMEAR = 0; SIGMA = 0.1\n"
            tmp += "#LREAL = .FALSE.\n"
            tmp += "#ADDGRID = .TRUE.\n"
            tmp += "#LWAVE = .FALSE.\n"
            tmp += "#LCHARG = .FALSE.\n"
            tmp += "\n"
        #
        return tmp
        

################################################################################################
# EMTO 
class KFCDFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the kfcd program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.jobnam = "default"
        self.kstrjobnam = "default"
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        filestring = ""
        tmpstring = "KFCD      MSGL..=  0"
        tmpstring = tmpstring.ljust(25)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam+"\n"
        filestring += tmpstring
        tmpstring = "STRNAM...="+self.kstrjobnam+"\n"
        filestring += tmpstring
        filestring += "DIR001=../kstr/smx/\n"
        filestring += "DIR002=../kgrn/chd/\n"
        filestring += "DIR003=../shape/shp/\n"
        filestring += "DIR004=../bmdl/mdl/\n"
        filestring += "DIR006=./\n"
        filestring += "Lmaxs.= 30 NTH..= 41 NFI..= 81 FPOT..= N\n"
        filestring += "OVCOR.=  Y UBG..=  N NPRN.=  0\n"
        return filestring

class KGRNFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the kgrn program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.jobnam = "default"
        self.kstrjobnam = "default"
        # To be put on the first line
        self.programdoc = ""
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.latticenr = 14
    def __str__(self):
        ed = ElementData()
        filestring = ""
        tmpstring = "KGRN"
        tmpstring = tmpstring.ljust(25)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM="+self.jobnam+"\n"
        filestring += tmpstring
        filestring += "STRT..=  A MSGL.=  0 EXPAN.= S FCD..=  Y FUNC..= SCA\n"
        tmpstring = "FOR001=../kstr/smx/"+self.kstrjobnam+".tfh\n"
        tmpstring += "FOR001=../kstr/smx/"+self.kstrjobnam+"10.tfh\n"
        filestring += tmpstring
        filestring += "DIR002=pot/\n"
        filestring += "DIR003=pot/\n"
        tmpstring = "FOR004=../bmdl/mdl/"+self.kstrjobnam+".mdl\n"
        filestring += tmpstring
        filestring += "DIR006=\n"
        filestring += "DIR009=pot/\n"
        filestring += "DIR010=chd/\n"
        # Use environment variable TMPDIR if possible
        tmpstring = "DIR011="
        if "TMPDIR" in os.environ:
            tmpstring += os.environ["TMPDIR"]
            # Make sure the string will end with a single /
            tmpstring = tmpstring.rstrip("/")
        else:
            # ...else check for /tmp
            if os.path.isdir("/tmp"):
                tmpstring += "/tmp"
            else:
                # ...and last resort is ./
                tmpstring += "."
        tmpstring += "/\n"
        filestring += tmpstring
        filestring += self.docstring.replace("\n"," ")+"\n"
        filestring += "Band: 10 lines\n"
        tmpstring = "NITER.= 50 NLIN.= 31 NPRN.=  0 NCPA.= 20 NT...=%3i"%len(self.cell.atomdata)+" MNTA.="
        # Work out maximal number of species occupying a site
        mnta = 1
        for a in self.cell.atomdata:
            for b in a:
                mnta = max(mnta,len(b.species))
        tmpstring += "%3i"%mnta+"\n"
        filestring += tmpstring
        filestring += "MODE..= 3D FRC..=  N DOS..=  N OPS..=  N AFM..=  P CRT..=  M\n"
        filestring += "Lmaxh.=  8 Lmaxt=  4 NFI..= 31 FIXG.=  2 SHF..=  0 SOFC.=  N\n"
        # Choose brillouin zone by lattice type
        # Output the smallest allowed n for each direction in this lattice type
        if self.latticenr == 1:
            nkx = 0
            nky = 2
            nkz = 0
        elif self.latticenr == 2:
            nkx = 0
            nky = 5
            nkz = 0
        elif self.latticenr == 3:
            nkx = 0
            nky = 3
            nkz = 0
        elif self.latticenr == 4:
            nkx = 0
            nky = 3
            nkz = 2
        elif self.latticenr == 5:
            nkx = 0
            nky = 2
            nkz = 2
        elif self.latticenr == 6:
            nkx = 0
            nky = 3
            nkz = 2
        elif self.latticenr == 7:
            nkx = 0
            nky = 3
            nkz = 3
        elif self.latticenr == 8:
            nkx = 2
            nky = 2
            nkz = 1
        elif self.latticenr == 9:
            nkx = 2
            nky = 2
            nkz = 1
        elif self.latticenr == 10:
            nkx = 2
            nky = 2
            nkz = 1
        elif self.latticenr == 11:
            nkx = 1
            nky = 1
            nkz = 1
        elif self.latticenr == 12:
            nkx = 1
            nky = 2
            nkz = 1
        elif self.latticenr == 13:
            nkx = 3
            nky = 3
            nkz = 0
        else:
            nkx = 2
            nky = 2
            nkz = 2
        filestring += "KMSH...= G IBZ..= %2i NKX..= %2i NKY..= %2i NKZ..= %2i FBZ..=  N\n"%(self.latticenr,nkx,nky,nkz)
        filestring += "KMSH2..= G IBZ2.=  1 NKX2.=  4 NKY2.=  0 NKZ2.= 51\n"
        filestring += "ZMSH...= C NZ1..= 16 NZ2..= 16 NZ3..=  8 NRES.=  4 NZD.= 500\n"
        filestring += "DEPTH..=  1.500 IMAGZ.=  0.020 EPS...=  0.200 ELIM..= -1.000\n"
        filestring += "AMIX...=  0.100 EFMIX.=  1.000 VMTZ..=  0.000 MMOM..=  0.000\n"
        filestring += "TOLE...= 1.d-05 TOLEF.= 1.d-05 TOLCPA= 1.d-05 TFERMI=  500.0 (K)\n"
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        # average wigner-seitz radius
        volume = abs(det3(self.cell.latticevectors))
        wsr = self.cell.lengthscale * 3*volume/(nosites * 4 * pi)**third
        filestring += "SWS......=%8f   NSWS.=  1 DSWS..=   0.05 ALPCPA= 0.9020\n"%wsr
        filestring += "Setup: 2 + NQ*NS lines\n"
        filestring += "EFGS...=  0.000 HX....=  0.100 NX...= 11 NZ0..=  6 STMP..= Y\n"
        # atom info
        filestring += "Symb   IQ IT ITA NZ  CONC   Sm(s)  S(ws) WS(wst) QTR SPLT\n"
        iq = 1
        it = 1
        for a in self.cell.atomdata:
            ita = 1
            # THIS MAKES ASSUMPTIONS ABOUT THE ORDERING OF ATOMDATA
            # But we're OK for all orderings implemented so far
            for comp in a[0].spcstring().split("/"):
                for b in a:
                    if comp in b.species:
                        tmpstring = comp.ljust(4)+"  "+"%3i%3i%3i"%(iq,it,ita)
                        tmpstring += "%4i"%ed.elementnr[comp]
                        tmpstring += "%7.3f%7.3f%7.3f%7.3f"%(a[0].species[comp],1,1,1)
                        tmpstring += "%5.2f%5.2f\n"%(0,0)
                        filestring += tmpstring
                        iq += 1
                ita += 1
                iq -= len(a)
            iq += len(a)
            it += 1
        filestring += "Atom:  4 lines + NT*NTA*6 lines\n"
        filestring += "IEX...=  4 NP..= 251 NES..= 15 NITER=100 IWAT.=  0 NPRNA=  0\n"
        filestring += "VMIX.....=  0.300000 RWAT....=  3.500000 RMAX....= 20.000000\n"
        filestring += "DX.......=  0.030000 DR1.....=  0.002000 TEST....=  1.00E-12\n"
        filestring += "TESTE....=  1.00E-12 TESTY...=  1.00E-12 TESTV...=  1.00E-12\n"
        for a in self.cell.atomdata:
            for comp in a[0].species:
                filestring += comp+"\n"
                try:
                    filestring += ed.emtoelements[comp]
                except KeyError:
                    filestring += "\n\n\n\n\n"
        return filestring

class ShapeFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the shape program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.jobnam = "default"
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        filestring = ""
        tmpstring = "SHAPE     HP......=N"
        tmpstring = tmpstring.ljust(25)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(10)+" MSGL.=  1\n"
        filestring += tmpstring
        filestring += "FOR001=../kstr/smx/"+self.jobnam+".tfh\n"
        filestring += "DIR002=shp/\n"
        filestring += "DIR006=./\n"
        filestring += "Lmax..= 30 NSR..=129 NFI..= 11\n"
        filestring += "NPRN..=  0 IVEF.=  3\n"
        return filestring

class BMDLFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the bmdl program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.jobnam = "default"
        self.latticenr = 14
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        lv = self.cell.latticevectors
        ed = ElementData()
        filestring = ""
        tmpstring = "BMDL      HP......=N"
        tmpstring = tmpstring.ljust(25)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(10)+" MSGL.=  1 NPRN.=  0\n"
        filestring += tmpstring
        filestring += "DIR001=mdl/\n"
        filestring += "DIR006=./\n"
        filestring += "Madelung potential, "+self.docstring.replace("\n"," ")+"\n"
        filestring += "NL.....= 7\n"
        filestring += "LAMDA....=    2.5000 AMAX....=    4.5000 BMAX....=    4.5000\n"
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        if self.latticenr == 0:
            tmpstring = "NQ3...=%3i LAT...= 0 IPRIM.= 0 NGHBP.=13 NQR2..= 0\n" % (nosites,self.latticenr)
        else:
            tmpstring = "NQ3...=%3i LAT...=%2i IPRIM.= 1 NGHBP.=13 NQR2..= 0\n" % (nosites,self.latticenr)
        filestring += tmpstring
        boa = self.b/self.a
        coa = self.c/self.a
        filestring += "A........= 1.0000000 B.......=%10f C.......=%10f\n"%(boa,coa)
        tmpstring = ""
        if self.latticenr == 0:
            for i in range(3):
                tmpstring += "BSX......=%10f BSY.....=%10f BSZ.....=%10f\n" % (lv[i][0],lv[i][1],lv[i][2])
        else:
            tmpstring +=  "ALPHA....=%10f BETA....=%10f GAMMA...=%10f\n" % (self.alpha, self.beta, self.gamma)
        filestring += tmpstring
        for a in self.cell.atomdata:
            for b in a:
                v = mvmult3(lv,b.position)
                filestring += "QX(IQ)...=%10f QY......=%10f QZ......=%10f" % (v[0],v[1],v[2])
                filestring += "      "+b.spcstring()+"\n"
        return filestring

class KSTRFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the kstr program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.jobnam = "default"
        self.latticenr = 14
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.hardsphere = 0.67
        self.iprim = 0
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        lv = self.cell.latticevectors
        ed = ElementData()
        filestring = ""
        tmpstring = "KSTR      HP......=N"
        tmpstring = tmpstring.ljust(25)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(10)+" MSGL.=  1 MODE...=B STORE..=Y HIGH...=Y\n"
        filestring += tmpstring
        filestring += "DIR001=smx/\n"
        filestring += "DIR006=./\n"
        filestring += "Slope matrices, "+self.docstring.replace("\n"," ")+"\n"
        # NL = maximal l from element blocks
        maxl = 1
        for a in self.cell.atomdata:
            for b in a:
                for i in b.species:
                    if ed.elementblock[i] == "p":
                        maxl = max(maxl,2)
                    elif ed.elementblock[i] == "d":
                        maxl = max(maxl,3)
                    elif ed.elementblock[i] == "f":
                        maxl = max(maxl,4)
        tmpstring = "NL.....= %1i NLH...=11 NLW...= 9 NDER..= 6 ITRANS= 3 NPRN..= 0\n" % maxl
        filestring += tmpstring
        # Setting the real space summation cutoff to 4.5*(wigner-seitz radius)
        volume = abs(det3(lv))
        wsr = (3*volume/(self.cell.natoms() * 4 * pi))**third
        tmpstring = "(K*W)^2..=  0.000000 DMAX....=%10f RWATS...=      0.10\n" % (wsr*4.5)
        filestring += tmpstring
        tmpstring = "NQ3...=%3i LAT...=%2i IPRIM.=%2i NGHBP.=13 NQR2..= 0\n" % (self.cell.natoms(),self.latticenr,self.iprim)
        filestring += tmpstring
        boa = self.b/self.a
        coa = self.c/self.a
        filestring += "A........= 1.0000000 B.......=%10f C.......=%10f\n"%(boa,coa)
        tmpstring = ""
        if self.iprim == 0:
            for i in range(3):
                tmpstring += "BSX......=%10f BSY.....=%10f BSZ.....=%10f\n" % (lv[i][0],lv[i][1],lv[i][2])
        else:
            tmpstring +=  "ALPHA....=%10f BETA....=%10f GAMMA...=%10f\n" % (self.alpha, self.beta, self.gamma)
        filestring += tmpstring
        for a in self.cell.atomdata:
            for b in a:
                v = mvmult3(lv,b.position)
                filestring += "QX(IQ)...=%10f QY......=%10f QZ......=%10f" % (v[0],v[1],v[2])
                filestring += "      "+b.spcstring()+"\n"
        for i in range(self.cell.natoms()):
            filestring += "a/w(IQ)..="
            for i in range(4):
                filestring += "%5.2f"%self.hardsphere
            filestring += "\n"
        filestring += "LAMDA....=    2.5000 AMAX....=    4.5000 BMAX....=    4.5000\n"
        return filestring

################################################################################################
# SPRKKR
class XBandSysFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].sys file for the xband program
    and the method __str__ that outputs the contents of the .sys file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.jobnam = "default"
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.minangmom = None
        self.filename = ""
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        ed = ElementData()
        # First identify any site which is not filled (concentrations add to 1.0)
        # and fill up with vacuum sphere (Vc)
        self.cell.fill_out_empty(label='Vc')
        # docstring on a single line
        filestring = deletenewline(self.docstring,replace=" ")
        filestring += "\n"+self.filename+"\n"
        filestring += "xband-version\n"
        filestring += "5.0\n"
        # It would be really cool to support lower dimensions...one day.
        filestring += "dimension\n"
        filestring += "3D\n"
        filestring += "Bravais lattice\n"
        if self.cell.crystal_system() == 'triclinic':
            filestring += "1  triclinic   primitive      -1     C_i \n"
        elif self.cell.crystal_system() == 'monoclinic':
            if self.cell.spacegroupsetting == 'P':
                filestring += "2  monoclinic  primitive      2/m    C_2h\n"
            elif self.cell.spacegroupsetting == 'A' or self.cell.spacegroupsetting == 'B' or \
                     self.cell.spacegroupsetting == 'C':
                filestring += "3  monoclinic  primitive      2/m    C_2h\n"
            else:
                print "xband only knows primitive and base-centered monoclinic settings!"
                sys.exit(43)
        elif self.cell.crystal_system() == 'orthorhombic':
            if self.cell.spacegroupsetting == 'P':
                filestring += "4  orthorombic primitive      mmm    D_2h\n"
            elif self.cell.spacegroupsetting == 'A' or self.cell.spacegroupsetting == 'B' or \
                     self.cell.spacegroupsetting == 'C':
                filestring += "5  orthorombic body-centered  mmm    D_2h\n"
            elif self.cell.spacegroupsetting == "I":
                filestring += "6  orthorombic body-centered  mmm    D_2h\n"
            elif self.cell.spacegroupsetting == "F":
                filestring += "7  orthorombic face-centered  mmm    D_2h\n"
            else:
                print "xband does not know %1s centering of an orthorhombic cell."%self.cell.spacegroupsetting
                sys.exit(43)
        elif self.cell.crystal_system() == "tetragonal":
            if self.cell.spacegroupsetting == "P":
                filestring += "8  tetragonal  primitive      4/mmm  D_4h\n"
            elif self.cell.spacegroupsetting == "I":
                filestring += "9  tetragonal  body-centered  4/mmm  D_4h\n"
            else:
                print "xband only knows primitive and body-centered tetragonal settings!"
                sys.exit(43)
        elif self.cell.crystal_system() == "trigonal":
            filestring += "10 trigonal    primitive      -3m    D_3d\n"
        elif self.cell.crystal_system() == "hexagonal":
            filestring += "11 hexagonal   primitive      6/mmm  D_6h\n"
        elif self.cell.crystal_system() == "cubic":
            if self.cell.spacegroupsetting == "P":
                filestring += "12 cubic       primitive      m3m    O_h \n"
            elif self.cell.spacegroupsetting == "F":
                filestring += "13 cubic       face-centered  m3m    O_h \n"
            elif self.cell.spacegroupsetting == "I":
                filestring += "14 cubic       body-centered  m3m    O_h \n"
            else:
                print "xband does not know %1s centering of a cubic cell."%self.cell.spacegroupsetting
                sys.exit(43)
        filestring += "space group number (ITXC and AP)\n"
        filestring += "%5i%5i"%(self.cell.spacegroupnr,Number2AP[self.cell.spacegroupnr])+"\n"
        filestring += "structure type\n"
        filestring += "UNKNOWN\n"
        filestring += "lattice parameter A  [a.u.]\n"
        filestring += "%18.12f\n"%self.cell.lengthscale
        filestring += "ratio of lattice parameters  b/a  c/a\n"
        filestring += "%18.12f%18.12f\n"%(self.cell.boa,self.cell.coa)
        filestring += "lattice parameters  a b c  [a.u.]\n"
        a = self.cell.lengthscale
        b = self.cell.b * self.cell.lengthscale / self.cell.a
        c = self.cell.c * self.cell.lengthscale / self.cell.a
        filestring += "%18.12f%18.12f%18.12f\n"%(a,b,c)
        filestring += "lattice angles  alpha beta gamma  [deg]\n"
        filestring += "%18.12f%18.12f%18.12f\n"%(self.cell.alpha,self.cell.beta,self.cell.gamma)
        filestring += "primitive vectors     (cart. coord.) [A]\n"
        for vec in self.cell.latticevectors:
            for p in vec:
                filestring += "%18.12f"%p
            filestring += "\n"
        # Get number of sites and fill out with empty spheres if the sites are not fully filled
        filestring += "number of sites NQ\n"
        nq = 0
        for a in self.cell.atomdata:
            nq += len(a)
        self.cell.fill_out_empty(label="Vc")
        filestring += "%3i\n"%nq
        filestring += " IQ ICL     basis vectors     (cart. coord.) [A]                      RWS [a.u.]  NLQ  NOQ ITOQ\n"
        # Average Wigner-Seitz radius
        rws = pow(3*self.cell.volume()/(4*pi*len(self.cell.atomset)),1.0/3.0)*self.cell.lengthscale
        iq = 0
        icl = 0
        itoq = 0
        for a in self.cell.atomdata:
            icl += 1
            itoqs = []
            for sp in a[0].species:
                itoq += 1
                itoqs.append(itoq)
            for b in a:
                iq += 1
                if self.minangmom:
                    angmom = max(max([ed.angularmomentum[ed.elementblock[spcs]] for spcs in b.species])+1,self.minangmom)
                else:
                    angmom = max([ed.angularmomentum[ed.elementblock[spcs]] for spcs in b.species])+1
                v = mvmult3(self.cell.latticevectors,b.position)
                filestring += "%3i%4i%18.12f%18.12f%18.12f  %18.12f%4i%5i "%(iq,icl,v[0],v[1],v[2],rws,angmom,len(a[0].species))
                for i in itoqs:
                    filestring += "%3i"%i
                filestring += "\n"
        filestring += "number of sites classes NCL\n"
        filestring += "%3i\n"%len(self.cell.atomdata)
        filestring += "ICL WYCK NQCL IQECL (equivalent sites)\n"
        iq = 0
        icl = 0
        for a in self.cell.atomdata:
            icl += 1
            filestring += "%3i   %1s%5i"%(icl,'-',len(a))
            for b in a:
                iq += 1
                filestring += "%3i"%iq
            filestring += "\n"
        filestring += "number of atom types NT\n"
        nt = 0
        for a in self.cell.atomdata:
            nt += len(a[0].species)
        filestring += "%3i\n"%nt
        filestring += " IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n"
        iq = 0
        it = 0
        for a in self.cell.atomdata:
            corr = 0
            for sp,conc in a[0].species.iteritems():
                it += 1
                filestring += " %2i%4i  %8s%5i%6.3f"%(it,ed.elementnr[sp],sp,len(a),conc)
                iq -= corr*len(a)
                for b in a:
                    iq += 1
                    filestring += "%3i"%iq
                corr = 1
                filestring += "\n"
        return filestring

class SPCFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the SPC program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.jobnam = "default"
        self.latticenr = 1
        self.compoundname = ""
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.iprim = 0
        #
        self.supercell=[1,1,1]
        self.pairs = 6
        self.triplets = 4
        self.quartets = 2
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        import datetime
        now = datetime.datetime.now()
        datestring = now.strftime("%d %b %y")
        filestring = "SPC       HP......=N                                         "+datestring+"\n"
        filestring += "JOBNAM...="+self.jobnam.ljust(10)+" MSGL.=  1 \n"
        if os.path.isdir('spc'):
            dirname = 'spc/'
        else:
            dirname = './'
        filestring += "FOR001="+dirname+"\n"
        filestring += "FOR002="+dirname+"\n"
        filestring += "FOR004="+dirname+"\n"
        filestring += "FOR006=\n"
        filestring += "FOR008="+dirname+"\n"
        filestring += "FOR009="+dirname+"\n"
        filestring += self.docstring
        filestring += "Supercell, "+self.compoundname+"\n"
        filestring += "NPRN..=  0 TEST.=  0 NCOL.=  0 STAT.=  1 TCLIM=  0 nsho.= 10\n"
        filestring += "NQ3...=%3i LAT..=%3i IPRIM=%3i HIGH.=  0 NSHC.= 10 NL...=  4 NLH..=  7\n"%(self.cell.natoms(),self.latticenr,self.iprim)
        self.a = self.cell.latticevectors[0].length()*self.cell.lengthscale
        self.b = self.cell.latticevectors[1].length()*self.cell.lengthscale
        self.c = self.cell.latticevectors[2].length()*self.cell.lengthscale
        # Renormalized lattice vectors
        lv = []
        for i in range(3):
            lv.append([])
            for j in range(3):
                lv[i].append(self.cell.latticevectors[i][j]*self.cell.lengthscale/self.a)
        filestring += "A........=%10f B.......=%10f C.......=%10f\n"%(self.a,self.b,self.c)
        tmpstring = ""
        if self.iprim == 0:
            for i in range(3):
                tmpstring += "BSX......=%10f BSY.....=%10f BSZ.....=%10f\n" % (lv[i][0],lv[i][1],lv[i][2])
        else:
            tmpstring +=  "ALPHA....=%10f BETA....=%10f GAMMA...=%10f\n" % (self.alpha, self.beta, self.gamma)
        filestring += tmpstring
        for a in self.cell.atomdata:
            for b in a:
                v = mvmult3(lv,b.position)
                filestring += "QX.......=%10f QY......=%10f QZ......=%10f" % (v[0],v[1],v[2])
                filestring += "      "+b.spcstring()+"\n"
        filestring += "LAMDA....=    2.5000 AMAX....=    4.5000 BMAX....=    4.5000\n"
        filestring += "Size of the super cell\n"
        filestring += "NA.......=%4i NB.......=%4i NC.......=%4i  Dmax     4.5\n"%(self.supercell[0],self.supercell[1],self.supercell[2])
        filestring += "NSDC.....=  20 NSDS.....=   1 NSDM.....=   1\n"
        filestring += "NMAXMX...=   3 TMLIM....= 1.0\n"
        filestring += "NT.......=%4i\n"%(len(self.cell.atomdata))
        filestring += "NTA(IQ)..="
        i = 0
        nat = 0
        for a in self.cell.atomdata:
            i += 1
            for b in a:
                nat += 1
                filestring += "%4i"%i
                if nat%15 == 0:
                    filestring += "\n          "
        filestring = filestring.rstrip(" ")
        if self.cell.natoms()%15 != 0:
            filestring += "\n"
        #
        filestring += "NTO......=%4i\n"%(len(self.cell.atomdata))
        filestring += "NTAO(IQ).="
        i = 0
        nat = 0
        for a in self.cell.atomdata:
            i += 1
            for b in a:
                nat += 1
                filestring += "%4i"%i
                if nat%15 == 0:
                    filestring += "\n          "
        filestring = filestring.rstrip(" ")
        if self.cell.natoms()%15 != 0:
            filestring += "\n"
        #
        filestring += "NQ3O.....=%4i\n"%(len(self.cell.atomdata))
        filestring += "IQO(IQ)..="
        i = 0
        nat = 0
        for a in self.cell.atomdata:
            i += 1
            for b in a:
                nat += 1
                filestring += "%4i"%i
                if nat%15 == 0:
                    filestring += "\n          "
        filestring = filestring.rstrip(" ")
        if self.cell.natoms()%15 != 0:
            filestring += "\n"
        conc = []
        ascii = string.ascii_uppercase+string.ascii_lowercase
        i = 0
        for b in self.cell.atomdata:
            natom = 0
            conc.append([])
            for a in self.cell.atomdata:
                for k,v in a[0].species.iteritems():
                    if a == b:
                        conc[i].append((ascii[natom],v))
                    else:
                        conc[i].append((ascii[natom],0.0))
                    natom += 1
            i += 1
        filestring += "NATOM....=%4i\n"%natom
        filestring += "SMB(IAT).="
        for a in conc[0]:
            filestring += "%4s"%a[0]
        filestring += "\n"
        filestring += "Concentrations on sublattices:\n"
        for a in conc:
            for c in a:
                filestring += "%f "%c[1]
            filestring += "\n"
        filestring += "Correlation functions and weights for each pairs of elem. (A-B, A-C, ... )\n"
        i = 0
        for a in self.cell.atomdata:
            i += 1
            filestring += "Sublattice\n"
            filestring += "%i\n"%i
            if len(a[0].species) == 1:
                continue
            filestring += "nc2  r_max\n"
            if len(a[0].species) == 1:
                filestring += "0     3.0\n"
                filestring += "i    alpha           weight\n"
            else:
                filestring += "%i     3.0\n"%self.pairs
                filestring += "i    alpha           weight\n"
                for p in range(self.pairs):
                    if p < 3:
                        filestring += "%i     0.0               1.0\n"%(p+1)
                    else:
                        filestring += "%i     0.0               0.0\n"%(p+1)
            filestring += "nc3\n"
            filestring += "0\n"
            filestring += "i   i1  i2  i3           <sss>          weight\n"
            filestring += "nc4\n"
            filestring += "0\n"
            filestring += "i   i1 i2 i3 i4 i5 i6    <ssss>         weight\n"
            filestring += "T_i,   T_f,  delt_T\n"
            filestring += "10.0   0.0   1.0\n"
            filestring += "100                   nsteps\n"
        return filestring

################################################################################################
# MOPAC FILE
class MOPACFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting a MOPAC file
    and the method __str__ that outputs the contents of the MOPAC file as a string.
    """
    def __init__(self,crystalstructure,string,setupall=False,firstline="",secondline="",thirdline="",freeze=-1):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit(newunit="angstrom")
        self.setupall = setupall
        self.firstline = firstline
        self.secondline = secondline
        self.thirdline = thirdline
        self.freeze = freeze
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        t = True
        for s in tmpstrings:
            t = t and (s[0] == '*')
        if t:
            self.docstring += "\n"
        else:
            self.docstring = ""
            for string in tmpstrings:
                string = string.lstrip("*")
                string = "*"+string+"\n"
                self.docstring += string
    def __str__(self):
        filestring = self.docstring
        if self.setupall:
            if self.firstline == "":
                filestring += " BZ \n"
            else:
                filestring += self.firstline.rstrip("\n")+"\n"
            filestring += self.secondline.rstrip("\n")+"\n"
            filestring += self.thirdline.rstrip("\n")+"\n"
        # Set up lattice vectors
        lv = []
        for i in range(3):
            lv.append(Vector([self.cell.lengthscale*self.cell.latticevectors[i][j] for j in range(3)]))
        # Print sites
        if self.freeze == 0:
            freezestring = " 0"
        elif self.freeze == 1:
            freezestring = " 1"
        else:
            freezestring = ""
        for a in self.cell.atomdata:
            for b in a:
                t = Vector(mvmult3(lv,b.position))
                filestring += str(b).split()[0]+"  "+str(t)+freezestring+"\n"
        # Print lattice vectors
        for l in lv:
            filestring += "Tv  "+str(l)+"\n"
        return filestring
################################################################################################
# AkaiKKR
class AkaiKKRFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].in file for the xband program
    and the method __str__ that outputs the contents of the .in file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.jobnam = "default"
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.minangmom = None
        self.filename = ""
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        ed = ElementData()
        # First identify any site which is not filled (concentrations add to 1.0)
        # and fill up with vacuum sphere (Vc)
        self.cell.fill_out_empty(label='Vc')
        # docstring on a single line
        filestring = deletenewline(self.docstring,replace=" ")
        filestring += "\n"+self.filename+"\n"
        # It would be really cool to support lower dimensions...one day.
        filestring += "Bravais lattice\n"
        if self.cell.crystal_system() == 'triclinic':
            filestring += "1  triclinic   primitive      -1     C_i \n"
            brvtyp = "trc"
        elif self.cell.crystal_system() == 'monoclinic':
            if self.cell.spacegroupsetting == 'P':
                filestring += "2  monoclinic  primitive      2/m    C_2h\n"
                brvtyp = "sm"
            elif self.cell.spacegroupsetting == 'A' or self.cell.spacegroupsetting == 'B' or \
                     self.cell.spacegroupsetting == 'C':
                filestring += "3  monoclinic  primitive      2/m    C_2h\n"
                brvtyp = "bsm"
            else:
                print "xband only knows primitive and base-centered monoclinic settings!"
                sys.exit(43)
        elif self.cell.crystal_system() == 'orthorhombic':
            if self.cell.spacegroupsetting == 'P':
                filestring += "4  orthorombic primitive      mmm    D_2h\n"
                brvtyp = "so"
            elif self.cell.spacegroupsetting == 'A' or self.cell.spacegroupsetting == 'B' or \
                     self.cell.spacegroupsetting == 'C':
                filestring += "5  orthorombic body-centered  mmm    D_2h\n"
                brvtyp = "bso"
            elif self.cell.spacegroupsetting == "I":
                filestring += "6  orthorombic body-centered  mmm    D_2h\n"
                brvtyp = "bco"
            elif self.cell.spacegroupsetting == "F":
                filestring += "7  orthorombic face-centered  mmm    D_2h\n"
                brvtyp = "fco"
            else:
                print "xband does not know %1s centering of an orthorhombic cell."%self.cell.spacegroupsetting
                sys.exit(43)
        elif self.cell.crystal_system() == "tetragonal":
            if self.cell.spacegroupsetting == "P":
                filestring += "8  tetragonal  primitive      4/mmm  D_4h\n"
                brvtyp = "st"
            elif self.cell.spacegroupsetting == "I":
                filestring += "9  tetragonal  body-centered  4/mmm  D_4h\n"
                brvtyp = "bct"
            else:
                print "xband only knows primitive and body-centered tetragonal settings!"
                sys.exit(43)
        elif self.cell.crystal_system() == "trigonal":
            filestring += "10 trigonal    primitive      -3m    D_3d\n"
            brvtyp = "trg"
        elif self.cell.crystal_system() == "hexagonal":
            filestring += "11 hexagonal   primitive      6/mmm  D_6h\n"
            brvtyp = "hcp"
        elif self.cell.crystal_system() == "cubic":
            if self.cell.spacegroupsetting == "P":
                filestring += "12 cubic       primitive      m3m    O_h \n"
                brvtyp = "sc"
            elif self.cell.spacegroupsetting == "F":
                filestring += "13 cubic       face-centered  m3m    O_h \n"
                brvtyp = "fcc"
            elif self.cell.spacegroupsetting == "I":
                filestring += "14 cubic       body-centered  m3m    O_h \n"
                brvtyp = "bcc"
            else:
                print "xband does not know %1s centering of a cubic cell."%self.cell.spacegroupsetting
                sys.exit(43)
        #
        # options
        run_level = self.run_level
        sbrvtyp = self.sbrvtyp
        collect_atoms = self.collect_atoms
        #
        # AkaiKKR file
        run_type = "go"
        bzqlty = 4
        filestring = "c "+self.docstring+"\n"
        filestring += "\n"
        filestring += "c------------------------------------------------------------\n"
        filestring += "     "+run_type+"   data/"+self.filename.replace(".in","")+"\n"
        filestring += "c------------------------------------------------------------\n"
        filestring += "c   brvtyp     a        c/a   b/a   alpha   beta   gamma\n"
        #
        if sbrvtyp == "brvtyp":
            #filestring += "lattice parameters  a b c  [a.u.]\n"
            a = self.cell.lengthscale
            b = self.cell.b * self.cell.lengthscale / self.cell.a
            c = self.cell.c * self.cell.lengthscale / self.cell.a
            filestring += "     %3s   %7.2f   %6.3f %6.3f   %5.2f   %5.2f   %5.2f\n"%(brvtyp,a,self.cell.coa,self.cell.boa,self.cell.alpha,self.cell.beta,self.cell.gamma)
        else:
             #filestring += "   aux\n"
             filestring += "   prv\n"
             #filestring += "primitive vectors     (cart. coord.) [A]\n"
             for vec in self.cell.latticevectors:
                 for p in vec:
                     filestring += "  %8.5f"%p
                 filestring += "\n"
             #filestring += "lattice parameter A  [a.u.]\n"
             filestring += " %10.5f\n"%(self.cell.lengthscale)
        #
        filestring += "c------------------------------------------------------------\n"
        filestring += "c   outtyp    bzqlty   maxitr   pmix\n"
        filestring += "    update    %3i       900     0.023\n"%bzqlty
        filestring += "c------------------------------------------------------------\n"
        #
        iq = 0
        it = 0
        iq_old = 0
        ncomp = []
        ncomp = 50*[1]
        type_ncomp = []
        type_ncomp = 50*["XX"]
        for a in self.cell.atomdata:
            corr = 0
            for sp,conc in a[0].species.iteritems():
                it += 1
                iq -= corr*len(a)
                for b in a:
                    iq += 1
                    if iq == iq_old:
                        ncomp[it] += 1
                        ncomp[it-1] = 0
                    iq_old = iq
                corr = 1
                type_ncomp[it] = sp
        it_max = it
        #
        filestring += "c    ntyp\n"
        ntyp = len(self.cell.atomdata)
        for i in range(it+1):
            natom_equ = 1
            for j in range(i+1,it+1):
                if type_ncomp[j] == type_ncomp[i] and ncomp[j] != -1:
                    natom_equ += 1
                    if collect_atoms == "yes":
                        ncomp[j] = -1
                        ntyp -= 1
                    else:
                        type_ncomp[j] += str(natom_equ)
                        type_ncomp[i] += "1"
        filestring += "    %3i\n"%ntyp                       
        #
        for i in range(it+1):
            if ncomp[it_max-i] == 0:
                ncomp[it_max-i] = ncomp[it_max-i+1]
                ncomp[it_max-i+1] = 0
                type_ncomp[it_max-i] += type_ncomp[it_max-i+1]
                type_ncomp[it_max-i+1] = " "
        for i in range(it+1):
            if type_ncomp[i] == " ":
                type_ncomp[i] = type_ncomp[i-1]
        filestring += "c------------------------------------------------------------\n"
        filestring += "c   type    ncmp    rmt    field   mxl  anclr   conc\n"
        it = 0
        for a in self.cell.atomdata:
            corr = 0
            for sp,conc in a[0].species.iteritems():
                it += 1
                if ncomp[it] == -1:
                    pass 
                elif ncomp[it] >= 1:      
                    filestring += "%8s"%(type_ncomp[it])
                    filestring += "  %4i   "%(ncomp[it])
                    filestring += "    1      0.0     2  "
                    filestring += "%4i    %6.3f\n"%(ed.elementnr[sp],conc)
                else:
                    filestring += "                                       %4i    %6.3f\n"%(ed.elementnr[sp],conc)
        filestring += "c------------------------------------------------------------\n"
        filestring += "c    natm\n"
        nq = 0
        for a in self.cell.atomdata:
            nq += len(a)
        self.cell.fill_out_empty(label="Vc")
        filestring += "    %3i\n"%nq
        filestring += "c------------------------------------------------------------\n"
        filestring += "c   atmicx                        atmtyp\n" 
        # Average Wigner-Seitz radius
        rws = pow(3*self.cell.volume()/(4*pi*len(self.cell.atomset)),1.0/3.0)*self.cell.lengthscale
        itoq = 0
        for a in self.cell.atomdata:
            itoqs = []
            for sp in a[0].species:
                itoq += 1
                itoqs.append(itoq)
            for b in a:
                v = mvmult3(self.cell.latticevectors,b.position)
                filestring += "  %8.5fa %8.5fb  %8.5fc %8s\n"%(v[0],v[1],v[2],type_ncomp[itoq])
                #filestring += "  %8.5f  %8.5f  %8.5f %8s\n"%(v[0],v[1],v[2],sp)
        filestring += "c------------------------------------------------------------\n"
        #
        return filestring 


class AkaiKKR_input_File(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].in file for the xband program
    and the method __str__ that outputs the contents of the .in file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.jobnam = "default"
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.minangmom = None
        self.filename = ""
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        ed = ElementData()
        # First identify any site which is not filled (concentrations add to 1.0)
        # and fill up with vacuum sphere (Vc)
        self.cell.fill_out_empty(label='Vc')
        # docstring on a single line
        filestring = deletenewline(self.docstring,replace=" ")
        filestring += "\n"+self.filename+"\n"
        # It would be really cool to support lower dimensions...one day.
        filestring += "Bravais lattice\n"
        if self.cell.crystal_system() == 'triclinic':
            filestring += "1  triclinic   primitive      -1     C_i \n"
            brvtyp = "trc"
        elif self.cell.crystal_system() == 'monoclinic':
            if self.cell.spacegroupsetting == 'P':
                filestring += "2  monoclinic  primitive      2/m    C_2h\n"
                brvtyp = "sm"
            elif self.cell.spacegroupsetting == 'A' or self.cell.spacegroupsetting == 'B' or \
                     self.cell.spacegroupsetting == 'C':
                filestring += "3  monoclinic  primitive      2/m    C_2h\n"
                brvtyp = "bsm"
            else:
                print "xband only knows primitive and base-centered monoclinic settings!"
                sys.exit(43)
        elif self.cell.crystal_system() == 'orthorhombic':
            if self.cell.spacegroupsetting == 'P':
                filestring += "4  orthorombic primitive      mmm    D_2h\n"
                brvtyp = "so"
            elif self.cell.spacegroupsetting == 'A' or self.cell.spacegroupsetting == 'B' or \
                     self.cell.spacegroupsetting == 'C':
                filestring += "5  orthorombic body-centered  mmm    D_2h\n"
                brvtyp = "bso"
            elif self.cell.spacegroupsetting == "I":
                filestring += "6  orthorombic body-centered  mmm    D_2h\n"
                brvtyp = "bco"
            elif self.cell.spacegroupsetting == "F":
                filestring += "7  orthorombic face-centered  mmm    D_2h\n"
                brvtyp = "fco"
            else:
                print "xband does not know %1s centering of an orthorhombic cell."%self.cell.spacegroupsetting
                sys.exit(43)
        elif self.cell.crystal_system() == "tetragonal":
            if self.cell.spacegroupsetting == "P":
                filestring += "8  tetragonal  primitive      4/mmm  D_4h\n"
                brvtyp = "st"
            elif self.cell.spacegroupsetting == "I":
                filestring += "9  tetragonal  body-centered  4/mmm  D_4h\n"
                brvtyp = "bct"
            else:
                print "xband only knows primitive and body-centered tetragonal settings!"
                sys.exit(43)
        elif self.cell.crystal_system() == "trigonal":
            filestring += "10 trigonal    primitive      -3m    D_3d\n"
            brvtyp = "trg"
        elif self.cell.crystal_system() == "hexagonal":
            filestring += "11 hexagonal   primitive      6/mmm  D_6h\n"
            brvtyp = "hcp"
        elif self.cell.crystal_system() == "cubic":
            if self.cell.spacegroupsetting == "P":
                filestring += "12 cubic       primitive      m3m    O_h \n"
                brvtyp = "sc"
            elif self.cell.spacegroupsetting == "F":
                filestring += "13 cubic       face-centered  m3m    O_h \n"
                brvtyp = "fcc"
            elif self.cell.spacegroupsetting == "I":
                filestring += "14 cubic       body-centered  m3m    O_h \n"
                brvtyp = "bcc"
            else:
                print "xband does not know %1s centering of a cubic cell."%self.cell.spacegroupsetting
                sys.exit(43)
        #
        # options
        run_level = self.run_level
        sbrvtyp = self.sbrvtyp
        collect_atoms = self.collect_atoms
        #
        # AkaiKKR file
        filestring = ""
        for option in range(run_level):
            if option == 0:
                run_type = "go"
                bzqlty = 4
            if option == 1:
                run_type = "dos"
                bzqlty = 8
            if option == 2:
                run_type = "spc"
                bzqlty = 8
                #
            filestring += "c------------------------------------------------------------\n"
            filestring += "     "+run_type+"   data/"+self.filename.replace(".in","")+"\n"
            filestring += "c------------------------------------------------------------\n"
            filestring += "c   brvtyp     a        c/a   b/a   alpha   beta   gamma\n"
            #
            if sbrvtyp == "brvtyp":
                #filestring += "lattice parameters  a b c  [a.u.]\n"
                a = self.cell.lengthscale
                b = self.cell.b * self.cell.lengthscale / self.cell.a
                c = self.cell.c * self.cell.lengthscale / self.cell.a
                filestring += "     %3s   %7.2f   %6.3f %6.3f   %3.2f   %3.2f   %3.2f\n"%(brvtyp,a,self.cell.coa,self.cell.boa,self.cell.alpha,self.cell.beta,self.cell.gamma)
            else:
                 #filestring += "   aux\n"
                 filestring += "   prv\n"
                 #filestring += "primitive vectors     (cart. coord.) [A]\n"
                 for vec in self.cell.latticevectors:
                     for p in vec:
                         filestring += "  %8.5f"%p
                     filestring += "\n"
                 #filestring += "lattice parameter A  [a.u.]\n"
                 filestring += " %10.5f\n"%(self.cell.lengthscale)
            #
            filestring += "c------------------------------------------------------------\n"
            filestring += "c   edelt    ewidth    reltyp   sdftyp   magtyp   record\n"
            filestring += "    0.001     1.0       sra      pbe      mag      2nd\n"
            filestring += "c------------------------------------------------------------\n"
            filestring += "c   outtyp    bzqlty   maxitr   pmix\n"
            filestring += "    update    %3i       900     0.023\n"%bzqlty
            filestring += "c------------------------------------------------------------\n"
            #
            iq = 0
            it = 0
            iq_old = 0
            ncomp = []
            ncomp = 50*[1]
            type_ncomp = []
            type_ncomp = 50*["XX"]
            for a in self.cell.atomdata:
                corr = 0
                for sp,conc in a[0].species.iteritems():
                    it += 1
                    iq -= corr*len(a)
                    for b in a:
                        iq += 1
                        if iq == iq_old:
                            ncomp[it] += 1
                            ncomp[it-1] = 0
                        iq_old = iq
                    corr = 1
                    type_ncomp[it] = sp
            it_max = it
            #
            filestring += "c    ntyp\n"
            ntyp = len(self.cell.atomdata)
            for i in range(it+1):
                natom_equ = 1
                for j in range(i+1,it+1):
                    if type_ncomp[j] == type_ncomp[i] and ncomp[j] != -1:
                        natom_equ += 1
                        if collect_atoms == "yes":
                            ncomp[j] = -1
                            ntyp -= 1
                        else:
                            type_ncomp[j] += str(natom_equ)
                            type_ncomp[i] += "1"
            filestring += "    %3i\n"%ntyp                       
            #
            for i in range(it+1):
                if ncomp[it_max-i] == 0:
                    ncomp[it_max-i] = ncomp[it_max-i+1]
                    ncomp[it_max-i+1] = 0
                    type_ncomp[it_max-i] += type_ncomp[it_max-i+1]
                    type_ncomp[it_max-i+1] = " "
            for i in range(it+1):
                if type_ncomp[i] == " ":
                    type_ncomp[i] = type_ncomp[i-1]
            filestring += "c------------------------------------------------------------\n"
            filestring += "c   type    ncmp    rmt    field   mxl  anclr   conc\n"
            it = 0
            for a in self.cell.atomdata:
                corr = 0
                for sp,conc in a[0].species.iteritems():
                    it += 1
                    if ncomp[it] == -1:
                        pass 
                    elif ncomp[it] >= 1:      
                        filestring += "%8s"%(type_ncomp[it])
                        filestring += "  %4i   "%(ncomp[it])
                        filestring += "    1      0.0     2  "
                        filestring += "%4i    %6.3f\n"%(ed.elementnr[sp],conc)
                    else:
                        filestring += "                                       %4i    %6.3f\n"%(ed.elementnr[sp],conc)
            filestring += "c------------------------------------------------------------\n"
            filestring += "c    natm\n"
            nq = 0
            for a in self.cell.atomdata:
                nq += len(a)
            self.cell.fill_out_empty(label="Vc")
            filestring += "    %3i\n"%nq
            filestring += "c------------------------------------------------------------\n"
            filestring += "c   atmicx                        atmtyp\n" 
            # Average Wigner-Seitz radius
            rws = pow(3*self.cell.volume()/(4*pi*len(self.cell.atomset)),1.0/3.0)*self.cell.lengthscale
            itoq = 0
            for a in self.cell.atomdata:
                itoqs = []
                for sp in a[0].species:
                    itoq += 1
                    itoqs.append(itoq)
                for b in a:
                    v = mvmult3(self.cell.latticevectors,b.position)
                    filestring += "  %8.5fa %8.5fb  %8.5fc %8s\n"%(v[0],v[1],v[2],type_ncomp[itoq])
                    #filestring += "  %8.5f  %8.5f  %8.5f %8s\n"%(v[0],v[1],v[2],sp)
            filestring += "c------------------------------------------------------------\n"
            #
            if option == 2:
                if brvtyp == "fcc":
                    filestring += "c--- k point as WIEN2k_FCC.klist_band, Conventional basis\n"
                    filestring += "220\n"
                    filestring += "c--- W-point\n"
                    filestring += "1 0.5 0\n"
                    filestring += "c--- L-point\n"
                    filestring += "0.5 0.5 0.5\n"
                    filestring += "c---  Gamma-point\n"
                    filestring += "0 0 0\n"
                    filestring += "c--- X-point\n"
                    filestring += "1 0 0\n"
                    filestring += "c--- W-point\n"
                    filestring += "1 0.5 0\n"
                    filestring += "c--- K-point\n"
                    filestring += "0.75 0.75 0\n"
                    filestring += "c---\n"
                    filestring += "end\n"
                elif brvtyp == "bcc":
                    filestring += "c--- k point as WIEN2k_BCC.klist_band, Conventional basis\n"
                    filestring += "240\n"
                    filestring += "c--- Gamma-point\n"
                    filestring += "0 0 0\n"
                    filestring += "c--- H-point\n"
                    filestring += "0 1 0\n"
                    filestring += "c--- N-point\n"
                    filestring += "0.5 0.5 0\n"
                    filestring += "c--- Gamma-point\n"
                    filestring += "0 0 0\n"
                    filestring += "c--- P-point\n"
                    filestring += "0.5 0.5 0.5\n"
                    filestring += "c---\n"
                    filestring += "end\n"
                elif brvtyp == "hcp":
                    filestring += "c--- k point as WIEN2k_HCP.klist_band, Conventional basis\n"
                    filestring += "240\n"
                    filestring += "c--- Gamma-point\n"
                    filestring += "0 0 0\n"
                    filestring += "c--- M-point\n"
                    filestring += "0.5 0 0\n"
                    filestring += "c--- K-point\n"
                    filestring += "0.333 0.333 0\n"
                    filestring += "c--- Gamma-point\n"
                    filestring += "0 0 0\n"
                    filestring += "c--- A-point\n"
                    filestring += "0 0 0.5\n"
                    filestring += "c---\n"
                    filestring += "end\n"
                else:
                    filestring += "c--- k point as WIEN2k_SC.klist_band, Conventional basis\n"
                    filestring += "300\n"
                    filestring += "c--- R-point\n"
                    filestring += "0.5 0.5 0.5\n"
                    filestring += "c--- Gamma-point\n"
                    filestring += "0 0 0\n"
                    filestring += "c--- X-point\n"
                    filestring += "0.5 0 0\n"
                    filestring += "c--- M-point\n"
                    filestring += "0.5 0.5 0\n"
                    filestring += "c--- Gamma-point\n"
                    filestring += "0 0 0\n"
                    filestring += "c---\n"
                    filestring += "end\n"
                filestring += "c------------------------------------------------------------\n"
            #
        return filestring 
################################################################################################
class OpenMXFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an OpenMX run and the method
    __str__ that outputs the contents of a OpenMX input file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell = crystalstructure
        #self.cell.newunit("bohr")
        self.cell.newunit("angstrom")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        self.printbraces = False
    def __str__(self):
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        #filestring = self.docstring
        # VASP needs lattice vector matrix to have positive triple product
        if det3(lattice) < 0:
            if lattice[0].length() == lattice[1].length() == lattice[2].length():
                # Shift the first and last for cubic lattices
                transmtx = [[0, 0, 1],
                            [0, 1, 0],
                            [1, 0, 0]]
            else:
                # Else shift the two shortest
                if lattice[0].length() > lattice[1].length() and lattice[0].length() > lattice[2].length():
                    transmtx = [[1, 0, 0],
                                [0, 0, 1],
                                [0, 1, 0]]
                elif lattice[1].length() > lattice[2].length() and lattice[1].length() > lattice[0].length():
                    transmtx = [[0, 0, 1],
                                [0, 1, 0],
                                [1, 0, 0]]
                else:
                    transmtx = [[0, 1, 0],
                                [1, 0, 0],
                                [0, 0, 1]]
        else:
            transmtx = [[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
        lattice = mmmult3(transmtx,lattice)
        # Print braces around values (or not)
        if self.printbraces:
            lbrace = "{"
            rbrace = "}"
        else:
            lbrace = " "
            rbrace = ""
        #
        # length scale and lattice
        filestring  = self.docstring
        filestring += "\n"
        #
        # The atom position and pseudopotential info
        alloy = False
        spcs = ""
        natom = 0
        ntypat = 0
        ntypatpp_ve = 0.0
        alloystring = ""
        xredstring = ""
        ppstring = ""
        self.dir = ""
        for a in self.cell.atomdata:
            for b in a:
                natom += 1
                if spcs != b.spcstring():
                    ntypat += 1
                    atname = b.spcstring()
                    if b.alloy():
                        znuclstring += "?? "
                        alloystring += atname+"   "+atname+"5.0-s1p1"+"         "+atname+pseudostring+"\n"
                        alloy = True
                    else:
                        radial_cutoff = ed.openmxelements[b.spcstring()].split("-")[0].replace(b.spcstring(),"").replace(" ","")
                        potname = atname+radial_cutoff+".pao"
                        potfile = self.dir+"/"+potname
                        if self.dir == "" or not os.path.exists(self.dir) or exists_potname == "":
                            ppstring += ed.openmxelements[b.spcstring()]+"\n"
                            if natom == 1:
                                print "  Please, you get valence.electron from pseudo-potential file.\nSet 'up' and 'down' spin at '<Atoms.SpeciesAndCoordinatesn'.\n 'up'+'down' spin = valence.electron. \n  pseudo-potential file, "
                            print "    "+potname
                v = Vector(mvmult3(transmtx,b.position))
                half_ve = ntypatpp_ve/2
                xredstring += "  %-3i   %-2s    %9.7f   %9.7f   %9.7f   %4.1f %4.1f\n"%(natom,str(b.spcstring()),v[0],v[1],v[2],half_ve,half_ve)
                spcs = b.spcstring()
        #
        filestring += "# Definition of Atomic Species\n"
        filestring += "Species.Number              %7i\n"%int(ntypat)
        filestring += "<Definition.of.Atomic.Species\n"
        if alloy:
            filestring += "    # "+alloystring
        filestring += ppstring
        filestring = filestring[:-1]+rbrace+" \n"
        filestring += "Definition.of.Atomic.Species>\n"
        filestring += "\n"
        #
        filestring += "# Atoms\n"
        filestring += "Atoms.Number                %7i\n"%int(natom)
        filestring += "Atoms.SpeciesAndCoordinates.Unit   FRAC # Ang|AU|FRAC\n"
        filestring += "<Atoms.SpeciesAndCoordinates\n"
        #
        if alloy:
            filestring += "    # "+alloystring
        filestring += xredstring
        filestring = filestring[:-1]+rbrace+" \n"
        #
        filestring += "Atoms.SpeciesAndCoordinates>\n"
        filestring += "\n"
        #
        filestring += "# Unit Vectors\n"
        filestring += "Atoms.UnitVectors.Unit             Ang # Ang|AU\n"
        filestring += "<Atoms.UnitVectors\n"
        #
        a = self.cell.lengthscale
        alattice = lattice
        for i in range(0,3):
            for j in range(0,3):
                alattice[i][j] = lattice[i][j]*a
        for vec in alattice:
            # self.cell.lengthscale
            filestring +="  %9.7f   %9.7f   %9.7f\n"%(Vector(vec)[0],Vector(vec)[1],Vector(vec)[2])
        #
        filestring += "Atoms.UnitVectors>\n"
        filestring += "\n"
        #
        # k-space information
        #kresolution = 0.2
        reclatvect = self.cell.reciprocal_latticevectors()
        for j in range(3):
            for i in range(3):
                reclatvect[j][i] = reclatvect[j][i] / self.cell.lengthscale
        # Lengths of reciprocal lattice vectors
        reclatvectlen = [elem.length() for elem in reclatvect]
        self.kgrid = [max(1,int(round(elem/self.kresolution))) for elem in reclatvectlen]
        times = 1
        #
        filestring += "scf.Kgrid                  "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"       # means n1 x n2 x n3\n"
        #
        return filestring


class OpenMX_input_File(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an OpenMX run and the method
    __str__ that outputs the contents of a OpenMX input file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell = crystalstructure
        #self.cell.newunit("bohr")
        self.cell.newunit("angstrom")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        self.printbraces = False
    def __str__(self):
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        #filestring = self.docstring
        # VASP needs lattice vector matrix to have positive triple product
        if det3(lattice) < 0:
            if lattice[0].length() == lattice[1].length() == lattice[2].length():
                # Shift the first and last for cubic lattices
                transmtx = [[0, 0, 1],
                            [0, 1, 0],
                            [1, 0, 0]]
            else:
                # Else shift the two shortest
                if lattice[0].length() > lattice[1].length() and lattice[0].length() > lattice[2].length():
                    transmtx = [[1, 0, 0],
                                [0, 0, 1],
                                [0, 1, 0]]
                elif lattice[1].length() > lattice[2].length() and lattice[1].length() > lattice[0].length():
                    transmtx = [[0, 0, 1],
                                [0, 1, 0],
                                [1, 0, 0]]
                else:
                    transmtx = [[0, 1, 0],
                                [1, 0, 0],
                                [0, 0, 1]]
        else:
            transmtx = [[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
        lattice = mmmult3(transmtx,lattice)
        # Print braces around values (or not)
        if self.printbraces:
            lbrace = "{"
            rbrace = "}"
        else:
            lbrace = " "
            rbrace = ""
        #
        # length scale and lattice
        filestring  = "# File settings\n"
        filestring += "System.CurrrentDirectory         ./  # default=./\n"
        filestring += "System.Name                  "+self.filename.replace(".dat","")+"\n"
        filestring += "level.of.stdout                   1  # default=1 (1-3)\n"
        filestring += "level.of.fileout                  1  # default=1 (0-2)\n"
        filestring += "\n"
        #
        filestring = filestring[:-1]
        filestring += rbrace+" \n"
        #
        # POTCAR library
        directory = self.pseudolib
        if directory != "":
            self.dir = directory
        else:
            try:
                self.dir = os.environ['OpenMX_PSEUDOLIB']
            except:
                try:
                    self.dir = os.environ['OpenMX_POTLIB']
                except:
                    self.dir = "../..//DFT_DATA13/PAO"
        # The atom position and pseudopotential info
        alloy = False
        spcs = ""
        natom = 0
        ntypat = 0
        alloystring = ""
        xredstring = ""
        ntypatpp_ve = 0.0
        ppstring = ""
        for a in self.cell.atomdata:
            for b in a:
                natom += 1
                if spcs != b.spcstring():
                    ntypat += 1
                    atname = b.spcstring()
                    if b.alloy():
                        znuclstring += "?? "
                        alloystring += atname+"   "+atname+"5.0-s1p1"+"         "+atname+pseudostring+"\n"
                        alloy = True
                    else:
                        radial_cutoff = ed.openmxelements[b.spcstring()].split("-")[0].replace(b.spcstring(),"").replace(" ","")
                        potname = atname+radial_cutoff+".pao"
                        potfile = self.dir+"/"+potname
                        if os.path.exists(potfile):
                            exists_potname = potname
                            if self.dfttype == "GGA-PBE":
                                ppstring += ed.openmxelements[b.spcstring()]+"\n"
                            else:
                                ppstring += ed.openmxelements[b.spcstring()].replace("PBE","CA")+"\n"
                            print potfile
                            f = open(potfile,'r')
                            lines = f.readlines()
                            f.close()
                            for line in lines:
                                #if line[0:16] == "valence.electron": # bug fix not use search() 
                                #    ntypatpp_ve = float(line[17:].replace("\n","").replace(" ",""))
                                if search("valence\.electron",line):
                                    ntypatpp_ve = float(line.split("valence.electron")[1].lstrip(" ").split()[0].strip(string.punctuation))
                        if self.dir == "" or not os.path.exists(self.dir) or exists_potname == "":
                            ppstring += ed.openmxelements[b.spcstring()]+"\n"
                            if natom == 1:
                                print "  Please, you get valence.electron from pseudo-potential file.\nSet 'up' and 'down' spin at '<Atoms.SpeciesAndCoordinatesn'.\n'up'+'down' spin = valence.electron. \n  pseudo-potential file, "
                            print "    "+potname
                v = Vector(mvmult3(transmtx,b.position))
                half_ve = ntypatpp_ve/2
                xredstring += "  %-3i   %-2s    %9.7f   %9.7f   %9.7f   %4.1f %4.1f\n"%(natom,str(b.spcstring()),v[0],v[1],v[2],half_ve,half_ve)
                spcs = b.spcstring()
        #
        filestring += "# Definition of Atomic Species\n"
        filestring += "Species.Number              %7i\n"%int(ntypat)
        filestring += "<Definition.of.Atomic.Species\n"
        if alloy:
            filestring += "    # "+alloystring
        filestring += ppstring
        filestring = filestring[:-1]+rbrace+" \n"
        filestring += "Definition.of.Atomic.Species>\n"
        filestring += "\n"
        #
        if self.openmx_dft_u:
            filestring += "<Hubbard.U.values                 #  eV\n" 
            natoms = 0
            spcstring = ""
            for a in self.cell.atomdata:
                for b in a:
                    spcs = b.spcstring()
                    if spcs != spcstring:
                        natoms += 1
                        atomic_number = int(ed.elementnr[b.spcstring()])
                        if 21 <= atomic_number and atomic_number <= 29:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.81 + 0.08*(atomic_number-21)
                            U = F0
                            Ueff = (U - J)*float(self.openmx_dftu_times)
                            filestring += "  "+str(spcs)+" 1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d %5.3f\n"%(Ueff)
                        elif 39 <= atomic_number and atomic_number <= 47:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.59 + 0.056*(atomic_number-39)
                            U = F0
                            Ueff = (U - J)*float(self.openmx_dftu_times)
                            filestring += "  "+str(spcs)+" 1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d %5.3f\n"%(Ueff)
                        elif 71 <= atomic_number and atomic_number <= 79:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.860 + 0.053*(atomic_number-71)
                            U = F0
                            Ueff = (U - J)*float(self.openmx_dftu_times)
                            filestring += "  "+str(spcs)+" 1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d %5.3f\n"%(Ueff)
                        elif 57 <= atomic_number and atomic_number <= 70:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.90 + 0.036*(atomic_number-57)
                            U = F0
                            Ueff = (U - J)*float(self.openmx_dftu_times)
                            filestring += "  "+str(spcs)+" 1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d 0.0 2d 0.0 1f %5.3f\n"%(Ueff)
                        elif 89 <= atomic_number and atomic_number <= 103:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.66 + 0.035*(atomic_number-57)
                            U = F0
                            Ueff = (U - J)*float(self.openmx_dftu_times)
                            filestring += "  "+str(spcs)+" 1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d 0.0 2d 0.0 1f %5.3f\n"%(Ueff)
                        else:
                            l = 0
                            U = 0.0 
                            J = 0.0 
                            Ueff = 0.0
                            filestring += "  "+str(spcs)+" 1s 0.0 2s 0.0 1p 0.0 2p 0.0\n"
                    spcstring = spcs 
            filestring += "Hubbard.U.values>>\n"
            filestring += "\n"
        #
        filestring += "# Atoms\n"
        filestring += "Atoms.Number                %7i\n"%int(natom)
        filestring += "Atoms.SpeciesAndCoordinates.Unit   FRAC # Ang|AU|FRAC\n"
        filestring += "<Atoms.SpeciesAndCoordinates\n"
        #
        if alloy:
            filestring += "    # "+alloystring
        filestring += xredstring
        filestring = filestring[:-1]+rbrace+" \n"
        #
        filestring += "Atoms.SpeciesAndCoordinates>\n"
        filestring += "\n"
        #
        filestring += "# Unit Vectors\n"
        filestring += "Atoms.UnitVectors.Unit             Ang # Ang|AU\n"
        filestring += "<Atoms.UnitVectors\n"
        #
        a = self.cell.lengthscale
        alattice = lattice
        for i in range(0,3):
            for j in range(0,3):
                alattice[i][j] = lattice[i][j]*a
        for vec in alattice:
            # self.cell.lengthscale
            filestring +="  %9.7f   %9.7f   %9.7f\n"%(Vector(vec)[0],Vector(vec)[1],Vector(vec)[2])
        #
        filestring += "Atoms.UnitVectors>\n"
        filestring += "\n"
        #
        tasks = self.run_type
        #
        filestring += "# SCF or Electronic System\n"
        if tasks == "nonscf+dos" or tasks == "nonscf+band":
            rescf = "on "
        else:
            rescf = "off"
        filestring += "scf.restart                "+rescf+"         # on|off,default=off\n"
        filestring += "scf.XcType                 "+self.dfttype+"     # LDA|LSDA-CA|LSDA-PW|GGA-PBE\n"
        if self.openmx_dft_u:
            filestring += "scf.Hubbard.U              off          # On|Off , default=off\n"
            filestring += "scf.Hubbard.Occupation     dual        # onsite|full|dual, default=dual\n"
        filestring += "scf.SpinPolarization       On          # On|Off|NC\n"
        for a in self.cell.atomdata:
            for b in a:
                atomic_number = int(ed.elementnr[b.spcstring()])
                if atomic_number >= 50:
                    soc = "On "
                else:
                    soc = "Off"
        filestring += "scf.SpinOrbit.Coupling     "+soc+"         # On|Off, default=off \n"
        filestring += "scf.ElectronicTemperature  300.0       # default=300 (K)\n"
        filestring += "scf.energycutoff           220.0       # default=150 (Ry)\n"
        #filestring += "scf.Ngrid                 32 32 32  \n" 
        filestring += "scf.maxIter                100         # default=40\n"
        filestring += "scf.EigenvalueSolver       band        # DC|GDC|Cluster|Band\n"
        #
        # k-space information
        #kresolution = 0.2
        reclatvect = self.cell.reciprocal_latticevectors()
        for j in range(3):
            for i in range(3):
                reclatvect[j][i] = reclatvect[j][i] / self.cell.lengthscale
        # Lengths of reciprocal lattice vectors
        reclatvectlen = [elem.length() for elem in reclatvect]
        self.kgrid = [max(1,int(round(elem/self.kresolution))) for elem in reclatvectlen]
        times = 1
        #
        filestring += "scf.Kgrid                  "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"       # means n1 x n2 x n3\n"
        filestring += "scf.Mixing.Type            rmm-diisk   # Simple|Rmm-Diis|Gr-Pulay|Kerker|Rmm-Diisk\n"
        filestring += "scf.Init.Mixing.Weight     0.05        # default=0.30 \n"
        filestring += "scf.Min.Mixing.Weight      0.01        # default=0.001 \n"
        filestring += "scf.Max.Mixing.Weight      0.30        # default=0.40 \n"
        filestring += "scf.Mixing.History         25          # default=5\n"
        filestring += "scf.Mixing.StartPulay      15          # default=6\n"
        filestring += "scf.criterion              %-4.1e     # default=1.0e-6 (Hartree) | cif2cell set ca. 1 meV/atom\n"%(3.675e-5*natom)
        filestring += "scf.lapack.dste            dstevx      # dstegr|dstedc|dstevx, default=dstevx\n"
        #
        if tasks == "stress":
            filestring += "scf.Kerker.factor         2.0\n"
            filestring += "scf.ProExpn.VNA              on        # default=on\n"
            filestring += "scf.stress.tensor            on        # default=off\n"
        #
        filestring += "\n"
        #
        if tasks == "1DFFT":
            filestring += "# 1D FFT\n"
            filestring += "1DFFT.NumGridK             900         # default=900\n"
            filestring += "1DFFT.NumGridR             900         # default=900\n"
            filestring += "1DFFT.EnergyCutoff        3600.0       # default=3DFFT.EnergyCutoff*3.0 (Ry)\n"
            filestring += "\n"
        #
        if tasks == "co":
            filestring += "# output of contracted orbitals\n"
            filestring += "CntOrb.fileout               off       # on|off, default=off\n"
            filestring += "Num.CntOrb.Atoms              1        # default=1\n"
            filestring += "<Atoms.Cont.Orbitals\n"
            filestring += " 1\n"
            filestring += "Atoms.Cont.Orbitals>\n"
            filestring += "\n"
        #
        if tasks == "Order-N":
            filestring += "# SCF Order-N\n"
            filestring += "orderN.HoppingRanges        6.8        # default=5.0 (Ang)\n" 
            filestring += "orderN.NumHoppings           2         # default=2\n"
            filestring += "orderN.KrylovH.order        350        # default=400\n"
            filestring += "#orderN.KrylovS.order       2048       # default=4*orderN.KrylovH.order\n"
            filestring += "orderN.Expand.Core           on\n"
            filestring += "orderN.Recalc.Buffer         on\n"
            filestring += "orderN.Exact.Inverse.S       on\n"
            filestring += "\n"
        #
        if tasks == "MD" or tasks == "opt":
            filestring += "# MD or Geometry Optimization\n"
            filestring += "MD.Type                    nomd        # NVE|NVT_VS|NVT_NH\n"
            filestring += "                                       # Nomd|Opt|DIIS|NVE|NVT_VS|NVT_NH\n"
            filestring += "MD.Opt.DIIS.History        4           # default=4\n"
            filestring += "MD.Opt.StartDIIS           5           # default=5\n"
            filestring += "MD.Opt.EveryDIIS           10000       # default=10\n"
            filestring += "MD.maxIter                 1           # default=1\n"
            filestring += "MD.TimeStep                1.0         # default=0.5 (fs)\n"
            filestring += "MD.Opt.criterion           0.0003      # default=1.0e-4 (Hartree/bohr)\n"
            filestring += "MD.Opt.DIIS.Mixing         0.1         # default=0.1\n"
            filestring += "MD.Initial.MaxStep         0.001       # defalut=0.02 (Ang) \n"
            filestring += "\n"
        #
        if tasks == "MO":
            filestring += "# MO output\n"
            filestring += "MO.fileout                 on          # on|off, default=off\n"
            filestring += "num.HOMOs                  1           # default=1\n"
            filestring += "num.LUMOs                  1           # default=1\n"
            filestring += "\n"
        #
        if tasks == "dos" or tasks == "nonscf+dos":
            filestring += "# DOS and PDOS\n"
            filestring += "Dos.fileout                on          # on|off, default=off\n"
            filestring += "Dos.Erange                -20.0  20.0  # default = -20 20\n"
            times = 2
            filestring += "Dos.Kgrid                  "+str(self.kgrid[0]*times)+" "+str(self.kgrid[1]*times)+" "+str(self.kgrid[2]*times)+"       # default = Kgrid1 Kgrid2 Kgrid3\n"
            filestring += "\n"
        if tasks == "band" or tasks == "nonscf+band":
            #filestring += "<Band.KPath.UnitCell\n"
            #filestring += " 5.76  0.00  0.00\n"
            #filestring += " 0.00  5.76  0.00\n"
            #filestring += " 0.00  0.00  5.76\n"
            #filestring += "Band.KPath.UnitCell>\n"
            #filestring += "# if <Band.KPath.UnitCell does not exist,\n"
            #filestring += "#     the reciprical lattice vector is employed.\n"
            #
            if self.cell.crystal_system() == "hexagonal":
                # hcp, primitive cell
                filestring += "# Band dispersion (hcp type)\n"
                filestring += "Band.dispersion              on        # on|off, default=off\n"
                filestring += "Band.Nkpath                4\n"
                filestring += "<Band.kpath                \n"
                filestring += "   20  0.0   0.0   0.0     0.5   0.0   0.0     g M\n"
                filestring += "   40  0.5   0.0   0.0     0.333 0.333 0.0     M K\n"
                filestring += "   40  0.333 0.333 0.0     0.0   0.0   0.0     K g\n"
                filestring += "   20  0.0   0.0   0.0     0.0   0.0   0.5     g A\n"
            elif self.cell.spacegroupsetting == "F":
                # fcc, primitive cell  
                filestring += "# Band dispersion (fcc type)\n"
                filestring += "Band.dispersion              on        # on|off, default=off\n"
                filestring += "Band.Nkpath                5\n"
                filestring += "<Band.kpath                \n"              
                filestring += "   20  0.5   0.25  0.75    0.5   0.0   0.0     W L\n"
                filestring += "   40  0.5   0.0   0.0     0.0   0.0   0.0     L g\n"
                filestring += "   40  0.0   0.0   0.0     0.5   0.5   0.0     g X\n"
                filestring += "   20  0.5   0.5   0.0     0.75  0.5   0.25    X W\n"
                filestring += "   10  0.75  0.5   0.25    0.75  0.375 0.375   W K\n"
            elif self.cell.spacegroupsetting == "I":
                # bcc, primitive cell
                filestring += "# Band dispersion (bcc type)\n"
                filestring += "Band.dispersion              on        # on|off, default=off\n"
                filestring += "Band.Nkpath                4\n"
                filestring += "<Band.kpath                \n"
                filestring += "   20  0.0   0.0   0.0     0.5  -0.5   0.5     g H\n"
                filestring += "   40  0.5  -0.5   0.5     0.0   0.0   0.5     H N\n"
                filestring += "   40  0.0   0.0   0.5     0.0   0.0   0.0     N g\n"
                filestring += "   20  0.0   0.0   0.0     0.25  0.25  0.25    g P\n"
            elif self.cell.spacegroupsetting == "P":
                # sc, primitive cell
                filestring += "# Band dispersion (simple cubic type)\n"
                filestring += "Band.dispersion              on        # on|off, default=off\n"
                filestring += "Band.Nkpath                4\n"
                filestring += "<Band.kpath                \n"
                filestring += "   20  0.5   0.5   0.5     0.0   0.0   0.0     R g\n"
                filestring += "   40  0.0   0.0   0.0     0.5   0.0   0.0     g X\n"
                filestring += "   40  0.5   0.0   0.0     0.5   0.5   0.0     X M\n"
                filestring += "   20  0.5   0.5   0.0     0.0   0.0   0.0     M g\n"
            else:
                # sc, primitive cell
                filestring += "# Band dispersion (simple cubic type)\n"
                filestring += "Band.dispersion              on        # on|off, default=off\n"
                filestring += "Band.Nkpath                4\n"
                filestring += "<Band.kpath                \n"
                filestring += "   20  0.5   0.5   0.5     0.0   0.0   0.0     R g\n"
                filestring += "   40  0.0   0.0   0.0     0.5   0.0   0.0     g X\n"
                filestring += "   40  0.5   0.0   0.0     0.5   0.5   0.0     X M\n"
                filestring += "   20  0.5   0.5   0.0     0.0   0.0   0.0     M g\n"
            #
            filestring += "Band.kpath>\n"
            filestring += "\n"
        #
        if tasks == "Hoverlap":
            filestring += "# output Hamiltonian and overlap\n"
            filestring += "HS.fileout                   off       # on|off, default=off\n"
            filestring += "\n" 
        #
        return filestring

################################################################################################
# XYZ FILE
class LAMMPSFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting an .data LAMMPS file
    and the method __str__ that outputs the contents of the .data file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # To be put on the second line
        self.programdoc = ""
    def __str__(self):
        filestring = ""
        filestring += "#"+self.docstring+"\n\n"
        filestring += "%i atoms\n" % sum([len(v) for v in self.cell.atomdata])
        atomTypes = {}
        nextAtomTypeId = 1

        for a in self.cell.atomdata:
            for b in a:
                atomType = str(b).split()[0]
                if not atomType in atomTypes:
                    atomTypes[atomType] = nextAtomTypeId
                    nextAtomTypeId += 1
        filestring += "%i atom types\n\n" % len(atomTypes)

        if self.cell.latticevectors[0][1]!=0:
            theta = math.atan2(-self.cell.latticevectors[0][1], self.cell.latticevectors[0][0])
            c = cos(theta)
            s = sin(theta)
            R = LatticeMatrix([[c, s, 0],
                                [-s, c, 0],
                                [0, 0, 1]])
            self.cell.latticevectors[0] = Vector(mvmult3(R,self.cell.latticevectors[0]))
            self.cell.latticevectors[1] = Vector(mvmult3(R,self.cell.latticevectors[1]))
            self.cell.latticevectors[2] = Vector(mvmult3(R,self.cell.latticevectors[2]))

        if self.cell.latticevectors[0][2]!=0:
            theta = math.atan2(-self.cell.latticevectors[0][2], self.cell.latticevectors[0][0])
            c = cos(theta)
            s = sin(theta)
            R = LatticeMatrix([[c, s, 0],
                                [0, 1, 0],
                                [-s, c, 0]])
            self.cell.latticevectors[0] = Vector(mvmult3(R,self.cell.latticevectors[0]))
            self.cell.latticevectors[1] = Vector(mvmult3(R,self.cell.latticevectors[1]))
            self.cell.latticevectors[2] = Vector(mvmult3(R,self.cell.latticevectors[2]))

        if self.cell.latticevectors[1][2]!=0:
            theta = math.atan2(-self.cell.latticevectors[1][2], self.cell.latticevectors[1][1])
            c = cos(theta)
            s = sin(theta)
            R = LatticeMatrix([[1, 0, 0],
                                [0, c, s],
                                [0, -s, c]])
            self.cell.latticevectors[0] = Vector(mvmult3(R,self.cell.latticevectors[0]))
            self.cell.latticevectors[1] = Vector(mvmult3(R,self.cell.latticevectors[1]))
            self.cell.latticevectors[2] = Vector(mvmult3(R,self.cell.latticevectors[2]))

        if self.cell.latticevectors[0][1]!=0 or self.cell.latticevectors[0][2] != 0 or self.cell.latticevectors[1][2]!=0 or self.cell.latticevectors[0][0] <= 0 or self.cell.latticevectors[1][1] <= 0 or self.cell.latticevectors[2][2] <= 0:
            print "Error in triclinic box. Vectors should follow these rules: http://lammps.sandia.gov/doc/Section_howto.html#howto-12"
            print "Ideally, this program should solve this, but it doesn't yet. You need to fix it."
            exit()

        xy = self.cell.lengthscale*self.cell.latticevectors[1][0]
        xz = self.cell.lengthscale*self.cell.latticevectors[2][0]
        yz = self.cell.lengthscale*self.cell.latticevectors[2][1]

        a = self.cell.latticevectors[0][0]*self.cell.lengthscale
        b = self.cell.latticevectors[1][1]*self.cell.lengthscale
        c = self.cell.latticevectors[2][2]*self.cell.lengthscale

        filestring += "0.0 %f xlo xhi\n" % a
        filestring += "0.0 %f ylo yhi\n" % b
        filestring += "0.0 %f zlo zhi\n" % c
        if xy!=0 or xz !=0 or yz != 0:
            filestring += str(xy) + " " + str(xz) + " " + str(yz) + " xy xz yz\n"
        
        filestring += "\n"
        filestring += "Masses\n\n"
        
        matomTypes = {}
        mnextAtomTypeId = 1
        for a in self.cell.atomdata:
            for b in a:
                sp_b = b.spcstring()
                matomType = str(b).split()[0]
                if not matomType in matomTypes:
                    filestring += str(mnextAtomTypeId)+" "+str(ed.elementweight[sp_b])+" # "+str(sp_b)+"\n"
                    matomTypes[matomType] = mnextAtomTypeId
                    mnextAtomTypeId += 1
        
        filestring += "\n"
        filestring += "Atoms\n\n"

        nextAtomId = 1

        #for b in [a for a in self.cell.atomdata]:
            #print str(b).split()[0]
        #atomTypes str(b).split()[0]

        lv = []
        for i in range(3):
            lv.append([])
            for j in range(3):
                lv[i].append(self.cell.lengthscale*self.cell.latticevectors[i][j])
        for a in self.cell.atomdata:
            for b in a:
                t = Vector(mvmult3(lv,b.position))
                atomType = str(b).split()[0]
                atomTypeId = atomTypes[atomType]
                if self.lammps_type == "charge":
                    filestring += str(nextAtomId)+" "+str(atomTypeId)+" 0.0 "+str(t)+"\n"
                else:
                    filestring += str(nextAtomId)+" "+str(atomTypeId)+" "+str(t)+"\n"
                nextAtomId += 1
        return filestring

################################################################################################
class FdmnesFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an *_inp.txt file and the method
    __str__ that outputs the contents of an *_inp.txt file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("angstrom")
        # Make sure the docstring has the form of a f90 comment
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("!")
            string = "!"+string+"\n"
            self.docstring += string
    def __str__(self):
        filestring = self.docstring
        filestring += "\n"
        # Lattice vectors
        filestring += " Crystal  ! Periodic material description (unit cell)\n"
        #filestring += "lattice parameters  a b c  [Angstrom]\n"
        a = self.cell.a
        b = self.cell.b
        c = self.cell.c
        filestring += " %7.3f %7.3f %7.3f  %5.2f  %5.2f  %5.2f  ! a, b, c, (Angstroem) alpha, beta, gamma (degree)\n"%(a,b,c,self.cell.alpha,self.cell.beta,self.cell.gamma)
        natom = 0
        for a in self.cell.atomdata:
            for b in a:
                if natom == 0:
                    filestring += " "+str(ed.elementnr[b.spcstring()])+" %8.5f"%float(b.position[0])+" %8.5f"%float(b.position[1])+" %8.5f"%float(b.position[2])+"  ! Z, x, y, z (unit cell unit)\n"
                else:
                    filestring += " "+str(ed.elementnr[b.spcstring()])+" %8.5f"%float(b.position[0])+" %8.5f"%float(b.position[1])+" %8.5f"%float(b.position[2])+"\n"
                natom += 1
        return filestring


class Fdmnes_input_File(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an fdmnes.in file and the method
    __str__ that outputs the contents of an fdmnes.in file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("angstrom")
        self.filename = ""
        # Make sure the docstring has the form of a f90 comment
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("!")
            string = "!"+string+"\n"
            self.docstring += string
    def __str__(self):
        #filestring = self.docstring
        filestring =  "! Fdmnes indata file\n"
        if  self.run_type == "tddft" and self.fdmnes_edge == "K":
            filestring += "! Calculation for the copper L23-edge in copper cfc\n"
        elif self.run_type == "circular" and self.fdmnes_edge == "K":
            filestring += "! Calculation for the copper L23-edge in copper cfc\n"
        else:
            filestring += "! Calculation for the copper "+str(self.fdmnes_edge)+"-edge in copper cfc\n"
        filestring += "! Finite difference method calculation with convolution\n"
        filestring += "\n"
        filestring += " Filout\n"
        filestring += "  "+self.filename+"\n"
        filestring += "\n"
        filestring += " Range   ! Energy range of calculation (eV)\n"
        if self.run_type == "xes":
            filestring += "  -60. 0.2  0. 0.5 10. 1. 80.  ! first energy, step, intermediary energy, step ..., last energy \n"
        elif self.run_type == "tddft":
            filestring += "  -1.  0.05 8. 0.5 30  ! first energy, step, intermediary energy, step ..., last energy \n"
        else:
            filestring += "  -15. 0.2  5. 0.5 20. 1. 50.  ! first energy, step, intermediary energy, step ..., last energy \n"
        filestring += "\n"
        filestring += " Radius  ! Radius of the cluster where final state calculation is performed\n"
        filestring += "  6.0    ! For a good calculation, this radius must be increased up to 6 or 7 Angstroems\n"
        filestring += "\n"
        # run type
        if self.fdmnes_dos == "DOS" or self.run_type.upper() == "DOS":
            filestring += " Density ! Density of state\n"
            filestring += "\n"
        elif self.fdmnes_dos_comp == "DOS_comp" or self.run_type.upper() == "DOS_COMP":
            filestring += " Density_comp ! output with the usual complex spherical harmonics\n"
            filestring += "\n"
        elif self.fdmnes_dos_all == "DOS_all" or self.run_type.upper() == "DOS_ALL":
            filestring += " Density_all ! DOS for every atom\n"
            filestring += "\n"
        if self.fdmnes_green == "Green" or self.run_type.lower() == "green":
            filestring += " Green   ! muffin-tin potential\n"
            filestring += "\n"
            if self.fdmnes_eimag != 0:
                filestring += " Eimag   ! Green setting case only.  It is good for Photoemission (XES)\n"
                filestring += "  "+str(self.fdmnes_eimag)+"    ! value of the uniform width (eV)\n"
                filestring += "\n"
        if self.fdmnes_spherical:
            filestring += " Spherical ! the spherical tensors (in number of electron)\n"
            filestring += "\n"
        if self.fdmnes_photon_energy == "photon_energy":
            filestring += " Energpho  ! show photon energy for x axis\n"
            filestring += "\n"
        #
        if self.dft_type.upper()  == "LDA":
            pass
        elif self.dft_type.upper()  == "GGA":
            filestring += " Perdew  ! GGA\n"
            filestring += "\n"
        elif self.dft_type.upper()  == "XA":
            filestring += " Xalpha\n"
            filestring += "  0.7\n"
            filestring += "\n"
        #
        if self.fdmnes_dft_u:
            filestring += "\n"
            filestring += " Hubbard  ! DFT+U\n"
            ldau_filestring = ""
            natoms = 0
            spcstring = ""
            for a in self.cell.atomdata:
                for b in a:
                    spcs = b.spcstring()
                    if spcs != spcstring:
                        natoms += 1
                        atomic_number = int(ed.elementnr[b.spcstring()])
                        if 21 <= atomic_number and atomic_number <= 29:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.81 + 0.08*(atomic_number-21)
                            U = F0
                            Ueff = U - J
                        elif 39 <= atomic_number and atomic_number <= 47:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.59 + 0.056*(atomic_number-39)
                            U = F0
                            Ueff = U - J
                        elif 71 <= atomic_number and atomic_number <= 79:
                            l = 2
                            F0 = 15.31 + 1.50*(atomic_number-21)
                            J = 0.860 + 0.053*(atomic_number-71)
                            U = F0
                            Ueff = U - J
                        elif 57 <= atomic_number and atomic_number <= 70:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.90 + 0.036*(atomic_number-57)
                            U = F0
                            Ueff = U - J
                        elif 89 <= atomic_number and atomic_number <= 103:
                            l = 3
                            F0 = 2.38 + 0.93*(atomic_number-57)
                            J = 0.66 + 0.035*(atomic_number-57)
                            U = F0
                            Ueff = U - J
                        else:
                            l = 0
                            U = 0.0
                            J = 0.0
                            Ueff = 0.0
                        #
                        Ueff = Ueff*float(self.fdmnes_dftu_times)
                        ldau_filestring += " %5.3f "%(Ueff)
                    spcstring = spcs 
            filestring += "  "+ldau_filestring+"\n"
        #
        if self.fdmnes_core_hole:
            filestring += " Screening  ! value of Screening < 1.0\n"
            filestring += "  "+str(self.fdmnes_core_hole)+"\n"
            filestring += "\n"
        #
        if self.fdmnes_edge == "K" and self.run_type == "tddft":
            filestring += " Edge   ! TDDFT calculation\n"
            filestring += "  L23\n"
            filestring += "\n"
        elif self.fdmnes_edge == "K" and self.run_type == "circular":
            filestring += " Edge   ! TDDFT calculation\n"
            filestring += "  L23\n"
            filestring += "\n"
        elif self.fdmnes_edge != "L23" and self.run_type == "tddft":
            pass
        elif self.fdmnes_edge != "L23" and self.run_type == "circular":
            pass
        else:
            filestring += " Edge\n"
            filestring += "  "+str(self.fdmnes_edge)+"\n"
            filestring += "\n"
        #
        if self.run_type == "tddft":
            filestring += " TDDFT  ! TDDFT calculation\n"
            filestring += "\n"
            spcstring = ""
            tddft_lmax = 2
            for a in self.cell.atomdata:
                for b in a:
                    spcs = b.spcstring()
                    if spcs != spcstring:
                        atomic_number = int(ed.elementnr[b.spcstring()])
                        if 57 <= atomic_number and atomic_number <= 70:
                            tddft_lmax = 3
                        elif 89 <= atomic_number and atomic_number <= 103:
                            tddft_lmax = 3
            filestring += " lmax   ! Limitation of the spherical harmonics to save some times\n"
            filestring += "  "+str(tddft_lmax)+"\n"
            filestring += "\n"
        #
        if self.fdmnes_scf == "SCF" or self.run_type.upper() == "SCF" or self.fdmnes_spin == "spin" or self.run_type == "spin" or self.run_type.upper() == "TDDFT" or self.run_type == "optic" or self.run_type == "dafs" or self.run_type == "nrixs":
            filestring += " SCF     ! SCF calculation\n"
            filestring += "\n"
        #
        if self.run_type == "xafs" or self.run_type == "fdmx":
            filestring += " FDMX  ! XAFS calculation (K-edge only for ver.2019/03/14)\n"
            filestring += "\n"
        #
        if self.fdmnes_spin == "spin" or self.run_type == "spin":
            filestring += " Magnetism  ! Spin Polarized calculation\n"
            filestring += "\n"
            filestring += " Polarisation\n"
            filestring += "  0. 0. 1.\n"
            filestring += "\n"
        #
        if self.run_type == "optic":
            filestring += " Optic\n"
            filestring += "\n"
        #
        if self.fdmnes_multipolar_type != "E1E1":
            filestring += " "+self.fdmnes_multipolar_type+"\n"
            filestring += "\n"
        #
        if self.run_type == "optic":
            if self.fdmnes_multipolar_type != "E1E2":
                filestring += " E1E2\n"
            if self.fdmnes_multipolar_type != "E2E2":
                filestring += " E2E2\n"
            if self.fdmnes_multipolar_type != "E1M1":
                filestring += " E1M1\n"
            if self.fdmnes_multipolar_type != "M1M1":
                filestring += " M1M1\n"     
            filestring += "\n"
            filestring += " Polarize\n"
            filestring += "  0.0 0.0 0.0   0.0 0.0 1.0\n"     
            filestring += "  0.0 0.0 0.0   0.0 1.0 0.0\n"
            filestring += "  0.0 0.0 0.0   1.0 0.0 0.0\n"
            filestring += "\n"
        #
        if self.run_type == "circular":
            filestring += " Full_Self_abs ! to get the self-absorption including the birefringence effect\n"
            filestring += "\n"
        #
        if self.run_type == "dichroism":
            filestring += " Polarize  ! circular dichroism\n"
            filestring += "  0.0 0.0 0.0   0.0 0.0 1.0\n"     
            filestring += "\n"
        #
        if self.run_type == "nrixs":
            filestring += " NRIXS\n"
            filestring += " 3. 6. 9.  ! the different values of the modulus of q in A^-1\n"
            filestring += "\n"
            filestring += " All_nrixs\n"
            filestring += "\n"
        #
        if self.run_type == "dafs":
            filestring += " DAFS\n"
            filestring += "  1 1 0  1 1 90.  ! h k l  Sigma Sigma Azimuth(degree)\n"
            filestring += "  ! h k l of reflection indices\n"
            filestring += "  ! Sigma=1, Pi=2\n"
            filestring += "  ! If the azimuth is not specified, a scan is performed\n"
            filestring += "\n"
        #
        if self.run_type == "rxs":
            if self.fdmnes_multipolar_type != "Quadrupole":
                filestring += " Quadrupole\n"
            filestring += "\n"
            filestring += " Self_abs      ! To perfom self-absorption correction\n"
            filestring += "\n"
            filestring += " Spherical\n"
            filestring += "\n"
            filestring += " RXS\n"
            filestring += "  1 1 1   1 1  0. ! Sigma-Sigma; 0. is the azimuth.\n"
            filestring += "  1 1 1   1 2     ! Sigma-Pi; Azimuthal scan is automaticaly performed because the azimuth is not specified \n"
            filestring += "  2 2 2   1 1  0. ! (2,2,2) is a not forbidden reflexion. \n"
            filestring += "\n"
        #
        Spinorbit_flag = 0
        atomic_number = 0
        spinorb = ".false."
        for a in self.cell.atomdata:
            for b in a:
                atomic_number = int(ed.elementnr[b.spcstring()])
                if atomic_number > 36:
                    spinorb = ".true."
        if spinorb == ".true.":
            spcstring = ""
            scalar_rel = "Relativism"
            for a in self.cell.atomdata:
                for b in a:
                    spcs = b.spcstring()
                    if spcs != spcstring:
                        atomic_number = int(ed.elementnr[b.spcstring()])
                        if 57 <= atomic_number and atomic_number <= 70:
                            soc_rel = "Spinorbit"
                            Spinorbit_flag = 1
                        elif 89 <= atomic_number and atomic_number <= 103:
                            soc_rel = "Spinorbit"
                            Spinorbit_flag = 1
            filestring += " "+scalar_rel+"  ! atomic number > 36\n"
            if  Spinorbit_flag == 1:
                filestring += " "+soc_rel+"  ! 4f electron near Fermi edge\n"
            filestring += "\n"
        #
        filestring += " Absorber\n"
        filestring += "  "+str(self.fdmnes_abs_no)+"\n"
        filestring += "\n"
        # Lattice vectors
        if self.fdmnes_molecular == "Molecule":
            filestring += " "+self.fdmnes_molecular+"  ! Periodic material description (unit cell)\n"
        else:
            filestring += " Crystal  ! Periodic material description (unit cell)\n"
        #filestring += "lattice parameters  a b c  [Angstrom]\n"
        a = self.cell.a
        b = self.cell.b
        c = self.cell.c
        filestring += "  "+"%7.4f %7.4f %7.4f  %3.2f  %3.2f  %3.2f  ! a, b, c, (Angstroem) alpha, beta, gamma (degree)\n"%(a,b,c,self.cell.alpha,self.cell.beta,self.cell.gamma)
        abs_sign = ""
        natom = 0
        for a in self.cell.atomdata:
            for b in a:
                if natom == (int(self.fdmnes_abs_no) - 1):
                    abs_sign = ": absorption atom"
                else:
                    abs_sign = ""
                #
                if natom == 0:
                    filestring += " "+str(ed.elementnr[b.spcstring()])+" %8.5f"%float(b.position[0])+" %8.5f"%float(b.position[1])+" %8.5f"%float(b.position[2])+"  ! %5i"%int(natom+1)+": Z, x, y, z (unit cell unit)"+"  "+abs_sign+"\n"
                else:
                    filestring += " "+str(ed.elementnr[b.spcstring()])+" %8.5f"%float(b.position[0])+" %8.5f"%float(b.position[1])+" %8.5f"%float(b.position[2])+"  "+"! %5i"%int(natom+1)+abs_sign+"\n"
                natom += 1
        #
        filestring += " \n"
        filestring += " Convolution\n"
        filestring += " \n"
        #
        if self.run_type == "tddft":
            filestring += " Estart\n"
            filestring += " -5.\n"
            filestring += "\n"
        #
        if self.run_type == "circular":
            filestring += " Circular    ! to get the circular - sigma and circular - left polarizations\n"
            filestring += " Double_cor  ! to get the double corrected intensity\n"
            filestring += " Estart\n"
            filestring += " -20.\n"
            filestring += "\n"
            if Spinorbit_flag != 1:
                filestring += " Spinorbite\n"
            if self.fdmnes_multipolar_type != "Dipmag":
                filestring += " Dipmag       ! to get the E1-M1 term and the M1-M1\n"
            filestring += "\n"
            filestring += " Axe_spin\n"
            filestring += " 0. 1. 0.\n"
            filestring += "\n"
        #
        if self.run_type == "optic":
            filestring += " Gamma_max\n"
            filestring += "  0.\n"
            filestring += "\n"
        #
        if self.run_type == "rxs":
            filestring += " Estart\n"
            filestring += " -20.\n"
            filestring += "\n"
        #
        if self.run_type == "fit":
            filestring += " Estart\n"
            filestring += " -10.\n"
            filestring += "\n"
            filestring += "! keywords for the comparison with the experiment and fit \n"
            filestring += " Experiment\n"
            filestring += "  "+self.filename+"_exp.txt\n"
            filestring += "\n"
            filestring += " Gen_shift\n"
            filestring += "  -6. -2. 21\n"
            filestring += "\n"
            filestring += " Parameter\n"
            filestring += "   Par_abc\n"
            filestring += "   -10. 0. 3\n"
            filestring += "\n"
            filestring += " Parameter\n"
            filestring += "   Par_Gamma_max\n"
            filestring += "   10. 15. 3\n"
            filestring += "\n"
        #
        filestring += " End\n"
        #
        f = open("fdmfile.txt", "w")
        fdmfile_filestring  = "! General indata file for FDMNES\n"
        fdmfile_filestring += "! with indata files examples\n"
        fdmfile_filestring += "\n"
        if self.run_type == "fit":
            fdmfile_filestring += " 1\n"
        else:
            fdmfile_filestring += " 2\n"
        fdmfile_filestring += self.filename+"_inp.txt\n"
        fdmfile_filestring += self.filename+"_conv_inp.txt\n"
        fdmfile_filestring += "\n"
        f.write(str(fdmfile_filestring))
        f.close()
        #
        f = open(self.filename+"_conv_inp.txt", "w")
        conv_filestring  = "! Fdmnes indata file\n"
        conv_filestring += "! Calculation of Chromium photoemission in Cr metal.\n"
        conv_filestring += "! The convolution is performed on a spectra already calculated in a previous step.\n"
        conv_filestring += "\n"
        conv_filestring += " Calculation\n"
        conv_filestring += "  "+self.filename+".txt\n"
        conv_filestring += "\n"
        conv_filestring += " Conv_out\n"
        conv_filestring += "  "+self.filename+"_photo_conv.txt\n"
        conv_filestring += "\n"
        conv_filestring += " Convolution\n"
        conv_filestring += "\n"
        if self.run_type == "xes":
            conv_filestring += " Photoemission  ! XES calculation\n"
            conv_filestring += "\n"
            conv_filestring += " Gaussian\n"
            conv_filestring += "  0.5  ! deltaE/E = about 5000-10000\n"
            conv_filestring += "\n"
        conv_filestring += " Gamma_hole\n"
        conv_filestring += "  0.7\n"
        conv_filestring += "\n"
        conv_filestring += "! To change the Fermi level (or energy of the first non occupied state)\n"
        conv_filestring += "! Efermi\n"
        conv_filestring += " !-2.62\n"
        conv_filestring += "\n"
        conv_filestring += "! To get the convoluted spectra starting at lower energy\n"
        conv_filestring += "! Estart\n"
        conv_filestring += " !-15.\n"
        conv_filestring += "\n"
        conv_filestring += "! To change the broadening width\n"
        conv_filestring += "! Gamma_max\n"      
        conv_filestring += " !7.\n"
        conv_filestring += "\n"
        conv_filestring += " End\n"
        f.write(str(conv_filestring))
        f.close()
        #
        f = open("xanes.plot", "w")
        plot_filestring  = "#!/usr/local/bin/gnuplot -persist\n"
        plot_filestring += "# Last modified: 2018/12/17\n"
        plot_filestring += "set terminal postscript eps enhanced 28 lw 2\n"
        plot_filestring += 'set output "'+self.filename+'.xanes.eps"\n'
        plot_filestring += "unset key\n"
        if self.fdmnes_photon_energy == "photon_energy":
            plot_filestring += "set xlabel 'Photon Energy (eV)'\n"     
        else:
            plot_filestring += "set xlabel 'Energy (eV)'\n"
        if self.run_type == "xes":
            plot_filestring += "set xtics 40\n"
        plot_filestring += "set mytics 5\n"
        plot_filestring += "set mxtics 5\n"
        plot_filestring += "set ylabel 'Absorption cross-section (Mbarn)'\n"
        if self.run_type == "xanes" or self.run_type == "dafs":
            plot_filestring += 'spec_filename = "'+self.filename+'_photo_conv.txt"\n'
        else:
            plot_filestring += 'spec_filename = "'+self.filename+'_conv.txt"\n'
        plot_filestring += "set yzeroaxis lt 2\n"
        plot_filestring += "plot spec_filename using 1:2 every 1:1:1 w l lt 1\n"
        plot_filestring += "\n"
        plot_filestring += '# bav_filename = "'+self.filename+'_bav.txt"\n'
        plot_filestring += '# ef = system("cat " . bav_filename . " | grep E_edge | head -1 |" . "'+" awk '{print $3}'"+'")\n'
        plot_filestring += "# set parametric\n"
        plot_filestring += "# set trange[0:1]\n"
        plot_filestring += "# const = abs(ef)\n"
        plot_filestring += "# plot const,t lt 2, spec_filename using 1:2 every 1:1:1 w l lt 1\n"
        f.write(str(plot_filestring))
        f.close()
        #
        if self.run_type == "tddft":
            f = open("tddft.plot", "w")
            tddft_plot_filestring  = "#!/usr/local/bin/gnuplot -persist\n"
            tddft_plot_filestring += "# Last modified: 2018/12/17\n"
            tddft_plot_filestring += "set terminal postscript eps enhanced 28 lw 2\n"
            tddft_plot_filestring += 'set output "'+self.filename+'.tddft.eps"\n'
            tddft_plot_filestring += "unset key\n"
            if self.fdmnes_photon_energy == "photon_energy":
                tddft_plot_filestring += "set xlabel 'Photon Energy (eV)'\n"       
            else:
                tddft_plot_filestring += "set xlabel 'Energy (eV)'\n"
            tddft_plot_filestring += "set mxtics 5\n"
            tddft_plot_filestring += "set mytics 5\n"
            tddft_plot_filestring += "set ylabel 'Absorption cross-section (Mbarn)'\n"
            tddft_plot_filestring += 'spec_filename = "'+self.filename+'_tddft.txt"\n'
            tddft_plot_filestring += "set yzeroaxis lt 2\n"
            tddft_plot_filestring += "plot spec_filename using 1:2 every 1:1:2 w l lt 1\n"
            f.write(str(tddft_plot_filestring))
            f.close()
        #
        if self.run_type == "xes":
            f = open("xes.plot", "w")
            xes_plot_filestring  = "#!/usr/local/bin/gnuplot -persist\n"
            xes_plot_filestring += "# Last modified: 2018/12/17\n"
            xes_plot_filestring += "set terminal postscript eps enhanced 28 lw 2\n"
            xes_plot_filestring += 'set output "'+self.filename+'.xes.eps"\n'
            xes_plot_filestring += "unset key\n"
            if self.fdmnes_photon_energy == "photon_energy":
                xes_plot_filestring += "set xlabel 'Photon Energy (eV)'\n"       
            else:
                xes_plot_filestring += "set xlabel 'Energy (eV)'\n"
            xes_plot_filestring += "set xtics 40\n"
            xes_plot_filestring += "set mxtics 5\n"
            xes_plot_filestring += "set mytics 5\n"
            xes_plot_filestring += "set ylabel 'Emission Intensity (arb.unit)'\n"
            xes_plot_filestring += 'spec_filename = "'+self.filename+'_photo_conv.txt"\n'
            xes_plot_filestring += "set yzeroaxis lt 2\n"
            xes_plot_filestring += "plot spec_filename using 1:2 every 1:1:1 w l lt 1\n"
            f.write(str(xes_plot_filestring))
            f.close()
        #
        if self.run_type == "circular":
            f = open("circular.plot", "w")
            circular_plot_filestring  = "#!/usr/local/bin/gnuplot -persist\n"
            circular_plot_filestring += "# Last modified: 2018/12/17\n"
            circular_plot_filestring += "set terminal postscript eps enhanced 28 lw 2\n"
            circular_plot_filestring += 'set output "'+self.filename+'.circular.eps"\n'
            circular_plot_filestring += "unset key\n"
            if self.fdmnes_photon_energy == "photon_energy":
                circular_plot_filestring += "set xlabel 'Photon Energy (eV)'\n"       
            else:
                circular_plot_filestring += "set xlabel 'Energy (eV)'\n"
            circular_plot_filestring += "set mxtics 5\n"
            circular_plot_filestring += "set mytics 5\n"
            circular_plot_filestring += "set ylabel 'Intensity (arb.unit)'\n"
            circular_plot_filestring += 'spec_filename = "'+self.filename+'_photo_conv.txt"\n'
            circular_plot_filestring += "set yzeroaxis lt 2\n"
            circular_plot_filestring += "plot spec_filename using 1:2 every 1:1:1 w l lt 1\n"
            f.write(str(circular_plot_filestring))
            f.close()
        #
        if self.run_type == "dichroism":
            f = open("dichroism.plot", "w")
            dichroism_plot_filestring  = "#!/usr/local/bin/gnuplot -persist\n"
            dichroism_plot_filestring += "# Last modified: 2018/12/17\n"
            dichroism_plot_filestring += "set terminal postscript eps enhanced 28 lw 2\n"
            dichroism_plot_filestring += 'set output "'+self.filename+'.dichroism.eps"\n'
            dichroism_plot_filestring += "unset key\n"
            if self.fdmnes_photon_energy == "photon_energy":
                dichroism_plot_filestring += "set xlabel 'Photon Energy (eV)'\n"       
            else:
                dichroism_plot_filestring += "set xlabel 'Energy (eV)'\n"
            dichroism_plot_filestring += "set mxtics 5\n"
            dichroism_plot_filestring += "set mytics 5\n"
            dichroism_plot_filestring += "set ylabel 'Absorption Intensity (arb.unit)'\n"
            dichroism_plot_filestring += "set y2label 'Circular Dichroism Intensity (arb.unit)'\n"
            dichroism_plot_filestring += 'spec_filename = "'+self.filename+'_photo_conv.txt"\n'
            dichroism_plot_filestring += "set yzeroaxis lt 2\n"
            dichroism_plot_filestring += "plot spec_filename using 1:2 every 1:1:1 axis x1y1 w l lt 1, spec_filename using 1:3 every 1:1:1 xis x1y2 w l lt 1\n"
            f.write(str(dichroism_plot_filestring))
            f.close()
        #
        if self.run_type == "optic":
            f = open("optic.plot", "w")
            optic_plot_filestring  = "#!/usr/local/bin/gnuplot -persist\n"
            optic_plot_filestring += "# Last modified: 2018/12/17\n"
            optic_plot_filestring += "set terminal postscript eps enhanced 28 lw 2\n"
            optic_plot_filestring += 'set output "'+self.filename+'.optic.eps"\n'
            optic_plot_filestring += "unset key\n"
            if self.fdmnes_photon_energy == "photon_energy":
                optic_plot_filestring += "set xlabel 'Photon Energy (eV)'\n"       
            else:
                optic_plot_filestring += "set xlabel 'Energy (eV)'\n"
            optic_plot_filestring += "set mytics 5\n"
            optic_plot_filestring += "set mytics 5\n"
            optic_plot_filestring += "set my2tics 5\n"
            optic_plot_filestring += "set ylabel 'Absorption Intensity (arb.unit)'\n"
            optic_plot_filestring += "set y2label 'Circular Dichroism Intensity (arb.unit)'\n"
            optic_plot_filestring += 'spec_filename = "'+self.filename+'_photo_conv.txt"\n'
            optic_plot_filestring += "set yzeroaxis lt 2\n"
            optic_plot_filestring += "plot spec_filename using 1:2 every 1:1:1 axis x1y1 w l lt 1 t 'Abs (001)', spec_filename using 1:3 every 1:1:1 axis x1y2 w l lt 1 t 'Dic(001)', spec_filename using 1:4 every 1:1:1 axis x1y1 w l lt 1 t 'Abs (010)', spec_filename using 1:5 every 1:1:1 axis x1y2 w l lt 1 t 'Dic(010)', spec_filename using 1:6 every 1:1:1 axis x1y1 w l lt 1 t 'Abs (100)', spec_filename using 1:7 every 1:1:1 axis x1y2 w l lt 1 t 'Dic(100)'\n"
            f.write(str(optic_plot_filestring))
            f.close()
        #
        if self.run_type == "dafs":
            f = open("dafs.plot", "w")
            dafs_plot_filestring  = "#!/usr/local/bin/gnuplot -persist\n"
            dafs_plot_filestring += "# Last modified: 2018/12/17\n"
            dafs_plot_filestring += "set terminal postscript eps enhanced 28 lw 2\n"
            dafs_plot_filestring += 'set output "'+self.filename+'.dafs.eps"\n'
            dafs_plot_filestring += "unset key\n"
            if self.fdmnes_photon_energy == "photon_energy":
                dafs_plot_filestring += "set xlabel 'Photon Energy (eV)'\n"       
            else:
                dafs_plot_filestring += "set xlabel 'Energy (eV)'\n"
            dafs_plot_filestring += "set mxtics 5\n"
            dafs_plot_filestring += "set mytics 5\n"
            dafs_plot_filestring += "set ylabel 'Intensity (arb.unit)'\n"
            dafs_plot_filestring += 'spec_filename = "'+self.filename+'_photo_conv.txt"\n'
            dafs_plot_filestring += "set yzeroaxis lt 2\n"
            dafs_plot_filestring += "plot spec_filename using 1:3 every 1:1:1 axis x1y1 w l lt 1\n"
            f.write(str(dafs_plot_filestring))
            f.close()
        #
        return filestring
    
