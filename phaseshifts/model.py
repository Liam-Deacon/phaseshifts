#!/usr/bin/env python
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Copyright: Copyright (C) 2013-2014 Liam Deacon                             #
#                                                                            #
# License: MIT License                                                       #
#                                                                            #
# Permission is hereby granted, free of charge, to any person obtaining a    #
# copy of this software and associated documentation files (the "Software"), #
# to deal in the Software without restriction, including without limitation  #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,   #
# and/or sell copies of the Software, and to permit persons to whom the      #
# Software is furnished to do so, subject to the following conditions:       #
#                                                                            #
# The above copyright notice and this permission notice shall be included in #
# all copies or substantial portions of the Software.                        #
#                                                                            #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    #
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        #
# DEALINGS IN THE SOFTWARE.                                                  #
#                                                                            #
##############################################################################
"""
**model.py**

Provides convenience functions for generating input and calculating
atomic charge densities for use with the Barbieri/Van Hove phase
shift calculation package.

"""

# cspell: ignore atmrad atorbs covrad malformatted nform nineq

from __future__ import print_function
from __future__ import division

import io
import os
from copy import deepcopy
from shutil import move
from glob import glob


from phaseshifts.core.atom import Atom
from phaseshifts.core.unit_cell import Unitcell
from phaseshifts import atorb
from phaseshifts.lib import libphsh


class MTZ_model(Model):
    """
    Muffin-tin potential Model subclass for producing input file
    for muffin-tin calculations in the Barbieri/Van Hove phase
    shift calculation package.
    """

    def __init__(self, unitcell, atoms, **kwargs):
        # type: (Unitcell, list[Atom], dict) -> None
        """
        Constructor for Model class.

        Parameters
        ----------
        unitcell : Unitcell
            An instance of the Unitcell class.
        atoms : list
            Array of Atom class instances which constitute the model.
        nh : int, optional
            Parameter for estimating the muffin-tin zero (default: 10).
        exchange : float, optional
            Hartree type exchange term alpha (default: 0.7200).
        c : float, optional
            The height of the slab in Angstroms - if this is much larger
            than the bulk c distance then there will be a large vacuum
            and therefore should be used when calculating a thin slab
            rather than a bulk muffin-tin potential. Default is to lookup
            the out-of-plane basis vector bulk value.
        nform : int or str, optional
            The phase shift calculation type, which can be 0 or 'cav' for
            using the cavpot subroutine, 1 or 'wil' for using the William's
            method, and 2 or 'rel' for using relativistic calculations suitable
            for the CLEED package.

        """
        self.atoms = []
        self.set_atoms(atoms)
        self.set_unitcell(unitcell)
        self.set_exchange(0.72)
        self.set_nh(10)
        self.mtz = None
        self.__dict__.update(kwargs)

    def set_nh(self, nh):
        """Sets the nh muffin-tin zero estimation parameter"""
        self.nh = int(nh)  # check this is not float

    def set_exchange(self, alpha):
        """Sets the alpha exchange term for muffin-tin calculation"""
        try:
            self.exchange = float(alpha)
        except:
            pass

    # set form of muffin-tin calculation: 0=cav, 1=wil, 2=rel
    def set_nform(self, nform):
        """
        Sets form of muffin-tin calculation

        Parameters
        ----------
        nform : int or str
            This governs the type of calculation, where nform can be:
            - **1.** "cav" or 0 - use Cavendish method
            - **2.** "wil" or 1 - use William's method
            - **3.** "rel" or 2 - Relativistic calculations

        """
        nform_lookup = {
            "cav": 0,
            "wil": 1,
            "rel": 2,
            0: 0,
            1: 1,
            2: 2,
        }
        self.nform = nform_lookup.get(nform, None) if nform in nform_lookup else int(nform)

    def set_slab_c(self, c):
        """
        Description
        -----------
        Set the maximum height of the slab in Angstroms - if this is
        much larger than the bulk c distance then there will be a large
        vacuum and therefore should be used when calculating a thin slab
        rather than a bulk muffin-tin potential.

        Examples
        --------
        For Re the bulk c distance is 2.76Å, whereas a possible slab c distance
        could be ~10Å.

        """
        try:
            self.c = float(c)
            self._c = self.c / 0.529
        except TypeError:
            pass

    def load_from_file(self, filename):
        # type: (str) -> None
        """
        Description
        -----------
        Load an cluster/slab input file and update the class instance variables

        Parameters
        ----------
        filename : str
            The path of the input file (e.g. cluster*.i or *slab*.i)

        Raises
        ------
        IOError : exception
            If the file cannot be read.
        TypeError : exception
            If a input line cannot be parsed correctly.

        """

        filename = glob(os.path.expanduser(os.path.expandvars(filename)))[0]

        def load_coordinates(line):
            """Load coordinates from line"""
            return list(map(float, line.split("#")[0].split()[:3]))

        try:
            with io.open(filename, mode="r", encoding="utf-8") as f:
                self.header = f.readline()
                a = float(f.readline().split("#")[0].split()[0]) * 0.529
                a1 = load_coordinates(f.readline())
                a2 = load_coordinates(f.readline())
                a3 = load_coordinates(f.readline())
                basis = [a1, a2, a3]
                c = float(a3[-1]) * 0.529  # change to Angstroms from Bohr
                self.set_unitcell(Unitcell(a, c, basis))
                self.nineq_atoms = int(f.readline().split("#")[0].split()[0])

                self.atoms = []
                line = f.readline().split("#")[0]
                while line.split()[0].isalpha():
                    n, protons, val, rad = [
                        t(s)
                        for (t, s) in zip(
                            (int, float, float, float),
                            f.readline().split("#")[0].split()[:4],
                        )
                    ]

                    self.atoms = [
                        Atom(
                            line.split()[0].capitalize(),
                            coordinates=load_coordinates(f.readline()),
                            tag=line.split()[1],
                            valence=val,
                            radius=rad,
                            protons=int(protons),
                        )
                        for _ in range(0, n)
                    ]
                    line = f.readline().split("#")[0]

                # print(line)
                self.set_nform(line.split()[0])
                self.set_exchange(f.readline().split("#")[0].split()[0])
                self.set_nh(f.readline().split("#")[0].split()[0])

        except IOError:
            raise IOError("Cannot read '%s'" % filename)

        except ValueError:
            raise ValueError("Malformatted input in '%s'" % filename)

    def create_atorbs(self, **kwargs):
        """
        Description
        -----------
        Create Atorb input files for each element in model.

        Arguments
        ---------
        output_dir : str
            Path to output directory for files
        library_dir : str
            Path to library of input files. Use this if you've previously
            created input files for a given element as it doesn't need to
            recalculate the radial charge density every time - note this
            is also a workaround to allow precalculated files if your machine's
            output is drastically different from what is expected.
        config : dict
            See help(atorb.gen_input) for list of keywords and values.

        Returns
        -------
        output_files : dict
            Dictionary list of atorb*.i input files for the Atorb class to
            calculate the charge density from.

        """

        # check output path
        output_dir = kwargs.get("output_dir") or kwargs.get("library_dir") or os.path.abspath(".")

        atorb_files = []
        at_files = []

        # generate input files for atomic orbitals and charge densities
        unique_elements = {str(atom.element.symbol) for atom in self.atoms}
        for element in unique_elements:
            atorb_dirpath = os.path.join(output_dir, "Atorb")
            os.makedirs(atorb_dirpath, exist_ok=True)
            atorb_file = os.path.join(atorb_dirpath, "atorb_%s.i" % element)
            if not os.path.isfile(atorb_file):  # create new atorb input file
                atorb_file = atorb.Atorb.gen_input(element, filename=atorb_file)

            at_file = os.path.join(output_dir, "Atorb", "at_" + element + ".i")
            if not os.path.isfile(at_file):
                at_file = atorb.Atorb.calculate_Q_density(
                    input=atorb_file, output_dir=os.path.join(output_dir, "Atorb")
                )

            # update lists
            atorb_files.append(atorb_file)
            at_files.append(at_file)

        return {"atorb_files": atorb_files, "at_files": at_files}

    def gen_atomic(self, **kwargs):
        """
        Description
        -----------

        Create atomic*.i input file for MTZ input based on model or
        list of files.

        Parameters
        ----------
        input_dir : str
            Input directory where at*.i files are kept.
        input_files : tuple
            List of input files to generate atomic input file from.
        output_file : str
            The filename of the resulting atomic*.i output file, which is
            simply a superimposed set of the radial charge densities from
            the individual input files.

        Returns
        -------
        output_file : str
            Returns the output file path string.

        Raises
        ------
        IOError : exception
            If either input directory or files do not exist.

        Notes
        -----
        If 'input_files' is not given then the default list of input files are
        inferred from the list of atoms in the model.

        """

        # input
        input_dir = (
            os.path.abspath(".")
            if "input_dir" not in kwargs
            else os.path.abspath(glob(os.path.expanduser(os.path.expandvars(kwargs["input_dir"])))[0])
        )

        if os.path.isfile(input_dir):
            input_dir = os.path.dirname(input_dir)

        if not os.path.isdir(input_dir):
            raise IOError("'%s' is not a valid input directory - does not exist!" % input_dir)

        # output filename
        output_file = str(os.path.abspath(kwargs.get("output_file") or "atomic.i"))

        # get list of input
        files = kwargs.get("input_files") or [  # assume using atoms from model
            os.path.join(input_dir, "at_" + atom.element.symbol + ".i") for atom in self.atoms
        ]

        # generate atomic.i input file by appending multiple at.i files
        with io.open(output_file, mode="w", encoding="ascii") as out_file_ptr:
            # loop through each atomic charge density file in list
            for input_file in files:
                if not os.path.isfile(str(input_file)) or input_file is None:
                    raise IOError("Radial charge density file " "'%s' does not exist!" % input_file)

                # append next input file to output
                with io.open(input_file, mode="r", encoding="utf-8") as input_file_ptr:
                    out_file_ptr.write(input_file_ptr.read())

        return output_file

    def get_MTZ(self, filename):
        """Retrieves muffin-tin potential from file"""
        try:
            with io.open(filename, mode="r", encoding="utf-8") as fptr:
                self.mtz = float(fptr.readline())  # read first line
        except IOError:
            raise IOError

    def calculate_muffin_tin_potential(self, mtz_string="", **kwargs):
        """
        Description
        -----------
        Calculate the muffin-tin potential (MTZ) for a given cluster file

        Parameters
        ----------
        atomic_file : str
            The path to the atomic input file. If this is omitted the default
            is generate one using the MTZ_model.gen_atomic() method.
        cluster_file : str
            The path to the cluster input file which may be a bulk or slab
            model.
        slab : int or bool
            Determines whether the MTZ calculation is for a slab model (True).
            The default is a bulk calculation.
        output : dict
            Dictionary output of 'mtz' - muffin-tin potential & 'output_file'
            - the path to the MTZ output file.


        Returns
        -------
        output_files : list(str)
            Paths to the MTZ output file.

        """

        # check to see if cluster input exists
        if "cluster_file" in kwargs:
            cluster_file = os.path.abspath(glob(os.path.expanduser(os.path.expandvars(kwargs["cluster_file"])))[0])
        else:
            cluster_file = os.path.abspath("cluster.i")

        if not os.path.isfile(cluster_file):
            raise IOError("MTZ cluster file '%s' does not exist!" % cluster_file)

        if not os.access(os.path.dirname(cluster_file), os.W_OK) and "atomic_file" not in kwargs:
            raise IOError("Do not have write access to '%s'" % os.path.dirname(cluster_file))

        # determine type of calculation - bulk or slab
        if "slab" in kwargs:
            slab = bool(kwargs["slab"])
            fid = "mtz"
        else:
            slab = False
            fid = "bmtz"
            mtz_string = ""

        if "atomic_file" in kwargs:
            atomic_file = os.path.abspath(kwargs["atomic_file"])
            if not os.path.isfile(atomic_file):
                raise IOError("Appended radial charge densities file " "'%s' does not exist!" % atomic_file)
        else:  # generate on the fly
            input_dir = os.path.abspath(os.path.dirname(cluster_file))
            self.create_atorbs(output_dir=input_dir)
            self.gen_atomic(input_dir=input_dir)

        if "output_file" in kwargs:
            output_file = os.path.abspath(kwargs["output_file"])
        else:
            output_file = "%s.i" % fid

        try:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
        except OSError:
            pass

        if os.path.isfile(output_file):
            move(output_file, output_file + ".bak")

        # create mufftin debug file
        mufftin_file = os.path.splitext(output_file)[0] + "_mufftin.d"
        info_file = os.path.splitext(output_file)[0] + "_info.txt"

        # call cavpot routine
        self.mtz = libphsh.cavpot(
            mtz_string,
            int(slab),
            atomic_file,
            cluster_file,
            mufftin_file,
            output_file,
            info_file,
        )

        # check to see if new file has been written
        if not os.path.isfile(output_file):
            raise IOError("Failed to write muffin-tin potential file '%s'" % output_file)

        return output_file

    def gen_input(self, **kwargs):
        """
        Description
        -----------
        Generate input file for use in the Barbieri/Van Hove phase shift
        calculation package (phsh1 or libphsh)

        Arguments
        ---------
        bulk : bool
            If True the c value is set to the bulk value...
        filename : str
            The path of the generated file (e.g. 'cluster.i')
        header : str
            A file header string.
        bulk : bool
            Specify whether generated file is a bulk or slab model
        pos_check : bool
            If true, the positions of each atom will be checked and
            any duplicates will cause a CoordinateError to be raised.
        nform : int or str
            Choose calculation type from either:
                - Cavendish (0)
                - Williams (1)
                - Relativistic (2)
        exchange : float
            Hartree-Fock exchange term alpha.
        nh : int
            Term used to estimate the muffin-tin zero.

        Returns
        -------
        filename on success

        Raises
        ------
        CoordinatesError : exception
            if the model atoms have duplicate coordinates and the 'pos_check'
            kwarg is set to True.

        """
        # create backup of atoms
        mtz_atoms = deepcopy(self.atoms)

        # get file header - should be single line (force if not)
        header = kwargs.get("header", "MTZ cluster input file").replace("\n", " ").replace("\r", "")

        if "bulk" in kwargs:
            if kwargs["bulk"]:
                fid = "_bulk"
            else:
                fid = "_slab"
        else:
            fid = ""

        self.set_nform(kwargs.get("nform", 2))

        if "exchange" in kwargs:
            self.set_exchange(kwargs["exchange"])
        elif not isinstance(self.exchange, float):
            self.set_exchange(0.72)

        if "nh" in kwargs:
            self.set_nh(kwargs["nh"])
        elif not isinstance(self.nh, int):
            self.set_nh(10)

        # get filename or make educated guess
        filename = kwargs.get("filename") or "cluster%s.i" % fid

        # determine if coordinate duplicates
        if "pos_check" in kwargs:
            if bool(kwargs["pos_check"]):
                self.check_coordinates()

        # auto-generate file
        a = float(self.unitcell._a) * 0.529  # pylint: disable=protected-access
        with io.open(filename, mode="w", encoding="ascii") as out_file_ptr:

            def write_line(data, comment="", justify=33, end="\n"):
                """Write line to file"""
                comment = ("# " + comment.lstrip("# ")) if comment else ""
                out_file_ptr.write(str(data).ljust(justify) + comment + end)

            out_file_ptr.write(header + "\n")
            write_line(data=" %7.4f" % self.unitcell._a, comment="a lattice parameter distance in Bohr radii")

            comment3 = "Notice the value %.2f (%s calculation)" % (self.unitcell.basis[2][2] / a, fid.replace("_", ""))
            for (
                i,
                comment,
            ) in zip(range(3), ["SPA coordinates of unit cell", "", comment3]):
                write_line(
                    data=" %7.4f %7.4f %7.4f" % tuple(self.unitcell.basis[i][j] / a for j in range(3)),
                    comment=comment,
                )

            # TODO: better nineq_atoms prediction
            try:
                nineq_atoms = int(kwargs.get("nineq_atoms", len(set(self.atoms))))
            except TypeError:
                nineq_atoms = self._nineq_atoms()[0]

            tags = []

            # check to see if nineq_atoms is estimated in code
            if isinstance(nineq_atoms, tuple):
                write_line(data="%4i" % nineq_atoms[0], comment="number of inequivalent atoms in this file (NINEQ)")

                # now loop through each inequivalent atom and add to file
                elements_dict = nineq_atoms[1]
                for element in elements_dict:
                    atoms = elements_dict.get(element)["atom_list"]

                    # Loop through each inequivalent atom type
                    for inequivalent_atom in set(atoms):
                        # get list of atoms of this type
                        # i.e. same element, radius & valence
                        inequivalent_atoms = [atom for atom in atoms if atom == inequivalent_atom]
                        inequivalent_tags = {atom.tag for atom in inequivalent_atoms if atom.tag not in tags}

                        # select first unused tag from list
                        for tag in inequivalent_tags:
                            if tag not in tags:
                                inequivalent_atom.tag = tag
                                break

                        # avoid duplicate tags for different atoms
                        while inequivalent_atom.tag in tags:
                            number = "".join([ch for ch in inequivalent_atom.tag if ch.isdigit()])
                            try:
                                number = int(number)
                                number += 1
                            except ValueError:
                                number = 1
                            inequivalent_atom.tag = "".join(
                                [ch for ch in inequivalent_atom.tag if ch.isalpha() or ch in ["_", "-", "+"]]
                            ) + "_{}".format(number)

                        write_line(
                            data="{0} {1}".format(inequivalent_atom.element.name.capitalize(), inequivalent_atom.tag),
                            comment="element, name tag",
                        )
                        write_line(
                            data="%4i %7.4f %7.4f %7.4f"
                            % (
                                len(inequivalent_atoms),
                                inequivalent_atom.protons,
                                inequivalent_atom.valence,
                                inequivalent_atom.radius / 0.529,
                            ),
                            comment="atoms in unit cell, Z, valence, Muffin-tin radius (Bohr radii)",
                        )

                        # write each inequivalent atom list to file
                        for atom in inequivalent_atoms:
                            write_line(
                                data=" %7.4f %7.4f %7.4f" % tuple(atom.coordinates[i] / a for i in range(3)),
                                comment="coordinates in SPA units",
                            )

            else:  # assume each element is an inequivalent atom
                write_line(data="%4i" % nineq_atoms, comment="number of inequivalent atoms in this file (NINEQ)")

                for atom in set(self.atoms):
                    while atom.tag in tags:
                        number_str = "".join([ch for ch in atom.tag if ch.isdigit()])
                        try:
                            number = int(number_str)
                            number += 1
                        except ValueError:
                            number = 1
                        atom.tag = (
                            "".join([ch for ch in atom.tag if ch.isalpha() or ch in ["_", "-", "+"]])
                            + "_"
                            + str(number)
                        )

                    tags.append(atom.tag)
                    write_line("%s %s" % (atom.element.name.capitalize(), atom.tag), comment="# element, name tag")
                    data = (tags.count(atom.tag), atom.protons, atom.valence, atom.radius / 0.529)
                    write_line(
                        data="%4i %7.4f %7.4f %7.4f" % data,
                        comment="# atoms in unit cell, Z, valence, Muffin-tin radius (Bohr radii)"
                    )
                    write_line(
                        data=" %7.4f %7.4f %7.4f" % tuple(atom.coordinates[i] / a for i in range(3)),
                        comment="# coordinates in SPA units\n"
                    )

            write_line(data="%4i" % self.nform, comment="# nform=2|1|0 (for rel, will or cav)")
            write_line(data=" %7.4f" % self.exchange, comment="# alpha (for Hartree type exchange term)")
            write_line("%4i" % self.nh, "# nh (for estimating muffin-tin zero)")

        # restore backup
        self.atoms = mtz_atoms

    def get_elements(self):
        """Return the unique elements in model"""
        return {atom.name for atom in self.atoms}


# #==============================================================================
# # Testing
# #==============================================================================
# at = Atom('C', [0, 0, 0])
# ab = Atom('Re', [0, 0, 0], tag='Re2')
# ac = Atom('Re', [0, 0, 0], tag='Re1')
# uc = Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
# print(uc)
# print(set([at, ab, ac]))
#
# mtz = MTZ_model(uc, atoms=[at, ab, ac])
# mtz.load_from_file('C:\\Users\\Liam\\Dropbox\\Programming\\Python\\LEED-PyV\\phaseshifts\\test\\Re0001\\cluster_Re_bulk.i')
# print(mtz.get_elements())
# mtz.load_from_file('C:\\Users\\kss07698\\Desktop\\test_cluster.bak.i')
# mtz.gen_input(filename='C:\\Users\\kss07698\\Desktop\\test_cluster.bak.i')
