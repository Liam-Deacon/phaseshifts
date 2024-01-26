"""Defines a basic model class."""

# cspell: ignore nineq

from ..errors import CoordinatesError
from .atom import Atom
from .unit_cell import Unitcell


class Model(object):
    """Generic model class."""

    def __init__(self, unitcell, atoms, **kwargs):  # type: ignore[misc,syntax]
        # type: (Unitcell, list[Atom], object) -> None
        """
        Constructor for Model class.

        Parameters
        ----------
        unitcell : Unitcell
            An instance of the Unitcell class.
        atoms : list
            Array of Atom class instances which constitute the model.

        """
        self.atoms = []  # type: list[Atom]
        self.set_atoms(atoms)
        self.unitcell = unitcell
        self.__dict__.update(kwargs)

    # checks if two models are equal
    def __eq__(self, other):
        is_equal = False
        if isinstance(other, Model):
            is_equal = self.atoms == other.atoms and self.unitcell == other.unitcell
        return is_equal

    # checks if two models are not equal
    def __neq__(self, other):
        return not self.__eq__(other)

    # reprinting of Atom object
    def __repr__(self):
        return "{}(atoms={!r}, unitcell={!r})".format(self.__class__.__name__, self.atoms, self.unitcell)

    def __hash__(self):
        """Redefine hash method for checking uniqueness of class instance."""
        return hash(self.__repr__())

    # estimate number of inequivalent atoms
    def _nineq_atoms(self):
        """
        Description
        -----------
        Internal method for estimating the number of inequivalent atoms

        Returns
        -------
        nineq_atoms, element_dict : tuple
            nineq_atoms : The estimated number of inequivalent atoms based on
                the valence and radius of each atom.
            element_dict : a dictionary of each element in the atom list where
                each element contains an atom dictionary of 'nineq_atoms',
                'n_atoms' and a complete 'atom_list'

        Example
        -------
        >>> C1 = Atom('C', [0, 0, 0])
        >>> Re1 = Atom('Re', [0, 0, 0], valence=2.0)
        >>> Re2 = Atom('Re', [0, 0, 0], radius=1)
        >>> uc = Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        >>> mtz = MTZ_model(uc, atoms=[C1, Re1, Re2])
        >>> print(mtz._nineq_atoms())
        (3, {'Carbon': {'n_atoms': 1, 'atom_list': [Atom(Carbon, tag='C',
        coordinates=[0, 0, 0], Z=6, radius=0.91, valence=0)], 'nineq_atoms':
        1,'nineq_atoms_list': set([Atom(Carbon, tag='C',
        coordinates=[0, 0, 0], Z=6, radius=0.91, valence=0)])}, 'Rhenium': {
        'n_atoms': 2, 'atom_list': [Atom(Rhenium, tag='Re',
        coordinates=[0, 0, 0], Z=75, radius=1.97, valence=2.0), Atom(Rhenium,
        tag='Re', coordinates=[0, 0, 0], Z=75, radius=1, valence=0)],
        'nineq_atoms': 2, 'nineq_atoms_list': set([Atom(Rhenium, tag='Re',
        coordinates=[0, 0, 0], Z=75, radius=1.97, valence=2.0),
        Atom(Rhenium, tag='Re', coordinates=[0, 0, 0], Z=75, radius=1,
        valence=0)])}})

        """
        nineq_atoms = 0
        element_dict = {}
        atom_dict = {}
        # loop through atom list, testing each element for duplicates
        # get list of elements
        elements = {atom.name for atom in self.atoms}
        for element in elements:
            atoms = [atom for atom in self.atoms if atom.name == element]
            n_atoms = len(set(atoms))
            nineq_atoms += n_atoms
            atom_dict = {
                "nineq_atoms": n_atoms,
                "n_atoms": len(atoms),
                "atom_list": atoms,
            }
            element_dict[element] = atom_dict
        return nineq_atoms, element_dict

    def add_element(self, element, position, **kwargs):    # type: ignore[misc,syntax]
        # type: (int|str, list[float], dict) -> None
        """
        Append an Atom instance to the model

        Parameters
        ----------
        element : str or int
            Either an element name, symbol or atomic number.
        position : list(float) or ndarray
            (1x3) array of the fractional coordinates of the atom
            within the unit cell in terms of the lattice vector a.

        """
        atom = Atom(element, position, **kwargs)
        self.atoms.append(atom)

    def check_coordinates(self):
        """
        Check for duplicate coordinates of different atoms in model.

        Raises
        ------
        CoordinateError : exception
            If duplicate positions found.

        """
        positions = [str(atom.coordinates) for atom in self.atoms]
        info = ""
        for position in {position for position in positions if positions.count(position) > 1}:
            for i, atom in enumerate([atom for atom in self.atoms if str(atom.coordinates) == position]):
                info += "%s, coordinates=%s, index=%i\n" % (
                    str(atom),
                    atom.coordinates,
                    i,
                )
        if len(set(positions)) < len(self.atoms):
            raise CoordinatesError("Not every atom position in model is unique!\n%s\n" % info)

    def set_atoms(self, atoms):
        # type: (list[Atom]) -> None
        """
        Set the atoms for the model.

        Parameters
        ----------
        atoms : list(Atoms)
            Array of Atom instances. Entries in the list which are
            not of type Atom will be ignored.

        Raises
        ------
        TypeError : exception
            If atoms parameter is not a list.

        """
        def _check_atom(atom):
            if not isinstance(atom, Atom):
                raise TypeError("atom {!r} must be an Atom instance, got '{}'".format(atom, type(atom)))
            return atom

        self.atoms = list(map(_check_atom, atoms))

    @property
    def unitcell(self):  # type: () -> Unitcell
        """Unit cell for the model."""
        return self._unitcell

    @unitcell.setter
    def unitcell(self, unitcell):
        # type: (Unitcell) -> None
        """
        Set the unitcell for the model

        Parameters
        ----------
        unitcell : Unitcell
            Instance of Unitcell class to set to model.

        Raises
        ------
        TypeError : exception
            If unitcell parameter is not a Unitcell.

        """
        if isinstance(unitcell, Unitcell):
            self._unitcell = unitcell
        else:
            raise TypeError("unitcell parameter must be of type Unitcell")
