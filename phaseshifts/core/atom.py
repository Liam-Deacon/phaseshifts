"""Define an Atom class."""
# cspell: ignore atmrad covrad mufftin

from __future__ import print_function
from __future__ import division

from phaseshifts.core import elements


class Atom(object):  # pylint: disable=too-many-instance-attributes
    """
    Atom class for input into cluster model for muffin-tin potential
    calculations.
    """

    def __init__(  # type: ignore[misc,syntax]
        self, element, coordinates=(0.0, 0.0, 0.0), radius=None, valence=None, occupancy=None, **kwargs
    ):
        # type: (int|str, tuple[float, float, float], float|None, float|None, float|None, object) -> None
        """
        Constructor for Atom class.

        Parameters
        ----------

        element : str or int
            This is either the elemental symbol, name or atomic number.
        coordinates : list[float, float, float] or :code:`ndarray`
            The fractional coordinates within the unitcell in terms of the
            basis vector a.
        tag : str, optional
            Add a name tag to this element (useful if dealing with multiple
            atoms of the same element in a given model). Default is the
            symbol for that element - numeric ids may be appended in the model
            class.
        radius : float, optional
            The muffin-tin radius of the atom in Angstroms (default is to
            lookup 'atmrad' in the element dictionary).
        valence : int, optional
            The valency of the atom (default is to assume neutral atom).
        occupancy : float, optional
            The fractional occupancy of the atom in the given site.

        """
        self.element = elements.ELEMENTS[element]
        self._coordinates = None  # dummy
        self.coordinates = coordinates
        self.name = self.element.name.title()  # type: str
        self.tag = self.element.symbol.title()  # type: str
        self.protons = self.element.protons  # type: int
        self.valence = valence or 0
        # NOTE: assume covrad for non-zero valency
        self.radius = self.element.covrad if valence and radius is None else self.element.atmrad  # type: float
        self.occupancy = occupancy or 1.0  # type: float
        self.__dict__.update(kwargs)

    def __eq__(self, other):
        """Checks whether two atoms are equal w.r.t. name, radius and valence."""
        is_equal = False
        if isinstance(other, Atom):
            is_equal = (
                self.name == other.name
                and self.radius == other.radius
                and self.valence == other.valence
            )
        return is_equal

    def __neq__(self, other):
        """Checks whether two atoms are not equal w.r.t. name, radius and valence."""
        return not self.__eq__(other)

    def __repr__(self):
        """Define Atom representation as string."""
        return "{}({!r}, tag={!r}, radius={!r}, " "valence={})".format(
            self.__class__.__name__,
            self.name,
            self.tag,
            self.radius,
            self.valence,
        )

    def __hash__(self):
        """Redefine hash method for checking uniqueness of class instance."""
        return hash(self.__repr__())

    @property
    def coordinates(self):
        """Fractional coordinates of the unit cell in terms of basic vector ``a``."""
        return self._coordinates

    @coordinates.setter
    def coordinates(self, coordinates):
        """Set coordinates of atom within unitcell in terms of ``a``."""
        self.coordinates = coordinates
        self._coordinates = tuple(r / 0.529 for r in map(float, coordinates))  # type: tuple[float, float, float]

    @property
    def valence(self):  # type: () -> float
        """The valency of the atom, e.g. 0 would indicate atomic, 2 would be an example of an ion."""
        return self._valence

    @valence.setter
    def valence(self, valency):
        # type: (float|int) -> None
        """Sets the valency of the atom."""
        self._valence = float(valency)

    def is_ion(self):  # type: () -> bool
        """Returns True if the atom is an ion."""
        return abs(self.valence) > 0

    def set_mufftin_radius(self, radius):
        # type: (float) -> None
        """
        Sets the muffin-tin radius of the atom in Angstroms.
        """
        self.radius = float(radius)
        self._radius = self.radius / 0.529  # in Bohr radii
