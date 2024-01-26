"""Define the Unitcell class."""
from .atom import Atom


class Unitcell(object):  # pylint: disable=too-many-instance-attributes
    """Unitcell class for representing a crystal structure."""

    def __init__(self, a, c, matrix_3x3, **kwargs):
        """
        Constructor for the Unitcell class

        Parameters
        ----------
        a : float
            The in-plane lattice vector in Angstroms
        c : float
            The out-of-plane lattice vector in Angstroms. For cubic systems this
            will be equal to a.
        matrix_3x3: ndarray
            A 3x3 matrix describing the x,y,z construction of a,b,c basis vectors
            of the unitcell. Units for x, y & z should be in terms of fractional
            coordinates.
        alpha : float, optional
            Angle alpha in degrees.
        beta : float, optional
            Angle beta in degrees.
        gamma : float, optional
            Angle gamma in degrees.

        """
        # Convert Angstrom input to Bohr radii
        self.set_a(a)
        self.set_c(c)

        # Set basis vectors
        self.set_vectors(matrix_3x3)

        # Set crystal system
        self.alpha = 90.0
        self.beta = 90.0
        self.gamma = 90.0

        # Update additional information
        self.__dict__.update(kwargs)

    # checks if two class instances are equal w.r.t. name, radius & valence
    def __eq__(self, other):
        return isinstance(other, Atom) and (
            self.a == other.a
            and self.c == other.c
            and self.alpha == other.alpha
            and self.beta == other.beta
            and self.gamma == other.gamma
        )

    # checks if two class instances are not equal w.r.t. name, radius & valence
    def __neq__(self, other):
        return not self.__eq__(other)

    # reprinting of class object
    def __repr__(self):
        return str("{}(a={}, c={}, alpha={}, beta={}, gamma={}, basis={})").format(
            self.__class__.__name__,
            self.a,
            self.c,
            self.alpha,
            self.beta,
            self.gamma,
            self.basis,
        )

    # redefine hash method for checking uniqueness of class instance
    def __hash__(self):
        return hash(self.__repr__())

    # set basis vectors from (3x3) matrix in fractional coordinates
    def set_vectors(self, m3x3):
        self.basis = m3x3

    # set a lattice parameter
    def set_a(self, a):
        """
        Description
        -----------
        Set the magnitude of the in-plane lattice vector a in Angstroms.

        Parameters
        ----------
        a: float
            The magnitude of the in-plane lattice vector in Angstroms

        Notes
        -----
        To retrieve a in terms of Angstroms use 'unitcell.a', whereas the
        internal parameter 'unitcell._a' converts a into Bohr radii
        (1 Bohr = 0.529Å), which is used for the muffin-tin potential
        calculations in libphsh (CAVPOT subroutine).

        """
        self.a = float(a)
        self._a = self.a / 0.529  # (1 Bohr = 0.529Å)

    # set c lattice parameter
    def set_c(self, c):
        """
        Description
        -----------
        Set the magnitude of the out-of-plane lattice vector a.

        Parameters
        ----------
        c : float
            The magnitude of the in-plane lattice vector in Angstroms

        Notes
        -----
        To retrieve c in terms of Angstroms use 'unitcell.c', whereas the
        internal parameter 'unitcell._c' converts c into Bohr radii
        (1 Bohr = 0.529Å), which is used for the muffin-tin potential
        calculations in libphsh (CAVPOT subroutine).

        """
        self.c = float(c)
        self._c = self.c / 0.529  # (1 Bohr = 0.529Å)

    @property
    def alpha(self):
        # type: () -> float
        """Angle alpha in degrees."""
        return self._alpha_degrees

    @alpha.setter
    def alpha(self, alpha):
        self._alpha_degrees = float(alpha) % 360.0

    @property
    def beta(self):
        # type: () -> float
        """Angle beta in degrees."""
        return self._beta_degrees

    @beta.setter
    def beta(self, beta):
        self._beta_degrees = float(beta) % 360.0

    @property
    def gamma(self):
        # type: () -> float
        """Angle gamma in degrees."""
        return self._gamma_degrees

    @gamma.setter
    def gamma(self, gamma):
        self._gamma_degrees = float(gamma) % 360.0
