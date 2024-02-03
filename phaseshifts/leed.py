"""Provides adaptors and validators when working with LEED files."""

import io

from phaseshifts.adaptors.cleed_adaptors import *  # noqa
from phaseshifts.validators.cleed_validators import *  # noqa


class CSearch:
    """class for csearch related data exchange"""

    def __init__(self, model_name, leed_command=None):
        self.set_model(model_name)
        self._get_results()

    def set_model(self, name):
        """sets the model name"""
        self.model = str(name)

    def get_iteration(self, iteration):
        try:
            with io.open(self.model + ".log", mode="r", encoding="ascii") as file_ptr:
                return [line for line in file_ptr if line.startswith("#") and "par" in line][
                    int(iteration)
                ]
        except (IndexError, ValueError, OSError):
            return None

    def get_r_factor(self, iteration):
        raise NotImplementedError

    def get_time_stamp(self, iteration):
        raise NotImplementedError

    def _get_results(self):
        raise NotImplementedError
