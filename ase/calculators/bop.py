"""Analytic Bond-Order Potential Calculator"""

from ase.calculators.calculator import Calculator, all_changes
from ase.neighborlist import NeighborList


class BOPcalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)
        self.neighborlist = None

    def get_hops(self, atoms):
        pass

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        self.neighborlist = NeighborList([5] * len(atoms), self_interaction=False)
