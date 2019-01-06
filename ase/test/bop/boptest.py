import numpy as np
from ase.bop.bopatoms import BOPAtoms
from ase.neighborlist import NeighborList


onsite_level = 0.5
struc = BOPAtoms(symbols='Fe', positions=[[0, 0, 0]], onsite_levels=[onsite_level], cell=1 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), pbc=(1, 1, 1))
assert struc[0].onsite_level == onsite_level

print(struc[0].number_valence_electrons)

cutoff = 1.01
nl = NeighborList([cutoff/2]*len(struc), skin=0.0)
nl.update(struc)

print(nl.get_neighbors(0))