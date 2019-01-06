import numpy as np
from ase.bop.bopatoms import BOPAtoms
from ase.neighborlist import NeighborList
from ase.bop.buildhop import TwoCenterHoppingIntegrals


onsite_level = 0.5
struc = BOPAtoms(symbols='Fe', positions=[[0., 0., 0.]], onsite_levels=[onsite_level],
                 cell=1 * np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]), pbc=(1, 1, 1))
assert struc[0].onsite_level == onsite_level

print(struc[0].number_valence_electrons)

cutoff = 0.6
cutofflist = [cutoff]*len(struc)

twocenterhops = TwoCenterHoppingIntegrals(bopatoms=struc, cutoffs=cutofflist, skin=0.0)
print(twocenterhops.nl.get_neighbors(0))

rotnum = 2
z_axis = np.array([0., 0., 1.])
z_axis = z_axis / np.linalg.norm(z_axis)
rel_pos = twocenterhops.get_relative_position(0, rotnum)
rotated_rel_pos = twocenterhops.get_rotation(0, rotnum, z_axis_global=z_axis).apply(rel_pos)
print(z_axis)
print(rotated_rel_pos)
