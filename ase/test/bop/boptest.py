import numpy as np
from ase.bop.bopatoms import BOPAtoms
from ase.neighborlist import NeighborList
from ase.bop.buildhop import TwoCenterHoppingIntegrals


lattice_constant = 2

onsite_level = 0.5
struc = BOPAtoms(symbols='Fe', positions=[[0., 0., 0.]], onsite_levels=[onsite_level],
                 cell=lattice_constant * np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]), pbc=(1, 1, 1))

assert struc[0].onsite_level == onsite_level


cutoff = lattice_constant / 2 + 0.1
cutofflist = [cutoff]*len(struc)

twocenterhops = TwoCenterHoppingIntegrals(bopatoms=struc, cutoffs=cutofflist, skin=0.0)
print(twocenterhops.nl.get_neighbors(0))

rotnum = 2
z_axis = np.array([0., 0., 1.])
z_axis = z_axis / np.linalg.norm(z_axis)
rel_pos = twocenterhops.get_relative_position(0, rotnum)
rotated_rel_pos = twocenterhops.get_rotation(0, rotnum, z_axis_global=z_axis).apply(rel_pos)

print(twocenterhops.get_single_hop_local(0, 3))
print(twocenterhops.get_single_hop_global(0, 3))

print(twocenterhops.get_dbond_rotation_matrix(np.pi, 0))
