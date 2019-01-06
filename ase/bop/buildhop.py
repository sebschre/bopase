import numpy as np
from typing import List
from ase.bop.bopatoms import BOPAtoms
from ase.neighborlist import NeighborList
from scipy.spatial.transform import Rotation


class TwoCenterHoppingIntegrals:
    def __init__(self, bopatoms: BOPAtoms, cutoffs: list, **kwargs):
        self.bopatoms = bopatoms
        self.nl = NeighborList(cutoffs=cutoffs, bothways=True, self_interaction=False, **kwargs)
        self.nl.update(bopatoms)

    def update_hops(self):
        raise NotImplemented

    def get_relative_position(self, index: int, jneigh: int):
        '''
        :param index: atom index
        :param jneigh: indexing neighboring atoms of atom index
        :return:
        '''
        # if index > len(self.bopatoms):
        #     raise IndexError
        # if jneigh > self.nl.nneighbors - 1:
        #     raise IndexError
        (index_list, relative_position_list) = self.nl.get_neighbors(index)
        return relative_position_list[jneigh]

    def get_rotation(self, index: int, jneigh:int, z_axis_global: np.array=np.array([0, 0, 1])) -> Rotation:
        rel_pos = self.get_relative_position(index, jneigh)
        v = np.cross(rel_pos, z_axis_global)
        sine   = np.abs(v)
        cosine = np.dot(rel_pos, z_axis_global)
        if cosine != -1:
            v_cross = [[0    , -v[2], v[1] ],
                       [v[2] , 0    , -v[0]],
                       [-v[1], v[0] , 0]]
            rotation_matrix = np.eye(3) + v_cross + np.dot(v_cross, v_cross) / (1 + cosine)
        else:
            rotation_matrix = np.diag([1, 1, -1])
        return Rotation.from_dcm(rotation_matrix)
