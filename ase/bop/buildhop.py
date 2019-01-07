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

    def get_single_hop_global(self, atomindex: int, jneigh: int):
        hop = self.get_single_hop_local(atomindex, jneigh)
        rot = self.get_rotation(atomindex, jneigh)
        return hop


    def get_single_hop_local(self, atomindex: int, jneigh: int):
        bopatom_i = self.bopatoms[atomindex]
        (bopatom_neighbors_i, scaled_pos_list) = self.nl.get_neighbors(atomindex)
        bopatom_j_neigh_i = self.bopatoms[bopatom_neighbors_i[jneigh]]
        pos_list = np.dot(scaled_pos_list, self.bopatoms.get_cell())
        rel_pos_ij = pos_list[jneigh]
        r_ij = np.linalg.norm(rel_pos_ij)

        ss_sigma = 0
        sp_sigma = 0
        sd_sigma = 0
        pp_sigma = 0
        pp_pi = 0
        pd_sigma = 0
        pd_pi = 0
        dd_sigma = 0
        dd_pi = 0
        dd_delta = 0
        slater_koster_matrix = np.zeros((9, 9))
        # initialize bond integrals
        if bopatom_i.number_valence_orbitals == 5 and bopatom_j_neigh_i.number_valence_orbitals == 5:
            dd_sigma = np.exp(-r_ij * 1)
            dd_pi = np.exp(-r_ij * 2)
            dd_delta = np.exp(-r_ij * 3)
            slater_koster_matrix[4, 4] = dd_delta
            slater_koster_matrix[5, 5] = dd_pi
            slater_koster_matrix[6, 6] = dd_pi
            slater_koster_matrix[7, 7] = dd_delta
            slater_koster_matrix[8, 8] = dd_sigma
        return slater_koster_matrix[4:, 4:]


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

    def get_rotation(self, atomindex: int, jneigh:int, z_axis_global: np.array=np.array([0, 0, 1])) -> Rotation:
        rel_pos = self.get_relative_position(atomindex, jneigh)
        v = np.cross(rel_pos, z_axis_global)
        cosine = np.dot(rel_pos, z_axis_global)
        if cosine != -1:
            v_cross = [[0    , -v[2], v[1] ],
                       [v[2] , 0    , -v[0]],
                       [-v[1], v[0] , 0]]
            rotation_matrix = np.eye(3) + v_cross + np.dot(v_cross, v_cross) / (1 + cosine)
        else:
            rotation_matrix = np.diag([1, 1, -1])
        return Rotation.from_dcm(rotation_matrix)
