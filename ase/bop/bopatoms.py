from ase.atoms import Atoms
from ase.bop.bopatom import BOPAtom


class BOPAtoms(Atoms):

    def __init__(self, symbols=None, *args, **kwargs):
        super().__init__(symbols, *args, **kwargs)
        if (isinstance(symbols, (list, tuple)) and
              len(symbols) > 0 and isinstance(symbols[0], BOPAtom)):
            # Get data from a list or tuple of Atom objects:
            data = [[atom.get_raw(name) for atom in symbols]
                    for name in
                    ['position', 'number', 'tag', 'momentum',
                     'mass', 'magmom', 'charge', 'onsite_level']]
            atoms = self.__class__(None, *data)
            symbols = None
