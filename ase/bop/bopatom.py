from ase.atom import Atom, atomproperty, names, chemical_symbols


names['onsite_level'] = ('onsite_levels', 0.0)


class BOPAtom(Atom):
    onsite_level = atomproperty('onsite_level', 'Atomic onsite level')

    def __init__(self, bopatoms=None, index=None):
        super().__init__(atoms=bopatoms, index=index)
