from ase.atom import Atom

from ase.atom import Atom, atomproperty, names


names['onsite_level'] = ('onsite_levels', 0.0)


class BOPAtom(Atom):
    onsite_level = atomproperty('onsite_level', 'Atomic onsite level')

    def __init__(self, onsite_level=None, atoms=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if atoms is None:
            self.data['onsite_level'] = onsite_level
