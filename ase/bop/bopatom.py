from ase.atom import Atom, atomproperty, names, chemical_symbols


names['onsite_level'] = ('onsite_levels', 0.0)


class BOPAtom(Atom):
    onsite_level = atomproperty('onsite_level', 'Atomic onsite level')

    def __init__(self, symbol='X', position=(0, 0, 0),
                 tag=None, momentum=None, mass=None,
                 magmom=None, charge=None,
                 atoms=None, index=None,
                 onsite_level=None):
        super().__init__(symbol=symbol, position=position,
                 tag=tag, momentum=momentum, mass=mass,
                 magmom=magmom, charge=charge,
                 atoms=atoms, index=index)

        self.data['onsite_level'] = onsite_level

    def get_raw(self, name):
        """Get name attribute, return None if not explicitely set."""
        if name == 'symbol':
            return chemical_symbols[self.get_raw('number')]

        if self.atoms is None:
            return self.data[name]

        plural = names[name][0]
        if plural in self.atoms.arrays:
            return self.atoms.arrays[plural][self.index]
        else:
            return None