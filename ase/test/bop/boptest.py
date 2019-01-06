from ase.bop.bopatom import BOPAtom
from ase.bop.bopatoms import BOPAtoms

from ase.atoms import Atoms
from ase.atom import Atom

#b = Atoms([Atom('Fe', [0,0,0])])

b = BOPAtoms(symbols='Fe', positions=[[0, 0, 0]], onsite_levels=[0.5])
print(b[0].onsite_level)