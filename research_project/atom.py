from mendeleev import element

element_mass = {}
def get_element_mass(atom):
    #atom.capitalize()
    if element_mass.get(atom) == None:
        element_mass[atom] = element(atom).atomic_weight
    return element_mass[atom]

class Atom:
    def __init__(self, symbol, atom_id, name, coords, resid, protein):
        self.symbol = symbol
        self.symbol.capitalize()
        self.atom_id = atom_id
        self.name = name
        self.atomic_mass = get_element_mass(self.symbol)
        self.coords = coords
        self.resid = resid
        self.protein = protein

    def get_symbol(self):
        return self.symbol

    def get_atom_id(self):
        return self.atom_id

    def get_mass(self):
        return self.atomic_mass

    def get_coords(self):
        return self.coords

    def get_resid(self):
        return self.resid

    def get_protein(self):
        return self.protein

def atom_main():
    atm = Atom("C", 1, "CA", (0, 0, 0), 2, "test")
    print(atm.get_symbol(), atm.get_atom_id(), atm.get_mass(), atm.get_coords(), atm.get_resid(), atm.get_protein())
