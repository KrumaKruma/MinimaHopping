
class Minimum():
    def __init__(self, atoms, n_visit, fingerprint, T, ediff, acc_rej):
        self.atoms = atoms
        self.fingerprint = fingerprint
        self.temperature = T
        self.ediff = ediff
        self.acc_rej = acc_rej
        self.n_visit = n_visit

    def __lt__(self, other):
        return self.atoms.get_potential_energy() < other.atoms.get_potential_energy()

    def __gt__(self, other):
        return self.atoms.get_potential_energy() > other.atoms.get_potential_energy()

