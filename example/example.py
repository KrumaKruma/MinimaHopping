
from ase.calculators.eam import EAM
from minimahopping.minhop import Minimahopping
from ase.cluster.wulff import wulff_construction


def main():
    # construct a chain with 13 atoms:
    # atoms = Icosahedron('Na', 2, latticeconstant=None)

    atoms = wulff_construction('Na', surfaces=[(1, 0, 0), (0, 1, 0),(0, 0, 1)], energies=[0.001, 0.001, 0.25],
                            size=55, # maximum number of atoms
                            structure='bcc', rounding='above')


    calculator = EAM(potential="Na_v2.eam.fs")
    atoms.calc = calculator
    mh = Minimahopping( atoms,
                        T0 = 2000,
                        beta_decrease = 1./1.05,
                        beta_increase = 1.05,
                        Ediff0 = .01,
                        alpha_a = 0.95,
                        alpha_r = 1.05,
                        n_soft = 20,
                        ns_orb = 1,
                        np_orb = 1,
                        width_cutoff = 3.5,
                        maxnatsphere = 100,
                        exclude = [],
                        dt = 0.1,
                        mdmin = 2,
                        fmax = 0.001,
                        enhanced_feedback = False,
                        energy_threshold = 0.1, #5 the noise
                        n_poslow = 30,
                        minima_threshold = 1e-3,
                        verbose = True,
                        new_start = False,
                        run_time = 'infinit',
                        use_intermediate_mechanism = False,
                        restart_interval = 2)
    mh(totalsteps=100)






if __name__ == '__main__':
    main()

