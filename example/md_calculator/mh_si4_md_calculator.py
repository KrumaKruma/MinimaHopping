#!/usr/bin/env python3
from ase.lattice.cubic import FaceCenteredCubic
from minimahopping.minhop import Minimahopping
from ase.calculators.espresso import Espresso
from ase.calculators.kim.kim import KIM
from ase.io import read
from ase.build import bulk
from ase.io import read, write

def main():
    #a0 = 3.#(15.0*4)**(1.0/3.0)
    #initial_configuration = bulk('Si')#, a=3.6, crystalstructure='fcc')
    #print(initial_configuration)
    initial_configuration = read("input.extxyz")
    pseudopotentials = {'Si': 'Si.pbe-nl-rrkjus_psl.1.0.0.UPF'}
    input_data = {
        'system': {
            'ecutwfc': 60.,
            'ecutrho': 125.,
            'ibrav'  : 0,
            'nosym' : True},
        'disk_io': 'none',
        'electrons': {
            'electron_maxstep': 500,
            'mixing_mode': 'local-TF',
            'mixing_beta': 0.7},
        'rsim':{
            'smear3d': 2.}}

    md_calculator = KIM("SW_StillingerWeber_1985_Si__MO_405512056662_005")
    
    calculator = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(3, 3, 3), koffset=(0, 0, 0),input_data=input_data)

    initial_configuration.calc = calculator
    with Minimahopping(initial_configuration, md_calculator=md_calculator, mdmin=2,verbose_output=True, T0=4000, dt0=0.01, use_MPI=False) as mh:
        mh(totalsteps=50)


if __name__ == '__main__':
    main()
    quit()

