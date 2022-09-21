from flare.ase.calculator import FLARE_Calculator
from flare.ase.atoms import FLARE_Atoms
from ase.io import read, write
from flare.gp import GaussianProcess
from flare.struc import Structure
from flare.utils.parameter_helper import ParameterHelper
from ase.calculators.espresso import Espresso
import numpy as np
from copy import deepcopy
from os.path import exists




class mh_otf():
    def __init__(self, gp_model, dft_calc, std_max=0.01, read=False, n_save = 10, n_train = 10):
        self.gp_model = gp_model
        self.dft_calc = dft_calc
        self.std_max = std_max
        self.flare_atoms = None
        self.read = read
        self.n_save = n_save
        self.n_train = n_train
        self.flare_calculator = FLARE_Calculator(self.gp_model,
                                            par = True,
                                            mgp_model = None,
                                            use_mapping = False)
        self.n_dft = 0
        self.dft_list = []

    def energyandforces(self, atoms):
        self.flare_atoms = FLARE_Atoms.from_ase_atoms(atoms)
        if self.read:
            self.flare_atoms.set_calculator(self.flare_calculator)
            self._read()
        if not self.gp_model.training_labels:
            energy, forces, stress = self.dft_update()
            is_dft = True
        else:
            energy, forces, stress, is_dft = self.otf_step()

        return energy, forces, stress, is_dft



    def _read(self):

        assert exists("otf.json"), "otf.json file does not exist. No restart possible. Set restart False!"
        assert exists("otf.extxyz"), "otf.extxyz file does not exist. No restart possible. Set restart False!"

        self.flare_atoms.calc.from_file("otf.json")
        structures_dft = read("otf.extxyz", index=":")
        for structure_dft in structures_dft:
            positions = structure_dft.get_positions()
            cell = structure_dft.get_cell()
            species = structure_dft.get_atomic_numbers()
            energy = structure_dft.get_potential_energy()
            forces = structure_dft.get_forces()
            dft_stress = structure_dft.get_stress()
            if dft_stress is not None:
                flare_stress = -np.array(
                    [
                        dft_stress[0],
                        dft_stress[5],
                        dft_stress[4],
                        dft_stress[1],
                        dft_stress[3],
                        dft_stress[2],
                    ]
                )

            structure = Structure(cell, species, positions)

            self.gp_model.update_db(struc=structure, forces=forces, energy=energy, stress=flare_stress)

        # print(len(self.gp_model.training_labels))



    def _write(self):
        for structure_dft in self.dft_list:
            write("otf.extxyz", structure_dft, append=True)
            self.flare_atoms.set_calculator(self.flare_calculator)
            self.flare_atoms.calc.write_model("otf")
            self.dft_list = []


    def otf_step(self):
        self.flare_atoms.set_calculator(self.flare_calculator)
        flare_energy = self.flare_atoms.get_potential_energy()
        flare_stds = self.flare_atoms.stds
        if np.max(flare_stds) > self.std_max:
            is_dft = True
            energy, forces, stress = self.dft_update()
            #print("DFT STEP:  ", energy)
            if self.n_dft%self.n_save == 0:
                self._write()
            self.n_dft += 1
        else:
            is_dft = False
            energy = self.flare_atoms.get_potential_energy()
            forces = self.flare_atoms.get_forces()
            stress = self.flare_atoms.get_stress()

            #print("FLARE STEP:   ", energy)

        return energy, forces, stress, is_dft


    def dft_update(self, ):
        self.flare_atoms.set_calculator(self.dft_calc)
        positions = self.flare_atoms.wrapped_positions
        cell = self.flare_atoms.get_cell()
        species = self.flare_atoms.get_atomic_numbers()
        energy = self.flare_atoms.get_potential_energy()
        forces = self.flare_atoms.get_forces()
        dft_stress = self.flare_atoms.get_stress()
        # self.dft_list.append([deepcopy(self.flare_atoms), energy, forces, dft_stress])
        self.dft_list.append(deepcopy(self.flare_atoms))
        # update_index = list(np.unique(np.where(self.flare_atoms.stds > 0.1)[0]))
        if dft_stress is not None:
            flare_stress = -np.array(
                [
                    dft_stress[0],
                    dft_stress[5],
                    dft_stress[4],
                    dft_stress[1],
                    dft_stress[3],
                    dft_stress[2],
                ]
            )

        structure = Structure(cell, species, positions)
        self.gp_model.update_db(struc=structure, forces=forces, energy=energy, stress=flare_stress)

        self.gp_model.set_L_alpha()
        #print(self.gp_model.likelihood)
        if self.n_dft%self.n_train == 0:
            self.gp_model.train(print_progress=False)
        #print(self.gp_model.likelihood)
        return energy, forces, dft_stress



def test():
    atoms = read("Si_in3.extxyz")
    flare_atoms = FLARE_Atoms.from_ase_atoms(atoms)
    # set up GP hyperparameters
    kernels = ['twobody', 'threebody'] # use 2+3 body kernel
    parameters = {'cutoff_twobody': 5.0,
                  'cutoff_threebody': 3.5}
    pm = ParameterHelper(
        kernels = kernels,
        random = True,
        parameters=parameters
    )

    hm = pm.as_dict()
    hyps = hm['hyps']
    cut = hm['cutoffs']
    print('hyps', hyps)

    gp_model = GaussianProcess(
        kernels = kernels,
        component = 'mc', # If you are using ASE, please set to "mc" no matter for single-component or multi-component
        hyps = hyps,
        cutoffs = cut,
        hyp_labels = ['sig2','ls2','sig3','ls3','noise'],
        opt_algorithm = 'L-BFGS-B',
        n_cpus = 1
    )


    pseudo_dir = '/home/marco/NequipMH/src/python/'
    pseudopotentials = {'Si': 'Si.UPF'}

    input_data = {
       'system': {
           'ecutwfc': 44,
           'ecutrho': 175},
       'disk_io': 'low',
       'electrons': {
           'mixing_beta' : 0.4
       }
    }  # automatically put into 'control'

    dft_calc = Espresso(pseudo_dir=pseudo_dir, pseudopotentials=pseudopotentials,
                   tstress=True, tprnfor=True, kpts=(3, 3, 3), input_data=input_data)




    OFT = mh_otf(gp_model, dft_calc)


    for i in range(2):

        energy, forces, stress, is_dft = OFT.energyandforces(atoms)
        x = atoms.get_positions()
        atoms.set_positions(np.random.normal(x, 0.005))
        print(is_dft)

if __name__=="__main__":
    test()

