import numpy as np
from ase.io import read, write
from ase.calculators.lj import LennardJones
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

from ase.optimize import BFGS
from ase.optimize import QuasiNewton
#from nequip.ase import NequIPCalculator



def md(atoms,dt):
    # MD which visits at least three max
    masses = atoms.get_masses()[:, np.newaxis] / atoms.get_masses()[:, np.newaxis]  # for the moment no masses
    epot_old = atoms.get_potential_energy()
    sign_old = -1
    n_change = 0
    i = 0
    while n_change < 4:
        forces = atoms.get_forces()
        velocities = atoms.get_momenta()
        positions = atoms.get_positions()

        #Update postions
        atoms.set_positions(positions + dt*velocities + 0.5*dt*dt*(forces/(masses)))

        new_forces = atoms.get_forces()
        atoms.set_momenta(velocities + 0.5*dt * ((forces + new_forces)/masses), apply_constraint=False)

        epot = atoms.get_potential_energy()
        sign = int(np.sign(epot_old-epot))
        if sign_old != sign:
            sign_old = sign
            n_change += 1
        #print(n_change, epot_old-epot)

        epot_old = epot
        #filename = "MD" + str(i).zfill(4) + ".xyz"
        #write(filename, atoms)
        i += 1

def optimizer(atoms, criterion):
    norm = 1.
    alpha = 0.001
    forces = atoms.get_forces()
    i = 0
    while norm > criterion:
        prev_forces = forces
        positions = atoms.get_positions()
        atoms.set_positions(positions + alpha*forces)
        forces = atoms.get_forces()
        cosangle = np.sum(forces*prev_forces)/(np.linalg.norm(forces)*np.linalg.norm(prev_forces))
        if cosangle < 0.5:
            alpha /= 2.
        else:
            alpha *= 1.05

        norm = np.linalg.norm(forces)
        #print(atoms.get_potential_energy(), norm)
        #filename = "OPTIM" + str(i).zfill(4) + ".xyz"
        #write(filename,atoms)
        i += 1


def escape_trial(atoms, dt, T):
    e_pot_curr = atoms.get_potential_energy()
    escape = 0.0
    beta_s = 1.01
    while escape < 1e-6:
        MaxwellBoltzmannDistribution(atoms, temperature_K=T)
        md(atoms, dt)
        #optimizer(atoms, 0.005)
        dyn = QuasiNewton(atoms)
        dyn.run(fmax = 1e-6)
        e_pot = atoms.get_potential_energy()
        escape = abs(e_pot-e_pot_curr)
        T *= beta_s

    return atoms

def in_history(e_pot_cur,history):
    i = 0
    for e_pot, n_v in history:
        e_diff = abs(e_pot - e_pot_cur)
        if e_diff < 1e-6:
            i += 1
    return i+1



def main():

    # Initial settings and paths
    path = "/home/marco/NequipMH/data/38/"
    filename = "acc-00010.xyz"

    #path = "/home/marco/NequipMH/data/"
    #filename = "LJ38.xyz"

    #model_path = " "
    dt = 0.01
    e_diff = .01
    alpha_a = 1./1.2
    alpha_r = 1.2
    n_acc = 0
    n_min = 0
    history = []
    T = 500
    beta_decrease = 1./1.1
    beta_increase = 1.1

    # Read local minimum input file
    atoms = read(path+filename)


    # Initialize calculator (Nequip)
    #atoms.calc = NequIPCalculator.from_deployed_model()

    #    model_path=model_path,
    #    species_to_type_name = {
    #        "C" : "NequIPTypeNameForCarbon",
    #        "O" : "NequIPTypeNameForOxygen",
    #        "H" : "NequIPTypeNameForHydrogen",
    #    }
    #)

    # TESTING: For now it will be a LJ potential
    calculator = LennardJones()
    calculator.parameters.epsilon = 1.0
    calculator.parameters.sigma = 1.0
    calculator.parameters.rc = 12.0

    atoms.calc = calculator

    pos_cur = atoms.get_positions()
    e_pot_cur = atoms.get_potential_energy()
    history.append((e_pot_cur, 1))
    for i in range(10000):
        atoms = escape_trial(atoms, dt, T)

        filename = "MIN" + str(n_min).zfill(5) + ".xyz"
        write(filename, atoms)

        e_pot = atoms.get_potential_energy()
        n_min += 1
        print("LOG:",n_min,  " Epot:  ",e_pot," Ediff:   ", e_diff," T:  ", T)
        if abs(e_pot_cur-e_pot) < e_diff:
            e_pot_cur = e_pot
            e_diff *= alpha_a
            n_acc += 1
            filename = "ACC" + str(n_acc).zfill(5) + ".xyz"
            write(filename, atoms)
            pos_cur = atoms.get_positions()
        else:
            e_diff *= alpha_r
            atoms.set_positions(pos_cur)
            print(atoms.get_potential_energy()-e_pot_cur)


        n_visits = in_history(e_pot, history)
        if n_visits > 1:
            T *= beta_increase * (1. + 4. * np.log(float(n_visits)))
        else:
            T *= beta_decrease

        history.append((e_pot, n_visits))

        if i%1 == 0:
            f = open("history.dat", "w")
            for e, n in history:
                f.write(str(e) + "   " + str(n)+"\n")
            f.close()





















if __name__ == '__main__':
    main()



