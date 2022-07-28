import numpy as np
from ase.io import read, write
from ase.calculators.lj import LennardJones
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import Atoms
from ase.md.verlet import VelocityVerlet
from ase import units
#import torque
from ase.optimize import BFGS
from ase.optimize import QuasiNewton
#from nequip.ase import NequIPCalculator
from ase.build import bulk
from ase.calculators.espresso import Espresso



def frac2cart(atoms, reduced_positions):
    cell = atoms.get_cell()
    positions = np.zeros(reduced_positions.shape)

    for i,at in enumerate(reduced_positions):
        positions[i,:] = np.matmul(cell,at)

    return positions

def cart2frac(atoms):
    positions = atoms.get_positions()
    cell = atoms.get_cell()
    inv_cell = np.linalg.inv(cell)
    reduced_positions = np.zeros(positions.shape)
    for i,at in enumerate(positions):
        reduced_positions[i,:] = np.matmul(inv_cell, at)


    return reduced_positions

def lattice_derivative(atoms):
    stress_tensor = atoms.get_stress(voigt=False)
    cell = atoms.get_cell(complete=False)
    inv_cell = np.linalg.inv(cell)
    prefact = np.linalg.det(cell)
    deralat = prefact * np.matmul(stress_tensor, inv_cell.T)
    return deralat



def get_moment(velocities):
    # eliminiation of momentum
    s = np.sum(velocities, axis=0)/velocities.shape[0]
    print(s)
    return None




def get_torque(positions, velocities, masses):
    total_mass = np.sum(masses)
    masses_3d = np.vstack([masses] * 3).T
    weighted_positions = positions * masses_3d
    cm = np.sum(weighted_positions, axis=0)
    cm /= total_mass

    tv = 0
    for at, v in zip(positions, velocities):
        tv += np.cross(at,v)
    print(tv)

    return None



def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
       return v
    return v / norm


def moment_of_inertia(masses, positions):
    inertia_tensor = np.zeros((3,3))
    for at, mass in zip(positions, masses):
        inertia_tensor[0,0] += mass * (at[1]**2 + at[2]**2)
        inertia_tensor[1,1] += mass * (at[0]**2 + at[2]**2)
        inertia_tensor[2,2] += mass * (at[0]**2 + at[1]**2)
        inertia_tensor[0,1] -= mass * (at[0] * at[1])
        inertia_tensor[0,2] -= mass * (at[0] * at[2])
        inertia_tensor[1,2] -= mass * (at[1] * at[2])

    inertia_tensor[1,0] = inertia_tensor[0,1]
    inertia_tensor[2,0] = inertia_tensor[0,2]
    inertia_tensor[2,1] = inertia_tensor[1,2]

    eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
    return eigenvalues, eigenvectors



def elim_moment(velocities):
    # eliminiation of momentum
    s = np.sum(velocities, axis=0)/velocities.shape[0]
    velocities -= s
    return velocities

def elim_torque(positions, velocities, masses):
    #elimination of torque
    masses = np.ones((positions.shape[0]))
    #calculate center of mass and subtracti it from positions
    total_mass = np.sum(masses)
    masses_3d = np.vstack([masses]*3).T
    weighted_positions = positions * masses_3d
    #get_torque(weighted_positions, velocities, masses)
    cm = np.sum(weighted_positions, axis=0)
    cm /= total_mass
    weighted_positions -= cm


    evaleria, teneria = moment_of_inertia(masses, positions)

    vrot = np.zeros((positions.shape[0], 3, 3))
    for iat, at in enumerate(positions):
        vrot[iat, :, 0] = np.cross(teneria[:,0], at)
        vrot[iat, :, 1] = np.cross(teneria[:,1], at)
        vrot[iat, :, 2] = np.cross(teneria[:,2], at)

    velocities = velocities.flatten()
    vrot = vrot.reshape((positions.shape[0]*3,3),order="C")
    #print("DEBI:   ", vrot[9,1])
    for i, vec in enumerate(vrot.T):
        vrot[:,i] = normalize(vec)

    #print(vrot.reshape((38*3,3)))

    weighted_positions += cm

    for i, eval in enumerate(evaleria):
        if abs(eval) > 1e-10:
            alpha = np.dot(vrot[:,i], velocities)
            velocities -= alpha * vrot[:,i]

    velocities = velocities.reshape((positions.shape[0],3))

    #get_torque(weighted_positions, velocities, masses)
    return velocities



def soften(atoms, nsoft):
    # Softening constants
    eps_dd = 1e-2
    alpha = 1e-3
    # Get initial energy
    e_pot_in = atoms.get_potential_energy()
    positions_in = atoms.get_positions()
    masses = atoms.get_masses()

    # Normalize initial guess
    velocities = atoms.get_velocities()
    norm_const = eps_dd/np.sqrt(np.sum(velocities**2))
    velocities *= norm_const

    # Softening cycle
    for it in range(nsoft):
        w_positions = positions_in + velocities
        atoms.set_positions(w_positions)

        e_pot = atoms.get_potential_energy()
        forces = atoms.get_forces()

        # fd2 is only a control variable
        fd2 = 2. * (e_pot - e_pot_in)/eps_dd**2

        sdf = np.sum(velocities * forces)
        sdd = np.sum(velocities * velocities)
        curve = -sdf/sdd
        if it == 0:
            curve0 = curve

        tt = np.sqrt(np.sum(forces*forces))
        forces += curve * velocities
        res = np.sqrt(np.sum(forces*forces))

        print(it, tt, res, curve, fd2, e_pot - e_pot_in)
        #return velocities
        #quit()
        w_positions = w_positions + alpha * forces
        velocities = w_positions-positions_in


        velocities = elim_moment(velocities)

        #get_torque(w_positions, velocities, masses)
        velocities = elim_torque(w_positions, velocities, masses)
        #velocities = torque.elim_torque_reza(w_positions.flatten(), velocities.flatten())
        velocities = velocities.reshape(w_positions.shape[0],3)
        #get_moment(velocities)
        #get_torque(w_positions, velocities, masses)

        sdd = eps_dd/np.sqrt(np.sum(velocities**2))
        if res < (curve*eps_dd*0.5):
            break
        velocities *= sdd
    velocities /= norm_const
    # Restore initial positions in atoms object
    atoms.set_positions(positions_in)

    return velocities


def md(atoms,dt):
    # MD which visits at least three max
    masses = atoms.get_masses()[:, np.newaxis] / atoms.get_masses()[:, np.newaxis]  # for the moment no masses
    epot_old = atoms.get_potential_energy()
    sign_old = -1
    n_change = 0
    i = 0
    while n_change < 4:
        forces = atoms.get_forces()
        velocities = atoms.get_velocities()
        positions = atoms.get_positions()

        #Update postions
        atoms.set_positions(positions + dt*velocities + 0.5*dt*dt*(forces/(masses)))

        new_forces = atoms.get_forces()
        atoms.set_velocities(velocities + 0.5 * dt * ((forces + new_forces)/masses))

        epot = atoms.get_potential_energy()
        #print("MD", str(epot))
        sign = int(np.sign(epot_old-epot))
        if sign_old != sign:
            sign_old = sign
            n_change += 1
        #print(n_change, epot_old-epot)

        epot_old = epot
        #filename = "MD" + str(i).zfill(4) + ".xyz"
        #write(filename, atoms)
        i += 1

    return None

def vcsmd(atoms,cell_atoms,dt):

    # MD which visits at least three max
    # Initializations
    masses = atoms.get_masses()[:, np.newaxis] / atoms.get_masses()[:, np.newaxis]  # for the moment no masses
    cell_masses = cell_atoms.get_masses()[:, np.newaxis]  # for the moment no masses
    epot_old = atoms.get_potential_energy()
    sign_old = -1
    n_change = 0
    i = 0

    forces = atoms.get_forces()
    cell_forces = -lattice_derivative(atoms)
    #while n_change < 4:
    for j in range(300):
        velocities = atoms.get_velocities()
        positions = atoms.get_positions()

        # Update postions
        atoms.set_positions(positions + dt*velocities + 0.5*dt*dt*(forces/(masses)))

        # Update lattice so that fractional coordinates remain invariant
        reduced_postitions = cart2frac(atoms)
        cell_positions = cell_atoms.get_positions()
        cell_velocities = cell_atoms.get_velocities()
        cell_atoms.set_positions(cell_positions + dt*cell_velocities + 0.5*dt*dt*(cell_forces/cell_masses))
        atoms.set_cell(cell_atoms.get_positions())
        positions = frac2cart(atoms,reduced_postitions)
        atoms.set_positions(positions)

        # Update velocities
        new_forces = atoms.get_forces()
        atoms.set_velocities(velocities + 0.5 * dt * ((forces + new_forces) / masses))

        # Update velocities of the cell atoms
        new_cell_forces = -lattice_derivative(atoms)
        cell_atoms.set_velocities(cell_velocities + 0.5 * dt * ((cell_forces + new_cell_forces) / cell_masses))

        forces = new_forces
        cell_forces = new_cell_forces

        epot = atoms.get_potential_energy()
        #print("MD", str(epot))
        sign = int(np.sign(epot_old-epot))
        if sign_old != sign:
            sign_old = sign
            n_change += 1
        #print(n_change, epot_old-epot)

        epot_old = epot
        e_kin = 0.5*np.sum(masses*atoms.get_velocities() * atoms.get_velocities())
        e_kin += 0.5 * np.sum(cell_masses * cell_atoms.get_velocities() * cell_atoms.get_velocities())
        print("STEP:   ", i,epot, e_kin, epot + e_kin)
        filename = "MD" + str(i).zfill(4) + ".ascii"
        write(filename, atoms)
        i += 1

    return None









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

        # Soften part for now commented...
        #print(np.sqrt(np.sum(atoms.get_momenta()*atoms.get_momenta())))
        #velocities = soften(atoms, 10)
        #atoms.set_velocities(velocities)
        #print(np.sqrt(np.sum(atoms.get_momenta() * atoms.get_momenta())))
        if True in atoms.pbc:
            cell_atoms = Atoms('Ar3', positions=atoms.get_cell())
            cell_atoms.set_masses([10,10,10])
            MaxwellBoltzmannDistribution(cell_atoms, temperature_K=20)
            #velocities = soften(atoms,10)
            #velocities = vcs_soften(atoms, cell_atoms, 10)
            vcsmd(atoms,cell_atoms,  dt)
            #md(atoms,dt)
        else:
            velocities = soften(atoms, 10)
            atoms.set_velocities(velocities)
            md(atoms,dt)
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
    #path = "/home/marco/NequipMH/data/38/"
    #filename = "acc-00000.xyz"
    path = "/home/marco/NequipMH/data/"
    filename = "MIN_IN.xyz"
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
    T = 300
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
    calculator.parameters.rc = 6.0

    # atoms = bulk('NaCl', crystalstructure='rocksalt', a=6.0)
    # pseudo_dir = '/home/marco/NequipMH/src/python/'
    # pseudopotentials = {'Na': 'Na.UPF',
    #                     'Cl': 'Cl.UPF'}
    #
    # input_data = {
    #     'system': {
    #         'ecutwfc': 100,
    #         'ecutrho': 500},
    #     'disk_io': 'low'}  # automatically put into 'control'
    #
    # calc = Espresso(pseudo_dir=pseudo_dir, pseudopotentials=pseudopotentials,
    #                 tstress=True, tprnfor=True, kpts=(1, 1, 1), input_data=input_data)

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



