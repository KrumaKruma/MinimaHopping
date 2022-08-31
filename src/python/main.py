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
from ase.constraints import UnitCellFilter
from ase import units
import reshapecell
import time
import bazant
import periodic_sqnm
import warnings



class cell_atom:
    def __init__(self, positions, mass=None, velocities=None):
        self.positions = positions
        self.masses = np.array([mass,mass,mass])
        self.velocities = velocities

    def set_velocities_boltzmann(self, temperature):
        xi = np.random.standard_normal((len(self.masses), 3))
        temp = units.kB * temperature
        self.velocities = xi * np.sqrt(self.masses * temp)[:,np.newaxis]
        self.velocities /= self.masses






def get_area(lattice):
    area = 0.
    area += np.linalg.norm(np.cross(lattice[:, 0], lattice[:, 1]))
    area += np.linalg.norm(np.cross(lattice[:, 0], lattice[:, 2]))
    area += np.linalg.norm(np.cross(lattice[:, 1], lattice[:, 2]))

    return area



def reshape_cell(atoms, imax):
    lattice_in = atoms.get_cell()

    lattice = np.zeros((3,3))
    lattice_min = np.zeros((3,3))
    area_min = 1e100

    volium_in = abs(np.linalg.det(lattice_in))

    for i13 in range(-imax, imax, 1):
        for i12 in range(-imax, imax, 1):
            lattice[0,:] = lattice_in[0,:] + i12 * lattice_in[1,:] + i13 * lattice_in[2,:]
            for i21 in range(-imax, imax, 1):
                for i23 in range(-imax, imax, 1):
                    lattice[1,:] = lattice_in[1,:] + i21 * lattice_in[0,:] + i23 * lattice_in[2,:]
                    for i31 in range(-imax, imax, 1):
                        for i32 in range(-imax, imax, 1):
                            lattice[2,:] = lattice_in[2,:] + i31 * lattice_in[0,:] + i32 * lattice_in[1,:]
                            volium = abs(np.linalg.det(lattice))
                            if abs(volium_in-volium) < 1e-8:
                                area = get_area(lattice)
                                if area < area_min:
                                    area_min = area
                                    lattice_min[:,:] = lattice[:,:]

    atoms.set_cell(lattice_min,scale_atoms=False, apply_constraint=False)





def energyandforces(atoms):
    positions = atoms.get_positions()
    cell = atoms.get_cell(complete=False).T
    nat = positions.shape[0]
    e_pot, forces, deralat, stress_tensor = bazant.energyandforces_bazant(cell,positions.T,nat)

    return e_pot, forces.T, deralat, stress_tensor






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
    #BAZANT TEST
    #______________________________________________________________________________________
    e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
    #stress_tensor = atoms.get_stress(voigt=False,apply_constraint=False)
    #______________________________________________________________________________________
    cell = atoms.get_cell(complete=False)

    inv_cell = np.linalg.inv(cell)
    prefact = np.linalg.det(cell)
    deralat = - prefact * np.matmul(stress_tensor, inv_cell)
    return deralat.T



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
    #masses = np.ones((positions.shape[0]))
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


def vcs_soften(atoms, cell_atoms, nsoft):
    # Softening constants
    eps_dd = 1e-2
    alpha_pos = 1e-3
    alpha_lat = 1e-3
    cell_masses = cell_atoms.masses  # for the moment no masses
    masses = atoms.get_masses()
    # Get initial energy
    #BAZANT TEST
    #___________________________________________________________________________________________________________
    # e_pot_in   = atoms.get_potential_energy()
    e_pot_in , forces, deralat, stress_tensor = energyandforces(atoms)
    #___________________________________________________________________________________________________________
    positions_in = atoms.get_positions()
    cell_positions_in = cell_atoms.positions


    # Normalize initial guess
    velocities = atoms.get_velocities()
    cell_velocities = cell_atoms.velocities
    norm_const = eps_dd/(np.sqrt(np.sum(velocities**2) + np.sum(cell_velocities**2)))
    velocities *= norm_const
    cell_velocities *= norm_const


    # Softening cycle
    for it in range(nsoft):
        w_positions = positions_in + velocities
        w_cell_positions = cell_positions_in + cell_velocities

        atoms.set_positions(w_positions)
        reduced_positions = cart2frac(atoms)
        atoms.set_cell(w_cell_positions, scale_atoms=False, apply_constraint=False)
        positions = frac2cart(atoms, reduced_positions)
        atoms.set_positions(positions)

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        # forces = atoms.get_forces()
        # e_pot = atoms.get_potential_energy()
        e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________



        deralat = lattice_derivative(atoms)

        # fd2 is only a control variable
        fd2 = 2. * (e_pot-e_pot_in)/(eps_dd**2)

        sdf = np.sum(velocities * forces) + np.sum(cell_velocities * deralat)
        sdd = np.sum(velocities * velocities) + np.sum(cell_velocities * cell_velocities)

        curve = -sdf/sdd
        if it == 0:
            curve0 = curve

        tt = np.sqrt(np.sum(forces*forces) + np.sum(deralat*deralat))
        forces += curve * velocities
        deralat += curve * cell_velocities
        res = np.sqrt(np.sum(forces*forces) + np.sum(deralat*deralat))


        #print("SOFTEN:   ", it, tt, res, curve/fd2, e_pot - e_pot_in)

        w_positions = w_positions + alpha_pos * forces
        w_cell_positions = w_cell_positions + alpha_lat * deralat
        velocities = w_positions-positions_in
        cell_velocities = w_cell_positions - cell_positions_in

        velocities = elim_moment(velocities)
        cell_velocities = elim_torque(w_cell_positions, cell_velocities, cell_masses)

        sdd = eps_dd/np.sqrt(np.sum(velocities**2)+np.sum(cell_velocities**2))
        if res < (curve*eps_dd*0.5):
            break
        velocities *= sdd
        cell_velocities *= sdd


    velocities /= norm_const
    cell_velocities /= norm_const
    # Restore initial positions in atoms object
    atoms.set_positions(positions_in)
    atoms.set_cell(cell_positions_in)
    cell_atoms.positions = cell_positions_in

    return velocities, cell_velocities

def soften(atoms, nsoft):
    # Softening constants
    eps_dd = 1e-2
    alpha = 1e-3
    # Get initial energy
    # BAZANT TEST
    # ___________________________________________________________________________________________________________
    # e_pot_in = atoms.get_potential_energy()
    e_pot_in, forces, deralat, stress_tensor = energyandforces(atoms)
    # ___________________________________________________________________________________________________________
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

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        #e_pot = atoms.get_potential_energy()
        #forces = atoms.get_forces()
        e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________


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

        #print(it, tt, res, curve/ fd2,e_pot - e_pot_in)
        #return velocities
        #quit()
        w_positions = w_positions + alpha * forces
        velocities = w_positions-positions_in


        velocities = elim_moment(velocities)

        #get_torque(w_positions, velocities, masses)
        velocities = elim_torque(w_positions, velocities, masses)
        #velocities = torque.elim_torque_reza(w_positions.flatten(), velocities.flatten())
        #velocities = velocities.reshape(w_positions.shape[0],3)
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

    #BAZANT TEST
    #___________________________________________________________________________________________________________
    # epot_old  = atoms.get_potential_energy()
    # forces = atoms.get_forces()
    epot_old, forces, deralat, stress_tensor = energyandforces(atoms)
    #___________________________________________________________________________________________________________

    sign_old = -1
    n_change = 0
    i = 0
    #while n_change < 4:
    for i in range(100000):
        velocities = atoms.get_velocities()
        positions = atoms.get_positions()

        #Update postions
        atoms.set_positions(positions + dt*velocities + 0.5*dt*dt*(forces/(masses)))

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        # epot_old  = atoms.get_potential_energy()
        # new_forces = atoms.get_forces()
        epot_old, new_forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________
        atoms.set_velocities(velocities + 0.5 * dt * ((forces + new_forces)/masses))

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        # epot  = atoms.get_potential_energy()
        epot, forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________

        #print("MD", str(epot))
        sign = int(np.sign(epot_old-epot))
        if sign_old != sign:
            sign_old = sign
            n_change += 1
        #print(n_change, epot_old-epot)

        epot_old = epot
        e_kin = 0.5*np.sum(masses*atoms.get_velocities() * atoms.get_velocities())
        print("STEP:   ", i,epot, e_kin, epot + e_kin)
        #filename = "MD" + str(i).zfill(4) + ".xyz"
        #write(filename, atoms)
        i += 1

    return None

def vcsmd(atoms,cell_atoms,dt):

    # MD which visits at least three max
    # Initializations
    masses = atoms.get_masses()[:, np.newaxis] / atoms.get_masses()[:, np.newaxis]  # for the moment no masses
    cell_masses = cell_atoms.masses[:, np.newaxis]

    #BAZANT TEST
    #___________________________________________________________________________________________________________
    #epot_old  = atoms.get_potential_energy()
    epot_old, forces, deralat, stress_tensor = energyandforces(atoms)
    #___________________________________________________________________________________________________________

    sign_old = -1
    n_change = 0
    i = 0

    #BAZANT TEST
    #___________________________________________________________________________________________________________
    #forces = atoms.get_forces()
    e_pot_cur, forces, deralat, stress_tensor = energyandforces(atoms)
    #___________________________________________________________________________________________________________
    cell_forces = lattice_derivative(atoms)
    while n_change < 4:
        e_kin = 0.5*np.sum(masses*atoms.get_velocities() * atoms.get_velocities())
        e_kin += 0.5 * np.sum(cell_masses * cell_atoms.velocities * cell_atoms.velocities)
        print("MD   STEP:   ", i,epot_old, e_kin, epot_old + e_kin)

        velocities = atoms.get_velocities()
        positions = atoms.get_positions()

        # Update postions
        atoms.set_positions(positions + dt*velocities + 0.5*dt*dt*(forces/(masses)))

        # Update lattice so that fractional coordinates remain invariant
        reduced_postitions = cart2frac(atoms)
        cell_positions = cell_atoms.positions
        cell_velocities = cell_atoms.velocities
        cell_atoms.positions = cell_positions + dt*cell_velocities + 0.5*dt*dt*(cell_forces/cell_masses)
        atoms.set_cell(cell_atoms.positions,scale_atoms=False, apply_constraint=False)
        positions = frac2cart(atoms,reduced_postitions)
        atoms.set_positions(positions)

        # Update velocities
        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        # new_forces = atoms.get_forces()
        e_pot_cur, new_forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________

        atoms.set_velocities(velocities + 0.5 * dt * ((forces + new_forces) / masses))

        # Update velocities of the cell atoms
        new_cell_forces = lattice_derivative(atoms)
        cell_atoms.velocities = cell_velocities + 0.5 * dt * ((cell_forces + new_cell_forces) / cell_masses)

        forces = new_forces
        cell_forces = new_cell_forces

        # Update velocities
        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        # epot = atoms.get_potential_energy()
        epot, new_forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________

        #print("MD", str(epot))
        sign = int(np.sign(epot_old-epot))
        if sign_old != sign:
            sign_old = sign
            n_change += 1
        #print(n_change, epot_old-epot)

        epot_old = epot
        #e_kin = 0.5*np.sum(masses*atoms.get_velocities() * atoms.get_velocities())
        #e_kin += 0.5 * np.sum(cell_masses * cell_atoms.get_velocities() * cell_atoms.get_velocities())
        #print("STEP:   ", i,epot, e_kin, epot + e_kin)
        filename = "MD" + str(i).zfill(4) + ".ascii"
        write(filename, atoms)
        i += 1
        #Time4 = time.time()
        #Tot_Time = Time4-Time1
        #DFT_Time = Time3-Time2
        #Rest_Time = Tot_Time - DFT_Time
        #print("Total Time:   ", Tot_Time, "Time DFT:   ", DFT_Time, "REST Time:  ", Rest_Time)

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
    #BAZANT TEST
    #___________________________________________________________________________________________________________
    #e_pot_curr = atoms.get_potential_energy()
    e_pot_curr, forces, deralat, stress_tensor = energyandforces(atoms)
    #___________________________________________________________________________________________________________
    escape = 0.0
    beta_s = 1.001
    while escape < 1e-3:

        MaxwellBoltzmannDistribution(atoms, temperature_K=T)

        if True in atoms.pbc:
            mass = .75 * np.sum(atoms.get_masses()) / 10.
            cell_atoms = cell_atom(mass = mass, positions = atoms.get_cell())
            cell_atoms.set_velocities_boltzmann(temperature=T)
            velocities, cell_velocities = vcs_soften(atoms, cell_atoms, 20)
            atoms.set_velocities(velocities)
            cell_atoms.velocities = cell_velocities
            vcsmd(atoms,cell_atoms,  dt)
            reshape_cell(atoms,3)
            reshape_cell(atoms,3)
            quit()
            #positions, cell_positions = reshapecell.reshapecell(atoms.get_cell().T/0.52917721067, atoms.get_positions().T/0.52917721067)
            #atoms.set_cell(cell_positions.T*0.52917721067)
            #atoms.set_positions(positions.T*0.52917721067)
        else:
            velocities = soften(atoms, 20)
            atoms.set_velocities(velocities)
            md(atoms,dt)
        #optimizer(atoms, 0.005)
        ucf = UnitCellFilter(atoms)
        dyn = QuasiNewton(ucf)
        dyn.run(fmax = 1e-2)
        e_pot = atoms.get_potential_energy()
        escape = abs(e_pot-e_pot_curr)
        T *= beta_s

    return atoms

def in_history(e_pot_cur,history):
    i = 0
    for e_pot, n_v in history:
        e_diff = abs(e_pot - e_pot_cur)
        if e_diff < 1e-2:
            i += 1
    return i+1



def main():

    # Initial settings and paths
    #path = "/home/marco/NequipMH/data/38/"
    #filename = "acc-00000.xyz"
    #path = "/home/marco/NequipMH/data/"
    #filename = "LJC.extxyz"
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
    T = 3000
    beta_decrease = 1./1.1
    beta_increase = 1.1

    # Read local minimum input file
    #atoms = read(filename)

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
    #calc = LennardJones()
    #calc.parameters.epsilon = 1.0
    #calc.parameters.sigma = 1.0
    #calc.parameters.rc = 12.0

    #atoms = bulk('NaCl', crystalstructure='rocksalt', a=6.0)
    atoms = read("Si_in2.extxyz")
    #atoms.pbc = True
    nat = atoms.get_positions().shape[0]
    init_lat = atoms.get_cell().T
    inital_step_size = 0.01
    nhist_max=10
    lattice_weight = 2
    alpha_min = 0.001
    eps_subsop = 1e-4

    optim = periodic_sqnm.periodic_sqnm(nat, init_lat, inital_step_size,nhist_max,lattice_weight,alpha_min,eps_subsop)

    max_force_comp = 100
    max_force_threshold = 0.002
    i = 0
    while max_force_comp > max_force_threshold:
        #BAZANT TEST
        #___________________________________________________________________________________________________________
        #energy = atoms.get_potential_energy() and also deralat
        energy, forces, deralat, stress_tensor = energyandforces(atoms)
        #___________________________________________________________________________________________________________

        print(i, energy, np.max(forces))#/51.42208619083232)
        i += 1
        max_force_comp = np.max(forces)
        positions = atoms.get_positions()
        lattice = atoms.get_cell().T

        new_positions, new_lattice = optim.optimizer_step(positions.T, lattice, energy, forces.T, deralat)
        atoms.set_positions(new_positions.T)
        atoms.set_cell(new_lattice.T, scale_atoms=False, apply_constraint=False)

        if i > 10:
            warning_msg = "Geometry did not converge in " + str(i) + " optimizations steps"
            warnings.warn(warning_msg, FutureWarning)
            break


    quit()
    reshape_cell(atoms,3)
    reshape_cell(atoms, 3)


    #pseudo_dir = '/home/marco/NequipMH/src/python/'
    #pseudopotentials = {'Si': 'Si.UPF'}

    #input_data = {
    #    'system': {
    #        'ecutwfc': 44,
    #        'ecutrho': 175},
    #    'disk_io': 'low',
    #    'electrons': {
    #        'mixing_beta' : 0.4
    #    }
    #}  # automatically put into 'control'

    #calc = Espresso(pseudo_dir=pseudo_dir, pseudopotentials=pseudopotentials,
    #                tstress=True, tprnfor=True, kpts=(3, 3, 3), input_data=input_data)

    #atoms.calc = calc
    #ucf = UnitCellFilter(atoms)
    #dyn = QuasiNewton(ucf)
    #dyn.run(fmax=1e-2)
    #write("LJC.extxyz", atoms)
    pos_cur = atoms.get_positions()

    #BAZANT TEST
    #___________________________________________________________________________________________________________
    #e_pot_cur = atoms.get_potential_energy()
    e_pot_cur, forces, deralat, stress_tensor = energyandforces(atoms)
    #___________________________________________________________________________________________________________

    history.append((e_pot_cur, 1))
    for i in range(10000):
        atoms = escape_trial(atoms, dt, T)

        filename = "MIN" + str(n_min).zfill(5) + ".ascii"
        write(filename, atoms)

        e_pot = atoms.get_potential_energy()
        n_min += 1
        print("LOG:",n_min,  " Epot:  ",e_pot," Ediff:   ", e_diff," T:  ", T)
        if abs(e_pot_cur-e_pot) < e_diff:
            e_pot_cur = e_pot
            e_diff *= alpha_a
            n_acc += 1
            filename = "ACC" + str(n_acc).zfill(5) + ".ascii"
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



