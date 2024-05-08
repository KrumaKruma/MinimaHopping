import numpy as np
import minimahopping.mh.lattice_operations as lat_opt
import minimahopping.biomode.biomode as split_forces 
import minimahopping.logging.logger as logging
import ase.atom
from minimahopping.mh.cell_atom import Cell_atom


def soften(atoms: ase.atom.Atom, calculator: ase.calculators.calculator.Calculator, cell_atoms: Cell_atom, nsoft: int, alpha_pos: float = 1e-3, alpha_lat: float = None):
    atoms = atoms.copy()
    atoms.calc = calculator
    eps_dd = 1e-2


    # Check if softenig steps is larger than zero
    if nsoft > 0:
        # initialization and normalization of velocities
        positions_in, normed_velocities, e_pot_in, norm_const = initialize(atoms = atoms, eps_dd = eps_dd)

        # first softening step to get the initial residum
        res_initial, curve, new_normed_velocities, = update_velocities(atoms = atoms, 
                                                                                                  positions_in = positions_in, 
                                                                                                  normed_velocities = normed_velocities, 
                                                                                                  e_pot_in = e_pot_in, 
                                                                                                  eps_dd = eps_dd, 
                                                                                                  alpha_pos = alpha_pos, 
                                                                                                  alpha_lat = alpha_lat)
        normed_velocities = new_normed_velocities

        # Soften loop
        for i in range(nsoft):
            res, curve, new_normed_velocities = update_velocities(atoms = atoms, 
                                                                                positions_in = positions_in, 
                                                                                normed_velocities = normed_velocities, 
                                                                                e_pot_in = e_pot_in, 
                                                                                eps_dd = eps_dd, 
                                                                                alpha_pos = alpha_pos, 
                                                                                alpha_lat = alpha_lat)
            # Criterion for early stopping of softening
            if res < (curve * eps_dd * 0.5):
                break
        
            # Update velocities
            normed_velocities = new_normed_velocities

        # Warning if softening has diverged and not converged
        if res_initial < res:
            warning_msg = "Softening diverged"
            logging.logger.warning(warning_msg)

        # Renormalize velocities
        velocities = normed_velocities / norm_const

        # Return output
    else:
        # if there is no softening remove momentum and torque
        velocities = elim_moment(velocities=atoms.get_velocities())
        velocities = elim_torque(velocities=velocities, positions=atoms.get_positions(), masses=atoms.get_masses())
    
    if cell_atoms is not None:
        cell_velocities = elim_moment(velocities = cell_atoms.velocities)
        cell_velocities = elim_torque(velocities = cell_velocities,positions = cell_atoms.positions,masses = cell_atoms.masses)
    else:
        cell_velocities = None 
    return velocities, cell_velocities


    

def initialize(atoms: ase.atom.Atom, eps_dd: float):
    '''
    Initialization before the iterative part of the softening
    '''

    # Get the initial positions and velocities
    positions_in = atoms.get_positions()
    e_pot_in = atoms.get_potential_energy()
    velocities = atoms.get_velocities()

    # Get the normalization constant
    normed_velocities = elim_moment(velocities = velocities)
    normed_velocities = elim_torque(velocities = velocities, positions = positions_in,masses = atoms.get_masses())
    norm_const = get_norm_constant(velocities, eps_dd)
    normed_velocities = velocities * norm_const

    return positions_in, normed_velocities, e_pot_in, norm_const


def get_norm_constant(velocities, eps_dd):
    '''
    get the normalization constant for the velocities
    '''
    norm_const = eps_dd / np.sqrt(np.sum(velocities ** 2))

    return norm_const


def update_velocities(atoms: ase.atom.Atom, 
                      positions_in: np.ndarray, 
                      normed_velocities: np.ndarray, 
                      e_pot_in: float, 
                      eps_dd: float, 
                      alpha_pos: float, 
                      alpha_lat: float):
    '''
    Performing one softening steps of the velocities
    '''
    positions = positions_in.copy() + normed_velocities
    atoms.set_positions(positions)

    e_pot = atoms.get_potential_energy()
    forces = atoms.get_forces()
    # Only a parameter for check
    fd2 = 2* (e_pot - e_pot_in) / eps_dd ** 2




    sdf = np.sum(normed_velocities * forces)
    sdd = np.sum(normed_velocities * normed_velocities)

    curve = -sdf/sdd

    tt = np.sum(forces * forces) # only for debugging reasons
    forces = forces + curve * normed_velocities
    res = np.sum(forces * forces)


    tt = np.sqrt(tt)
    res = np.sqrt(res)

    debug_msg = "SOFTEN:   {:1.5f}    {:1.5f}    {:1.5f}    {:1.5f}    {:1.5f}".format(tt, res, curve , fd2, e_pot - e_pot_in)
    logging.logger.debug(debug_msg)

    print(debug_msg)

    positions = positions + alpha_pos * forces
    normed_velocities = positions - positions_in

    normed_velocities = elim_moment(velocities = normed_velocities)
    normed_velocities = elim_torque(velocities = normed_velocities, positions = positions,masses = atoms.get_masses())
    divisor = np.sqrt(np.sum(normed_velocities ** 2))

    sdd = eps_dd / divisor

    normed_velocities = normed_velocities * sdd

    return res, curve, normed_velocities


def elim_moment(velocities: np.ndarray):
    """
    Elimination of the momentum in the velocities
    """
    # eliminiation of momentum
    _s = np.sum(velocities, axis=0) / velocities.shape[0]
    no_moment_velocities = velocities - _s
    return no_moment_velocities


def elim_torque(velocities: np.ndarray, positions: np.ndarray, masses: np.ndarray):
    """
    Elimination of the torque in the velocites
    """

    # Calculate center of mass and subtract it from positions
    total_mass = np.sum(masses)
    masses_3d = np.vstack([masses] * 3).T
    weighted_positions = positions * masses_3d

    cm = np.sum(weighted_positions, axis=0) / total_mass
    weighted_positions -= cm

    # Calculate moments of inertia with both functions
    evaleria, teneria = moment_of_inertia(positions = positions,masses = masses)

    # New vrot calculation: Vectorized operation replacing the loop
    teneria_reshaped = teneria.T.reshape(3, 1, 3)
    vrot = np.cross(teneria_reshaped, positions[None, :, :]).transpose(1, 2, 0)

    # flatten velocities and reshape vrot to match dimensions
    velocities = velocities.flatten()
    vrot = vrot.reshape((positions.shape[0] * 3, 3), order="C")

    # normalize vrot using a loop
    for i, vec in enumerate(vrot.T):
        vrot[:, i] = normalize(v = vec)
        weighted_positions += cm

    # New Implementation: Vectorized operation replacing the above for loop
    # mask for elements of evaleria that are greater than 1e-10
    mask = np.abs(evaleria) > 1e-10

    # calculate alpha and update velocities using np.einsum for dot product and updating velocities
    alpha = np.einsum('ij,i->j', vrot[:, mask], velocities)
    velocities -= np.einsum('ij,j->i', vrot[:, mask], alpha)

    # reshape velocities back to the original shape
    velocities = velocities.reshape((positions.shape[0], 3))


    return velocities

def moment_of_inertia(positions: np.ndarray, masses: np.ndarray):
    '''
    Calcualtion of the eigenvalues and eigenvectors of the inertia tensor
    '''
    inertia_tensor = np.zeros((3, 3))
    for at, mass in zip(positions, masses):
        inertia_tensor[0, 0] += mass * (at[1] ** 2 + at[2] ** 2)
        inertia_tensor[1, 1] += mass * (at[0] ** 2 + at[2] ** 2)
        inertia_tensor[2, 2] += mass * (at[0] ** 2 + at[1] ** 2)
        inertia_tensor[0, 1] -= mass * (at[0] * at[1])
        inertia_tensor[0, 2] -= mass * (at[0] * at[2])
        inertia_tensor[1, 2] -= mass * (at[1] * at[2])

    inertia_tensor[1, 0] = inertia_tensor[0, 1]
    inertia_tensor[2, 0] = inertia_tensor[0, 2]
    inertia_tensor[2, 1] = inertia_tensor[1, 2]

    eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
    return eigenvalues, eigenvectors


def normalize(v: np.ndarray):
    """
    Function that normalized a vector of arbitrary length
    """
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm



#         subroutine soften(nsoft,nat,rxyz,dd,curv0,curv,res,it,padd,pweight, &
#                           count_soft)
#         implicit real*8 (a-h,o-z)
#         dimension rxyz(3*nat),dd(3*nat)
#         real*8 , allocatable :: wpos(:),fxyz(:)
#         logical padd

#         eps_dd=1.d-2
# !        alpha=4.d-2   ! Au gold
# !        alpha=0.8d-3  ! step size for LJ systems
# !        alpha=3.d-2  ! step size for 
#          ! alpha=40.d0   !  Carbon in eV/A**2
#           alpha=0.5d0   !  Carbon in atomic units

#         allocate(wpos(3*nat),fxyz(3*nat))

#         call  energyandforces(nat,rxyz,fxyz,etot0,padd,pweight,count_soft)



# ! normalize initial guess
#         sdd=0.d0
#         do i=1,3*nat
#         sdd=sdd+dd(i)**2
#         enddo
#         sdd=eps_dd/sqrt(sdd)
#         do i=1,3*nat
#         dd(i)=sdd*dd(i)
#         enddo


#     do it=1,nsoft

#         do i=1,3*nat
#         wpos(i)=rxyz(i)+dd(i)
#         enddo
#         call  energyandforces(nat,wpos,fxyz,etot,padd,pweight,count_soft)

#         fd2=2.d0*(etot-etot0)/eps_dd**2

#         sdf=0.d0
#         sdd=0.d0
#         do i=1,3*nat
#         sdf=sdf+dd(i)*fxyz(i)
#         sdd=sdd+dd(i)*dd(i)
#         enddo

#         curv=-sdf/sdd
#         if (it.eq.1) curv0=curv
#         res=0.d0
#         tt=0.d0
#         do i=1,3*nat
#             tt=tt+fxyz(i)**2
#             fxyz(i)=fxyz(i)+curv*dd(i)
#             res=res+fxyz(i)**2
#         enddo
#         res=sqrt(res)
#         tt=sqrt(tt)
#         write(100,'(a,i5,4(e13.5),e18.10)') 'it,fnrm,res,curv,fd2,etot ',it,tt,res,curv,fd2,etot-etot0

#         do i=1,3*nat
#             wpos(i)=wpos(i)+alpha*fxyz(i)
#         enddo
#         do i=1,3*nat
#             dd(i)=wpos(i)-rxyz(i)
#         enddo
#         call  elim_moment(nat,dd)
#         call  elim_torque(nat,rxyz,dd)
#         sdd=0.d0
#         do i=1,3*nat
#         sdd=sdd+dd(i)*dd(i)
#         enddo
#         if (res.le.curv*eps_dd*5.d-1) goto 1000
#         sdd=eps_dd/sqrt(sdd)
#         do i=1,3*nat
#         dd(i)=dd(i)*sdd
#         enddo

    #   enddo





