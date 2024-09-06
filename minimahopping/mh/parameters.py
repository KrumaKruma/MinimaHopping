from dataclasses import dataclass, field
from dataclasses_json import dataclass_json
import logging

@dataclass_json
@dataclass
class minimaHoppingParameters:
    T0: float = 1000
    """T0 (float): Initial temperature to start MH"""
    _T: float = -1.0
    """"_T (float): Temperature for restart (if -1 then last temperature is read)"""
    Ediff0: float = .1
    """"Ediff0 (float): Initial energy difference for accepting minima"""
    _eDiff: float = -1.0
    """"_eDiff, Private (float): Energy difference for accepting minima after restarts."""
    beta_decrease: float = 1./1.05
    """beta_decrease (float): Factor for decreasing the temperature"""
    beta_increase: float = 1.05
    """beta_increase (float): Factor for increasing the temperature"""
    alpha_accept: float = 1/1.05
    """alpha_accept (float): Factor for decreasing Ediff"""
    alpha_reject: float = 1.05
    """alpha_reject (float): Factor for increasing Ediff"""
    fixed_cell_simulation: bool = False
    """fixed_cell_simulation(bool): If True the simulation cell is fixed and no variable cell shape MD and optimization are performed"""
    n_soft: int = 20
    """n_soft (int): Number of softening steps"""
    soften_positions: float = 1e-3
    """soften_positions (float): Step size for softening the positions"""
    soften_lattice: float = 1e-3
    """soften_lattice (float): Step size for softening the lattice (only relevant for periodic boundary conditions)"""
    n_S_orbitals: int = 1
    """n_S_orbitals (int): Number of s orbitals for constructing the OMFP"""
    n_P_orbitals: int = 0
    """n_P_orbitals (int): Number of p orbitals for constructing the OMFP"""
    width_cutoff: int = 4
    """width_cutoff (float): Cutoff for the OMFP"""
    exclude: list = field(default_factory=list)
    """exclude (list): List of elements to exclude in the OMFP"""
    dt0: float = 0.08
    """dt0 (float): Initial time step for the MD"""
    _dt: float = -1.0
    """_dt, Private (float): Time step of the MD for restart."""
    dt_min: float = 0.0001
    """dt_min, (float): Minimal time step of the MD."""
    mdmin: int = 2
    """mdmin (int): Number of minima visited before stopping the MD"""
    md_max_steps: int = 10000
    """md_max_steps (int): Maximum number of MD steps. If reached MH continuous"""
    collect_md_data: bool = False
    """collect_md_data (bool): Flag to collect MD data which later could be used e.g. in machine learning"""
    margin: float = 0.3
    """margin (float): Margin for fixing of fragmentation, if particle closer than (2*max_rcov + 2*margin*max_rcov)"""
    symprec: float = 1e-5
    """symprec (float): Distance tolerance in Cartesian coordinates to find crystal symmetry for reshape cell operation (see spglib documentation for more info)"""
    fmax_pre_optimization : float = 0.1
    """fmax_pre_optimization (float): Maximal force component. Used as stopping criterion for pre geometry optimization (only used if second calculator is available)"""
    fmax: float = 0.01
    """fmax (float): Maximal force component. Used as a stopping criterion for the geometry optimization."""
    initial_step_size: float = None
    """initial_step_size (float): Initial step size of the geometry optimizer. If None the initial step size is estimated"""
    nhist_max: int = 10
    """nhist_max (int): Maximal length of history list in sqnm geometry optimizer"""
    lattice_weight: float = 2.
    """lattice_weight (float): Weight / size of the supercell that is used to transform lattice derivatives. Use a value between 1 and 2. Default is 2."""
    alpha_min: float = 1e-3
    """alpha_min (float): Lower limit on the step size. 1.e-2 is the default."""
    eps_subsp: float = 1e-3
    """esp_subsp (float): Lower limit on linear dependencies of basis vectors in history list."""
    enhanced_feedback:bool = False
    """enhanced_feedback (bool): Enhanced feedback (rise temperature according to T = T * beta_increase * (1. + 1. * ln(n_visits)))."""
    energy_threshold: float = 0.001
    """energy_threshold (float): If the energy difference of two structures is below the fingerprint is compared"""
    output_n_lowest_minima: int = 20
    """output_n_lowest_minima (int): Outputs the n lowest minima"""
    fingerprint_threshold: float = 5e-2
    """fingerprint_threshold (float): OMFP distance threshold for distinguishing minima"""
    verbose_output: bool = False
    """verbose_output (bool): If True MD and OPT logs are written."""
    new_start: bool = False
    """new_start (bool): Start from scratch even if restart files are present (deprecated)"""
    run_time: str = "infinite"
    """run-time (str): Runtime in the format (d-hh:mm:ss) or infinite for infinite runtime."""
    use_intermediate_mechanism: bool = False
    """use_intermediate_mechanism (bool): Sets if intermediate minimas will be stored and accepted if necessary."""
    write_graph_output: bool = True
    """write_graph_output (bool): Determines wether graph is calculated and written to file. Creates rather large files when hundred of thousands structures are found."""
    use_MPI: bool = False
    """use_MPI (bool): Sets if MPI with a master client model should be used"""
    logLevel: int = logging.INFO
    """logLevel (int): Sets the logging level."""
    _n_accepted: int = 0
    """_n_accepted (int): Private counting variable for number of accepted minima."""
    _n_rejected: int = 0
    """_n_rejected (int): Private counting variable for number of rejected minima."""
    _n_same: int = 0
    """_n_same (int): Private counting variable for number of unsuccessful escape attempts."""
    maxNumberOfMinima: int = 0
    """maxNumberOfMinima (int): Maximal number of minima that will be stored in database. Only the maxNumberOfMinima lowest energy structures will be stored. If a structure is higher in energy than the 
    maxNumberOfMinima it will considered as a new structure. If the number is 0 or negative, it will be considered as infinite."""
    logToStdOut: bool = True
    """logToStdOut (bool): If True log messages are written to stdout. If False log messages are written to a file. When MPI is used the log messages are always written to a file."""


    def getFixedParameterList(self):
        return ['n_S_orbitals', 'n_P_orbitals', 'width_cutoff', 'exclude', 'fingerprint_threshold', 'use_MPI']

if __name__ == '__main__':
    params = minimaHoppingParameters(exclude=["H", "He", 'Li'], dt0=1)
    print(params.getFixedParameterList())
    print(params.to_dict())