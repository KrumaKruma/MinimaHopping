from dataclasses import dataclass, field
from dataclasses_json import dataclass_json
import logging

@dataclass_json
@dataclass
class minimaHoppingParameters:
    T0: float = 1000                                  # Initial temperature to start MH
    _T: float = -1.0                                  # Temperature for restart (if -1 then last temperature is read)
    _eDiff: float = -1.0                              # Ediff for restart (if -1 then last Ediff is read)
    Ediff0: float = .1                                # Initital Ediff
    beta_decrease: float = 1./1.02                    # factor for decreaseing the temperature
    beta_increase: float = 1.02                       # factor for increaseing the temperature
    alpha_accept: float = 1/1.02                      # factor for decreaseing Ediff
    alpha_reject: float = 1.02                        # factor for increaseing Ediff
    n_soft: int = 20                                  # number of softening steps
    soften_positions: float = 1e-2                    # step size for softening the positions
    soften_lattice: float = 1e-3                      # step size for softening the lattice (only relevant for periodic boundary conditions)
    n_S_orbitals: int = 1                             # number of s orbitals for constructing the OMFP
    n_P_orbitals: int = 1                             # number of p orbitals for constructing the OMFP 
    width_cutoff: int = 3.5                           # cutoff for the OMFP
    maxnatsphere: int = 50                           # Truncation of the OMFP length 
    exclude: list = field(default_factory=list)       # List of elements to exclude in the OMFP
    dt0: float = 0.01                                 # Initial dt for the MD
    _dt: float = -1.0                                 # dt of the MD for restart (if -1 then last dt is read)
    mdmin: int = 2                                    # number of minima visited before stopping the MD
    collect_md_data: bool = False                     # flag to collect MD data which later could be used e.g. in machine learning
    fmax: float = 0.01                               # optimization stop criterion 
    enhanced_feedback:bool = False                    # enhanced feedback (rise temperature according to T = T * beta_increase * (1. + 1. * ln(n_visits)))
    energy_threshold: float = 0.01                    # if the energy difference of two structures is below the fingerprint is compared
    output_n_lowest_minima: int = 20                  # outputs the n lowest minima
    fingerprint_threshold: float = 5e-3               # OMFP distance threshold for distinguishing minima
    verbose_output: bool = False                      # if True MD and OPT logs are written
    new_start: bool = False                           
    run_time: str = "infinite"                        # runtime in the format (d-hh:mm:ss) or inifinite for infinite runtime
    use_intermediate_mechanism: bool = False          
    write_graph_output: bool = True                    
    use_MPI: bool = False
    logLevel: int = logging.INFO
    _n_accepted: int = 0
    _n_rejected: int = 0
    _n_same: int = 0

    def getFixedParameterList(self):
        return ['n_S_orbitals', 'n_P_orbitals', 'width_cutoff', 'exclude', 'fingerprint_threshold', 'use_MPI', 'maxnatsphere']

if __name__ == '__main__':
    params = minimaHoppingParameters(exclude=["H", "He", 'Li'], dt0=1)
    print(params.getFixedParameterList())
    print(params.to_dict())
