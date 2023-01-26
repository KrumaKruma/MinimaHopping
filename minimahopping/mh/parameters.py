from dataclasses import dataclass, field
from dataclasses_json import dataclass_json

@dataclass_json
@dataclass
class minimaHoppingParameters:
    T0: float = 1000
    _T: float = -1.0
    _eDiff: float = -1.0
    beta_decrease: float = 1./1.02
    beta_increase: float = 1.02
    Ediff0: float = .1
    alpha_accept: float = 1/1.02
    alpha_reject: float = 1.02
    n_soft: int = 20
    soften_positions: float = 1e-2
    soften_lattice: float = 1e-3
    n_S_orbitals: int = 1
    n_P_orbitals: int = 1 
    width_cutoff: int = 3.5
    maxnatsphere: int = 100
    exclude: list = field(default_factory=list)
    dt0: float = 0.01
    _dt: float = -1.0
    mdmin: int = 2
    collect_md_data: bool = False
    fmax: float = 0.001
    enhanced_feedback:bool = False
    energy_threshold: float = 0.0001
    output_n_lowest_minima: int = 20
    fingerprint_threshold: float = 1e-3
    verbose_output: bool = True
    new_start: bool = False
    run_time: str = "infinite"
    use_intermediate_mechanism: bool = False
    write_graph_output: bool = True
    use_MPI: bool = False
    totalWorkers: int = 1

if __name__ == '__main__':
    params = minimaHoppingParameters(exclude=["H", "He", 'Li'], dt0=1)
    print(params.to_dict())
