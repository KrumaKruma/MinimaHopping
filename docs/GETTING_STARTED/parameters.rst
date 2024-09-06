Parameters
++++++++++

Minima Hopping Parameters
-------------------------

The minima hopping algorithm requires some input parameters. All the parameters are set to default but may be changed
to make the global optimization more efficient or more precise. Below you can find a list of all the input parameters
and a description of them.



.. csv-table:: Parameters
   :header: Parameter, Default, Description
   :widths: 15 10 60

   T0, 1000, Initial temperature to start MH
   _T, -1.0, Temperature for restart (if -1 then last temperature is read)
   Ediff0, 0.1, Initial energy difference for accepting minima
   _eDiff, -1.0, Energy difference for accepting minima after restarts.
   beta_decrease, 1./1.05, Factor for decreasing the temperature
   beta_increase, 1.05, Factor for increasing the temperature
   alpha_accept, 1./1.05, Factor for decreasing Ediff
   alpha_reject, 1.05, Factor for increasing Ediff
   fixed_cell_simulation, False, If True the simulation cell is fixed and no variable cell shape MD and optimization are performed.
   n_soft, 20, Number of softening steps
   soften_positions, 1e-3, Step size for softening the positions
   soften_lattice, 1e-3, Step size for softening the lattice (only relevant for periodic boundary conditions)
   n_S_orbitals, 1, Number of s orbitals for constructing the OMFP
   n_P_orbitals, 0, Number of p orbitals for constructing the OMFP
   width_cutoff, 4.0, Cutoff for the OMFP
   exclude, [], List of elements to exclude in the OMFP
   dt0, 0.08, Initial time step for the MD
   _dt, -1.0, Time step of the MD for restart.
   dt_min, 0.0001, Minimal time step of the MD.
   mdmin, 2, Number of minima visited before stopping the MD.
   md_max_steps, 10000, Maximum number of MD steps. If reached MH continuous
   collect_md_data, False, Flag to collect MD data which later could be used e.g. in machine learning
   margin, 0.3, "Margin for fixing of fragmentation, if particle closer than (2*max_rcov + 2*margin*max_rcov)"
   symprec, 1e-05, Distance tolerance in Cartesian coordinates to find crystal symmetry for reshape cell operation (see spglib documentation for more info).
   fmax_pre_optimization, 0.1, Maximal force component. Used as stopping criterion for pre geometry optimization (only used if second calculator is available).
   fmax, 0.01, Maximal force component. Used as a stopping criterion for the geometry optimization.
   initial_step_size, None, Initial step size of the geometry optimizer. If None the initial step size is estimated.
   nhist_max, 10, Maximal length of history list in sqnm geometry optimizer
   lattice_weight, 2.0, Weight / size of the supercell that is used to transform lattice derivatives. Use a value between 1 and 2. Default is 2.
   alpha_min, 1e-3, Lower limit on the step size.
   eps_subsp, 1e-3, Lower limit on linear dependencies of basis vectors in history list.
   enhanced_feedback, False, Enhanced feedback (rise temperature according to T = T * beta_increase * (1. + 1. * ln(n_visits))).
   energy_threshold, 0.001, If the energy difference of two structures is below the fingerprint is compared
   output_n_lowest_minima, 20, Outputs the n lowest minima
   fingerprint_threshold, 5e-2, OMFP distance threshold for distinguishing minima
   verbose_output, False, If True MD and OPT logs are written.
   new_start, False, Start from scratch even if restart files are present (deprecated).
   run_time, infinite, Runtime in the format (d-hh:mm:ss) or infinite for infinite runtime.
   use_intermediate_mechanism, False, Sets if intermediate minimas will be stored and accepted if necessary.
   write_graph_output, True, Determines wether graph is calculated and written to file. Creates rather large files when hundred of thousands structures are found.
   use_MPI, False, Sets if MPI with a master client model should be used.
   logLevel, logging.INFO, Sets the logging level.
   _n_accepted, 0, Private counting variable for number of accepted minima.
   _n_rejected, 0, Private counting variable for number of rejected minima.
   _n_same, 0, Private counting variable for number of unsuccessful escape attempts.
   maxNumberOfMinima, 0, "Maximal number of minima that will be stored in database. Only the maxNumberOfMinima lowest energy structures will be stored. If a structure is higher in energy than the maxNumberOfMinima it will considered as a new structure. If the number is 0 or negative, it will be considered as infinite."
   logToStdOut, True, "logToStdOut (bool): If True log messages are written to stdout. If False log messages are written to a file. When MPI is used the log messages are always written to a file."
    


Restart
~~~~~~~
The python minimahopping algorithm automatically checks if restart files are available. If new_start is False as set by default the minimahopping is restarting from the last accepted minimum.
If a temperature, dt or Ediff different to the initial given parameters is required these parameters can be adjusted by changing _T, _dt or _eDiff, respectively.


Enhanced Feedback
~~~~~~~~~~~~~~~~~
Usually if a minimum is found more than once the temperature is increased the a fixed factor. However, there is build in mechanism
where the temperature is increased according to 

.. math:: 
   T = T * \beta_{increase} * (0.2 + ln(n_{visits}))

where ``n_visits`` is the number of times a particular minimum has already been visited. The same formula is also used as a feedback to 
increase ``Ediff``:

.. math::
   E_{diff} = E_{diff} * \alpha_{rejected} * (0.2 + ln(n_{visits}))

Intermediate Mechanism
~~~~~~~~~~~~~~~~~~~~~~
First note, that in this implementation a minima hopping step is considered
done once a minimum is accepted. If a minima hopping step takes more than one
escape step and the intermediate mechanism is used, the minimum of lowest
potential energy of all escape steps is compared to the current minium and
accepted once it is smaller than ``Ediff``. If no intermediate mechanism is
used, the minimum of the last escape step is compared to the current minimum
and accepted once smaller as ``Ediff``. 

Critical Parameters
~~~~~~~~~~~~~~~~~~~

.. caution::
   The parameters ``fmax``, ``fingerprint_threshold`` as well as all other fingerprint parameters are crucial for
   distinguishing different minima. A tutorial how to adjust the ``fingerprint_threshold`` parameter 
   for clusters can be found :ref:`here <clusters adjust_fp>` and for bulk systems :ref:`here <bulk adjust_fp>`. 


Fingerprint Adjustment
----------------------
In order to adjust the critical parameters ``fingerprint_threshold`` and ``fmax`` as well as ``energy_threshold`` we strongly suggest to use the 
fingerprint adjustment tool.

.. csv-table:: Parameters Fingerprint Adjustment
   :header: Parameter, Default, Description
   :widths: 15 10 60

    fmax, 0.01, max force component for the local geometry optimization
    iteration, 10, number of md and optimizations performed
    temperature, 500, Temperature in Kelvin
    dt, 0.08, timestep for the MD
    md_min, 1, criteria to stop the MD trajectory (no. of minima)
    n_S_orbitals, 1, number of s orbitals in the OMFP fingerprint
    n_P_orbitals, 0, number of p orbitals in the OMFP fingerprint
    width_cutoff, 4.0, width cutoff for the OMFP fingerprint
    exclude, [], List of elements to exclude in the OMFP

It is important to keep the temperature, the timestep and the ``md_min`` low, so that after the optimization converges to the same minimum. 

