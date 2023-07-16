
Parameters
++++++++++

Minima Hopping Parameters
-------------------------

The minima hopping algorithm requires some input paramters. All the parameters are set to defaul but may be changed
to make the global optimization more efficient or more precise. Below you can find a list of all the input parameters
and a description of them.



.. csv-table:: Parameters
   :header: Parameter, Default, Description
   :widths: 15 10 60

    T0, 1000, Initial temperature to start the MH
    _T, -1.0, Temperature for restart (if -1 then last temperature is read)
    Ediff0, 0.1, Initital energy difference for accepting minima
    _eDiff, -1.0, Energy difference for accepting minima after restarts (if -1 then last Ediff is read)
    beta_decrease, 0.91, Factor for decreaseing the temperature
    beta_increase, 1.1, Factor for increaseing the temperature
    alpha_accept, 0.95, factor for decreasing Ediff
    alpha_reject, 1.05, factor for increasing Ediff
    n_soft, 20, number of softening steps
    soften_positions, 1e-2, step size for softening the positions
    soften_lattice, 1e-3, step size for softening the lattice (only relevant for periodic boundary conditions)
    n_S_orbitals, 1, number of s orbitals for constructing the OMFP
    n_P_orbitals, 0, number of p orbitals for constructing the OMFP
    width_cutoff, 3.5, cutoff for the OMFP
    exclude, [], List of elements to exclude in the OMFP
    dt0, 0.01, Initial time step for the MD
    _dt, -1.0, Time step of the MD for restart (if -1 then last dt is read)
    mdmin, 2, Number of minima visited before stopping the MD
    collect_md_data, False, flag to collect MD data which later could be used e.g. in machine learning
    symprec, 1e-5, Distance tolerance in Cartesian coordinates to find crystal symmetry for reshape cell operation (see spglib documentation for more info)
    fmax, 0.01, Maximal force component. Used as a stopping criterion for the geometry optimization.
    initial_step_size, None, Inital step size of the geometry optimizer. If None the initial step size is estimated
    nhist_max, 10, Maximal length of history list in sqnm geometry optimizer
    lattice_weight, 2, Weight / size of the supercell that is used to transform lattice derivatives. Use a value between 1 and 2. Default is 2.
    alpha_min,  1e-3, Lower limit on the step size.
    eps_subsp, 1e-3, Lower limit on linear dependencies of basis vectors in history list.
    enhanced_feedback, False, enhanced_feedback (bool):Enhanced feedback (rise temperature according to T = T * beta_increase * (1. + 1. * ln(n_visits))).
    energy_threshold, 0.001, if the energy difference of two structures is below the fingerprint is compared
    output_n_lowest_minima, 20, Outputs the n lowest minima
    fingerprint_threshold, 1e-3, OMFP distance threshold for distinguishing minima
    verbose_output, True, If True MD and OPT logs are written
    new_start, False, Start from scratch even if restart files are present (deporecated)
    run_time, 'infinity', Runtime in the format (d-hh:mm:ss) or inifinite for infinite runtime.
    use_intermediate_mechanism, False, Sets if intermediate minimas will be stored and accepted if necessary.
    write_graph_output, True, Determines wether graph is calculated and written to file. 
    use_MPI, False, Sets if MPI with a master client model should be used
    logLevel, logging.INFO, Sets the logging level
    _n_accepted, 0, Private counting variable for number of accepted minima
    _n_rejected, 0, Private counting variable for number of rejected minima
    _n_same, 0, Private conting variable for number of unsuccessful escape attempts
    maxNumberOfMinima, 0, Maximal number of minima that will be stored in database. If the number is 0 or negative it will be considered as infinite.



Restart
~~~~~~~
The python minimahopping algorithm automatically checks if restart files are aviable. If new_start is False as set by default the minimahopping is restarting from the last accepted minimum.
If a temperature, dt or Ediff different to the initial given parameters is required these parameters can be adjusted by changing _T, _dt or _eDiff respectivly.


Enhanced Feedback
~~~~~~~~~~~~~~~~~
Usually if a minimum is found more than once the temperature is increased the a fixed factor. However, there is build in mechanism
where the temperature is increased according to 

.. math:: 
   T = T * \beta_{increase} * (0.2 + ln(n_{visits}))

where n_visits is the number of times a particular minimum has already been visited. The same formula is also used as a feedback to 
increase Ediff:

.. math::
   E_{diff} = E_{diff} * \alpha_{rejected} * (0.2 + ln(n_{visits}))


Critical Parameters
~~~~~~~~~~~~~~~~~~~

.. caution::
   The paramter `minima_threshold` and the `fmax` as well as all the fingerprint paramters are crucial for
   distinguishing different minima. A tutorial how to adjust the `minima_threshold`  parameter can be adjusted to `fmax`
   can be found here (LINK TO THE TUTORIAL)



Fingerprint Adjustment
----------------------
In order to adjust the critical paramters `minima_threshold` and `fmax` as well as `energy_threshold` we strongly suggest to use the 
fingerprint adjustment tool.

.. csv-table:: Parameters Fingerprint Adjustment
   :header: Parameter, Default, Description
   :widths: 15 10 60

    fmax, 0.005, max force component for the local geometry optimization
    iteration, 100, number of md and optimizations performed
    temperature, 500, Temperature in Kelvin
    dt, 0.01, timestep for the MD
    md_min, 1, criteria to stop the MD trajectory (no. of minima)
    ns_orb, 1, number of s orbitals in the OMFP fingerprint
    np_orb, 1, number of p orbitals in the OMFP fingerprint
    width_cutoff, 3.5, width cutoff for the OMFP fingerprint
    maxnatsphere, 100, maximal number of atoms in one local atomic environment

It is important to keep the temperature, the timestep and the md_min low, so that after the optimization converges to the same minimum. 

