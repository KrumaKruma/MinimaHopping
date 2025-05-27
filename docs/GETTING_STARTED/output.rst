Output
++++++
Minima hopping produces several output files. 
Depending on the parameters set the output can vary.
This section is explaining how the output of Minima Hopping is structured and which output files are generated.
Minima Hopping is generating two output folders one is called ``minima`` and the other one ``output``.

Minima folder
-------------
In the ``minima`` folder the lowest minima are output ordered from lowest to highest energy.
The number of how many lowest minima are written is given by the parameter ``output_n_lowest_minima``.
In the case of periodic simulations the files are in ``.ascii`` format in the non-periodic case the files are written in the ``.xyz`` format.
The same folder also contains tow more files, namely ``all_minima.extxyz`` and ``all_minima_no_duplicates.extxyz``.
As the filenames indicate contains former all minima found by Minima Hopping and latter contains only all unique minima.


Output folder
-------------

Standard output
~~~~~~~~~~~~~~~
The minimal output in the ``output`` folder contains the following files:

.. csv-table:: Minimal output 
   :header: Filename, Description
   :widths: 15 60

    accepted_minima.extxyz, Contains all accepted minima
    all_minima.extxyz, Contains all minima
    all_minima_no_duplicates.extxyz, Contains all unique minima (accepted/rejected) 
    history.dat, Information summary of all minima found
    minimahopping.log, Log file of the current minimahopping run


.. caution::
    Be aware that in the case of a restart the ``minimahopping.log`` file is overwritten.
    If you require any information contained in this file copy it to a different folder before the restart.
    The file is not needed for performing the restart.

If the parameter ``verbose_output`` is set to True each step of the MD and the geometry optimization is written.
Hence, the following files are written:

.. csv-table:: verbose output 
   :header: Filename, Description
   :widths: 15 60

    MD_log.dat, Contains information about the current MD run e.g. energy conservation
    MD.extxyz, Contains the structure of each MD step
    geometry_optimization_log, Contains information about each geometry optimization step 
    geometry_optimization_trajectory.extxyz, Contains the structure of each geometry optimization step

.. caution::
    Be aware that in each MD-geometry optimization cycle these files are overwritten in order to keep the output minimal.

If it is desired to collect the structures generated during each MD e.g. for the construction of machine learning potentials the parameter ``collect_md_data`` can be set to true.
This generates a file named ``MD_collection.extxyz`` where the structure is added if two structures in the MD differ on average by 0.2 Angstroem.
This file is not overwritten during a restart and can get very large during a Minima Hopping run.

The ``output`` folder also contains a ``restart`` subfolder containing all the restart files required for restarting Minima Hopping.
It contains the ``params.json`` file where all the current parameters are stored and the ``poscur.extxyz`` file which is the last accepted minimum.
In the ``minima.pickle.shelve.dat`` the database is stored including also fingerprints as well as the structures.
If the parameter ``write_graph_output`` is set to True a ``trajectory.dat`` and ``graph.dat`` file is written containing information about the paths the algorithm is taking.

.. caution::
   Changing files in the ``restart`` folder can lead to corrupted restarts or even makes restarts of the Minima Hopping impossible. 
   Therefore, no changes in the files of the restart files should be made.  


MPI parallelization
~~~~~~~~~~~~~~~~~~~
In the case of MPI parallelized Minima Hopping runs the ``output`` directory contains subdirectories for each worker and a master directory.
The same files as described above are written, however, for each worker seperately.
The graph output is written in the ``restart`` subdirectory in the master directory.
Furthermore contains the master directory also a ``minimahopping.log`` file where the communication between the workers is written. 
This information can be very helpful if a MPI-parallelized Minima Hopping run is crashing and detecting the worker which is crashed. 
   


