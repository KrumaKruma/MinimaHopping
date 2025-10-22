# Example: Foundation Models
In this example a foundation model machine learning potential is used to perform Minima Hopping on a $\mathrm{CaTiO_3}$ crystal. 
Here we are using the [GRACE Potential](https://doi.org/10.1103/PhysRevX.14.021036), but many other implementations are available which should be easy to use as long as they provide an ASE calculator. 

For more details refer to the [documentation](https://python-minima-hopping.readthedocs.io/en/latest/EXAMPLES/foundation_models.html)

After installing the GRACE potential ([see here for instructions](https://gracemaker.readthedocs.io/en/latest/gracemaker/install/)), you can simply run the provided script using `python main.py`.

The script also supports MPI parallelization which can be enabled by setting `USE_MPI=True`. 
In this case, GPUs will only be assigned to the slave tasks while the master task will not require a GPU since it only handles the database. 

To run Minima Hopping with 4 slaves use a SLURM submit script that requests 5 tasks but only 4 GPUs.
```bash
#!/usr/bin/bash -l
#SBATCH --job-name=MH
#SBATCH --output=out.txt
#SBATCH --partition=my_partition
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --ntasks-per-node=5
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00

srun python3 main.py
```
