# xtb_md_scripts
tools for setup and analysis of molecular dynamics simulations using xtb

# Example usage
First, get a .xyz file of the molecule you want to simulate. Perhaps 5CB; and we'll generate this using our smi2xyz.py script:<br>
```python ~/code/hpc_tools/smi2xyz.py -i "CCCCCc1ccc(c2ccc(C#N)cc2)cc1" -o 5cb.xyz```
<br>
Then feed this into the orca or xtb builder (respective scripts):<br>
```python ../rjm_orca_md.py -i 5cb.xyz -n 3 -ns 20000```<br>
```python ../rjm_xtb_md.py -i 5cb.xyz -n 3 -ns 20000 ```
<br><br>
There are a few options for things like number of molecules ```-n```, number of steps ```-ns```, temperature, using SHAKE... all sorts.
<br><br>
To use xtb_replicates.sh:<br>
```sh ~/xtb_md/xtb_replicates.sh -i 5cb.xyz -nrep 500 -ns 100000 -enviro hpc -n 2```<br>
This creates 500 replica simulations, each with different starting configurations, and runs with default options for 100 ps (```-ns 100000```) for two molecules (```-n 2```) configured for execution on hpc (```-enviro hpc```). The written .sh files are configured for submission to slurm on the AIRE HPC at UoL.

# TO DO
* Analysis code is fairly boilerplate, it probably isn't the final version. It needs tweaking to handle the replica simulations.
* We need to run some benchmarks so we can assess the sort of time we expect for using multiple cores, for larger numbers of molecules and so on.
