# xtb_md_scripts
tools for setup and analysis of molecular dynamics simulations using xtb

# Example usage
First, get a .xyz file of the molecule you want to simulate. Perhaps 5CB; and we'll generate this using our smi2xyz.py script:<br>
```python ~/code/hpc_tools/smi2xyz.py -i "CCCCCc1ccc(c2ccc(C#N)cc2)cc1" -o 5cb.xyz```
<br>
Then feed this into the orca or xtb builder (respective scripts):<br>
```python ../rjm_orca_md.py -i 5cb.xyz -n 3 -ns 20000```
```python ../rjm_xtb_md.py -i 5cb.xyz -n 3 -ns 20000 ```
<br><br>
There are a few options for things like number of molecules ```-n```, number of steps ```-ns```, temperature, using SHAKE... all sorts.

# TO DO
* need to write code for analysis
* need to write some bash or something for setting up many tens or hundreds of replicas of the simulation (each a unique starting point) so we get better statistics
* need to write scripts for submitting to slurm
* And need to benchmark on the HPC systems also (all just run on a laptop, single core).
