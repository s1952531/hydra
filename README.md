# Hydra package

```
     __              __          
    / /_  __  ______/ /________ _
   / __ \/ / / / __  / ___/ __ `/
  / / / / /_/ / /_/ / /  / /_/ / 
 /_/ /_/\__, /\__,_/_/   \__,_/  
       /____/
```

The following steps should help to install and run this code library on practically any Linux machine:

1) Clone the hydra repository into your home directory, e.g. via

   ```git clone git@github.com:HennesHajduk/hydra.git```

2) Use of csh-shells is recommended for working with hydra. To do so, add ```~/hydra/scripts``` as well as the current directory to your search paths. In the file ```~/.cshrc```, this can be done as follows:

   ```setenv PATH <path1>:<path2>:/home/<user>/hydra/scripts:.```

   where ```/home/<user>``` is your home directory and ```<path1>``` and ```<path2>``` are two default search paths (any other number is possible here, just separate each path with a colon).

3) Run the script ```configure``` in ```hydra/scripts``` to define the fortran compiler and parent directory for all individual hydra runs. Make sure to specify the absolute path here.

4) Install the required linear algebra libraries. On Linux-based OS, this should be possible via

   ```sudo apt-get install libblas-dev liblapack-dev```

5) Add the following lines to the file ```~/.cshrc```:

   ```
   alias bat '(  nice +19 nohup time \!:2-$ ) >&\! \!^ &'   
   setenv GFORTRAN_UNBUFFERED_ALL 1
   ```

To update all scripts, execute the command ```source ~/.cshrc```. Otherwise, start such a shell by typing ```csh```.

You should now be able to use hydra, simply execute the command ```hydra``` and set up a simulation by specifying the input parameters inquired from you.
Subsequently, the code for this setup is compiled and moved from a temporary directory to the one where execution takes place and output is dumped.
This task is executed in the background and can be monitored via the ```log```-file.
There will also be some postprocessing scripts within each run-directory for visualization and data analysis purposes.

Hydra is hierarchical and modular in the sense that first the method, geometry, and system of equations are specified and then the physical parameters are read in. If you are interested in seeing the overall package structure, type ```hydra --h```.


## ADJUSTMENTS FOR RUNNING HYDRA ON COMPUTE SERVERS

If you are running hydra on some server with a job scheduling system (e.g., slurm, see [slurm.schedmd.com/documentation](https://slurm.schedmd.com/documentation.html)), some of the above steps require slight modifications. In our opinion, the easiest way to use hydra as on your local machine is by adapting the following of the above steps:

4) On a cluster, you typically do not have the rights required to install software, such as OpenBLAS, which instead is often pre-installed and may need to be activated. In this case, use commands akin to

   ```module spider OpenBLAS```

to find out, which <version> is available, then load it via

   ```module load OpenBLAS/<version>```

5) To avoid running computationally demanding processes on login nodes (the admins might go so far as to revoke your access to the cluster), you can simply remove the ```bat```-command from the script ```flow-setup``` that is executed for the run. As a consequence, hydra will only be compiled on a login node but not yet executed. As above, make sure that these changes take effect in your shell via ```source```.

To then start a simulation, type ```hydra``` once more, to setup a run as above. After specifying the final run-directory, ```cd``` there. The ```log``` file will now be empty. To start the simulation, add a job to the queue via a job script that only runs the main executable and writes the output to the ```log```-file. The script ```fox.slurm``` does precisely this and can be found in

```hydra/ca/basin/qgml/scripts/```

Jobs are scheduled via ```sbatch fox.slurm```

## Citations

Please cite with:
```c
@misc{hydra,
  key = {hydra},
  title = {Hydra [Software]},
  howpublished = {\url{https://github.com/daviddritschel/hydra}}
}
```

## Contributors

- [David Dritschel](https://www.st-andrews.ac.uk/mathematics-statistics/people/dgd1)
- [Hennes Hajduk](https://www.mn.uio.no/geo/english/people/aca/metos/hennesh/index.html)

## License

- [BSD](https://github.com/daviddritschel/hydra/LICENSE)

<br>
