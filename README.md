# Hydra package

```
     __              __          
    / /_  __  ______/ /________ _
   / __ \/ / / / __  / ___/ __ `/
  / / / / /_/ / /_/ / /  / /_/ / 
 /_/ /_/\__, /\__,_/_/   \__,_/  
       /____/
```

The following steps should help to install and run this code library on a Linux machine:

1) Install the required linear algebra libraries. On most operating systems, this can be done via

   ```sudo apt-get install libblas-dev liblapack-dev```

   On a cluster OpenBLAS is often pre-installed and needs just to be loaded. In this case, use

   ```module spider OpenBLAS```

   to find out, which <version> is available, then activate this via

   ```module load OpenBLAS/<version>```.

2) Clone the hydra repository into your home directory, e.g. via

   ```git clone git@github.com:HennesHajduk/hydra.git```.

3) Modify the scripts ```machine``` and ```workdir``` in the directory ```hydra/scripts``` as follows:

   a) Specify the FORTRAN-compiler to be used in ```machine```.

   b) Specify the directory where hydra output is dumped in ```workdir```.

4) Use of csh-shells is recommended for working with hydra. In that case, add the following line to the file ```~/.csh_aliases``` (create such a file if not yet existent):

   ```alias bat '(  nice +19 nohup time \!:2-$ ) >&\! \!^ &'```.

5) Add ```~/hydra```, ```~/hydra/scripts```, and ```.``` to the list of search paths. In your ```.cshrc```, this can be done as follows:

   ```setenv PATH <path1>:<path2>:/home/<user>/hydra:/home/<user>/hydra/scripts:.```.

   where ```/home/<user>``` is your home directory and ```<path1>``` and ```<path2>``` are two default search paths (any number other than two is possible here, just separate each path with a colon).

6) Modify your ```.cshrc``` file by including the following content:

   ```
   if( -f ~/.csh_aliases && -r ~/.csh_aliases && -o ~/.csh_aliases) then
      source ~/.csh_aliases
   endif
   
   setenv GFORTRAN_UNBUFFERED_ALL 1
   ```

7) If your shell already is of csh-type, then execute the command ```source ~/.cshrc``` otherwise, start such a shell by typing ```tcsh```.

8) You should now be able to use hydra, simply type ```hydra``` and set up a simulation by specifying input parameters.

Hydra is hierarchical and modular in the sense that first the method, geometry, and equation model are specified and then the physical parameters are read in.
Subsequently, the code for this setup is compiled and moved from a temporary directory to the one, where execution takes place and output is dumped.
This task is executed in the background and can be monitored via the ```log```-file.
There will also be some postprocessing scripts within each run-directory for visualization and data analysis purposes.


## IMPORTANT ADJUSTMENT FOR RUNNING HYDRA ON COMPUTE SERVERS

If you are running hydra on some server with a job scheduling system (e.g., Slurm, see [slurm.schedmd.com/documentation](https://slurm.schedmd.com/documentation.html)), the easiest way to use hydra as on your local machine is as follows:

To avoid running computationally demanding processes on login nodes (the admins of the cluster might react badly to you doing this), you can simply remove the ```bat```-command from the script ```hydra/scripts/flow-setup```.
The code will then be compiled as usual but not yet executed.
To perform this step, add a job that only runs the main executable and writes output to the ```log```-file to the queue. The script ```fox.slurm``` does precisely this and can be found in

```hydra/ca/basin/qgml/scripts/```.
