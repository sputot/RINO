# ![RINO](https://github.com/cosynus-lix/RINO/blob/master/RINO2.jpg)


# WARNING - Documentation currently being updated 


This is a library to compute guaranteed inner and outer approximations of reachable sets for uncertain discrete-time or continous-time dynamical systems, with (possibly time-varying) perturbations and control inputs, where some of the control inputs can be specified as outputs of a neural network.

For continuous-time systems, it relies on Taylor model based reachability analysis to compute outer envelopes of all possible trajectories of an uncertain system, as implemented in other reachability tools (but with the specificity to rely on affine arithmetic for the evaluation of the Taylor models). Additionally, it uses a generalized mean-value theorem to deduce inner tubes, that contain only states guaranteed to be reached. Finally, it also studies robust versions of these tubes, when there can be control inputs and perturbations to the system.

## Dependencies and Installation

### Using docker

Get the RINO directory and run ```docker build .```.
An image ```shaxxx...``` is built which you can run by ```docker run -it --name rino shaxxx...```.
You can then execute RINO as indicated below.

### Building from sources

You need g++, LAPACK and BLAS installed.

Install the FILIB++ Interval Library, available from http://www2.math.uni-wuppertal.de/wrswt/software/filib.html (we used Version 3.0.2), and set variable $FILIBHOME

Get and unzip the FADBAD++ automatic diffentiation package, available from http://www.fadbad.com/fadbad.html (we used FADBAD++ 2.1), and set variable $FADBADHOME

A slightly modified of the third party package for Affine Arithmetic aaflib-0.1 (http://aaflib.sourceforge.net) has been included in the current directory. 
Future plans include separating more cleanly the initial version and our modifications...

Go to directory aaflib-0.1 within the current package and compile by "make static". 

Returning to the main directory, you can now compile by "make" and obtain the "main" executable. 
The installation has been mostly tested on MacOS, but should also work on Ubuntu. 

## Running the analysis

### Running existing examples

For now, the models of systems to analyse are defined in ode_def.h/ode_def.cpp for continuous-time systems, discrete_system.h/discrete_system.cpp for discrete-time systems, and given some fixed ids.
Running an existing example is then performed at command line, by 
```
./rino [-systype system_type -syschoice system_id] [-nnfile-sfx filename.sfx] [config_file.txt]
```
where 
- system_type is either ode (for a system of ODEs - Ordinary Differential Equations) or dde (for a system of DDEs - Delay Differential Equations) or discrete (for a discrete-time dynamical system)
- system_id is an integer specifying the predefined system identifier (matching variable syschoice in file ```ode_def.h```)
- optional filename.sfx is the name of a file containing a neural network in the Sherlock sfx format (https://github.com/souradeep-111/sherlock/blob/master/sherlock-network-format.pdf) 
- optional configuration file allows to specify analysis parameters, inputs, parameters and initial conditions of the system, and the visualized outputs (all these can also be set in the code, but if both are specified, the configuration file overrides the code). 
- configuration files for some examples are available in directory ```Examples/ConfigFiles```. 
- at command line, either systype and syschoice should be specified, or a configuration file containing this information should be provided (if both are provided, config file information overrides command-line options)

For example:
  - the Brusselator example of Reference [HSCC 2017] below is run by 
  ```./rino -systype ode -syschoice 2```
  or to modify the parameter and initial conditions, by 
  ```./rino -configfile Examples/ConfigFiles/cfg_ode_2.txt```
  - the self-driving car example of Reference [HSCC 2019]  is run by ```./rino -systype ode -syschoice 6``` or ```./rino -configfile Examples/ConfigFiles/cfg_ode_6.txt```
  - the running example of [CAV 2018] is run  by ```./rino -systype dde -syschoice 1``` or ```./rino -configfile Examples/ConfigFiles/cfg_dde_1.txt```
  - Example 10 of  [CAV 2018] is run  by ```./rino -systype dde -syschoice 3``` or ```./rino -configfile Examples/ConfigFiles/cfg_dde_3.txt```
  - Example 9 (self-driving car with uncertain PID coefficients) of  [CAV 2018] is run  by ```./rino -systype dde -syschoice 8``` or ```./rino -configfile Examples/ConfigFiles/cfg_dde_8.txt```
  - The platton examples of [CAV 2018] are run  by ```./rino -systype dde -syschoice 10 [Examples/ConfigFiles/cfg_dde_10.txt]``` (5 vehicles) or ```./rino -systype dde -syschoice 11 [Examples/ConfigFiles/cfg_dde_11.txt]``` (10 vehicles)
<!---  - the crazyflie model of Reference [HSCC 2019]  is run by "./rino -systype ode -syschoice 18 [Examples/ConfigFiles/cfg_ode_18.txt]" )-->
- the test model of [ADHS 2021] with Algorithm 1 is run by ```./rino -systype discrete -syschoice 15 -nbsteps 25 [-iter_method 1 -AEextension_order 1 -skew 1]``` or ```./rino -configfile Examples/ConfigFiles/cfg_discrete_15.txt```
- the SIR epidemic model of [ADHS 2021] with Algorithm 1 (Fig. 3) is run by ```./rino -configfile Examples/ConfigFiles/cfg_discrete_16_1.txt``` and the same model with slightly modified initial conditions and Algorithm 2 (Fig. 4) is run by  ```./rino -configfile Examples/ConfigFiles/cfg_discrete_16_2.txt```.
- the Honeybees Site Choice Model [ADHS 2021] with Algorithm 2 is run by ```./rino -configfile Examples/ConfigFiles/cfg_discrete_17.txt```

The corresponding predefined systems are defined:
  - for ODEs and DDEs in ```ode_def.h``` (system and constant parameters) and ```ode_def.cpp``` (parameters, initial conditions and input ranges)
  - for discrete-time systems in ```discrete_system.h``` and ```discrete_system.cpp``` 


#### Sample configuration file: parameters common to all system types

```
# system type (ode, dde, discrete)
systype = ode
# system ID
syschoice = 2

time-horizon = 5.
# time step: only for ODEs
integration-step = 0.02

# order for Taylor models in time
order = 3

# refined mean-value theorem for inner-approximation (default is 0: classical mean value)
refined-mean-value = 1

# ranges for initial conditions
initial-values = [-0.1,0.1] [0,0.1] [1.9,2.1] [2.9,3.1]
# accuracy can be refined by partitioning the initial domain, you can specify the component you whish to subdivide and the number of partitions by:
# for example for 2 subdivisions for the last component; for the time being only one component at at time can be partitioned
initial-values = [-0.1,0.1] [0,0.1] [1.9,2.1] ([2.9,3.1],2) 

# ranges for (constant or piecewise constant) inputs
# if all inputs are constant you can simply write inputs = [-0.1,0.1] [-0.1,0.1]
# for piecewise constant parameters, you can also specify number of steps (default is 1 for constant, x for piecewise constant with x different pieces)
inputs = ([-0.1,0.1],1) ([-0.1,0.1],10) # or inputs = [-0.1,0.1] ([-0.1,0.1],10)

# ranges for uncertain but constant parameters
params = [2,2] [3,4]

# which dimensions of the inputs (numbered starting from 1) are disturbances (the others are control inputs)
uncontrolled = 1
# which dimensions are time-varying inputs or parameters
variable =

# number of samples in each dimension taken for reachset estimation by trajectory samples
samples-per-dim = 20

# for the visualization: if 0, the visualization script is not called
create-png = 0

# for the visualization: if 0, only .png files are produced (if create-png is 1), if 1 figures are plot on screen 
interactive-visualization = 1

# number of (time) points saved and printed on graphs 
points-per-graph = 50

# the dimensions we wish to visualize graphically (output files are produced for all dimensions in any case)
# to print all (default value when the option is not set in the config file): set variables-to-display = all
variables-to-display = 1 2
```

#### Sample configuration file: parameters specific to ODEs (when systype is ode)

```
time-horizon = 5.

# time step
integration-step = 0.02

# order for Taylor models in time
order = 3

# refined mean-value theorem for inner-approximation (default is 0: classical mean value)
refined-mean-value = 1


# ranges for (constant or piecewise constant) inputs
# if all inputs are constant you can simply write inputs = [-0.1,0.1] [-0.1,0.1]
# for piecewise constant parameters, you can also specify number of steps (default is 1 for constant, x for piecewise constant with x different pieces)
inputs = ([-0.1,0.1],1) ([-0.1,0.1],10) # or inputs = [-0.1,0.1] ([-0.1,0.1],10)

# which dimensions of the inputs (numbered starting from 1) are disturbances (the others are control inputs)
uncontrolled = 1

# for the visualization: if 0, the visualization script is not called
create-png = 0
# for the visualization: if 0, only .png files are produced (if create-png is 1), if 1 figures are plot on screen 
interactive-visualization = 1
# the dimensions we wish to visualize graphically (output files are produced for all dimensions in any case)
# to print all (default value when the option is not set in the config file): set variables-to-display = all
variables-to-display = 1 2
```

#### Sample configuration file for DDEs

```
time-horizon = 5.

# only for DDEs (for ODEs, starting-time is set to 0)
starting-time = -0.3

# only for DDEs
delay = 0.3 

# defines the time-step by delay/nb-time-subdivisions
nb-time-subdivisions = 3

# order for Taylor models in time
order = 3

# refined mean-value theorem for inner-approximation (default is 0: classical mean value)
refined-mean-value = 1
```

#### Sample configuration file for discrete-time systems

```
# number of discrete time steps
nbsteps = 25

# iterating algorithm (Algorithm 1 or 2 in [ADHS2021])
iter-method = 1

# computing skewed box joint approximation or regular boxes
skew = 1

# Order of the AE extension (1 for mean-value, 2 for Taylor-based higher-order extension)
AEextension-order = 1
```

#### New: Running neural network controlled dynamical systems

```
TORA
./rino -configfile Examples/ConfigFiles/cfg_tora_tanh.txt
./rino -configfile Examples/ConfigFiles/cfg_tora_sigmoid.txt

B1
./rino -configfile Examples/ConfigFiles/cfg_B1_tanh.txt
./rino -configfile Examples/ConfigFiles/cfg_B1_sigmoid.txt

B2
./rino -configfile Examples/ConfigFiles/cfg_B2_sigmoid.txt

B3
./rino -configfile Examples/ConfigFiles/cfg_B3_tanh.txt
./rino -configfile Examples/ConfigFiles/cfg_B3_sigmoid.txt

B4:
./rino -configfile Examples/ConfigFiles/cfg_B4_tanh.txt
./rino -configfile Examples/ConfigFiles/cfg_B4_sigmoid.txt

B5:
./rino -configfile Examples/ConfigFiles/cfg_B5_tanh.txt
./rino -configfile Examples/ConfigFiles/cfg_B5_sigmoid.txt
  
ACC:
./rino -configfile Examples/ConfigFiles/cfg_acc_tanh.txt

Continuous-time Mountain Car
./rino -configfile Examples/ConfigFiles/cfg_MC_sigmoid.txt

Discrete-time Mountain Car
./rino -configfile Examples/ConfigFiles/cfg_discrete_mc.txt
```



### Visualizing results (ODEs)

After running an example, all results are in the subdirectory ‘output’. The data are provided in the following files : 
- ```sumup.txt```: summary of configuration, running time and ranges at the end of the analysis
- ```samplesreachset.yaml```: sampled trajectories (used to assess accuracy of reachability results)
- ```approxreachset.yaml```: over and under-approximated reachset (projected, robust, joint ranges) and accuracy measures (eta, gamma) at each time step 

A python visualization file Visu_output.py is available in the GUI directory. 
It is run from the analyzer but you can also run it from the GUI directory, for example (for an interactive analysis and to produce .png files only for variables x[1] and x[2]) by:
```
cd GUI; python3 Visu_output.py --interactive=1 --printvar=-1-2; cd ..
```
When the script is run by analyzer, the options set above in command line can be set in the configuration file by: 
```
interactive-visualization = 1
variables-to-display = 1 2
```
In particular, for the k ranging from 1 to system dimension, the following results files print the projected ranges on dimension k as function of time:
- ```xk_max.png``` and ```xk_max_sample.png```: the maximal inner and outer-approximations, with and without sampled trajectories (e.g. ```x1_max.png```)
- ```xk_min_max.png```: when minimal approximations are computed, the minimal, maximal inner and outer-approximations on a same graph
- ```xk_rob_max.png```, ```xk_rob_min_max.png```: when robust approximations are computed, the minimal, robust and maximal inner and outer-approximations

For any couple (k,l) we also print 2-dimensional projections :
- ```xkxl.png```:  maximal (and when relevant minimal, robust) under and over approximations of (xk,xl) 
Three-dimensional projections when relevant are also printed, only the corners of boxes are printed for more lisibility.

We also provide 
- ```xi_max.png``` and ```xi_subplots_min_max.png```: reachset for all variables on one graph
- ```eta.png```, ```gamma.png```: error measures (eta_o = (width of sampled set)/(width of over_approx) ; eta_i = (width of under_approx)/(width of sampled set);
gamma = (width of under_approx)/(width of over_approx)

These different type of inner and outer approximations are those described in "Inner and Outer Reachability for the Analysis of Control Systems" ([HSCC2019] in References below).

Note that the files produced can slightly vary depending on the system type (ode, dde, discrete-time)

### Modifying / adding one's own example

For now, please take inspiration for the existing examples.

## Authors and References

This package, written by [Sylvie Putot](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/), implements the ideas presented in:
- [ADHS 2021] Tractable higher-order under-approximating AE extensions for non-linear systems, In IFAC Conference on Analysis and Design of Hybrid Systems, ADHS 2021 [ [DOI](https://www.sciencedirect.com/science/article/pii/S2405896321012799) ]
- [CDC-LCSS 2020]  Robust under-approximations and application to reachability of non-linear control systems with disturbance, In IEEE Conference on Decision and Control, CDC 2020 and IEEE Control Systems Letters (Vol 4, Oct. 2020) [ [DOI](https://ieeexplore.ieee.org/document/9099222) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/csl2020.pdf) ]
-  [HSCC 2019] Inner and Outer Reachability for the Analysis of Control Systems, by Eric Goubault and Sylvie Putot, Proceedings of the 22th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2019, Montreal [ [DOI](https://dl.acm.org/citation.cfm?id=3311794) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc19.pdf) ]
-  [CAV 2018] Inner and Outer Approximating Flowpipes for Delay Differential Equations, by Eric Goubault, Sylvie Putot and Lorenz Sahlmann, Proceedings of 30th International Conference on Computer Aided Verification, CAV 2018, Springer LNCS volume 10982 [ [DOI](https://www.springer.com/us/book/9783319961415) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/cav18.pdf) ]  
-  [HSCC 2017] Forward inner-approximated reachability of non-linear continuous systems, by Eric Goubault and Sylvie Putot, Proceedings of the 20th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2017 [ [DOI](https://dl.acm.org/citation.cfm?id=3049811) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc17.pdf) ]

Versions until 2019 were focused on uncertain continuous-time dynamics. More recent versions include discrete-time dynamics and neural-network controlled systems. 
Please contact putot@lix.polytechnique.fr for suggestions or difficulties with the package.


## License

This project is licensed under the GNU LGPLv3 license - see the [LICENSE](LICENSE) file for details

