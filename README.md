# ![RINO](https://github.com/cosynus-lix/RINO/blob/master/RINO2.jpg)



This is a library to compute guaranteed inner and outer approximations of reachable sets for uncertain discrete-time or continous-time dynamical systems, with (possibly time-varying) perturbations and control inputs, where some of the control inputs can be specified as outputs of a neural network.

For continuous-time systems, it relies on Taylor expansions in time and affine arithmetic (i.e. zonotopes) in space based reachability analysis to compute outer envelopes of all possible trajectories of an uncertain system. Additionally, it uses a generalized mean-value theorem to deduce inner tubes, that contain only states guaranteed to be reached. It also studies robust versions of these tubes, when there can be control inputs and perturbations to the system. Finally, the control can be specified as the output of a neural network which inputs are the system state. 

# Dependencies and Installation

### Using docker

Get the RINO directory and run ```docker build .```.
An image ```shaxxx...``` is built which you can run by ```docker run -it --name rino shaxxx...```.
You can then execute RINO  from directory /home/RINO as indicated in the "Running the analysis" Section.

### Building from sources

You need g++, LAPACK and BLAS installed. Python visualization was tested with Python 3.8.8. 

Install the FILIB++ Interval Library, available from http://www2.math.uni-wuppertal.de/wrswt/software/filib.html (we used Version 3.0.2), and set variable $FILIBHOME

Get and unzip the FADBAD++ automatic diffentiation package, available from http://www.fadbad.com/fadbad.html (we used FADBAD++ 2.1), and set variable $FADBADHOME.
Copy files fadbad.h and fadiff.h from RINO/FADBAD_Modified/ into your FADBAD++ distribution (we modified these files to add differentiation of activation functions).

A slightly modified of the third party package for Affine Arithmetic aaflib-0.1 (http://aaflib.sourceforge.net) has been included in the current directory. 
Future plans include separating more cleanly the initial version and our modifications...
Go to directory aaflib-0.1 within the current package and compile by "make static". 

Returning to the main directory, you can now compile by "make" and obtain the "main" executable. 
The installation has been mostly tested on MacOS, but should also work on Ubuntu. 

# Running the analysis


For now, the dynamics of systems on which to perform reachability analysis are defined as C++ functions in ode_def.h/ode_def.cpp for continuous-time systems, discrete_system.h/discrete_system.cpp for discrete-time systems, and given some fixed ids.
Running an example is then performed at command line, in directory /home/RINO, by 
```
./rino [-systype system_type -syschoice system_id] [-nnfile-sfx filename.sfx] [config_file.txt]
```
where 
- system_type is either ode (for a system of ODEs - Ordinary Differential Equations) or dde (for a system of DDEs - Delay Differential Equations) or discrete (for a discrete-time dynamical system)
- system_id is an integer specifying the predefined system identifier (matching variable syschoice in file ```ode_def.h``` for ODEs and DDEs and  ```discrete_system.h``` for discrete-time systems )
- optional filename.sfx is the name of a file containing a neural network in the Sherlock sfx format (https://github.com/souradeep-111/sherlock/blob/master/sherlock-network-format.pdf) 
- optional configuration file allows to specify analysis parameters, inputs, parameters and initial conditions of the system, and the visualized outputs (all these can also be set in the code, but if both are specified, the configuration file overrides the code). 
- configuration files for some examples are available in directory ```Examples/ConfigFiles```. 
- at command line, either systype and syschoice should be specified, or a configuration file containing this information should be provided (if both are provided, config file information overrides command-line options)

The predefined systems (see running existing examples Section) are defined as C++ code:
  - for ODEs and DDEs in ```ode_def.h``` (system and constant parameters) and ```ode_def.cpp``` (parameters, initial conditions and input ranges)
  - for discrete-time systems in ```discrete_system.h``` and ```discrete_system.cpp``` 
The examples of neural network files are in directory ```Examples/Networks```.


## Running existing systems

### Continuous-time differential systems (ODEs)

  - the Brusselator example of Reference [HSCC 2017] below is run by 
  ```./rino -systype ode -syschoice 2```
  or if you want to use a configuration file to modify the parameter and initial conditions, by 
  ```./rino -configfile Examples/ConfigFiles/cfg_ode_2.txt```
  - the self-driving car example of Reference [HSCC 2019]  is run by ```./rino -systype ode -syschoice 6``` or ```./rino -configfile Examples/ConfigFiles/cfg_ode_6.txt```
  <!---  - the crazyflie model of Reference [HSCC 2019]  is run by "./rino -systype ode -syschoice 18 [Examples/ConfigFiles/cfg_ode_18.txt]" )-->
  
  ### Continuous-time delay differential systems with constant delays (DDEs)
  
  - the running example of [CAV 2018] is run  by ```./rino -systype dde -syschoice 1``` or ```./rino -configfile Examples/ConfigFiles/cfg_dde_1.txt```
  - Example 10 of  [CAV 2018] is run  by ```./rino -systype dde -syschoice 3``` or ```./rino -configfile Examples/ConfigFiles/cfg_dde_3.txt```
  - Example 9 (self-driving car with uncertain PID coefficients) of  [CAV 2018] is run  by ```./rino -systype dde -syschoice 8``` or ```./rino -configfile Examples/ConfigFiles/cfg_dde_8.txt```
  - The platoon examples of [CAV 2018] are run  by ```./rino -systype dde -syschoice 10 [Examples/ConfigFiles/cfg_dde_10.txt]``` (5 vehicles) or ```./rino -systype dde -syschoice 11 [Examples/ConfigFiles/cfg_dde_11.txt]``` (10 vehicles)

### Discrete-time dynamical systems

- the test model of [ADHS 2021] with Algorithm 1 is run by ```./rino -systype discrete -syschoice 15 -nbsteps 25 [-iter_method 1 -AEextension_order 1 -skew 1]``` or ```./rino -configfile Examples/ConfigFiles/cfg_discrete_15.txt```
- the SIR epidemic model of [ADHS 2021] with Algorithm 1 (Fig. 3) is run by ```./rino -configfile Examples/ConfigFiles/cfg_discrete_16_1.txt``` and the same model with slightly modified initial conditions and Algorithm 2 (Fig. 4) is run by  ```./rino -configfile Examples/ConfigFiles/cfg_discrete_16_2.txt```.
- the Honeybees Site Choice Model [ADHS 2021] with Algorithm 2 is run by ```./rino -configfile Examples/ConfigFiles/cfg_discrete_17.txt```

### Neural network controlled dynamical systems (continuous-time or discrete-time)

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


## Modifying / adding one's own example

Please take inspiration for the existing examples. In short :
- for ODEs, the dynamical system should be given an (unused) syschoice ID and the dynamics defined in file ```ode_def.h```, in class OdeFunc, by the vector yp  of time derivatives as function of state vector y and possibly parameters.  
The dimension of the state vector (sysdim) needs to be set in file ```ode_def.cpp``` in function define_system_dim(). The parameters and initial values (if not given at command line in a configuration file) need to be set in function init_system() in the same file. 
- for discrete-time systems, the process is very similar but in a different location. Give the system  an (unused) syschoice ID and define the dynamics in file ```discrete_system.h```, in class DiscreteFunc, by the vector z of next discrete state as a function of the current state x and possibly parameters. 
The dimension of the state vector (sysdim) and the initial conditions need to be set in file ```discrete_system.h``` in function init_discrete_system().  



# Parameters and Configuration File

## Sample configuration file: parameters common to all system types

```
# system type (ode, dde, discrete)
systype = ode

# system ID
syschoice = 2

# ranges for initial conditions
initial-values = [-0.1,0.1] [0,0.1] [1.9,2.1] [2.9,3.1]
# accuracy can be refined by partitioning the initial domain, you can specify the component you whish to subdivide and the number of partitions by:
# for example for 2 subdivisions for the last component (for the time being only one component at at time can be partitioned): 
# initial-values = [-0.1,0.1] [0,0.1] [1.9,2.1] ([2.9,3.1],2) 

# ranges for (constant or piecewise constant) inputs
# if all inputs are constant you can simply write inputs = [-0.1,0.1] [-0.1,0.1]
# for piecewise constant parameters, you can also specify number of steps over the time horizon (default is 1 for constant, x for piecewise constant with x different pieces)
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

## Parameters specific to ODEs (when systype is ode)

```
time-horizon = 5.

# time step
integration-step = 0.02

# order for Taylor models in time
order = 3

# refined mean-value theorem for inner-approximation (default is 0: classical mean value)
refined-mean-value = 1
```

## Parameters specific to DDEs (when systype is ode)

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

## Parameters specific to discrete-time systems (when systype is discrete)

```
# number of discrete time steps
nbsteps = 25

# iterating algorithm (Algorithm 1 or 2 in [ADHS2021])
iter-method = 1

# computing skewed box joint approximation or regular boxes (influences the precision for Algorithm 1, mostly for visualization for Algorithm 2)
skew = 1

# Order of the AE extension (1 for mean-value, 2 for Taylor-based higher-order extension)
AEextension-order = 1
```

## Parameters specific to neural network controlled dynamical systems (systype can be either ode or discrete)

```
# file containing the neural network in Sherlock sfx format (https://github.com/souradeep-111/sherlock/blob/master/sherlock-network-format.pdf) 
nnfile-sfx = Examples/Networks/tora_tanh.sfx

# when relevant, final offset and scaling of the output of the network:
nn-offset = 0. # also read from nn file (erasing this one)
nn-scaling = 11.  # also read from nn file (erasing this one)

# control time step
control-step = 0.1
```


# Visualizing results

After running an example, all results are in the subdirectory ‘output’. They are provided in the following files : 
- ```sumup.txt```: summary of configuration, running time and ranges at the final state of the analysis (part of this information can also be found with more significant digits in ```sumup.yaml```)
- ```samplesreachset.yaml```: sampled trajectories (used to assess accuracy of reachability results)
- ```approxreachset.yaml```: over and under-approximated reachset (projected, robust, joint ranges) and accuracy measures (eta, gamma) at each time step 

A python visualization file Visu_output.py is available in the GUI directory. It can be run from the analyzer (if variable create-png is set to 1 in the configuration file) but you can also run it separately, 
provided that the above data files are presen tn the output subdirectory of RINO. 
For example, for an interactive analysis (prints the figures on screen, otherwise the files are simply saved in the output directory) and to produce figures only for variables x[1] and x[2]), by:
```
cd GUI; python3 Visu_output.py --interactive=1 --printvar=-1-2; cd ..
```
When the script is run by analyzer, the options set above in command line can be set in the configuration file by: 
```
interactive-visualization = 1
variables-to-display = 1 2
```
In particular, for the k ranging from 1 to system dimension, the following results files print the projected ranges on dimension k as function of time:
- ```xk_max.png```  (e.g. ```x1_max.png```)  and ```xk_max_sample.png```: the maximal inner and outer-approximations, with and without sampled trajectories
- ```xk.png```, ```xk_sample.png```: additionally to the maximal inner and outer-approximations, the robust approximations when relevant, with and without sampled trajectories

For any couple (k,l) we also print 2-dimensional projections :
- ```xkxl.png```:  maximal (and when relevant robust) inner and outer-approximations of the joint range (xk,xl) as skewed boxes (see e.g. [CDC-LCSS 2020])
- ```xkxl_sample.png```: sampled trajectories for (xk,xl) 
- ```xkxl_approx_sample.png```: on the same graph the inner and outer-approximations of the joint range (xk,xl) as skewed boxes and sampled trajectories
- ```xkxl_box_sample.png```: same as above but the approximations are printed as boxes (useful in a few cases where the skewed boxes have a bad behavior)
- ```xkxl_finalstate.png```: box and skewed box inner and outer-approximations, robust when relevant, and sampled points at the final state of the analysis

Three-dimensional projections when relevant are also printed, only the corners of boxes are printed for more lisibility.

We also provide 
- ```xi_max.png``` and ```xi_subplots_min_max.png```: reachset for all variables on one graph
- ```eta.png```, ```gamma.png```: error measures (eta_o = (width of sampled set)/(width of over_approx) ; eta_i = (width of under_approx)/(width of sampled set);
gamma = (width of under_approx)/(width of over_approx)

These different type of inner and outer approximations are those described in "Inner and Outer Reachability for the Analysis of Control Systems" ([HSCC2019] in References below).

Note that the files produced can slightly vary depending on the system type (ode, dde, discrete-time)



# Authors and References

This package, written by [Sylvie Putot](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/), implements the ideas presented in:
- [ADHS 2021] Tractable higher-order under-approximating AE extensions for non-linear systems, In IFAC Conference on Analysis and Design of Hybrid Systems, ADHS 2021 [ [DOI](https://www.sciencedirect.com/science/article/pii/S2405896321012799) ]
- [CDC-LCSS 2020]  Robust under-approximations and application to reachability of non-linear control systems with disturbance, In IEEE Conference on Decision and Control, CDC 2020 and IEEE Control Systems Letters (Vol 4, Oct. 2020) [ [DOI](https://ieeexplore.ieee.org/document/9099222) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/csl2020.pdf) ]
-  [HSCC 2019] Inner and Outer Reachability for the Analysis of Control Systems, by Eric Goubault and Sylvie Putot, Proceedings of the 22th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2019, Montreal [ [DOI](https://dl.acm.org/citation.cfm?id=3311794) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc19.pdf) ]
-  [CAV 2018] Inner and Outer Approximating Flowpipes for Delay Differential Equations, by Eric Goubault, Sylvie Putot and Lorenz Sahlmann, Proceedings of 30th International Conference on Computer Aided Verification, CAV 2018, Springer LNCS volume 10982 [ [DOI](https://www.springer.com/us/book/9783319961415) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/cav18.pdf) ]  
-  [HSCC 2017] Forward inner-approximated reachability of non-linear continuous systems, by Eric Goubault and Sylvie Putot, Proceedings of the 20th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2017 [ [DOI](https://dl.acm.org/citation.cfm?id=3049811) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc17.pdf) ]

Versions until 2019 were focused on uncertain continuous-time dynamics. More recent versions include discrete-time dynamics and neural-network controlled systems. 
Please contact putot@lix.polytechnique.fr for suggestions or difficulties with the package.


# License

This project is licensed under the GNU LGPLv3 license - see the [LICENSE](LICENSE) file for details

