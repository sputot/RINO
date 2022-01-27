# ![RINO](https://github.com/cosynus-lix/RINO/blob/master/RINO2.jpg)


# WARNING (January 22) - Documentation currently being updated 


This is a library to compute guaranteed inner and outer approximations of reachable sets for uncertain discrete-time or continous-time dynamical systems, with (possibly time-varying) perturbations and control inputs, where some of the control inputs can be specified as outputs of a neural network.

For continuous-time systems, it relies on Taylor model based reachability analysis to compute outer envelopes of all possible trajectories of an uncertain system, as implemented in other reachability tools (but with the specificity to rely on affine arithmetic for the evaluation of the Taylor models). Additionally, it uses a generalized mean-value theorem to deduce inner tubes, that contain only states guaranteed to be reached. Finally, it also studies robust versions of these tubes, when there can be control inputs and perturbations to the system.

## Dependencies

You need g++, LAPACK and BLAS installed.

Install the FILIB++ Interval Library, available from http://www2.math.uni-wuppertal.de/wrswt/software/filib.html (we used Version 3.0.2), and set variable $FILIBHOME

Get and unzip the FADBAD++ automatic diffentiation package, available from http://www.fadbad.com/fadbad.html (we used FADBAD++ 2.1), and set variable $FADBADHOME

A modified of the third party package for Affine Arithmetic aaflib-0.1 (http://aaflib.sourceforge.net) has been included in the current directory, because some modifications were needed (additional functions and modifications of trigonometric functions). 
Future plans include separating more cleanly the initial version and our modifications...

## Installing

Go to directory aaflib-0.1 within the current package and compile by "make static". 

Returning to the main directory, you can now compile by "make" and obtain the "main" executable. 
The installation has been mostly tested on MacOS, but should also work on Ubuntu. 

## Running the analysis

### Running existing examples

For now, the models of systems to analyse are defined in ode_def.h/ode_def.cpp, and given some fixed ids.
Running an existing example is then performed at command line, by 
```
./rino [-systype system_type -syschoice system_id] [-nnfile-sfx filename.sfx] [config_file.txt]
```
where 
- system_type is either ode (for a system of ODEs - Ordinary Differential Equations) or dde (for a system of DDEs - Delay Differential Equations) or discrete (for a discrete-time dynamical system)
- system_id is an integer specifying the predefined system identifier (matching variable syschoice in file ode_def.h)
- optional filename.sfx is the name of a file containing a neural network in the Sherlock sfx format (https://github.com/souradeep-111/sherlock/blob/master/sherlock-network-format.pdf) 
- optional configuration file allows to specify analysis parameters, inputs, parameters and initial conditions of the system, and the visualized outputs (all these can also be set in the code, but if both are specified, the configuration file overrides the code). 
- configuration files for some examples are available in directory Examples/ConfigFiles. 
- at command line, either systype and syschoice should be specified, or a configuration file containing this information should be provided (if both are provided, config file information overrides command-line options)

In praticular:
  - the Brusselator example of Reference [HSCC 2017] below is run by "./rino -systype ode -syschoice 2 [Examples/ConfigFiles/cfg_ode_2.txt] "
  - the self-driving car example of Reference [HSCC 2019]  is run by "./rino -systype ode -syschoice 6 [Examples/ConfigFiles/cfg_ode_6.txt]"
  -   - the crazyflie model of Reference [HSCC 2019]  is run by "./rino -systype ode -syschoice 18 [Examples/ConfigFiles/cfg_ode_18.txt]" 
  - the self-driving car example of Reference [CAV 2018]  is run by "./rino -systype dde -syschoice 7 [Examples/ConfigFiles/cfg_dde_7.txt]" (here the model is with delays, hence the system_type 1)
  - the platoon examples  of Reference [CAV 2018] are run  by "./rino -systype dde -syschoice 10 [Examples/ConfigFiles/cfg_dde_10.txt]" (5 vehicles) or "./rino -systype dde -syschoice 11 [Examples/ConfigFiles/cfg_dde_11.txt]" (10 vehicles) 

The corresponding systems (both for ODEs and DDES) are defined in ode_def.h (system and constant parameters) and ode_def.cpp (dimensions of the system, initial conditions, uncertain control inputs and perturbations, whether they are constant or time-varying, and control inputs or perturbations, and finally the integration settings - order of Taylor models, initial final time, time step etc). 
More documentation on how to use these (and better input mechanisms) should come.

#### Sample configuration file (ODE)

```
# dimension of the ODE
system-dimension = 4
# dimension of the uncertain parameter and inputs
inputs-dimension = 2

time-horizon = 5.
# time step: only for ODEs
sampling-time = 0.02

# only for DDEs
delay = 0.3 
# only for DDEs (for ODEs, starting-time is set to 0)
starting-time = -0.3
# only for DDEs: defines the time step by delay/nb-time-subdivisions
nb-time-subdivisions = 3

# order for Taylor models
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

# which dimensions of the inputs (numbered starting from 1) are disturbances (the others are control inputs)
uncontrolled = 1
# which dimensions are time-varying inputs or parameters
variable =


# for the visualization: if 0, only .png files are produced
interactive-visualization = 1
# the dimensions we wish to visualize graphically (output files are produced for all dimensions in any case)
# to print all (default value when the option is not set in the config file): set variables-to-display = all
variables-to-display = 1 2
```

### Visualizing results

After running an example, all results are in the subdirectory ‘output’. For the i-th variable, there are basically 6 interesting results files produced:  
- xiinner.out (maximal inner-approximation function of time)
- xiinner_robust.out (robust inner-approximation function of time) 
- xiinner_minimal.out (minimal inner-approximation function of time) 
- xiouter.out (maximal outer-approximation function of time)
- xiouter_robust.out (robust outer-approximation function of time)
- xiouter_minimal.out (minimal outer-approximation function of time)

A python visualization file Visu_output.py is available in the GUI directory. 
It is run from the analyzer but you can also run it from the GUI directory, for example (for an interactive analysis and to produce .png files only for variables x[1] and x[2]) by
```
python3 Visu_output.py --interactive=1 --printvar=-1-2
```
When the script is run by analyzer, the options set above in command line can be set in the configuration file by: 
```
interactive-visualization = 1
variables-to-display = 1 2
```
In particular, at least figure per variable is produced, printing the evolution with time of
    - the maximal inner and outer-approximations (e.g. x1_max.png)
    - when minimal approximations are computed, the minimal, maximal inner and outer-approximations on a same graph (e.g. x1_min_max.png)
    - when robust approximations are computed, the minimal, robust and maximal inner and outer-approximations
We also provide graphs embedding all variables, as well as some error measures (see the files in the output directory).

These different type of inner and outer approximations are those described in "Inner and Outer Reachability for the Analysis of Control Systems" ([HSCC2019] in References below)

### Modifying / adding one's own example

For now, please take inspiration for the existing examples.

## Authors and References

This package, written by [Sylvie Putot](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/), implements the ideas presented in:
- [ADHS 2021] Tractable higher-order under-approximating AE extensions for non-linear systems, In IFAC Conference on Analysis and Design of Hybrid Systems, ADHS 2021 [ [pdf](https://arxiv.org/pdf/2101.11536.pdf)]
- [CDC-LCSS 2020]  Robust under-approximations and application to reachability of non-linear control systems with disturbance, In IEEE Conference on Decision and Control, CDC 2020 and IEEE Control Systems Letters (Vol 4, Oct. 2020) [ [DOI](https://ieeexplore.ieee.org/document/9099222) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/csl2020.pdf)]
-  [HSCC 2019] Inner and Outer Reachability for the Analysis of Control Systems, by Eric Goubault and Sylvie Putot, Proceedings of the 22th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2019, Montreal [ [DOI](https://dl.acm.org/citation.cfm?id=3311794) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc19.pdf) ]
-  [CAV 2018] Inner and Outer Approximating Flowpipes for Delay Differential Equations, by Eric Goubault, Sylvie Putot and Lorenz Sahlmann, Proceedings of 30th International Conference on Computer Aided Verification, CAV 2018, Springer LNCS volume 10982 [ [DOI](https://www.springer.com/us/book/9783319961415) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/cav18.pdf) ]  
-  [HSCC 2017] Forward inner-approximated reachability of non-linear continuous systems, by Eric Goubault and Sylvie Putot, Proceedings of the 20th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2017 [ [DOI](https://dl.acm.org/citation.cfm?id=3049811) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc17.pdf) ]

Please contact me (putot@lix.polytechnique.fr) for suggestions or difficulties with the package.

Versions until 2019 were focused on uncertain continuous-time dynamics. More recent versions include discrete-time dynamics and neural-network controlled systems. Input-output formatting has been updated, please contact me in case of difficulties running older examples. 

## License

This project is licensed under the GNU LGPLv3 license - see the [LICENSE](LICENSE) file for details

