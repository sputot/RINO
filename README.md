# RINO
Robust INner and Outer Approximated Reachability

This is a verification library to compute inner and outer approximations of reachable sets for uncertain continous-time dynamical systems, with (possibly time-varying) perturbations and control inputs.

An early version is available from http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/software.html, cleaner version coming soon on this repository.

## Dependencies

You need g++, LAPACK and BLAS installed.

Install the FILIB++ Interval Library, available from http://www2.math.uni-wuppertal.de/~xsc/software/filib.html (we used Version 3.0.2), and set variable $FILIBHOME

Get and unzip the FADBAD++ automatic diffentiation package, available from http://www.fadbad.com/fadbad.html (we used FADBAD++ 2.1), and set variable $FADBADHOME

A modified of the third party package for Affine Arithmetic aaflib-0.1 (http://aaflib.sourceforge.net) has been included in the current directory, because some modifications were needed (additional functions and modifications of trigonometric functions). 
Future plans include separating more cleanly the initial version and our modifications...

## Installing

Go to directory aaflib-0.1 within the current package and compile by "make static". 

Returning to the main directory, you can now compile by "make" and obtain the "main" executable. 
The installation has been mostly tested on MacOS, but should also work on Ubuntu. 

## Running the analysis

### Running existing examples

Running an existing example is performed at command line, by 
```
./main type_system nb_system
```
where 
- type_system is either 0 (for a system of ODEs - Ordinary Differential Equations) or 1 (for a system of DDEs - Delay Differential Equations)
- nb_system is an integer specifying the predefined system identifier.

The corresponding systems (both for ODEs and DDES) are defined in ode_def.h (system and constant parameters) and ode_def.cpp (dimensions of the system, initial conditions, uncertain control inputs and perturbations, whether they are constant or time-varying, and control inputs or perturbations, and finally the integration settings - order of Taylor models, initial final time, time step etc). 
More documentation on how to use these (and better input mechanisms) should come...

### Visualizing results

After running an example, all results are in the subdirectory ‘output’. For the i-th variable, there are basically 6 interesting results files produced:  
- xiinner.out (maximal inner-approximation function of time)
- xiinner_robust.out (robust inner-approximation function of time) 
- xiinner_minimal.out (minimal inner-approximation function of time) 
- xiouter.out (maximal outer-approximation function of time)
- xiouter_robust.out (robust outer-approximation function of time)
- xiouter_minimal.out (minimal outer-approximation function of time)

A gnuplot file is also produced, that can be run from the source code repository by
```
cd output; gnuplot -p 'gnuplot_script.gp' ; cd ..
```
In particular, one figure per variable, xi.png, is produced (simply named x.png when the system is 1-dim as the running example) printing all these inner and outer-approximations with respect to time. 

These different type of inner and outer approximations are those described in "Inner and Outer Reachability for the Analysis of Control Systems" (see References below)

### Modifying / adding one's own example

...

## Authors and References

This package, written by [Sylvie Putot](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/), implements the ideas presented in:
-  Inner and Outer Reachability for the Analysis of Control Systems, by Eric Goubault and Sylvie Putot, Proceedings of the 22th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2019, Montreal [ [DOI](https://dl.acm.org/citation.cfm?id=3311794) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc19.pdf) ]
-  Inner and Outer Approximating Flowpipes for Delay Differential Equations, by Eric Goubault, Sylvie Putot and Lorenz Sahlmann, Proceedings of 30th International Conference on Computer Aided Verification, CAV 2018, Springer LNCS volume 10982 [ [DOI](https://www.springer.com/us/book/9783319961415) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/cav18.pdf) ]  
-  Forward inner-approximated reachability of non-linear continuous systems, by Eric Goubault and Sylvie Putot, Proceedings of the 20th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2017 [ [DOI](https://dl.acm.org/citation.cfm?id=3049811) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc17.pdf) ]

## Versions

Plans to propose cleaner versions in the future include (contributions and suggestions welcome, naturally) :
  - cleanly separate aaflib from our modifications
  - propose a parser for input systems (in SpaceEx's xml format)
  - better mechanism for specifying and taking into account time-varying inputs and perturbations
  - improve output format
 We also plan to improve precision and functionalities.

## License

This project is licensed under the GNU LGPLv3 license - see the [LICENSE](LICENSE) file for details

## Acknowledgments

Thanks to Franck Djeumou for helpful contributions, among which 
- the crazyflie quadcoptor system example (full non linear model of the dynamics with realistic control loop)
- modifications on the aaflib library for trigonometric functions
- better compatibility with Ubuntu
