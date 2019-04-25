# RINO
Robust INner and Outer Approximated Reachability

This is a verification library to compute inner and outer approximations of reachable sets for uncertain continous-time dynamical systems, with (possibly time-varying) perturbations and control inputs.

An early version is available from http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/software.html, cleaner version coming soon on this repository.

## Dependencies

You need g++, LAPACK and BLAS installed.

Install the FILIB++ Interval Library, available from http://www2.math.uni-wuppertal.de/~xsc/software/filib.html (we used Version 3.0.2)

Get and unzip the FADBAD++ automatic diffentiation package, available from http://www.fadbad.com/fadbad.html (we used FADBAD++ 2.1), and set variable $FADBADHOME to the FADBAD++ directory

A modified of the third party package for Affine Arithmetic aaflib-0.1 (http://aaflib.sourceforge.net) has been included in the current directory, because some modifications were needed (additional functions and modifications of trigonometric functions). 
Future plans include separating more cleanly the initial version and our modifications...

## Installing

Go to directory aaflib-0.1 within the current package and compile by "make static". 

Returning to the main directory, you can now compile by "make" and obtain the "main" executable. 

## Running the analysis

### Running existing examples

### Visualizing results

### Adding one's own example

## References

This package implements the ideas presented in:
-  Inner and Outer Reachability for the Analysis of Control Systems, by Eric Goubault and Sylvie Putot, Proceedings of the 22th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2019, Montreal [ [DOI](https://dl.acm.org/citation.cfm?id=3311794) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc19.pdf) ]
-  Inner and Outer Approximating Flowpipes for Delay Differential Equations, by Eric Goubault, Sylvie Putot and Lorenz Sahlmann, Proceedings of 30th International Conference on Computer Aided Verification, CAV 2018, Springer LNCS volume 10982 [ [DOI](https://www.springer.com/us/book/9783319961415) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/cav18.pdf) ]  
-  Forward inner-approximated reachability of non-linear continuous systems, by Eric Goubault and Sylvie Putot, Proceedings of the 20th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2017 [ [DOI](https://dl.acm.org/citation.cfm?id=3049811) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc17.pdf) ]

## Versions

Plans to propose a cleaner version in the future include (contributions and suggestions welcome, naturally) :
  - cleanly separate aaflib from our modifications
  - propose a parser for input systems (in SpaceEx's xml format)
  - better mechanism for specifying and taking into account time-varying inputs and perturbations
  - improve output format

## License

This project is licensed under the GNU LGPLv3 license - see the [LICENSE](LICENSE) file for details

## Acknowledgments

Thanks to Franck Djeumou for his contribution with the crazyflie quadcoptor system example and the modifications on the aaflib library for trigonometric functions. 
