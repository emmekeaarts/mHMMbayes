## Test environments
* local OS X install, R 4.2.1, Apple silicon arm64 build
* win-builder (devel and release), results can be obtained from:  and 
* R-hub (Ubuntu Linux 20.04.1 LTS, R-release, GCC and Fedora Linux, R-devel, clang, gfortran and Debian Linux, R-devel, GCC ASAN/UBSAN), results can be obtained from: and and
 
## R CMD check results
There were no ERRORs or WARNINGs. 
There are 2 NOTEs:

* ONLY on win-builder: checking CRAN incoming feasibility ... NOTE
Maintainer: 'Emmeke Aarts <e.aarts@uu.nl>'. I believe this is because the Authors in DESCRIPTION is changed. 
* ONLY on Fedora Linux (R-hub): checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found. I cannot change that Tidy is not on the path, or update Tidy on the external Fedora Linux server. 
