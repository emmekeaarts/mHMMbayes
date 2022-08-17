## Test environments
* local OS X install, R 4.2.1, Apple silicon arm64 build
* win-builder (devel and release), results can be obtained from: https://win-builder.r-project.org/sd3lG6iWxxKj/00check.log and https://win-builder.r-project.org/dMv5pngtYyjA/00check.log
* R-hub (Ubuntu Linux 20.04.1 LTS, R-release, GCC and Fedora Linux, R-devel, clang, gfortran and Debian Linux, R-devel, GCC ASAN/UBSAN), results can be obtained from: https://artifacts.r-hub.io/mHMMbayes_0.2.0.tar.gz-ea64d34144cd49678e0bf162d46d238a/mHMMbayes.Rcheck/00check.log and https://artifacts.r-hub.io/mHMMbayes_0.2.0.tar.gz-b0f0a64d6efe41edb7de02886800a090/mHMMbayes.Rcheck/00check.log and https://builder.r-hub.io/status/original/mHMMbayes_0.2.0.tar.gz-0c1f7589ed704ec1b07ad5a01547382d
 
## R CMD check results
There were no ERRORs or WARNINGs. 
There are 2 NOTEs:

* ONLY on win-builder: checking CRAN incoming feasibility ... NOTE
Maintainer: 'Emmeke Aarts <e.aarts@uu.nl>'. I believe this is because the Authors in DESCRIPTION is changed. 
* ONLY on Fedora Linux (R-hub): checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found. I cannot change that Tidy is not on the path, or update Tidy on the external Fedora Linux server. 
