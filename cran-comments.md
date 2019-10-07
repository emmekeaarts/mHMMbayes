## Resubmission
This is a resubmission. In this version I have:

* Reduced the size of the package to about 2 MB by adjusting the datasets used in te vignette.

## Test environments
* local OS X install, R 3.6.1
* win-builder (devel and release), results can be obtained from (updated for resubmission): 
 https://win-builder.r-project.org/C9LsZA4aQZ2k/00check.log and 
 https://win-builder.r-project.org/16oo7v0wHeXp/00check.log

## R CMD check results
There were no ERRORs or WARNINGs. 

There are 3 NOTEs:

* Possibly mis-spelled words in DESCRIPTION: Viterbi -> This is not mis-spelled,
it refers to the Viterbi algorithm. 
* ONLY on the release version of R on win-builder, one of the examples takes a little over 10 seconds (11.4 s). The example is scraped to the minimum running time and most examples are marked with "dontrun". I would appreciate it very much if this example could remain unchanged.  
* New submission (this is my FIRST package submission)
