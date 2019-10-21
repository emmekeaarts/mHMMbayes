## Resubmission
This is a resubmission. In this version I have:

* changed the way information messages are given to the console, that is, print() inside functions is replaced by either message() or if(x)cat()
* user options wrt par() are resest in vignettes and examples as well 

## Test environments
* local OS X install, R 3.6.1
* win-builder (devel and release), results can be obtained from (updated for resubmission): https://win-builder.r-project.org/U8IN5uUZR0qo/00check.log and https://win-builder.r-project.org/tzysOvVEcuUD/00check.log

## R CMD check results
There were no ERRORs or WARNINGs. 

There are 2 NOTEs:

* Possibly mis-spelled words in DESCRIPTION: 
    Haan (15:55) -> name of an author
    Rabiner (13:38) -> name of an author
    Rietdijk (15:60) -> name of an author
    Viterbi (25:44) -> refers to the Viterbi algorithm 
    al (15:72) -> part of reference
    de (15:52) -> (part of) name of an author
    et (15:69) -> part of reference
* New submission (this is my FIRST package submission)



# OLD

## Resubmission
This is a resubmission. In this version I have:

* Added references to the Description field that describe the general method that is used. As yet, there are no references that describe the exact algorithm that is used. When this is published, I will add it to the Description field. 
* \dontrun{} is replaced by \donttest{}, and I wrote unit tests for each (external) function using testthat. 
* An immediate call of on.exit() is added to the plot functions such that the user's par() settings are reset.

## Test environments
* local OS X install, R 3.6.1
* win-builder (devel and release), results can be obtained from (updated for resubmission): https://win-builder.r-project.org/i1Lv3T784hIt/00check.log and https://win-builder.r-project.org/5bt84GS7v4tZ/00check.log

## R CMD check results
There were no ERRORs or WARNINGs. 

There are 2 NOTEs:

* Possibly mis-spelled words in DESCRIPTION: 
    Haan (15:55) -> name of an author
    Rabiner (13:38) -> name of an author
    Rietdijk (15:60) -> name of an author
    Viterbi (25:44) -> refers to the Viterbi algorithm 
    al (15:72) -> part of reference
    de (15:52) -> (part of) name of an author
    et (15:69) -> part of reference
* New submission (this is my FIRST package submission)


## Resubmission
This is a resubmission. In this version I have:

* Reduced the size of the package to about 2 MB by adjusting the datasets used in te vignette.

## Test environments
* local OS X install, R 3.6.1
* win-builder (devel and release), results can be obtained from (updated for resubmission): https://win-builder.r-project.org/xQe4gl2nljot/00check.log and https://win-builder.r-project.org/RKrvaBM1K64m/00check.log
 

## R CMD check results
There were no ERRORs or WARNINGs. 

There are 2 NOTEs:

* Possibly mis-spelled words in DESCRIPTION: Viterbi -> This is not mis-spelled, it refers to the Viterbi algorithm. 
* New submission (this is my FIRST package submission)
