# ewrefxn: useful helper R fxns used in my daily work

* how to create such a package:  
  this package was created following the tutorial [here](http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/)  
     1. install and load devtools&roxygen2 packages with command:  
     	` install.packages('devtools');`  
	    ` devtools::install_github('klutometis/roxygen');`  
    	` library(devtools);library(roxygen2) `  
     2. change workdir to a father directory like f_wd:  
     	`setwd('f_wd')`
     3. create the project with the fxn `create` in devtools package
     4. add fxns files to the R subdirectory of the project
     5. run `document` at the root of project dir
     6. cd to the f_wd and run `install` to test the package  
* to install: `library(devtools); devtools::install_github('htc502/ewrefxn')`
