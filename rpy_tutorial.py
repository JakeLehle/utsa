#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 12:33:44 2021

@author: jake
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:59:24 2021

@author: jakel
"""
import os
os.chdir("/home/jake/Documents/SC_RNA_SEQ")
#%%
import rpy2
import rpy2.situation
for row in rpy2.situation.iter_info():
    print(row)
#rpy2 version 3.4.4
import rpy2.robjects as robjects #you want this version of the package to use R like you normally would
#### Importing Packages
from rpy2.robjects.packages import importr
# import R's "base" package
base = importr('base')
# import R's "utils" package
utils = importr('utils')
# There is a twist though. R object names can contain a “.” (dot) while in Python the dot means “attribute in a namespace”. Because of this, importr is trying to translate “.” into “_”. The details will not be necessary in most of the cases, but when they do the documentation for R packages should be consulted.


#%%
### Installing R Packages
# import rpy2's package module
import rpy2.robjects.packages as rpackages
# import R's utility package
utils = rpackages.importr('utils')
# select a mirror for R packages
utils.chooseCRANmirror(ind=1) # select the first mirror in the list
# R package names
packnames = ('ggplot2', 'hexbin', 'lazyeval')
# R vector of strings
from rpy2.robjects.vectors import StrVector
# Selectively install what needs to be install.
# We are fancy, just because we can.
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))
   
    
print('Out of loop')
#%%

#Exsample Script
pi = robjects.r('pi')
pi[0]
# This next part is defining a new function f 
robjects.r('''
        # create a function `f`
        f <- function(r, verbose=FALSE) {
            if (verbose) {
                cat("I am calling f().\n")
            }
            2 * pi * r
        }
        # call the function `f` with argument value 3
        f(3)
        ''')
# You can see that the product of f(3) is 18.49 or 2*pi*3
#Since that function f is now present in the R Global Environment, and can be output by assigning it to a value in the python environment and printing that value as shown:
r_f = robjects.globalenv['f']
print(r_f.r_repr())
#In addition not the function r_f is callable and can be used like a normal python function.
res = r_f(3)
print(res)
# you can also see all the variables in the global env like this
print(robjects.r.ls())
#%%
#In R, data are mostly represented by vectors, even when looking like scalars.#
#When looking closely at the R object pi used previously, we can observe that this is in fact a vector of length 1.
len(robjects.r['pi'])
#Accessing the one value in that vector has to be stated explicitly
robjects.r['pi'][0]
#Creating R vectors can be achieved simply
#Vectors containing strings
res = robjects.StrVector(['abc', 'def'])
print(res.r_repr())
#Vectors containing whole integers 
res = robjects.IntVector([1, 2, 3])
print(res.r_repr())
#Vectors with factions
res = robjects.FloatVector([1.1, 2.2, 3.3])
print(res.r_repr())
#Converting these numerical vectors into matricies 
v = robjects.FloatVector([1.1, 2.2, 3.3, 4.4, 5.5, 6.6])
m = robjects.r['matrix'](v, nrow = 2)
print(m)
#%%
#Calling R functions. Here I want to bring in the function sum and then use it like I would a python function but assigning R objects into it.
rsum = robjects.r['sum']
rsum(robjects.IntVector([1,2,3]))[0]
#Lets see that again with another function like sorting 
rsort = robjects.r['sort']
res = rsort(robjects.IntVector([1,2,3]), decreasing=True)
print(res.r_repr())
#%%
#Graphics and plots
# This code sets up the graphics plots so they will first create a file and then save all of the code into that file.
# This is the best way I have found to do this because R is set up as a threaded interative process and expects to be modified while python is set up to be less interactive and will crash the kernel if you give it multiple interactive x11 processes 
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
#This is what you will need to make non interactive plots
grdevices = importr('grDevices')

r = robjects.r

x = robjects.IntVector(range(10))
y = r.rnorm(10)

grdevices.png(file="/home/jake/Documents/SC_RNA_SEQ/file_1.png", width=512, height=512)
#These two lines are what is actually going to be plotted. The console will interprit these as In and Out R type arguments. 
r.layout(r.matrix(robjects.IntVector([1,2,3,2]), nrow=2, ncol=2))
r.plot(r.runif(10), y, xlab="runif", ylab="foo/bar", col="red")
#This closes the file
grdevices.dev_off()
