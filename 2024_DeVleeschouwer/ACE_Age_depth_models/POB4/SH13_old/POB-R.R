#sets working directory on mac
setwd("/Users/Steve/Dropbox/BAS/Data/C14/MacBacon_2.2")
source("Bacon.R")

#POB4_bomb - test runs 
Bacon("POB4bomb",depths.file=FALSE,thick=10,postbomb=4, rotate.axes=TRUE,acc.mean=20,
      mem.mean=0.7,mem.strength=20,yr.max=10000,d.min=0,d.max=400,title.location='bottomright',dark=0.7)

