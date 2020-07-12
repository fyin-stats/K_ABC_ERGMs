##################### main function

# set the path as the local path of the folder
setwd()

#####################
# karate club network modeling

# K-ABC-IS (20 reps)
source("karate.kbergm.IS.R")

# K-ABC-AIS (20 reps)
source("karate.kbergm.AIS.R")

# AEA (20 reps)
source("karate.bergm.foreach.R")

# AEA (high quality)
source("karate.bergm.R")

# K-ABC vs AEA 
source("karate.bergm.kbergm.comp.R") # requires the results of the above scripts

#########################
# faux mesa high school network
# K-ABC-AIS (20 reps)
source("faux.kbergm.AIS.R")

# AEA (20 reps)
source("faux.bergm.foreach.R")

# AEA (high quality)
source("faux.bergm.R")

# K-ABC vs AEA
source("faux.bergm.kbergm.comp.R") # requires the results of the above scripts

#########################
# runtime comparison
#########################
source("kbergm.speed.test.main.R")

########################
source("relative.time.plot.R")
