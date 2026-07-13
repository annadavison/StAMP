#####The Sturgeon ABM for Migration Prediction (StAMP) #########################
# The following is an agent based model built to examine beluga sturgeon (Huso huso) 
# migration and population in the Danube River up to the Iron Gates dams. The 
# parameters included in this model are accumulated as far as possible from a 
# mixture of literature review and consultation with experts. The sampling of the 
# YoY simulated in this model reflects the sampling conducted by DDNI each year.

# Agents are individual adult and subadult sturgeon as well as collectives of eggs 
# and YoY from the same mother. Spatial information was created using ArcGIS Pro
# using real spawning, wintering and feeding sites.

# The output is several csv files which record key dates in migration and 
# spawning, the sturgeon present in the model, overall population numbers and
# the simulation of YoY sampling.

################################################################################

# First load in the packages used
library(lubridate)
library(tidyverse)
library(dplyr)
library(leaflet)
library(magrittr)
library(sf)
library(ggplot2)
#library(TrackReconstruction)
library(stringr)

# Then set the working directory - this contains the function script and the spatial information
# It also has a folder within it to store the output csv files and a folder for the environmental data
#setwd(
# Set your working directory
#"..."
#)
# Load in the script with the custom functions for the model
source("experiment_sturgeonfunctions.R")

################## Load in Data ################################################

##### Spatial Data #####

# The spatial data covers the Lower Danube accessible to the Beluga sturgeon 
# (up to the Iron Gates Dams) and Black Sea surrounding the Romanian coast. 
# Ensure that all has crs = 4326

# This model works with moving sturgeon following polyline routes from their
# locations to their destination. Polylines have been plotted over the Danube
# river which can be combined to produce these routes. One essential component
# for this route construction to work is to know the index of the vertex of
# the polyline which lies just before the spawning, wintering and feeding sites.

# To this end, we have an index and route column in the spawning and wintering
# site csv files which are loaded below

spawningsites <- read.csv("Spawning_sites.csv")
spawningsites <- st_as_sf(spawningsites,
                          coords = c("Long", "Lat"),
                          crs = 4326,
                          na.fail = FALSE)

winteringsites <- read.csv("Wintering_sites.csv")
winteringsites <- st_as_sf(winteringsites,
                           coords = c("Long", "Lat"),
                           crs = 4326,
                           na.fail = FALSE)

# Add the Polylines for Route Construction
SG <- read.csv("Sfantu_Gheorghe_route.csv")
sulina <- read.csv("Sulina_route.csv")
tulcea <- read.csv("Tulcea_stretch.csv")
K1 <- read.csv("Killia_entry_1.csv")
K2 <- read.csv("Killia_entry_2.csv")
K1_2 <- read.csv("Killia_entry_1-2.csv")
K3 <- read.csv("Killia_entry_3.csv")
K <- read.csv("Killia_route.csv")
mid <- read.csv("mid-bit.csv")
spawbranch <- read.csv("spawn_branch.csv")
wintbranch <- read.csv("winter_branch.csv")
danube <- read.csv("danube_proper.csv")

# Now use these polylines to create routes that are split into a start and end
# section. These can then be combined in the model to accommodate the different
# routes sturgeon will take according to their starting position and destination.

# Start route options (travel from the Black Sea)
sg_route <- add_routes(option = "route", SG, tulcea)
sulina_route <- add_routes(option = "route", sulina, tulcea)
killia_route1 <- add_routes(option = "route", K1, K1_2, K)
killia_route2 <- add_routes(option = "route", K2, K1_2, K)
killia_route3 <- add_routes(option = "route", K3, K)

# Combine all the routes to choose from for the first part of the sturgeon route.
routes <- c(sg_route, sulina_route, killia_route1, killia_route2, killia_route3)
# Create a data frame with all the startpoints
startpoints <- data.frame(X_coord = NA,
                          Y_coord = NA,
                          Routes = c("sg_route","sulina_route","killia_route1",
                                     "killia_route2","killia_route3"))
    
for (i in 1:length(routes)) {
  route <- routes[[i]]
  startpoint_X <- route[1, "X"]
  startpoint_Y <- route[1, "Y"]
  startpoints$X_coord[i] <- startpoint_X
  startpoints$Y_coord[i] <- startpoint_Y
}
startpoints <- st_as_sf(startpoints,
                        coords = c("X_coord", "Y_coord"),
                        crs = 4326)

# Now have the end route options (depending on the destination)
spawnbranchend_route <- add_routes(option = "route", mid, spawbranch, danube)
winterbranchend_route <- add_routes(option = "route", mid, wintbranch, danube)

# Since the sturgeon only have a predetermined spawning site, the wintering
# sturgeon also need to be able to select where they winter which is assumed to
# be related to the distance they need to travel to the spawning site.

# Simple assignment of wintering sites used by those travelling 
# to each spawning site
spawntowinter <- data.frame(spawnID = c(1:23), winterID = c(1,1,2,2,3,3,4,4,5,6,0,0,0,0,0,0,0,0,0,0,0,0,0))

# Now load in the spatial data to determine where the sturgeon are when they are
# feeding in the Black Sea. This is simplified and meant to represent where 
# sturgeon ready to spawn that year may be clustered more around the river in
# preparation for spawning.

# Load a csv with pre-sampled random points (10000)
feed <- read.csv("blackseapointsrough.csv")
feed_sf <- st_as_sf(feed,
                    coords = c("POINT_X", "POINT_Y"),
                    crs = 4326,
                    na.fail = FALSE)
    
# YoY Feeding site
# Coordinates and index are removed to obscure exact sensitive location
YoYfeedingsite <- c(..., ...)
YoYfeedingsite <- st_point(YoYfeedingsite)
# This is the index on the mid-bit route which is just before the feeding site
YoYfeedindex <- ...

# Set destination of the first part of the YoY travel to the YoY feeding site
destYoY <- YoYfeedingsite

# Set the end destination of the YoY travel which is around the mouth of the 
# Danube
desty <- st_coordinates(sg_route)
desty <- desty[1,1:2]
desty <- st_point(desty)

################################################################################
##### Agent Information #####
# Adults can be:
# Feeding, Migrating, Waiting, Spawning, ReturnMigrating, Wintering or Dead

# YoY can be:
# Incubating, Waiting, Firstmigrating, Riverfeeding, Secondmigrating, Feeding

# These are all used to control their migration steps.
#-------------------------------------------------------------------------------
# Extract the parameters from the SLURM script
#-------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
R <- as.numeric(args[1])
rep <- as.numeric(args[2])
# The underneath were used for testing
#R <- 1
#rep <- 1

if (R == 1){
  # This combo was seen to be the best for an increasing adult population in 
  # previous trials
  meanadultage <- 13140 # Increased the max age for sturgeon to 36
  fishcatch <- FALSE
} else {
  # For both baseline and decreasing there is an average age of 24 for adult death
  meanadultage <- 8760
  if (R == 3){
    fishcatch <- TRUE
  } else {
    fishcatch <- FALSE
  }
}

# Then use this to establish the mortality rate
# SUBADULTS
# For a stable population, there must be 2 adults surviving from the young
# produced by a single female over her lifetime. Females are assumed to spawn 5
# times in their life with an average starting age of 14 and death at 24 (from
# Paraschiv et al., 2006) and spawning every 2 or 3 years.
# The period of time which subadults exist (between YoY and adult) is taken as
# a period from 1 to 13 (10-13 for males and 13-15 for females) which is 12 
# years
# 100 is the number of YoY surviving to 1 year old
deathsubadult <- 1 - ((2/5)/100)^(1/4563)
# The last value reflects the average age for becoming a subadult at 12.5

# Adult sturgeon speed
# These speeds are calculated by taking the mid-point of the speed ranges in 
# Hont et al., 2019 and then converting km/h into km/day on the advice of the 
# DDNI team. These are then converted into m/day since all calculations will be 
# in m.
# Speed travelling upstream for adults
upspeed <- 37200  #m/day
# Speed travelling downstream for adults
downspeed <- 85920 # m/day

# YoY speed
# The YoY speed is based on the calculations used by DDNI when predicting YoY 
# arrival at the feeding site from the closest spawning site. This is estimated 
# to take 45 days to travel 188 km. 188/45 = 4.17777777 km/day
youngspeed <- 4178

# Water peak detection
# This was the method determined to be the best after trial and sensitivity analysis
detectpeak <- function(data, timestep) {
  # Not calculated in the first four days of the data, but not a problem since
  # the data starts in January
  if (timestep > 4) {
    one <- timestep - 1
    two <- timestep - 2
    three <- timestep - 3
    dip <- diff(data[one:timestep])
    trend1 <- diff(data[two:one])
    trend2 <- diff(data[three:two])
    if (dip <= 0 & trend1 > 0 & trend2 > 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}

#Spawning temperature trigger
spawntemp <- 7.0

#################### REFRESH ZONE ##############################################
# All code within here will be refreshed at the start of each run after it has
# been saved as a csv file to record the results of the run. These are tracking
# dataframes.

########################## Initialisation ####################################

##### Load Environmental Variables #####
# These have been retrieved for the region from 2015 01 01 until 2024 05 10
# The data is taken from a station at Isaccea at rkm 103

EnvVarHist <- read.csv("Full_Isaccea_EnvData_2024.csv")

# Assign the date as POSIXct so that it can be easily manipulated later
EnvVarHist$Date <- as.POSIXct(EnvVarHist$Date, format = "%d/%m/%Y")
# Add month and year columns for ease of checking values later
EnvVarHist$Month <- month(EnvVarHist$Date)
EnvVarHist$Year <- year(EnvVarHist$Date)
row.names(EnvVarHist) <- 1:nrow(EnvVarHist)

##### Load in agents #####
# Load in the pre-made, standard sturgeon csv files to be used as the starting 
# population.

# Use this for the baseline
sturgeon <- read.csv("newmort_baseline_sturpop.csv")

# When running for the increasing population, the AorD values need to be
# switched in the adult sturgeon to reflect the changes made to the mortality
# rate. Replace with standard AorD values so it is reproducible
if (R == 1){
  newAd <- read.csv("newadults36.csv")
  sturgeon[sturgeon$Class == "adult",] <- newAd
}
spawningtrack <- read.csv("baseline_standardspawnstart.csv", sep = ";")

# Convert to sf data frame
sturgeon_sf <- st_as_sf(sturgeon,
                        coords = c("X_Coord", "Y_Coord"),
                        crs = 4326,
                        na.fail = FALSE)
sturgeon_sf$Spawnsite <- as.numeric(sturgeon_sf$Spawnsite)

spawningtrack$Date_left <- as.POSIXct(spawningtrack$Date_left, format = "%d-%m-%Y")
spawningtrack$Date_spawned <- as.POSIXct(spawningtrack$Date_spawned, format = "%d-%m-%Y")
##############################################################################
# Load in the data frame for the existing YoY and the Eggs to be spawned into
younggroups <- read.csv("Start_younggroups.csv", row.names = NULL, sep = ";")
colnames(younggroups) <- c("ID", "Class", "Number", "Age", "Strategy", "State", "Spawnsite", "Leaving_feed", "X", "Y")
younggroups$X <- str_replace_all(younggroups$X, "([c(X = ])", "")
younggroups$Y <- str_replace_all(younggroups$Y, "([Y = )])", "")
younggroups <- st_as_sf(younggroups, coords = c("X", "Y"), crs = 4326)

# These two data frames contain information on all the agents in the model

##### Create the Data frames needed to record and run the model #####

# This data frame is paramount to recording and calculating the movement of the
# adult sturgeon
moveinfo <- data.frame(SturID = numeric(),
                       Total_Route = numeric(),
                       Days_travelled = numeric())
moveinfo <- st_sf(moveinfo, geometry = st_sfc())

# A separate movement data frame for the YoY since they have their mothers' IDs
moveYoY <- data.frame(YoYID = numeric(),
                      Total_Route = numeric(),
                      Days_travelled = numeric())

moveYoY <- st_sf(moveYoY, geometry = st_sfc())

# A data frame to track when the winter and spring migrations start each year
migrationtrack <- data.frame(Strategy = character(),
                             Year = numeric(),
                             Adults_migrating = numeric(),
                             Date_Start = as.POSIXct(character()))

# A data frame to record when eggs hatch and how long it took them
hatchingtrack <- data.frame(eggID = numeric(), Date = as.POSIXct(character()), Incubation_Time = numeric())

# A data frame to record the caught YoY and the number at the feeding site
# Currently this just has their IDs which can be matched to the spawning sites
track_catch <- data.frame(Date = as.POSIXct(character()), Pop_at_site = numeric(), Catch = numeric(), Mother_ID = numeric(), SpawnSite = numeric())

# Track the overall population in each of the life stages
track_pop <- data.frame(Time_Step = numeric(),
                        Date = as.POSIXct(character()),
                        Adults = numeric(),
                        SubAdults = numeric(),
                        YoY = numeric(),
                        Eggs = numeric())

# This part is to overcome an issue with having empty coordinates in row 1 of 
# sturgeon_sf
napoint <- NA

##### For Loop - Model Run #####
# Establish the time to see how long the model takes to run
start_time <- Sys.time()

# The model runs for each day in the environmental data
for (day in 1:nrow(EnvVarHist)) {
  print(EnvVarHist$Date[day])
  # Firstly mortality and aging are addressed
  # For each of the eggs and YoY
  if (nrow(younggroups) != 0) {
    for (group in 1:nrow(younggroups)) {
      ################### Eggs and YoY Mortality #############################
      # First comes the mortality for living through the day and then add 
      # the day to their age.
      # Mortality is currently established within the range given for the
      # first year in Jaric et al. (2009) ending at 100 subadults when
      # starting with the average number of eggs.
      
      # Mortality is currently a custom created exponential decrease 
      # function which gives an annual mortality between the .9996 and 1.0
      # as taken from Jaric et al. (2009) for this class
      survivingYoung <- survivingYoY(younggroups$Number[group], asYoYnum)
      # Replace the number with the new number of living individuals
      younggroups$Number[group] <- survivingYoung
      # Add a day to their age
      younggroups$Age[group] <- younggroups$Age[group] + 1
      
      #################### Egg to YoY ########################################
      if (younggroups$Class[group] == "eggs") {
        # If the class is egg, the chance of the egg hatching (changing to YoY) 
        # is dependent on the temperature.
        # Find number of days the egg has been alive and calculate average 
        # temperature for that period.
        eggdays <- younggroups$Age[group]
        temp <- mean(EnvVarHist[((day - eggdays):day), "Water_temp"])
        # Using a graph pulled from Ivanov 1987, a function for the
        # relationship between temperature and incubation time was created.
        hours <- incubationtime(temp)
        # If the time to hatch is less than the age of the eggs, then hatch
        if ((hours / 24) <= eggdays) {
          # Record the date of hatching in hatching track and assign the
          # correct class and state.
          eggdate <- data.frame(eggID = younggroups$ID[group], Date = EnvVarHist$Date[day], Incubation_Time = hours)
          hatchingtrack <- rbind(hatchingtrack, eggdate)
          younggroups$Class[group] <- "YoY"
          younggroups$State[group] <- "Waiting"
        }
      }
      # Now just processes for the YoY
      if (younggroups$Class[group] == "YoY") {
        ########################## MOVEMENT FOR YOY ##########################
        # First section is to ensure it doesn't include YoY spawned last year
        filteryear <- subset(hatchingtrack, year(hatchingtrack$Date) == EnvVarHist$Year[day])
        yearhatch <- filteryear$Incubation_Time[filteryear$eggID == younggroups$ID[group]]
        # Then select the YoY which have been waiting for two weeks already
        # This value was set by DDNI consultation
        if (((younggroups$Age[group] - (yearhatch[1]/24)) > 14) & younggroups$State[group] == "Waiting"){
          # Firstly, if the YoY are spawned at site 23 (downstream of the feeding
          # site), then they will simply continue down to the Black Sea
          if (younggroups$Spawnsite[group] == 23){
            # Secondmigrating is the term for travelling down to the Black Sea
            younggroups$State[group] <- "Secondmigrating"
            # Now give them routes to follow
            returncoords <- st_coordinates(spawnbranchend_route)
            returncoords <- returncoords[1:9, 1:2]
            # Now create a polyline again from the coordinates
            full_route <- st_linestring(returncoords)
            # Use the custom function to plot a course from the previous route 
            # to the YoY starting position
            startfeed <- asthesturswims(st_point(returncoords[9,]), younggroups$geometry[group])
            # Make the full route by combining with the Sf Gheorghe route
            # which all YoY are using for simplicity. The order is reversed
            # since this route is created and then run through backwards using
            # the day travel function
            # Split in half because of issues with the st_union function
            half <- st_union(sg_route, full_route)
            full_route <- st_union(half, startfeed)
            # Find the total length of the route to be used in daytravel
            coords_full <- st_coordinates(full_route)
            total <- totallength(coords_full)
            # Add the new data to the moveYoY data frame, essential to moving
            newdat <- data.frame(YoYID = younggroups$ID[group],
                                 Total_Route = total,
                                 Days_travelled = 0)
            newro <- st_sf(newdat, geometry = st_sfc(full_route))
            moveYoY <- rbind(moveYoY, newro)
          } else {
            # For all the other YoY spawned further upstream, they will travel
            # to the feeding site. Therefore they are firstmigrating
            younggroups$State[group] <- "Firstmigrating"
            # All spawning sites are located on spawnbranchend_route so take 
            # this and remove one index above the feeding site
            co <- st_coordinates(spawnbranchend_route)
            findex <- YoYfeedindex + 1
            # Now find the spawning site ID for each YoY group
            stuff <- younggroups$Spawnsite[group]
            # Use this to find the index of the spawning site to crop the route
            spawnstuff <- spawningsites$Index[spawningsites$ID == stuff]
            # For all sites upstream of site 22 this is where the index is
            # located
            if (stuff <= 21){
              yindex <- nrow(mid) + nrow(spawbranch) + spawnstuff
            } else {
              # When the site is 22 (23 is dealt with earlier) then the route
              # is on the spawbranch segment
              yindex <- nrow(mid) + spawnstuff
            }
            # Now use these indexes to crop the route to between the feeding
            # and the origin spawning site
            coo <- co[findex:yindex, 1:2]
            # Create a polyline again from these coordinates
            full_route <- st_linestring(coo)
            # Create the small starting and ending connections to make routes
            # from the start and end locations to the route constructed before
            # Make sure to create them as though the route is travelling from
            # the Black Sea inland and with the full_route creation
            ending <- asthesturswims(st_point(full_route[nrow(full_route),]), younggroups$geometry[group])
            starting <- asthesturswims(destYoY, st_point(full_route[1,]))
            # Split in half because of issues with st_union
            half <- st_union(starting, full_route)
            full_route <- st_union(half, ending)
            # Create a total length for the route
            coords_full <- st_coordinates(full_route)
            total <- totallength(coords_full)
            # Add the new data to the moveYoY data frame
            newdat <- data.frame(YoYID = younggroups$ID[group],
                                 Total_Route = total,
                                 Days_travelled = 0)
            newro <- st_sf(newdat, geometry = st_sfc(full_route))
            moveYoY <- rbind(moveYoY, newro)
            # No need to travel here since they will be compelled to do so in 
            # the next section
          }
        }
        if (younggroups$State[group] == "Firstmigrating") {
          # Use the daytravel function to proceed downstream
          youngcoord <- daytravel(youngspeed, moveYoY[moveYoY$YoYID == younggroups$ID[group],], destYoY, "down")
          # Update their coordinate
          younggroups$geometry[group] <- youngcoord
          # Add another day to the days travelled
          moveYoY$Days_travelled[moveYoY$YoYID == younggroups$ID[group]] <-
            moveYoY$Days_travelled[moveYoY$YoYID == younggroups$ID[group]] + 1
          # Once they reach the destination (feeding site), change their
          # state to riverfeeding to start that process.
          if (younggroups$geometry[group] == st_sfc(destYoY)){
            younggroups$State[group] <- "Riverfeeding"
            # Extract the day they will leave from a normal distribution based
            # on the assumption that YoY stay for an average of three weeks
            younggroups$Leaving_feed[group] <- younggroups$Age[group] + round(rnorm(n = 1, mean = 21, sd = 2.3))
            # Remove their moveYoY info so they are ready for the next part
            moveYoY <- moveYoY[moveYoY$YoYID != younggroups$ID[group], ]
          }
        }
        if (younggroups$State[group] == "Riverfeeding"){
          # Indicate that the YoY group arrived at the feeding site in the 
          # console.
          print("chilling")
          print(younggroups$ID[group])
          # When they have stayed for their determined time, then they will
          # leave
          if (younggroups$Age[group] >= younggroups$Leaving_feed[group]){
            # Set their state to Secondmigrating to indicate their travelling
            # downstream to the Black Sea now
            younggroups$State[group] <- "Secondmigrating"
            # Print information to the console for tracking their progress
            print("I'm leaving")
            print(younggroups$ID[group])
            # Now give them routes to follow
            # The first part is cropping the first spawnbranchend_route using
            # up to the index before the feeding site
            returncoords <- st_coordinates(spawnbranchend_route)
            returncoords <- returncoords[1:YoYfeedindex, 1:2]
            full_route <- st_linestring(returncoords)
            # Then create a small join from the feeding site to the start of
            # the previous route
            startfeed <- asthesturswims(st_point(returncoords[YoYfeedindex,]), younggroups$geometry[group])
            # Combine these with the Sf Gheorghe route which all the YoY take
            # split in half because of issues with the st_union function
            half <- st_union(sg_route, full_route)
            full_route <- st_union(half, startfeed)
            # Calculate the total route length
            coords_full <- st_coordinates(full_route)
            total <- totallength(coords_full)
            # Add the new data to the moveinfo data frame
            newdat <- data.frame(YoYID = younggroups$ID[group],
                                 Total_Route = total,
                                 Days_travelled = 0)
            newro <- st_sf(newdat, geometry = st_sfc(full_route))
            moveYoY <- rbind(moveYoY, newro)
            # No movement here as they will be compelled to move by the next 
            # section
          }
        }
        if (younggroups$State[group] == "Secondmigrating"){
          # Now use daytravel to proceed along the route until they reach the
          # last point of sg_route which is the mouth of the river.
          youngcoord <- daytravel(youngspeed, moveYoY[moveYoY$YoYID == younggroups$ID[group],], desty, "down")
          # Assign the new coordinate
          younggroups$geometry[group] <- youngcoord
          # Add another day travelled to moveYoY
          moveYoY$Days_travelled[moveYoY$YoYID == younggroups$ID[group]] <-
            moveYoY$Days_travelled[moveYoY$YoYID == younggroups$ID[group]] + 1
          # When they reach their destination, remove all their movement info
          # and assign them to feeding like the adults and subadults.
          if (younggroups$geometry[group] == st_sfc(desty)){
            younggroups$State[group] <- "Feeding"
            moveYoY <- moveYoY[moveYoY$YoYID != younggroups$ID[group], ]
          }
        }
        
        ############## YoY to Subadult #######################################
        # When the YoY become 1 year old, they become subadults and are added
        # to the sturgeon_sf data frame
        if (younggroups$Age[group] == 365) {
          # Use the convtosub function to generate the information for the new
          # subadults
          newsubs <- convtosub(younggroups[group,])
          sturgeon_sf <- rbind(sturgeon_sf, newsubs)
          # Set their number to zero so this group is removed at the end
          younggroups$Number[group] <- 0
        }
      }
    }
  }
  # At the end, all groups either completely dead or converted to subadult are
  # removed from the younggroups data frame. Done here because otherwise the 
  # row numbers get confused when removed halfway through the for loop.
  younggroups <- subset(younggroups, younggroups$Number != 0)
  
  ########################### Adults and Subadults ###########################
  for (i in 1:nrow(sturgeon_sf)) {
    # 
    if (sturgeon_sf$State[i] != "Dead") {
      ############# Subadult Aging and Mortality #############################
      if (sturgeon_sf$Class[i] == "subadult") {
        # Make sure the new YoY do not experience a double mortality on the day
        # they become subadults
        if (sturgeon_sf$Age[i] > 364) {
          # Mortality rate is calibrated for a stable population
          sdead <- runif(1) < deathsubadult
          # If true, set them as dead and remove their geometry
          if (sdead == TRUE) {
            sturgeon_sf$State[i] <- "Dead"
            sturgeon_sf$geometry[i] <- NA
          }
        }
        # Then add a day onto their age
        sturgeon_sf$Age[i] <- sturgeon_sf$Age[i] + 1
        # Then see if they become an adult if they are above the lower bracket
        # within which they mature to an adult - often cited in literature
        if (sturgeon_sf$Age[i] >= sturgeon_sf$AdultOrDead[i]) {
          sturgeon_sf$Class[i] <- "adult"
          # Assign them the day they will die from a normal distribution
          # based on a death at a mean age defined at the start of the simulation
          assigndead <- setAorD(sturgeon_sf$Sex[i], sturgeon_sf$Age[i], meanadultage)
          sturgeon_sf$AdultOrDead[i] <- assigndead
        }
      }
      ################### Adult Aging and Mortality ##########################
      if (sturgeon_sf$Class[i] == "adult") {
        if (sturgeon_sf$Age[i] >= sturgeon_sf$AdultOrDead[i]) {
          sturgeon_sf$State[i] <- "Dead"
          # Assign an empty point to dead sturgeon. Using the value napoint
          sturgeon_sf$geometry[i] <- napoint
          # This section was made to overcome an issue which occurs when the
          # sturgeon in row 1 dies. It will only not work if that sturgeon
          # dies first
          if (is.na(napoint)){
            napoint <- sturgeon_sf$geometry[i]
          }
          # Remove dead sturgeon from the moveinfo data frame
          if (any(moveinfo$SturID == sturgeon_sf$ID[i])){
            moveinfo <- moveinfo[moveinfo$SturID != sturgeon_sf$ID[i], ]
          }
        }
        # Then add a day onto their age
        sturgeon_sf$Age[i] <- sturgeon_sf$Age[i] + 1
      }
      # Finally, set any sturgeon which have just spawned in the previous day
      # to returnmigrating so they will start moving.
      if (sturgeon_sf$State[i] == "Spawning") {
        sturgeon_sf$State[i] <- "ReturnMigrating"
      }
    }
  }
  # Remove the dead sturgeon from the data frame after the for loop for
  # increased efficiency
  sturgeon_sf <- sturgeon_sf[sturgeon_sf$State != "Dead", ]
  rownames(sturgeon_sf) <- NULL
  ############# Environmental Conditions for Spring Migration ################
  # Currently the requirements are between January and June, over 6 degrees
  # and under 21 degrees.
  if (EnvVarHist$Month[day] %in% 1:6 &
      EnvVarHist$Water_temp[day] > 6.0 &
      EnvVarHist$Water_temp[day] < 21.0) {
    print("Spring migrate conditions")
    # Now assign the start of the migration date to migration track which also
    # indicates for the sturgeon to decide whether they are migrating this
    # year or not.
    if (!(nareplace(any(migrationtrack$Year == EnvVarHist$Year[day] 
                        & migrationtrack$Strategy == "S")) 
          | nareplace(any(EnvVarHist$Date[day] == migrationtrack$Date_Start)))) {
      # Now add that the migration has started
      print("First day spring migration")
      new_migr <- data.frame(Strategy = "S",
                             Year = EnvVarHist$Year[day],
                             Adults_migrating = NA,
                             Date_Start = as.POSIXct(EnvVarHist$Date[day], format = "%d/%m/%Y"))
      migrationtrack <- rbind(migrationtrack, new_migr)
      
      # Now establish the migrating population for this year
      for (a in 1:nrow(sturgeon_sf)) {
        # The wintering sturgeon have already established whether they are 
        # migrating this year so all wintering sturgeon are included. Spring
        # sturgeon are tested using readytospawn
        if ((readytospawn(sturgeon_sf[a, ], EnvVarHist$Date[day], spawningtrack) == TRUE &
             sturgeon_sf$State[a] == "Feeding" &
             sturgeon_sf$Strategy[a] == "S" &
             sturgeon_sf$Class[a] == "adult"
        ) | sturgeon_sf$State[a] == "Wintering") {
          sturgeon_sf$State[a] <- "Migrating"
          # Record which sturgeon are travelling for spawning this year
          spdate <- data.frame(SturID = sturgeon_sf$ID[a], 
                               Spawning_Site = sturgeon_sf$Spawnsite[a], 
                               Date_left = EnvVarHist$Date[day], 
                               Date_spawned = as.POSIXct(NA))
          spawningtrack <- rbind(spawningtrack, spdate)
        }
      }
      migrationtrack$Adults_migrating[migrationtrack$Date_Start == EnvVarHist$Date[day]] <- nrow(sturgeon_sf[sturgeon_sf$State == "Migrating",])
    }
    ## Spring Migration Movement ##
    for (b in 1:nrow(sturgeon_sf)) {
      # All sturgeon with the migrating state will now move
      if (sturgeon_sf$State[b] == "Migrating") {
        # Set the destination to their spawning site
        dest <- spawningsites$geometry[sturgeon_sf$Spawnsite[b]]
        # If the sturgeon is not already in the destination area
        if (!st_equals(sturgeon_sf$geometry[b], dest, sparse = FALSE)) {
          # And if the sturgeon does not already have a route it is following
          # then create a route
          if (!any(moveinfo$SturID == sturgeon_sf$ID[b])) {
            # For wintering sturgeon
            if (sturgeon_sf$Strategy[b] == "W") {
              # Simplify the process, all travel via the winterbranchend_route
              coords_first <- st_coordinates(winterbranchend_route)
              # First establish the index where the route needs to be cut off
              # before the spawning site they are travelling to
              windex <- nrow(wintbranch) + nrow(mid) + spawningsites$Index[sturgeon_sf$Spawnsite[b]]
              # Establish the index we start from with the wintering site
              int <- spawntowinter$winterID[spawntowinter$spawnID == sturgeon_sf$Spawnsite[b]]
              # The first section is for the wintering site which is on the 
              # wintbranch section
              if (sturgeon_sf$Spawnsite[b] %in% 9:10){
                wintdex <- nrow(mid) + winteringsites$Index[int]
              } else {
                wintdex <- nrow(wintbranch) + nrow(mid) + winteringsites$Index[int]
              }
              # Now subset the coords so they are just between the destination
              # and start wintering site
              coords_second <- coords_first[wintdex:windex, 1:2]
              full_route <- st_linestring(coords_second)
              # Now create a small join from either end of the route to the 
              # start and destination
              starting <- asthesturswims(sturgeon_sf$geometry[b], st_point(full_route[1,]))
              ending <- asthesturswims(st_point(full_route[nrow(full_route),]), dest)
              # Create the full route and calculate the total length
              # split in half due to issues with the st_union function
              half <- st_union(starting, full_route)
              full_route <- st_union(half, ending)
              coords_full <- st_coordinates(full_route)
              total <- totallength(coords_full)
            } else{
              # For all the spring sturgeon, first they must identify which 
              # route they will take which is dependent on which they are
              # closest to.
              # First establish the distance from the sturgeon to the start
              # points.
              routedistances <- st_distance(startpoints, sturgeon_sf$geometry[b])
              # Then find which of them are the closest
              nearest_route_index <- which.min(routedistances)
              nearest_route <- routes[[nearest_route_index]]
              # Convert to a linestring
              nearest_route <- st_linestring(nearest_route)
              # Record the startpoint as well for the nearest_route
              start_point <- startpoints[nearest_route_index, ]
              # This is then used to plot a direct course to the start of the 
              # route
              starting <- asthesturswims(sturgeon_sf$geometry[b], start_point)
              # Now the route needs to be cropped at the spawning site which
              # changes for the two sites which don't include further
              # stretches of the route.
              if (sturgeon_sf$Spawnsite[b] %in% 1:21){
                # This is the index recorded plus the total for the start of the route
                sindex <- nrow(mid) + nrow(spawbranch) + spawningsites$Index[sturgeon_sf$Spawnsite[b]]
              } else {
                if (sturgeon_sf$Spawnsite[b] == 22){
                  sindex <- nrow(mid) + spawningsites$Index[as.numeric(sturgeon_sf$Spawnsite[b])]
                } else {
                  sindex <- spawningsites$Index[sturgeon_sf$Spawnsite[b]]
                }
              }
              # Now use the index from the last part to crop the route
              endbit <- st_coordinates(spawnbranchend_route)
              endbit <- endbit[1:sindex, 1:2]
              endbit <- st_linestring(endbit)
              # Create the last final bit from the end of the path to the 
              # destination
              pointthing <- st_point(endbit[sindex,])
              ending <- asthesturswims(pointthing, dest)
              # We can then construct the full route the sturgeon will follow
              # Split in two steps because it can't handle in one go
              half <- st_union(starting, nearest_route)
              half2 <- st_union(endbit, ending)
              full_route <- st_union(half, half2)
              # and calculate the total length of the route
              coords_full <- st_coordinates(full_route)
              total <- totallength(coords_full)
            }
            # Now the moveinfo df can be filled with this information
            newdata <- data.frame(SturID = sturgeon_sf$ID[b],
                                  Total_Route = total,
                                  Days_travelled = 0)
            newrow <- st_sf(newdata, geometry = st_sfc(full_route))
            moveinfo <- rbind(moveinfo, newrow)
          }
          # The daytravel function is then used to move the sturgeon on for
          # all migrating sturgeon. Now the ones who needed moveinfo rows have
          # them
          coord <- daytravel(upspeed, moveinfo[moveinfo$SturID == sturgeon_sf$ID[b],], dest, "up")
          # Add another day travelled to the record for the moveinfo table
          moveinfo$Days_travelled[moveinfo$SturID == sturgeon_sf$ID[b]] <-
            moveinfo$Days_travelled[moveinfo$SturID == sturgeon_sf$ID[b]] + 1
          # Assign the new coord to the sturgeon
          sturgeon_sf$geometry[b] <- coord
        } else {
          # If the sturgeon is already at the site, then they will now wait to 
          # spawn
          print("I've Arrived!")
          sturgeon_sf$State[b] <- "Waiting"
          # All spring sturgeon retain their routes as they will use this to
          # travel back but the wintering sturgeon need a new route prepared
          # which will take them all the way back to the Black Sea.
          if (sturgeon_sf$Strategy[b] == "W") {
            # Remove their moveinfo data
            moveinfo <- moveinfo[moveinfo$SturID != sturgeon_sf$ID[b], ]
            # Take a destination from the feeding site randomly
            coord <- sample(feed$OID, size = 1)
            dest <- feed_sf$geometry[coord]
            # Go through the same process now as with the spring sturgeon
            # starting the spring migration
            routedistances <- st_distance(startpoints, dest)
            nearest_route_index <- which.min(routedistances)
            nearest_route <- routes[[nearest_route_index]]
            nearest_route <- st_linestring(nearest_route)
            start_point <- startpoints[nearest_route_index, ]
            starting <- asthesturswims(dest, start_point)
            rindex <- nrow(mid) + nrow(spawbranch) + spawningsites$Index[sturgeon_sf$Spawnsite[b]]
            endbit <- st_coordinates(spawnbranchend_route)
            endbit <- endbit[1:rindex, 1:2]
            endbit <- st_linestring(endbit)
            ending <- asthesturswims(st_point(endbit[rindex,]), sturgeon_sf$geometry[b])
            half <- st_union(starting, nearest_route)
            half2 <- st_union(endbit, ending)
            full_route <- st_union(half, half2)
            coords_full <- st_coordinates(full_route)
            total <- totallength(coords_full)
            # Add all this new information to the moveinfo data frame
            newdata <- data.frame(SturID = sturgeon_sf$ID[b],
                                  Total_Route = total,
                                  Days_travelled = 0)
            newrow <- st_sf(newdata, geometry = st_sfc(full_route))
            moveinfo <- rbind(moveinfo, newrow)
          } else {
            # For the spring sturgeon, no complexity, just return the days
            # travelled to 0 for calculating the return migration.
            moveinfo$Days_travelled[moveinfo$SturID == sturgeon_sf$ID[b]] <- 0
          }
        }
      }
    }
  }
  ################## SPAWNING ################################################
  # If it is March to May, the temperature is over 7 degrees and there is a 
  # peak in water level detected by the custom function, then all waiting
  # sturgeon will spawn
  if (EnvVarHist$Month[day] %in% 3:5 &
      EnvVarHist$Water_temp[day] > spawntemp &
      detectpeak(EnvVarHist$Water_level, day) == TRUE) {
    print("SPAWN")
    for (c in 1:nrow(sturgeon_sf)) {
      if (sturgeon_sf$State[c] == "Waiting") {
        sturgeon_sf$State[c] <- "Spawning"
        # If sex is female then use the custom geneggs function to create eggs
        if (sturgeon_sf$Sex[c] == "F") {
          neweggs <- geneggs(sturgeon_sf[c, ])
          print(neweggs$Number)
          print(neweggs$ID)
          # Add the groups of eggs to the younggroups data frame as a new
          # agent
          younggroups <- rbind(younggroups, neweggs)
          # Add the date the sturgeon spawned to the correct part of the
          # spawningtrack data frame.
          spawningtrack$Date_spawned[spawningtrack$SturID == sturgeon_sf$ID[c] & year(spawningtrack$Date_left) == EnvVarHist$Year[day]] <- as.POSIXct(EnvVarHist$Date[day])
          
        }
      }
    }
  }
  ################ Winter Migration ##########################################
  # In the months August to December, over 6 degrees and under 21 degrees, the
  # winter strategists will migrate to winter in the Danube.
  if (EnvVarHist$Month[day] %in% 8:12 &
      EnvVarHist$Water_temp[day] > 6.0 &
      EnvVarHist$Water_temp[day] < 21.0) {
    print("Winter migrate conditions")
    # Extra step if the migration has just begun to assign information to the
    # migrationtrack data frame
    if (!(nareplace(any(migrationtrack$Year == EnvVarHist$Year[day] & migrationtrack$Strategy == "W")) 
          | nareplace(any(EnvVarHist$Date[day] == migrationtrack$Date_Start)))) {
      # Now add that the migration has started
      new_migr <- data.frame(Strategy = "W",
                             Year = EnvVarHist$Year[day],
                             Adults_migrating = NA,
                             Date_Start = EnvVarHist$Date[day])
      migrationtrack <- rbind(migrationtrack, new_migr)
      # Now establish the migrating population from adult winter strategists
      # using the custom readytospawn function.
      for (d in 1:nrow(sturgeon_sf)) {
        if (sturgeon_sf$Strategy[d] == "W" &
            sturgeon_sf$Class[d] == "adult" &
            !(sturgeon_sf$State[d] == "Dead") &
            readytospawn(sturgeon_sf[d, ], EnvVarHist$Date[day], spawningtrack) == TRUE){
          sturgeon_sf$State[d] <- "Migrating"
        }
      }
      migrationtrack$Adults_migrating[migrationtrack$Date_Start == EnvVarHist$Date[day]] <- nrow(sturgeon_sf[sturgeon_sf$State == "Migrating",])
    }
    ## Winter Migration Movement ##
    # Now take all the sturgeon which will migrate this winter
    for (e in 1:nrow(sturgeon_sf)) {
      if (sturgeon_sf$State[e] == "Migrating") {
        # Establish their wintering site based on their spawning destination
        deststep <- spawntowinter$winterID[spawntowinter$spawnID == sturgeon_sf$Spawnsite[e]]
        dest <- winteringsites$geometry[deststep]
        # If the sturgeon is not already at the destination area
        if (!st_equals(sturgeon_sf$geometry[e], dest, sparse = FALSE)) {
          # and if the sturgeon does not already have a route it is following
          if (!any(moveinfo$SturID == sturgeon_sf$ID[e])) {
            # Then first find the distances to the startpoints of all the 
            # routes from the sturgeon
            routedistances <- st_distance(startpoints, sturgeon_sf$geometry[e])
            # Then find which of them are the closest
            nearest_route_index <- which.min(routedistances)
            nearest_route <- routes[[nearest_route_index]]
            # Convert to a linestring
            nearest_route <- st_linestring(nearest_route)
            # Record the start point of the nearest_route
            start_point <- startpoints[nearest_route_index, ]
            # This is then used to plot a direct course to the start of the 
            # route
            starting <- asthesturswims(sturgeon_sf$geometry[e], start_point)
            # All travel via the winter branch with some on the branch and 
            # therefore don't need the rest of the route
            if (dest %in% winteringsites$geometry[5:6]){
              wmindex <- nrow(mid) + winteringsites$Index[deststep]
            } else {
              wmindex <- nrow(mid) + nrow(wintbranch) + winteringsites$Index[deststep]
            }
            # Now use the index to crop the route to the wintering site
            endsec <- st_coordinates(winterbranchend_route)
            endsec <- endsec[1:wmindex, 1:2]
            endsec <- st_linestring(endsec)
            # Create the final section of the route to the destination
            ending <- asthesturswims(st_point(endsec[wmindex,]), dest)
            # We can then construct the full route the sturgeon will follow
            # and the total distance of the route
            half <- st_union(starting, nearest_route)
            half2 <- st_union(endsec, ending)
            full_route <- st_union(half, half2) 
            coords_full <- st_coordinates(full_route)
            total <- totallength(coords_full)
            # Now the moveinfo df can be filled with this information
            newdata <- data.frame(SturID = sturgeon_sf$ID[e],
                                  Total_Route = total,
                                  Days_travelled = 0)
            newrow <- st_sf(newdata, geometry = st_sfc(full_route))
            moveinfo <- rbind(moveinfo, newrow)
          } 
          # Then all the migrating sturgeon will move with the daytravel
          # function including those starting and continuing.
          coord <- daytravel(upspeed, moveinfo[moveinfo$SturID == sturgeon_sf$ID[e],], dest, "up")
          # Add another day to the days travelled
          moveinfo$Days_travelled[moveinfo$SturID == sturgeon_sf$ID[e]] <-
            moveinfo$Days_travelled[moveinfo$SturID == sturgeon_sf$ID[e]] + 1
          # Assign the new coord for the sturgeon
          sturgeon_sf$geometry[e] <- coord
        } else {
          # If the sturgeon is at the site, change their state to wintering
          # and they will wait there for spring.
          print("I'm wintering now!")
          sturgeon_sf$State[e] <- "Wintering"
        }
      }
    }
  }
  
  # Final considerations for each of the sturgeon
  for (f in 1:nrow(sturgeon_sf)) {
    ################## Return Migration ######################################
    if (sturgeon_sf$State[f] == "ReturnMigrating") {
      # The destination is now the first point of the original route
      dest <- st_coordinates(moveinfo$geometry[moveinfo$SturID == sturgeon_sf$ID[f]])
      dest <- dest[1, 1:2]
      dest <- st_point(dest)
      # Use daytravel but using the downspeed and the the direction argument
      # set to down in order to reverse the direction of travel
      coord <- daytravel(downspeed, moveinfo[moveinfo$SturID == sturgeon_sf$ID[f],], dest, "down")
      # Assign the new coord for the sturgeon
      sturgeon_sf$geometry[f] <- coord
      # Add another day travelled to the record for the moveinfo data frame
      moveinfo$Days_travelled[moveinfo$SturID == sturgeon_sf$ID[f]] <-
        moveinfo$Days_travelled[moveinfo$SturID == sturgeon_sf$ID[f]] + 1
      if (dest == coord) {
        print("I'm back!")
        # Return the sturgeon to their feeding state
        sturgeon_sf$State[f] <- "Feeding"
        # Remove the sturgeon entry from the moveinfo dataframe as return
        # journey is complete
        moveinfo <- moveinfo[moveinfo$SturID != sturgeon_sf$ID[f], ]
      }
    }
  }
  ################### Fishing ################################################
  # This is imposed during simulations where a reducing adult population is 
  # desired in order to reduce the population without directly altering the
  # mortality rates for the YoY in a similar way to the pressures imposed by
  # overfishing
  if (fishcatch == TRUE){
    fcatch <- runif(1, min=0, max = 1)
    # Calculated a reasonable rate of decline at 15 fished per year which is a 
    # probability of 0.04109589 per day that one will die.
    if (fcatch >= 0.9589041){
      aliveadults <- sturgeon_sf$ID[sturgeon_sf$Class == "adult"]
      caught <- sample(aliveadults, 1)
      sturgeon_sf$State[sturgeon_sf$ID == caught] <- "Dead"
      # Assign an empty point to dead sturgeon. Using the value napoint
      sturgeon_sf$geometry[sturgeon_sf$ID == caught] <- napoint
      # This section was made to overcome an issue which occurs when the
      # sturgeon in row 1 dies. Small chance of it not working, only if the first
      # sturgeon dies first but still work around
      if (is.na(napoint)){
        napoint <- sturgeon_sf$geometry[sturgeon_sf$ID == caught]
      }
      # Remove dead sturgeon from the moveinfo data frame
      if (any(moveinfo$SturID == caught)){
        moveinfo <- moveinfo[moveinfo$SturID != caught, ]
      }
      sturgeon_sf <- sturgeon_sf[sturgeon_sf$State != "Dead", ]
      rownames(sturgeon_sf) <- NULL
    }
  }
  ######################### YoY Sampling #####################################
  if (nrow(younggroups) > 0){
    # Subset the YoY to include those at the feeding site
    catchableYoY <- subset(younggroups, younggroups$geometry == st_sfc(destYoY))
    # If this subset is over 0 then we proceed
    if (nrow(catchableYoY) > 0){
      # Establish the positions object to accumulate over the for loop
      YoY_positions <- NA
      for (g in 1:nrow(catchableYoY)){
        # Now calculate the positions of the YoY groups for sampling.
        # Currently this is achieved using a basic system.  I have defined a 
        # stretch of river (500m) with a 2m sampling zone in the middle. YoY
        # groups have a random mean location assigned on the stretch and a SD
        # which increases with the number of days since hatching. Then the 
        # positions of each of the individuals in the group are given with
        # a normal distribution.
        
        # Randomly pick where the group is at the site
        mean_location <- sample(1:500, 1)
        # Use the days since hatching to calculate the standard deviation
        days_since_hatching <- EnvVarHist$Date[day] - max(hatchingtrack$Date[hatchingtrack$eggID == catchableYoY$ID[g]])
        sd <- 4 * (1 + 0.1 * days_since_hatching)
        # Now calculate the positions of the YoY in this group and add that to
        # the positions of all the other YoY already simulated.
        add_positions <- rnorm(catchableYoY$Number[g], mean = mean_location, sd = sd)
        if (any(is.na(YoY_positions))){
          YoY_positions <- add_positions
        } else {
          YoY_positions <- c(YoY_positions, add_positions)
        }
      }
      # Now sample the number of YoY in the sampling area (in the middle)
      start_position <- 245
      end_position <- 255
      catch <- sum(YoY_positions >= start_position & YoY_positions <= end_position)
      # Print the number caught
      cat("Number of YoY caught", catch, "\n")
      # Add the information to the tracking df for the YoY catch
      catch_row <- data.frame(Date = EnvVarHist$Date[day], 
                              Pop_at_site = sum(catchableYoY$Number), 
                              Catch = catch, 
                              Mother_ID = paste(catchableYoY$ID, collapse = ", "), 
                              SpawnSite = paste(catchableYoY$Spawnsite, collapse = ", "))
      track_catch <- rbind(track_catch, catch_row)
    }
  }
  # reset catch so it isn't double recorded
  catch <- 0
  
  ########### Record Information #############################################
  # Here the population for each class are added up and added to the track_pop 
  # data frame to track over the time steps
  countad <- sturgeon_sf %>%
    filter(Class == "adult", State != "Dead") %>%
    nrow()
  countsub <- sturgeon_sf %>%
    filter(Class == "subadult", State != "Dead") %>%
    nrow()
  steppop <- data.frame(Time_Step = day,
                        Date = EnvVarHist$Date[day],
                        Adults = countad,
                        SubAdults = countsub,
                        YoY = sum(younggroups$Number[younggroups$Class == "YoY"]),
                        Eggs = sum(younggroups$Number[younggroups$Class == "eggs"]))
  track_pop <- rbind(track_pop, steppop)
}
# Find the time of the full model to run over the whole EnvVarHist
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

# Save all the data from the run
# Hatchingtrack records the ID of the eggs (same as the mother), when they
# hatched and how long incubation took
filename <- paste0("hatch", R, "_", rep, ".csv")
write.csv(hatchingtrack, file = filename, row.names = FALSE)

# Migrationtrack records the start dates of the migration for spring and
# winter each year
filename <- paste0("migrate", R, "_", rep, ".csv")
write.csv(migrationtrack, file = filename, row.names = FALSE)

# Spawningtrack records which sturgeon migrated to spawn each year, which
# spawning sites they headed for and when they left
filename <- paste0("spawn", R, "_", rep, ".csv")
write.csv(spawningtrack, file = filename, row.names = FALSE)

# Sturgeon_sf contains the information on the full population of sturgeon
# alive in the model at the end of the simulation
filename <- paste0("sturgeon", R, "_", rep, ".csv")
write.csv(sturgeon_sf, file = filename, row.names = FALSE)

# Track_pop shows the numbers of sturgeon alive in each of the age classes
# for each day in the model
filename <- paste0("population", R, "_", rep, ".csv")
write.csv(track_pop, file = filename, row.names = FALSE)

# Track_catch records the number of YoY "caught" at the feeding site and the
# dates of the catch alongside the IDs of the agents.
filename <- paste0("catch", R, "_", rep, ".csv")
write.csv(track_catch, file = filename, row.names = FALSE)

# Print the time it took for the full number of runs
overallend_time <- Sys.time()
overallelapsed_time <- overallend_time - start_time
print(overallelapsed_time)