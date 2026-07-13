################### StAMP Functions ############################################
# These functions have been developed for use in the Sturgeon ABM for Migration
# Prediction (StAMP)

########################### ROUTES AND MOVEMENT ################################

#################### Copied from TrackReconstruction for Anunna ################
# To avoid requiring additional package downloads, this function was copied from
# the TrackReconstruction package v 1.3 (Battaile, 2021)
CalcBearing <-
  function(initialLat,initialLong,finalLat,finalLong)
  {
    initialLat<-initialLat/360*2*pi
    initialLong<-initialLong/360*2*pi
    finalLat<-finalLat/360*2*pi
    finalLong<-finalLong/360*2*pi
    atan2(sin(finalLong-initialLong)*cos(finalLat),cos(initialLat)*sin(finalLat)-sin(initialLat)*cos(finalLat)*cos(finalLong-initialLong))
  }

########################## Adding Routes Functon ###############################
# This function will take csv files for polylines produced using ArcGIS Pro and
# combine them into an sf linestring. This is used for creating the routes for
# the sturgeon to follow when migrating
add_routes <- function(option, ...){
  line <- rbind(...)
  line <- line[ ,5:6]
  
  if (option == "route"){
    line <- as.matrix(line)
    route <- st_linestring(line, dim = "XY")
    return(route)
  }
  
  if (option == "point"){
    points <- st_as_sf(line, coords = c("X", "Y"), crs = 4326)
    return(points)
  }
}

############################## Haversine Formula ###############################
# This is used to determine the distance between two coordinates in metres given
# that the earth is a sphere.

haversine_distance <- function(lat1, lon1, lat2, lon2) {
  # Convert from degrees to radians
  lat1 <- lat1 * pi / 180
  lon1 <- lon1 * pi / 180
  lat2 <- lat2 * pi / 180
  lon2 <- lon2 * pi / 180
  
  # Approx radius of the Earth in meters
  earth_radius <- 6371000
  
  # Haversine formula
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  distance <- earth_radius * c
  
  return(distance)
}

###################### Create a Route Between Two Points #######################
# A route plotted "as the sturgeon swims" between two points. This is used to 
# connect the destination or sturgeon location to the end of a linestring to 
# complete a sturgeon route.

asthesturswims <- function(start, end) {
  line <- st_linestring(rbind(
    st_coordinates(start),
    st_coordinates(end)
  ))
  return(line)
}

########################## Total Lengths Function ##############################
# Used to calculate the total length of the routes created and used by the 
# sturgeon agents and then used in the daytravel function. The points argument
# needs to be of the format st_coordinates() produces.

# Reminder lat = Y and long = X!!

totallength <- function(points){
  # Establish the length as 0 so it can be added to
  length <- 0
  num <- 1:(nrow(points)-1)
  # Calculate the distances between each of the indexes using the haversine
  # distance function
  for (i in num){
    # Index 2 is Y and 1 is X for each of these atomic vectors
    lat1 <- points[i, "Y"]
    lon1 <- points[i, "X"]
    lat2 <- points[(i+1), "Y"]
    lon2 <- points[(i+1), "X"]
    distance <- haversine_distance(lat1, lon1, lat2, lon2)
    length <- length + distance
  }
  return(length)
}

######################### Sturgeon Movement Function ###########################
# This is the function by which agents move around and uses four arguments.
# The distance argument is the sturgeon speed, the movedata argument is the row
# of the moveInfo data frame for that sturgeon, dest is the destination as an sf
# point and the direction is given as up for migration upriver and down for the 
# return

daytravel <- function(distance, movedata, dest, direction) {
  # Extract elements from the moveInfo row to use later
  route <- movedata$geometry
  route_coords <- st_coordinates(route)
  daystrav <- movedata$Days_travelled
  total_length <- movedata$Total_Route
  
  # Calculate the distance covered at this point by multipling the days already
  # travelled by the speed and then add another day's speed (m/day)
  distcovered <- (daystrav * distance) + distance
  # If the total route length is smaller or equal to the distance covered, then
  # the sturgeon has reached its destination
  if (total_length <= distcovered){
    new_coord <- dest
    # If the sturgeon is still travelling...
  } else {
    # If they are travelling upstream calculate the distances between each
    # of the points in the route
    if (direction == "up"){
      for (i in 1:(nrow(route_coords)-1)){
        lat1 <- route_coords[[i, "Y"]]
        lon1 <- route_coords[[i, "X"]]
        lat2 <- route_coords[[(i + 1), "Y"]]
        lon2 <- route_coords[[(i + 1), "X"]]
        difference <- haversine_distance(lat1, lon1, lat2, lon2)
        # Then take away this calculated difference between the two indices from
        # the distance to be travelled until it reaches 0 or difference is
        # greater than the distance yet to be travelled.
        if (distcovered >= difference){
          distcovered <- distcovered - difference
          # In the unlikely event that it exactly ends at one of the indices,
          # then set that as the new coord.
          if (distcovered == 0){
            new_coord <- st_point(route_coords[i + 1, 1:2])
            break # Exit the for loop
          }
          # If the sturgeon won't make it to the next index, calculate the new
          # location between the two points
        } else {
          # This bearing function from the TrackReconstruction package gives radians
          bearing <- CalcBearing(lat1, lon1, lat2, lon2)
          # Approximate radius of the Earth in meters
          earth_radius <- 6371000
          # convert degrees to radians
          lat1 <- lat1 * pi / 180
          lon1 <- lon1 * pi / 180
          # Calculate the new latitude and longitudes
          newlat <- asin(sin(lat1) * cos(distcovered / earth_radius) + cos(lat1) * sin(distcovered / earth_radius) * cos(bearing))
          newlon <- lon1 + atan2(sin(bearing) * sin(distcovered / earth_radius) * cos(lat1), cos(distcovered / earth_radius) - sin(lat1) * sin(newlat))
          # convert radians back to degrees
          newlat <- newlat * 180 / pi
          newlon <- newlon * 180 / pi
          # convert to an sf object
          new_coord <- st_point(c(newlon, newlat))
          # No more distance to cover in this day
          distcovered <- 0
          break # Exit the for loop
        }
      }
    }
    if (direction == "down"){
      # If the sturgeon is travelling downstream, then do a similar process with
      # a few notable exceptions. 
      # The for loop now counts backwards through the coordinate indices.
      for (i in (nrow(route_coords)):2){
        lat1 <- route_coords[[i, "Y"]]
        lon1 <- route_coords[[i, "X"]]
        lat2 <- route_coords[[(i - 1), "Y"]]
        lon2 <- route_coords[[(i - 1), "X"]]
        difference <- haversine_distance(lat1, lon1, lat2, lon2)
        if (distcovered > difference){
          distcovered <- distcovered - difference
          if (distcovered == 0){
            # Changed to minus 1 for the next index downstream
            new_coord <- route_coords[i - 1, 1:2]
            new_coord <- st_point(new_coord)
            break # Exit the for loop
          }
        } else {
          bearing <- CalcBearing(lat1, lon1, lat2, lon2)
          earth_radius <- 6371000
          lat1 <- lat1 * pi / 180
          lon1 <- lon1 * pi / 180
          newlat <- asin(sin(lat1) * cos(distcovered / earth_radius) + cos(lat1) * sin(distcovered / earth_radius) * cos(bearing))
          newlon <- lon1 + atan2(sin(bearing) * sin(distcovered / earth_radius) * cos(lat1), cos(distcovered / earth_radius) - sin(lat1) * sin(newlat))
          newlat <- newlat * 180 / pi
          newlon <- newlon * 180 / pi
          new_coord <- st_point(c(newlon, newlat))
          distcovered <- 0
          break # Exit the for loop
        }
      } 
    }
  }
  # Return the new coordinate
  return(new_coord) 
}

############################### AGENT CREATION #################################

######################## Convert YoY to Subadults ##############################
# When the YoY reach a year old, they will be converted to subadults and added
# to the sturgeon_sf data frame. The input is the row in the younggroups data
# frame.

convtosub <- function(YoYgroup){
  # Where they spawn is handed down to them from their mother since there is 
  # some evidence to suggest they express site fidelity
  spawning <- as.numeric(YoYgroup$Spawnsite)
  sub <- data.frame(ID = 1:YoYgroup$Number,
                    Class = rep("subadult", YoYgroup$Number),
                    Sex = NA,
                    Age = rep(364, YoYgroup$Number), # Remove one day as they will age a day as a subadult
                    Strategy = rep(YoYgroup$Strategy, YoYgroup$Number),
                    State = rep("Feeding", YoYgroup$Number),
                    Spawnsite = rep(spawning, YoYgroup$Number),
                    AdultOrDead = NA,
                    X_Coord = NA,
                    Y_Coord = NA)
  # The same probability of sex is 50:50 and they are all assigned locations in
  # the Black Sea.
  for (s in sub$ID){
    sex <- sample(c("M", "F"), size = 1, prob = c(0.5, 0.5))
    if (sex == "F"){
      AorD <- round(rnorm(n=1, mean = 5110, sd = 183))
    } else {
      AorD <- round(rnorm(n=1, mean = 4198, sd = 274))
    }
    coord <- sample(feed$OID, size = 1)
    x <- feed$POINT_X[coord]
    y <- feed$POINT_Y[coord]
    # enter the values into the dataframe
    sub$Sex[s] <- sex
    sub$AdultOrDead[s] <- AorD
    sub$X_Coord[s] <- x
    sub$Y_Coord[s] <- y
  }
  # Find the last ID number in sturgeon_sf data frame in order to add on the
  # new subadult IDs
  lastID <- tail(sturgeon_sf$ID, 1)
  sub$ID <- sub$ID + lastID
  sub <-
    st_as_sf(
      sub,
      coords = c("X_Coord", "Y_Coord"),
      crs = 4326,
      na.fail = FALSE
    )
  return(sub)
}

############################## Generate Eggs ###################################
# This function is to create the collective agents for the eggs. The
# input argument is the line from sturgeon_sf for the mother of the eggs. 

# The average stated in the literature for the number of eggs produced is 574400
# and the range is 228400 to 964800. A normal distribution was created to
# approximate this distribution with 574400 as the mean and 90000 as the SD
# which produces a range approximately within the one quoted in the literature.

geneggs <- function(mother){
  clutch_size <- rnorm(n=1, mean = 574400, sd = 90000)
  clutch_size <- round(clutch_size, digits = 0)
  # Lay out in the order of the data frame
  eggs <- data.frame(ID = NA, Class = NA, Number = NA, Age = NA, Strategy = NA, State = NA, Spawnsite = NA, Leaving_feed = NA)
  new_row_values <- c(mother$ID, "eggs", clutch_size, 0, mother$Strategy, "Incubating", as.numeric(mother$Spawnsite), Leaving_feed = NA)
  eggs[1,] <- new_row_values
  eggs <- st_sf(eggs, geometry = st_sfc(mother$geometry))
  eggs$Age <- as.numeric(eggs$Age)
  eggs$Number <- as.numeric(eggs$Number)
  
  return(eggs)
}

############################ ASSORTED EXTRAS! ##################################

############################ NA and Logic Function #############################
# I couldn't find any existing function for this so I just made one. This is to
# fix the instance where the logic returns an NA which I want to mean a FALSE
# but otherwise the logic should return the correct answer.
nareplace <- function(logic){
  if (is.na(logic)){
    end <- FALSE
  } else {
    end <- logic
  }
  return(end)
}

############################### Ready to Spawn? ################################
# This function is used to determine whether the sturgeon are going to migrate
# this year depending on when they last migrated. The input arguments are the
# current day in the model and the row from sturgeon_sf for the sturgeon.

readytospawn <- function(adultstur, date, spawntrack){
  # First check to see if they've spawned in the model before
  if (adultstur$ID %in% spawntrack$SturID){
    # Then take the last time they spawned
    lastspawn <- max(spawntrack$Date_left[spawntrack$SturID == adultstur$ID])
    lastspawn <- difftime(date, lastspawn, units = "days")
    # Check if they spawned over a year and a half ago (not the full two because
    # there is some variation in the months they spawn in).
    # If they spawned last year then they don't spawn this year
    if (lastspawn < 548) {
      doispawn <- FALSE
    }
    # Sturgeon will either spawn in year 2 or 3 after they last spawned so
    # they have a 50% chance of spawning in year 2
    if (lastspawn > 548){
      doispawn <- runif(1) < 0.5
    } 
    # If it is now year 3 after they spawned, then they will certainly spawn
    if (lastspawn > 913){
      doispawn <- TRUE
    }
  }
  # If they haven't spawned yet, then they will since it is their first year.
  else {
    doispawn <- TRUE
  }
  return(doispawn)
}

############################ Detect Water Peak #################################
# I have moved this to the ABM script page to avoid confusion if running multiple 
# experiment runs at the same time testing the different water peak methods.

###################### Calculate the egg incubation time #######################
# Based on a graph from Ivanov, 1987.
incubationtime <- function(x){
  return(1432.81902952 * (sqrt(202/280))^x)
}

############################# Mortality Functions ##############################
## YoY Mortality ##
# The average annual mortality for the YoY was determined to be 0.9996 to 1.0 in
# the Joric et al., 2009 paper. However, this is a very large range. Working 
# with the starting average egg of 574400 and an ending number of YoY at 100,
# the daily survival rate is established

survivingYoY <- function(number, endingnum){
  survivalrate <- (endingnum/574400)^(1/365)
  living <- number*survivalrate
  return(round(living, digits = 0))
}

################# Create Functions for Adult Mortality #########################
# This death rate reflects an average age at death of 24 and an exponential rate
# using the information on maturation and the maximum age of sturgeon

avmaturefem <- 14*365
avmaturemal <- 11.5*365
maxage <- 100*365

setAorD <- function(sex, matureage, meandeathage){
  if (sex == "F"){
    femrate <- 1 / (meandeathage - avmaturefem)
    f_ageatdeath <- matureage + rexp(1, rate = femrate)
    f_ageatdeath <- pmin(f_ageatdeath, maxage)
    AorD <- f_ageatdeath
  } else {
    malrate <- 1 / (meandeathage - avmaturemal)
    m_ageatdeath <- matureage + rexp(1, rate = malrate)
    m_ageatdeath <- pmin(m_ageatdeath, maxage)
    AorD <- m_ageatdeath
  }
  return(round(AorD))
}

