
library(proxy)
library(tidyr)
library(celltrackR)
library(ggplot2)
library(celltrackR)

################################################################################

xmax<-600 # Maximum x coordinate
ymax<-600 # Maximum y coordinate
xmin<--600 # Minimum x coordinate
ymin<--600 # Minimum y coordinate

################################################################################

pred_characteristics <- c("location_x", "location_y", "status") # Names of state variables of agents (foragers)
pdtr.character <- length(pred_characteristics) 

predpop <- array(NA, dim=c(1, pdtr.character)) 
colnames(predpop) <- c(pred_characteristics) # Array with one row and column names corresponding to state variables

#####

predpop[, "location_x"] <- 0 # Initial forager x location set to zero
predpop[, "location_y"] <- 0 # Initial forager y location set to zero


################################################################################
# Functions: 
################################################################################

# 1. Directed foraging movement function (Move.Pred)

Move.Pred<-function(coordinates, xloc.prey, yloc.prey,SLf.dist,TAf.dist,SLt.dist,TAt.dist){
  
  x <-as.numeric(coordinates[,1]) # X location of agent
  y <-as.numeric(coordinates[,2]) # Y location of agent
  
  a.sq<-(x-xloc.prey)^2 
  b.sq<-(y-yloc.prey)^2
  c<-sqrt(a.sq+b.sq) # Pairwise agent-patch distances
  dists_df<-data.frame(c) 
  foraging<-which(dists_df<preypop[,"radius"]) # Foraging agents are those within patch radius
  
  status<-ifelse(length(foraging)>0, 1, 0) 
  pred.status<-status # Change status: 1 = foraging, 0 = travelling
  
  if(pred.status==1){ # Change movement mode to "foraging" mode if within prey patch
    Step.length<-sample(SLf.dist,1)
    Turn.angle<-sample(TAf.dist,1) # Draws step lengths and turning angles from "foraging" distribution
  }
  else{
    Step.length<-sample(SLt.dist,1)
    Turn.angle<-sample(TAt.dist,1)
  }
  
  if(Turn.angle < 0){ # Translate turning anlge to value between 0 and 2*pi
    Turn.angle<-abs(Turn.angle)+pi
  } else{}
  
  x.next<-cos(Turn.angle)*Step.length # Derive change in x location from selected TA and SL
  y.next<-sin(Turn.angle)*Step.length # Derive change in y location from selected TA and SL
  
  xloc.new<-x.next+x # Add x.next to current x location to produce new x location
  yloc.new<-y.next+y # Add y.next to current y location to produce new y location
  
  # Generate output for this function, containing updated x and y locations, SL, TA, and status
  move<-matrix(0,1,5)
  colnames(move)<-c("x","y","SL","TA","Status")
  move[,1]<-xloc.new
  move[,2]<-yloc.new
  move[,3]<-Step.length
  move[,4]<-Turn.angle
  move[,5]<-pred.status
  
  # Return the output generated above
  return(move)
}

# 2. Brownian motion function (Move.Brownian)

Move.Brownian <- function(xloc.prey,yloc.prey,ts){
  
  # Use 'simulateTracks' function in package 'cellTrackR' to automatically generate 2D Brownian movement trajectory of 100 timesteps
  brownian.tracks <- simulateTracks(1, brownianTrack(nsteps = ts, dim=2, mean=0, sd=20))
  brownian_df<-data.frame(brownian.tracks)
  brownian_df<-brownian_df[,(3:4)]
  
  # Plot trajectory:
  plot(brownian.tracks)
  
  # Return the output generated above
  return(brownian_df)
}

################################################################################
# 1. Directed foraging (DF) model:
################################################################################


SLf.dist<-rpois(c(1:1000), lambda=5) # Step length distribution for agents in DF simulation with status "foraging"
TAf.dist<-runif(1000,min=0, max=2*pi) # Turning angle distribution for agents in DF simulation with status "foraging"
hist(TAf.dist)

simulations<-1000 # Number of simulations
n.sim<-1:simulations 
ts<-100 # Number of timesteps in each simulation
time<-2:ts
pred.tracks<-array(data=0, dim=c(ts,5)) # Create array to store output of each simulation of DF model; 
# Number of rows of  array = number of timesetps (ts); 5 columns; all values set to zero
colnames(pred.tracks)<-c("x","y" ,"SL","TA","Status") # Set output column names
pred.tracks.df<-data.frame(pred.tracks) # Convert output array to dataframe
head(pred.tracks.df)

# Create output for all simulations in DF model (length = number of simulations, 10 columns)
movement.directed<-array(0, dim=c(simulations, 10)) 
colnames(movement.directed)<-c("Simulation no.","N patches","Mean patch dist","Min patch dist","Max patch dist","Total dist moved","Displacement","Mean TA","Mean SL","% ts foraging")
movement.directed.df<-data.frame(movement.directed)
movement.directed.df[,1]<-c(1:simulations)

for(s in n.sim){ # Loop the following through 1000 simulations: 
  
  n.prey<-sample(1:20,1) # Randomly select  number of resource patches between 1 and 20
  prey_characteristics <- c("location_x", "location_y", "radius","area") # State variables of patches 
  prey.chr <- length(prey_characteristics)
  
  preypop <- array(NA, dim=c(n.prey, prey.chr))
  colnames(preypop) <- c(prey_characteristics) # Array  ontaining patch state variables
  
  preypop[,"location_x"]<-sample(c(-500:500), n.prey, replace = TRUE) # Assign random x location within xmax and xmin - this is the x coordinate for the patch centre
  preypop[,"location_y"]<-sample(c(-500:500), n.prey, replace = TRUE) # Assign random y location within ymax and ymin - this is the y coordinate for the patch centre
  
  preypop[,"area"]<-as.integer(sample(runif(n=1000, min=2000, max=20000), n.prey)) # Assign random area from uniform distribution between 2000 and 20000 
  
  preypop[,"radius"]<-sqrt(as.numeric(preypop[,"area"])/pi) # Derive patch radius using area
  
  preypop.df<-data.frame(preypop)
  n.prey<-dim(preypop.df)[1] #Length of preypop = number of prey patches
  
  xloc.prey<-as.numeric(preypop.df[,1]) 
  yloc.prey<-as.numeric(preypop.df[,2])
  
  dc<-sqrt((xloc.prey^2) + (yloc.prey^2)) # Calculate the distance of each patch centre from the centre of the landscape (coordinate 0,0)
  
  plot(xloc.prey, yloc.prey, xlab="x coordinate", ylab="y coordinate", pch=c(""), xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  symbols(x=xloc.prey, y=yloc.prey,circles=preypop.df[,3], add=T, inches=F) # Plot patches as circles on landscape
  
  m<-sample(runif(n=1000,min=0, max=2*pi),1) # Generate a  heading from uniform distribution between 0 and 2*pi along which agent will initially move in simulation
  
for(t in time){ # Run this loop for 100 timesteps (the following programs the movement of the forager in the landscape generated previously)
  SLt.dist<-rnorm(1000,50,5) # Generate step length distribution for 'travelling' agents
  hist(SLt.dist)
  TAt.dist1<-rnorm(1000, mean=m, sd=2) # Draw turning angle from 'travelling' distribution
  crop1<-unique(which(TAt.dist1<0)) 
  crop2<-unique(which(TAt.dist1>2*pi)) # Crop1 and Crop2 contain values of TAt.dist1 which fall outside of range 0-2*pi (i.e. 360 degrees)
  TAt.dist<-TAt.dist1[-c(crop1, crop2)] # Remove these values
  
  # Predator takes a step of given step length and turning angle (sampled from SLt.dist and TAt.dist if travelling or SLf.dist and TAf.dist if foraging)
  # using Move.Pred function
  step<-Move.Pred(pred.tracks.df[c(t-1),1:2], xloc.prey, yloc.prey,SLf.dist,TAf.dist,SLt.dist,TAt.dist)
 
  #The following three 'if' loops reorientate the trajectory of the agent if its location exceeds the x and y boundary limits
  # i.e. they deflect agents away from boundary
   if(step[,"x"]>xmax){
    step[,"x"]<-xmax
    m<-sample(runif(n=1000,min=0, max=2*pi),1)
  }
  else{}
  if(step[,"x"]<xmin){
    step[,"x"]<-xmin
    m<-sample(runif(n=1000,min=0, max=2*pi),1)
  }
  else{}
  if(step[,"y"]>ymax){
    step[,"y"]<-ymax
    m<-sample(runif(n=1000,min=0, max=2*pi),1)
  }
  else{}
  if(step[,"y"]<ymin){
    step[,"y"]<-ymin
    m<-sample(runif(n=1000,min=0, max=2*pi),1)
  }
  else{}
  
  pred.tracks.df[t,]<-step # Assign new xy location and foraging status to corresponding timesetp of output "pred.tracks.df" generated earlier
  
} 
  # Plot agent trajectories as defined by pred.tracks.df
  plot(pred.tracks.df$x,pred.tracks.df$y, pch=c(""),col='blue', xlab="x coordinate", ylab="y coordinate", xlim=c(xmin,xmax), ylim=c(ymin,ymax), main = "B")
  points(pred.tracks.df$x[1], pred.tracks.df$y[1], pch=3, col='red', cex=2) # Red cross to mark start point
  points(pred.tracks.df$x[100], pred.tracks.df$y[100], pch=3, col='blue', cex=2) # Blue cross to mark end point
  lines(pred.tracks.df$x,pred.tracks.df$y) # Plot trajectory
  points.prey<-points(xloc.prey, yloc.prey, pch="")
  symbols(x=xloc.prey, y=yloc.prey,circles=preypop.df[,3], add=T, inches=F) # Overlay patches on movement plot
  
  # Assign data  from each simulation to movement.directed.df output generated earlier: this is the overall output for this model
  movement.directed.df[s,"N.patches"]<-dim(preypop)[1]
  movement.directed.df[s,"Mean.patch.dist"]<-mean(dc)
  movement.directed.df[s,"Min.patch.dist"]<-min(dc)
  movement.directed.df[s,"Max.patch.dist"]<-max(dc)
  movement.directed.df[s,"Total.dist.moved"]<-sum(pred.tracks.df[,"SL"])
  movement.directed.df[s,"Displacement"]<-sqrt((pred.tracks.df[ts,"x"]^2)+(pred.tracks.df[ts,"y"]^2))
  movement.directed.df[s,"Mean.TA"]<-mean(pred.tracks.df[,"TA"])
  movement.directed.df[s,"Mean.SL"]<-mean(pred.tracks.df[,"SL"])
  movement.directed.df[s,"X..ts.foraging"]<-sum(pred.tracks.df[,"Status"])
  
}

# Save data to CSV file (insert file pathway)
write.csv(movement.directed.df,"-insert file pathway here -",row.names = TRUE)


###############################################################################
# 2.Brownian motion (BM) model
###############################################################################


# As before, create output for each simulation:
pred.tracks<-array(data=0, dim=c(ts,5))
colnames(pred.tracks)<-c("x","y" ,"SL","TA","Status")
brownian.tracks.df<-data.frame(pred.tracks)
head(brownian.tracks.df)

# As before, create output for whole model:
movement.data<-array(0, dim=c(simulations, 10))
colnames(movement.data)<-c("Simulation no.","N patches","Mean patch dist","Min patch dist","Max patch dist","Total dist moved","Displacement","Mean TA","Mean SL","% ts foraging")
movement.data.df<-data.frame(movement.data)
movement.data.df[,1]<-c(1:simulations)

for(s in n.sim){ 
  n.prey<-sample(1:20,1)
  prey_characteristics <- c("location_x", "location_y", "radius","area")
  prey.chr <- length(prey_characteristics)
  
  preypop <- array(NA, dim=c(n.prey, prey.chr))
  colnames(preypop) <- c(prey_characteristics)  
  
  preypop[,"location_x"]<-sample(c(-500:500), n.prey, replace = TRUE)
  preypop[,"location_y"]<-sample(c(-500:500), n.prey, replace = TRUE)
  
  preypop[,"area"]<-as.integer(sample(runif(n=1000, min=2000, max=20000), n.prey))
  
  preypop[,"radius"]<-sqrt(as.numeric(preypop[,"area"])/pi)
  
  preypop.df<-data.frame(preypop)
  n.prey<-dim(preypop.df)[1]
  
  xloc.prey<-as.numeric(preypop.df[,1]) 
  yloc.prey<-as.numeric(preypop.df[,2])
  
  dc<-matrix(sqrt((xloc.prey^2) + (yloc.prey^2)))
  n.dc<-dim(dc)[1]
  
  for(d in 1:n.dc){ # The following loop assigns an initial foraging status to agents (1 = foraging, 0 = travelling) 
                    # based on whether initial location (0,0) overlaps with prey patch
    if(dc[d,]<preypop[d,"radius"]){
      brownian.tracks.df[1,5]<-1
     }
    else{
      brownian.tracks.df[1,5]<-0
    }
  }
  
  plot(xloc.prey, yloc.prey, xlab="x coordinate", ylab="y coordinate", pch=c(""), xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  symbols(x=xloc.prey, y=yloc.prey,circles=preypop.df[,3], add=T, inches=F) # Plot patches 
  
  # Simulate Brownian motion for 100 timesteps using Move.Brownian function:
  step.b<-Move.Brownian(xloc.prey, yloc.prey,ts)
  step.df<-data.frame(step.b[-101,])
  brownian.tracks.df[,(1:2)]<-step.df
  brownian.length<-dim(brownian.tracks.df)[1]
  
  
  for(i in 2:brownian.length){
    x<-brownian.tracks.df[i,1]-brownian.tracks.df[(i-1),1]
    y<-brownian.tracks.df[i,2]-brownian.tracks.df[(i-1),2]
    SL<-abs(sqrt((x^2)+(y^2))) # Calculate step length of timestep t based on difference in x and y location between timesteps t and t-1
    brownian.tracks.df[i,3]<-SL # Assign step length to output
    
    brownian.tracks.df[i,4]<-asin(x/SL) # Calculate the turning angle between t and t-1
    
    a.prey<-abs(brownian.tracks.df[i,1]-xloc.prey)
    b.prey<-abs(brownian.tracks.df[i,2]-yloc.prey) 
    c.prey.df<-data.frame(sqrt((a.prey^2)+(b.prey^2))) # Pairwise distances between agent and each patch centre
    foraging<-which(c.prey.df[,1]<preypop[,3]) # Is the forager within a patch (i.e. is distance from patch centre less than patch radius?) 
    n.forage<-unique(foraging)
    
    if(length(n.forage)>0){ # Change foraging status depending on whether agent is within prey patch or not (1=foraging, 0 = travelling) 
      brownian.tracks.df[i,5]<-1
    } else{
      brownian.tracks.df[i,5]<-0
    }
  }
  
  # As before, plot agent trajectories and overlay patches
  plot(brownian.tracks.df$x,brownian.tracks.df$y, pch=c(""),col='blue', xlab="x coordinate", ylab="y coordinate", xlim=c(xmin,xmax), ylim=c(ymin,ymax), main="A")
  points(brownian.tracks.df$x[1], brownian.tracks.df$y[1], pch=3, col='red', cex=2)
  points(brownian.tracks.df$x[100], brownian.tracks.df$y[100], pch=3, col='blue', cex=2)
  lines(brownian.tracks.df$x,brownian.tracks.df$y)
  points.prey<-points(xloc.prey, yloc.prey, pch="")
  symbols(x=xloc.prey, y=yloc.prey,circles=preypop.df[,3], add=T, inches=F)
  
  # Assign data  from each simulation to movement.directed.df output generated earlier
  movement.data.df[s,"N.patches"]<-dim(preypop)[1]
  movement.data.df[s,"Mean.patch.dist"]<-mean(dc)
  movement.data.df[s,"Min.patch.dist"]<-min(dc)
  movement.data.df[s,"Max.patch.dist"]<-max(dc)
  movement.data.df[s,"Total.dist.moved"]<-sum(brownian.tracks.df[,"SL"])
  movement.data.df[s,"Displacement"]<-sqrt((brownian.tracks.df[ts,"x"]^2)+(brownian.tracks.df[ts,"y"]^2))
  movement.data.df[s,"Mean.TA"]<-mean(brownian.tracks.df[,"TA"])
  movement.data.df[s,"Mean.SL"]<-mean(brownian.tracks.df[,"SL"])
  movement.data.df[s,"X..ts.foraging"]<-sum(brownian.tracks.df[,"Status"])
  
}
 # Save movement.data.df to CSV
write.csv(movement.data.df,"-insert file pathway here -",row.names = TRUE)

