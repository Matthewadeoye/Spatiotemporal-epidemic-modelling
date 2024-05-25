#Design matrices

#model 0
DesignMatrixModel0<- function(y, adjacencymatrix){
  ndept<- nrow(y)
  time<- ncol(y)
z_it<- matrix(0, ndept, time)
z_it2<- matrix(0, ndept, time)
return(list(z_it, z_it2))
}

#Model 1
DesignMatrixModel1<- function(y, adjacencymatrix){
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)
  for(i in 1:ndept){
    for(t in 2:time){
      if(y[i, t-1] > 0)
        z_it[i, t]<- 1
    }
  }
  z_it2<- matrix(0, ndept, time)
  return(list(z_it, z_it2))
}

#Model 2
DesignMatrixModel2<- function(y, adjacencymatrix){
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)
  
  for(i in 1:ndept){
    indexes = which(adjacencymatrix[i, ] > 0 & 1:ndept != i)
    
    z_it[i, 1] = 0
    for(t in 2:time){
      for(b in 1:length(indexes)){
        neighbor_index <- indexes[b]
        if(y[i, t-1] > 0 || y[neighbor_index, t-1] > 0){
          z_it[i, t] = 1
      }else{
         z_it[i, t] = 0
       }
      }
    }
  }
  z_it2<- matrix(0, ndept, time)
  return(list(z_it, z_it2))
}

#Model 3
DesignMatrixModel3<- function(y, adjacencymatrix){
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)
  z_it2<- matrix(0, ndept, time)
  
  for(i in 1:ndept){
    z_it[i, 1] = 0
    for(t in 2:time){
      if(y[i, t-1] > 0){
        z_it[i, t]<- 1
      }else{
        z_it[i, t]<- 0
      }
    }
  }
  
  for(i in 1:ndept){
    indexes = which(adjacencymatrix[i, ] > 0 & 1:ndept != i)
    
    z_it2[i, 1] = 0
    for(t in 2:time){
      flag = 0
      for(b in 1:length(indexes)){
        neighbor_index <- indexes[b]
        if(y[neighbor_index, t-1] > 0){
          flag = 1
          break;
        }
      }
      z_it2[i, t] = flag
    }
  }
  return(list(z_it, z_it2))
}

#model 4
DesignMatrixModel4<- function(y, adjacencymatrix){
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)
  for(i in 1:ndept){
    for(t in 2:time){
      z_it[i, t]<- log(y[i, t-1] + 1)
    }
  }
  z_it2<- matrix(0, ndept, time)
  return(list(z_it, z_it2))
}

#model 5
DesignMatrixModel5<- function(y, adjacencymatrix){
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)
  for(i in 1:ndept){
    indexes <- which(adjacencymatrix[i, ] > 0 & 1:ndept != i)
    
    for(t in 2:time){
      sum_neighbors <- 0
      for(b in indexes){
       sum_neighbors<- sum_neighbors + y[b, t-1]  
      }
        z_it[i, t]<- sum_neighbors
    }
  }
  for(i in 1:ndept){
    for(t in 2:time){
      z_it[i, t] = log(y[i, t-1] + z_it[i, t] + 1)
    }
  }
  z_it2<- matrix(0, ndept, time)
  return(list(z_it, z_it2))
}

#model 6
DesignMatrixModel6 <- function(y, adjacencymatrix) {
  ndept <- nrow(y)
  time <- ncol(y)
  z_it <- matrix(0, ndept, time)
  z_it2 <- matrix(0, ndept, time)
  
  # Compute z_it
  for (i in 1:ndept) {
    for (t in 2:time) {
      z_it[i, t] <- log(y[i, t-1] + 1)
    }
  }
  
  # Compute z_it2
  for(i in 1:ndept){
    indexes <- which(adjacencymatrix[i, ] > 0 & 1:ndept != i)
    
    for(t in 2:time){
      sum_neighbors <- 0
      for(b in indexes){
        sum_neighbors<- sum_neighbors + y[b, t-1]  
      }
      z_it2[i, t]<- log(sum_neighbors + 1)
    }
  }
  return(list(z_it, z_it2))
}