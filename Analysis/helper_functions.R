
# Import the data from the model output -----------------------------------

# * Get the simulated data in a exploitable format ------------------------

make_list_from_output <- function (file_name){
  data <- read.table(file_name, sep = " ")
  my_list <- list()
  for(i in 1:length(data[,2])){
    if(data[i,2]=="sim"){
      if(data[i,3] == "end"){
        break
      }
      for(j in 1:(length(data[,2])-i)){
        if(data[i+j,2]=="sim"){
          x <- as.data.frame(data[(i+2):(i+j-1),2:6])
          x[,2] <- as.numeric(x[,2])
          x[,3] <- as.numeric(x[,3])
          x[,4] <- as.numeric(x[,4])
          x[,5] <- as.numeric(x[,5])
          my_list[[length(my_list)+1]] <- x
          colnames(my_list[[length(my_list)]]) <- data[2,2:6]
          rownames(my_list[[length(my_list)]]) <- NULL
          break
        }
      }
    }
  }
  return(my_list)
}


# Calculations from the outputs -------------------------------------------

# * Individuals measures --------------------------------------------------

# * * Get the latency, rank, day means per individual ---------------------

get_means <- function(list_name){
  for(i in 1:length(list_name)){
    if(i == 1){
      z <- list_name[[1]]
    }
    else{
      z <- rbind(z, list_name[[i]])
    }
  }
res_latency <- tapply(z$total_time, z$ind_id, mean)
res_rank <- tapply(z$rank, z$ind_id, mean)
res_day <- tapply(z$day, z$ind_id, mean)
my_list <- list(latency = res_latency, rank = res_rank, day = res_day)
return(my_list)
}

# * * Get the individual probability to get infected ----------------------

get_proba_inf <- function(list_name, table_name){
  x <- data.frame(matrix(ncol = 1, nrow = nrow(table_name)))
  colnames(x) <- c("proba_infection")
  for(i in 1:nrow(table_name)){
    x2 <- 0
    for(j in 1:length(list_name)){
      for(k in 1:nrow(list_name[[j]])){
        if(list_name[[j]][,1][k] == table_name$id_ind[i]){
          x2 <- x2 + 1
        }
      }
    }
    x[i,1] <- x2 / length(list_name)
  }
  return(x)
}

# * * Get the mean transmission events duration, mean proportion o --------

get_outbreak_dyn <- function(list_name, table_name, threshold){
  v <- data.frame(matrix(ncol = 4, nrow = nrow(table_name)))
  colnames(v) <- c("duration", "outbreak_size", "duration_complete_transmission", "duration_incomplete_transmission")
  for(i in 1:nrow(table_name)){
    v2 <- 0
    v3 <- 0
    v4 <- 0
    v5 <- 0
    for(j in 1:length(list_name)){
      if(list_name[[j]][1,1] == table_name$id_ind[i]){
        if(v2[1] == 0){
          v2[1] <- max(list_name[[j]]$total_time)
        }
        else{
          v2 <- c(v2, max(list_name[[j]]$total_time))
        }
        if(v3[1] == 0){
          v3[1] <- max(list_name[[j]]$rank)
        }
        else{
          v3 <- c(v3, max(list_name[[j]]$rank))
        }
        if(nrow(list_name[[j]]) >= threshold * nrow(table_name)){
          if(v4[1] == 0){
            v4[1] <- max(list_name[[j]]$total_time)
          }
          else{
            v4 <- c(v4, max(list_name[[j]]$total_time))
          }
        }
        else{
          if(v5[1] == 0){
            v5[1] <- max(list_name[[j]]$total_time)
          }
          else{
            v5 <- c(v5, max(list_name[[j]]$total_time))
          }
        }
      }
    }
    v[i,1] <- mean(v2)
    v[i,2] <- mean(v3) / nrow(table_name)
    if(mean(v4) == 0){
      v[i,3] <- NA
    }
    else{
      v[i,3] <- mean(v4)
    }
    if(mean(v5) == 0){
      v[i,4] <- NA
    }
    else{
      v[i,4] <- mean(v5)
    }
  }
  return(v)
}

# * * Individual probability to be the first infected via sociality (2nd infected individual) -------

get_proba_inf_1st <- function(list_name, table_name){
  x <- data.frame(matrix(ncol = 1, nrow = nrow(table_name)))
  colnames(x) <- c("proba_1st_infection")
  for(i in 1:nrow(table_name)){
    x2 <- 0
    for(j in 1:length(list_name)){
      if(nrow(list_name[[j]]) > 1){
        if(list_name[[j]][2,1] == table_name$id_ind[i]){
          x2 <- x2 + 1
        }
      }
    }
    x[i,1] <- x2 / length(list_name)
  }
  return(x)
}

# * * Outbreak proportion depending on the first infected individu --------

get_outbreak_proportion_time_constrained <- function(list_name, table_name, nb_days){
  v <- data.frame(matrix(ncol = 1, nrow = nrow(table_name)))
  for(i in 1:nrow(table_name)){
    v2 <- 0
    for(j in 1:length(list_name)){
      if(list_name[[j]][1,1] == table_name$id_ind[i]){
        if(max(list_name[[j]][["day"]]) > nb_days){
          for(k in 1:nrow(list_name[[j]])){
            if(list_name[[j]][["day"]][k] > nb_days){
              if(v2[1] == 0){
                v2[1] <- list_name[[j]][["rank"]][k-1]
                break
              }
              else{
                v2 <- c(v2, list_name[[j]][["rank"]][k-1])
                break
              }
            }
          }
        }
        else{
          if(v2[1] == 0){
            v2[1] <- max(list_name[[j]][["rank"]])
          }
          else{
            v2 <- c(v2, max(list_name[[j]][["rank"]]))
          }
        }
      }
    }
    v[i,1] <- mean(v2) / nrow(table_name)
  }
  return(v)
}

# * * Individual probability to get infected during a given number --------

get_proba_inf_time_constrained <- function(list_name, table_name, nb_days){
  x <- data.frame(matrix(ncol = 1, nrow = nrow(table_name)))
  colnames(x) <- c("proba_infection")
  for(i in 1:nrow(table_name)){
    x2 <- 0
    for(j in 1:length(list_name)){
      for(k in 1:nrow(list_name[[j]])){
        if(list_name[[j]][["day"]][k] > nb_days){
          break
        }
        else{
          if(list_name[[j]][,1][k] == table_name$id_ind[i]){
            x2 <- x2 + 1
          }
        }
      }
    }
    x[i,1] <- x2 / length(list_name)
  }
  return(x)
}

# * * Get the mean latency between two consecutive infections depe --------

get_consecutive_inf_latency <- function(list_name, table_name){
  v <- data.frame(matrix(ncol = 1, nrow = nrow(table_name)))
  colnames(v) <- "consecutive_infection_latency"
  for(i in 1:nrow(table_name)){
    v2 <- 0
    for(j in 1:length(list_name)){
      if(list_name[[j]][1,1] == table_name$id_ind[i]){
        if(v2[1] == 0){
          v2[1] <- mean(list_name[[j]]$latency)
        }
        else{
          v2 <- c(v2, mean(list_name[[j]]$latency))
        }
      }
    }
    v[i,1] <- mean(v2)
  }
  return(v)
}

# * Global measures -------------------------------------------------------

# * * Mean total number of individuals infected ---------------------------

nb_inf <- function (list_name){
  w <- rep(0, length(list_name))
  for (i in 1:length(list_name)){
    w[i] <- nrow(list_name[[i]])
  }
  return(w)
}

# * * Mean duration of simulation -----------------------------------------

sim_time <- function(list_name){
  y <- rep(0, length(list_name))
  for(i in 1:length(list_name)){
    y[i] <- max(list_name[[i]]$total_time)
  }
  return(y)
}

# * * Mean duration to get all individuals infected -----------------------

all_inf_time <- function(list_name, trait_table){
  y2 <- 0
  for(i in 1:length(list_name)){
    if(length(list_name[[i]]$total_time) == nrow(trait_table)){
      if(y2[1] == 0){
        y2[1] <- max(list_name[[i]]$total_time)
      }
      else{
        y2 <- c(y2, max(list_name[[i]]$total_time))
      }
    }
  }
  return(y2)
}