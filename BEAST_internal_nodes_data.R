setwd("~/tanager_spec")
load("tanager_spec/BEAST.trees.subset.RData")
library(splines2)
library(ggplot2)
library(tidyverse)
library(tidylog)
library(pavo)
library(photobiology)
library(ggspectra)

bspline.matrix <- bSpline(300:700, df = 38)


# Extracting root node ASR coefficients for 100 trees
root.reflectance.ASR <- NULL
for (i in 1:dim(root.spline.coef)[1]){
  
  x <- i
  
  root.coef <- root.spline.coef[x,]
  
  root.refl <- bspline.matrix %*% matrix(root.coef[2:39])
  
  root.refl <- root.refl + root.coef[1]
  
  root.reflectance.ASR <- cbind(root.reflectance.ASR, root.refl)
  
}


matplot(x = 300:700, y = root.reflectance.ASR[,-1], lwd = 2, type = "l", xlab = "Wavelength (nm)", ylab = "Reflectance (%)", main = "Root ASR BEAST")

root.reflectance.ASR <- data.frame(wl = 300:700, root.reflectance.ASR)

root.reflectance.ASR <- as.rspec(root.reflectance.ASR)

root.reflectance.ASR <- procspec(root.reflectance.ASR, fixneg = "zero")



matplot(x = root.reflectance.ASR$wl, y = root.reflectance.ASR[,-1], pch = 19, cex = 0.5)



# Repeating the same process but in the tidyverse

# Creates a column for run number and the intercept coefficient of each run for the root
y <- tibble(run = 1:100, root_intercept = root.spline.coef[,1])

# Splits the rows (non-intercept spline coefficients) for each run of the root
z <- split(root.spline.coef[,-1], f = 1:100, 38)

# Adds the spline coefficient lists to our tibble
y <- y %>% mutate(spline_coefs = z)

# Takes the spline coefficients our basis matrix and creates reflectance values (stored as lists) into our tibble
y <- y %>% group_by(run) %>% mutate(reflectance = list(bspline.matrix %*% matrix(unlist(spline_coefs), ncol = 1))) %>% ungroup()

# Rescales the reflectance values by adding the intercept coefficient to those values
y <- y %>% group_by(run) %>% mutate(rescaled_reflectance = list(unlist(reflectance) + root_intercept)) %>% ungroup()

# Adds a column containing the rspec dataframes to the tibble  
y <- y %>% group_by(run) %>% mutate(test = list(procspec(as.rspec(data.frame(wl = 300:700, reflectance = unlist(rescaled_reflectance))), fixneg = "zero")))

# Adds a column to the tibble that has the rgb color from the rspec object
y <- y %>% group_by(run) %>% mutate(color = spec2rgb(data.frame(test)))
  
# Creates a separate tibble that unlists the rescaled_reflectance (stored as a singular vector)
x <- y %>% unnest(cols = rescaled_reflectance)

# Manually 'procspec' our rescaled_reflectance
x$rescaled_reflectance[x$rescaled_reflectance < 0] <- 0

# Adds a wavelength column
x <- x %>% mutate(wl = 300:700)

# Plots each predicted reflectance curve from the BEAST runs
x %>% group_by(x$run) %>% ggplot(aes(x = wl, y = rescaled_reflectance)) + geom_point(alpha = 0.25, size = 0.75, color = as.character(rep(y$color, each = 401))) + theme() + theme_bw() + stat_wl_strip(ymin = -1.5, ymax = -0.5) + scale_fill_identity() + labs(x = "Wavelength (nm)", y = "Reflectance (%)", title = "Root Reflectance Predicted by BEAST")



# Applying what we did but with the internal nodes

a <- tibble(tree_name = rep(names(BEAST.trees.subset), each = length(53:101)))

a <- a %>% mutate(nodes_list = rep(c(1,2,3,5,8,9,10,11,12,14,17,20,21,24,26,29,31,33,34,35,38,40,44,47,49,51,54,55,56,58,60,64,65,67,68,69,70,72,76,77,81,83,85,88,89,90,92,94,97), 100), nodes_internal = rep(c(53:101), 100)) %>% mutate(join_list = 1:length(nodes_list))

extraction <- function(multiPhylo, internal_node_values){
  
  internal_node_coefs <- NULL
  
  for (run in names(multiPhylo)){
    
    print(run)
    
    spline_coefs <- multiPhylo[[run]]$annotations
    
    node_coefs <- NULL
    
    for (i in 1:length(spline_coefs)){
      
      node <- internal_node_values[i]
      
      print(node)
      
      node_coefs[i] <- spline_coefs[[node]]
      
    }
  
  internal_node_coefs <- c(internal_node_coefs, node_coefs)
  
  }
  return(tibble(internal_node_coefs))
}

r <- extraction(BEAST.trees.subset, internal_node_values = a$nodes_list[1:49])

r <- r %>% rowwise() %>% filter(class(internal_node_coefs) == "list") %>% ungroup() %>% mutate(nodes_list = rep(c(1,2,3,5,8,9,10,11,12,14,17,20,21,24,26,29,31,33,34,35,38,40,44,47,49,51,54,55,56,58,60,64,65,67,68,69,70,72,76,77,81,83,85,88,89,90,92,94,97), 100)) %>% mutate(join_list = 1:length(nodes_list))

a <- a %>% left_join(r, by = "join_list") %>% select(-c(nodes_list.y, join_list)) %>% rename(nodes_list = nodes_list.x)


# test case with a specific node

test.node <- a %>% filter(nodes_list == 1)

test.node <- test.node %>% rowwise() %>% mutate(intercept_coef = unlist(internal_node_coefs)[1]) %>% ungroup()

test.node <- test.node %>% rowwise() %>% mutate(spline_coefs = list(unlist(internal_node_coefs)[-1])) %>% ungroup()

test.node <- test.node %>% select(-internal_node_coefs)

test.node <- test.node %>% rowwise() %>% mutate(reflectance = list(bspline.matrix %*% matrix(unlist(spline_coefs), ncol = 1))) %>% ungroup()

test.node <- test.node %>% rowwise() %>% mutate(rescaled_reflectance = list(unlist(reflectance) + intercept_coef)) %>% ungroup()

test.node <- test.node %>% rowwise() %>% mutate(test = list(procspec(as.rspec(data.frame(wl = 300:700, reflectance = unlist(rescaled_reflectance))), fixneg = "zero")))

test.node <- test.node %>% rowwise() %>% mutate(color = spec2rgb(data.frame(test)))


# Creates a separate tibble that unlists the rescaled_reflectance (stored as a singular vector)
plot.test.node <- test.node %>% unnest(cols = rescaled_reflectance)

# Manually 'procspec' our rescaled_reflectance
plot.test.node$rescaled_reflectance[plot.test.node$rescaled_reflectance < 0] <- 0

# Adds a wavelength column
plot.test.node <- plot.test.node %>% mutate(wl = rep(300:700, 100))

# Plots each predicted reflectance curve from the BEAST runs
plot.test.node %>% rowwise() %>% ggplot(aes(x = wl, y = rescaled_reflectance)) + geom_point(alpha = 0.25, size = 0.75, pch = 19, color = as.character(rep(y$color, each = 401))) + theme() + theme_bw() + stat_wl_strip(ymin = -1.5, ymax = -0.5) + scale_fill_identity() + labs(x = "Wavelength (nm)", y = "Reflectance (%)", title = "Node 53 Reflectance Predicted by BEAST")

# function to plot all nodes

ASR_plotter <- function(dataframe, nodes){
  
  if (sum(names(dataframe) == "nodes_list") == 0){
    stop("Dataframe needs a column titled 'nodes_list'.")
  }
  
  if (sum(names(dataframe) == "internal_node_coefs") == 0){
    stop("Dataframe needs a column titled 'internal_node_coefs'.")
  }
  
  if (class(dataframe$internal_node_coefs) != "list"){
    stop("Spline coefficient column needs to be a list.")
  }
    
  if (length(nodes) == 1){
    
    
    test.node <- dataframe %>% filter(nodes_list == nodes)
    
    test.node <- test.node %>% rowwise() %>% mutate(intercept_coef = unlist(internal_node_coefs)[1]) %>% ungroup()
    
    test.node <- test.node %>% rowwise() %>% mutate(spline_coefs = list(unlist(internal_node_coefs)[-1])) %>% ungroup()
    
    test.node <- test.node %>% select(-internal_node_coefs)
    
    test.node <- test.node %>% rowwise() %>% mutate(reflectance = list(bspline.matrix %*% matrix(unlist(spline_coefs), ncol = 1))) %>% ungroup()
    
    test.node <- test.node %>% rowwise() %>% mutate(rescaled_reflectance = list(unlist(reflectance) + intercept_coef)) %>% ungroup()
    
    test.node <- test.node %>% rowwise() %>% mutate(test = list(procspec(as.rspec(data.frame(wl = 300:700, reflectance = unlist(rescaled_reflectance))), fixneg = "zero")))
    
    test.node <- test.node %>% rowwise() %>% mutate(color = spec2rgb(data.frame(test)))
    
    
    # Creates a separate tibble that unlists the rescaled_reflectance (stored as a singular vector)
    plot.test.node <- test.node %>% unnest(cols = rescaled_reflectance)
    
    # Manually 'procspec' our rescaled_reflectance
    plot.test.node$rescaled_reflectance[plot.test.node$rescaled_reflectance < 0] <- 0
    
    # Adds a wavelength column
    plot.test.node <- plot.test.node %>% mutate(wl = rep(300:700, 100))
    
    # Plots each predicted reflectance curve from the BEAST runs
    plot.test.node %>% rowwise() %>% ggplot(aes(x = wl, y = rescaled_reflectance)) + geom_point(alpha = 0.25, size = 0.75, pch = 19, color = as.character(rep(y$color, each = 401))) + theme() + theme_bw() + stat_wl_strip(ymin = -1.5, ymax = -0.5) + scale_fill_identity() + labs(x = "Wavelength (nm)", y = "Reflectance (%)", title = paste("Node", nodes, "Reflectance Predicted by BEAST"))
    
  }
  
  
}



# Spline coefficient densities at each node for 100 runs

test <- a %>% unnest(cols = internal_node_coefs)

test <- test %>% mutate(coef_name = rep(c("Int.", "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38"), 4900))

test2 <- test %>% filter(coef_name == "Int." & nodes_internal == 53) %>%  mutate(join_list = 1:length(coef_name))

test2 %>% ggplot(aes(x = unlist(internal_node_coefs))) + geom_density(fill = "blue")

test3 <- test %>% filter(coef_name == "Int." & nodes_internal == 54) %>% mutate(join_list = 1:length(coef_name))

test3 %>% ggplot(aes(x = unlist(internal_node_coefs))) + geom_density(fill = "red")

test4 <- test2 %>% left_join(test3, by = "join_list")

test4 %>% ggplot() + geom_density(aes(x = unlist(internal_node_coefs.x)), fill = "black", alpha = 0.5) + geom_density(aes(x = unlist(internal_node_coefs.y)), fill = "yellow2", alpha = 0.5) + guides(fill = guide_legend(title = "Node")) + theme() + scale_fill_manual(values = c("53", "54"))
  
test5 <- test4 %>% pivot_longer(cols = c(internal_node_coefs.x, internal_node_coefs.y))

test5 %>% ggplot(aes(x = unlist(value), fill = name)) + geom_density(alpha = 0.25) + scale_fill_manual(labels = c(53,54), name = "Node", values = c("black", "yellow2"))
#####

# Numerical values in annotation list corresponding to the internal nodes (not including tip labels) for our tanager tree
internal.node.labels <- c(1,2,3,5,8,9,10,11,12,14,17,20,21,24,26,29,31,33,34,35,38,40,44,47,49,51,54,55,56,58,60,64,65,67,68,69,70,72,76,77,81,83,85,88,89,90,92,94,97)

# Extracting all internal node ASR coefficients for 100 trees

internal.node.data <- NULL

for (i in 1:length(BEAST.trees.subset)){
  
  x <- i
  
  # Takes 1 of the subset of trees
  tree <- BEAST.trees.subset[[x]]
  
  internal.node.data.subset <- NULL

  
  # for loop extracting the spline coefficients for each internal node of our tree
    for (j in 1:length(internal.node.labels)){
  
      internal.node <- internal.node.labels[j]
      
      internal.coef <- tree$annotations[[internal.node]]
  
      internal.coef <- c(unlist(internal.coef$spline_coefficients_BEAST))
  
      internal.node.data.subset <- rbind(internal.node.data.subset, internal.coef)
    }
  
  internal.node.data <- rbind(internal.node.data, internal.node.data.subset)
  
}

# Creating a reflectance data set for each internal node

internal.node.reflectance <- NULL

for (i in 1:49){
  
  x <- i

  internal.node.reflectance <- NULL
    
  column.seq <- seq(i, 4900-(49-x), length.out = 100)
  
  node.coefs <- internal.node.data[column.seq,]
  
  node.data <- NULL
  
    for (j in 1:dim(node.coefs)[1]){
      
      y <- j
      
      node.coefficients <- node.coefs[y,]
      
      refl <- bspline.matrix %*% matrix(node.coefficients[2:39])
      
      refl <- refl + node.coefficients[1]
      
      node.data <- cbind(node.data, refl)
      
    }

  internal.node.reflectance <- cbind(internal.node.reflectance, node.data)


plot(x = 300:700, y = rep(100, 401), col = "white", ylim = c(min(internal.node.reflectance),max(internal.node.reflectance)), main = i+52)
for (i in 1:100){
  
  x <- i
  
  lines(x = 300:700, y = c(internal.node.reflectance[,x]), col = i)
  
  
}

}

# to get color of line, need pavo still
# library(photobiology)
# library(ggspectra)
# test <- raw_spct(w.length = 300:700, counts = internal.node.reflectance[,1])
# ggplot(test) + geom_line() + stat_wl_strip(ymin = 4.5, ymax = 5) + scale_fill_identity() + ylim(c(4.5, max(test$counts))) + theme_bw()
# rgb_spct(test)
# node.number <- NULL
#   
# for (k in 1:length(internal.node.labels)){
#    
#     z <- rep(k+52, 100)
#     
#     node.number <- c(node.number, z)
#     
#   }
#   
# internal.node.reflectance <- rbind(node.number, internal.node.reflectance)
#   



