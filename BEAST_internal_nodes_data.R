library(splines2)
library(ggplot2)
library(tidyverse)
library(tidylog)
library(pavo)
library(photobiology)
library(ggspectra)
library(ape)
library(tictoc)
library(here)
here()
load("BEAST.trees.subset.RData")

# Creates a basis matrix with our respective domain (300-700 nm)
bspline.matrix <- bSpline(300:700, df = 38)

# for loop method of extracting our coefficients and generating the predicted reflectance of root node
# Extracting root node ASR coefficients for 100 trees
root.reflectance.ASR <- NULL
for (i in 1:dim(root.spline.coef)[1]){
  
  x <- i
  
  root.coef <- root.spline.coef[x,]
  
  root.refl <- bspline.matrix %*% matrix(root.coef[2:39])
  
  root.refl <- root.refl + root.coef[1]
  
  root.reflectance.ASR <- cbind(root.reflectance.ASR, root.refl)
  
}

# Plots all of the reflectance curves prior to correcting the negative values
matplot(x = 300:700, y = root.reflectance.ASR[,-1], lwd = 2, type = "l", xlab = "Wavelength (nm)", ylab = "Reflectance (%)", main = "Root ASR BEAST")

# Turning our reflectance data into a dataframe
root.reflectance.ASR <- data.frame(wl = 300:700, root.reflectance.ASR)

# Turning our reflectance dataframe into an rspec object
root.reflectance.ASR <- as.rspec(root.reflectance.ASR)

# Correcting our negative reflectance values into 0
root.reflectance.ASR <- procspec(root.reflectance.ASR, fixneg = "zero")

# Plotting our corrected reflectance curves
matplot(x = root.reflectance.ASR$wl, y = root.reflectance.ASR[,-1], pch = 19, cex = 0.5)


#####
#Repeating the same process but in the tidyverse

# Creates a column for run number and the intercept coefficient of each run for the root
root.reflectance.BEAST <- tibble(run = 1:100, root_intercept = root.spline.coef[,1])

# Splits the rows (non-intercept spline coefficients) for each run of the root and adds the spline coefficient lists to our tibble 
root.reflectance.BEAST <- root.reflectance.BEAST %>% mutate(spline_coefs = split(root.spline.coef[,-1], f = 1:100, 38))

# Takes the spline coefficients our basis matrix and creates reflectance values (stored as lists) into our tibble
root.reflectance.BEAST <- root.reflectance.BEAST %>% group_by(run) %>% mutate(reflectance = list(bspline.matrix %*% matrix(unlist(spline_coefs), ncol = 1))) %>% ungroup()

# Rescales the reflectance values by adding the intercept coefficient to those values
root.reflectance.BEAST <- root.reflectance.BEAST %>% group_by(run) %>% mutate(rescaled_reflectance = list(unlist(reflectance) + root_intercept)) %>% ungroup()

# Adds a column containing the rspec dataframes to the tibble  
root.reflectance.BEAST <- root.reflectance.BEAST %>% group_by(run) %>% mutate(rspec = list(procspec(as.rspec(data.frame(wl = 300:700, reflectance = unlist(rescaled_reflectance))), fixneg = "zero"))) %>% ungroup()

# Adds a column to the tibble that has the rgb color from the rspec object
root.reflectance.BEAST <- root.reflectance.BEAST %>% group_by(run) %>% mutate(color = spec2rgb(data.frame(rspec))) %>% ungroup()

# Creates a separate tibble that unlists the rescaled_reflectance (stored as a singular vector)
plot.root.reflectance <- root.reflectance.BEAST %>% unnest(cols = rescaled_reflectance)

# Removing any unnecessary columns in our plotting tibble
plot.root.reflectance <- plot.root.reflectance %>% select(-c(root_intercept, spline_coefs, reflectance,rspec, color))

# Manually 'procspec' our rescaled_reflectance
plot.root.reflectance$rescaled_reflectance[plot.root.reflectance$rescaled_reflectance < 0] <- 0

# Adds a wavelength column
plot.root.reflectance <- plot.root.reflectance %>% mutate(wl = rep(300:700, 100))

# Plots each predicted reflectance curve from the BEAST runs
plot.root.reflectance %>% group_by(run) %>% ggplot(aes(x = wl, y = rescaled_reflectance)) + geom_point(alpha = 0.25, size = 0.75, color = as.character(rep(root.reflectance.BEAST$color, each = 401))) + theme() + theme_bw() + stat_wl_strip(ymin = -1.5, ymax = -0.5) + scale_fill_identity() + labs(x = "Wavelength (nm)", y = "Reflectance (%)", title = "Root Reflectance Predicted by BEAST")

#####
# Plotting mean root reflectance using a 1 sd interval

mean_root <- root.reflectance.ASR %>% column_to_rownames(var = "wl")

means <- rowMeans(mean_root)

mean_root <- mean_root %>% mutate(stdev = apply(mean_root, 1, sd)) %>% select(c(stdev, c(-stdev))) %>%  mutate(wl = c(300:700), wl_means = means)

mean_root %>% ggplot(aes(x = wl, y = means)) + geom_ribbon(aes(ymin = wl_means - stdev, ymax = wl_means + stdev), fill = "lightgray", alpha = 0.5) + geom_line(lwd = 2, color = as.character(spec2rgb(procspec(as.rspec(data.frame(wl = mean_root$wl, reflectance = mean_root$wl_means)), fixneg = "zero")))) + theme_bw() + labs(x = "Wavelength (%)", y = "Reflectance (%)", title = "Root Reflectance Predicted by BEAST", subtitle = "Includes 1 standard deviation around mean")

#####
# Applying what we did but with the internal nodes

# Creates a tibble with the trees drawn from our MCMC
intnode.reflectance.BEAST <- tibble(tree_name = rep(names(BEAST.trees.subset), each = length(53:101)))

# Adds 3 new columns
# 1. a series of node numbers used in BEAST
# 2. a series of node numbers according to the tree 
# NOTE: Tree is labeled 1:101, with the first 51 points being the tips and the rest are labeled working from the root through the tree to tips, starting from the bottom and working its way up to the top of tree. BEAST, however, only does the latter, incorporating the tips within the flow of the tree. For example, node 1 in BEAST is actually node 53 on the tree, and node 4 in BEAST would be the first tip (TanCyo).
# 3. a series in order to join our tibble with another tibble holding the coefficient information. Need to join using a column that has unique values in each row otherwise it will create unnecessary duplicates
intnode.reflectance.BEAST <- intnode.reflectance.BEAST %>% mutate(nodes_list = rep(c(1,2,3,5,8,9,10,11,12,14,17,20,21,24,26,29,31,33,34,35,38,40,44,47,49,51,54,55,56,58,60,64,65,67,68,69,70,72,76,77,81,83,85,88,89,90,92,94,97), 100), nodes_internal = rep(c(53:101), 100)) %>% mutate(join_list = 1:length(nodes_list))

# Function used to extract the internal node data from the BEAST output
extraction <- function(multiPhylo, internal_node_values){
  
  if (class(multiPhylo) != "multiPhylo"){
    stop("You are not using a 'multiPhylo' object!")
  }
  
  # Storage vector
  internal_node_coefs <- NULL
  
  # for loop that takes the annotations from the respective tree and stores it into a vector
  for (run in names(multiPhylo)){
    
    print(run)
    
    spline_coefs <- multiPhylo[[run]]$annotations
    
    node_coefs <- NULL
    
    # for loop that extracts the spline coefficients from an internal node from each tree 
    for (i in 1:length(spline_coefs)){
      
      node <- internal_node_values[i]
      
      print(node)
      
      node_coefs[i] <- spline_coefs[[node]]
      
    }
    
    internal_node_coefs <- c(internal_node_coefs, node_coefs)
    
  }
  return(tibble(internal_node_coefs))
}

# Applying our function to the BEAST trees data
intnode.data <- extraction(BEAST.trees.subset, internal_node_values = intnode.reflectance.BEAST$nodes_list[1:49])

# Using the 'extraction' function creates a bunch of NULL values for the nodes that were not selected for
# We remove them here, add a nodes list column (using the BEAST node naming conventions), and add a join_list column to add this tibble to our reflectance tibble for the internal nodes
intnode.data <- intnode.data %>% rowwise() %>% filter(class(internal_node_coefs) == "list") %>% ungroup() %>% mutate(nodes_list = rep(c(1,2,3,5,8,9,10,11,12,14,17,20,21,24,26,29,31,33,34,35,38,40,44,47,49,51,54,55,56,58,60,64,65,67,68,69,70,72,76,77,81,83,85,88,89,90,92,94,97), 100)) %>% mutate(join_list = 1:length(nodes_list))

# Combine the two tibbles and compare the node_list columns to ensure the data joined correctly
intnode.reflectance.BEAST <- intnode.reflectance.BEAST %>% left_join(intnode.data, by = "join_list")

# Removes columns not needed in the future
intnode.reflectance.BEAST <- intnode.reflectance.BEAST %>% select(-c(nodes_list.y, join_list)) %>% rename(nodes_list = nodes_list.x)




# The code below follows a near-identical process to that of what we used with the root node. The only minor change is that there has to be an initial line of code to filter which node that you want to specifically look at
# Case study using node 65
test.node <- intnode.reflectance.BEAST %>% filter(nodes_internal == 65)

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
plot.test.node %>% rowwise() %>% ggplot(aes(x = wl, y = rescaled_reflectance)) + geom_point(alpha = 0.25, size = 0.75, pch = 19, color = as.character(rep(test.node$color, each = 401))) + theme() + theme_bw() + stat_wl_strip(ymin = -1.5, ymax = -0.5) + scale_fill_identity() + labs(x = "Wavelength (nm)", y = "Reflectance (%)", title = paste("Node", as.character(unique(test.node$nodes_internal)), "Reflectance Predicted by BEAST"))

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
    plot.test.node <- plot.test.node %>% rowwise() %>% ggplot(aes(x = wl, y = rescaled_reflectance)) + geom_point(alpha = 0.25, size = 0.75, pch = 19, color = as.character(rep(test.node$color, each = 401))) + theme() + theme_bw() + stat_wl_strip(ymin = -1.5, ymax = -0.5) + scale_fill_identity() + labs(x = "Wavelength (nm)", y = "Reflectance (%)", title = paste("Node", as.character(test.node$nodes_internal), "Reflectance Predicted by BEAST"))
    
    return(plot.test.node)
  }
  
  if (length(nodes) > 1){
    
    for (i in 1:length(nodes)){
      
      node <- nodes[i]
      
      test.node <- dataframe %>% filter(nodes_list == node)
      
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
      plot.test.node <- plot.test.node %>% rowwise() %>% ggplot(aes(x = wl, y = rescaled_reflectance)) + geom_point(alpha = 0.25, size = 0.75, pch = 19, color = as.character(rep(test.node$color, each = 401))) + theme() + theme_bw() + stat_wl_strip(ymin = -1.5, ymax = -0.5) + scale_fill_identity() + labs(x = "Wavelength (nm)", y = "Reflectance (%)", title = paste("Node", as.character(test.node$nodes_internal), "Reflectance Predicted by BEAST"))
      
      print(plot.test.node)
    }
    
  }
  
}

tanagerTree <- read.nexus("tanagerTree.NEXUS")

tic()
pdf("BEAST_ASR_no_root.pdf", width = 8, height = 6)
plot(tanagerTree)
nodelabels(cex = 0.75, frame = "circle", bg = "yellow")
suppressMessages({ASR_plotter(dataframe = intnode.reflectance.BEAST, nodes = c(unique(intnode.reflectance.BEAST$nodes_list)))})
dev.off()
toc()

#####
# Spline coefficient densities at each node for 100 runs

test <- intnode.reflectance.BEAST %>% unnest(cols = internal_node_coefs)

test <- test %>% mutate(coef_name = rep(c("Int.", "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38"), 4900))

test2 <- test %>% filter(coef_name == "Int." & nodes_internal == 53) %>%  mutate(join_list = 1:length(coef_name))

test2 %>% ggplot(aes(x = unlist(internal_node_coefs))) + geom_density(fill = "blue")

test3 <- test %>% filter(coef_name == "Int." & nodes_internal == 54) %>% mutate(join_list = 1:length(coef_name))

test3 %>% ggplot(aes(x = unlist(internal_node_coefs))) + geom_density(fill = "red")

test4 <- test2 %>% left_join(test3, by = "join_list")

test4 %>% ggplot() + geom_density(aes(x = unlist(internal_node_coefs.x)), fill = "black", alpha = 0.5) + geom_density(aes(x = unlist(internal_node_coefs.y)), fill = "yellow2", alpha = 0.5) + guides(fill = guide_legend(title = "Node")) + theme() + scale_fill_manual(values = c("53", "54"))

test5 <- test4 %>% pivot_longer(cols = c(internal_node_coefs.x, internal_node_coefs.y))

test5 %>% ggplot(aes(x = unlist(value), fill = name)) + geom_density(alpha = 0.25) + scale_fill_manual(labels = c(53,54), name = "Node", values = c("black", "yellow2"))
######

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




test <- tibble(reflectance = test.node$rescaled_reflectance)%>% unnest(cols = reflectance) %>% mutate(wl = rep(c(300:700), length(reflectance)/length(300:700)), tree = rep(unique(test.node$tree_name), each = 401*49), node = rep(c(rep(unique(test.node$nodes_internal), each = 401)), 100)) %>% pivot_wider(names_from = "wl", values_from = "reflectance")

# test plot
plot <- test %>% filter(node == 65)

plot <- plot %>% select(-c(tree, node))

plot <- data.frame(wl = 300:700, t(plot)) %>% as_tibble() %>% column_to_rownames(var = "wl")

means <- rowMeans(plot)


plot <- plot %>% mutate(stdev = apply(plot, 1, sd)) %>% select(c(stdev, c(-stdev))) %>%  mutate(wl = c(300:700), wl_means = means)

plot %>% ggplot(aes(x = wl, y = wl_means)) + geom_ribbon(aes(ymin = wl_means - stdev, ymax = wl_means + stdev), fill = "lightgray", alpha = 0.5) + geom_line(lwd = 2, color = as.character(spec2rgb(procspec(as.rspec(data.frame(wl = plot$wl, reflectance = plot$wl_means)), fixneg = "zero")))) + theme_bw() + labs(x = "Wavelength (%)", y = "Reflectance (%)", title = "Node 65 Reflectance Predicted by BEAST", subtitle = "Includes 1 standard deviation around mean")

# Generalized plotter
ASR_plotter_v2 <- function(dataframe, nodes){
  
  if (sum(names(dataframe) == "reflectance") == 0){
    stop("Dataframe needs a column titled 'reflectance'.")
  }
  
  if (sum(names(dataframe) == "") == 0){
    stop("Dataframe needs a column titled 'internal_node_coefs'.")
  }
  
  if (class(dataframe$internal_node_coefs) != "list"){
    stop("Spline coefficient column needs to be a list.")
  }
  
  if (length(nodes) == 1){
  
  
  
  
  
  
}}
               





# not fixed topology internal node extractor

x <- BEAST.trees.subset[["STATE_5604"]]

edges <- data.frame(x$edge)

edges <- data.frame(edges, BEAST_annotation_number = c(1:length(c(edges[,1]))))

root <- edges[1,1]

edges <- as_tibble(edges) %>% filter(X2 > root)

edges <- edges %>% filter(X2 > root)

y <- BEAST.trees.subset[["STATE_14334"]]

BEAST_annotation_finder_for_internal_nodes <- function(tree){
  
  if (class(tree) == "phylo" | length(tree) == 1){
  # Extracts the edges data frame from the tree
  edges <- data.frame(tree$edge)
  
  # Adds the internal notation values to the edges data frame
  edges <- data.frame(edges, BEAST_annotation_number = c(1:length(c(edges[,1]))))
  
  # Extracts the root node value
  root <- edges[1,1]
  
  # Filters the data frame to show which values from the BEAST internal node values correspond to the actual tree internal nodes
  edges <- as_tibble(edges) %>% filter(X2 > root)
  
  # Returns a vector of the internal nodes that should be drawn from the tree
  return(edges$BEAST_annotation_number)
  }
  
  if (class(tree) == "multiPhylo"){
    
    tree_annotations <- NULL
    
    tree_names <- NULL
    
    for (i in 1:length(tree)){
      
      tree_subset <- tree[[i]]
      
      # Extracts the edges data frame from the tree
      edges <- data.frame(tree_subset$edge)
      
      # Adds the internal notation values to the edges data frame
      edges <- data.frame(edges, BEAST_annotation_number = c(1:length(c(edges[,1]))))
      
      # Extracts the root node value
      root <- edges[1,1]
      
      # Filters the data frame to show which values from the BEAST internal node values correspond to the          actual tree internal nodes
      edges <- as_tibble(edges) %>% filter(X2 > root)
      
      # Combines data into data frame
      tree_annotations <- cbind(tree_annotations, edges$BEAST_annotation_number)
      
    }
    
    return(tree_annotations)
  }
}

