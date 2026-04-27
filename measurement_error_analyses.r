library(shapes)
library(geomorph)

## Merge the two repeated measures datasets, while keeping track of the origin of each data point
## Assumes that procGPA is already available; originally used with the data saved in 
# the workspace mlmorph_dataset_wmargins_241220.RData
# labeled_data1 and labeled_data2 are the two dataframes to be merged, while procGPA1 and procGPA2 
# are the corresponding procGPA arrays of k x 2 x n dimensions (k is usually 12)
merge_repeated_measurements <- function(labeled_data1, labeled_data2, procGPA1, procGPA2){
    # Add a column to each dataframe to keep track of the origin of each data point
    labeled_data1$origin <- 1
    labeled_data2$origin <- 2
    
    # Merge the two dataframes
    merged_data <- rbind(labeled_data1, labeled_data2)
    
    # Merge the two procGPA arrays
    array_dims <- dim(procGPA1)
    procGPA_semiraw <- array(0, dim = c(array_dims[1], array_dims[2], array_dims[3] * 2))
    procGPA_semiraw[,,1:array_dims[3]] <- procGPA1
    procGPA_semiraw[,,(array_dims[3] + 1):(array_dims[3] *2)] <- procGPA2

    procGPA_joint <- plotGPA(procGPA_semiraw, rtrn = T)
    
    # Create a list of the joint data and procGPA array, name the two list elements
    joint_list <- list(labeled_data = merged_data, semiraw_coordinates = procGPA_semiraw, procGPA = procGPA_joint)

    return(joint_list)
}


# Function to calculate the differences between two sets of landmark coordinates
# arr1 and arr2 are 12 x 2 x n arrays of landmark coordinates
landmark_differences <- function(arr1, arr2){
    len <- dim(arr1)[3]
    diff <- array(0, dim = c(12, 2, len))
    for(i in 1:len){
        diff[,,i] <- arr2[,,i]  - arr1[,,i]
    }
    diff_mod <- sqrt(apply(diff^2, 3, sum))
    return(diff)
}

# Function to calculate the euclidean distance between two sets of landmark coordinates, for each
# landmark consisting of a pair of coordinates
# arr1 and arr2 are 12 x 2 x n arrays of landmark coordinates
# Returns a vector of length 12 with the euclidean distances for each landmark

# Function to calculate the mean Euclidean distance for each landmark
# arr1 and arr2 are 12 x 2 x n arrays of landmark coordinates
mean_landmark_distances <- function(arr1, arr2) {
    len <- dim(arr1)[3]
    diff <- array(0, dim = c(12, 2, len))
    
    for (i in 1:len) {
        diff[,,i] <- arr2[,,i] - arr1[,,i]
    }
    
    # Calculate the Euclidean distance for each landmark
    euclidean_distances <- sqrt(apply(diff^2, c(1, 2), sum))
    
    # Calculate the mean Euclidean distance for each landmark
    mean_distances <- apply(euclidean_distances, 1, mean)
    
    return(mean_distances)
}

# Given two [ 12 x 2 x n ] arrays of landmark coordinates (where n is equal between them), 
# calculates the Procrustes distance between each pair of corresponding configurations

pair_procDist <- function(arr1, arr2){
    len <- dim(arr1)[3]
    procDist <- array(0, dim = c(1, len))
    for(i in 1:len){
        procDist[1,i] <- procdist(arr1[,,i], arr2[,,i])
    }
    return(procDist)
}

# Given two [ 12 x 2 x n ] arrays of landmark coordinates (where n is equal between them), 
# calculates the Euclidean distance between each pair of corresponding coordinates, for each corresponding
# specimen (i.e. for each n)

pair_euclideanDist <- function(arr1, arr2){
    len <- dim(arr1)[3]
    euclideanDist <- array(0, dim = c(12, len))
    for(i in 1:len){
        euclideanDist[,i] <- sqrt(rowSums((arr1[,,i] - arr2[,,i])^2))
    }
    euclideanDist <- as.data.frame(t(euclideanDist))
    return(euclideanDist)
}

mean_shape_euclideanDist <- function(mean_shape, shapes_array) {
    len <- dim(shapes_array)[3]
    euclideanDist <- array(0, dim = c(12, len))
    for(i in 1:len){
        euclideanDist[,i] <- sqrt(rowSums((shapes_array[,,i] - mean_shape)^2))
    }
    euclideanDist <- as.data.frame(t(euclideanDist))
    colnames(euclideanDist) <- sprintf("%02d", 1:12)
    return(euclideanDist)
}

# make a ggplot boxplot of paired Euclidean distances 
# for an example of the input dataframe, see euclidean_distance_measurement_error_rescaled_final (n * 12 
# dataframe, basically just the output of pair_euclideanDist with colnames named 01, 02, ..., 12)
library(tidyr)
library(ggplot2)
plot_euclideanDist_boxplot <- function(euclideanDistances, title = "Euclidean Distances per Landmark", ylimit = c(0,0.11)) {
    colnames(euclideanDistances) <- sprintf("%02d", 1:12)
    df_long <- pivot_longer(euclideanDistances, cols = everything(), names_to = "Landmark", values_to = "Distance")

    plot <- ggplot(df_long, aes(x = Landmark, y = Distance)) +
        geom_boxplot() +
        labs(title = title, x = "Landmark", y = "Euclidean Distance") +
        #theme_cowplot(12) +
        ylim(ylimit)

    return(plot)
}




# Compute mean shapes for each repeat
mean_shape1 <- apply(full_GPA_mesErr1_afterjoint[,,1], c(1,2), mean)
mean_shape2 <- apply(full_GPA_mesErr2_afterjoint[,,2], c(1,2), mean)

# Plot displacement arrows
plot(mean_shape1, type = "n", asp = 1, main = "Landmark Displacement Due to Error")
arrows(mean_shape1[,1], mean_shape1[,2], mean_shape2[,1], mean_shape2[,2],
       length = 0.1, col = "red")
points(mean_shape1, pch = 16, col = "blue")
points(mean_shape2, pch = 16, col = "green")
legend("topright", legend = c("Repeat 1", "Repeat 2", "Displacement"),
       col = c("blue", "green", "red"), pch = c(16, 16, NA), lty = c(NA, NA, 1))

# Example mean shapes


# Calculate displacement vectors
displacement_vectors <- mean_shape2 - mean_shape1

# Scale factor
scale_factor <- 2

# Rescale displacement vectors
scaled_displacement_vectors <- displacement_vectors * scale_factor

# Plot displacement arrows
plot(mean_shape1, type = "n", asp = 1, main = "Landmark Displacement Due to Error")
arrows(mean_shape1[,1], mean_shape1[,2], 
       mean_shape1[,1] + scaled_displacement_vectors[,1], 
       mean_shape1[,2] + scaled_displacement_vectors[,2],
       length = 0.1, col = "red")
points(mean_shape1, pch = 16, col = "blue")
points(mean_shape2, pch = 16, col = "green")
legend("topright", 
       col = c("blue", "green", "red"), pch = c(16, 16, NA), lty = c(NA, NA, 1))


# function for saving procD.lm output as a csv file
library(xtable)
procD.lm_output <- function(procD.lm_obj, filename = NULL) {
    
    table <-  xtable(procD.lm_obj$aov.table)
    
    # Create a data frame with the residuals, distances, and coordinates
    #output_df <- data.frame(residuals = residuals, distances = distances, coords)

    if(is.null(filename)) {
        filename <- "procD_lm_output.csv"
    }
    
    # Write the data frame to a CSV file
    write.csv(table, file = filename, row.names = T)
}

# Calculate the mean value and standard deviation of each column in a dataframe, and then 
# the overall mean of the column means and the overall standard deviation of the column means

mean_std_df <- function(df) {
    col_means <- apply(df, 2, mean)
    col_stds <- apply(df, 2, sd)
    
    overall_mean <- mean(col_means)
    overall_std <- sd(col_means)
    
    return(list(col_means = col_means, col_stds = col_stds, overall_mean = overall_mean, overall_std = overall_std))
}

mean_shape <- apply(full_GPA_mesErr_afterjoint, c(1,2), mean)