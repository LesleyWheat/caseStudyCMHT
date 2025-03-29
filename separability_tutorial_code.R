# Description:
# This program demonstrates the behaviour of several measures from the ECoL
# library at https://github.com/lpfgarcia/ECoL .
# The measures are used on a generated ball inside a hypersphere with 
# variable variance, number of features, number of samples per class and number 
# of iterations.
# For the noise, the features are assumed to be independent and have the same
# variance value.
# Note: this program is not multi-threaded.

#------------------------------------------------------------------------------
#### Setup ####

# Clear environment data/variables
rm(list=ls()) 

# Enable libraries
library(gena) # Hypersphere
library(MASS) # Multivariate normal distribution
library(ECoL) # Complexity measures
library(ggplot2) # For plotting graphs
library(investr) # Confidence levels for graphs

#------------------------------------------------------------------------------
# Parameters

# Output directory name
output_folder = 'output'

# Set seed for reproducible results
set.seed(123)

# Samples per class
# Note large numbers of samples will significantly increase runtime
sample_size_list <- c(100, 500)

# Data dimensions that will be used
# Number of features
NOF_list <- c(1, 10, 20) 

# Generate a range of Variances that will be used
variance_list <- append(seq(30, 40, by=2),seq(10, 30, by=1))
variance_list <- append(variance_list, seq(5, 10, by=1))
variance_list <- append(variance_list, seq(1, 5, by=0.5))

# Number of randomly generated datasets generated for each parameter combination
runs <- 3

#### end ####
#------------------------------------------------------------------------------
#### Setup ####
# Create output directory
dir.create(file.path(getwd(), output_folder), showWarnings = FALSE)

#### Generate Data ####
# This section generates the simulated data based on the parameters
# and calculates the complexity measure results. These results are saved 
# in a dataframe and the data is discarded.

# Starting index for keeping track of iterations
index_total = 1

# Create dataframe for saving measure data
results <- data.frame(matrix(ncol = 14, nrow = 0))
column_names <- c("index_total", "index_run", "distribution_type",
                  "points_per_class","NOF","variance", "runtime",
                  "N1","N2", "N3", "T1", "LSC", "F1v", "F2")
colnames(results) <- column_names

# What measures the key codes stand for:
#   N1: Fraction of points lying on the class boundary
#   N2: Average intra/inter class nearest neighbor distances
#   N3: Leave-one-out error rate of the 1-nearest neighbor algorithm
#   N4: Nonlinearity of the one-nearest neighbor classifier
#   T1: Fraction of maximum covering spheres on data
#   LSC: Local-Set cardinality average
#   F1v: The directional-vector Fisher's discriminant ratio
#   F2: Overlapping of the per-class bounding boxes

# Run through all dataset combinations to create datasets and measure them
for (NOF in NOF_list){
  print(paste("Starting NOF: ", NOF)) # Inform user where program is
  
  for (sample_size in sample_size_list){
    print(paste("Starting sample_size: ", sample_size)) # Inform user
    
    for (i in 1:runs){
      for (variance in variance_list){
        # Save current time
        start_time <- Sys.time()
          
        # Create diagonal covariance matrix, assuming each feature is independent
        # Scaling the variance
        svar <- variance*(1/sqrt(NOF))
        covariance <- diag(svar, NOF, NOF)
        
        # Generate class A from a multivariate normal distribution
        classA <- mvrnorm(n = sample_size,
                          mu = numeric(NOF), 
                          Sigma = covariance)
        
        # Generate class B from a sphere
        pointsB <- rhypersphere(n = sample_size, 
                                dim = NOF, 
                                radius = 10, 
                                center = numeric(NOF),
                                type = "boundary")
        # Add noise
        noise <- mvrnorm(n = sample_size,
                         mu = numeric(NOF), 
                         Sigma = covariance)
        
        classB <- pointsB + noise
        
        # Format data for measures
        class_id <- c(rep(1, sample_size), rep(2, sample_size))
        
        normalize = function(v) (v - min(v, na.rm = TRUE)) / diff(range(v, na.rm = TRUE))
        class_data <- data.frame(lapply(data.frame(rbind(classA, classB)), normalize))
        
        # Get measures
        neighbourhood_measures <- data.frame(neighborhood(class_data, class_id))
        feature_measures <- data.frame(overlapping(class_data, class_id,
                                                   measures = c("F1v","F2"), summary=c("max")))
        
        # Save run time
        runtime <- Sys.time() - start_time
        
        # Compile run results
        run_data = list(index_total, i, "ball in sphere", sample_size, 
                        NOF, variance, runtime,
                        neighbourhood_measures["mean","N1"],
                        neighbourhood_measures["mean","N2"],
                        neighbourhood_measures["mean","N3"],
                        neighbourhood_measures["mean","T1"], 
                        neighbourhood_measures["mean","LSC"],
                        feature_measures["max","F1v"])
        
        # Save results
        results[nrow(results) + 1,] <- run_data
        
        # Increment index
        index_total <- index_total + 1
      }
    }
  }
}

print("Completed data generation")
#### end ####
#------------------------------------------------------------------------------
#### Generate Graphs ####

# Display all the graphs for a specific set of features

list_of_y_labels = c("N1","N2", "N3", "T1", "LSC", "F1v") # Measures to display

# Iterate through graphs to make
for (graph_NOF in NOF_list){
for (graph_samples in sample_size_list){
for(y_index in 1:length(list_of_y_labels)){
  # Setup variables
  subtitle_short <- paste(paste(graph_NOF, "NOF,"), # Setup graph subtitle
                          paste(graph_samples, "samples per class"))
  
  x_label <- sym('variance')
  y_label <- sym(list_of_y_labels[y_index])
  
  # Let's skip curve fitting for F1v as it won't work
  if (list_of_y_labels[y_index] != 'F1v'){
  # Formula for curve fitting
  fit_formula <- as.formula(paste(paste(y_label, ' ~ SSlogis('), 
                                  paste(x_label, ', Asym, xmid, scal) + y_offset')))
  
  # Set initial guesses for regression
  y_min <- min(results[results$NOF == graph_NOF & results$points_per_class == graph_samples,
                       list_of_y_labels[y_index]])
  y_max <- max(results[results$NOF == graph_NOF & results$points_per_class == graph_samples,
                       list_of_y_labels[y_index]])
  x_min <- min(results[results$NOF == graph_NOF & results$points_per_class == graph_samples,
                       'variance'])
  x_max <- max(results[results$NOF == graph_NOF & results$points_per_class == graph_samples,
                       'variance'])
  
  # These starters are set up for scaled variance and a radius of 100
  IG_up <- list(Asym = y_max - y_min, xmid = 30, scal = 5, y_offset = y_min)
  IG_up_2 <- list(Asym = y_max - y_min, xmid = 12, scal = 5, y_offset = y_min)
  IG_up_3 <- list(Asym = y_max - y_min, xmid = 3, scal = 5, y_offset = y_min)
  IG_down <- list(Asym = y_min - y_max, xmid = 30, scal = 5, y_offset = y_max)
  IG_down_2 <- list(Asym = y_min - y_max, xmid = 12, scal = 5, y_offset = y_max)
  IG_down_3 <- list(Asym = y_min - y_max, xmid = 3, scal = 5, y_offset = y_max)
  IG_flat <- list(Asym = y_min - y_max, xmid = x_max - x_min, scal = 1, y_offset = y_max)
  IG_flat2 <- list(Asym = y_max - y_min, xmid = x_max - x_min, scal = 3, y_offset = y_min)
  IG_downSmall <- list(Asym = y_max - y_min, xmid = 0.03, scal = 0.02, y_offset = y_max)
  
  IG_list = list(IG_up, IG_up_3, IG_up_2,
                 IG_down, IG_down_2, IG_down_3,
                 IG_flat, IG_flat2, IG_downSmall)
  
  # Try different guesses for curve fitting until one works
  for(guess in IG_list){
    # Setup variables and options for curve fitting
    nls_options <- nls.control(tol = 0.01, minFactor = 1/2048)
    nls_success <- FALSE
    
    # Catch errors from fit function
    tryCatch({
      # Attempt to fit the logistic curve to the data
      fit <- nls(formula = fit_formula, 
                 results[results$NOF == graph_NOF & results$points_per_class == graph_samples, ],
                 control = nls_options,
                 start = guess)
      
      # If no error, then fitting worked and stop trying
      nls_success <- TRUE
      break
    }, error = function(E1){
      # Error capture
      # Ignore error
    })
  }
  
  # If first fit did not work, try again with different parameters
  if (!nls_success){
    # Inform the user
    print(paste(paste("NLS fitting failed for:", list_of_y_labels[y_index]), 
                ". Retrying with lower tolerance"))
    
    # Change fitting options to increase chance of success
    nls_options <- nls.control(tol = 0.05, minFactor = 1/(2*2048))
    for(guess in IG_list){
      tryCatch({
        fit <- nls(formula = fit_formula, 
                   results[results$NOF == graph_NOF & results$points_per_class == graph_samples, ],
                   control = nls_options,
                   start = guess)
        nls_success <- TRUE
        break
      }, error = function(E1){
        # Error capture
        # Ignore error
      })
    }
    # Inform user
    if (!nls_success){
      print(paste("NLS fitting again failed for: ", list_of_y_labels[y_index]))
    }else{
      print(paste("NLS fitting success on retry for: ", list_of_y_labels[y_index]))
    }
  }
  }else{
    nls_success = FALSE
  }
  
  # Calculate prediction intervals if fitting is successful
  if (nls_success){
    tryCatch({
    cimat1 <- predFit(fit,
                      interval = "prediction",
                      level = 0.95,
                      se.fit = FALSE)
    # Select data for the plot
    plot <- ggplot(data = cbind(results[results$NOF == graph_NOF & results$points_per_class == graph_samples, ], cimat1), 
                 aes(x=!!x_label))
    }, error = function(E1){
      # Select data for the plot
      plot <- ggplot(data = cbind(results[results$NOF == graph_NOF & results$points_per_class == graph_samples, ]), 
                     aes(x=!!x_label))
    })
  }else{
    # Select data for the plot
    plot <- ggplot(data = cbind(results[results$NOF == graph_NOF & results$points_per_class == graph_samples, ]), 
                   aes(x=!!x_label))
  }
  
  # Add points to plot
  plot <- plot + geom_point(mapping = aes(x=!!x_label, y=!!y_label))
  
  # If curve fitting was successful, then add it to the plot
  if (nls_success){
    plot <- plot + geom_ribbon(aes(ymin = lwr, ymax = upr),
                               color = "grey70",
                               alpha = 0.25,
                               fill = "grey50")
    plot <- plot + geom_line(aes(y = fit), color = "black")
  }
  
  # Set plot style
  plot <- plot + theme_light()
  plot <- plot + labs(x = 'Scaled variance', y = y_label, 
              subtitle = paste(subtitle_short, "- approximate 95%-prediction intervals"), 
              title = paste(y_label, paste("change with", x_label)))
  
  # Save plot to output directory
  ggplot2::ggsave(
    paste(y_label, subtitle_short,".png"),
    plot = last_plot(),
    device = 'png',
    width = 6, 
    height = 4,
    path = file.path(getwd(), output_folder)
  )

}
}
}

#### end ####
