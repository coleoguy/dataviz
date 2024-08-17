
# Function
SeqDef <- function(tree, df, data.col = 2, invert = F, scale = T){
  # This checks to make sure that the number of tips and 
  # dataframe rows are the same.
  if(length(tree$tip.label) != nrow(df)){
    stop("Tree and Data should have same Length")
  } 
  # This line finds the max branching time of our phylo object (tree)
  td <- max(branching.times(tree))
  #  Computes the pairwise distances between the pairs of tips from a
  #  phylo tree using its branch lengths. 
  dist.matrix <- cophenetic.phylo(tree)
  # This line finds the individual concordance of each tip and divides 
  # it by the max branch length to find a proportion ?
  dist.prop <- 1 - ((dist.matrix/2) / td)
  # Create a vector with length equal to the number of labels on the tree. 
  synscores <- vector(length = length(tree$tip.label))
  # Assign the tip labels of our tree to the names of the synscores vector. 
  names(synscores) <- tree$tip.label
  # This loop begins by stating that we will be iterating 
  # through values of 1: length of tree tips, and use these values to 
  # create a 'focal.tip' that we will use for our calculations. 
  # We then make a vector (x) which is equal to the length of the number 
  # of tree tips. 
  # we then create a nested loop that will iterate in values of 1:length 
  # of tree tips, and use these values to create a current tip ('cur.tip'). 
  for(i in 1:length(tree$tip.label)){
    focal.tip <- tree$tip.label[i] 
    x <- vector(length = length(tree$tip.label))  
    for(j in 1:length(tree$tip.label)){
      cur.tip <- tree$tip.label[j]
      # Now that we have our loop to look at our two tips, we can do calculation. 
      # From our matrix dist.prop, we select the value corresponding to the 
      # two tips of interest and multiply it by a specific data column that 
      # corresponds to our current tip.
      # We assign each of these values for a given current tip to a vector (x). 
      x[j] <- dist.prop[row.names(dist.prop) == cur.tip, 
                        colnames(dist.prop) == focal.tip] * 
        df[, data.col][df[, 1] == cur.tip]
    }
    # If true, we do 1 - the mean of (x). otherwise just do mean x. 
    # store this value into synscores. 
    if(invert){
      synscores[i] <- 1 - mean(x) 
    }else{
      synscores[i] <- mean(x) 
    }
  }
  # If scale = T, we will scale the data in such a way that the smallest value 
  # becomes 0, and the largest value becomes 1. 
  if(scale){
    synscores <- (synscores - min(synscores)) /
      (max(synscores) - min(synscores))
  }
  # Prepare the results
  # Below we assign a specific column of our data frame to 'empscores'. 
  empscores <- df[, data.col]
  # we create a list (results) that consists of our phylo object, synscores, 
  # and empscores. 
  results <- list(tree, synscores, empscores)
  names(results) <- c("tree","seqdef","empdata")
  # we assign the class of our results to 'seqdef'
  class(results) <- "seqdef"
  return(results)
}

