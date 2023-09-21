
kingsize <- function(dataframe, windowsize) {
  tic("time for trimming: ")
  # Calculate the midpoint of the original start and end coordinates
  midpoint <- floor(dataframe$start + (dataframe$stop - dataframe$start) / 2)
  
  # Calculate the new start and end coordinates based on the midpoint and windowsize
  dataframe$start <- floor(midpoint - (windowsize / 2))
  dataframe$stop <- floor(midpoint + (windowsize / 2))
  
  toc()
  # Return the updated dataframe
  return(dataframe)
}
