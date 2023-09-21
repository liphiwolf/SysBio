queensize <- function(folder_path, window_size, overwrite = FALSE) {
  tic("time for trimming: ")
  
  if (!overwrite) {
    # List to store the modified data frames
    modified_data_frames <- list()
  }
  # Check if folder_path is a directory or a single file
  if (file.info(folder_path)$isdir) {
    # If folder_path is a directory, get all the .bed files in it
    bed_files <- list.files(path = folder_path, pattern = "\\.bed$", full.names = TRUE)
  } else {
    # If folder_path is a single file, process it directly
    bed_files <- folder_path
  }
  
  # Create a progress bar
  pb <- progress_bar$new(total = length(bed_files), format = "[:bar] :percent :eta")
  
  # Loop through each BED file
  for (bed_file in bed_files) {
    # Read the BED file into a data frame
    bed_data <- read.table(bed_file, header = FALSE, sep = "\t")
    
    # Perform the operations on each row
    bed_data[, 2] <- floor((bed_data[, 2] + bed_data[, 10]) - (window_size/2))
    bed_data[, 3] <- bed_data[, 2] + window_size
    bed_data[, 10] <- round(window_size/2)
    
    if (overwrite) {
      # Overwrite the original BED file with the updated data
      write.table(bed_data, file = bed_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    } else {
      # Save the modified data frame to the list
      modified_data_frames[[basename(bed_file)]] <- bed_data
    }
    
    # Increment the progress bar
    pb$tick()
  }
  toc()
  # Return a message indicating the operation is complete
  return("Trimming complete. Original BED files have been updated.")
}
