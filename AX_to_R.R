# script to convert matlab .m files to .r files

library(matconv) # is this equivalent to import? 
library(stringr)
filelist <- list.files(pattern="*.m")
for (file in filelist){
  filename <- paste(file[-2])
  print(filename)
  filename <- str_sub(filename, end=-3)
  new_file <- paste(filename, ".R", sep="")
  print(new_file)
    tryCatch( {
        suppressWarnings(mat2r(file, new_file, verbose = 1))
        },
        error = function(cond) {
          message(paste("Problem with file:", file))
          message("Here's the original error message:")
          message(conditionMessage(cond))
          # Choose a return value in case of error
          NA
        }
    )

}
