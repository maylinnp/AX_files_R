# Function to read all the lines in a file
# returns a variable containing all lines

read_file = function(filepath) {
  lines = c()
  con = file(filepath, "r")
  BWD <- FALSE
  line_no <- 0
  while ( TRUE ) {
    line_no <- line_no +1
    line = readLines(con, n = 1)
    # stop at the first empty line
    if ( length(line) == 0 ) {
      break
    }
    # capture if there is back titration data
    if (grepl("BWD", line)){
      BWD <- TRUE
      eoFWD <- line_no-1
    } 
    lines = append(lines, line)
    if (length(lines) == 1 ){
      first_line = line
    }
  }
  if (eoFWD == 0 ){
    eoFWD <- length(lines)
  }
  close(con)
  return (list(first_line = first_line, lines = lines, BWD = BWD, eoFWD = eoFWD))
}