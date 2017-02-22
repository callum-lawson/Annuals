# Trial R script
# This is to see if a new R script can be added

newfunction <- function() 
  print("hello world")

# Some more changes that will be deleted

dir.create(file.path(getwd(), "Output"), showWarnings = FALSE)
  # creates directory if it doesn't exist already
write.csv(1:10, file = "Output/MyData.csv")
