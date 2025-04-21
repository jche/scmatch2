# practice code for S3 class learning

# constructor
new_Date <- function(x = double()) {
  stopifnot(is.double(x))
  structure(x, class = "Date")
}
new_Date(c( 1))

new_difftime <- function(x = double(), units = "secs") {
  stopifnot(is.double(x))
  units <- match.arg(units, c("secs", "mins", "hours", "days", "weeks"))

  structure(x,
            class = "difftime",
            units = units
  )
}
difftime(c(1, 10, 3600), units = "secs")

tmp1 <- new_difftime(c(1, 10, 3600), "secs")
class(tmp1)

new_difftime(52, "weeks")

difftime(350, 30, units = "secs" )


# Create object
my_s3_object <- list(name = "John", age = 30)

# Set class to something unique
class(my_s3_object) <- "my_custom_person"

# Define print method for this custom class
print.my_custom_person <- function(x) {
  cat("Name:", x$name, "\n")
  cat("Age:", x$age, "\n")
}

# Now this should work
print(my_s3_object)
