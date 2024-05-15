### Installing a list of packages from a file (and their dependencies)###

# Read the package names from the text file
package_names <- readLines("package_list.txt")

# Install packages that are not already installed
new_packages <- package_names[!(package_names %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
lapply(package_names, require, character.only = TRUE)



