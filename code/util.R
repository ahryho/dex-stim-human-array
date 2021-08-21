# Function to load libraries

LoadPackages <- function(pkg.list){
  new.pkg <- pkg.list[!(pkg.list %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  suppressMessages(sapply(pkg.list, library, character.only = TRUE))
}

LoadPackagesBio <- function(pkg.list){
  new.pkg <- pkg.list[!(pkg.list %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg, dependencies = TRUE)
  suppressMessages(sapply(pkg.list, library, character.only = TRUE))
}
