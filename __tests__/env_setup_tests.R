rm(list = ls())

# Check if the code is running on GitHub Actions
if (Sys.getenv("GITHUB_ACTIONS") == "true") {
  # Set the working directory for GitHub Actions
  setwd("/home/runner/work/SBARTProject/SBARTProject")
} else {
  # Set the working directory for your local machine
  setwd("/home/paolo/Bayesian_Statistics/Project23-24/SBARTProject") 
}