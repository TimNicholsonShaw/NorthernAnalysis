library('tidyverse')
library('ggplot2')

read_in_and_preprocessor <- function(file_loc) {
    df <- read.csv(file_loc) # read as CSV

    df$Replicate <- as.factor(df$Replicate) # make sure that replicate isn't treated as numbers
    df$Treatment <- as.factor(df$Treatment) # factorize everything else
    df$Substrate <- as.factor(df$Substrate)

    # Normalize all columns to 100% abundance
    df$X120 <- df$X120/df$X0
    df$X240 <- df$X240/df$X0
    df$X360 <- df$X360/df$X0
    df$X0 <- df$X0/df$X0

    # Pivot to a version that's more amenable to ggplot
    df <- pivot_longer(df, cols=c("X0", "X120", "X240", "X360"), names_to="Time", values_to="Abundance")

    # Remove the X added to Time column and convert to an integer
    df$Time <- sub(".", "", df$Time)
    df$Time <- as.integer(df$Time)

    return(df)
}

