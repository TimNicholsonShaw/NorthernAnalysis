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

find_half_life <- function(df) { # returns half life in minutes
    # really good explanation of non-linear fitting https://rpubs.com/mengxu/exponential-model
    # trying to fit to model: y = alpha*e^(beta*x) + theta
    # need to find good starting values before fitting can start

    # select approximate theta. Needs be lower than the minimum of Abbundance, but greater than 0
    theta.0 <- min(df$Abundance) * 0.5 

    #Estimate the other parameters using a linear model
    model.0 <- lm(log(Abundance - theta.0) ~ Time, data=df)
    alpha.0 <- exp(coef(model.0)[1])
    beta.0 <- coef(model.0)[2]

    # Starting paramteters
    start <- list(alpha=alpha.0, beta=beta.0, theta=theta.0)

    # make the model using starting values
    model <- nls(Abundance ~ alpha * exp(beta * Time) + theta, data=df, start=start)

    beta = coef(model)[2]
    return(log(0.5)/beta)
}