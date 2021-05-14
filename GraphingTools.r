library('tidyverse')
library('ggplot2')

read_in_and_preprocessor <- function(file_loc) {
    df <- read.csv(file_loc) # read as CSV

    df$Replicate <- as.factor(df$Replicate) # make sure that replicate isn't treated as numbers
    df$Treatment <- as.factor(df$Treatment) # factorize everything else
    df$Treatment <- relevel(df$Treatment, "siLuc") #siLuc is always the first factor
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


find_half_life <- function(df){ # returns half life in minutes
    df_model = lm(log(Abundance) ~ Time, data=df)
    return(log(0.5)/coef(df_model)[2])
}

find_se <- function(df){
    df_summary <- summarySE(df, 
                            measurevar="Abundance", 
                            groupvars=c("Time", "Substrate", "Treatment"), 
                            na.rm=T)
    return(df_summary)
}

find_R2 <- function(df){
    return(cor(df$Time, df$Abundance)^2)
}

plot_overview <- function(df) { # plots all data in df as facet wrapped table
    ggplot(df, aes(x=Time, y=log(Abundance), color=Treatment)) +
        geom_point() + 
        facet_wrap(~ Substrate + Replicate) + 
        geom_smooth(method="lm", formula=(y~x), se=F)
}

plot_paper_curve_by_treatment <- function(df) { # plot curve formatted for paper
    df_summary <- summarySE(df, 
                            measurevar="Abundance", 
                            groupvars=c("Time", "Substrate", "Treatment"), 
                            na.rm=T)
    ggplot(df_summary, aes(x=Time, y=log(Abundance), color=Treatment)) +
        geom_point(size=3) + 
        geom_smooth(method="lm", formula=(y~x), size=1.3, se=F) + 
        geom_errorbar(aes(ymin=log(Abundance-se), ymax=log(Abundance+se), width=5)) +
        theme_minimal() +
        theme(panel.grid.major.y=element_blank(), 
                panel.grid.minor.x=element_blank(),
                panel.grid.major.x=element_blank(),
                panel.grid.minor.y=element_blank(),
                axis.ticks=element_line(),
                panel.border=element_rect(size=1, fill=NA, color='lightgrey'),
                text=element_text(size=10)) +
        scale_y_continuous(breaks=c(0.000001,-0.105360516, -0.223143551, -0.356674944, -0.510825624, -0.693147181, -0.916290732, -1.203972804, -1.609437912, -2.302585093),
                            label=(c(100,90,80,70,60,50,40,30,20,10))) +
        scale_x_continuous(breaks=c(0,120,240,360)) +
        scale_color_manual(values=c( "#A4A4A5","#000000")) +
        ylab("mRNA Abundance (%)") +
        xlab("Time (min)")
}
