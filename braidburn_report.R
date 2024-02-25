#### Environmental Surveying and Monitoring
### Freshwater quality and community composition report - Braid Burn

rm(list = ls())

install.packages("devtools")
install.packages("ggdist")
devtools::install_github("gavinsimpson/ggvegan")
library(ggplot2)
library(vegan)
library(ggvegan)
library(ggdist)

braidburn_data <- read.csv("braidburn_biota.csv")
#maybe dont have to do the below. could just repeat analysis using only [, 1:4] and [,5:7]
LNR <- c(rep("0", 4), rep("1", 3)) #sites 1-4 are not part of a local nature reserve (LNR). Adding this to the dataset. 
braidburn_data.2 <- data.frame(
  braidburn_data[, 1], #copying first existing column
  LNR = LNR, #adding in the new column
  braidburn_data[ , 2:ncol(braidburn_data)] #adding the rest of the columns
) 

View(braidburn_data)


################################
#Biodiversity indices
burn.simp <- diversity(braidburn_data,index = "simpson")  #Simpson's index 
burn.simp
par(mfrow = c(1,1))
hist(burn.simp, 
     main = paste("Simpson's Biodiversity Index"),
     xlab = "Simpson's index value",
     ylab = "Frequency",
     col = "#006600")
#for just sites not included in the LNR
burn.simp.nlnr <- diversity(braidburn_data[1:4, ],index = "simpson")  #Simpson's index 
burn.simp.nlnr
hist(burn.simp.nlnr, 
     main = paste("Simpson's Biodiversity Index (Not LNR)"),
     xlab = "Simpson's index value",
     ylab = "Frequency",
     col = "#006600")
#for just sites included in the LNR
burn.simp.lnr <- diversity(braidburn_data[5:7, ],index = "simpson")  #Simpson's index 
burn.simp.lnr
hist(burn.simp.lnr, 
     main = paste("Simpson's Biodiversity Index (LNR)"),
     xlab = "Simpson's index value",
     ylab = "Frequency",
     col = "#006600")
#probably better to display the data in a single violin plot rather than histograms 
#creating a group dataframe for these 3 sets of indices
# Create a data frame
simp.diversity_data <- data.frame(
  Group = c(rep("LNR", length(burn.simp.lnr)),
            rep("Non-LNR", length(burn.simp.nlnr)),
            rep("Whole Dataset", length(burn.simp))),
  Simpson = c(burn.simp.lnr, burn.simp.nlnr, burn.simp)
)
#creating a violin plot with these indices
simpsons_plot <- ggplot(simp.diversity_data, aes(x = Group, y = Simpson, fill = Group)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.7) +
  geom_jitter(aes(color = Group), width = 0.2, size = 0.7, alpha = 0.7) +
  labs(x = "Group",
       y = "Simpson's Diversity Index") +
  theme_minimal() +
  scale_fill_manual(values = c("LNR" = "#006600", "Non-LNR" = "#33FF00", "Whole Dataset" = "#669966")) +
  scale_color_manual(values = c("LNR" = "#006600", "Non-LNR" = "#33FF00", "Whole Dataset" = "#669966"))



#Shannon's index 
burn.shan <- diversity(braidburn_data,index = "shannon")
burn.shan
hist(burn.shan,
     main = paste("Shannon's Diversity Histogram"),
     xlab = "Shannon's index value",
     ylab = "Frequency",
     col = "#006600")
#just the sites not included in the LNR
burn.shan.nlnr <- diversity(braidburn_data[1:4, ],index = "shannon")
burn.shan.nlnr
hist(burn.shan.nlnr,
     main = paste("Shannon's Diversity Histogram (Not LNR)"),
     xlab = "Shannon's index value",
     ylab = "Frequency",
     col = "#006600")
#jsut the sites included in the LNR
burn.shan.lnr <- diversity(braidburn_data[5:7,],index = "shannon")
burn.shan.lnr
hist(burn.shan.lnr,
     main = paste("Shannon's Diversity Histogram (LNR)"),
     xlab = "Shannon's index value",
     ylab = "Frequency",
     col = "#006600")

#creating a group dataframe for these 3 sets of indices
# Create a data frame
shan.diversity_data <- data.frame(
  Group = c(rep("LNR", length(burn.shan.lnr)),
            rep("Non-LNR", length(burn.shan.nlnr)),
            rep("Whole Dataset", length(burn.shan))),
  Shannon = c(burn.shan.lnr, burn.shan.nlnr, burn.shan)
)
#creating a violin plot with these indices
shannon.plot <- ggplot(shan.diversity_data, aes(x = Group, y = Shannon, fill = Group)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.7) +
  geom_jitter(aes(color = Group), width = 0.2, size = 0.7, alpha = 0.7) +
  labs(x = "Group",
       y = "Shannon's Diversity Index") +
  theme_minimal() +
  scale_fill_manual(values = c("LNR" = "#006600", "Non-LNR" = "#33FF00", "Whole Dataset" = "#669966")) +
  scale_color_manual(values = c("LNR" = "#006600", "Non-LNR" = "#33FF00", "Whole Dataset" = "#669966"))

# Look at the differences between the two indices using Bray Curtis dissimilarity 
bray = vegdist(braidburn_data, "bray") 
hist(bray, xlim = range(0.0,1.0),
     main = paste("Bray Curtis Dissimilarity of Plots"),
     xlab = "Dissimilarity Value",
     ylab = "Frequency",
     col = "#33FF66")

#displaying both plots using patchwork package
install.packages("patchwork")
library("patchwork")

combined_plots <- shannon.plot/simpsons_plot
combined_plots

# Rarefaction Curves
rarecurve(braidburn_data, col = "#33FF66")

## Count the number of individuals in each sample
burn.abund<- rowSums(braidburn_data)  #gives the number of individuals found in each plot
burn.abund

##############################
## Species Accumulation Curve 
#Looking at the entire dataset and the LNR inclusions/exclusions
# Create Not LNR subset (rows 1-4)
not_lnr_data <- braidburn_data[1:4, ]

# Create LNR subset (rows 5-7)
lnr_data <- braidburn_data[5:7, ]

whole_dataset_accum <- specaccum(braidburn_data, method = "random")

# Run specaccum for Not LNR
not_lnr_accum <- specaccum(not_lnr_data, method = "random")

# Run specaccum for LNR
lnr_accum <- specaccum(lnr_data, method = "random")

# Combine the data for plotting
combined_data <- rbind(
  data.frame(sites = 1:length(whole_dataset_accum$richness),
             richness = whole_dataset_accum$richness,
             sd = whole_dataset_accum$sd,
             group = "All Sites"),
  data.frame(sites = 1:length(not_lnr_accum$richness),
             richness = not_lnr_accum$richness,
             sd = not_lnr_accum$sd,
             group = "Not LNR"),
  data.frame(sites = 1:length(lnr_accum$richness),
             richness = lnr_accum$richness,
             sd = lnr_accum$sd,
             group = "LNR")
)

# Create the ggplot
ggplot(combined_data, aes(x = sites, y = richness, group = group, color = group, fill = group)) +
  geom_line(size = 2) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), alpha = 0.5) +
  labs(x = "Number of Sites",
       y = "Number of Species") +
  scale_color_manual(values = c("#FF6600", "#FFCC00", "#CC3333")) +
  scale_fill_manual(values = c("#FF6600", "#FFCC00", "#CC3333")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),   #axis labels
        axis.title = element_text(size = 14),  #axis titles
        axis.text.x = element_text(size = 10), #x-axis tick labels
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) 



#####################################
## Beta diversity indices
install.packages("adespatial")
library(adespatial)

## Jaccard dissimilarity
burn.div.jac <- beta.div.comp(braidburn_data, coef = "J", quant = T)
burn.div.jac$part

## Sorensen dissimilarity 
burn.div.sor <- beta.div.comp(braidburn_data, coef = "S", quant = T)
burn.div.sor$part

#### Individual Locations
## Identify the local contribution to Beta-diversity 
# For Jaccard replacement/turnover
local.burn.repl <- LCBD.comp(burn.div.jac$repl, 
                             sqrt.D = T)
local.burn.repl

plot(local.burn.repl$LCBD, 
     type = "p", 
     ylab = "Contribution to beta-diversity",
     xlab = "Site Number")


# For Jaccard richness/turnover
local.burn.rich <- LCBD.comp(burn.div.jac$rich, 
                             sqrt.D = T)
local.burn.rich

plot(local.burn.rich$LCBD, 
     type = "p", 
     ylab = "Contribution to beta-diversity",
     xlab = "Site Number")

## Species contribution to Beta-diversity #
burn.scBD <- beta.div(braidburn_data, method = "hellinger")
burn.scBD

barplot(burn.scBD$SCBD, 
        horiz = TRUE,
        col="#FF66CC")


#########################################
## Non-Metric Multidimensional Scaling (NMDS)
## Perform the initial analysis ##
burn.nmds <- metaMDS(braidburn_data,
                     k = 2,     
                     distance = "bray",
                     maxit = 999,
                     trymax = 500,
                     wascores = TRUE,
)

burn.nmds

## Shepard's stressplot
stressplot(burn.nmds)
#NMDS plot
autoplot(burn.nmds, geom = "text", label.size = 4)+
  scale_color_manual(values = c("sites" = "#FF3399", "species" = "#330066")) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),   #axis labels
        axis.title = element_text(size = 14),  #axis titles
        axis.text.x = element_text(size = 10), #x-axis tick labels
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))


#############
# Separating the datapoints by the site area - in an LNR or not. 
#don't think this is actually useful information - will omit from report
# Create a data frame for plotting
plot_data <- data.frame(NMDS1 = burn.nmds$points[, 1],
                        NMDS2 = burn.nmds$points[, 2],
                        LNR = factor(braidburn_data.2$LNR, labels = c("Not LNR", "LNR")))

# Create the ggplot
ggplot(plot_data, aes(x = NMDS1, y = NMDS2, color = LNR, shape = LNR)) +
  geom_point(size = 4) +
  geom_segment(aes(x = NMDS1, y = NMDS2,
                   xend = lead(NMDS1),
                   yend = lead(NMDS2)),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c("#0073C2", "#CC99FF"), name = "LNR") +
  scale_shape_manual(values = c(18, 20), name = "LNR") +
  labs(title = "Groupings by LNR inclusion",
       x = "NMDS1",
       y = "NMDS2") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),   #axis labels
        axis.title = element_text(size = 14),  #axis titles
        axis.text.x = element_text(size = 10), #x-axis tick labels
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))



## Test for differences in those groups ##

LNR.diff <- adonis(braidburn_data~LNR,
                     data = braidburn_data.2,
                     permutations = 1000,
                     method = "bray")
LNR.diff

## Need to check for equal variances ##
dist.burn <- vegdist(braidburn_data)
burn.aov <- anova(betadisper                # Are there differences between the sites surveyed within the ANOVA?
                  (dist.burn, braidburn_data.2$LNR))
burn.aov  #OK to run the Adonis test - looking for non-significance


############
#WHPT score 


