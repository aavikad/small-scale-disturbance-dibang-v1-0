
#### pre and post disturbance ############
# DIBANG #

library(dplyr)
library(tidyverse)
library(devtools)
library(ggplot2)
library(ggpubr)

### exploration######
#### pre disturbance richness ###
bird_dist<- read.csv("disturbance_2021_2023.csv", header=TRUE)
bird_dist<- na.omit(bird_dist)
bird_dist<- bird_dist %>% distinct()
head(bird_dist)


## how many species are there in total ###
bird_dist_names = subset(bird_dist, select=c(disturbance, scientific_name))
unique(bird_dist_names$scientific_name) #73 species in total 

## how many species are there in pre disturbance vs post disturbance forests?##

pre_dist<- bird_dist %>% filter(disturbance== "pre-disturbance")
pre_dist_bird = subset(pre_dist, select=c(plot, scientific_name))

post_dist<- bird_dist %>% filter(disturbance== "post-disturbance")
post_dist_bird = subset(post_dist, select=c(plot, scientific_name))

aggregate(data=bird_dist_names,
          scientific_name ~ disturbance, function(x) length(unique(x))) # pre-disturbance 57, post-disturbance 61
intersect(pre_dist$scientific_name, post_dist$scientific_name) # 45 common 

setdiff(pre_dist$scientific_name, post_dist$scientific_name) #11 unique

setdiff(post_dist$scientific_name, pre_dist$scientific_name) #17 unique

### alpha diversity #######
###### species by sites matrix- incidence_raw data with 0s and 1s
library(iNEXT)
pre_dist_bird_matrix<- pre_dist_bird %>%
  group_by(plot, scientific_name) %>%
  summarise(n=n())

pre_dist_bird_matrix_wide<- pre_dist_bird_matrix %>%
  pivot_wider(names_from = plot, values_from = n, values_fill=0) %>%
  column_to_rownames("scientific_name")

post_dist_bird_matrix<- post_dist_bird %>%
  group_by(plot, scientific_name) %>%
  summarise(n=n())

post_dist_bird_matrix_wide<- post_dist_bird_matrix %>%
  pivot_wider(names_from = plot, values_from = n, values_fill=0) %>%
  column_to_rownames("scientific_name")

dist_list<- list(pre_disturbance=pre_dist_bird_matrix_wide, post_disturbance=post_dist_bird_matrix_wide)
richness.dist.raw<- iNEXT(dist_list, q=0, datatype = "incidence_raw", size=NULL, endpoint=NULL, se= TRUE, knots=40, conf=0.95, nboot=50)
richness.dist.raw$iNextEst$coverage_based
plot(richness.dist.raw, type=1, cex.axis = 1.5) # sample size-based R/E curve
plot(richness.dist.raw, type=2, cex.axis = 1.5) # completeness curve
plot(richness.dist.raw, type=3, cex.axis = 1.5) # coverage-based rarefaction/extrapolation curve


#rarefied_richness_coverage<- richness.dist.raw2$iNextEst$coverage_based # To Extract asymptotic richness estimates (order = 0)

## some more inext plots (optional)
ggiNEXT(richness.dist.raw, type=1)+
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 16),
        axis.text = element_text(size = 16)) # sample size-based R/E curve

ggiNEXT(richness.dist.raw, type=2)+ 
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 16),
        axis.text = element_text(size = 16)) # completeness curve

ggiNEXT(richness.dist.raw, type=3)+ theme_bw() +
 theme(legend.position = "right",
       text = element_text(size = 16),
       axis.text = element_text(size = 16)) # coverage-based rarefaction/extrapolation curve

#lapply(richness.dist.raw, function(x) write.table( data.frame(x), 'richness_rarefied.csv'  , append= T, sep=',' )) # to export
rarefied_richness_coverage<- read.csv("rarefied_richness_coverage.csv") # cleaned

pre_dist_richness<- rarefied_richness_coverage %>%
  filter(Assemblage== "pre_dist")

post_dist_richness<- rarefied_richness_coverage %>%
  filter(Assemblage=="post_dist")

ggqqplot(pre_dist_richness$qD)
shapiro.test(pre_dist_richness$qD) # near normal

ggqqplot(post_dist_richness$qD)
shapiro.test(post_dist_richness$qD) # not normal

wilcox.test(pre_dist_richness$qD, post_dist_richness$qD) #W = 151, p-value = 0.1918 not significantly different

### alpha richness of one trophic group ## sub-canopy invertivore ########
library(vegan)

sub_canopy_pre<- pre_dist %>%
  filter(trophic_niche_stratification== "sub-canopy-invertivore")

sub_canopy_post<- post_dist %>%
  filter(trophic_niche_stratification== "sub-canopy-invertivore")
  
pre_dist_sub_canopy<- sub_canopy_pre %>%
     group_by(scientific_name, plot) %>%
   summarise(n=n())

pre_dist_sub_canopy_wide<- pre_dist_sub_canopy %>%
   pivot_wider(names_from = scientific_name, values_from = n, values_fill=0) %>%
   column_to_rownames("plot")

post_dist_sub_canopy<- sub_canopy_post %>%
  group_by(scientific_name, plot) %>%
  summarise(n=n())

post_dist_sub_canopy_wide<- post_dist_sub_canopy %>%
  pivot_wider(names_from = scientific_name, values_from = n, values_fill=0) %>%
  column_to_rownames("plot")


# First Define a function to calculate alpha richness with first-order Jackknife for each site
species_richness <- function(site_data) {
  sum(site_data > 0)
}


alpha_richness_jackknife_first_order <- function(dataset) {
  num_sites <- nrow(dataset)
  alpha_richness_jackknife_first_order <- numeric(num_sites)
  
  for (i in 1:num_sites) {
    # Exclude the current site
    reduced_data <- dataset[-i, , drop = FALSE]
    
    # Calculate species richness for the remaining sites
    site_richness <- apply(reduced_data, 1, species_richness)
    
    # Calculate the first-order Jackknife estimate for the current site
    alpha_richness_jackknife_first_order[i] <- mean(site_richness)
  }
  
  return(alpha_richness_jackknife_first_order)
}

#Now Define a function to calculate alpha richness with second-order Jackknife for each site
alpha_richness_jackknife_second_order <- function(dataset) {
  num_sites <- nrow(dataset)
  alpha_richness_jackknife_first_order <- alpha_richness_jackknife_first_order(dataset)
  
  # Calculate the mean of the first-order Jackknife estimates
  mean_first_order <- mean(alpha_richness_jackknife_first_order)
  
  # Calculate the second-order Jackknife estimate
  alpha_richness_jackknife_second_order <- (2 * mean_first_order) - alpha_richness_jackknife_first_order
  
  return(alpha_richness_jackknife_second_order)
}

# Calculate alpha richness with second-order Jackknife for pre_dist and post_dist sampling periods
sub_canopy_pre_dist <- alpha_richness_jackknife_second_order(pre_dist_sub_canopy_wide)
sub_canopy_post_dist <- alpha_richness_jackknife_second_order(post_dist_sub_canopy_wide)

# Print the second-order Jackknife estimates for pre_dist and post_dist sampling periods
print(sub_canopy_pre_dist)
print(sub_canopy_post_dist)

sub_canopy<- read.csv("sub_canopy.csv") # has combined jackknife richness and vegetation variables
sub_canopy$plot<- factor(sub_canopy$plot)
sub_canopy$sampling_stage_two<- factor(sub_canopy$sampling_stage, levels= c("pre-disturbance", "post-disturbance"))

ggplot(sub_canopy, aes(x = plot, y = jackknife, color=sampling_stage_two)) +
  geom_point(size=3.5, alpha=0.5) + 
  labs(x = "Forest Plots", y = "Richness", color="Sampling stage") + 
  ggtitle("Richness of sub canopy invertivores across sampling stages") +
  theme_bw()+
  scale_color_manual(values = c("#DF536B", "#28E2E5"))+
  theme(text = element_text(size = 14))+
  theme(
    text = element_text(size = 14),  
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 14) 
  )

### beta- diversity (community) #########
library(betapart)

all_bird_dist<- bird_dist %>%
  group_by(scientific_name, plot) %>%
  summarise(n=n())

all_bird_dist_wide<- all_bird_dist %>%
  pivot_wider(names_from = scientific_name, values_from = n, values_fill=0) %>%
  column_to_rownames("plot")

all_birds_dist_mat<- vegdist(all_bird_dist_wide, method="jaccard") ## calculate distance matrix - dissimilarity ##
all_birds_dist_mat1 <- as.matrix(all_birds_dist_mat, labels =T) ## create matrix

nmds_dist<- metaMDS(all_birds_dist_mat1, distance= "jaccard", k=2, maxit= 999, trymax=500, wascores= TRUE, autotransform = FALSE)
goodness(nmds_dist) ### goodness of fit ###
stressplot(nmds_dist)

nmds_scores<- as.data.frame(scores(nmds_dist)) # extract scores of nmds
nmds_scores$sites<- rownames(nmds_scores)

plot(nmds_dist, type= "t") # display empty ordination space
orditorp(nmds_dist, display = "sites", cex=1)


nmds_scores <- nmds_scores %>%
  mutate(grp = ifelse(grepl("5$", sites), "post-disturbance",
                      ifelse(grepl("6$", sites), "post-disturbance",
                             ifelse(grepl("7$", sites), "post-disturbance",
                                    ifelse(grepl("8$", sites), "post-disturbance", 
                                           ifelse(grepl("1$", sites), "pre-disturbance",
                                                  ifelse(grepl("2$", sites), "pre-disturbance",
                                                         ifelse(grepl("3$", sites), "pre-disturbance",
                                                                ifelse(grepl("4$", sites), "pre-disturbance", NA)))))))))

grp.a <- nmds_scores[nmds_scores$grp == "post-disturbance", ][chull(nmds_scores[nmds_scores$grp == 
                                                                                  "post-disturbance", c("NMDS1", "NMDS2")]), ]  # hull values for grp A

grp.b <- nmds_scores[nmds_scores$grp == "pre-disturbance", ][chull(nmds_scores[nmds_scores$grp == 
                                                                                 "pre-disturbance", c("NMDS1", "NMDS2")]), ]

hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data


ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) + # add the convex hulls
  geom_point(data=nmds_scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=4) + # add the point markers
  scale_colour_manual(values=c("pre-disturbance" = "#DF536B", "post-disturbance" = "#28E2E5")) +
  scale_shape_manual(values=c(1, 2))+
  scale_fill_manual(values=c("pre-disturbance" = "#DF536B" , "post-disturbance" = "#28E2E5")) +
  coord_equal() +
  theme_bw()


# Run adonis on jaccard dissimilarity matrix
adonis_result <- adonis2(all_birds_dist_mat ~ nmds_scores$grp) # 0.054  not significant
print(adonis_result)


####### trait analysis######
#### visualsing ### 
### box plots ###
traits<- read.csv("traits.csv")
traits<- na.omit(traits)
traits$disturbance<- factor(traits$disturbance, levels=c("pre-disturbance", "post-disturbance"))
traits$logmass<- log(traits$mass)

ggplot(traits, aes(x=disturbance, y=logmass, color=disturbance), col="white")+
  geom_boxplot()+
  geom_jitter()+ theme_bw() + xlab("disturbance periods")

library(xfun)
library(hrbrthemes)

### density plot ###
suppressWarnings({ggplot(traits, aes(x=logmass, fill= disturbance, color=disturbance))+
    geom_density(alpha=0.3)+
    theme_ipsum()+ 
    labs(x="log(mass)", y= "density")+
    theme(
      axis.title.x = element_text(vjust = 0),  # Adjust vjust to center horizontally
      axis.title.y = element_text(vjust = 0.5, angle = 90)  # Adjust vjust to center vertically
    )
})

traits_pre<- traits %>%
  filter(disturbance== "pre-disturbance")

traits_post<- traits %>%
  filter(disturbance== "post-disturbance")

ggqqplot(traits_pre$logmass)
shapiro.test(traits_pre$logmass) # W = 0.92131, p-value = 6.112e-05 not normal

ggqqplot(traits_post$logmass)
shapiro.test(traits_post$logmass) # W = 0.91813, p-value = 1.352e-05 not normal

wilcox.test(traits_pre$logmass, traits_post$logmass) #W = 4551.5, p-value = 0.3497 not significant
fligner.test(list(traits_pre$logmass, traits_post$logmass)) #Fligner-Killeen:med chi-squared = 1.7697, df = 1, p-value = 0.1834 # not significant

all_trophic<- union(traits_pre$trophic_niche_stratification, traits_post$trophic_niche_stratification)
table_pre <- table(factor(traits_pre$trophic_niche_stratification, levels = all_trophic))
table_post <- table(factor(traits_post$trophic_niche_stratification, levels = all_trophic))
combined_trophic <- cbind(table_pre, table_post)
fisher_trophic<- fisher.test(combined_trophic, simulate.p.value = TRUE, B = 10000)
print(fisher_trophic) # not significant p-value =  0.3057 not significant
summary(fisher_trophic)


# Calculate the test statistic manually
n <- sum(combined_trophic)  # Total number of observations
a <- combined_trophic[1, 1]  # Frequency in cell 1,1
b <- combined_trophic[1, 2]  # Frequency in cell 1,2
c <- combined_trophic[2, 1]  # Frequency in cell 2,1
d <- combined_trophic[2, 2]  # Frequency in cell 2,2

# Calculate Fisher's exact test statistic
test_statistic <- (a * d - b * c) / sqrt((a + b) * (a + c) * (b + d) * (c + d))

trophic_table<- table(traits$trophic_niche_stratification,traits$disturbance)
fisher.test(trophic_table) # p-value = 0.8776


library(RColorBrewer)

pastel_palette <- c("#FAC08F", "#FFDE8A", "#FFF1AA", "#C0F8D1", "#A4E3D4",
                    "#B9D1F1", "#B3B7F5", "#E3B4D1", "#F9C4B2", "#FFD2AA")
par(cex.axis = 1.0)
barplot(prop.table(trophic_table, 2),
        legend.text = TRUE,
        args.legend = list(cex=0.80,x="topright"),
        col = pastel_palette,
        ylab = "Proportion of Trophic groups",
        xlab = "Sampling stages",
        cex.axis=1.2,
        cex.lab = 1.2)


### indicator species analysis ##### ======

library(indicspecies)

## pre-disturbance period ###

pre_dist_bird_matrix2<- pre_dist_bird %>%
  group_by(scientific_name, plot) %>%
  summarise(n=n())

pre_dist_bird_wide2<- pre_dist_bird_matrix2 %>%
  pivot_wider(names_from = scientific_name, values_from = n, values_fill=0) %>%
  column_to_rownames("plot")

pre_mat<- as.data.frame(ifelse(pre_dist_bird_wide2 >0,1,0))

# Move row names to the last column
pre_mat2 <- cbind(pre_mat[, -1, drop = FALSE], ID = rownames(pre_mat))
colnames(pre_mat2)[ncol(pre_mat2)]<- "Group"
rownames(pre_mat2)<- NULL

# Perform IndVal analysis
indval_result_pre <- multipatt(pre_mat2[, -ncol(pre_mat2)], pre_mat2$Group, func= "r",
                               control = how(nperm = 999))

summary(indval_result_pre)

# Extract indicator values and their significance
indicator_values_pre <- indval_result_pre$sign
p_values_pre<- indicator_values_pre$p.value

# Display results
ind_pre <- data.frame(
  Species = colnames(pre_mat2)[1:(ncol(pre_mat2) - 1)],
  Indicator_Value = indicator_values_pre,
  P_value= p_values_pre
)

View(ind_pre) # ubiquitous species have NA values

#write.csv(ind_pre, "ind_pre.csv")

### post-disturbance period ##

post_dist_bird_matrix2<- post_dist_bird %>%
  group_by(scientific_name, plot) %>%
  summarise(n=n())

post_dist_bird_wide2<- post_dist_bird_matrix2 %>%
  pivot_wider(names_from = scientific_name, values_from = n, values_fill=0) %>%
  column_to_rownames("plot")

post_mat<- as.data.frame(ifelse(post_dist_bird_wide2 >0,1,0))

# Move row names to the last column
post_mat2 <- cbind(post_mat[, -1, drop = FALSE], ID = rownames(post_mat))
colnames(post_mat2)[ncol(post_mat2)]<- "Group"
rownames(post_mat2)<- NULL

# Perform IndVal analysis
indval_result_post <- multipatt(post_mat2[, -ncol(post_mat2)], func= "r", post_mat2$Group, control = how(nperm = 999))

View(indval_result_post)

# Extract indicator values and their significance
indicator_values_post <- indval_result_post$sign
p_values_post<- indicator_values_post$p.value

# Display results
ind_post <- data.frame(
  Species = colnames(post_mat2)[1:(ncol(post_mat2) - 1)],
  Indicator_Value = indicator_values_post,
  P_value= p_values_post
)
View(ind_post)
#write.csv(ind_post, "ind_post.csv")

#### vegetation assessment ########
birds_veg<- read.csv("vegetation.csv")

veg_pre<- birds_veg %>%
  filter(sampling_period== "pre-disturbance")

veg_post<- birds_veg %>%
  filter(sampling_period == "post-disturbance")

foliage_columns_pre <- veg_pre[, 11:13]
foliage_columns_post <- veg_post[, 11:13]

shannon_wiener <- function(proportions) {
  -sum(proportions * log(proportions))
}

veg_pre$Shannon_Wiener_Index <- apply(foliage_columns_pre, 1, shannon_wiener)
veg_post$Shannon_Wiener_Index <- apply(foliage_columns_post, 1, shannon_wiener)


min_max_scaling <- function(x) {
   (x - min(x)) / (max(x) - min(x))
 }

veg_pre$FHD <- min_max_scaling(veg_pre$Shannon_Wiener_Index)
veg_post$FHD <- min_max_scaling(veg_post$Shannon_Wiener_Index)
veg_both<- rbind(veg_pre, veg_post)
View(veg_both)

#write.csv(veg_both, "veg_both.csv") # to export

## vegetation data exploration ##
t.test(veg_pre$canopy_total, veg_post$canopy_total) # t = 10.842, df = 34.557, p-value = 1.152e-12 significant
t.test(veg_pre$leaf_litter_total, veg_post$leaf_litter_total) # t = 10.348, df = 37.874, p-value = 1.36e-12 significant
ks.test(veg_pre$FHD, veg_post$FHD) # D = 0.4, p-value = 0.07342 not significant

### visualising data ###

veg_both$sampling_period<- factor(veg_both$sampling_period, levels=c("pre-disturbance", "post-disturbance"))

ggplot(veg_both, aes(x=sampling_period, y=canopy_total, color=sampling_period), col="white")+
  geom_boxplot()+
  geom_jitter()+ theme_bw() + xlab("Sampling_period")

ggplot(veg_both, aes(x=sampling_period, y=leaf_litter_total, color=sampling_period), col="white")+
  geom_boxplot()+
  geom_jitter()+ theme_bw() + xlab("Sampling_period")


ggplot(veg_both, aes(x = leaf_litter_total, y = canopy_total)) +
  geom_point() +
  labs(title = "Leaf litter % and canopy cover %",
       x = "Leaf Litter Total",
       y = "canopy_total")


##### GLMMs#########
library(Matrix)
library(lme4)
library(lmerTest)
library(DHARMa)
library(AICcmodavg)
library(ggpubr)
library(MASS)

## checking for correlation ##
veg_all<- read.csv("veg_both_rich.csv") # vegetation and species richness
corr1 <- cor.test(veg_all$canopy_total, veg_all$leaf_litter_total)
print(corr1) #83%, p-value = 1.936e-11

corr2<- cor.test(veg_all$canopy_total, veg_all$FHD)
print(corr2) #t = -1.0121, df = 38, p-value = 0.5867 # not significant

corr3<- cor.test(veg_all$leaf_litter_total, veg_all$FHD)
print(corr3) #0.08%,  p-value = 0.9698 ## not significant

library(reshape2)
corr_matrix<- cor(veg_all[c("canopy_total", "leaf_litter_total", "FHD")])
corr_data <- melt(corr_matrix)

# Plot the heatmap with correlation values
ggplot(corr_data, aes(Var1, Var2, fill = value, label = ifelse(abs(round(value, 2)) >= 0.5, paste0(round(value, 2), "*"), ""))) +
  geom_tile(color = "white") +
  geom_text(size = 6, color = "black") +  # Add correlation values
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limit = c(-1,1)) +
  theme_minimal() +
  labs(title = "Correlation Chart of Vegetation Variables",
       x = "Variables",
       y = "Variables")+
  theme(axis.text.x = element_text(size = 14, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 14, margin = margin(t = 40, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(size = 16, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 16, margin = margin(t = 30, r = 0, b = 0, l = 0)),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14))


### using alpha richness (alpha_rarefied) as response variable ###
veg_all$sampling_stage<- factor(veg_all$sampling_stage)
veg_all$sampling_stage<- relevel(veg_all$sampling_stage, ref = "pre-disturbance") #pre-disturbance as a reference 
veg_all$plots<- factor(veg_all$plot)

### null model ###

model_0<- lmerTest::lmer(alpha_rarefied ~ 1+ (1|plot), data = veg_all)
summary(model_0)

model_1<- lmerTest::lmer(alpha_rarefied ~ (FHD+canopy_total)*sampling_stage+ (1|plot),
                         data=veg_all) # FHD + canopy 
summary(model_1)


model_2 <-lmerTest::lmer(alpha_rarefied ~ (FHD + leaf_litter_total)*sampling_stage+ (1|plot),
                data = veg_all) # FHD + leaf litter
summary(model_2)

model_2_res<- residuals(model_2)
qqnorm(model_2_res)
qqline(model_2_res)

residuals_2 <- simulateResiduals(model_2) #null hypothesis is that there is no overdispersion
plot(residuals_2)
testDispersion(residuals_2) #dispersion = 0.9521, p-value = 0.904

model_3 <-lmerTest::lmer(alpha_rarefied ~ (FHD)*sampling_stage+ (1|plot),
                         data = veg_all) #FHD
summary(model_3)

model_4 <-lmerTest::lmer(alpha_rarefied ~ (canopy_total)*sampling_stage+ (1|plot),
                         data = veg_all) # canopy 
summary(model_4)

model_5 <-lmerTest::lmer(alpha_rarefied ~ (leaf_litter_total)*sampling_stage+ (1|plot),
                         data = veg_all) # leaf litter
summary(model_5)

model_set1 <- list(model_0, model_1, model_2, model_3, model_4, model_5)
mod.names1<- c('null', 'FHD+canopy', 'FHD+leaf litter', 'FHD', 'canopy', 'leaf litter')
aictab(cand.set = model_set1, modnames = mod.names1) # model_1 with FHD performs best

# Model selection based on AICc:
#   
#   K   AICc Delta_AICc AICcWt Cum.Wt  Res.LL
# FHD+leaf litter 8 265.23       0.00   0.98   0.98 -122.29
# FHD+canopy      8 273.90       8.66   0.01   0.99 -126.63
# FHD             6 275.05       9.82   0.01   1.00 -130.25
# leaf litter     6 280.50      15.27   0.00   1.00 -132.98
# canopy          6 298.87      33.64   0.00   1.00 -142.16
# null            3 307.56      42.33   0.00   1.00 -150.45


# model_2 plot #
ggplot(veg_all, aes(x = FHD, y = alpha_rarefied, color = factor(sampling_stage))) +
  geom_smooth(method = "lm", se = TRUE, size = 1.5, alpha = 0.1) +  # Add best fit line
  labs(x = "FHD", y = "Alpha Rarefied",  color = "Sampling Stage") +
  theme_minimal()+
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),  
        legend.title = element_text(size = 14),
  legend.text= element_text(size=14))+
  ggtitle("Species richness (alpha rarefied) and FHD across sampling stages")

### using alpha richness (sub_canopy) as response variable ###
sub_canopy$sampling_stage<- factor(sub_canopy$sampling_stage)
sub_canopy$sampling_stage<- relevel(sub_canopy$sampling_stage, ref = "pre-disturbance")

## null model
library(glmmTMB)


model_00 <- glmmTMB(jackknife ~ 1+ (1|plot),
                          data = sub_canopy)
summary(model_00)

model_6<- glmmTMB(jackknife ~ (FHD + leaf_litter_total)* sampling_stage + (1|plot), 
                  data=sub_canopy) #FHD + leaf litter
summary(model_6)

model_6_res<- residuals(model_6)
qqnorm(model_6_res)
qqline(model_6_res)

residuals_6 <- simulateResiduals(model_6) #null hypothesis is that there is no overdispersion
plot(residuals_6)
testDispersion(residuals_6) #dispersion = 0.92693, p-value = 0.776

model_7<- glmmTMB(jackknife ~ (FHD+canopy_total)*sampling_stage + (1|plot), 
               data=sub_canopy) # FHD + canopy
summary(model_7)


model_8<- glmmTMB(jackknife ~ (FHD)* sampling_stage + (1|plot), 
                         data=sub_canopy) #FHD
summary(model_8)


model_9<- glmmTMB(jackknife ~ (canopy_total)* sampling_stage + (1|plot), 
                         data=sub_canopy) #canopy cover
summary(model_9)

model_10<- glmmTMB(jackknife ~ (leaf_litter_total)* sampling_stage + (1|plot), 
                         data=sub_canopy) #leaf litter
summary(model_10)


model_set2 <- list(model_00, model_6, model_7, model_8, model_9, model_10)
mod.names2<- c('null', 'FHD+leaf litter', 'FHD+canopy', 'FHD', 'canopy', 'leaf litter')
aictab(cand.set = model_set2, modnames = mod.names2) # model 6 with FHD + leaf litter performs the best

# Model selection based on AICc:
#   
#   K   AICc Delta_AICc AICcWt Cum.Wt     LL
# FHD+leaf litter 8 -57.84       0.00   0.95   0.95  39.24
# leaf litter     6 -50.51       7.33   0.02   0.97  32.53
# FHD             6 -50.28       7.56   0.02   1.00  32.41
# FHD+canopy      8 -46.82      11.01   0.00   1.00  33.73
# canopy          6 -42.14      15.69   0.00   1.00  28.34
# null            3 100.87     158.71   0.00   1.00 -47.10


# model_4 plot #
ggplot(sub_canopy, aes(x = FHD, y = jackknife, color = factor(sampling_stage))) +
  geom_smooth(method = "lm", se = TRUE, size = 1.5, alpha = 0.1) +  # Add best fit line
  labs(x = "FHD", y = "Alpha Rarefied",  color = "Sampling Stage") +
  theme_minimal()+
  theme(axis.text = element_text(size = 14),   # Increase axis label font size
        axis.title = element_text(size = 14),  # Increase axis title font size
        legend.title = element_text(size = 14))+
  ggtitle("Species richness (alpha sub canopy) and FHD across sampling stages")

