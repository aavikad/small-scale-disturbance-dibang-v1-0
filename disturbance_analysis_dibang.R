
## Aavika Dhanda, Sonya Clegg
## The effects of small-scale disturbance on avian communities of an eastern Himalayan tropical forest ##

library(dplyr)
library(tidyverse)
library(devtools)
library(ggplot2)
library(ggpubr)
library(git2r)
library(usethis)

### PART I ###
####### trait analysis###### ==========

#### visualsing ### 
### box plots ###
traits<- read.csv("traits.csv")
traits<- na.omit(traits)
traits$disturbance<- factor(traits$disturbance, levels=c("pre-disturbance", "post-disturbance"))
traits$logmass<- log(traits$mass)

ggplot(traits, aes(x=disturbance, y=logmass, color=disturbance), col="white")+
  geom_boxplot()+
  geom_jitter()+ theme_bw() + xlab("disturbance periods") # exploration

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
print(fisher_trophic) # not significant
summary(fisher_trophic)

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


### indicator species analysis ##### 

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

#write.csv(ind_pre, "ind_pre.csv") # to export

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

View(indval_result_post) # ubiquitous species have NA values

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
#write.csv(ind_post, "ind_post.csv")# to export

### PART II ###
#### pre disturbance richness ###
bird_dist<- read.csv("sp_disturbance_dibang_2021_2023.csv", header=TRUE)
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

setdiff(post_dist$scientific_name, pre_dist$scientific_name) #16 unique

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
richness.dist.raw$iNextEst$coverage_based # alpha estimates
plot(richness.dist.raw, type=1, cex.axis = 1.5) # sample size-based R/E curve
plot(richness.dist.raw, type=2, cex.axis = 1.5) # completeness curve
plot(richness.dist.raw, type=3, cex.axis = 1.5) # coverage-based rarefaction/extrapolation curve

#rarefied_richness_coverage<- richness.dist.raw2$iNextEst$coverage_based # To Extract asymptotic richness estimates (order q= 0)

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
rarefied_richness_coverage<- read.csv("rarefied_richness_coverage.csv") # cleaned ?

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
sub_canopy$sampling_stage<- factor(sub_canopy$sampling_stage, levels= c("pre-disturbance", "post-disturbance"))

ggplot(sub_canopy, aes(x = plot, y = jackknife, color=sampling_stage)) +
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


#ordiplot(nmds_dist, type= "n", main="ellipses")

#ordiellipse(nmds_dist, groups = bird_dist$dis, draw= "polygon")

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


plot(veg_both$canopy_total, veg_both$leaf_litter_total, xlab = "Canopy Cover", ylab = "Leaf Litter", main = "Scatter Plot")

install.packages("vcd")
library(vcd)
assocstats(table(veg_all$canopy_total, veg_all$sampling_stage))


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

## checking for correlation ##
veg_all<- read.csv("veg_both_rich.csv")
corr1 <- cor.test(veg_all$canopy_total, veg_all$leaf_litter_total)
print(corr1) #t = 9.3845, df = 38, p-value = 1.936e-11

corr2<- cor.test(veg_all$canopy_total, veg_all$FHD)
print(corr2) #t = -1.0121, df = 38, p-value = 0.3179 # not significant

corr3<- cor.test(veg_all$leaf_litter_total, veg_all$FHD)
print(corr3) #t = -3.7138, df = 38, p-value = 0.0006534


### using alpha richness (alpha_rarefied) as response variable ###
veg_all$sampling_stage<- factor(veg_all$sampling_stage)
veg_all$sampling_stage<- relevel(veg_all$sampling_stage, ref = "pre-disturbance") #pre-disturbance as a reference 
veg_all$plots<- factor(veg_all$plot)

### null model ###
model_0<- lmer(alpha_rarefied ~ 1*sampling_stage + (1|plot), data = veg_all)
summary(model_0)

model_1 <-lmerTest::lmer(alpha_rarefied ~ FHD *sampling_stage+ (1|plot),
                data = veg_all)
summary(model_1)


model_2<- lmerTest::lmer(alpha_rarefied ~ canopy_total*sampling_stage+ (1|plot),
               data=veg_all)
summary(model_2)


model_3<- lmerTest::lmer(alpha_rarefied ~ leaf_litter_total * sampling_stage+ (1|plot),
                         data=veg_all)
summary(model_3)

model_set1 <- list(model_1, model_2, model_3, model_0)
mod.names1<- c('FHD', 'canopy', 'leaf', 'null')
aictab(cand.set = model_set1, modnames = mod.names1) # model_1 with FHD performs best

#Model selection based on AICc:

# K   AICc Delta_AICc AICcWt Cum.Wt  Res.LL
# FHD    6 275.05       0.00   0.94   0.94 -130.25
# leaf   6 280.50       5.46   0.06   1.00 -132.98
# canopy 6 298.87      23.83   0.00   1.00 -142.16
# null   3 307.56      32.52   0.00   1.00 -150.45

residuals_1 <- simulateResiduals(model_1) #null hypothesis is that there is no overdispersion
plot(residuals_1)
testDispersion(residuals_1) #dispersion = 1.0984, p-value = 0.672

# model_1 plot #
ggplot(veg_all, aes(x = FHD, y = alpha_rarefied, color = factor(sampling_stage))) +
  geom_smooth(method = "lm", se = TRUE, size = 1.5, alpha = 0.1) +  # Add best fit line
  labs(x = "FHD", y = "Alpha Rarefied",  color = "Sampling Stage") +
  theme_minimal()+
  theme(axis.text = element_text(size = 14),   # Increase axis label font size
        axis.title = element_text(size = 14),  # Increase axis title font size
        legend.title = element_text(size = 14))+
  ggtitle("Species richness (alpha rarefied) and FHD across sampling stages")

### using alpha richness (sub_canopy) as response variable ###
sub_canopy$sampling_stage<- relevel(sub_canopy$sampling_stage, ref = "pre-disturbance")

## null model
library(glmmTMB)
#model_00 <- glmmTMB(jackknife ~ 1 *sampling_stage+ (1|plot), data = sub_canopy, family = gaussian())
#summary(model_00)

model_00 <-lmerTest::lmer(jackknife~ sampling_stage+ (1|plot),
                          data = sub_canopy)
summary(model_00)

model_4<- lmerTest::lmer(jackknife ~ FHD*sampling_stage + (1|plot), 
               data=sub_canopy)
summary(model_4)

model_5<- lmerTest::lmer(jackknife ~ canopy_total* sampling_stage + (1|plot), 
               data=sub_canopy)
summary(model_5)

model_6<- lmerTest::lmer(jackknife ~ leaf_litter_total* sampling_stage + (1|plot), 
                         data=sub_canopy)
summary(model_6)


model_set2 <- list(model_4, model_5, model_6, model_00)
mod.names2<- c('FHD', 'canopy', 'leaf', 'null')
aictab(cand.set = model_set2, modnames = mod.names2)

# Model selection based on AICc:
#   
#   K   AICc Delta_AICc AICcWt Cum.Wt Res.LL
# null   4 -35.15       0.00   0.69   0.69  22.14
# FHD    6 -33.58       1.56   0.31   1.00  24.06
# leaf   6 -21.08      14.06   0.00   1.00  17.81
# canopy 6 -19.02      16.13   0.00   1.00  16.78

residuals_4 <- simulateResiduals(model_4) #null hypothesis is that there is no overdispersion
plot(residuals_4)
testDispersion(residuals_4) #dispersion = 0.92693, p-value = 0.776

# model_4 plot #
ggplot(sub_canopy, aes(x = FHD, y = jackknife, color = factor(sampling_stage))) +
  geom_smooth(method = "lm", se = TRUE, size = 1.5, alpha = 0.1) +  # Add best fit line
  labs(x = "FHD", y = "Alpha Rarefied",  color = "Sampling Stage") +
  theme_minimal()+
  theme(axis.text = element_text(size = 14),   # Increase axis label font size
        axis.title = element_text(size = 14),  # Increase axis title font size
        legend.title = element_text(size = 14))+
  ggtitle("Species richness (alpha sub canopy) and FHD across sampling stages")



###======= null model ######
model_0<- lmer(coverage_based.qD ~ 1 + (1|plot), data = veg_all)
summary(model_0)

model_1a <-lmer(coverage_based.qD ~ FHD + (1|plot),
                data = veg_all)
summary(model_1a)

delta_AIC_1a<- AIC(model_0)- AIC(model_1a)
print(delta_AIC_1a)

residuals_1a <- simulateResiduals(model_1a) #null hypothesis is that there is no overdispersion
plot(residuals_1a)
testDispersion(residuals_1a) #p-value greater than 0.05- cannot reject null hypothesis



# make a linear model first to check residual spread ##
model_1 <- lm(coverage_based.qD ~ FHD, data = veg_all)
plot(model_1, which=1)
qqnorm(residuals(model_1))
qqline(residuals(model_1))
shapiro.test(residuals(model_1)) # W = 0.91256, p-value = 0.004528 (<0.05) not normal
hist(residuals(model_1)) # negatively skewed

## data exploration ##
# ggplot(veg_all, aes(x = plot, y = alpha_rarefied, color=sampling_period)) +
#   geom_point() + 
#   geom_smooth(method = "lm", se = FALSE) +  # Add a smooth line
#   #labs(x = "Plot", y = "Predicted coverage_based.qD") +  # Add axis labels
#   theme_minimal()  # Apply a minimal theme



# make a linear model first to check residual spread ##
model_2 <- lm(coverage_based.qD ~ canopy_total, data = veg_all)
plot(model_2, which=1)
qqnorm(residuals(model_2))
qqline(residuals(model_2))
shapiro.test(residuals(model_2)) # not normal <0.05
hist(residuals(model_2)) # negatively skewed

model_2a <-lmer(coverage_based.qD ~ canopy_total + (1|plot),
                 data = veg_all)

summary(model_2a)

delta_AIC_2a<- AIC(model_0)- AIC(model_2a)
print(delta_AIC_2a)

residuals_2a <- simulateResiduals(model_2a)
plot(residuals_2a)
testDispersion(residuals_2a)


# make a linear model first to check residual spread ##
model_3 <- lm(coverage_based.qD ~ leaf_litter_total, data = veg_all)
plot(model_3, which=1)
qqnorm(residuals(model_3))
qqline(residuals(model_3))
shapiro.test(residuals(model_3)) # not normal <0.05
hist(residuals(model_3)) # negatively skewed

model_3a <-lmer(coverage_based.qD ~ leaf_litter_total + (1|plot),
                 data = veg_all)
summary(model_3a)

delta_AIC_3a<- AIC(model_0)- AIC(model_3a)
print(delta_AIC_3a)

residuals_3a <- simulateResiduals(model_3a)
plot(residuals_3a)
testDispersion(residuals_3a)


### using alpha sub-canopy richness as response variables ###


# > colnames(sub_canopy)
# [1] "PlotID"                 "plot"                   "sampling_period"        "season"                 "jackknife"             
# [6] "canopy_total"           "leaf_litter_total"      "avg_high_canopy_height" "avg_mid_canopy_height"  "avg_understory_height" 
# [11] "FHD"




# make a linear model first to check residual spread ##
model_11 <- lm(jackknife ~ FHD , data = sub_canopy)
plot(model_11, which=1)
qqnorm(residuals(model_11))
qqline(residuals(model_11))
shapiro.test(residuals(model_11)) # W = 0.91256, p-value = 0.004528 (<0.05) not normal
hist(residuals(model_11)) #highly skewed


model_1b<-lmer(jackknife ~ FHD + (1|sampling_period) + (1|plot),
                 data = sub_canopy)
summary(model_1b)

delta_AIC_1b<- AIC(model_00)- AIC(model_1b)
print(delta_AIC_1b)

residuals_1b <- simulateResiduals(model_1b)
plot(residuals_1b)
testDispersion(residuals_1b)


# make a linear model first to check residual spread ##
model_22 <- lm(jackknife ~ canopy_total, data = sub_canopy)
plot(model_22, which=1)
qqnorm(residuals(model_22))
qqline(residuals(model_22))
shapiro.test(residuals(model_22)) # normal >0.05
hist(residuals(model_22)) # normal distribution

model_2b <-lmer(jackknife ~ canopy_total + (1|sampling_period) + (1|plot),
                 data = sub_canopy)
summary(model_2b)

delta_AIC_2b<- AIC(model_00)- AIC(model_2b)
print(delta_AIC_2b)

residuals_2b <- simulateResiduals(model_2b)
plot(residuals_2b)
testDispersion(residuals_2b)

# make a linear model first to check residual spread ##
model_33 <- lm(jackknife ~ leaf_litter_total, data = sub_canopy)
plot(model_33, which=1)
qqnorm(residuals(model_33))
qqline(residuals(model_33))
shapiro.test(residuals(model_33)) 
hist(residuals(model_33)) # near normal

model_3b <-lmer(jackknife ~ leaf_litter_total+(1|sampling_period) + (1|plot),
                 data = sub_canopy)

summary(model_3b)

delta_AIC_3b<- AIC(model_00)- AIC(model_3b)
print(delta_AIC_3b)

residuals_3b <- simulateResiduals(model_3b)
plot(residuals_3b, nsim = 2000)
testDispersion(residuals_3b)