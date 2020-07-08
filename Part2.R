library(dplyr)
library(tidyverse)
library(ggplot2)
library(MCMCpack)
library(mclust)
library(factoextra)
wine_data <- read.csv('C:/Users/Rohit/Desktop/winemag-data-130k-v2.csv')
wine_data_USA <- wine_data%>%filter(
  country == 'US'
)%>%select(points,price)
wine_data_USA <- na.omit(wine_data_USA)
wine_data_USA <- wine_data_USA%>%filter(
  price < 800
)
wine_data_USA$variety <- factor(wine_data_USA$variety)
wine_data_USA <- scale(wine_data_USA)
fit <- Mclust(wine_data_USA)
summary(fit)
plot(fit, what = "classification", log = 'y')
plot(fit, what = "uncertainty", log = 'y')
plot(fit, what = "BIC")
fit$BIC
fit2 <- Mclust(wine_data_USA, G = 8, modelNames = "VVV")
plot(fit2, what = "classification")
plot(fit2, what = "uncertainty")
fviz_mclust(fit, "classification", geom = "point", 
            pointsize = 1.5, palette = "jco")
fviz_mclust(fit, "BIC", palette = "jco")
fviz_mclust(fit, "uncertainty", palette = "jco")
cluster2 <- wine_data_USA %>%
  scale() %>%
  as.data.frame() %>%
  mutate(cluster = fit$classification) %>%
  filter(cluster == 8) %>%
  select(-cluster)
wine_data_USA$cluster <- fit$classification
cluster2<- wine_data_USA%>%
  select(cluster)%>%
  group_by(cluster)%>%
  summarise(counts = n())
wine_data_USA$cluster <- factor(wine_data_USA$cluster)
ggplot(wine_data_USA, aes(cluster))+
  geom_bar()

probabilities <- fit$z 
colnames(probabilities) <- paste0('C', 1:8)

probabilities <- probabilities %>%
  as.data.frame() %>%
  mutate(id = row_number()) %>%
  tidyr::gather(cluster, probability, -id)

ggplot(probabilities, aes(probability)) +
  geom_histogram() +
  facet_wrap(~ cluster, nrow = 2)

uncertainty <- data.frame(
  id = 1:nrow(wine_data_USA),
  cluster = fit$classification,
  uncertainty = fit$uncertainty
)

uncertainty %>%
  group_by(cluster) %>%
  filter(uncertainty > 0.45) %>%
  ggplot(aes(uncertainty, reorder(id, uncertainty))) +
  geom_point() +
  facet_wrap(~ cluster, scales = 'free_y', nrow = 3)+
  theme(axis.text.y = element_blank())