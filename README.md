# Diversity-and-Multispecies-Invasion
Data and example code from *Rare, common, alien and native species follow different rules in an understorey plant community*.

# Summary
This repo holds the data and R code to run the analysis from the paper *Rare, common, alien and native species follow different rules in an understorey plant community* by S.Reeve, D.Deane, C. McGrannachan, G.Horner, C.Hui and M.McGeoch, published in Ecology and Evolution in 2022. In this study, we look at the understorey plant community in a set of national parks nestled in an agricultural landscape matrix in the central north of Victoria, Australia.   

![figA1_1.png](./figA1_1.png)   

## Abstract
Biological invasions are a leading threat to biodiversity globally. Increasingly, ecosystems experience multiple introductions, which can have significant effects on patterns of diversity. The way these communities assemble will depend partly on whether rare and common alien species respond to environmental predictors in the same manner as rare and common native species, but this is not well understood. To examine this question across four national parks in south-eastern Australia, we sampled the understorey plant community of eucalypt-dominated dry forest subject to multiple plant introductions. The drivers of diversity and turnover in alien and native species of contrasting frequency of occurrence (low, intermediate, and high) were each tested individually. We found alien species diversity and turnover were both strongly associated with abiotic conditions (e.g., soil pH), while distance had little influence because of the greater extent of occurrence and more homogeneous composition of common aliens. In contrast, native species diversity was not associated with abiotic conditions and their turnover was as strongly influenced by distance as by abiotic conditions. In both alien and native species, however, the most important predictors of turnover changed with frequency of occurrence. Although local coexistence appears to be facilitated by life history trade-offs, species richness of aliens and natives was negatively correlated and native species might face greater competition in areas with more neutral soils (e.g, pH > ~5.5) where alien richness and relative frequency were both highest. We conclude that diversity and turnover in the generally more widespread alien species is mainly driven by species sorting along an environmental gradient associated with pH and nutrient availability, whereas turnover of native species is driven by more neutral processes associated with dispersal limitation. We show alien and native plant species respond to different environmental factors, as do rare and common species within each component.

## Data  
There are three R data objects 'plants.RData', 'traits.RData' and 'enviro.RData'. These contain respectively: 
1. 'plants' a list of two elements ('native' and 'alien'), which give the plant frequencies (integer count between 1 and 25 for the number of adjacent sub quadrats in which it was present) for each of 50 plots for 178 native species and 73 alien species observed at the sites. 
2. 'traits' a list of two elements giving the traits data used for native and alien species, which were life history ('life') and growth form ('form'). 
3. the environmental predictors (whether used in the analysis or not) as shown in Fig. A1.6 in the text.

```
load('plants.RData')
dim(plants$native)
dim(plants$alien)
```

```
load('traits.RData')
head(traits$native)
head(traits$alien)
```

```
load('enviro.RData')
names(env)
```



