# import cleaned data ####
load(file="plants.RData")
native.ab <- plants$natives
alien.ab <- plants$aliens

dim(alien.ab)   # 50 73
dim(native.ab) # 50 178
# 251 total spp

# ofdall[match(rar20, names(ofdall))]

load(file="traits.RData")
names(traits)
trt.nat <- traits$natives
trt.ale <- traits$aliens

load(file="enviro.RData")

xydat<- env[,c(23,24,21,22)] 
names(xydat)# "long","lat", "north","east"; 

names(env)
#  [1] "Slope"                 "Altitude"              "Leaf.area.index"       "Live.basal.area"       "Phosphorus"            "Nitrate.Nitrogen"      "Ammonium.Nitrogen"    
#  [8] "pH"                    "Conductivity"          "Organic.Matter"        "Calcium"               "Magnesium"             "Potassium"             "Sodium"               
# [15] "Hydration"             "ECEC"                  "Carbon"                "Nitrogen"              "Carbon.Nitrogen.ratio" "pk"                    "east"                 
# [22] "north"                 "long"                  "lat"                  "acid"  

# pre-process ####
# 01. SR, freq omni ####
# frequency in plots
natab <- apply(native.ab, 1, function(x) sum(x/25)/sum(x>0))
aleab <- apply(alien.ab, 1, function(x) sum(x/25)/sum(x>0))

summary(apply(native.ab,1,sum))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 151.0   359.0   443.5   428.8   510.0   675.0 
summary(apply(alien.ab,1,sum))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 40.0   125.2   222.5   219.5   312.0   432.0 

# occupancy freq dist
a.ofd <- apply(alien.ab>0,2,sum)
n.ofd <- apply(native.ab>0,2,sum)

ave.ofda <- apply(alien.ab/25, 2, mean)
ave.odfn <- apply(native.ab/25,2,mean)

# SR srich
# rare natives subsetting 
n.rar.sr <- apply(native.ab[,which(n.ofd<10)]>0,1,sum)
a.rar.sr <- apply(alien.ab[,which(a.ofd<10)]>0,1,sum)
n.rar.ab <- apply(native.ab[,which(n.ofd<10)]/25,1,mean)
a.rar.ab <- apply(alien.ab[,which(a.ofd<10)]/25,1,mean)


ale.sr <- apply(alien.ab>0,1,sum)
nat.sr <- apply(native.ab > 0, 1, sum)

mean(ale.sr) # 19.42
sd(ale.sr) # 8.246434

quantile(ale.sr, c(0.025,0.975)) # 6.450 31.775
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.00   11.25   19.50   19.42   27.00   34.00 

mean(nat.sr) # 40.26
sd(nat.sr) # 8.630298
summary(nat.sr)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 22.00   33.25   42.00   40.26   46.75   57.00 

quantile(nat.sr, c(0.025,0.975)) # 24.00 54.55 

t.test(nat.sr, ale.sr, paired=TRUE)
# t = 10.213, df = 49, p-value = 9.944e-14
# 20.84 

boxplot(nat.sr, ale.sr)

# simple linear regression native SR ~ alien SR 
lmsr <- lm(nat.sr ~ ale.sr)
summary(lmsr)

# post hoc tests ####
# role of soil pH
scatter.smooth(aleab ~ env$pH)
scatter.smooth(natab ~ env$pH)

# rare native Species richness
boxplot(n.rar.sr ~ env$acid)
kruskal.test(n.rar.sr ~ env$acid)
# Kruskal-Wallis chi-squared = 15.462, df = 1, p-value = 8.417e-05

tapply(n.rar.sr, env$acid, median)
#  <5.5    > 5.5 
#    8     5

# and frequency
boxplot(n.rar.ab ~ env$acid)
tapply(n.rar.ab, env$acid, mean)
kruskal.test(n.rar.ab ~ env$acid)
# Kruskal-Wallis chi-squared = 7.6331, df = 1, p-value = 0.005731


# Aliens ~ pH
boxplot(a.rar.sr ~ env$acid)
boxplot(a.rar.ab ~ env$acid)
# rare aliens opposite to rare natives...

# role of the C:N ratio 
scatter.smooth(ale.sr ~ env$CNratio)
abline(v = 21, lty=2)

# subset into high and low CN ratio
hic <- which(env$CNratio>21)
lop <- which(env$pH<5)
alien.lop <- alien.ab[lop,]
alien.lop <- alien.lop[,which(apply(alien.lop,2,sum)>0)]; dim(alien.lop)

alien.hic <- alien.ab[hic,]
alien.hic <- alien.hic[,which(apply(alien.hic,2,sum)>0)]; dim(alien.hic)
ale.tol <- names(alien.hic)[names(alien.hic) %in% names(a.ofd)[a.ofd>25]]

# names of the tolerant aliens
baddies <- names(alien.lop)[which(names(alien.lop) %in% ale.tol)]
# [1] "Aira.elegantissima"     "Anthoxanthum.aristatum" "Briza.spp."             "Centaurium.spp."        "Cirsium.vulgare"        "Galium.murale"         
# [7] "Hypericum.perforatum"   "Hypochaeris.glabra"     "Hypochaeris.radicata"   "Lysimachia.arvensis"    "Sonchus.asper"          "Sonchus.oleraceus"     
# [13] "Trifolium.spp."         "Vulpia.spp."           

# subset tolerant aliens, compare their frequency at low and high C:N ratio
alebad <- alien.ab[,which(names(alien.ab) %in% baddies)]
names(alebad)
bad.ab <- apply(alebad/25,1,mean)
boxplot(bad.ab~acid, data=env)

tapply(bad.ab, env$acid, mean)
# 0.2934694 0.5402956 
kruskal.test(bad.ab ~ env$acid)

# Extent of occurrence (EOO) ###
# calc polyhull 
library(sp)
xydf <- data.frame(y=env$north, x =env$east)
x1 = env$east
y1 = env$north
hpts <- chull(x=x1, y=y1 )
hpts <- c(hpts, hpts[1])
xy <- cbind(env$east,env$north)
chull.coords <- xy[hpts,]
chull.poly <- Polygon(chull.coords, hole=F)
chull.area <- chull.poly@area
tot.ext <- chull.area/(1000*1000)

# is the EOO of common native vs alien spp different?
a50 <- which(a.ofd>20) # 20 as a threhsol
a75 <- which(a.ofd> 37)
aspp50 <- alien.ab[,a50]/37
asr5 <- apply(aspp50>0,1,sum)/length(a50)
asri <- which(asr5>0.75)

xy50a <- env[asri,c(21,22)]

# EOO calcs
hpta5 <- chull(xy50a)
hpta5 <- c(hpta5, hpta5[1])
xya5 <- cbind(xy50a[,1], xy50a[,2])
chull.corda5 <- xya5[hpta5,]
chull.polya5 <- Polygon(chull.corda5, hole=F)
ale5EOO <- chull.polya5@area
ale5EOO/(1000*1000)


hptn5 <- chull(xy50n)
hptn5 <- c(hptn5, hptn5[1])
xyn5 <- cbind(xy50n[,1], xy50n[,2])
chull.cordn5 <- xyn5[hptn5,]
chull.polyn5 <- Polygon(chull.cordn5, hole=F)
nat5EOO <- chull.polyn5@area
nat5EOO/(1000*1000)

plot(env$east,env$north)
points(xy50a,col=2)
lines(chull.corda5,col=2)
lines(chull.cordn5,col=4)

eoo.fn <- function(dat, xy, thresh = 5){
  eoo <- numeric()
  for(i in 1:ncol(dat)){
    idx <- which(dat[,i]>0)
    if(length(idx)< thresh) next
    xyi <- xy[idx,]
    pts <- chull(xyi)
    pts <- c(pts, pts[1])
    xya <- cbind(xyi[,1], xyi[,2])
    cords <- xya[pts,]
    polyh <- Polygon(cords, hole=F)
    eoo[i] <- polyh@area/(1000*1000)
  }
  return(eoo)
}

a.eoo <- eoo.fn(dat=alien.ab, xy=xydf, thresh=3)
n.eoo <- eoo.fn(dat=native.ab, xy=xydf, thresh=3)

boxplot(a.eoo[which(a.ofd>25)], n.eoo[which(n.ofd>25)]  )

gtr50 <-matrix(c(table(a.ofd>25),table(n.ofd>25)),nrow=2)
fisher.test(gtr50)
boxplot(a.eoo ~ a.ofd>25)
boxplot(n.eoo ~ n.ofd>25)
eoo.lst <- list(alien=a.eoo,native=n.eoo)
boxplot(eoo.lst, main="extent of occurrence (polyhull area)", ylab="EOO (km2)")

plot(a.eoo ~ a.ofd)
plot(n.eoo ~ n.ofd)
kruskal.test(eoo.lst)

# eoo Fig 2 ####
eoovec <- c(a.eoo, n.eoo)
spp <- c(rep("Alien", length(a.eoo)), rep("Native",length(n.eoo)))
pldat <- data.frame(occ = eoovec, comp = as.factor(spp))

# makes sure you specify the x-predictor as a factor in the dataframe...
library(ggplot2)
bp1 <- ggplot(aes(y = occ,  x = spp), data = pldat) + 
  theme_minimal()+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = c(0.87, 0.2),
        legend.title=element_blank(),
        panel.grid = element_blank())+
  geom_boxplot()+
  #labs(title ="Macroinvertebrate FFG responses")+
  xlab("")+
  ylab(expression('Extent of occurrence ( km' ^ '2'~')'))+
  theme(legend.position="none") # suppress the legend
print(bp1)

# just the tolerant ones
boxplot(a.eoo[-which(names(a.ofd)%in% baddies)], a.eoo[which(names(a.ofd)%in% baddies)])
length(a.eoo[-which(names(a.ofd)%in% baddies)])  #59
length(a.eoo[which(names(a.ofd)%in% baddies)]) # 14

eoovec2 <- c(a.eoo[-which(names(a.ofd)%in% baddies)], a.eoo[which(names(a.ofd)%in% baddies)])
spp2 <- c(rep("Other", 59), rep("Tolerant",14))
pldat2 <- data.frame(occ = eoovec2, comp = as.factor(spp2))
pldat2$comp = factor(pldat2$comp, levels = c("Tolerant", "Other"))

bp2 <- ggplot(aes(y = occ,  x = comp), data = pldat2) + 
  theme_minimal()+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = c(0.87, 0.2),
        legend.title=element_blank(),
        panel.grid = element_blank())+
  geom_boxplot()+
  #labs(title ="Macroinvertebrate FFG responses")+
  xlab("")+
  ylab(expression('Extent of occurrence ( km' ^ '2'~')'))+
  theme(legend.position="none"); print(bp2)


# pH and/or C/N ratio fig 
pldat3 <- data.frame(sr = c(ale.sr,nat.sr), cnrat = rep(env$CNratio,2), comp = c(rep("Alien",50),rep("Native",50)))
# pldat3 <- data.frame(sr = c(ale.sr,nat.sr), cnrat = rep(env$pH,2), comp = c(rep("Alien",50),rep("Native",50)))

labs <- data.frame(x= c(15,15), y = c(40,70), lab=c("A. Alien", "B. Native"), comp= c("Alien","Native"))
bp3 <- ggplot(aes(y = sr,  x = cnrat), data = pldat3) + 
  theme_minimal()+
  theme_bw()+
  facet_wrap(~comp,ncol=1, scales="free_y")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = c(0.87, 0.2),
        legend.title=element_blank(),
        panel.grid = element_blank())+
  geom_smooth(method="lm", col=1)+
  geom_point(pch=21,fill="white",size = 3)+
  #geom_smooth(method="loess",span=0.7,col=1)+
  #xlim(c(12,27))+
  xlab("Carbon-Nitrogen ratio")+
  ylab(expression('Species richness'))+
  geom_text(aes(x=x,y=y,label=lab),data=labs)+
  theme(legend.position="none"); print(bp3)

#***********************************************************************************************
# analysis code #### 
# glms ####
# standardise data 
envstd <- env
envstd[,-c(20,25)] <- vegan::decostand(env[,-c(20,25)], method="standardize")
xyd <- dist(env[,c(21,22)]) # euc distance between sites 

# richness GLMs ####
p01 <- MASS::glm.nb(nat.sr ~pH + LiveBasalArea + ECEC + OrgMatter, data=envstd);summary(p01)
p02 <- MASS::glm.nb(ale.sr ~ pH + LiveBasalArea + ECEC + OrgMatter, data=envstd); summary(p02)

# spatial autocor
mantel(xdis= xyd, dist(resid(p01))) # Mantel statistic r: 0.01277; Significance: 0.378 
mantel(xdis= xyd, dist(resid(p02))) # Mantel statistic r: 0.08296; Significance: 0.132 

cor(predict(p01), nat.sr)^2 # 0.2917144
cor(predict(p02), ale.sr)^2 #  0.5331785 

# output summary tables 
sjPlot::tab_model(p01, show.se= TRUE,show.stat=TRUE, digits = 3, transform=NULL, string.est = "Estimate", show.ci = FALSE)#, file="natsrNB.doc")

# manyglm ####
# caution- can take many hours  to run 
library(mvabund)
alien_abu <- mvabund(alien.ab)
native_abu <- mvabund(native.ab)

alefr <- alien.ab/25

meanvar.plot(alefr)
names(envstd)
alien.abm <- as.matrix(alien.ab)
native.abm <- as.matrix(native.ab)

# mg.aleab <- manyglm(alien_abu ~ OrgMatter + LiveBasalArea + ECEC + pH, data=env, family="negative binomial")

# mg.natab <- manyglm(native_abu ~OrgMatter + LiveBasalArea + ECEC + pH, data=env, family="negative binomial")

# ale.nb <- summary(mg.aleab, test="LR", nBoot=9999)
# nat.nb <- summary(mg.natab, test="LR", nBoot=9999)

# extract median coefficient values 
apply(mg.aleab$coef,1,median)
# (Intercept)  Organic.Matter Live.basal.area            ECEC              pH 
# -2.23294373     -0.05800055     -0.30608536      0.27312159      0.42711322 

apply(mg.natab$coef,1,median)
# Intercept)  Organic.Matter Live.basal.area            ECEC              pH 
# 7.96723631     -0.03471003      0.02908231      0.01074266     -1.11053153 

# zetadiv ####
library(zetadiv)
native.pa <- as.data.frame(ifelse(native.ab>0,1,0))
alien.pa <- as.data.frame(ifelse(alien.ab>0,1,0))

# zeta decline ####  
set.seed(11)
# raw zeta random combinations
zdnatAllRaw <- Zeta.decline.mc(data.frame(native.pa), orders= 1:50, normalize = FALSE, sam=10000, NON= FALSE, xy = xydf)

# Raw zeta near-neighbour
zdnatNNRaw <- Zeta.decline.mc(data.frame(native.pa), orders= 1:50, normalize = FALSE, sam=10000, NON= TRUE, xy = xydf)

# Simpson random 
zdnatAllSim <- Zeta.decline.mc(data.frame(native.pa), orders= 1:50, normalize = "Simpson", sam=10000, NON= FALSE, xy = xydf)

# Simpson near-neighbour
zdnatNNSim <- Zeta.decline.mc(data.frame(native.pa), orders= 1:50, normalize = "Simpson", sam=10000, NON= TRUE, xy = xydf)



# msgdm ####
# alien eg
ale.g2 <- Zeta.msgdm(data.spec = alien.pa, data.env = env[,c(16,4,10,8)], xy = xydf, data.spec.pred = NULL,
                     order = 2,  sam = 10000,  reg.type = "ispline",  family = stats::gaussian(),  method.glm = "glm.fit.cons",
                     cons = -1,   cons.inter = 1,   confint.level = 0.95,   bs = "mpd",   kn = -1,  order.ispline = 2,
                     kn.ispline = 1,  distance.type = "Euclidean",  dist.custom = NULL,  rescale = FALSE,
                     rescale.pred = TRUE,  method = "mean",  normalize = FALSE,  silent = FALSE,  empty.row = 0,
                     control = list(),  glm.init = FALSE)

# repeat changing for order 2-10 and for alien.pa to native.pa 

# plotting (e.g., Fig 5), repeat for higher orders, aliens and natives
pag2 <- Plot.ispline(msgdm=ale.g2, data.env = env[,c(16,4,10,8)], distance=TRUE)#, biotic=TRUE)

# variation partitioning #####
vpa2 <- Zeta.varpart(ale.g2, reg.type="ispline")#[4:7,]

