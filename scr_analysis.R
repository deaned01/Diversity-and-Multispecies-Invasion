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
a50 <- which(a.ofd>20)
a75 <- which(a.ofd> 37)
aspp50 <- alien.ab[,a50]/37
asr5 <- apply(aspp50>0,1,sum)/length(a50)
asri <- which(asr5>0.75)
a50i <- which(apply(aspp50,1,mean)>0.25)
aspp75 <- alien.ab[,a75]/25
a75i <- which(apply(aspp75,1,mean)>0.5)

n50 <- which(n.ofd>20)
n75 <- which(n.ofd>37)
nspp50 <- native.ab[,n50]/37
n50i <- which(apply(nspp50,1,mean)>0.25)
nsr5 <- apply(nspp50>0,1,sum)/length(n50)
nsri <- which(nsr5>0.75)
nspp75 <- native.ab[,n75]/25
n75i <- which(apply(nspp75,1,mean)>0.5)

xy50a <- env[asri,c(21,22)]
xy75a <- env[a75i,c(21,22)]
xy50n <- env[nsri,c(21,22)]
xy75n <- env[n75i,c(21,22)]

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


# analysis code 
# 05. glms ####
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
sjPlot::tab_model(p02, show.se= TRUE, show.stat=TRUE, digits = 3, transform=NULL, string.est = "Estimate", show.ci = FALSE, file="alesrNB.doc")

# MV frequency = manyglm ####
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

boxplot(t(mg.natab$coef[2:5,]),ylim=c(-10,10)); abline(0,0,lty=2)
boxplot(t(mg.aleab$coef[2:5,]),ylim=c(-10,10)); abline(0,0,lty=2)

# aa.nb <- anova(mg.aleab)
# an.nb <- anova(mg.natab)

# mg.aleabp <- manyglm(alien.abm ~ pH + LiveBasalArea + ECEC + Organic.Matter, data=envstd, family="poisson")
# mg.natabp <- manyglm(native.abm ~ pH + LiveBasalArea + ECEC + Organic.Matter, data=envstd, family="poisson")
# 
# ale.pois<- summary(mg.aleabp)
# nat.pois <- summary(mg.natabp)

# compare with frequency instead of count 
ale.cofmg <- t(mg.aleab$coeff)
boxplot(ale.cofmg[,-1], ylim=c(10,-10),lab="Coefficient value", main="Alien");abline(0,0,lty=2)

nat.cofmg <- t(mg.natab$coeff)
boxplot(nat.cofmg[,-1], ylim=c(10,-10),lab="Coefficient value",main="Native");abline(0,0,lty=2)

mgsum.ale.ns <- summary(mg.aleab)
mgsum.nat.ns <- summary(mg.natab)

plot(mg.aleab$coef[2,],ale.ab)
anova(mg.aleab)
hist(apply(native.ab,2,mean))

umg.aleab <- manyglm(alien.abm ~ pH + Nitrogen + Slope + Sodium, data=envstd, family="negative binomial")
umg.natab <- manyglm(native.abm ~ pH + Nitrogen + Slope + Sodium, data=envstd, family="negative binomial")

umgsum.ale <- summary(umg.aleab)
umgsum.nat <- summary(umg.natab)



head(alien.abm)
plot(mg.aleab)
summary(mg.aleab)
anova(mg.aleab)

str(mg.aleab, max.level=1)

# zetadiv ####
library(zetadiv)
native.pa <- as.data.frame(ifelse(native.ab>0,1,0))
alien.pa <- as.data.frame(ifelse(alien.ab>0,1,0))

# zeta decline 
set.seed(11)
# raw zeta random combinations
zdnatAllRaw <- Zeta.decline.mc(data.frame(native.pa), orders= 1:50, normalize = FALSE, sam=10000, NON= FALSE, xy = xydf)
zdaleAllRaw <- Zeta.decline.mc(data.frame(alien.pa), orders= 1:50, normalize = FALSE, sam=10000, NON= FALSE, xy = xydf)

# Raw zeta near-neighbour
zdnatNNRaw <- Zeta.decline.mc(data.frame(native.pa), orders= 1:50, normalize = FALSE, sam=10000, NON= TRUE, xy = xydf)
zdaleNNRaw <- Zeta.decline.mc(data.frame(alien.pa), orders= 1:50, normalize = FALSE, sam=10000, NON= TRUE, xy = xydf)

# Simpson random 
zdnatAllSim <- Zeta.decline.mc(data.frame(native.pa), orders= 1:50, normalize = "Simpson", sam=10000, NON= FALSE, xy = xydf)
zdaleAllSim <- Zeta.decline.mc(data.frame(alien.pa), orders= 1:50, normalize = "Simpson", sam=10000, NON= FALSE, xy = xydf)

# Simpson near-neighbour
zdnatNNSim <- Zeta.decline.mc(data.frame(native.pa), orders= 1:50, normalize = "Simpson", sam=10000, NON= TRUE, xy = xydf)
zdaleNNSim <- Zeta.decline.mc(data.frame(alien.pa), orders= 1:50, normalize = "Simpson", sam=10000, NON= TRUE, xy = xydf)

zdnat <- list(nAllRaw = zdnatAllRaw, nNNRaw = zdnatNNRaw, aAllRaw = zdaleAllRaw, aNNRaw = zdaleNNRaw,
              nAllSim = zdnatAllSim, nNNSim = zdnatNNSim, aAllSim = zdaleAllSim, aNNSim = zdaleNNSim)

names(zdnat)zdnatAllSi
save(zdnat, file="zetaDeclineRawFits.RData")
load(file="zetaDeclineRawFits.RData")


nA <- zdnatAllRaw # zdnat$nAllRaw
aA <- zdaleAllRaw # zdnat$aAllRaw
nX <- zdnatNNRaw # zdnat$nNNRaw
aX <- zdaleNNRaw # zdnat$aNNRaw

nAs <- zdnatAllSim # zdnat$nAllRaw
aAs <- zdaleAllSim # zdnat$aAllRaw
nXs <- zdnatNNSim # zdnat$nNNRaw
aXs <- zdaleNNSim # zdnat$aNNRaw


zval <- list(nA$zeta.val, nX$zeta.val, aA$zeta.val,  aX$zeta.val,
          nAs$zeta.val, aAs$zeta.val, nXs$zeta.val, aXs$zeta.val)
zsd <- list(nA$zeta.val.sd, nX$zeta.val.sd, aA$zeta.val.sd, aX$zeta.val.sd,
         nAs$zeta.val.sd, aAs$zeta.val.sd, nXs$zeta.val.sd, aXs$zeta.val.sd)

zeta <- unlist(lapply(zval, function(x) x))
zlo <- lapply(1:length(zval), function(x, y, i) x[[i]] - y[[i]],x=zval,  y = zsd)
zhi <- lapply(1:length(zval), function(x, y, i) x[[i]] + y[[i]],x=zval,  y = zsd)
zlov <- unlist(zlo)
zhiv <- unlist(zhi)

ords <- rep(1:50, 8)
crv <- rep(letters[1:8],each=50)

# ggplot zeta decline 
pan <- c(rep("a", 100), rep("b",100), rep("c",100),rep("d",100))
tfm <- c(rep("raw",200),rep("simpson",200))
sch <- c(rep("All",50), rep("NN",50),rep("All",50), rep("NN",50),rep("All",100), rep("NN",100))

alenat <- c(rep("Native",100),rep("Alien",100),rep("Native",50),rep("Alien",50),rep("Native",50),rep("Alien",50))
ptyp <- c(rep("A",100),rep("C",100),rep("B",100),rep("D",100))
#crv <- c(rep("a1",50),rep("a2",50),rep("a3",50),rep("a4",50))
zet <- data.frame(zeta = zeta, zhi = zhiv, zlo = zlov, ords = ords, pan=pan, tfm = tfm, ptyp=ptyp, 
                  orig = alenat, sch=sch,crv=crv)

zet2 <- zet#[1:200,1:3] <- .ln <- 
zet2[1:200,1:3] <- log2(zet[1:200,1:3])

# correct CI > 1 for simpson
names(zet)
zet[which(zet$tfm=="simpson"& zet$zhi>1),1:2] <- 1
zet$int <- interaction(zet$sch,zet$orig)
levels(zet$int)
par(mfrow=c(1,1))
library(scales)

zet$int <- as.factor(zet$int)
zet$int <- factor(zet$int, levels = c("All.Alien",  "NN.Alien",   "All.Native", "NN.Native"),
                  labels = c("Alien, random","Alien near-neighbour", "Native, random","Native near-neighbour"))

labs <- data.frame(y=c(6,1.02,5,1.02),x=c(1,1,1,1),lab=c("(a)","(c)","(b)","(d)"), 
                   int = c("Alien, random","Alien near-neighbour", "Native, random","Native near-neighbour"), ptyp=c("A","B","C","D"),
                   crv = rep("a",8), tfm = c("raw","simpson","raw","simpson"))
zetpl.dat <- zet2
zplt  <- ggplot(data=zetpl.dat, aes(x=ords, y = zeta, group= int,col=int))+
  geom_point(aes(x = ords, y = zeta, group=crv, pch = int),size=2)+
  geom_line(aes(x = ords, y = zlo, group=crv), linetype=2)+
  geom_line(aes(x = ords, y = zhi, group=crv), linetype=2)+
  theme_minimal()+
  theme(strip.background = element_blank(),
        panel.background = element_rect(),
        strip.text = element_blank(),
        #legend.position="none",
        legend.position = c(0.85, 0.88),
        legend.title=element_blank(),
        panel.grid = element_blank())+
  xlab("Order of zeta")+
  ylab("Zeta diversity")+
  scale_color_manual(values=c("#0072B2", "#009E73","#E69F00", "#999999"))+
  scale_shape_manual(values = c(1,2,5,6))+#,labels=c("Modelled", "Observed"))+
  # scale_linetype_discrete(c("dashed", "dotted"))+#,labels=c("Modelled", "Observed"))+
#  scale_y_continuous(trans = 'log2', expand = c(0,0) , breaks = trans_breaks("log2", function(x) 2^x))+
  geom_text(aes(x=x,y=y,label = lab),fontface="bold",data=labs,col=1)+
  facet_wrap(ptyp~tfm, ncol=2, scales="free_y")

ggsave(zplt, file="fig4.tiff", width = 150, height = 150, units = "mm", dpi = 800)
ggsave(zplt, file="fig4.png", width = 150, height = 150, units = "mm", dpi = 100)
ggsave(zplt, file="fig4.eps", width = 150, height = 150, units = "mm", dpi = 800)

save(zetpl.dat, file="zetaDeclinePlotDat.RData")

envgdm <- envstd[,c(4,8, 10,16)]#dataframe()


# msgdm ####
# custom fn to extract quantiles (hacked from zetadiv::Plot.ispline() )

pts.hack <- function(isplines, nq = 11){
  # isplines is the output of zetadiv::Return.ispline() 
  # nq= number of quantiles for the ispline (num.quantiles arg from orig fn)
  
  env.resc <- isplines$env
  Isplines.pred <- isplines$Ispline
  data.env.num <- isplines$env
  distance <- isplines$distance
  my.order <- isplines$my.order
  num.quantiles <- nq
  
  outx <- outy <- matrix(nrow=num.quantiles, ncol = ncol(data.env.num))
  
  ind.points<- numeric()
  
  for(i in 1:ncol(data.env.num)){
    for (j in 1:num.quantiles) {
      ind.points[j] <- which.min(abs(stats::quantile(env.resc[, 
                                                              i], seq(0, 1, 1/(num.quantiles - 1)))[j] - 
                                       sort(env.resc[, i])))
    }
    
    outx[,i] <- sort(env.resc[, i])[ind.points]
    outy[,i] <- Isplines.pred[order(env.resc[,i]), i][ind.points]
  }
  out <- list(x = outx, y = outy)
  return(out)
}

plot(pp$x[,1], pp$y[,1],ylim=c(0,8))
lines(pp$x[,1],pp$y[,1])
names(env)

# alien ####
ale.g2 <- Zeta.msgdm(data.spec = alien.pa, data.env = env[,c(16,4,10,8)], xy = xydf, data.spec.pred = NULL,
                     order = 2,  sam = 10000,  reg.type = "ispline",  family = stats::gaussian(),  method.glm = "glm.fit.cons",
                     cons = -1,   cons.inter = 1,   confint.level = 0.95,   bs = "mpd",   kn = -1,  order.ispline = 2,
                     kn.ispline = 1,  distance.type = "Euclidean",  dist.custom = NULL,  rescale = FALSE,
                     rescale.pred = TRUE,  method = "mean",  normalize = FALSE,  silent = FALSE,  empty.row = 0,
                     control = list(),  glm.init = FALSE)


ale.g4 <- Zeta.msgdm(data.spec = alien.pa, data.env = env[,c(16,4,10,8)], xy = xydf, data.spec.pred = NULL,
                     order = 4,  sam = 10000,  reg.type = "ispline",  family = stats::gaussian(),  method.glm = "glm.fit.cons",
                     cons = -1,   cons.inter = 1,   confint.level = 0.95,   bs = "mpd",   kn = -1,  order.ispline = 2,
                     kn.ispline = 1,  distance.type = "Euclidean",  dist.custom = NULL,  rescale = FALSE,
                     rescale.pred = TRUE,  method = "mean",  normalize = FALSE,  silent = FALSE,  empty.row = 0,
                     control = list(),  glm.init = FALSE)


ale.g5 <- Zeta.msgdm(data.spec = alien.pa, data.env = env[,c(16,4,10,8)], xy = xydf, data.spec.pred = NULL,
                     order = 5,  sam = 10000,  reg.type = "ispline",  family = stats::gaussian(),  method.glm = "glm.fit.cons",
                     cons = -1,   cons.inter = 1,   confint.level = 0.95,   bs = "mpd",   kn = -1,  order.ispline = 2,
                     kn.ispline = 1,  distance.type = "Euclidean",  dist.custom = NULL,  rescale = FALSE,
                     rescale.pred = TRUE,  method = "mean",  normalize = FALSE,  silent = FALSE,  empty.row = 0,
                     control = list(),  glm.init = FALSE)

ale.g10 <- Zeta.msgdm(data.spec = alien.pa, data.env = env[,c(16,4,10,8)], xy = xydf, data.spec.pred = NULL,
                     order = 10,  sam = 10000,  reg.type = "ispline",  family = stats::gaussian(),  method.glm = "glm.fit.cons",
                     cons = -1,   cons.inter = 1,   confint.level = 0.95,   bs = "mpd",   kn = -1,  order.ispline = 2,
                     kn.ispline = 1,  distance.type = "Euclidean",  dist.custom = NULL,  rescale = FALSE,
                     rescale.pred = TRUE,  method = "mean",  normalize = FALSE,  silent = FALSE,  empty.row = 0,
                     control = list(),  glm.init = FALSE)


# native ####
nat.g2 <- Zeta.msgdm(data.spec = native.pa, data.env = env[,c(16,4,10,8)], xy = xydf, data.spec.pred = NULL,
                     order = 2,  sam = 10000,  reg.type = "ispline",  family = stats::gaussian(),  method.glm = "glm.fit.cons",
                     cons = -1,   cons.inter = 1,   confint.level = 0.95,   bs = "mpd",   kn = -1,  order.ispline = 2,
                     kn.ispline = 1,  distance.type = "Euclidean",  dist.custom = NULL,  rescale = FALSE,
                     rescale.pred = TRUE,  method = "mean",  normalize = FALSE,  silent = FALSE,  empty.row = 0,
                     control = list(),  glm.init = FALSE)


nat.g4 <- Zeta.msgdm(data.spec = native.pa, data.env = env[,c(16,4,10,8)], xy =  xydf, data.spec.pred = NULL,
                     order = 4,  sam = 10000,  reg.type = "ispline",  family = stats::gaussian(),  method.glm = "glm.fit.cons",
                     cons = -1,   cons.inter = 1,   confint.level = 0.95,   bs = "mpd",   kn = -1,  order.ispline = 2,
                     kn.ispline = 1,  distance.type = "Euclidean",  dist.custom = NULL,  rescale = FALSE,
                     rescale.pred = TRUE,  method = "mean",  normalize = FALSE,  silent = FALSE,  empty.row = 0,
                     control = list(),  glm.init = FALSE)


nat.g5 <- Zeta.msgdm(data.spec = native.pa, data.env = env[,c(16,4,10,8)], xy =  xydf, data.spec.pred = NULL, #alien.pa,
                     order = 5,  sam = 10000,  reg.type = "ispline",  family = stats::gaussian(),  method.glm = "glm.fit.cons",
                     cons = -1,   cons.inter = 1,   confint.level = 0.95,   bs = "mpd",   kn = -1,  order.ispline = 2,
                     kn.ispline = 1,  distance.type = "Euclidean",  dist.custom = NULL,  rescale = FALSE,
                     rescale.pred = TRUE,  method = "mean",  normalize = FALSE,  silent = FALSE,  empty.row = 0,
                     control = list(),  glm.init = FALSE)

nat.g10 <- Zeta.msgdm(data.spec = native.pa, data.env = env[,c(16,4,10,8)], xy =  xydf,# data.spec.pred = alien.pa,
                      order = 10,  sam = 10000,  reg.type = "ispline",  family = stats::gaussian(),  method.glm = "glm.fit.cons",
                      cons = -1,   cons.inter = 1,   confint.level = 0.95,   bs = "mpd",   kn = -1,  order.ispline = 2,
                      kn.ispline = 1,  distance.type = "Euclidean",  dist.custom = NULL,  rescale = FALSE,
                      rescale.pred = TRUE,  method = "mean",  normalize = FALSE,  silent = FALSE,  empty.row = 0,
                      control = list(),  glm.init = FALSE)

zetmsgdm <- list(nat.g2 = nat.g2, nat.g4 = nat.g4, nat.g5 = nat.g5, nat.g10 = nat.g10,
                 ale.g2 = ale.g2, ale.g4 = ale.g4, ale.g5 = ale.g5, ale.g10 = ale.g10)
save(zetmsgdm, file=file.path("zetamsgdm.RData"))

par(mfrow=c(2,2))
png2 <- Plot.ispline(msgdm=nat.g2, data.env = env[,c(16,4,10,8)], distance=TRUE, num.quantiles=11)#, biotic=TRUE)
png4 <-Plot.ispline(msgdm=nat.g4, data.env=env[,c(16,4,10,8)], distance=TRUE)
png5 <- Plot.ispline(msgdm=nat.g5, data.env=env[,c(16,4,10,8)], distance=TRUE)
png10 <- Plot.ispline(msgdm=nat.g10, data.env=env[,c(16,4,10,8)], distance=TRUE)
# par(mfrow=c(2,2))
pag2 <- Plot.ispline(msgdm=ale.g2, data.env = env[,c(16,4,10,8)], distance=TRUE)#, biotic=TRUE)
pag4 <- Plot.ispline(msgdm=ale.g4, data.env=env[,c(16,4,10,8)], distance=TRUE)
pag5 <- Plot.ispline(msgdm=ale.g5, data.env=env[,c(16,4,10,8)], distance=TRUE)
pag10 <- Plot.ispline(msgdm=ale.g10, data.env=env[,c(16,4,10,8)], distance=TRUE)


# zeta varpart #####
vpa2 <- Zeta.varpart(ale.g2, reg.type="ispline")#[4:7,]
vpa4 <- Zeta.varpart(ale.g4, reg.type="ispline")#[4:7,]
vpa5 <- Zeta.varpart(ale.g5, reg.type="ispline")#[4:7,]
vpa10 <- Zeta.varpart(ale.g10, reg.type="ispline")#[4:7,]
vpn2 <- Zeta.varpart(nat.g2, reg.type="ispline")#[4:7,]
vpn4 <- Zeta.varpart(nat.g4, reg.type="ispline")#[4:7,]
vpn5 <- Zeta.varpart(nat.g5, reg.type="ispline")#[4:7,]
vpn10 <- Zeta.varpart(nat.g10, reg.type="ispline")#[4:7,]

vap <-cbind(vpa2, vpa4, vpa5, vpa10, vpn2, vpn4, vpn5, vpn10)
vapR2 <- round(vap[4:7,],3)
row.names(vapR2) <- c("xy","xy:env","env","unexpl")
names(vapR2) <- c("Alien Z2","Alien Z4","Alien Z5","Alien Z10",
                  "Native Z2","Native Z4","Native Z5","Native Z10")
write.csv(vapR2, file="varpartGDM.csv")

# ispline ~ raw ####
res.x <- function(x, xraw){
  xbar = x*(max(xraw)-min(xraw)) + min(xraw)
  return(xbar)
  }

iale2 <- Return.ispline(ale.g2, data.env=env[,c(16,4,10,8)], distance=TRUE)
inat2 <- Return.ispline(nat.g2, data.env=env[,c(16,4,10,8)], distance=TRUE)

iale4 <- Return.ispline(ale.g4, data.env=env[,c(16,4,10,8)], distance=TRUE)
inat4 <- Return.ispline(nat.g4, data.env=env[,c(16,4,10,8)], distance=TRUE)

iale5 <- Return.ispline(ale.g5, data.env=env[,c(16,4,10,8)], distance=TRUE)
inat5 <- Return.ispline(nat.g5, data.env=env[,c(16,4,10,8)], distance=TRUE)

iale10 <- Return.ispline(ale.g10, data.env=env[,c(16,4,10,8)], distance=TRUE)
inat10 <- Return.ispline(nat.g10, data.env=env[,c(16,4,10,8)], distance=TRUE)


IAxy.2 <- res.x(x = iale2$env[,5], xraw= dist(xydf)/1000)
INxy.2 <- res.x(x = inat2$env[,5], xraw= dist(xydf)/1000)

IAxy.4 <- res.x(x = iale4$env[,5], xraw= dist(xydf)/1000)
INxy.4 <- res.x(x = inat4$env[,5], xraw= dist(xydf)/1000)

IAxy.5 <- res.x(x = iale5$env[,5], xraw= dist(xydf)/1000)
INxy.5 <- res.x(x = inat5$env[,5], xraw= dist(xydf)/1000)

IAxy.10 <- res.x(x = iale10$env[,5], xraw= dist(xydf)/1000)
INxy.10 <- res.x(x = inat10$env[,5], xraw= dist(xydf)/1000)


IAecec.2 <- res.x(x = iale2$env[,1], xraw= env$ECEC)
INecec.2 <- res.x(x = inat2$env[,1], xraw= env$ECEC)

IAecec.4 <- res.x(x = iale4$env[,1], xraw= env$ECEC)
INecec.4 <- res.x(x = inat4$env[,1], xraw= env$ECEC)

IAecec.5 <- res.x(x = iale5$env[,1], xraw= env$ECEC)
INecec.5 <- res.x(x = inat5$env[,1], xraw= env$ECEC)

IAecec.10 <- res.x(x = iale10$env[,1], xraw= env$ECEC)
INecec.10 <- res.x(x = inat10$env[,1], xraw= env$ECEC)

IAlba.2 <- res.x(x = iale2$env[,2], xraw= env[,4])
INlba.2 <- res.x(x = inat2$env[,2], xraw= env[,4])

IAlba.4 <- res.x(x = iale4$env[,2], xraw= env[,4])
INlba.4 <- res.x(x = inat4$env[,2], xraw= env[,4])

IAlba.5 <- res.x(x = iale5$env[,2], xraw= env[,4])
INlba.5 <- res.x(x = inat5$env[,2], xraw= env[,4])

IAlba.10 <- res.x(x = iale10$env[,2], xraw= env[,4])
INlba.10 <- res.x(x = inat10$env[,2], xraw= env[,4])

IAom.2 <- res.x(x = iale2$env[,3], xraw= env[,10])
INom.2 <- res.x(x = inat2$env[,3], xraw= env[,10])
IAom.4 <- res.x(x = iale4$env[,3], xraw= env[,10])
INom.4 <- res.x(x = inat4$env[,3], xraw= env[,10])
IAom.5 <- res.x(x = iale5$env[,3], xraw= env[,10])
INom.5 <- res.x(x = inat5$env[,3], xraw= env[,10])
IAom.10 <- res.x(x = iale10$env[,3], xraw= env[,10])
INom.10 <- res.x(x = inat10$env[,3], xraw= env[,10])

IApH.2 <- res.x(x = iale2$env[,4], xraw= env[,8])
INpH.2 <- res.x(x = inat2$env[,4], xraw= env[,8])
IApH.4 <- res.x(x = iale4$env[,4], xraw= env[,8])
INpH.4 <- res.x(x = inat4$env[,4], xraw= env[,8])
IApH.5 <- res.x(x = iale5$env[,4], xraw= env[,8])
INpH.5 <- res.x(x = inat5$env[,4], xraw= env[,8])
IApH.10 <- res.x(x = iale10$env[,4], xraw= env[,8])
INpH.10 <- res.x(x = inat10$env[,4], xraw= env[,8])

# ggplot raw ispline ####
library(ggplot2)
dxy <- c(IAxy.2, IAxy.4, IAxy.5, IAxy.10, INxy.2, INxy.4, INxy.5, INxy.10)
isxy <- c(iale2$Ispline[,5],iale4$Ispline[,5],iale5$Ispline[,5],iale10$Ispline[,5],
          inat2$Ispline[,5],inat4$Ispline[,5],inat5$Ispline[,5],inat10$Ispline[,5])

isxy <- c(iale2$Ispline[,5],iale4$Ispline[,5],iale5$Ispline[,5],iale10$Ispline[,5],
          inat2$Ispline[,5],inat4$Ispline[,5],inat5$Ispline[,5],inat10$Ispline[,5])


ph <- c(IApH.2,IApH.4,IApH.5,IApH.10, INpH.2,INpH.4,INpH.5,INpH.10)
isph <- c(iale2$Ispline[,4],iale4$Ispline[,4],iale5$Ispline[,4],iale10$Ispline[,4],
          inat2$Ispline[,4],inat4$Ispline[,4],inat5$Ispline[,4],inat10$Ispline[,4])

ec <- c(IAecec.2,IAecec.4,IAecec.5,IAecec.10, INecec.2,INecec.4,INecec.5,INecec.10)
isec <- c(iale2$Ispline[,1],iale4$Ispline[,1],iale5$Ispline[,1],iale10$Ispline[,1],
          inat2$Ispline[,1],inat4$Ispline[,1],inat5$Ispline[,1],inat10$Ispline[,1])

om <- c(IAom.2,IAom.4,IAom.5,IAom.10, INom.2,INom.4,INom.5,INom.10)
isom <- c(iale2$Ispline[,3],iale4$Ispline[,3],iale5$Ispline[,3],iale10$Ispline[,3],
          inat2$Ispline[,3],inat4$Ispline[,3],inat5$Ispline[,3],inat10$Ispline[,3])

lba <- c(IAlba.2, IAlba.4, IAlba.5, IAlba.10, INlba.2, INlba.4, INlba.5, INlba.10)
islba <- c(iale2$Ispline[,2], iale4$Ispline[,2], iale5$Ispline[,2], iale10$Ispline[,2],
           inat2$Ispline[,2], inat4$Ispline[,2], inat5$Ispline[,2], inat10$Ispline[,2])

ords <- c(rep(2,50),rep(4,50),rep(5,50),rep(10,50))
comps <- c(rep("alien", 200),rep("native", 200))
var <- c(rep("ec", 400), rep("lba", 400), rep("om", 400), rep("ph", 400), rep("xy", 400))
plis <- data.frame(var=var, zeta=ords, veg=comps, env = c(ec,lba,om,ph,dxy), mod = c(isec, islba, isom,isph,isxy))
plis$var <- factor(plis$var, labels=c("ECEC", "Live basal area", "Organic matter","Soil pH", "Distance"))
plis$zeta <- factor(plis$zeta, labels=c("Zeta 2", "Zeta 4", "Zeta 5", "Zeta 10"))


isplt <- ggplot(aes(x = env, y= mod, group = veg),data=plis)+
  geom_line(aes(linetype=veg))+
  theme_minimal()+
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        panel.background = element_rect())+
  ylab("Modelled value (I spline)")+
  xlab("")+
  facet_grid(zeta ~ var, scales="free"); print(isplt)
  
library(grid)
g <- ggplotGrob(isplt)
xax <- which(g$layout$name=="xlab-b")
g[["grobs"]][[xax]]$children[[1]]$label <- c("ECEC (cmol^2/kg)","Area (m2)", "% OM", "pH","Distance (km)")
g[["grobs"]][[xax]]$children[[1]]$x <- grid::unit(c(0.1, 0.3, 0.5, 0.7,0.9), "npc")
#g[["grobs"]][[xax]]$children[[1]]$x <- grid::unit(c(0.12,0.37,0.62,0.87), "npc")
grid.draw(g)
ggsave(g, file="isplot.png", height=130,width=170,units="mm",dpi=200)

# ispline Fig 5 ####
isxy.x <- c(iale2$env[,5],iale4$env[,5],iale5$env[,5],iale10$env[,5],
          inat2$env[,5],inat4$env[,5],inat5$env[,5],inat10$env[,5])

isph.x <- c(iale2$env[,4],iale4$env[,4],iale5$env[,4],iale10$env[,4],
            inat2$env[,4],inat4$env[,4],inat5$env[,4],inat10$env[,4])

isec.x <- c(iale2$env[,1],iale4$env[,1],iale5$env[,1],iale10$env[,1],
            inat2$env[,1],inat4$env[,1],inat5$env[,1],inat10$env[,1])

isom.x <- c(iale2$env[,3],iale4$env[,3],iale5$env[,3],iale10$env[,3],
            inat2$env[,3],inat4$env[,3],inat5$env[,3],inat10$env[,3])

islba.x <- c(iale2$env[,2],iale4$env[,2],iale5$env[,2],iale10$env[,2],
            inat2$env[,2],inat4$env[,2],inat5$env[,2],inat10$env[,2])

ords <- c(rep(2,50),rep(4,50),rep(5,50),rep(10,50))
comps2 <- c(rep("alien", 200),rep("native", 200))
var2 <- c(rep("xy", 400), rep("ph", 400), rep("ec", 400), rep("om", 400), rep("lba", 400))


ispandf <- data.frame(x = c(isxy.x, isph.x, isec.x, isom.x, islba.x), y =c(isxy, isph, isec, isom, islba),
                      zeta = ords, comp = comps2, var = var2)
ispandf$var <- factor(ispandf$var, levels = c("xy", "ph", "ec","om","lba"),
                                              labels=c("Distance","Soil pH", "ECEC", "Organic matter", "Live basal area"))
ispandf$zeta <- factor(ispandf$zeta, labels=c("Zeta 2", "Zeta 4", "Zeta 5", "Zeta 10"))
ispandf$comp <- factor(ispandf$comp,levels=c("native","alien"), labels=c("Native", "Alien"))

# pts data 
# pp <- pts.hack(isplines=iale2)
zlst <- list(iale2, iale4, iale5, iale10, inat2, inat4, inat5, inat10)
zpt <- lapply(zlst, function(x) pts.hack(x))

xlst <-(lapply(zpt, function(z) rbind(z$x)))
ylst <- (lapply(zpt, function(z) rbind(z$y)))

xpt <- cbind(unlist(xlst))
ypt <- cbind(unlist(ylst))
noms <- rep(names(iale2$env),each=11)
ordz <- c(rep(2,55),rep(4,55),rep(5,55),rep(10,55))
anom <- rep(c("alien","native"), each=55*4)

noms[noms== "ECEC"] <- "ec" 
noms[noms== "Live.basal.area"] <- "lba" 
noms[noms== "Organic.Matter"] <- "om"
noms[noms== "distance"] <- "xy"
noms[noms== "pH"] <- "ph"

pdat <- data.frame(xpt = xpt, ypt=ypt, var = noms, comp = anom, zeta = rep(ordz,2))
pdat$var <- factor(pdat$var, levels = c("xy", "ph", "ec","om","lba"),
                   labels=c("Distance","Soil pH", "ECEC", "Organic matter", "Live basal area"))      
pdat$zeta <- factor(pdat$zeta, labels=c("Zeta 2", "Zeta 4", "Zeta 5", "Zeta 10"))
pdat$comp <- factor(pdat$comp,levels=c("native","alien"), labels=c("Native", "Alien"))

# ispandf$var <- factor(ispandf$var, levels = c("xy", "ph", "ec","om","lba"),
#                       labels=c("Distance","Soil pH", "ECEC", "Organic matter", "Live basal area"))
# ispandf$zeta <- factor(ispandf$zeta, labels=c("Zeta 2", "Zeta 4", "Zeta 5", "Zeta 10"))
# ispandf$comp <- factor(ispandf$comp,levels=c("native","alien"), labels=c("Native", "Alien"))

labs <- data.frame(x = rep(0.02,8), y = c(9,9,10,10,8,8,6,6), lab=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"),
                   comp = rep(c("native","alien"), 4), zeta = rep(c("Zeta 2","Zeta 4","Zeta 5","Zeta 10"), each=2), 
                   var=rep("Distance",8))

labs$zeta <- factor(labs$zeta, levels = c("Zeta 2","Zeta 4","Zeta 5","Zeta 10"), 
                    labels=c("Zeta 2", "Zeta 4", "Zeta 5", "Zeta 10"))
labs$comp <- factor(labs$comp,levels=c("native","alien"), labels=c("Native", "Alien"))

cbPalette <- c("#0072B2", "#009E73","#E69F00", "#999999", "#000000", "#F0E442",  "#D55E00", "#CC79A7","#56B4E9")

ispan.fig <- ggplot(aes(x=x, y = y, group=interaction(comp,var,zeta)), data=ispandf)+
  geom_line(aes(colour=var,linetype=var))+
  geom_point(aes(x=xpt,y=ypt,group=interaction(comp,var,zeta), colour=var, pch=var), fill="white", data=pdat)+ #, fill="white", size=2
  theme_minimal()+
  scale_colour_manual(values = c("#0072B2", "#009E73","#E69F00", "#999999", "#000000"))+
  scale_shape_manual(values=c(21,22,23,24,25))+
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        panel.background = element_rect())+
  ylab("Modelled value (I spline)")+
  xlab("Rescaled value")+
  geom_text(aes(x=x,y=y,label = lab),data=labs)+
  facet_grid(zeta ~ comp, scales="free_y"); print(ispan.fig)

ggsave(ispan.fig, file="ispanel.png", height=170,width=140,units="mm",dpi=200)

# Native community #

# Using data file plant_community 
native_plant_community.2 <- plant_community

# Data frame indexing to just look at the alien community
native_plant_community.2 <- native_plant_community.2[native_plant_community.2$Origin == 'native',]

# Removes the species ID and origin columns
native_plant_community.2 <- native_plant_community.2[-c(1:2) ]


# Sums the columns in the biodiversity matrix to get species richness per plot #
native_plant_community_site_rich <- as.data.frame(colSums(native_plant_community.2))
native_plant_community_site_rich
names(native_plant_community_site_rich) <- c("siterichness")


#add column
native_plant_community_site_rich$Origin <- "native"



# Combine the alien and native site richness dataframes
combined_site_rich <- rbind(alien_plant_community_site_rich,native_plant_community_site_rich)

# Interleaved histograms #
separated.plot <- ggplot(combined_site_rich, aes(x=siterichness, fill=Origin)) +
  geom_histogram(breaks=seq(0, 60, by = 5), position = "dodge", col = "black") +
  theme_bw(base_family = "sans") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position =  "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.line = element_line(colour = "black", size = .5)) +
  scale_fill_manual(values=c("alien"="lightcoral", "native"="mediumseagreen")) +
  scale_x_continuous(name = "Species Richness", expand = c(0, 0)) + 
  scale_y_continuous(name = "Frequency",expand = c(0, 0), breaks = seq(0, 20, 5), limits = c(0, 20))

separated.plot

ggsave(filename=("Outputs/Figures/Exploratory/colour_Interleavered plant community_site rich.png"),
       width = 5, height = 4.5, dpi = 700)

theme_set(theme_cowplot(font_size=14, font_family = "sans"))
plot_grid(combined.plot,separated.plot,ncol=2,
          align='hv',labels=c('a)','b)'), hjust = -1.5, label_size = 16)

ggsave(filename=("Outputs/Figures/Exploratory/To use/Colour combined and interleavered plant community_site rich.png"))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#    #~# Rank occupancy plots for the alien and native communities
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Alien community #

# Using data file plant_community 
alien_plant_community.3 <- plant_community

# Data frame indexing to just look at the alien community
alien_plant_community.3 <- alien_plant_community.3[alien_plant_community.3$Origin == 'alien',]
str(alien_plant_community.3)

# Converting the data from integer to numeric 
alien_plant_community.3[3:51] <- lapply(alien_plant_community.3[3:51], as.numeric)
str(alien_plant_community.3)

# Removes the species ID and origin columns and then the column names
alien_plant_community.3 <- alien_plant_community.3[-c(1:2) ]
colnames(alien_plant_community.3) <-NULL
alien_plant_community.3

# Sums the rows to get the frequency of occurance for each species across the plots
alien_rank <- as.data.frame(rowSums(alien_plant_community.3))

names(alien_rank) <- c("rank")

#add columns and sort teh data into descending values
alien_rank$sp.identity <- 1:nrow(alien_rank) 
alien_rank <- alien_rank[rev(order(alien_rank$rank)),]
alien_rank$sp.identity <- 1:nrow(alien_rank) 
str(alien_rank)

alien_rank[,2] <- as.numeric(as.character(alien_rank[,2]))

# Plotting alien rank frequency
alien_rank_plot <- ggplot(alien_rank, aes(x = reorder(sp.identity, -rank), y = rank)) +
  geom_bar(stat= "identity", colour="white", fill="lightcoral", size = 0.5) +
  theme_bw(base_family = "sans") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.line = element_line(colour = "black", size = .5)) +
  scale_x_discrete(name = "Rank", expand = c(0, 0), breaks = seq(0, 80, 10)) + 
  scale_y_continuous(name = "Frequency",expand = c(0, 0), breaks = seq(0, 60, 10), limits = c(0, 60))




# Native community #

# Using data file plant_community 
native_plant_community.3 <- plant_community

# Data frame indexing to just look at the native community
native_plant_community.3 <- native_plant_community.3[native_plant_community.3$Origin == 'native',]
str(native_plant_community.3)

# Converting the data from integer to numeric 
native_plant_community.3[3:51] <- lapply(native_plant_community.3[3:51], as.numeric)
str(native_plant_community.3)

# Removes the species ID and origin columns and then the column names
native_plant_community.3 <- native_plant_community.3[-c(1:2) ]
colnames(native_plant_community.3) <-NULL
native_plant_community.3

# Sums the rows to get the frequency of occurance for each species across the plots
native_rank <- as.data.frame(rowSums(native_plant_community.3))

names(native_rank) <- c("rank")

#add columns and sort teh data into descending values
native_rank$sp.identity <- 1:nrow(native_rank) 
native_rank <- native_rank[rev(order(native_rank$rank)),]
native_rank$sp.identity <- 1:nrow(native_rank) 
str(native_rank)

native_rank[,2] <- as.numeric(as.character(native_rank[,2]))

# Plotting native rank frequency
native_rank_plot <- ggplot(native_rank, aes(x = reorder(sp.identity, -rank), y = rank)) +
  geom_bar(stat= "identity", colour="white", fill="mediumseagreen") +
  theme_bw(base_family = "sans") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "none",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.line = element_line(colour = "black", size = .5)) +
  scale_x_discrete(name = "Rank", expand = c(0, 0), breaks = seq(0, 80, 10)) + 
  scale_y_continuous(name = "Frequency",expand = c(0, 0), breaks = seq(0, 60, 10), limits = c(0, 60))

# Combine the alien and native rank frequency plots

# add origin column
alien_rank$Origin <- "alien"

#add origin column
native_rank$Origin <- "native"


# Adding rows to make the alien and native 'ranks' the same number
# Create dataframe 
rank <- c(rep(0,103))
Origin <- c(rep("alien", 103))
sp.identity <- c(seq(72,174, by = 1))

alien.additions <- data.frame(rank,Origin, sp.identity)

# Combine dataframes by row
alien_rank <- rbind(alien_rank,alien.additions) 


# Combine alien and native dataframes by row
combined_rank <- rbind(alien_rank,native_rank)

str(combined_rank)


# Interleaved histograms #
interleavered_rank_plot <- ggplot(combined_rank, aes(x = sp.identity, y = rank, fill=Origin)) +
  geom_bar(stat= "identity", position = "dodge", width = 0.9, colour="white", size = 0.1) +
  theme_bw(base_family = "sans") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.line = element_line(colour = "black", size = .5)) +
  scale_fill_manual(values=c("alien"="lightcoral", "native"="mediumseagreen")) +
  scale_x_continuous(name = "Individual species rank", expand = c(0, 0), breaks = seq(0, 80, 20), limits = c(0, 80)) + 
  scale_y_continuous(name = "Frequency",expand = c(0, 0), breaks = seq(0, 50, 10), limits = c(0, 50))
interleavered_rank_plot


ggsave(filename=("Outputs/Figures/Exploratory/50 plots averaged HS7/interleavered rank plot_colour.2.png"),
       width = 5, height = 4.5, dpi = 700)



#~# Combining all exploratory plots (so far) for one figure
theme_set(theme_cowplot(font_size=12, font_family = "sans"))
p1 <- plot_grid(combined.plot,separated.plot,combined, interleavered, ncol=2,
                align='hv',labels=c('a)','b)', 'c)', 'd)'), hjust = -1, label_size = 14)+
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

p1


p2 <- plot_grid(native.rich.alien.rich,interleavered_rank_plot, ncol=2, align='hv',labels=c('c)','f)'), 
                hjust = -1, label_size = 12) + 
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

#p3 <- plot_grid(native.rich.alien.rich, ncol=2, align='hv',labels=c('f)'), hjust = -1, label_size = 14)

plot_grid(p1, p2,  ncol = 1, align = 'hv',rel_heights = c(2,1,1)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

ggsave(filename=("Outputs/Figures/Exploratory/To use/V3newest_Combined site richness, frequency occupancy and rank frequency plots.png"), 
       width = 9, height = 11, dpi = 500)

#~# Combining all exploratory plots for manuscript figure
theme_set(theme_cowplot(font_size=12, font_family = "sans"))
p1 <- plot_grid(combined.plot, interleavered_rank_plot, separated.plot, combined, ncol=2, axis = "bt",
                align='hv',labels=c('a)','d)', 'b)', 'e)'), hjust = -1, label_size = 14)+
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

p1

p2 <- plot_grid(native.rich.alien.rich,percent.stacked , ncol=2, align='hv',labels=c('c)','f)'), 
                axis = "bt", hjust = -1, label_size = 12) + 
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

p2

p3 <- plot_grid(p1, p2,  ncol = 1, align = 'hv',rel_heights = c(2,1,1)) +
                theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))


legend_b <- get_legend(percent.stacked + theme(legend.position="bottom"))
# 
# # add the legend underneath the row we made earlier. Give it 10% of the height
# # of one plot (via rel_heights).
p_with_legend <- plot_grid(legend_b, p3,  ncol = 1, rel_heights = c(0.05, 1))
p_with_legend

ggsave(filename=("Outputs/Figures/Exploratory/To use/V1_for manuscript_Combined site richness, frequency occupancy and rank frequency plots.png"), 
       width = 11, height = 12, dpi = 700)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#     Standard rarefaction curve analysis                                            
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Entire community #

# Using the data file plant_community
# Rows are plant species, columns are plot IDs

entire_plant_community.4 <- plant_community
str(entire_plant_community.4)

# Removes the species ID and origin columns and column names
entire_plant_community.4 <- entire_plant_community.4[-c(1:2) ]
colnames(entire_plant_community.4) <- NULL

# Transposing the table #
entire_plant_community.4.mat <- t(as.matrix(entire_plant_community.4))
entire_plant_community.4.mat

# Species accumulation curve # 
entire_plant_community.4_sac <- specaccum(entire_plant_community.4.mat, method ="random", permutations = 1000)
plot(entire_plant_community.4_sac, ci.type="line", xlab =" Number of sites", ylab = "Number of species", cex.lab=1.5)
# Note: Curve does not reach asymptote

# Extrapolated richness #
specpool(entire_plant_community.4.mat)
245/290


# Alien community #

# Using data file plant_community 
alien_plant_community.4 <- plant_community

# Data frame indexing to just look at the alien community
alien_plant_community.4 <- alien_plant_community.4[alien_plant_community.4$Origin == 'alien',]
str(alien_plant_community.4)

# removing species identity and plotID - now we have just incidence data #
alien_plant_community.4 <- alien_plant_community.4[-c(1:2) ]
colnames(alien_plant_community.4) <- NULL


# Transposing the table #
alien_plant_community.4.mat <- t(as.matrix(alien_plant_community.4))
alien_plant_community.4.mat

# Species accumulation curve # - default is 100 permutations
alien_plant_community.4_sac <- specaccum(alien_plant_community.4.mat, method ="random", permutations = 1000)
plot(alien_plant_community.4_sac, ci.type="line", xlab =" Number of sites", ylab = "Number of species", cex.lab=1.5)
summary(alien_plant_community.4_sac)
# Note: Curve does not reach asymptote

# Extrapolated richness #
specpool(alien_plant_community.4.mat) 
71/84



# Native community #

# Using data file plant_community 
native_plant_community.4 <- plant_community

# Data frame indexing to just look at the alien community
native_plant_community.4 <- native_plant_community.4[native_plant_community.4$Origin == 'native',]
str(native_plant_community.4)

# removing species identity and plotID - now we have just incidence data #
native_plant_community.4 <- native_plant_community.4[-c(1:2) ]
colnames(native_plant_community.4) <- NULL


# Transposing the table #
native_plant_community.4.mat <- t(as.matrix(native_plant_community.4))
native_plant_community.4.mat

# Species accumulation curve # - default is 100 permutations
native_plant_community.4_sac <- specaccum(native_plant_community.4.mat, method ="random", permutations = 1000)
plot(native_plant_community.4_sac, ci.type="line", xlab =" Number of sites", ylab = "Number of species", cex.lab=1.5)
summary(native_plant_community.4_sac)
# Note: Curve does not reach asymptote

# Extrapolated richness #
specpool(native_plant_community.4.mat) 
174/206


#~# Plotting native and alien rarefaction curves together
par(mfrow=c(1,1), family = "sans",mai = c(1, 1, 1, 1))
plot(entire_plant_community.4_sac,  xlab =" Number of sites", ylab = "Number of species (richness)")
plot(alien_plant_community.4_sac, add = TRUE, col = 2) # col = 2 is green
plot(native_plant_community.4_sac, add = TRUE, col = 3) # col = 3 is red

# Save the plot to: Outputs/1. Plant Community/Figures/Rarefaction/Core combined standard rarefaction.png









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Inverse (switched absences and presences)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Entire community #

# Using the data file inverse_plant_community
# Rows are inverse_plant species, columns are plot IDs

entire_inverse_plant_community.4 <- inverse_plant_community
str(entire_inverse_plant_community.4)

# Removes the species ID and origin columns and column names
entire_inverse_plant_community.4 <- entire_inverse_plant_community.4[-c(1:2) ]
colnames(entire_inverse_plant_community.4) <- NULL

# Transposing the table #
entire_inverse_plant_community.4.mat <- t(as.matrix(entire_inverse_plant_community.4))
entire_inverse_plant_community.4.mat

# Species accumulation curve # 
entire_inverse_plant_community.4_sac <- specaccum(entire_inverse_plant_community.4.mat, method ="random", permutations = 1000)
plot(entire_inverse_plant_community.4_sac, ci.type="line", xlab =" Number of sites", ylab = "Number of species", cex.lab=1.5)
# Note: Curve does not reach asymptote

# Extrapolated richness #
specpool(entire_inverse_plant_community.4.mat)
245/290


# Alien community #

# Using data file inverse_plant_community 
alien_inverse_plant_community.4 <- inverse_plant_community

# Data frame indexing to just look at the alien community
alien_inverse_plant_community.4 <- alien_inverse_plant_community.4[alien_inverse_plant_community.4$Origin == 'alien',]
str(alien_inverse_plant_community.4)

# removing species identity and plotID - now we have just incidence data #
alien_inverse_plant_community.4 <- alien_inverse_plant_community.4[-c(1:2) ]
colnames(alien_inverse_plant_community.4) <- NULL


# Transposing the table #
alien_inverse_plant_community.4.mat <- t(as.matrix(alien_inverse_plant_community.4))
alien_inverse_plant_community.4.mat

# Species accumulation curve # - default is 100 permutations
alien_inverse_plant_community.4_sac <- specaccum(alien_inverse_plant_community.4.mat, method ="random", permutations = 1000)
plot(alien_inverse_plant_community.4_sac, ci.type="line", xlab =" Number of sites", ylab = "Number of species", cex.lab=1.5)
summary(alien_inverse_plant_community.4_sac)
# Note: Curve does not reach asymptote

# Extrapolated richness #
specpool(alien_inverse_plant_community.4.mat) 
71/84



# Native community #

# Using data file inverse_plant_community 
native_inverse_plant_community.4 <- inverse_plant_community

# Data frame indexing to just look at the alien community
native_inverse_plant_community.4 <- native_inverse_plant_community.4[native_inverse_plant_community.4$Origin == 'native',]
str(native_inverse_plant_community.4)

# removing species identity and plotID - now we have just incidence data #
native_inverse_plant_community.4 <- native_inverse_plant_community.4[-c(1:2) ]
colnames(native_inverse_plant_community.4) <- NULL


# Transposing the table #
native_inverse_plant_community.4.mat <- t(as.matrix(native_inverse_plant_community.4))
native_inverse_plant_community.4.mat

# Species accumulation curve # - default is 100 permutations
native_inverse_plant_community.4_sac <- specaccum(native_inverse_plant_community.4.mat, method ="random", permutations = 1000)
plot(native_inverse_plant_community.4_sac, ci.type="line", xlab =" Number of sites", ylab = "Number of species", cex.lab=1.5)
summary(native_inverse_plant_community.4_sac)
# Note: Curve does not reach asymptote

# Extrapolated richness #
specpool(native_inverse_plant_community.4.mat) 
174/206

#~# Plotting native and alien rarefaction curves together
par(mfrow=c(1,1), family = "sans",mai = c(1, 1, 1, 1))
plot(entire_inverse_plant_community.4_sac,  xlab="Number of sites", 
     ylab="Number of species (richness)")
plot(entire_plant_community.4_sac, add = TRUE)
plot(alien_inverse_plant_community.4_sac, add = TRUE, col = 2) # col = 2 is green
plot(alien_plant_community.4_sac, add = TRUE, col = 2) # col = 2 is green
plot(native_inverse_plant_community.4_sac, add = TRUE, col = 3) # col = 3 is red
plot(native_plant_community.4_sac, add = TRUE, col = 3) # col = 3 is red

plot(native_plant_community.4_sac, add = TRUE, col = 3) # col = 3 is red

# Save the plot to: Outputs/1. inverse_plant Community/Figures/Rarefaction/Core combined standard rarefaction.png



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#     Spatially constrained rarefaction curve analysis                                            
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Entire community #

# Using the data file plant_community
# Rows are plant species, columns are plot IDs

entire_plant_community.5 <- plant_community
str(entire_plant_community.5)

# Removes the species ID and origin columns and column names
entire_plant_community.5 <- entire_plant_community.5[-c(2) ]
#colnames(entire_plant_community.5) <- NULL

# Change the shape of the community data frame
entire_plant_community.5.t <- as.data.frame(t(entire_plant_community.5[,-c(1)]))  
entire_plant_community.5.t <- cbind(ID=row.names(entire_plant_community.5.t),entire_plant_community.5.t)
head(entire_plant_community.5.t)

# Pointpattern from Bacaro et al 2012 method
core_mxp <- pointpattern(plot_coords, longlat=TRUE)

# Spatial rarefaction from Bacaro et al 2012 method
SCR_curve_all <- SCR(entire_plant_community.5.t, core_mxp)

# the mxp matrix represents the output of the pointpattern function, see Appendix 1 non_spatially_explicit<-specaccum(community_data[2:ncol(community_data)], method="exact")
# The following code can be used to generate Figure 3 in Bacaro et al. (2012).
# In this example the curves are not different from each other since the community
# data being used are synthetic (see community_data)
# Refer to Figure 3 in Bacaro et al. (2012) and to Palmer et al. (2007)
# for a published dataset on vascular plants from a North Carolina piedmont forest. 

plot(seq(1:50),SCR_curve_all[,1], xlim=c(0,50), ylim=c(0,300), xlab="Number of plots", ylab="Species Richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_all[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,3], lty=1, lwd=3, col="dark grey") 


# Alien community #

# Using data file plant_community 
alien_plant_community.5 <- plant_community

# Data frame indexing to just look at the alien community
alien_plant_community.5 <- alien_plant_community.5[alien_plant_community.5$Origin == 'alien',]
str(alien_plant_community.5)

# Removes the origin column
alien_plant_community.5 <- alien_plant_community.5[-c(2) ]

# Change the shape of the community data frame
alien_plant_community.5.t <- as.data.frame(t(alien_plant_community.5[,-c(1)]))  
alien_plant_community.5.t <- cbind(ID=row.names(alien_plant_community.5.t),alien_plant_community.5.t)
head(alien_plant_community.5.t)

# Pointpattern from Bacaro et al 2012 method
core_mxp <- pointpattern(plot_coords, longlat=TRUE)

# Spatial rarefaction from Bacaro et al 2012 method
SCR_curve_alien <- SCR(alien_plant_community.5.t, core_mxp)


# the mxp matrix represents the output of the pointpattern function, see 
# Appendix 1 non_spatially_explicit<-specaccum(community_data[2:ncol(community_data)], method="exact")
# The following code can be used to generate Figure 3 in Bacaro et al. (2012).
# In this example the curves are not different from each other since the community
# data being used are synthetic (see community_data)
# Refer to Figure 3 in Bacaro et al. (2012) and to Palmer et al. (2007)
# for a published dataset on vascular plants from a North Carolina piedmont forest. 

plot(seq(1:50),SCR_curve_alien[,1], xlim=c(0,50), ylim=c(0,80), xlab="Number of plots", 
     ylab="Species Richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_alien[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_alien[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_alien[,3], lty=1, lwd=3, col="dark grey")

# Native community #

# Using data file plant_community 
native_plant_community.5 <- plant_community

# Data frame indexing to just look at the alien community
native_plant_community.5 <- native_plant_community.5[native_plant_community.5$Origin == 'native',]
str(native_plant_community.5)

# Removes the origin column
native_plant_community.5 <- native_plant_community.5[-c(2) ]

# Change the shape of the community data frame
native_plant_community.5.t <- as.data.frame(t(native_plant_community.5[,-c(1)]))  
native_plant_community.5.t <- cbind(ID=row.names(native_plant_community.5.t),native_plant_community.5.t)
head(native_plant_community.5.t)

# Pointpattern from Bacaro et al 2012 method
#core_mxp <- pointpattern(core_plot_coords, longlat=TRUE)

# Spatial rarefaction from Bacaro et al 2012 method
SCR_curve_nat <- SCR(native_plant_community.5.t, core_mxp)

# the mxp matrix represents the output of the pointpattern function, 
# see Appendix 1 non_spatially_explicit<-specaccum(community_data[2:ncol(community_data)], method="exact")
# The following code can be used to generate Figure 3 in Bacaro et al. (2012).
# In this example the curves are not different from each other since the community
# data being used are synthetic (see community_data)
# Refer to Figure 3 in Bacaro et al. (2012) and to Palmer et al. (2007)
# for a published dataset on vascular plants from a North Carolina piedmont forest. 

plot(seq(1:50),SCR_curve_nat[,1], xlim=c(0,50), ylim=c(0,200), xlab="Number of plots", 
     ylab="Species Richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_nat[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_nat[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_nat[,3], lty=1, lwd=3, col="dark grey")



# Combining the three graphs
par(mfrow=c(1,1))
plot(seq(1:50),SCR_curve_all[,1], xlim=c(0,50), ylim=c(0,300), 
     xlab="Number of sites", ylab="Number of species (richness)", type="p", pch=10, ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_all[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,3], lty=1, lwd=3, col="dark grey") 


lines(seq(1:50),SCR_curve_alien[,1], lwd=3, col="red") 
lines(seq(1:50),SCR_curve_alien[,2], lty=1 ,lwd=3, col="red") 
lines(seq(1:50),SCR_curve_alien[,3], lty=1, lwd=3, col="red")

lines(seq(1:50),SCR_curve_nat[,1], lwd=3, col="green") 
lines(seq(1:50),SCR_curve_nat[,2], lty=1 ,lwd=3, col="green") 
lines(seq(1:50),SCR_curve_nat[,3], lty=1, lwd=3, col="green")



#~# Combining graphs according to their community type 
par(mfrow=c(3,1), family = "sans",mai = c(0.7, 0.75, 0.15, 0.5))
#~# Entire plant community comparison
plot(seq(1:50),SCR_curve_all[,1], xlim=c(0,50), ylim=c(0,300), cex.axis = 1.5, cex.lab = 2,
     xlab="Number of plots", ylab="Species richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_all[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,3], lty=1, lwd=3, col="dark grey") 
plot(entire_plant_community.4_sac, ci.type="line", lwd= 3, add = TRUE, xlab ="Number of sites", ylab = "Number of species", cex.lab=1.5)

#~# Alien plant community comparison
plot(seq(1:50),SCR_curve_alien[,1], xlim=c(0,50), ylim=c(0,80), cex.axis = 1.5, cex.lab = 2,
     xlab="Number of plots", ylab="Species richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_alien[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_alien[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_alien[,3], lty=1, lwd=3, col="dark grey") 
plot(alien_plant_community.4_sac, ci.type="line", lwd= 3, add = TRUE, xlab ="Number of sites", 
     ylab = "Number of species", cex.lab=1.5)

#~# Native plant community comparison
plot(seq(1:50),SCR_curve_nat[,1], xlim=c(0,50), ylim=c(0,200), cex.axis = 1.5, cex.lab = 2,
     xlab="Number of plots", ylab="Species richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_nat[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_nat[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_nat[,3], lty=1, lwd=3, col="dark grey") 
plot(native_plant_community.4_sac, ci.type="line", lwd= 3, add = TRUE, xlab ="Number of sites", ylab = "Number of species", cex.lab=1.5)

# Save to (filename=("Outputs/Figures/Exploratory/To use/Combined standard and SC rarefaction curves.png"))







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Inverse (switched absences and presences)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Entire community #

# Using the data file inverse_plant_community
# Rows are inverse_plant species, columns are plot IDs

entire_inverse_plant_community.5 <- inverse_plant_community
str(entire_inverse_plant_community.5)

# Removes the species ID and origin columns and column names
entire_inverse_plant_community.5 <- entire_inverse_plant_community.5[-c(2) ]
#colnames(entire_inverse_plant_community.5) <- NULL

# Change the shape of the community data frame
entire_inverse_plant_community.5.t <- as.data.frame(t(entire_inverse_plant_community.5[,-c(1)]))  
entire_inverse_plant_community.5.t <- cbind(ID=row.names(entire_inverse_plant_community.5.t),entire_inverse_plant_community.5.t)
head(entire_inverse_plant_community.5.t)

# Pointpattern from Bacaro et al 2012 method
core_mxp <- pointpattern(plot_coords, longlat=TRUE)

# Spatial rarefaction from Bacaro et al 2012 method
SCR_curve_all <- SCR(entire_inverse_plant_community.5.t, core_mxp)

# the mxp matrix represents the output of the pointpattern function, see Appendix 1 non_spatially_explicit<-specaccum(community_data[2:ncol(community_data)], method="exact")
# The following code can be used to generate Figure 3 in Bacaro et al. (2012).
# In this example the curves are not different from each other since the community
# data being used are synthetic (see community_data)
# Refer to Figure 3 in Bacaro et al. (2012) and to Palmer et al. (2007)
# for a published dataset on vascular inverse_plants from a North Carolina piedmont forest. 

plot(seq(1:50),SCR_curve_all[,1], xlim=c(0,50), ylim=c(0,300), xlab="Number of plots", ylab="Species Richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_all[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,3], lty=1, lwd=3, col="dark grey") 


# Alien community #

# Using data file inverse_plant_community 
alien_inverse_plant_community.5 <- inverse_plant_community

# Data frame indexing to just look at the alien community
alien_inverse_plant_community.5 <- alien_inverse_plant_community.5[alien_inverse_plant_community.5$Origin == 'alien',]
str(alien_inverse_plant_community.5)

# Removes the origin column
alien_inverse_plant_community.5 <- alien_inverse_plant_community.5[-c(2) ]

# Change the shape of the community data frame
alien_inverse_plant_community.5.t <- as.data.frame(t(alien_inverse_plant_community.5[,-c(1)]))  
alien_inverse_plant_community.5.t <- cbind(ID=row.names(alien_inverse_plant_community.5.t),alien_inverse_plant_community.5.t)
head(alien_inverse_plant_community.5.t)

# Pointpattern from Bacaro et al 2012 method
core_mxp <- pointpattern(plot_coords, longlat=TRUE)

# Spatial rarefaction from Bacaro et al 2012 method
SCR_curve_alien <- SCR(alien_inverse_plant_community.5.t, core_mxp)


# the mxp matrix represents the output of the pointpattern function, see 
# Appendix 1 non_spatially_explicit<-specaccum(community_data[2:ncol(community_data)], method="exact")
# The following code can be used to generate Figure 3 in Bacaro et al. (2012).
# In this example the curves are not different from each other since the community
# data being used are synthetic (see community_data)
# Refer to Figure 3 in Bacaro et al. (2012) and to Palmer et al. (2007)
# for a published dataset on vascular inverse_plants from a North Carolina piedmont forest. 

plot(seq(1:50),SCR_curve_alien[,1], xlim=c(0,50), ylim=c(0,80), xlab="Number of plots", 
     ylab="Species Richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_alien[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_alien[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_alien[,3], lty=1, lwd=3, col="dark grey")

# Native community #

# Using data file inverse_plant_community 
native_inverse_plant_community.5 <- inverse_plant_community

# Data frame indexing to just look at the alien community
native_inverse_plant_community.5 <- native_inverse_plant_community.5[native_inverse_plant_community.5$Origin == 'native',]
str(native_inverse_plant_community.5)

# Removes the origin column
native_inverse_plant_community.5 <- native_inverse_plant_community.5[-c(2) ]

# Change the shape of the community data frame
native_inverse_plant_community.5.t <- as.data.frame(t(native_inverse_plant_community.5[,-c(1)]))  
native_inverse_plant_community.5.t <- cbind(ID=row.names(native_inverse_plant_community.5.t),native_inverse_plant_community.5.t)
head(native_inverse_plant_community.5.t)

# Pointpattern from Bacaro et al 2012 method
#core_mxp <- pointpattern(core_plot_coords, longlat=TRUE)

# Spatial rarefaction from Bacaro et al 2012 method
SCR_curve_nat <- SCR(native_inverse_plant_community.5.t, core_mxp)

# the mxp matrix represents the output of the pointpattern function, 
# see Appendix 1 non_spatially_explicit<-specaccum(community_data[2:ncol(community_data)], method="exact")
# The following code can be used to generate Figure 3 in Bacaro et al. (2012).
# In this example the curves are not different from each other since the community
# data being used are synthetic (see community_data)
# Refer to Figure 3 in Bacaro et al. (2012) and to Palmer et al. (2007)
# for a published dataset on vascular inverse_plants from a North Carolina piedmont forest. 

plot(seq(1:50),SCR_curve_nat[,1], xlim=c(0,50), ylim=c(0,200), xlab="Number of plots", 
     ylab="Species Richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_nat[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_nat[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_nat[,3], lty=1, lwd=3, col="dark grey")



# Combining the three graphs
par(mfrow=c(1,1), family = "sans")
plot(seq(1:50),SCR_curve_all[,1], xlim=c(0,50), ylim=c(0,300), xlab="Number of sites", ylab="Number of species", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_all[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,3], lty=1, lwd=3, col="dark grey") 


lines(seq(1:50),SCR_curve_alien[,1], lwd=3, col="red") 
lines(seq(1:50),SCR_curve_alien[,2], lty=1 ,lwd=3, col="red") 
lines(seq(1:50),SCR_curve_alien[,3], lty=1, lwd=3, col="red")

lines(seq(1:50),SCR_curve_nat[,1], lwd=3, col="green") 
lines(seq(1:50),SCR_curve_nat[,2], lty=1 ,lwd=3, col="green") 
lines(seq(1:50),SCR_curve_nat[,3], lty=1, lwd=3, col="green")



#~# Combining graphs according to their community type 
par(mfrow=c(3,1), family = "sans",mai = c(0.7, 0.75, 0.15, 0.5))
#~# Entire inverse_plant community comparison
plot(seq(1:50),SCR_curve_all[,1], xlim=c(0,50), ylim=c(0,300), cex.axis = 1.5, cex.lab = 2,
     xlab="Number of plots", ylab="Species Richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_all[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_all[,3], lty=1, lwd=3, col="dark grey") 
plot(entire_inverse_plant_community.4_sac, ci.type="line", lwd= 3, add = TRUE, xlab ="Number of sites", ylab = "Number of species", cex.lab=1.5)

#~# Alien inverse_plant community comparison
plot(seq(1:50),SCR_curve_alien[,1], xlim=c(0,50), ylim=c(0,80), cex.axis = 1.5, cex.lab = 2,
     xlab="Number of plots", ylab="Species Richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_alien[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_alien[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_alien[,3], lty=1, lwd=3, col="dark grey") 
plot(alien_inverse_plant_community.4_sac, ci.type="line", lwd= 3, add = TRUE, xlab ="Number of sites", 
     ylab = "Number of species", cex.lab=1.5)

#~# Native inverse_plant community comparison
plot(seq(1:50),SCR_curve_nat[,1], xlim=c(0,50), ylim=c(0,200), cex.axis = 1.5, cex.lab = 2,
     xlab="Number of plots", ylab="Species Richness", type="p", pch=10,ann=TRUE,col="white") 
lines(seq(1:50),SCR_curve_nat[,1], lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_nat[,2], lty=1 ,lwd=3, col="dark grey") 
lines(seq(1:50),SCR_curve_nat[,3], lty=1, lwd=3, col="dark grey") 
plot(native_inverse_plant_community.4_sac, ci.type="line", lwd= 3, add = TRUE, xlab ="Number of sites", ylab = "Number of species", cex.lab=1.5)

# Save to (filename=("Outputs/Figures/Exploratory/To use/Combined standard and SC rarefaction curves.png"))












# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#     Percent occupancy graphs   
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Using the data file percent_occ
percent_occ

# Ordering the dataframe by the occupancy of alien
percent_occ_order <- percent_occ[order(percent_occ$alien.occupancy),]
percent_occ_order <- cbind(id=1, nrow(percent_occ_order), percent_occ_order)
percent_occ_order1 <- melt(percent_occ_order, id.vars=c("id", "Plot"),variable.name= "alien.occupancy",value.name="percent")




percent_occ_order1 <- melt(percent_occ, id.var="Plot")


# Plotting the graph
o <- percent_occ_order1 %>% filter(variable == "alien.occupancy") %>% arrange(value) %>% extract2("Plot")

percent_occ_order1 %>% 
  mutate(Plot = factor(Plot, o)) %>% 
  ggplot() +
  aes(x = Plot, y = value, fill = variable) +
  geom_bar(position = position_stack(reverse = TRUE), stat = "identity", width = .75, col = "black") + 
  xlab("Plot") +
  ylab("Percentage of species") +
  theme_bw(base_family = "sans") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size = 16),
        axis.text.y = element_text(vjust=0.5, hjust=1, size = 16),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "none") + 
  scale_fill_manual(values=c("alien.occupancy"="white", "native.occupancy"="grey")) +
  scale_y_continuous(breaks=c(seq(0,100,10)), limits = c(0,100), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

ggsave(filename=("Outputs/Figures/Exploratory/To use/Plant percent_occ.png"), width = 11, height = 9)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#    #~# The relationships between the core native and alien communities #~#                                          
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install.packages("ggiraphExtra")
library(ggiraphExtra )

# Using the data file alien_native_rel
alien_native_rel

# Checking the dataframe
head(alien_native_rel)
str(alien_native_rel) 

# Changing the abundance data to numeric not integers
str(alien_native_rel)
alien_native_rel[3:6] <- lapply(alien_native_rel[3:6], as.numeric)
alien_native_rel[1:2] <- lapply(alien_native_rel[1:2], as.factor)

install.packages("ggpmisc")
install.packages("broom")

my.formula <- alien.richness~native.richness

compare_means(alien.richness~native.richness + Park, data=alien_native_rel)

# Plotting native richness against alien richess 

fit1=lm(alien.richness~native.richness + Park, data=alien_native_rel)
summary(fit1)

native.rich.alien.rich <- ggplotRegression1(fit1) +
  scale_x_continuous("Native richness", breaks=c(seq(20,60,10)),
                     limits = c(20,60), expand = c(0, 0)) + 
  scale_y_continuous("Alien richness", breaks=c(seq(0,40,10)),
                     limits = c(0,40),expand = c(0, 0)) 

native.rich.alien.rich

native.rich.alien.rich + annotate("text", x = 1:2, y = 1:2, label = "ljihvgutfib")

native.rich.alien.rich <- ggplot(alien_native_rel, aes(x=native.richness, y=alien.richness, 
                                                       shape = Park)) + geom_point() +
                          geom_smooth(method="lm", se=FALSE, color="black") +
                              theme_bw(base_family = "sans") +
                              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), panel.border = element_blank(),
                                    plot.margin = unit(c(1,1,1,1), "cm"),
                                    axis.text = element_text(size = 10),
                                    axis.title = element_text(size = 12),
                                    axis.line = element_line(colour = "black", size = .5),
                                    legend.position = c(1, 1),
                                    legend.justification = c("right", "top")) +
                          scale_x_continuous("Native Richness", breaks=c(seq(20,60,10)),
                                             limits = c(20,60), expand = c(0, 0)) + 
                          scale_y_continuous("Alien Richness", breaks=c(seq(0,40,10)),
                                             limits = c(0,40),expand = c(0, 0)) 

native.rich.alien.rich

ggsave(filename=("Outputs/Figures/Exploratory/50 plots averaged HS7/alien and native richness relationships.png"), width = 16, height = 8)



native.rich.alien.rich <- ggplot(alien_native_rel, aes(x=native.richness, y=alien.richness, 
                                                       shape = Park)) + geom_point() +
  geom_smooth(method="lm", se=FALSE, color="black") +






fit <- lm(native.richness ~ alien.richness, data = alien_native_rel)
summary(fit)


# Plotting native occupancy against alien occupancy

fit <- lm(alien.occupancy ~ native.occupancy, data = alien_native_rel)
summary(fit)
Anova(fit)

native.occ.alien.occ <- ggplotRegression1(fit) +
                                scale_x_continuous(name = "Native frequency of occurrence", breaks=c(seq(100,700,100)),limits = c(100,700), expand = c(0, 0)) + 
                                scale_y_continuous(name = "Alien frequency of occurrence", breaks=c(seq(0,500,100)),limits = c(0,500),expand = c(0, 0)) 

native.occ.alien.occ

ggsave(filename=("Outputs/Figures/Exploratory/50 plots averaged HS7/alien and native occupancy relationship.png"), width = 16, height = 8)


# Plotting native richness against native occupancy 
fit <- lm(native.occupancy ~ native.richness, data = alien_native_rel)
summary(fit)

native.rich.native.occ <- ggplotRegression1(fit) +
                                  scale_x_continuous("Native richness", breaks=c(seq(20,60,10)),limits = c(20,60), expand = c(0, 0)) + 
                                  scale_y_continuous("Native frequency of occurrence", breaks=c(seq(0,800,200)),limits = c(0,800),expand = c(0, 0))
                                
ggsave(filename=("Outputs/Figures/Exploratory/50 plots averaged HS7/native richness against native occupancy.png"), width = 16, height = 8)


# NEW CODE FOR ANOVA (occupancy)
test <- ggplot(alien_native_occ, aes(x = Origin, y = Occupancy, colour = Origin))  + 
  geom_boxplot()
test

#One way anova
fit.1 <- aov(Occupancy ~ Origin, data = alien_native_occ)
summary(fit.1)

#Checking the data to see if an anova is ok
plot(fit.1, 1)
leveneTest(Occupancy ~ Origin, data = alien_native_occ) # it is

# NEW CODE FOR ANOVA (richness)
test <- ggplot(alien_native_ric, aes(x = Origin, y = Richness, colour = Origin))  + 
  geom_boxplot()
test

#One way anova
fit.2 <- aov(Richness ~ Origin, data = alien_native_ric)
summary(fit.2)

# Plotting alien richness against alien occupancy 
fit <- lm(alien.occupancy ~ alien.richness, data = alien_native_rel)
summary(fit)


alien.rich.alien.occ <- ggplotRegression1(fit) +
                                  scale_x_continuous("Alien richness", breaks=c(seq(0,40,10)),limits = c(0,40), expand = c(0, 0)) + 
                                  scale_y_continuous("Alien frequency of occurrence", breaks=c(seq(0,500,100)),limits = c(0,500),expand = c(0, 0))

ggsave(filename=("Outputs/Figures/Exploratory/50 plots averaged HS7/alien richness against alien occupancy.png"), width = 16, height = 8)


# Adding two columns for the entire community richness and occupancy
# Changing the abundance data to numeric not integers
str(alien_native_rel)
alien_native_rel[3:6] <- lapply(alien_native_rel[3:6], as.numeric)

#add column for entire richness
alien_native_rel$entire.richness <- alien_native_rel$alien.richness + alien_native_rel$native.richness

#add column for entire richness
alien_native_rel$entire.occupancy <- alien_native_rel$alien.occupancy + alien_native_rel$native.occupancy

alien_native_rel

# Plotting entire richness against entire occupancy 
fit <- lm(entire.occupancy ~ entire.richness, data = alien_native_rel)
summary(fit)

entire.rich.entire.occ <- ggplotRegression1(fit) +
                                  scale_x_continuous("Whole community richness", breaks=c(seq(30,80,10)),limits = c(30,80), expand = c(0, 0)) + 
                                  scale_y_continuous("Whole community frequency of occurrence", breaks=c(seq(200,1000,200)),limits = c(200,1000),expand = c(0, 0))

ggsave(filename=("Outputs/Figures/Exploratory/50 plots averaged HS7/entire richness against entire occupancy.png"), width = 16, height = 8)


# Figure with each richness comparion and occupancy comparison on the same plot
# Plotting native richness against alien richess 
fit <- lm(native.occupancy ~ native.richness, data = alien_native_rel)
summary(fit)

fit1 <- lm(alien.occupancy ~ alien.richness, data = alien_native_rel)
summary(fit1)

native.alien.richness.occ <- ggplot() + 
              geom_point(data=alien_native_rel, aes(x=native.richness, y=native.occupancy), size = 2) +
              geom_smooth(data=alien_native_rel, aes(x=native.richness, y=native.occupancy),
                          method="lm", se=F, colour = "black", size = 0.75) + 
              geom_point(data=alien_native_rel, aes(x=alien.richness, y=alien.occupancy), size = 2, colour = "dark grey") +
              geom_smooth(data=alien_native_rel, aes(x=alien.richness, y=alien.occupancy),
                          method="lm", se=F, colour = "dark grey", size = 0.75) + 
              theme_bw(base_family = "sans") +
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                    axis.text.x = element_text(vjust=0.5, size = 14),
                    axis.text.y = element_text(vjust=0.5, hjust=1, size = 14),
                    axis.title.x = element_text(size=16),
                    axis.title.y = element_text(size=16),
                    plot.margin = unit(c(1,1.5,0.5,0.5), "cm")) +
              scale_x_continuous("Richness", breaks=c(seq(0,60,10)),limits = c(0,60), expand = c(0, 0)) + 
              scale_y_continuous("Frequency of occurrence", breaks=c(seq(0,700,100)),limits = c(0,700),expand = c(0, 0))

native.alien.richness.occ

ggplotRegression2(fit, fit1) +
  scale_x_continuous("Richness", breaks=c(seq(0,60,10)),limits = c(0,60), expand = c(0, 0)) + 
  scale_y_continuous("Whole community frequency of occurrence", breaks=c(seq(0,700,100)),limits = c(0,700),expand = c(0, 0))

ggsave(filename=("Outputs/Figures/Exploratory/50 plots averaged HS7/alien and native comparison.png"), width = 16, height = 8)


combined.1 <- plot_grid(native.rich.alien.rich, native.occ.alien.occ, native.rich.native.occ,
                        alien.rich.alien.occ,entire.rich.entire.occ, native.alien.richness.occ, ncol=2,
                        align='hv', hjust = 0.01, labels = "auto", label_size = 16)
combined.1

ggsave(filename=("Outputs/Figures/Exploratory/50 plots averaged HS7/V3_Combined.png"), width = 11, height = 14, dpi = 500)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#    #~# GLMs run on overall, native and alien richness with environmental variables as predictors #~#          
#        - Also testing the variance inflation factor
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Using the following datasets: plantcom_env stand_envir

# Environmental variables to use (chosen through inspecting the Pearsons cross correlation coefficients 
# (see Environmental variables.csv), and considering which are important for teasing apart the PCA components)
#   First component
#     - OM, pH, LBA, 
#   Second component
#     - Altitude, Mg

# Now using the variables OM, pH, ECEC and C_N
# (08.06.20) Ru-running using OM, LBA, ECEC and pH

#~# Merge coordinates into plantcom data
plantcom_env <- merge(plot_coords,plantcom_env,by.x="ID",by.y="Site_id",all=TRUE)
str(plantcom_env)

# Plots
par(mfrow=c(1,1))

hist(plantcom_env$combined.richness)
hist(plantcom_env$OM)
hist(plantcom_env$LBA)
hist(plantcom_env$ECEC)
hist(plantcom_env$pH)

pairs(plantcom_env[,c("OM","LBA","ECEC","pH")], lower.panel = panel.smooth, upper.panel = panel.cor)


# Checking for spatial autocorrelation within the richness data

# Tools to measure SA in residuals
# Create near-neighbour object for moran's I
siteDat <- plantcom_env[,c("X","Y")]
coordinates(siteDat) <- ~X+Y
siteDat.nb<-knn2nb(knearneigh(siteDat))


# #~#~# Poisson regression - combined.richness - on the entire community
# #       (Count data that aren't necessarily bounded by an upper limit)
# 
 entire_plantcom_envvar.glm.1 <- glm(combined.richness ~ . , data=plantcom_env[,c("combined.richness", 
                                                                                  "OM", "LBA", "ECEC", 
                                                                                  "pH")], 
                                     family= poisson (link="log"))
 summary(entire_plantcom_envvar.glm.1)
# 
# # Psuedo R2: Calculated by 100 * (null - residual)/ null
# # = 100 * (64.664 - 62.354)/ 64.664
 100 * (64.664 - 62.354)/ 64.664
# 
 plot(entire_plantcom_envvar.glm.1)
# dev.off()
 hist(residuals(entire_plantcom_envvar.glm.1))
# # Model looks good, so just need to check for spatial correlation issues
# 
# #~# SA checks
 par(mfrow=c(2,2))
 SA.bubble.plot(entire_plantcom_envvar.glm.1,plantcom_env)
 lm.morantest(entire_plantcom_envvar.glm.1, nb2listw(siteDat.nb, style="W")) 
# # no significant spatial autocorrelation in regression residuals
# # * Doesn't look like much of a problem


#~#~# Poisson regression - alien.richness

alien_plantcom_envvar.glm <- glm(alien.richness ~ ., data=plantcom_env[,c("alien.richness", 
                                                                          "OM", "LBA", "ECEC", "pH")], 
                                 family= poisson (link="log"))

plot(alien_plantcom_envvar.glm)

summary(alien_plantcom_envvar.glm) 
Anova(alien_plantcom_envvar.glm)

# Psuedo R2: Calculated by 100 * (null - residual)/ null
# = 100 * ((179.764  - 85.262  )/ 179.764)
100 * (179.764 - 85.262)/ 179.764

par(mfrow=c(2,2))
plot(alien_plantcom_envvar.glm)
dev.off()
hist(residuals(alien_plantcom_envvar.glm))
# Model looks good, so just need to check for spatial correlation issues

str(alien_plantcom_envvar.glm)
#~# SA checks
SA.bubble.plot(alien_plantcom_envvar.glm,plantcom_env)
lm.morantest(alien_plantcom_envvar.glm, nb2listw(siteDat.nb, style="W")) 
# no significant spatial autocorrelation in regression residuals
# * Doesn't look like much of a problem


#~#~# Poisson regression - native.richness

native_plantcom_envvar.glm <- glm(native.richness ~ . , data=plantcom_env[,c("native.richness", 
                                                                             "OM", "LBA", "ECEC", "pH")], 
                                  family= poisson (link="log"))
summary(native_plantcom_envvar.glm) 

# Psuedo R2: Calculated by 100 * (null - residual)/ null
# = 100 * (90.681 - 65.040)/ 90.681
100 * (90.681 - 65.040)/ 90.681

par(mfrow=c(2,2))
plot(native_plantcom_envvar.glm)
dev.off()
hist(residuals(native_plantcom_envvar.glm))
# Model looks good, so just need to check for spatial correlation issues

#~#~# SA checks
SA.bubble.plot(native_plantcom_envvar.glm,plantcom_env)
lm.morantest(native_plantcom_envvar.glm, nb2listw(siteDat.nb, style="W")) # no significant spatial autocorrelation in regression residuals
# * Doesn't look like much of a problem


# Checking the Variance Inflation Factors
vif(alien_plantcom_envvar.glm)

# OM 2.248695 Slope 1.028073  ECEC 1.469692  C.N 1.862390 
# OM 4.493243  LBA 1.241582   ECEC 4.293055   pH 5.004422  


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#    #~# MVABUND: MANYGLMs run on overall, native and alien richness with environmental variables as predictors #~#          
#        - Also testing the variance inflation factor
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Environmental variables to use (chosen through inspecting the Pearsons cross correlation coefficients 
# (see Environmental variables.csv), and considering which are important for teasing apart the PCA components)
#   Using the variables OM, pH, ECEC and C_N


# Set so can see all residual plots on one screen
par(mfrow= c(2,2))

#~# Five datasets are used within this analysis: species abundance data and environmental variables data
# entire_occ, alien_occ, native_occ and raw_envir


head(entire_occ)
str(entire_occ)

head(alien_occ)
str(alien_occ)

head(native_occ)
str(native_occ)

head(raw_envir)
str(raw_envir)


# #~# Entire plant community #~# 
# entire_oc <- entire_occ
# 
# # Changing the abundance data to numeric not integers
# entire_oc[2:258] <- lapply(entire_occ[2:258], as.numeric)
# 
# 
# #~# Create a mvabund object
# abu <- mvabund(entire_oc[,-1])
# 
# #~# Having a look at the spread of the data
# boxplot(entire_oc[,2:11],horizontal = FALSE,las=2, main="Abundance") # Extremely zero inflated
# 
# #~# Having a quick look at the abundances
# plot(abu ~ entire_oc$Site.id, cex.axis = 0.8, cex = 0.8, transformation = "no")
# 
# #~# Fit a Poisson GLM
# entire_oc_env.mglm1 <- manyglm(abu ~ Organic.Matter + pH + ECEC + Carbon.Nitrogen.ratio, 
#                                data = raw_envir, family = 'poisson')
# 
# # Does the model fit the data?
# plot(entire_oc_env.mglm1)
# 
# #~# Fit a Negative Binomial GLM
# entire_oc_env.mglm2 <- manyglm(abu ~ Organic.Matter + pH + ECEC + Carbon.Nitrogen.ratio, 
#                                data = raw_envir, family = 'negative_binomial')
# 
# # Does the model fit the data?
# plot(entire_oc_env.mglm2)
# 
# 
# #~# Interpreting the results
# summary(entire_oc_env.mglm2)
# anova(entire_oc_env.mglm2)


#~# Alien plant community #~# 

# Changing the abundance data to numeric not integers
alien_occ[2:76] <- lapply(alien_occ[2:76], as.numeric)




foo <- mvformula(alien_abu~ raw.en)

best.r.sq( foo, n.xvars= 3)

#~# Create a mvabund object
alien_abu <- mvabund(alien_occ[,-1])

str(alien_abu)
raw.en <- as.matrix(raw_envir[2:20])
is.matrix(raw.en)
str(raw.en)

best.r.sq( alien_abu~ raw.en, n.xvars= 4)

raw.en[, 2]


#~# Having a look at the spread of the data
boxplot(alien_occ[,2:11],horizontal = FALSE,las=2, main="Abundance") # Extremely zero inflated

#~# Having a quick look at the abundances
plot(alien_abu ~ alien_occ$Site.id, cex.axis = 0.8, cex = 0.8, transformation = "no")

#~# Fit a Poisson GLM
alien_occ_env.mglm1 <- manyglm(alien_abu ~ Organic.Matter + Live.basal.area + ECEC + pH, 
                               data = raw_envir, family = 'poisson')

# Does the model fit the data?
plot(alien_occ_env.mglm1)

data(spider)

spiddat <- mvabund(spider$abund)
X <- spider$x
str(X)
best.r.sq( spiddat~X )


#~# Fit a Negative Binomial GLM
alien_occ_env.mglm2 <- manyglm(alien_abu ~ Organic.Matter + Live.basal.area + ECEC + pH, 
                               data = raw_envir, family = 'negative_binomial')

# Does the model fit the data?
plot(alien_occ_env.mglm2)

#~# Interpreting the results
summary(alien_occ_env.mglm2)
anova(alien_occ_env.mglm2)


#~# Native plant community #~# 

# Changing the abundance data to numeric not integers
native_occ[2:183] <- lapply(native_occ[2:183], as.numeric)


#~# Create a mvabund object
native_abu <- mvabund(native_occ[,-1])

#~# Having a look at the spread of the data
boxplot(native_occ[,2:11],horizontal = FALSE,las=2, main="Abundance") # Extremely zero inflated

#~# Having a quick look at the abundances
plot(native_abu ~ native_occ$Site.id, cex.axis = 0.8, cex = 0.8, transformation = "no")

#~# Fit a Poisson GLM
native_occ_env.mglm1 <- manyglm(native_abu ~ Organic.Matter + Live.basal.area + ECEC + pH, 
                                data = raw_envir, family = 'poisson')

# Does the model fit the data?
plot(native_occ_env.mglm1)


#~# Fit a Negative Binomial GLM
native_occ_env.mglm2 <- manyglm(native_abu ~ Organic.Matter + Live.basal.area + ECEC + pH, 
                                data = raw_envir, family = 'negative_binomial')

# Does the model fit the data?
plot(native_occ_env.mglm2)

#~# Interpreting the results
summary(native_occ_env.mglm2)
anova(native_occ_env.mglm2)


# ANOVAS run #

anova(entire_occ_env.mglm2)
anova(alien_occ_env.mglm2)
anova(native_occ_env.mglm2)








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Re-running with OM, Slope, C:N and LBA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~# Alien plant community #~# 

# Changing the abundance data to numeric not integers
alien_occ[2:76] <- lapply(alien_occ[2:76], as.numeric)




foo <- mvformula(alien_abu~ raw.en)

best.r.sq( foo, n.xvars= 3)

#~# Create a mvabund object
alien_abu <- mvabund(alien_occ[,-1])

str(alien_abu)
raw.en <- as.matrix(raw_envir[2:20])
is.matrix(raw.en)
str(raw.en)

best.r.sq( alien_abu~ raw.en, n.xvars= 4)

raw.en[, 2]


#~# Having a look at the spread of the data
boxplot(alien_occ[,2:11],horizontal = FALSE,las=2, main="Abundance") # Extremely zero inflated

#~# Having a quick look at the abundances
plot(alien_abu ~ alien_occ$Site.id, cex.axis = 0.8, cex = 0.8, transformation = "no")

#~# Fit a Negative Binomial GLM
alien_occ_env.mglm2 <- manyglm(alien_abu ~ Organic.Matter + Live.basal.area + ECEC + Carbon.Nitrogen.ratio, 
                               data = raw_envir, family = 'negative_binomial')

# Does the model fit the data?
plot(alien_occ_env.mglm2)

#~# Interpreting the results
summary(alien_occ_env.mglm2)
anova(alien_occ_env.mglm2)


#~# Native plant community #~# 

# Changing the abundance data to numeric not integers
native_occ[2:183] <- lapply(native_occ[2:183], as.numeric)


#~# Create a mvabund object
native_abu <- mvabund(native_occ[,-1])

#~# Having a look at the spread of the data
boxplot(native_occ[,2:11],horizontal = FALSE,las=2, main="Abundance") # Extremely zero inflated

#~# Having a quick look at the abundances
plot(native_abu ~ native_occ$Site.id, cex.axis = 0.8, cex = 0.8, transformation = "no")

#~# Fit a Poisson GLM
native_occ_env.mglm1 <- manyglm(native_abu ~ Organic.Matter + Live.basal.area + ECEC + Carbon.Nitrogen.ratio, 
                                data = raw_envir, family = 'poisson')

# Does the model fit the data?
plot(native_occ_env.mglm1)

#~# Fit a Negative Binomial GLM
native_occ_env.mglm2 <- manyglm(native_abu ~ Organic.Matter + Live.basal.area + ECEC + Carbon.Nitrogen.ratio, 
                                data = raw_envir, family = 'negative_binomial')

# Does the model fit the data?
plot(native_occ_env.mglm2)

#~# Interpreting the results
summary(native_occ_env.mglm2)
anova(native_occ_env.mglm2)


# ANOVAS run #

anova(entire_occ_env.mglm2)
anova(alien_occ_env.mglm2)
anova(native_occ_env.mglm2)