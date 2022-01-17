#Hyloxalus
#
##### Packages #####
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Call packages
packages <- c("ggplot2", #figs
              "dplyr","tidyr","Rmisc",
              "geiger","phytools","ggtree") #phylogenetic signal
ipak(packages)

#Data####
data<-read.csv("acoustics.csv")
head(data)
data$RR<-log(data$Value.Hv/data$Mean.Value.Hspp)
mFun<-subset(data,Variable=="Mean Fundamental frequency (Hz)")
meanFun<-mean(mFun$RR, na.rm = T)
mNote<-subset(data,Variable=="Mean Note duration (s)")
meanNote<-mean(mNote$RR, na.rm = T)
mPeak<-subset(data,Variable=="Mean Dominant frequency (Hz)")
meanPeak<-mean(mPeak$RR, na.rm = T)

#Graficar ####

#B1: Overall effect size RR of peak frequency###
mPeak<-subset(mPeak,RR!="NA")
R<-1000
B1<-rep(0,R);
for(kk in 1:R){
##Sample 1 data per study
Overall<-ddply(mPeak,.(Species))
  boot.sample <- sample(mPeak$RR, replace = TRUE)
  B1[kk] <- mean(boot.sample)
  print(B1)
}
boxplot(B1)
quantile(B1,c(0.025,0.975))
hist(B1, breaks = 30)
dB1<-c(Bmean=mean(B1),quantile(B1,c(0.025,0.975)),n=nrow(mPeak))
dB1

#B2: Overall effect size RR of Fund Peak###
mFun<-subset(mFun,RR!="NA")
R<-1000
B2<-rep(0,R);
for(kk in 1:R){
##Sample 1 data per study
Overall<-ddply(mFun,.(Species))
  boot.sample <- sample(mFun$RR, replace = TRUE)
  B2[kk] <- mean(boot.sample)
  print(B2)
}
boxplot(B2)
quantile(B2,c(0.025,0.975))
hist(B2, breaks = 30)
dB2<-c(Bmean=mean(B2),quantile(B2,c(0.025,0.975)),n=nrow(mFun))
dB2

#B3: Overall effect size RR of Note duration###
mNote<-subset(mNote,RR!="NA")
R<-1000
B3<-rep(0,R);
for(kk in 1:R){
##Sample 1 data per study
Overall<-ddply(mNote,.(Species))
  boot.sample <- sample(mNote$RR, replace = TRUE)
  B3[kk] <- mean(boot.sample)
  print(B3)
}
boxplot(B3)
quantile(B3,c(0.025,0.975))
hist(B3, breaks = 30)
dB3<-c(Bmean=mean(B3),quantile(B3,c(0.025,0.975)),n=nrow(mNote))
dB3

#(podria sacar los overall para los tres subgrupos: bocagei, subpunctatus y other)

#Join results
#Resume data of Bootstrap (Bootstrap mean, confidence limits[2.5%-97.5%],n)
dbR=data.frame(dB1,dB2,dB3)
t.dbR<-t(dbR)
datboots<-as.data.frame(t.dbR)
head(datboots)
datboots$Variable<-c("Mean Dominant frequency (Hz)",
                     "Mean Fundamental frequency (Hz)",
                     "Mean Note duration (s)")
head(datboots)
colnames(datboots)<-c("RR","Lower","Upper","n","Variable")
datboots$Spp<-c("OVERALL","OVERALL","OVERALL")

meanGraph<-data.frame(Variable = c('Mean Dominant frequency (Hz)',
                                   'Mean Fundamental frequency (Hz)',
                                   'Mean Note duration (s)'), 
                      meanValues = c(meanPeak,meanFun,meanNote))

meanGraph$Variable = factor(meanGraph$Variable, levels=c('Mean Dominant frequency (Hz)',
                                             'Mean Fundamental frequency (Hz)',
                                             'Mean Note duration (s)'))
data$Variable = factor(data$Variable, levels=c('Mean Dominant frequency (Hz)',
                                             'Mean Fundamental frequency (Hz)',
                                             'Mean Note duration (s)'))
datboots$Variable = factor(datboots$Variable, levels=c('Mean Dominant frequency (Hz)',
                                             'Mean Fundamental frequency (Hz)',
                                             'Mean Note duration (s)'))
data$Spp<- factor(data$Spp,levels = c("vertebralis '","toachi '","sanctamariensis '",
                                      "pulchellus '",
                                      "patitae '","nexipus '","fuliginosus \"", 
                                      "fascianigrus \"",
                                      "elachyhistus '","azureiventris '","awa '","arliensis '",
                                      "abditaurantius '", "subpunctatus '","felixcoperari '",
                                      "cepedai '","yasuni '","sauli '","maculosus \"",
                                      "italoi '",
                                      "bocagei \"",
                                      "OVERALL"))
my_x_title <- expression(paste("Species of ", italic("Hyloxalus")))

ggplot(data, aes(x=Spp, y=RR)) + 
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_hline(aes(yintercept=meanValues),meanGraph, linetype = "dashed",colour="gray")+
  geom_point(alpha=0.5) +
  geom_pointrange(aes(ymin=Lower,ymax=Upper),datboots, size=0.7)+
  labs(title="", x=my_x_title, y = "Response ratio")+  #Para anadir titulo y cambiar los nombres de los ejes
  coord_flip()+
  theme_bw()+ #Para eliminar el color gris del fondo
  theme(panel.grid.major = element_blank(), #elimina las lineas grandes dentro del plot
        panel.grid.minor = element_blank(), #elimina las lineas pequenas dentro del plot
        axis.text.y = element_text(face="italic"), #hace que los nombres de las variables esten en italico y en angulo
        legend.position = "none",
        strip.background = element_blank())+  #Posicion de la leyenda
  facet_grid(.~Variable, scales = "free_x")
#Guardar esta figura en buena resolución .tiff a 300 dpi 
ggsave("RR-Hyloxalus.tiff", units="in", width=8, height=4.5, dpi=300, compression = 'lzw')



#Phylogenetic tree (Not this time) ####
an_trees=read.nexus("files2_consensus.nex")

# Generating the consensus tree (when the position agree in >%50 of times) - just 1 time
#an_tree=consensus.edges(an_trees, p = 0.5,check.labels = FALSE)
# Exporting data
#write.tree(an_tree,file="Hyloxalus_tree_majorR.tre")#save the consensus tree
# Reading data
an_tree <- read.tree("Hyloxalus_tree_majorR.tre")
# Ploting the tree
plot(an_tree,edge.width=2) #type="fan"
spp<-as.data.frame(table(data$Species,data$Orden))
spp<-subset(spp,Freq!=0)
rownames(spp)=spp$Var1  
spp<-spp[order(match(spp$Var1, an_tree$tip.label)), ]

tZ<-ggtree(an_tree, layout='fan',open.angle=180,ladderize = F,size=.5) %<+% spp+
  geom_tiplab2(font="italic", align=F,offset=200,aes(color=Var2),hjust=0.5,size=1)+
  scale_color_manual(name="Order",
                     values = c("Anura" = "#80d043", "Urodela"="#2271d8"))+
  theme(legend.position = "none")  #Posicion de la leyenda

tz<-rotate_tree(tZ,270)
tz

p_g_ageC<-summarySE(data=data,measurevar="g_massC",groupvars=c("Species"),na.rm = TRUE)
p_g_ageC
rownames(p_g_ageC)<-p_g_ageC$Species
#Phylogenetic signal?
phylosig(an_tree, p_g_ageC$g_massC, method="K", 
         test=TRUE, nsim=1000, 
         se=NULL, start=NULL,control=list())
phylosig(an_tree, p_g_ageC$g_massC, method = "lambda" , 
         test=TRUE, nsim=1000, 
         se=NULL, start=NULL,control=list())

