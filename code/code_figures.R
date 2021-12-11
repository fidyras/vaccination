library(tidyverse)
library(scales)
library(patchwork)
library(sf)
library(ggmap)




mdg2<-read_sf("data/shapefile/MDG_ADM1.shp")
POPDATA <- read.csv("data/age_dist.csv")
VACDATA <- read.csv("data/vacc_alloc2.csv")
REG <- unique(POPDATA$key)
TOTALPOP <- sum(POPDATA$sum_pop)
TOTALVAC <-  ceiling(c(TOTALPOP*0.2)) 

## map of vaccine allocation by strategy

POPDATA%>%group_by(Region)%>%summarise(population=sum(sum_pop))->a

mdg2%>%
  left_join(a,by=c("NAME"="Region"))%>%
  right_join(VACDATA,by=c("NAME"="Region"))%>%
  mutate("vo1"=ceiling(TOTALVAC*freq_pop),"vo2"=ceiling(TOTALVAC*freq_60),"vo3"=ceiling(TOTALVAC*freqcases),"vo4"=ceiling(TOTALVAC*freqdeaths))->mdg2
mdg2<-rmapshaper::ms_simplify(mdg2, keep = 0.15, keep_shapes = T)
mdg2<-st_as_sf(mdg2)


ggplot(mdg2)+
  geom_sf(aes(fill=freq_pop))+
  scale_fill_viridis_b(option="B",breaks=c(seq(0.025,0.125,0.025)),label=percent,name="")+
  #guides(fill=guide_colorbar(nbin=12))+
  #labs(title="Pro-rata")+
  theme_void()->vo1

ggplot(mdg2)+
  geom_sf(aes(fill=freq_60))+
  scale_fill_viridis_b(option="B",breaks=c(seq(0.025,0.125,0.025)),label=percent,name="")+
  #guides(fill=guide_colorbar(nbin=12))+
  #labs(title="Distribution of elderly")+
  theme_void()->vo2

my_breaks=c(0.01,0.03,0.05,0.1,0.3,0.5,1)

ggplot(mdg2)+
  geom_sf(aes(fill=freqcases))+
  scale_fill_viridis_b(name="",breaks=my_breaks,trans="log2",labels=percent,option="B")+
  #labs(title="Distribution of cases")+
  theme_void()->vo3

ggplot(mdg2)+
  geom_sf(aes(fill=freqdeaths))+
  scale_fill_viridis_b(name="",breaks=my_breaks,trans="log2",labels=percent,option="B")+
  #guides(fill=guide_colorbar(nbin=12))+
  #labs(title="Distribution of cases")+
  theme_void()->vo4


############
#### FIGURE 2
vo1|vo2|vo3|vo4|plot_annotation(tag_levels = "A")




#####################
############# Baseline
#### reduction of deaths as a function of strategy

basedata<-read.csv("Simdata/basedata2_test.csv")%>%
  select(-1)

basedata_r2%>%mutate("PropR_dist"="uniform")%>%
  mutate("PropR"=spG)%>%filter(spG==0)%>%
  group_by(vo,acceptV,totalvac,Simul,effV,startV,PropR,PropR_dist,testing,onlyS,vacByAge,r0)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:12,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->baseline


baseline$percent_red<-ifelse(baseline$percent_red<0,0,baseline$percent_red)
baseline$acceptance<-as.numeric(baseline$acceptance)
baseline$Simul<-as.factor(baseline$Simul)
baseline$vo<-as.factor(baseline$vo)
baseline$r0<-as.factor(baseline$r0)
levels(baseline$vo) <- c("Population", "age", "cases", "deaths", "uniform")



base_fig<-baseline%>%
  filter(acceptance==0.7)%>%
  filter(onlyS==FALSE & testing==FALSE)%>%
  ggplot(aes(x=vo,y=percent_red,group=vo,fill=vo))+
  geom_violin(trim=FALSE)+
  #geom_boxplot(width=0.1)+
  #stat_summary(fun=median, geom="point", size=2)+
  theme_classic(base_size = 14)+
  scale_y_continuous(name="Reduction in mortality",label = percent)+
  scale_x_discrete(name="",breaks=NULL)+
  coord_cartesian(ylim=c(0.25,0.5))+
  facet_wrap(~r0)+
  scale_fill_brewer(palette= "Dark2",name="Allocation strategy",label=c("Pro-rata", "Age", "Cases", "Deaths", "Uniform"))

base_fig+plot_layout(guides = "collect")+plot_annotation(tag_levels = "A")
##########

################### SIMULATIONS
###########   SIMULATION 1
#### reduction of deaths as a function of doses available 
# vary total vaccine supply with acceptance rate of 70% and 50% of HCW doing 20 vaccines/day

Simdata1<-read.csv("Simdata/simulation1.csv")%>%select(-1)

Simdata1$Simul<-as.factor(Simdata1$Simul)
Simdata1%>%
  group_by(vo,acceptV,totalvac,Simul,effV,startV)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:6,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->Simulation1

Simulation1$acceptance<-as.numeric(Simulation1$acceptance)
Simulation1$vo<-as.factor(Simulation1$vo)
levels(Simulation1$vo) <- c("Population", "Age", "Cases", "Deaths", "Uniform")
Simulation1$percent_red<-ifelse(Simulation1$percent_red<0,0,Simulation1$percent_red)

Simulation1%>%filter(acceptance==0.70)%>%group_by(nvax,vo)%>%summarise(median=median(percent_red,na.rm=TRUE))->median_sim1


Simulation1%>%
  filter(acceptance==0.7)%>%
  ggplot()+
  geom_line(aes(x=nvax,y=percent_red,group=interaction(vo,Simul),color=vo),size=0.1,alpha=0.2)+
  geom_smooth(data=median_sim1,aes(x=nvax,y=median,group=vo,color=vo),se=FALSE,size=1)+
  theme_classic(base_size = 12)+
  scale_x_continuous(name="Total vaccine supply",
                     label=percent)+
  scale_y_continuous(name="Reduction in mortality",label = percent,
                     limits=c(0,1))+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))->sim1_fig

sim1_fig


####### Simulation 2
#####vary vaccine acceptance assuming enough doses to cover 20% of the population and 50% of HCW performing 20 vaccines/day
############## vaccinating only seronegative individuals (testing prior to vaccination, at vaccination site)

Simdata2<-read.csv("Simdata/simulation2_data.csv")%>%select(-1)
Simdata2$acceptV<-as.factor(Simdata2$acceptV)
Simdata2%>%
  group_by(vo,acceptV,totalvac,Simul,effV,startV)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:6,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->Simulation2

Simulation2$acceptance<-as.numeric(Simulation2$acceptance)

Simulation2$vo<-as.factor(Simulation2$vo)
levels(Simulation2$vo) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")


Simulation2$percent_red<-ifelse(Simulation2$percent_red<0,0,Simulation2$percent_red)
Simulation2%>%
  #filter(nvax==0.2)%>%
  group_by(acceptance,vo)%>%
  summarise(median=mean(percent_red))->median_sim2

Simulation2%>% 
  #filter(nvax==0.2)%>% #filter baseline number of vaccines available
  ggplot()+
  geom_smooth(aes(x=acceptance,y=percent_red,group=interaction(Simul,vo),color=vo),
              size=0.05,se=F)+
  geom_smooth(data=median_sim2,aes(x=acceptance, y=median, group=interaction(vo),color=vo),
              size=1,se=F)+
  theme_classic(base_size = 12)+
  scale_x_continuous(name="Vaccine acceptance",
                     label=percent)+
  scale_y_continuous(name="Reduction in mortality",label = percent)+
  scale_color_brewer(palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"),guide=NULL)->sim2_fig



sim2_fig






#######simulation 3

#### reduction of deaths as a function of the rollout speed 
#vary number of HCW doing 20 vaccines/day assuming enough vaccines for 20% of the population and 70% acceptance

simData3<-read.csv("Simdata/Simulation3.csv")%>%select(-1)
simData3$vo<-as.factor(simData3$vo)
levels(simData3$vo) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")

simData3%>%
  group_by(vo,acceptV,vpd,Simul)%>%
  summarise(deaths=sum(Number))%>%
  mutate(staff=vpd/20)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:4,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->Simulation3

Simulation3%>%filter(acceptance==0.7)%>%
  ggplot()+
  geom_smooth(aes(x=staff,y=percent_red,group=vo,color=vo),position=position_dodge(width=0.1),size=1.5,se=FALSE)+
  theme_classic(base_size = 12)+
  geom_line(aes(x=staff,y=percent_red,group=interaction(Simul,vo),color=vo),
            size=0.1,se=F,alpha=0.5)+
  scale_x_continuous(name="Speed of vaccine rollout",
                     label=percent)+
  scale_y_continuous(name="Reduction in mortality",label = percent,limits=c(0,0.5))+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))->sim3_fig



################### SIMULATION 5
###########   vary start date of vaccination

Simdata5<-read.csv("Simdata/Simulation5.csv")
Simdata5$Simul<-as.factor(Simdata5$Simul)
Simdata5%>%
  group_by(vo,acceptV,totalvac,Simul,effV,vpd,startV)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:6,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->Simulation5

Simulation5$vo<-as.factor(Simulation5$vo)
levels(Simulation5$vo) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")

Simulation5%>%filter(acceptance==0.70)%>%group_by(startV,vo)%>%summarise(median=median(percent_red))->median_sim5
Simulation5%>%filter(acceptance==0.7)%>%
  ggplot()+
  geom_line(aes(x=startV,y=percent_red,group=interaction(vo,Simul),color=vo),size=0.1,alpha=0.5)+
  geom_line(data=median_sim5,aes(x=startV,y=median,group=vo,color=vo),se=FALSE,size=1)+
  theme_classic(base_size = 12)+
  scale_x_continuous(name="Start of vaccination (days)")+
  scale_y_continuous(name="Reduction in mortality",label = percent,
                     limits=c(0,0.5))+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))->sim5_fig




######
############# FIGURE 3

(base_fig|(sim1_fig+sim2_fig)/(sim3_fig+sim5_fig))+
  plot_annotation(tag_levels = "A")+plot_layout(guides="collect")& theme(legend.position = 'bottom')




###########################################################################
####################### vaccinating only seronegative individuals 

###########################################################################
##### testing prior to vaccination, at vaccination site: ie no replacement)

basedata<-read.csv("Simdata/basedata2_test.csv")%>%
  select(-1)

basedata%>%
  filter(onlyS==FALSE)%>%
  mutate("PropR"=spG)%>%
  group_by(vo,acceptV,totalvac,Simul,effV,startV,PropR,PropR_dist,testing,vacByAge,onlyS)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:11,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->baseline


baseline$percent_red<-ifelse(baseline$percent_red<0,0,baseline$percent_red)
baseline$acceptance<-as.numeric(baseline$acceptance)
baseline$Simul<-as.factor(baseline$Simul)
baseline$vo<-as.factor(baseline$vo)
levels(baseline$vo) <- c("Population", "Age", "Cases", "Deaths", "Uniform")

baseline%>%filter(acceptance==0.70)%>%group_by(testing,vo,PropR,onlyS)%>%summarise(median=median(percent_red,na.rm=TRUE))->median_test
baseline%>%
  filter(acceptance==0.7)%>%
  ggplot(aes(x=PropR,y=percent_red,group=vo,fill=vo))+
  geom_line(data=median_test,aes(x=PropR,y=median,group=interaction(vo,testing,onlyS),color=vo,lty=testing),se=FALSE,size=1)+
  geom_line(aes(group=interaction(vo,Simul,testing),color=vo,lty=testing),alpha=0.05)+
  theme_classic(base_size = 14)+ theme(legend.position = "none")+
  scale_y_continuous(name="Reduction in mortality",label = percent)+
  scale_x_continuous(name="Underlying seroprevalence",label=percent)+
  coord_cartesian(ylim=c(0.1,0.8),xlim=c(0,0.4))+
  scale_linetype_manual(values=c(1,3),name="Ab-testing prior to vaccination (on site)")+
  scale_color_brewer(palette= "Dark2",name="Allocation strategy",label=c("Pro-rata", "Age", "Cases", "Deaths", "Uniform"))->base_test_fig


##############################################################################
###### testing prior to vaccination, before vaccination site: ie only Susceptibles come to vccination)


basedata%>%
  filter(testing==FALSE)%>%
  mutate("PropR"=spG)%>%
  group_by(vo,acceptV,totalvac,Simul,effV,startV,PropR,PropR_dist,testing,vacByAge,onlyS)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:11,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->baseline


baseline$percent_red<-ifelse(baseline$percent_red<0,0,baseline$percent_red)
baseline$acceptance<-as.numeric(baseline$acceptance)
baseline$Simul<-as.factor(baseline$Simul)
baseline$vo<-as.factor(baseline$vo)
levels(baseline$vo) <- c("Population", "age", "cases", "deaths", "uniform")

baseline%>%filter(acceptance==0.70)%>%group_by(testing,vo,PropR,onlyS)%>%summarise(median=median(percent_red,na.rm=TRUE))->median_test
baseline%>%
  filter(acceptance==0.7)%>%
  ggplot(aes(x=PropR,y=percent_red,group=vo,fill=vo))+
  geom_line(data=median_test,aes(x=PropR,y=median,group=interaction(vo,testing,onlyS),color=vo,lty=onlyS),se=FALSE,size=1)+
  geom_line(aes(group=interaction(vo,Simul,testing,onlyS),color=vo,lty=onlyS),alpha=0.05)+
  theme_classic(base_size = 14) + theme(legend.position = "none")+
  scale_y_continuous(name="Reduction in mortality",label = percent)+
  scale_x_continuous(name="Underlying seroprevalence",label=percent)+
  coord_cartesian(ylim=c(0.1,0.8),xlim=c(0,0.4))+
  scale_linetype_manual(values=c(1,3),name="Only Susceptibles come to vaccination")+
  scale_color_brewer(palette= "Dark2",name="Allocation strategy",label=c("Pro-rata", "Age", "Cases", "Deaths", "Uniform"))->onlyS_fig




######################################################

############################ distribution of seropositives


basedata_spG<-read.csv("Simdata/basedata_spG.csv")%>%select(-c(1,2))


basedata_spG%>%
  group_by(vo,acceptV,totalvac,Simul,effV,startV,PropR,PropR_dist,testing,vacByAge,onlyS)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:11,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->baseline

baseline$acceptance<-as.numeric(baseline$acceptance)
baseline$PropR_dist<-as.factor(baseline$PropR_dist)
baseline$vo<-as.factor(baseline$vo)
levels(baseline$vo)<-c("pro-rata", "age", "cases", "deaths", "uniform")

baseline%>%filter(acceptance==0.70)%>%filter(PropR_dist!="deaths")%>%
  group_by(PropR_dist,vo,PropR)%>%
  summarise(median=median(percent_red,na.rm=TRUE))->median_PropR

baseline%>%
  filter(acceptance==0.7)%>%
  filter(PropR_dist!="deaths")%>%
  ggplot()+
  geom_line(aes(x=PropR,y=percent_red,color=vo,group=interaction(Simul,vo,PropR_dist),lty=PropR_dist),alpha=0.1)+
  geom_smooth(data=median_PropR,aes(x=PropR,y=median,color=vo,lty=PropR_dist,group=interaction(vo,PropR_dist)),size=1.5,se=F)+
  theme_classic(base_size = 14)+
  scale_linetype_manual(name="Distribution of seropositives",values=c(2,1),labels=c("cases","uniform"),guide=NULL)+
  scale_y_continuous(name="Reduction of mortality",labels=percent)+
  scale_x_continuous(name="Underlying seroprevalence",labels=percent)+
  coord_cartesian(xlim=c(0,0.39))+
  scale_color_brewer(palette= "Dark2",name="allocation strategy")->distr_sero

#############
#### FIGURE 4

distr_sero/(onlyS_fig+base_test_fig)+plot_annotation(tag_levels = "A")


