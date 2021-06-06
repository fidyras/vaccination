library(tidyverse)
library(scales)
library(patchwork)

#####################
############# Baseline
#### reduction of deaths as a function of strategy and facet by underlying seroprevalence

basedata<-read.csv("Simdata/baseline.csv")%>%select(-1)%>%mutate("PropR"=0.0,"PropR_dist"="Uniform")
basedata_PropR<-rbind(basedata,basedata0.1,basedata0.2)
write.csv(basedata_PropR,"Simdata/basedata_PropR.csv")
summaryD%>%select(-1)%>%mutate("PropR"=0.0,"PropR_dist"="Uniform")%>%
  group_by(vo,acceptV,totalvac,Simul,effV,startV,PropR,PropR_dist)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:8,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->baseline


baseline$percent_red<-ifelse(baseline$percent_red<0,0,baseline$percent_red)
baseline$acceptance<-as.numeric(baseline$acceptance)
baseline$Simul<-as.factor(baseline$Simul)
baseline$vo<-as.factor(baseline$vo)
levels(baseline$vo) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")

baseline%>%filter(acceptance==0.7)%>%
  ggplot(aes(x=vo,y=percent_red,group=vo,fill=vo))+
  geom_boxplot()+
  geom_boxplot()+
  #stat_summary(fun=median, geom="point", size=2)+
  theme_classic(base_size = 14)+
  scale_y_continuous(name="Reduction in mortality",label = percent)+
  scale_x_discrete(name="Prioritization strategy",breaks=NULL)+
  scale_fill_brewer(palette= "Dark2",name="Allocation strategy",label=c("Pro-rata", "Distribution of >60", "Distribution of Cases", "Distribution of deaths", "Uniform"))+
  labs(title="Reduction in Mortality based on allocation strategy",
       subtitle="Baseline conditions 0-40% seroprevalence")+facet_wrap(~PropR)

baseline%>%group_by(vo,acceptance)%>%summarise(median(percent_red))



###################
###########   SIMULATION 1
#### reduction of deaths as a function of doses available 
# vary total vaccine supply with acceptance rate of 70% and 50% of HCW doing 20 vaccines/day

Simdata1<-read.csv("Simdata/simulation1.csv")%>%select(-1)
#Simdata1<-rbind(Simdata1,summaryD)
#write.csv(Simdata1,"Simdata/simulation1.csv")
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
levels(Simulation1$vo) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")
Simulation1$percent_red<-ifelse(Simulation1$percent_red<0,0,Simulation1$percent_red)

Simulation1%>%filter(acceptance==0.70)%>%group_by(nvax,vo)%>%summarise(median=median(percent_red,na.rm=TRUE))->median_sim1


Simulation1%>%filter(acceptance==0.7)%>%
  ggplot()+
  geom_line(aes(x=nvax,y=percent_red,group=interaction(vo,Simul),color=vo),size=0.1,alpha=0.2)+
  geom_smooth(data=median_sim1,aes(x=nvax,y=median,group=vo,color=vo),se=FALSE,size=1)+
  theme_classic(base_size = 10)+
  scale_x_continuous(name="Total vaccine supply (% of population)",
                     label=percent)+
  scale_y_continuous(name="Reduction in mortality",label = percent,
                     limits=c(0,1))+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))+
  labs(
       subtitle="Total supply")->sim1_fig

sim1_fig

####### Simulation 2
#####vary vaccine acceptance assuming enough doses to cover 20% of the population and 50% of HCW performing 20 vaccines/day
Simdata2<-read.csv("Simdata/simulation2.csv")%>%select(-1)
# Simdata2<-rbind(Simdata2,summaryD)
Simdata2$acceptV<-as.factor(Simdata2$acceptV)
Simdata2%>%
  group_by(vo,acceptV,totalvac,Simul)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:4,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->Simulation2

Simulation2$acceptance<-as.numeric(Simulation2$acceptance)

Simulation2$vo<-as.factor(Simulation2$vo)
levels(Simulation2$vo) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")

#Simulation2$age<-as.factor(Simulation2$age)
#levels(Simulation2$age)<- c("0-9","10-19","20-29","30-39","40-49","50-59","60+")

Simulation2$percent_red<-ifelse(Simulation2$percent_red<0,0,Simulation2$percent_red)
Simulation2%>%
  #filter(nvax==0.2)%>%
  group_by(acceptance,vo)%>%summarise(median=median(percent_red))->median_sim2

Simulation2%>% 
  #filter(nvax==0.2)%>% #filter baseline number of vaccines available
  ggplot()+
  geom_line(aes(x=acceptance,y=percent_red,group=interaction(Simul,vo),color=vo),
            size=0.1,se=F,alpha=0.2)+
  geom_smooth(data=median_sim2,aes(x=acceptance,y=median,group=vo,color=vo),
            size=1,se=F)+
  theme_classic(base_size = 10)+
  scale_x_continuous(name="Vaccine acceptance (% of population)",
                     label=percent)+
  scale_y_continuous(name="Reduction in mortality",label = percent,limits=c(0,0.5))+
  scale_color_brewer(palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"),guide=NULL)+
  labs(subtitle="Acceptance rate ")->sim2_fig



sim2_fig

###############################################

#######simulation 3

#### reduction of deaths as a function of the rollout speed 
#vary number of HCW doing 20 vaccines/day assuming enough vaccines for 20% of the population and 70% acceptance

simData3<-read.csv("Simdata/Simulation3.csv")
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
  theme_classic(base_size = 10)+
  geom_line(aes(x=staff,y=percent_red,group=interaction(Simul,vo),color=vo),
            size=0.1,se=F,alpha=0.5)+
  scale_x_continuous(name="Speed of vaccine rollout",
                     label=percent)+
  scale_y_continuous(name="Reduction in mortality",label = percent,limits=c(0,0.5))+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))+
  labs(subtitle="Rollout speed")->sim3_fig

####### Simulation 4
#### vary vaccine efficacy

Simdata4<-read.csv("Simdata/Simulation4.csv")%>%select(-1)
Simdata4%>%
  group_by(vo,acceptV,totalvac,Simul,effV,startV)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:6,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->Simulation4

Simulation4$vo<-as.factor(Simulation4$vo)
levels(Simulation4$vo) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")
Simulation4%>%#filter(acceptance==0.70&nvax==0.2)%>%
  group_by(effV,vo)%>%
  summarise(median=median(percent_red))->median_sim4

Simulation4$acceptance<-as.numeric(Simulation4$acceptance)
Simulation4%>%#filter(acceptance==0.70&nvax==0.2)%>%
  ggplot()+
  geom_line(aes(x=effV,y=percent_red,group=interaction(Simul),color=vo),
            alpha=0.1,size=0.5)+
  geom_line(data=median_sim4,aes(x=effV,y=median,group=vo,color=vo),
            size=1,se=F)+
  theme_classic(base_size = 10)+
  scale_x_continuous(name="Vaccine efficacy (%)",
                     label=percent)+
  scale_y_continuous(name="Reduction in mortality",label = percent,limits = c(0,0.8))+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))+
  labs(subtitle="Vaccine efficacy")->sim4_fig


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
  geom_smooth(data=median_sim5,aes(x=startV,y=median,group=vo,color=vo),se=FALSE,size=1)+
  theme_classic(base_size = 10)+
  scale_x_continuous(name="Start of vaccination (days)")+
  scale_y_continuous(name="Reduction in mortality",label = percent,
                     limits=c(0,0.5))+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))+
  labs(subtitle="Start day")->sim5_fig

######
library(patchwork)
(base_fig|(sim1_fig+sim2_fig)/(sim3_fig+sim5_fig))+
  plot_annotation(tag_levels = "A")+plot_layout(guides="collect")


###### visualize reduction as function of rollout speed and acceptance or total supply and speed of rollout

visual<-read.csv("Simdata/vpd_accept.csv")
visual%>%
  group_by(vo,acceptV,totalvac,Simul,vpd,effV,startV)%>%
  summarise(deaths=sum(Number),vPerD=sum(vPerD))%>%
  mutate(staff=vpd/20)%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:9,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->visdata

visdata$vo<-as.factor(visdata$vo)
levels(visdata$vo) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")
str(visdata)
visdata$acceptance<-as.numeric(visdata$acceptance)
visdata%>%#filter(nvax==0.2)%>%
  ggplot()+
  geom_tile(aes(x=vPerD,y=acceptance,fill=percent_red))+
  theme_classic(base_size = 14)+
  scale_x_continuous(name="rollout speed (number of vaccines per day)")+
  scale_y_continuous(name="Total Supply (% population)",label = percent)+
  scale_fill_viridis_c(option="A",label=percent, name="Reduction in Mortality")+
  facet_wrap(~vo)





## map of vaccine allocation by strategy
library(sf)
library(ggmap)
mdg2<-read_sf("data/shapefile/MDG_ADM1.shp")
POPDATA%>%group_by(Region)%>%summarise(population=sum(sum_pop))->a
mdg2%>%
  left_join(a,by=c("NAME"="Region"))%>%
  right_join(VACDATA,by=c("NAME"="Region"))%>%
  mutate("vo1"=ceiling(totalvac*freq_pop),"vo2"=ceiling(totalvac*freq_60),"vo3"=ceiling(totalvac*freqcases),"vo4"=ceiling(totalvac*freqdeaths))->mdg2
mdg2<-rmapshaper::ms_simplify(mdg2, keep = 0.15, keep_shapes = T)
mdg2<-st_as_sf(mdg2)
ggplot(mdg2)+
  geom_sf(aes(fill=freq_pop))+
  scale_fill_viridis_c(option="B",breaks=c(seq(0.025,0.125,0.025)),label=c("2.5%","5%","7.5%","10%","12.5%"),name="Proportion of doses \n to allocate")+
  #guides(fill=guide_colorbar(nbin=12))+
  #labs(title="Pro-rata")+
  theme_void()->vo1

ggplot(mdg2)+
  geom_sf(aes(fill=freq_60))+
  scale_fill_viridis_c(option="B",breaks=c(seq(0,0.20,0.025)),guide="none")+
  #guides(fill=guide_colorbar(nbin=12))+
  #labs(title="Distribution of elderly")+
  theme_void()->vo2

ggplot(mdg2)+
  geom_sf(aes(fill=freqcases))+
  scale_fill_viridis_b(option="D",breaks=c(seq(0,0.7,0.05)))+
  #guides(fill=guide_colorbar(nbin=12))+
  #labs(title="Distribution of cases")+
  theme_void()->vo3

ggplot(mdg2)+
  geom_sf(aes(fill=freqdeaths))+
  scale_fill_viridis_b(option="E",breaks=c(seq(0,0.70,0.05)))+
  #guides(fill=guide_colorbar(nbin=12))+
  #labs(title="Distribution of cases")+
  theme_void()->vo4

vo1+vo2+vo3+vo4+plot_annotation(tag_levels = "A")+plot_layout(guides="collect")




###### SEIR timeseries



ts<-read.csv("Simdata/SEIR_accept_0.7_start_10_.csv")  

ts$allocV_ts<-as.factor(ts$allocV_ts)
levels(ts$allocV_ts) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")
ts$age<-as.factor(ts$age)
levels(ts$age)<- c("0-9","10-19","20-29","30-39","40-49","50-59","60+")
options(scipen=999)
ts$acceptV<-as.factor(ts$acceptV)
library(scales)

ts%>%
  group_by(time,status,allocV_ts,regID_ts,acceptV,Simul)%>%
  summarise(Number=sum(Number))%>%ungroup()%>%
  pivot_wider(names_from = status, values_from=Number)%>%
  mutate("R"=U+R,"I"=A+I)%>%
  select(time,S,E,A,I,R,D,nVac,allocV_ts,regID_ts,acceptV,Simul)%>%
  pivot_longer(cols=-c(time,allocV_ts,regID_ts,acceptV,Simul),names_to="status", values_to="Number")%>%
  filter(status=="R"|status=="S"|status=="I"|status=="D"|status=="E"|status=="nVac")%>%
  filter(allocV_ts=="Population")->SEIR_ts

SEIR_ts$acceptV<-as.factor(SEIR_ts$acceptV)
SEIR_ts%>%
  filter(status=="nVac")%>%
  ggplot(aes(x=time,y=Number,group=interaction(allocV_ts,status,regID_ts,acceptV,Simul),color=regID_ts,lty=acceptV))+
  geom_line(size=0.3)+
  scale_y_continuous(labels= comma)+
  #scale_color_brewer(palette="Dark2",name="Status")+
  theme_bw()+
  geom_vline(aes(xintercept=10),lty=2, color= "grey30")+
  labs(title = "Vaccination per day")

library(scales)


ts%>%#filter(regID_ts=="AN")%>%
  group_by(time,status,allocV_ts,acceptV,Simul,age)%>%
  summarise(Number=sum(Number))%>%
  pivot_wider(names_from = status, values_from=Number)%>%
  mutate("R"=U+R)%>%
  select(time,S,E,A,I,R,D,nVac,allocV_ts,acceptV,Simul,age)%>%
  pivot_longer(cols=-c(time,allocV_ts,acceptV,Simul,age),names_to="status", values_to="Number")%>%
  # filter(status=="R")%>%
  ggplot(aes(x=time,y=Number,group=interaction(allocV_ts,status),color=status))+
  geom_line(size=0.75,alpha=0.5)+
  # scale_linetype_discrete(name="Allocation strategy")+
  scale_y_continuous(labels= comma)+
  scale_color_brewer(palette="Dark2",name="age category")+
  theme_bw()+
  labs(title = "vaccination accross age categories through time",subtitle = "baseline conditions")+
  facet_wrap(~acceptV)
library(tidyverse)
ts%>%filter(Simul==1)%>%
  group_by(status,allocV_ts,acceptV,Simul,time)%>%
  #filter(regID_ts=="AN")%>%
  summarise(Number=sum(Number))%>%
  pivot_wider(names_from = status, values_from=Number)%>%
  mutate("R"=U+R)%>%
  select(S,E,A,I,R,D,nVac,allocV_ts,acceptV,Simul,time)%>%
  pivot_longer(cols=-c(allocV_ts,acceptV,time,Simul),names_to="status", values_to="Number")%>%
  #filter(status=="nVac"|status=="D"|status=="R")%>%
  filter(allocV_ts=="1")->c

c%>%
  # filter(status=="R"|status=="S")%>%
  ggplot()+
  geom_line(aes(x=time,y=Number,color=status,group=interaction(status,Simul)))+
  scale_y_continuous(label=comma)+
  facet_wrap(~acceptV)
c%>%filter(time==0)
vpd+vpd_age

ts_age_novax+coord_cartesian(ylim=c(0,500))->ts_age_novax
ts_age_vax+coord_cartesian(ylim=c(0,500))->ts_age_vax
ts_age_novax+ts_age_vax+plot_layout(guides="collect")

timeseries_D_novax+coord_cartesian(ylim=c(0,1500))->timeseries_D_novax
timeseries_D_vax+coord_cartesian(ylim=c(0,1500))->timeseries_D_vax
(timeseries_novax+timeseries_vax)/(timeseries_D_novax+timeseries_D_vax)+plot_layout(guides="collect")

