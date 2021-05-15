library(tidyverse)
library(scales)
library(patchwork)
#####################
############# Baseline
#### reduction of deaths as a function of strategy
baseline<-read.csv("Simdata/baseline.csv")


baseline%>%
  group_by(vo,acceptV,totalvac,Simul,effV,startV)%>%
  summarise(deaths=sum(Number))%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:6,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->basedata


basedata$percent_red<-ifelse(basedata$percent_red<0,0,basedata$percent_red)
basedata$acceptance<-as.numeric(basedata$acceptance)
basedata$Simul<-as.factor(basedata$Simul)
basedata$vo<-as.factor(basedata$vo)
levels(basedata$vo) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")
basedata%>%filter(acceptance==0.7&nvax==0.2)%>%
  ggplot()+
  geom_violin(aes(x=vo,y=percent_red,group=vo,fill=vo))+
  theme_classic(base_size = 14)+
  scale_y_continuous(name="Reduction in mortality",label = percent)+
  scale_x_discrete(name="Prioritization strategy",breaks=NULL)+
  scale_fill_brewer(palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))+
  labs(title="Reduction in Mortality based on allocation strategy",
       subtitle="Baseline conditions")->base_fig





###################
###########   SIMULATION 1
#### reduction of deaths as a function of doses available 
# vary total vaccine supply with acceptance rate of 70% and 50% of HCW doing 20 vaccines/day

Simdata1<-read.csv("Simdata/simulation1.csv")
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

Simulation1%>%filter(acceptance==0.70)%>%group_by(nvax,vo)%>%summarise(median=median(percent_red))->median_sim1


Simulation1%>%filter(acceptance==0.7)%>%
  ggplot()+
  geom_line(aes(x=nvax,y=percent_red,group=interaction(vo,Simul),color=vo),size=0.1,alpha=0.5)+
  geom_smooth(data=median_sim1,aes(x=nvax,y=median,group=vo,color=vo),se=FALSE,size=1)+
  theme_classic(base_size = 14)+
  scale_x_continuous(name="Total vaccine supply (% of population)",
                     label=percent)+
  scale_y_continuous(name="Reduction in mortality",label = percent,
                     limits=c(0,1))+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))+
  labs(
       subtitle="Variable supply")->sim1_fig

####### Simulation 2
#####vary vaccine acceptance assuming enough doses to cover 20% of the population and 50% of HCW performing 20 vaccines/day
Simdata2<-read.csv("Simdata/simulation2.csv")
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
Simulation2$percent_red<-ifelse(Simulation2$percent_red<0,0,Simulation2$percent_red)
Simulation2%>%filter(nvax==0.2)%>%group_by(acceptance,vo)%>%summarise(median=median(percent_red))->median_sim2

Simulation2%>% filter(nvax==0.2)%>% #filter baseline number of vaccines available
  ggplot()+
  geom_line(aes(x=acceptance,y=percent_red,group=interaction(Simul,vo),color=vo),
            size=0.1,se=F,alpha=0.5)+
  geom_line(data=median_sim2,aes(x=acceptance,y=median,group=vo,color=vo),
            size=1,se=F)+
  theme_classic(base_size = 14)+
  scale_x_continuous(name="Vaccine acceptance (% of population)",
                     label=percent)+
  scale_y_continuous(name="Reduction in mortality",label = percent)+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))+
  labs(
       subtitle="Variable acceptance rate ")->sim2_fig



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
  theme_classic(base_size = 14)+
  geom_line(aes(x=staff,y=percent_red,group=interaction(Simul,vo),color=vo),
            size=0.1,se=F,alpha=0.5)+
  scale_x_continuous(name="Speed of vaccine rollout (% of healthcare workers)",
                     label=percent)+
  scale_y_continuous(name="Reduction in mortality",label = percent)+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))+
  labs(subtitle="Rollout speed")->sim3_fig

####### Simulation 4
#### vary vaccine efficacy

Simdata4<-read.csv("Simdata/Simulation4.csv")
Simdata4$Simul<-as.factor(Simdata4$Simul)
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
Simulation4%>%filter(acceptance==0.70&nvax==0.2)%>%
  group_by(effV,vo)%>%
  summarise(median=median(percent_red))->median_sim4


Simulation4%>%filter(acceptance==0.70&nvax==0.2)%>%
  ggplot()+
  geom_line(aes(x=effV,y=percent_red,group=interaction(Simul,vo),color=vo),
            alpha=0.1,size=0.5)+
  geom_line(data=median_sim4,aes(x=effV,y=median,group=vo,color=vo),
            size=1,se=F)+
  theme_classic(base_size = 14)+
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
  theme_classic(base_size = 14)+
  scale_x_continuous(name="Start of vaccination after begining of epidemic (days)")+
  scale_y_continuous(name="Reduction in mortality",label = percent,
                     limits=c(0,1))+
  scale_color_brewer(guide="none",palette= "Dark2",name="Allocation strategy",label=c("Population", "Proportion of Older people", "Proportion of Cases", "Proportion of deaths", "Uniform"))+
  labs(subtitle="Vaccination start day")->sim5_fig

######
library(patchwork)
(sim1_fig+sim2_fig+sim3_fig)/(sim4_fig+sim5_fig)+
  plot_annotation(tag_levels = "A")


###### visualize reduction as function of total supply and acceptance or total supply and speed of rollout

visual<-read.csv("Simdata/visual.csv")
visual%>%
  group_by(vo,acceptV,totalvac,Simul,vpd,effV,startV)%>%
  summarise(deaths=sum(Number))%>%
  mutate(staff=vpd/20)%>%
  mutate(nvax=totalvac/TOTALPOP)%>%
  pivot_wider(names_from=acceptV,values_from=deaths)%>%
  mutate(novax=`0`)%>%
  pivot_longer(cols=-c(1:8,novax),names_to="acceptance",values_to="deaths")%>%
  ungroup()%>%
  mutate(percent_red=1-(deaths/novax))->visdata

visdata$vo<-as.factor(visdata$vo)
levels(visdata$vo) <- c("Population", "Elderly", "Number of cases", "Number of deaths", "Uniform allocation")

visdata%>%filter(acceptance==0.7)%>%
  ggplot()+
  geom_tile(aes(x=staff,y=nvax,fill=percent_red))+
  theme_classic(base_size = 14)+
  scale_x_continuous(name="rollout speed (% of HCW)",
                     label=percent)+
  scale_y_continuous(name="Total Supply (% population)",label = percent)+
  scale_fill_viridis_c(option="A",label=percent, name="Reduction in Mortality")+
  facet_wrap(~vo)
