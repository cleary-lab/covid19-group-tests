library(tidyverse)
waning_dat <- read_csv("data/waning.csv")
waning_dat <- waning_dat %>% group_by(Patients) %>% mutate(time1 = time - min(time,na.rm=TRUE)) %>% drop_na()
greb <- waning_dat %>% group_by(time) %>% summarise(n_positive=sum(positive),
                                                    n=n(),
                                                    prop_positive=n_positive/n,
                                                    lower_confint=prop.test(n_positive, n)$conf.int[1],
                                                    upper_confint=prop.test(n_positive, n)$conf.int[2]
                                                    
                                                    )
greb <- greb %>% drop_na()
ggplot(greb) +
  geom_ribbon(aes(x=time,ymin=lower_confint,ymax=upper_confint),alpha=0.25) +
  geom_line(aes(x=time,y=prop_positive))

waning_dat1 <- waning_dat %>% 
  drop_na() %>% 
  ungroup() %>% 
  select(Patients, positive, time) %>% 
  unique() %>% 
  arrange(time,Patients)

waning_dat1$Patients <- as.factor(waning_dat1$Patients)
waning_dat1$time1 <- paste0("day",waning_dat1$time1)
waning_dat1$positive <- as.factor(waning_dat1$positive)
waning_dat1 %>%
  pivot_wider(names_from=time1, values_from=positive) %>% View

final <- matrix(nrow=nrow(greb),ncol=2)
for(i in 1:nrow(greb)){
  n_positive <- greb$n_positive[i]
  n <- greb$n[i]
  res <- prop.test(n_positive,n)
  final[i,1] <- res$conf.int[1]
  final[i,2] <- res$conf.int[2]
  
}
colnames(final) <- c("lower","upper")
greb <- cbind(greb, final)


ggplot(greb) + 
  #geom_ribbon(aes(x=time,ymin=lower,ymax=upper),fill="blue",alpha=0.5) +
  geom_jitter(data=waning_dat, aes(x=time1,y=positive),height=0.01,width=0.01) +
  geom_line(aes(x=time1,y=prop_positive)) 

model <- glm(data=waning_dat, positive~time, family="binomial")

predict(model)



ggplot() + 
  geom_line(data=waning_dat,aes(x=time1,y=positive,group=Patients)) + facet_wrap(~Patients)

