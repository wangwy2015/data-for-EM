library(vegan) 
library(spaa)
library (psych)
library(ggplot2)
setwd("D:/703/wwy/0C-李飞数据/分长期人工草地/生态位宽度指数")

dune<-read.csv('OTU.csv', row.names = 1, check.names = FALSE)
niche_width <- niche.width(dune, method = 'levins')
niche_width
boxplot(unlist(niche_width),xlab="", ylab = 'niche breadth')
barplot(unlist(niche_width),xlab="", ylab = 'niche breadth')
my_niche_width<-t(niche_width) 
write.csv(my_niche_width,"生态位宽度指数_bacteria.csv")


dune<-read.csv('OTU_fungi.csv', row.names = 1, check.names = FALSE)
niche_width <- niche.width(dune, method = 'levins')
niche_width
boxplot(unlist(niche_width),xlab="", ylab = 'niche breadth')
barplot(unlist(niche_width),xlab="", ylab = 'niche breadth')
my_niche_width<-t(niche_width) #先对数据进行转置
write.csv(my_niche_width,"生态位宽度指数_fungi.csv")


dune<-read.csv('OTU_筛选18.csv', row.names = 1, check.names = FALSE)
niche_width <- niche.width(dune, method = 'levins')
niche_width
boxplot(unlist(niche_width),xlab="", ylab = 'niche breadth')
barplot(unlist(niche_width),xlab="", ylab = 'niche breadth')
my_niche_width<-t(niche_width) 
write.csv(my_niche_width,"生态位宽度指数_bacteria_筛选18.csv")


setwd("D:/703/wwy/0C-李飞数据/分长期人工草地/生态位宽度指数")
data<-read.csv("生态位宽度指数统计分析.csv",row.names = 1,header = T)
library(dplyr)
library(ggplot2)
data_bac<-data[1:5,]
data_bac$Type<- factor(data_bac$Type, levels=c('ASM','AKM','SDM', 'ASR', 'ALR'))
p_bac<-ggplot(data=data_bac,aes(x=Type,y= mean,fill = Type))+
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.25, size = 0.3, position = position_dodge(0.8)) +
  geom_text(aes(label=sig,y= mean + se ),size=6,
            position = position_dodge(width = 0.8), 
            vjust=-0.3)+
  labs(x='',color="")+
  ylab("Nich breadth")+
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = rep(c('#0072B2','#F0E442','#009E73','#56B4E9','#CC79A7'),3))
p_bac  

data_fun<-data[6:10,]
data_fun$Type<- factor(data_fun$Type, levels=c('ASM','AKM','SDM', 'ASR', 'ALR'))
p_fun<-ggplot(data=data_fun,aes(x=Type,y= mean, fill = Type))+
  geom_bar(stat = 'identity', position = 'dodge')+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.25, size = 0.3, position = position_dodge(0.8)) +
  geom_text(aes(label=sig,y= mean + se),size=6,
            position = position_dodge(width = 0.8), 
            vjust=-0.3)+
  labs(x='',color="")+
  ylab("Nich breadth")+
  theme_bw(base_size = 18)+  
  theme(axis.text = element_text(colour = 'black'),
        legend.position = "")+
  scale_fill_manual(values = rep(c('#0072B2','#F0E442','#009E73','#56B4E9','#CC79A7'),3))
p_fun  


data$Type<- factor(data$Type, levels=c('ASM','AKM','SDM', 'ASR', 'ALR'))
p <- ggplot(data, aes(Type, mean)) +
  geom_bar(stat = 'identity', position = 'dodge', fill="grey")+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.25, size = 0.3, position = position_dodge(0.8)) +
  geom_text(aes(label=sig,y= mean + se),size=6,
            position = position_dodge(width = 0.8), 
            vjust=-0.3)+
  labs(x='',color="")+
  ylab("Nich breadth")+
  theme_bw(base_size = 18)+
  theme(panel.grid=element_blank())+
  theme(axis.text = element_text(colour = 'black'),
        legend.position = "")+
  #scale_fill_manual(values = rep(c('#0072B2','#F0E442','#009E73','#56B4E9','#CC79A7'),3))+
  facet_wrap(~Experiment, scales = "free_y")
p

ggsave('生态位宽度指数.pdf', p, width = 8, height = 5)
ggsave('生态位宽度指数.png', p, width = 8, height = 5)






####随机森林模型预测
library(randomForest)
library(caret)
library(tidyverse)
setwd("D:/703/wwy/0C-李飞数据/分长期人工草地/随机森林模型")
data <- read.csv("EMF随机森林预测.csv")
data=na.omit(data)  
# 拆分数据集为训练集和测试集
set.seed(123) 
train_index <- sample(1:nrow(data), 0.7 * nrow(data))  # 70%的数据作为训练集
train_data <- data[train_index, ]
test_data <- data[-train_index, ]
# 构建随机森林模型
model <- randomForest(EMF ~ ., data = train_data, ntree = 100)
# 查看模型概览
print(model)
# 在测试集上进行预测
predictions <- predict(model, test_data)
# 评估模型性能
accuracy <- sum(predictions == test_data$EMF) / length(test_data$EMF)

accuracy <- sum(predictions) / length(test_data$EMF)

print(paste("模型准确率：", accuracy))


importance_ozo <- model$importance
importance_ozo
importance_plot <- tibble(var = rownames(importance_ozo), 
                          IncMSE = importance_ozo[,1],
                          IncNodePurity = importance_ozo[,2])


library(randomForest)
library(caret)
library(tidyverse)
library(ranger)
library(rfPermute)
library(ggsci)
setwd("D:/703/wwy/0C-李飞数据/分长期人工草地/随机森林模型")
data <- read.csv("EMF随机森林预测.csv")
RFdata <- scale(data)
RFdata <- na.omit(RFdata)
set.seed(100)
richness_rf <- randomForest(EMF ~ ., data= RFdata,
                            importance=TRUE,proximity=TRUE)
richness_rf

#set.seed(100)
#richness_perm <- rf.significance(richness_rf, RFdata[,-1], nperm=99, ntree=501) 


set.seed(100)
richness_rfP<- rfPermute(EMF ~ ., data = RFdata, ntree = 500,
                         na.action = na.omit, nrep = 100,num.cores = 1)
richness_dat <- importance(richness_rfP, sort.by = NULL, decreasing = TRUE)

Fig.1 = richness_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*",""))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>%
  ggplot(aes(x = names, y = X.IncMSE))+
  geom_bar(aes(fill = label),stat = "identity")+
  scale_fill_npg()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),legend.title =element_blank() )+
  theme(legend.position = 'none') +
  geom_text(aes(y = X.IncMSE + 1,label = label),size=9)+
  labs(x = "Variables", y = "Increase in MSE (%)")+
  annotate("text", x = 1.5, y = 8, label = 'EMF' , colour = "black",size=5) +
  geom_label( x = 4.5, y = 6.5, label = expression('R'^2==.19~~'p<0.05' ), colour = "black",label.size=0,size=5) +
  coord_flip()
Fig.1


data <- read.csv("EMF随机森林预测_3.csv")
RFdata <- scale(data)
RFdata <- na.omit(RFdata)
set.seed(100)
richness_rf <- randomForest(EMF ~ ., data= RFdata,
                            importance=TRUE,proximity=TRUE)
richness_rf

#set.seed(100)
#richness_perm <- rf.significance(richness_rf, RFdata[,-1], nperm=99, ntree=501) 


set.seed(100)
richness_rfP<- rfPermute(EMF ~ ., data = RFdata, ntree = 500,
                         na.action = na.omit, nrep = 100,num.cores = 1)
richness_dat <- importance(richness_rfP, sort.by = NULL, decreasing = TRUE)
write.csv(richness_dat,"随机森林模型预测指标贡献.csv")

Fig.1 = richness_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*",""))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>%
  ggplot(aes(y = names, x = X.IncMSE))+
  geom_bar(aes(fill = label),stat = "identity")+
  scale_fill_npg()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),legend.title =element_blank() )+
  theme(axis.text.x = element_text(color='black',size=15),
        axis.text.y = element_text(color='black',size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.text= element_text(size=15))+
  theme(legend.position = 'none') +
  geom_text(aes(x = X.IncMSE + 1,label = label),size=9)+
  labs(y = "", x = "Increase in MSE (%)")+
  #annotate("text", x = 2, y = 4.5, label = 'EMF' , colour = "black",size=5) +
  geom_label( x = 2, y = 8, label = expression('R'^2==.70~~', p<0.05' ), colour = "black",label.size=0,size=5) +
  coord_flip()
Fig.1

ggsave('随机森林模型指标贡献.png', Fig.1, width = 5, height = 4)
ggsave('随机森林模型指标贡献.pdf', Fig.1, width = 5, height = 4)
