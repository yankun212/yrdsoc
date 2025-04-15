########Mixed effects model
rm(list = ls())

library(lmerTest)
library(readxl)
library(lme4)
library(vegan)
library(ade4)
library(glmm.hp)
library(lmerTest)
library(modelr)          
library(broom)
library(broom.mixed)          
library(ggsci)          
library(glmm.hp)

#Load data
data<-read.table("metadata.txt")

#Set the factor order and color
col<-c("#f2401c","#ffc000","#ffff01","#91d04f","#00af50","#01b0f1","#0071c0","#7030a0")
data$group<-factor(data$group,levels = c("S1","S2","S3","S4","S5","S6","S7","S8","S9"))
data$plant<-factor(data$plant,levels = c("Suaeda-glauca","Tamarix-chinensis","Salix-matsudana","Robinia-pseudoacacia","Fraxinus-chinensis","Populus"))

#Build model
f1<-lmer(SOC~Salt+(1|plant),data=data)
summary(f1)
coef(f1)
r.squaredGLMM(f1)

#visualization
a %>% add_predictions(f1) %>% ggplot(aes(x = Salt, y =SOC, group = data$plant,color =data$plant))+
  geom_point(size=3)+
  geom_line(aes(x = Salt, y = pred),size=1)+
  theme_bw()+theme(text = element_text(size = 15,face="bold",family = "serif"))+
  scale_color_manual(values = col)


#########Random Forest models
rm(list = ls())

library(rfPermute)
library(ggplot2)
library(A3)

#Load data
data<-read.table("metadata.txt")

#Build model
set.seed(123)
rfP <- rfPermute(MAOC~., data = data, importance = TRUE, ntree = 500, nrep = 1000, num.cores = 1)
rfP$rf

#Importance scores of extracted predictors (standardized scores)
importance_scale <- data.frame(importance(rfP, scale = TRUE), check.names = FALSE)
importance_scale

#Extract the significance of the importance score of the predictor (take the standardized score as an example)

importance_scale.pval <- (rfP$pval)[ , , 2]
importance_scale.pval


#Rank predictor variables by importance score, for example according to "%IncMSE"
importance_scale <- importance_scale[order(importance_scale$'%IncMSE', decreasing = TRUE), ]

#Simply plot the %IncMSE value of the predictor
importance_scale$name <- rownames(importance_scale)
importance_scale$name <- factor(importance_scale$name, levels = importance_scale$name)

#visualization
p <- ggplot() +
  geom_col(data = importance_scale, aes(x = name, y = `%IncMSE`), width = 0.5, fill = '#FFC068', color = NA) +
  labs(title = NULL, x = NULL, y = 'Increase in MSE (%)', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme_bw()+theme(text = element_text(size = 15,face="bold",family = "serif"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(expand = c(0, 0), limit = c(-2.5, 12))

p

# Label the significance information of the predictor
# By default, p<0.05 is *, p<0.01 is **, p<0.001 is ***
for (I in rownames(importance_scale)) {
  importance_scale[I,'%IncMSE.pval'] <- importance_scale.pval[I,'%IncMSE']
  if (importance_scale[I,'%IncMSE.pval'] >= 0.05) importance_scale[I,'%IncMSE.sig'] <- ''
  else if (importance_scale[I,'%IncMSE.pval'] >= 0.01 & importance_scale[I,'%IncMSE.pval'] < 0.05) importance_scale[I,'%IncMSE.sig'] <- '*'
  else if (importance_scale[I,'%IncMSE.pval'] >= 0.001 & importance_scale[I,'%IncMSE.pval'] < 0.01) importance_scale[I,'%IncMSE.sig'] <- '**'
  else if (importance_scale[I,'%IncMSE.pval'] < 0.001) importance_scale[I,'%IncMSE.sig'] <- '***'
}

p <- p +
  geom_text(data = importance_scale, aes(x = name, y = `%IncMSE`, label = `%IncMSE.sig`), nudge_y = 0.5)

p

# Top right note the known interpretation rate of the model
p <- p +
  annotate('text', label = 'MAOC', x = 9, y = 15, size = 4) +
  annotate('text', label = sprintf('italic(R^2) == %.2f', 41.82), x = 9, y = 8, size = 3, parse = TRUE)

p

#The significance of the entire model is calculated
#model.fn=randomForest Call the method of random forest to calculate
#p.acc=0.001 Represents an estimate of the p-value based on 1000 random permutations
set.seed(123)
otu_forest.pval <- a3(MAOC~., data = data, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 50))
otu_forest.pval

##########PCA
rm(list = ls())

library(factoextra)

pca <- prcomp(data, scale. = TRUE)

fviz_pca_biplot(pca2,geom = c("point"),col.var = "black",
                col.ind  = a$Salt,palette="npg",pointsize = 3, 
                labelsize = 4,gradient.cols = c("#01b0f1","#ffff01","#FC4E07"))+
  theme_bw()+theme(text = element_text(size = 15,face="bold",family = "serif"))


##########SEM
rm(list = ls())
library(lme4)
library(piecewiseSEM)

##MNC

mod1 <- psem(
  lmer(TN~Salt+mineral+(1|plant),data = data),
  lmer(MNC~interaction+growth+qCO2+TN+Salt+(1|plant),data = data),
  lmer(growth~Salt+TN+interaction+(1|plant),data = data),
  lmer(qCO2~Salt+TN+interaction+(1|plant),data = data),
  lmer(interaction~Salt+TN+mineral+(1|plant),data = data),
  qCO2 %~~% growth,
  data = data
)
summary(mod1)


##MNC/SOC
mod <- psem(
  lmer(TN~Salt+mineral+(1|plant),data = data),
  lmer(MNC.SOC~interaction+growth+qCO2+TN+Salt+(1|plant),data = data),
  lmer(growth~Salt+TN+interaction+(1|plant),data = data),
  lmer(qCO2~Salt+TN+interaction+(1|plant),data = data),
  lmer(interaction~Salt+TN+mineral+(1|plant),data = data),
  qCO2 %~~% growth,
  data = data
)
summary(mod)


##MAOC

modh <- psem(
  lmer(TN~Salt+mineral+(1|plant),data = data),
  lmer(MAOC~Salt+TN+MNC+mineral+(1|plant),data = data),
  lmer(MNC~interaction+growth+qCO2+TN+Salt+(1|plant),data = data),
  lmer(growth~Salt+TN+interaction+(1|plant),data = data),
  lmer(qCO2~Salt+TN+interaction+(1|plant),data = data),
  lmer(interaction~Salt+TN+mineral+(1|plant),data = data),
  qCO2 %~~% growth,
  data = data
)
summary(modh)


###MAOC/SOC

modh <- psem(
  lmer(TN~Salt+mineral+(1|plant),data = data),
  lmer(MAOC.SOC~Salt+TN+MNC+mineral+(1|plant),data = data),
  lmer(MNC~interaction+growth+qCO2+TN+Salt+(1|plant),data = data),
  lmer(growth~Salt+TN+interaction+(1|plant),data = data),
  lmer(qCO2~Salt+TN+interaction+(1|plant),data = data),
  lmer(interaction~Salt+TN+mineral+(1|plant),data = data),
  qCO2 %~~% growth,
  data = data
)

summary(modh)


