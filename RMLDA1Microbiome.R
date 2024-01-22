setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/Thesis_NVP/Chapter2")
library(rstan)
library(cmdstanr)
DataComb<- read.csv("~/Documents/Documents - Namitha’s MacBook Pro/Thesis_NVP/Chapter1/MicrobiomeCombined.csv")
summary(DataComb)
DataComb1<-DataComb[-which(DataComb$SubjectID %in% c("69-001")),]
DataComb<-DataComb1
Data1=DataComb[,10:118]
N=max(rowSums(Data1))
V=ncol(Data1)
#Data Pre-processing
Data1$SubjectID=as.factor(DataComb$SubjectID)
n=length(unique(Data1$SubjectID))
M=length(unique(Data1$SubjectID))
Data=data.frame(matrix(ncol=(ncol(Data1)-1), nrow=n))
Data_n=data.frame(matrix(ncol=(ncol(Data1)-1), nrow=n))
num_sub<-c()
num_sub_n<-c()
for(i in 1:n)
{
  Data[i,1:(ncol(Data1)-1)]=floor(colMeans(Data1[which(Data1$SubjectID==unique(Data1$SubjectID)[i]),1:109]))
  num_sub[i]<-nrow(Data1[which(Data1$SubjectID==unique(Data1$SubjectID)[i]),1:109])
  Data_n[i,1:(ncol(Data1)-1)]=floor(colMeans(Data1[which(Data1$SubjectID==unique(Data1$SubjectID)[i]& DataComb$CL4=="Healthy"),1:109]))
  num_sub_n[i]<-nrow(Data1[which(Data1$SubjectID==unique(Data1$SubjectID)[i] & DataComb$CL4=="Healthy"),1:109])
}

subid<-unique(Data1$Subject)

DataComb1<-DataComb[which(DataComb$SubjectID%in% subid[which(num_sub>=5 & num_sub<=25)]),]
DataComb<-DataComb1
Data1<-DataComb
Data1=DataComb[,10:118]
Data1$SubjectID=as.factor(DataComb$SubjectID)
DataComb$SubjectID<-DataComb1$SubjectID
subid<-unique(Data1$SubjectID)
n=length(unique(Data1$SubjectID))
num_sub<-c()
for(i in 1:n)
{
  num_sub[i]<-nrow(Data1[which(Data1$SubjectID==unique(Data1$SubjectID)[i]),1:109])
}
n_d<-max(num_sub)
M<-length(unique(Data1$SubjectID))
#J
V
M
n_d
N

Data<-Data1[,-110]

pred_name=colnames(Data1)[-110]

summary(as.factor(DataComb$CL4))
DataComb$CL4_regroup=c()
DataComb$CL4_regroup[which(DataComb$CL4=="Healthy")]="Healthy"
DataComb$CL4_regroup[which(DataComb$CL4=="Infection_L")]="Infection"
DataComb$CL4_regroup[which(DataComb$CL4=="Infection/Fiber")]="Infection"
DataComb$CL4_regroup[which(is.na(DataComb$CL4_regroup)==TRUE)]="Stress"
lev1<-c("Healthy","Infection","Stress")
DataComb$CL4_regroup<-factor(DataComb$CL4_regroup,levels=lev1)
summary(DataComb$CL4_regroup)

DataComb$CL4_regroup=DataComb$CL4
DataComb$CL4_regroup[which(DataComb$CL4=="Infection_L")]="Infection"
DataComb$CL4_regroup[which(DataComb$CL4=="Ant_L")]="Ant"
DataComb$CL4_regroup[which(DataComb$CL4=="Colonoscopy_L")]="Colonoscopy"
DataComb$CL4_regroup[which(DataComb$CL4=="Imz_L")]="Imz"
DataComb$CL4_regroup[which(DataComb$CL4=="Post-Travel_L")]="Post-Travel"

summary(as.factor(DataComb$CL4_regroup))/nrow(DataComb)

########
library(tm)
Data1$SubjectID=as.factor(DataComb1$SubjectID)

V<-ncol(Data1)-1
w<-array(,dim = c(M, n_d, V))
for( d in 1:M)
{
  w[d,1:length(which(Data1$SubjectID==subid[d])),]<-as.matrix(as.DocumentTermMatrix(Data1[which(Data1$SubjectID==subid[d]),-c(ncol(Data1))], 
                                                                                    weighting = weightTf))
}


for( d in 1:M)
{
  for( r in 1:n_d)
  {
    for(v in 1:V)
    {
      if(is.na(w[d,r,v]==TRUE))
      {
        w[d,r,v]=0
      }
    }
  }
}


summary(as.factor(DataComb$CL4_regroup))
DataComb$CL4_regroup<-as.factor(DataComb$CL4_regroup)


###RM-LDA #####
tab<-list()
acc1<-c()
for(J in 2:15)
  #{
J=11
datar <- list(J=J, V=V, M=M, n_d=n_d,w=w,mu=rep(0,J),num_sub=num_sub)
datar
model<- cmdstan_model(stan_file = "RM_LDA_dtm.stan")
model$print()
set.seed(5)
fitr <- model$variational(data = datar,
                          seed=7,output_samples = 5000,
                          eval_elbo=50,grad_samples =20,
                          elbo_samples = 20,
                          algorithm = "meanfield",
                          output_dir = NULL,
                          iter = 3000,
                          adapt_iter = 50,
                          save_latent_dynamics=TRUE,
                          tol_rel_obj = 10^-1)

predtheta<-array(fitr$summary("theta")$mean,dim=c(M,n_d,J))
beta1<-matrix(fitr$summary("beta")$mean,nrow=J,byrow=FALSE)
thetaest<-c()

for(d in 1:M)
{
  for( r in 1:num_sub[d])
  {
    thetaest<-rbind(thetaest,predtheta[d,r,])
  }
}

library(nnet)
DataC<-data.frame(matrix(,nrow=nrow(thetaest),ncol=J+2))
DataC[,1:J]<-thetaest
DataC[,(J+1)]<-as.factor(DataComb$CL4_regroup)
DataC[,(J+2)]<-as.factor(DataComb$SubjectID)


colnames(DataC)<-c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","CL4_regroup","SubjectID")
a1<-names(which(summary(DataC$CL4_regroup)<20))
#a1<-names(which(summary(as.factor(DataComb$CL4_regroup))<20));a1
summary(DataC)
w1<-c()
w1<-rep(1,nrow(DataComb))
w1[which(DataComb$CL4_regroup %in% a1)]<-4
m1 <- multinom(CL4_regroup ~T1+T2+T3+T4+T5+T6+T7+T8+T9+T10+T11+SubjectID-1,data=DataC,weights=w1)
t1<-predict(m1, newdata = DataC, "class")

tab1<-table(t1,DataComb$CL4_regroup);tab1
sum(diag(tab1))/sum(tab1)


Data<-DataComb[,10:118]
tot_sub<-rowSums(Data)

summary(Data)
DataRM<-Data
DataRM$CL4_regroup<-as.factor(DataComb$CL4_regroup)
DataRM$SubjectID<-as.factor(DataComb$SubjectID)
dim(DataRM)

Datatr<-c()
for(i in 1:length(num_sub))
{
  Datatr<-rbind(Datatr,DataRM[which(DataRM$SubjectID==subid[i])[1:(num_sub[i]-1)],])
}


Datatst<-c()
for(i in 1:length(num_sub))
{
  Datatst<-rbind(Datatst,DataRM[which(DataRM$SubjectID==subid[i])[num_sub[i]],])
}

dim(Datatr)
library(mclogit)
actTr<-Datatr$CL4_regroup
actTst<-Datatst$CL4_regroup
subtr<-Datatr$SubjectID
subtst<-Datatst$SubjectID

thetatr<-c()
for(d in 1:M)
{
  for( r in 1:(num_sub[d]-1))
  {
    thetatr<-rbind(thetatr,predtheta[d,r,])
  }
}

thetatst<-c()
for(d in 1:M)
{
  for( r in num_sub[d])
  {
    thetatst<-rbind(thetatst,predtheta[d,r,])
  }
}

Datatr<-data.frame(matrix(,nrow=nrow(thetatr),ncol=J+2))
Datatr[,1:J]<-thetatr

Datatr[,(J+1)]<-as.factor(actTr)
Datatr[,(J+2)]<-as.factor(subtr)
colnames(Datatr)<-c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","CL4_regroup","SubjectID")
summary(Datatr)
w1<-c()
w1<-rep(1,nrow(Datatr))
w1[which(actTr %in% a1)]<-4
m1 <- multinom(CL4_regroup ~T1+T2+T3+T4+T5+T6+T7+T8+T9+T10+T11+SubjectID-1,data=Datatr,weights=w1)


t1<-predict(m1, newdata = Datatr, "class")

tabtr<-table(t1,actTr);tabtr
sum(diag(tabtr))/sum(tabtr)
Datatst<-data.frame(matrix(,nrow=nrow(thetatst),ncol=J+2))
Datatst[,1:J]<-thetatst
Datatst[,(J+1)]<-as.factor(actTst)
Datatst[,(J+2)]<-as.factor(subtst)
colnames(Datatst)<-c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","CL4_regroup","SubjectID")
summary(Datatst)
t1<-predict(m1, newdata = Datatst, "class")

tabtst<-table(t1,actTst);tabtst
sum(diag(tabtst))/sum(tabtst)
summary(actTst)
tab1<-table(t,actTr);tab1
sum(diag(tab1))/sum(tab1)
sigma<-array(fitr$summary("sigma")$mean,dim=c(M,J,J))
summary(m1)

######### RM-Multinomial ######
Data<-DataComb[,10:118]
tot_sub<-rowSums(Data)

summary(Data)
DataRM<-Data
DataRM$CL4_regroup<-as.factor(DataComb$CL4_regroup)
DataRM$SubjectID<-as.factor(DataComb$SubjectID)
dim(DataRM)

Datatr<-c()
for(i in 1:length(num_sub))
{
  Datatr<-rbind(Datatr,DataRM[which(DataRM$SubjectID==subid[i])[1:(num_sub[i]-1)],])
}


Datatst<-c()
for(i in 1:length(num_sub))
{
  Datatst<-rbind(Datatst,DataRM[which(DataRM$SubjectID==subid[i])[num_sub[i]],])
}

dim(Datatr)
library(mclogit)
actTr<-Datatr$CL4_regroup
actTst<-Datatst$CL4_regroup
subtr<-Datatr$SubjectID
subtst<-Datatst$SubjectID
model <- mblogit(CL4_regroup ~ Bacteroides + Faecalibacterium + Lachnospiracea_incertae_sedis + unclassified_Lachnospiraceae + Prevotella + unclassified_Ruminococcaceae + Alistipes + unclassified_Bacteria + Blautia + Parabacteroides + unclassified_Clostridiales + Roseburia + Oscillibacter + Ruminococcus + unclassified_Proteobacteria + Phascolarctobacterium + Barnesiella + Paraprevotella + Akkermansia + unclassified_Alphaproteobacteria + Clostridium.IV + Parasutterella + unclassified_Porphyromonadaceae + Sutterella + unclassified_Bacteroidetes + Clostridium.XI + Collinsella + Escherichia.Shigella + Coprococcus + unclassified_Erysipelotrichaceae + unclassified_Firmicutes + Dorea + Clostridium.XlVa + Odoribacter + unclassified_Bacteroidales + unclassified_Fusobacteriaceae + Klebsiella + unclassified_Prevotellaceae + Clostridium.sensu.stricto + Butyricicoccus + Erysipelotrichaceae_incertae_sedis + Clostridium.XlVb + Streptococcus + Staphylococcus + Butyricimonas + Flavonifractor + unclassified_Coriobacteriaceae + unclassified_Burkholderiales + Turicibacter + Victivallis + Clostridium.XVIII + unclassified_Acidaminococcaceae + Acidaminococcus + Veillonella + Streptophyta + Lactobacillus + Bilophila + unclassified_Clostridiales_Incertae.Sedis.XIII + Anaerotruncus + Eggerthella + Oxalobacter + Gordonibacter + Fusobacterium + Slackia + Enterobacter + Desulfovibrio + Anaerovorax + Bifidobacterium + Pseudoflavonifractor + Cloacibacillus + unclassified_Desulfovibrionaceae + Holdemania + unclassified_Clostridia + Enterococcus + Anaerofilum + Actinomyces + Lactococcus + Coraliomargarita + TM7_genera_incertae_sedis + Peptococcus + Coprobacillus + Granulicatella + Campylobacter + Subdoligranulum + Lactonifactor + Corynebacterium + Gemella, 
                 random=~1|SubjectID,
                 data = Datatr,weights=w1)


t2<-predict(model,type="response",data=Datatr)
t3<-apply(t2,1,which.max)
table(t3)
c1<-colnames(t2)
t4<-c()
for(i in 1:length(t3))
{
  t4[which(t3==i)]<-c1[i]
}

t4<-factor(t4,levels=names(table(DataComb$CL4_regroup)))
tab4<-table(t4,Datatr$CL4_regroup);tab4

sum(diag(tab4))/sum(tab4)


t2<-predict(model,type="response",newdata=Datatst)
t3<-apply(t2,1,which.max)
table(t3)
c1<-colnames(t2)
t4<-c()
for(i in 1:length(t3))
{
  t4[which(t3==i)]<-c1[i]
}

t4<-factor(t4,levels=names(table(DataComb$CL4_regroup)))
tab4<-table(t4,Datatst$CL4_regroup);tab4

sum(diag(tab4))/sum(tab4)


b1<-beta1
library(lessR)
actual<-DataComb$CL4_regroup
summary(actual)
BarChart(actual,ylab =" Count-Health Status",xlab="Health Status")

topics<- c("Topic 1", "Topic 2", "Topic 3","Topic 4","Topic 5", "Topic 6", "Topic 7","Topic 8","Topic 9", "Topic 10", "Topic 11")
pred_name=colnames(Data1)[-110]
num_top_words <- 6
top_words1<-c()
prop<-c()
for(i in 1:4)
{
  
  colnames(b1)<-pred_name
  # Get the indices of the top words for each topic
  top_words_indices <- apply(b1, 2, order, decreasing = TRUE)[i, 1:num_top_words]
  # Get the actual top words for each topic
  top_words<-names(top_words_indices)
  print(i)
  print(b1[i,which(colnames(b1)%in% top_words)])
  prop<-as.vector(c(prop,b1[i,which(colnames(b1)%in% top_words)]))
  top_words1<- c(top_words1,top_words)
}
# Define the number of top words to consider for each topic


# Convert data to long format
df <- data.frame(topic = rep(topics[1:4], each = num_top_words),
                 Word = top_words1,
                 proportion = prop)
library(RColorBrewer)
color_palette <- c(
  brewer.pal(8, "Set2"),
  brewer.pal(2, "Dark2")
)  # Change the palette and number of colors as desired

# Barplot with proportions
ggplot(df, aes(x = topic, y = proportion, fill = Word)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +ylim(0,1)+
  geom_text(aes(label = round(proportion*100,2)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 1.5) +
  labs(x = "Topics", y = "Proportions", title = "Distribution of words over J topics") +
  scale_fill_manual(values = c("#FFD9D9", "#D9FFD9", "#D9D9FF", "#FFFFCC", "#FFCCFF", "#CCFFFF")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))



pred_name=colnames(Data1)[-110]
num_top_words <- 6
top_words1<-c()
prop<-c()
for(i in 5:8)
{
  
  colnames(b1)<-pred_name
  # Get the indices of the top words for each topic
  top_words_indices <- apply(b1, 2, order, decreasing = TRUE)[i, 1:num_top_words]
  # Get the actual top words for each topic
  top_words<-names(top_words_indices)
  print(i)
  print(b1[i,which(colnames(b1)%in% top_words)])
  prop<-as.vector(c(prop,b1[i,which(colnames(b1)%in% top_words)]))
  top_words1<- c(top_words1,top_words)
}
# Define the number of top words to consider for each 
df <- data.frame(topic = rep(topics[5:8], each = num_top_words),
                 Word = top_words1,
                 proportion = prop)
library(RColorBrewer)
color_palette <- c(
  brewer.pal(8, "Set2"),
  brewer.pal(2, "Dark2")
)  # Change the palette and number of colors as desired

# Barplot with proportions
ggplot(df, aes(x = topic, y = proportion, fill = Word)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +ylim(0,1)+
  geom_text(aes(label = round(proportion*100,2)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 1.5) +
  labs(x = "Topics", y = "Proportions", title = "Distribution of words over J topics") +
  scale_fill_manual(values = c("#FFD9D9", "#D9FFD9", "#D9D9FF", "#FFFFCC", "#FFCCFF", "#CCFFFF")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))


pred_name=colnames(Data1)[-110]
num_top_words <- 6
top_words1<-c()
prop<-c()

for(i in 9:11)
{
  colnames(b1)<-pred_name
  # Get the indices of the top words for each topic
  top_words_indices <- apply(b1, 2, order, decreasing = TRUE)[i, 1:num_top_words]
  # Get the actual top words for each topic
  top_words<-names(top_words_indices)
  print(i)
  print(b1[i,which(colnames(b1)%in% top_words)])
  prop<-as.vector(c(prop,b1[i,which(colnames(b1)%in% top_words)]))
  top_words1<- c(top_words1,top_words)
}

# Define the number of top words to consider for each 
df <- data.frame(topic = rep(topics[11:9], each = num_top_words),
                 Word = top_words1,
                 proportion = prop)
library(RColorBrewer)
color_palette <- c(
  brewer.pal(8, "Set2"),
  brewer.pal(2, "Dark2")
)  # Change the palette and number of colors as desired

# Barplot with proportions
ggplot(df, aes(x = topic, y = proportion, fill = Word)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +ylim(0,1)+
  geom_text(aes(label = round(proportion*100,2)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 1.5) +
  labs(x = "Topics", y = "Proportions", title = "Distribution of words over J topics") +
  scale_fill_manual(values = c("#FFD9D9", "#D9FFD9", "#D9D9FF", "#FFFFCC", "#FFCCFF", "#CCFFFF")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))


library(lessR)
actual<-DataComb$CL4_regroup
summary(actual)
BarChart(actual,ylab =" Count-Health Status",xlab="Health Status",values_cut=0,values_position="out")
