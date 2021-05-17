library(unix)
unix::rlimit_as(100*10^9)
library(readxl)
library(readr)
library(dplyr)
library(reshape2)
library(glmmLasso)
library(pbapply)
library(reshape2)
library(ggplot2)
library("robCompositions")
library("lme4")
library("faraway")
library("lmerTest")
library("plot.matrix")
library("viridis")
library("glmnet")
options(stringsAsFactors = FALSE)

###################Microbial data
DOM <- read_excel("MicrobialEcobyClass.xlsx", 
                                  col_types = c("text", "numeric", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric"))
#formatting the data
melted_data <-melt(DOM,
                   id=c("Class"),
                   variable.name = "seasons",
                  value.name = "class_value")
data_micro <-dcast(melted_data, seasons~ Class, value.var = "class_value")
##########Chemical Ecology Data
#PCA_NEG
PCA_NEG <- read_excel("ChemicalEcology.xlsx", 
                      sheet = "DOM_Neg")

#Formatting the data
 melted_data_chem <-melt(PCA_NEG,
                    id=c("Properties"),
                    variable.name = "seasons",
                    value.name = "chemical_value")
 
 data_chem <-dcast(melted_data_chem, seasons~ Properties, value.var = "chemical_value")

 #preparing the data as a composition
 prop_data1 <-t(apply(data_chem[,-1], 1, function(x) x/sum(x)))
 data_chemical <- data.frame(seasons =data_chem[,1],prop_data1)
 
 join_df <- merge(data_micro,  data_chemical, by="seasons")
 
 join_df$seasons <- as.character(join_df$seasons)
 
 join_df_new <- join_df %>%
         mutate(
                 splits = strsplit(seasons, "_")
         ) %>% 
         rowwise() %>% 
         mutate(
                 seasons = splits[1],
                 depth = splits[2])%>%
         select(-c(22))

 ####################
 # Microbial vrs PCA- NEG
 #############


 #Dimensions of the data
 D1 <- length(dimnames(data_micro)[[2]])-1
 D2 <- length(dimnames( data_chemical)[[2]])-1

 #Microbial data (father = DOM)
 father <- join_df_new[,2:(D1+1)]+ 0.00001
 names(father) <- as.character(seq(1, ncol(father)))

 # Chemical data (mother = DOM)
 mother <- join_df_new[,(D1+2):(ncol(join_df_new)-1)]+ 0.00001
 names(mother) <- as.character(seq(1, ncol(mother)))
 
 ###########################################
 # Effect of chemical on Microbial ecology
 ###########################################
 pval <- coef <- matrix(NA, ncol =D2, nrow = D1)
for(i in 1:D1){
         for(j in  1:D2){
                 zfath <- as.matrix(pivotCoord(cbind(father[, i], father[, -i])))
                 zmoth <- (pivotCoord(cbind(mother[, j], mother[, -j])))
                 names(zmoth) <- paste0("moth",c(j,seq(1,ncol(zmoth))[-j]))[-D2]
                 data_model <- data.frame(father = zfath[,1],zmoth, seasons= join_df_new$seasons, depth= join_df_new$depth)
                 formula <- as.formula(paste0("father ~ ",paste(names(zmoth),collapse = "+",sep = "")))
                 res <- summary(lm(formula, data=data_model))$coefficients
                 pval[i,j] <- res[2, 4] # entry of the p-value in the matrix
                 coef[i,j] <- res[2, 1] # entry of the coefficient in the matrix
                 print(c(i,j))
         }
}

 ##plot of the pvalues
 
 longData<-melt(pval)
 longData <- longData%>%
         mutate(adj_pval = ifelse(value <0.05, 0,1))
 
 ggplot(longData, aes(x = Var2, y = Var1)) + 
         geom_tile(aes(fill=adj_pval)) + 
         scale_fill_gradient(low="grey90", high="red") +
         labs(x="Microial Ecology", y="Chemical Ecology", title="P-values") +
         theme_bw() + 
         theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
               axis.text.y=element_text(size=9),
               plot.title=element_text(size=11))+
         geom_text(aes(label=round(value,1)))
 
 
 #Plot of the coefficients
 longData1<-melt(coef)
 
 ggplot(longData1, aes(x = Var2, y = Var1)) + 
         geom_tile(aes(fill=value)) + 
         scale_fill_gradient(low="grey90", high="red") +
         labs(x="Microbial Ecology", y="Chemical Ecology", title="Coefficients") +
         theme_bw() + 
         theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
               axis.text.y=element_text(size=9),
               plot.title=element_text(size=11))+
         geom_text(aes(label=round(value,1)))
 

###########################################
# Effect of micro_bio on chemical ecology
###########################################
pval <- coef <- matrix(NA, ncol =D1, nrow = D2)
for(i in 1:D1){
        for(j in  1:D2){
                zfath <- (pivotCoord(cbind(father[, i], father[, -i])))
                zmoth <- (pivotCoord(cbind(mother[, j], mother[, -j])))
                names(zfath) <- paste0("fath",c(i,seq(1,ncol(zfath))[-i]))[-D1]
                data_model <- data.frame(mother = zmoth[,1],zfath, seasons= join_df_new$seasons, depth= join_df_new$depth)
                # formula <- as.formula(paste0("father ~ ",paste("moth",seq(1,ncol(zmoth)-2),collapse = "+",sep = ""),"+",paste0("moth",ncol(zmoth)-1)))
                formula <- as.formula(paste0("mother ~ ",paste(names(zfath),collapse = "+",sep = "")))
                res <- cv.glmnet(x=as.matrix(zfath), y= zmoth[,1], alpha=1)
                res_new <- glmnet(x=as.matrix(zfath), y= zmoth[,1], lambda = res$lambda.min, alpha=1)
                #pval[i,j] <- res[2, 4] # entry of the p-value in the matrix
                coef[j,i] <- coef(res_new)[2, 1] # entry of the coefficient in the matrix
                if(coef[j,i] ==0){
                        pval[j,i] = 1
                }else{
                   formula <-  as.formula(paste0("mother ~ ",paste(rownames(coef(res_new))[which(coef(res_new)!=0)][-1],collapse = "+",sep = "")))   
                linmodel <- summary(lm(formula, data = data_model))$coefficients
                pval[j,i] = linmodel[2,4]
                coef[j,i] = linmodel[2,1]
                   }
                print(c(i,j))
        }
}

##plot of the pvalues

longData<-melt(pval)
longData <- longData%>%
        mutate(adj_pval = ifelse(value <0.05, 0,1))

ggplot(longData, aes(x = Var2, y = Var1)) + 
        geom_tile(aes(fill=adj_pval)) + 
        scale_fill_gradient(low="grey90", high="red") +
        labs(x="Chemical Ecology", y="Microbial Ecology", title="P-values") +
        theme_bw() + 
        theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                           axis.text.y=element_text(size=9),
                           plot.title=element_text(size=11))+
        geom_text(aes(label=round(value,1)))


#Plot of the coefficients
longData1<-melt(coef)

ggplot(longData1, aes(x = Var2, y = Var1)) + 
        geom_tile(aes(fill=value)) + 
        scale_fill_gradient(low="grey90", high="red") +
        labs(x="Chemical Ecology", y="Microbial Ecology", title="P-values") +
        theme_bw() + 
        theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
              axis.text.y=element_text(size=9),
              plot.title=element_text(size=11))+
        geom_text(aes(label=round(value,1)))




########################
# MICROBIAL AND PCA POSITIVE
#######################

#PCA_POS
PCA_POS <- read_excel("ChemicalEcology.xlsx", 
                      sheet = "DOM_Pos")

#Formatting the data
melted_data_chem <-melt(PCA_POS,
                        id=c("Properties"),
                        variable.name = "seasons",
                        value.name = "chemical_value")

data_chem <-dcast(melted_data_chem, seasons~ Properties, value.var = "chemical_value")

#preparing the data as a composition
prop_data1 <-t(apply(data_chem[,-1], 1, function(x) x/sum(x)))
data_chemical <- data.frame(seasons =data_chem[,1],prop_data1)

join_df <- merge(data_micro,  data_chemical, by="seasons")

join_df$seasons <- as.character(join_df$seasons)

join_df_new <- join_df %>%
        mutate(
                splits = strsplit(seasons, "_")
        ) %>% 
        rowwise() %>% 
        mutate(
                seasons = splits[1],
                depth = splits[2])%>%
        select(-c(22))


#Dimensions of the data
D1 <- length(dimnames(data_micro)[[2]])-1
D2 <- length(dimnames( data_chemical)[[2]])-1

#Microbial data (father = DOM)
father <- join_df_new[,2:(D1+1)]+ 0.00001
names(father) <- as.character(seq(1, ncol(father)))

# Chemical data (mother = DOM)
mother <- join_df_new[,(D1+2):(ncol(join_df_new)-1)]+ 0.00001
names(mother) <- as.character(seq(1, ncol(mother)))

###########################################
# Effect of chemical on Microbial ecology
###########################################
pval <- coef <- matrix(NA, ncol =D2, nrow = D1)
for(i in 1:D1){
        for(j in  1:D2){
                zfath <- as.matrix(pivotCoord(cbind(father[, i], father[, -i])))
                zmoth <- (pivotCoord(cbind(mother[, j], mother[, -j])))
                names(zmoth) <- paste0("moth",c(j,seq(1,ncol(zmoth))[-j]))[-D2]
                data_model <- data.frame(father = zfath[,1],zmoth, seasons= join_df_new$seasons, depth= join_df_new$depth)
                formula <- as.formula(paste0("father ~ ",paste(names(zmoth),collapse = "+",sep = "")))
                res <- summary(lm(formula, data=data_model))$coefficients
                pval[i,j] <- res[2, 4] # entry of the p-value in the matrix
                coef[i,j] <- res[2, 1] # entry of the coefficient in the matrix
                print(c(i,j))
        }
}

##plot of the pvalues

longData<-melt(pval)
longData <- longData%>%
        mutate(adj_pval = ifelse(value <0.05, 0,1))

ggplot(longData, aes(x = Var2, y = Var1)) + 
        geom_tile(aes(fill=adj_pval)) + 
        scale_fill_gradient(low="grey90", high="red") +
        labs(x="Microial Ecology", y="Chemical Ecology", title="P-values") +
        theme_bw() + 
        theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
              axis.text.y=element_text(size=9),
              plot.title=element_text(size=11))+
        geom_text(aes(label=round(value,1)))


#Plot of the coefficients
longData1<-melt(coef)

ggplot(longData1, aes(x = Var2, y = Var1)) + 
        geom_tile(aes(fill=value)) + 
        scale_fill_gradient(low="grey90", high="red") +
        labs(x="Microbial Ecology", y="Chemical Ecology", title="Coefficients") +
        theme_bw() + 
        theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
              axis.text.y=element_text(size=9),
              plot.title=element_text(size=11))+
        geom_text(aes(label=round(value,1)))


###########################################
# Effect of micro_bio on chemical ecology
###########################################
pval <- coef <- matrix(NA, ncol =D1, nrow = D2)
for(i in 1:D1){
        for(j in  1:D2){
                zfath <- (pivotCoord(cbind(father[, i], father[, -i])))
                zmoth <- (pivotCoord(cbind(mother[, j], mother[, -j])))
                names(zfath) <- paste0("fath",c(i,seq(1,ncol(zfath))[-i]))[-D1]
                data_model <- data.frame(mother = zmoth[,1],zfath, seasons= join_df_new$seasons, depth= join_df_new$depth)
                # formula <- as.formula(paste0("father ~ ",paste("moth",seq(1,ncol(zmoth)-2),collapse = "+",sep = ""),"+",paste0("moth",ncol(zmoth)-1)))
                formula <- as.formula(paste0("mother ~ ",paste(names(zfath),collapse = "+",sep = "")))
                res <- cv.glmnet(x=as.matrix(zfath), y= zmoth[,1], alpha=1)
                res_new <- glmnet(x=as.matrix(zfath), y= zmoth[,1], lambda = res$lambda.min, alpha=1)
                #pval[i,j] <- res[2, 4] # entry of the p-value in the matrix
                coef[j,i] <- coef(res_new)[2, 1] # entry of the coefficient in the matrix
                if(coef[j,i] ==0){
                        pval[j,i] = 1
                }else{
                        formula <-  as.formula(paste0("mother ~ ",paste(rownames(coef(res_new))[which(coef(res_new)!=0)][-1],collapse = "+",sep = "")))   
                        linmodel <- summary(lm(formula, data = data_model))$coefficients
                        pval[j,i] = linmodel[2,4]
                        coef[j,i] = linmodel[2,1]
                }
                print(c(i,j))
        }
}

##plot of the pvalues

longData<-melt(pval)
longData <- longData%>%
        mutate(adj_pval = ifelse(value <0.05, 0,1))

ggplot(longData, aes(x = Var2, y = Var1)) + 
        geom_tile(aes(fill=adj_pval)) + 
        scale_fill_gradient(low="grey90", high="red") +
        labs(x="Chemical Ecology", y="Microbial Ecology", title="P-values") +
        theme_bw() + 
        theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
              axis.text.y=element_text(size=9),
              plot.title=element_text(size=11))+
        geom_text(aes(label=round(value,1)))


#Plot of the coefficients
longData1<-melt(coef)

ggplot(longData1, aes(x = Var2, y = Var1)) + 
        geom_tile(aes(fill=value)) + 
        scale_fill_gradient(low="grey90", high="red") +
        labs(x="Chemical Ecology", y="Microbial Ecology", title="P-values") +
        theme_bw() + 
        theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
              axis.text.y=element_text(size=9),
              plot.title=element_text(size=11))+
        geom_text(aes(label=round(value,1)))


 