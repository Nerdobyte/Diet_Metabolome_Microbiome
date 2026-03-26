rm(list=ls())

library(missForest)
library(ggfortify)
library(tidyverse)
library(stringr)

GUT2 = read.csv(file="input_for_missforest_T20_ALLmol_nodup_20250816.csv", header = T, sep=",")
GUT2 = GUT2[,-c(1,6:19)]
GUT2 = GUT2 %>% filter(!(Location == "Oral"))
GUT = GUT2
GUT$Location = as.factor(GUT$Location)
GUT$Diet = as.factor(GUT$Diet)
GUT$PatientID = as.factor(GUT$PatientID)
summary(GUT)

GUT.n.Butyrate = GUT[,c(1:4,5)]
GUT.Valine = GUT[,c(1:4,6)]
GUT.Propionate = GUT[,c(1:4,7)]
GUT.Alanine = GUT[,c(1:4,8)]
GUT.Acetate = GUT[,c(1:4,9)]
GUT.Glutamine = GUT[,c(1:4,10)]
GUT.Methionine = GUT[,c(1:4,11)]
GUT.Aspartate = GUT[,c(1:4,12)]
GUT.Asparagine = GUT[,c(1:4,13)]
GUT.Tyrosine = GUT[,c(1:4,14)]
GUT.Phenylalanine = GUT[,c(1:4,15)]
GUT.SR = GUT[,c(1:4,16)]
GUT.Trigonelline = GUT[,c(1:4,17)]
GUT.Fumarate = GUT[,c(1:4,18)]
GUT.Tryptophan = GUT[,c(1:4,19)]
GUT.Histidine = GUT[,c(1:4,20)]
GUT.Formate = GUT[,c(1:4,21)]
GUT.alpha.Glucose_alpha.Maltose_5.24d = GUT[,c(1:4,22)]
GUT.beta.Glucose_3.25dd = GUT[,c(1:4,23)]
GUT.beta.Maltose_3.29dd = GUT[,c(1:4,24)]
GUT.TCDCA = GUT[,c(1:4,25)]
GUT.GCDCA = GUT[,c(1:4,26)]
GUT.TDCA = GUT[,c(1:4,27)]
GUT.GDCA = GUT[,c(1:4,28)]
GUT.TCA = GUT[,c(1:4,29)]
GUT.GCA = GUT[,c(1:4,30)]
GUT.Ciceritol = GUT[,c(1:4,31)]
GUT.Sucrose = GUT[,c(1:4,40)]

GUT.n.Butyrate.imp = missForest(GUT.n.Butyrate)
GUT.Valine.imp = missForest(GUT.Valine)
GUT.Propionate.imp = missForest(GUT.Propionate)
GUT.Alanine.imp = missForest(GUT.Alanine)
GUT.Acetate.imp = missForest(GUT.Acetate)
GUT.Glutamine.imp = missForest(GUT.Glutamine)
GUT.Methionine.imp = missForest(GUT.Methionine)
GUT.Aspartate.imp = missForest(GUT.Aspartate)
GUT.Asparagine.imp = missForest(GUT.Asparagine)
GUT.Tyrosine.imp = missForest(GUT.Tyrosine)
GUT.Phenylalanine.imp = missForest(GUT.Phenylalanine)
GUT.SR.imp = missForest(GUT.SR)
GUT.Trigonelline.imp = missForest(GUT.Trigonelline)
GUT.Fumarate.imp = missForest(GUT.Fumarate)
GUT.Tryptophan.imp = missForest(GUT.Tryptophan)
GUT.Histidine.imp = missForest(GUT.Histidine)
GUT.Formate.imp = missForest(GUT.Formate)
GUT.alpha.Glucose_alpha.Maltose_5.24d.imp = missForest(GUT.alpha.Glucose_alpha.Maltose_5.24d)
GUT.beta.Glucose_3.25dd.imp = missForest(GUT.beta.Glucose_3.25dd)
GUT.beta.Maltose_3.29dd.imp = missForest(GUT.beta.Maltose_3.29dd)
GUT.TCDCA.imp = missForest(GUT.TCDCA)
GUT.GCDCA.imp = missForest(GUT.GCDCA)
GUT.TDCA.imp = missForest(GUT.TDCA)
GUT.GDCA.imp = missForest(GUT.GDCA)
GUT.TCA.imp = missForest(GUT.TCA)
GUT.GCA.imp = missForest(GUT.GCA)
GUT.Ciceritol.imp = missForest(GUT.Ciceritol)
GUT.Sucrose.imp = missForest(GUT.Sucrose)

GUT.imputed = 
  cbind(GUT.n.Butyrate.imp$ximp, 
        GUT.Valine.imp$ximp[,5], 
        GUT.Propionate.imp$ximp[,5],
        GUT.Alanine.imp$ximp[,5], 
        GUT.Acetate.imp$ximp[,5],
        GUT.Glutamine.imp$ximp[,5], 
        GUT.Methionine.imp$ximp[,5], 
        GUT.Aspartate.imp$ximp[,5],
        GUT.Asparagine.imp$ximp[,5], 
        GUT.Tyrosine.imp$ximp[,5],
        GUT.Phenylalanine.imp$ximp[,5],
        GUT.SR.imp$ximp[,5], 
        GUT.Trigonelline.imp$ximp[,5], 
        GUT.Fumarate.imp$ximp[,5],
      GUT.Tryptophan.imp$ximp[,5], 
      GUT.Histidine.imp$ximp[,5], 
      GUT.Formate.imp$ximp[,5],
      GUT.alpha.Glucose_alpha.Maltose_5.24d.imp$ximp[,5],
      GUT.beta.Glucose_3.25dd.imp$ximp[,5],
      GUT.beta.Maltose_3.29dd.imp$ximp[,5],
      GUT.TCDCA.imp$ximp[,5],
      GUT.GCDCA.imp$ximp[,5],
      GUT.TDCA.imp$ximp[,5],
      GUT.GDCA.imp$ximp[,5],
      GUT.TCA.imp$ximp[,5],
      GUT.GCA.imp$ximp[,5],
      GUT.Ciceritol.imp$ximp[,5],
      GUT.Sucrose.imp$ximp[,5])
colnames(GUT.imputed) = colnames(GUT[, c(1:31,40)])
write.csv(GUT.imputed, file="output_for_missforest_T21pt1_ALLmol_nodup.csv")

rm(list = ls(pattern = "^GUT\\..*\\.imp$"))

GUT.imputed = read.csv(file="moremol_analysis/20241123/output_for_missforest_T22pt1_ALLmol_nodup.csv", header = T, sep=",")
GUT.imputed$Location_Timepoint = str_c(GUT.imputed$Location,"_",GUT.imputed$Timepoint)
GUT.imputed = GUT.imputed %>% relocate(Location_Timepoint, .after=Timepoint)
GUT.imputed = GUT.imputed[,-c(1,4:5)]
gutimpute_gather <- GUT.imputed %>% gather("Molecule", "Conc" , 4:39)
gutimpute_spread <- gutimpute_gather %>% spread("PatientID", "Conc")
write.csv(gutimpute_spread, file="output_for_missforest_T21pt2_ALLmol_nodup.csv")
write.csv(gutimpute_spread, file="output_for_missforest_T22pt2_ALLmol_nodup.csv")

########################################################################### PCAs

gutimpute_spread$Diet = as.factor(gutimpute_spread$Diet)
gutimpute_spread$Location_Timepoint = as.factor(gutimpute_spread$Location_Timepoint)

GUT.imputed.old = GUT.imputed
GUT.imputed = read.csv(file="output_for_missforest_T21pt1_ALLmol_nodup.csv", header = T, sep=",")
GUT.imputed = GUT.imputed[,-1]
df1 <- GUT.imputed[,c(5:31)]
pca_res <- prcomp(df1, center = T, scale. = T)
GUT.imputed$Diet = as.factor(GUT.imputed$Diet)
GUT.imputed$Location = as.factor(GUT.imputed$Location)
GUT.imputed$Timepoint = as.factor(GUT.imputed$Timepoint)
GUT.imputed$PatientID = as.factor(GUT.imputed$PatientID)
autoplot(pca_res, data = GUT.imputed, colour = 'PatientID')
autoplot(pca_res, data = GUT.imputed, colour = 'Location')
GUT.imputed$Diet = gsub("0","BC",GUT.imputed$Diet)
GUT.imputed$Diet = gsub("1","SC",GUT.imputed$Diet)
GUT.imputed$Diet = gsub("2","CC",GUT.imputed$Diet)
autoplot(pca_res, data = GUT.imputed, colour = 'Diet', )
autoplot(pca_res, data = GUT.imputed, colour = 'Timepoint')

################################################################################ ORIGINLAB PREP1

aa1 = dplyr::filter(gutimpute_spread, grepl("^n.Butyrate$",Molecule))
aa2 = dplyr::filter(gutimpute_spread, grepl("^Valine$",Molecule))
aa3 = dplyr::filter(gutimpute_spread, grepl("^Propionate$",Molecule))
aa4 = dplyr::filter(gutimpute_spread, grepl("^Alanine$",Molecule))
aa5 = dplyr::filter(gutimpute_spread, grepl("^Acetate$",Molecule))
aa6 = dplyr::filter(gutimpute_spread, grepl("^Glutamine$",Molecule))
aa7 = dplyr::filter(gutimpute_spread, grepl("^Methionine$",Molecule))
aa8 = dplyr::filter(gutimpute_spread, grepl("^Aspartate$",Molecule))
aa9 = dplyr::filter(gutimpute_spread, grepl("^Asparagine$",Molecule))
aa10 = dplyr::filter(gutimpute_spread, grepl("^Tyrosine$",Molecule))
aa11 = dplyr::filter(gutimpute_spread, grepl("^Phenylalanine$",Molecule))
aa12 = dplyr::filter(gutimpute_spread, grepl("^SR$",Molecule))
aa13 = dplyr::filter(gutimpute_spread, grepl("^Trigonelline$",Molecule))
aa14 = dplyr::filter(gutimpute_spread, grepl("^Fumarate$",Molecule))
aa15 = dplyr::filter(gutimpute_spread, grepl("^Tryptophan$",Molecule))
aa16 = dplyr::filter(gutimpute_spread, grepl("^Histidine$",Molecule))
aa17 = dplyr::filter(gutimpute_spread, grepl("^Formate$",Molecule))
aa18 = dplyr::filter(gutimpute_spread, grepl("^alpha.Glucose_alpha.Maltose_5.24d$",Molecule))
aa19 = dplyr::filter(gutimpute_spread, grepl("^beta.Glucose_3.25dd$",Molecule))
aa20 = dplyr::filter(gutimpute_spread, grepl("^beta.Maltose_3.29dd$",Molecule))
aa21 = dplyr::filter(gutimpute_spread, grepl("^TCDCA$",Molecule))
aa22 = dplyr::filter(gutimpute_spread, grepl("^GCDCA$",Molecule))
aa23 = dplyr::filter(gutimpute_spread, grepl("^TDCA$",Molecule))
aa24 = dplyr::filter(gutimpute_spread, grepl("^GDCA$",Molecule))
aa25 = dplyr::filter(gutimpute_spread, grepl("^TCA$",Molecule))
aa26 = dplyr::filter(gutimpute_spread, grepl("^GCA$",Molecule))
aa27 = dplyr::filter(gutimpute_spread, grepl("^Ciceritol_new..gastric.mean.eretic.correction.$",Molecule))
aa28 = dplyr::filter(gutimpute_spread, grepl("^Sucrose..correct.peak._.gastric.mean.eretic.correction.$",Molecule))
aa29 = dplyr::filter(gutimpute_spread, grepl("^Oral_Raffinose$",Molecule))
aa30 = dplyr::filter(gutimpute_spread, grepl("^Oral_Stachyose$",Molecule))
aa31 = dplyr::filter(gutimpute_spread, grepl("^Isoleucine$",Molecule))
aa32 = dplyr::filter(gutimpute_spread, grepl("^Taurine$",Molecule))
aa33 = dplyr::filter(gutimpute_spread, grepl("^Lactate$",Molecule))
aa34 = dplyr::filter(gutimpute_spread, grepl("^Leucine$",Molecule))
aa35 = dplyr::filter(gutimpute_spread, grepl("^Lysine$",Molecule))
aa36 = dplyr::filter(gutimpute_spread, grepl("^Choline$",Molecule))

################################################################################ ORIGINLAB PREP2
######### prep for 3D originlab graphs

query2 = aa36
name_query2 = "aa36"

aa_data = query2
aa_data$means = rowMeans(aa_data[,c(4:16)], na.rm = T)
aa_data = aa_data[,c(1,2,3,17)]
aa_data = aa_data %>% 
  separate(Location_Timepoint, into = c("Location", "Timepoint"), sep="_(?=[^_]+$)")
aa_data$Diet_Location = str_c(aa_data$Diet,"_",aa_data$Location)
aa_data = aa_data %>% relocate(Diet_Location, .after=Location)
aa_data = aa_data[,-c(1,2,5)]
aa_data = aa_data %>% 
  pivot_wider(names_from = Diet_Location, values_from = means)
aa_data = aa_data[order(as.integer(aa_data$Timepoint)),]                
aa_data = aa_data[,c(1,5,9,13,
                     3,7,11,
                     2,6,10,
                     4,8,12)]
aa_data[1,c(5:7)] = aa_data[1,c(2:4)]
aa_data = aa_data[,-c(2:4)]
aa_data$Timepoint = factor(aa_data$Timepoint, levels = c("0",   "15",  "30",  "60",  "120", "135", "150", "165", "180", "240", "300", "360", "480"))
colnames(aa_data) = c("Timepoint","BC_Gastric","SC_Gastric","CC_Gastric","BC_Duodenal","SC_Duodenal","CC_Duodenal","BC_Ileal","SC_Ileal","CC_Ileal")
aa_data = aa_data[-1,]                                                #check oral here
write.csv(aa_data, file=paste0("originlab_inputs_20250213/",name_query2,"_data_frameshiftcorrected.csv"),
          row.names = F, na = "")

ggplot(aa_data, aes(x=Timepoint, y=as.numeric(BC_Gastric), colour=BC_Gastric, group=BC_Gastric)) + geom_point() + geom_line()

############################################################################### BOXPLOTS

library(ggplot2)
library(dplyr)

a = 40

### ORAL ###

GUT.imputed = read.csv(file="moremol_analysis/20241123/output_for_missforest_T15pt1_ALLmol_nodup.csv", header = T, sep=",")
GUT.imputed = GUT.imputed[,-1]

boxwhisk_data = GUT.imputed[,c(1:4,a)]
molecule_name = colnames(boxwhisk_data)[5]
colnames(boxwhisk_data)[5] = "Concentration"
boxwhisk_data = dplyr::filter(boxwhisk_data, grepl("^Oral$",Location))
boxwhisk_data$Timepoint_Diet = str_c(boxwhisk_data$Timepoint,"_",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("0","BC",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("1","SC",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("2","CC",boxwhisk_data$Diet)
boxwhisk_data$Timepoint = factor(boxwhisk_data$Timepoint)

ggplot(boxwhisk_data, aes(x=Timepoint, y=Concentration, fill=Diet)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38","#619CFF")) +
  xlab("Time (minutes)") +
  ylab(paste0("Oral ",molecule_name," Concentration (mmol/L)")) +
  geom_point(aes(shape = Diet, group = Diet),
             position = position_dodge(width = 0.75), size = 2
  ) +
  scale_shape_manual(values = c(15, 16, 17))

### GASTRIC ###

boxwhisk_data = GUT.imputed[,c(1:4,a)]
molecule_name = colnames(boxwhisk_data)[5]
colnames(boxwhisk_data)[5] = "Concentration"
boxwhisk_data = dplyr::filter(boxwhisk_data, grepl("^Gastric$",Location))
boxwhisk_data$Timepoint_Diet = str_c(boxwhisk_data$Timepoint,"_",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("0","BC",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("1","SC",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("2","CC",boxwhisk_data$Diet)
boxwhisk_data$Timepoint = factor(boxwhisk_data$Timepoint)

ggplot(boxwhisk_data, aes(x=Timepoint, y=Concentration, fill=Diet)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38","#619CFF")) +
  xlab("Time (minutes)") +
  ylab(paste0("Gastric ",molecule_name," Concentration (mmol/L)")) +
  geom_point(aes(shape = Diet, group = Diet),
             position = position_dodge(width = 0.75), size = 2
  ) +
    scale_shape_manual(values = c(15, 16, 17))

### DUODENAL ###

boxwhisk_data = GUT.imputed[,c(1:4,a)]
molecule_name = colnames(boxwhisk_data)[5]
colnames(boxwhisk_data)[5] = "Concentration"
boxwhisk_data = dplyr::filter(boxwhisk_data, grepl("^Duodenal$",Location))
boxwhisk_data$Timepoint_Diet = str_c(boxwhisk_data$Timepoint,"_",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("0","BC",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("1","SC",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("2","CC",boxwhisk_data$Diet)
boxwhisk_data$Timepoint = factor(boxwhisk_data$Timepoint)

ggplot(boxwhisk_data, aes(x=Timepoint, y=Concentration, fill=Diet)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38","#619CFF")) +
  xlab("Time (minutes)") +
  ylab(paste0("Duodenal ",molecule_name," Concentration (mmol/L)")) +
  geom_point(aes(shape = Diet, group = Diet),
             position = position_dodge(width = 0.75), size = 2
  ) +
  scale_shape_manual(values = c(15, 16, 17))

### ILEAL ###

boxwhisk_data = GUT.imputed[,c(1:4,a)]
molecule_name = colnames(boxwhisk_data)[5]
colnames(boxwhisk_data)[5] = "Concentration"
boxwhisk_data = dplyr::filter(boxwhisk_data, grepl("^Ileal$",Location))
boxwhisk_data$Timepoint_Diet = str_c(boxwhisk_data$Timepoint,"_",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("0","BC",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("1","SC",boxwhisk_data$Diet)
boxwhisk_data$Diet = gsub("2","CC",boxwhisk_data$Diet)
boxwhisk_data$Timepoint = factor(boxwhisk_data$Timepoint)

ggplot(boxwhisk_data, aes(x=Timepoint, y=Concentration, fill=Diet)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38","#619CFF")) +
  xlab("Time (minutes)") +
  ylab(paste0("Ileal ",molecule_name," Concentration (mmol/L)")) +
  geom_point(aes(shape = Diet, group = Diet),
             position = position_dodge(width = 0.75), size = 2
  ) +
  scale_shape_manual(values = c(15, 16, 17))

################################################################################
################################################################################

########################## FOR HALLA

GUT_halla = gutimpute_spread
GUT_halla = GUT_halla[grepl(paste(c("Ileal","Oral"), collapse="|"), GUT_halla$Location),]
GUT_halla = GUT_halla %>% 
  separate(Location_Timepoint, into = c("Location", "Timepoint"), sep="_(?=[^_]+$)")
#GUT_halla = GUT_halla[,-2]
GUT_halla = GUT_halla %>% gather("PatientID", "Conc" , 5:17)
GUT_halla = GUT_halla %>% spread("Molecule", "Conc")
GUT_halla$Locator_ICL <- paste0(GUT_halla$PatientID,"_",
                               GUT_halla$Diet,"_",
                               GUT_halla$Timepoint)
GUT_halla = GUT_halla %>% relocate(Locator_ICL, .before=Diet)
rownames(GUT_halla) = GUT_halla$Locator_ICL

GUT.micro.metabo = read.csv(file="metadata_25072024_for_use_with_metabolites_matrix.csv", 
                            sep=",",
                            header = T)
rownames(GUT.micro.metabo) = GUT.micro.metabo$Locator_ICL
GUT_halla = merge(GUT_halla, GUT.micro.metabo, by=0, all=TRUE )
GUT_halla = GUT_halla[!is.na(GUT_halla$sample), ]
rownames(GUT_halla) = GUT_halla$sample
GUT_halla = GUT_halla[,-c(1,2,4,5,57:74)]
colnames(GUT_halla)[1] = "Diet"

GUT_halla_D0 = GUT_halla[grepl("0", GUT_halla$Diet),]
GUT_halla_D0 = GUT_halla_D0[,-1]
GUT_halla_D0 = t(GUT_halla_D0)
write.csv(GUT_halla_D0, file="Halla_20240623/GUT_halla_D0.csv")

GUT_halla_D1 = GUT_halla[grepl("1", GUT_halla$Diet),]
GUT_halla_D1 = GUT_halla_D1[,-1]
GUT_halla_D1 = t(GUT_halla_D1)
write.csv(GUT_halla_D1, file="Halla_20240623/GUT_halla_D1.csv")

GUT_halla_D2 = GUT_halla[grepl("2", GUT_halla$Diet),]
GUT_halla_D2 = GUT_halla_D2[,-1]
GUT_halla_D2 = t(GUT_halla_D2)
write.csv(GUT_halla_D2, file="Halla_20240623/GUT_halla_D2.csv")

################################################################################
################################################################################
