library(factoextra)
library(ggplot2)
library(pracma)
library(readr)

# pca in WGS
RMD_WGS = read_csv('results/ssm_remove/RMD/RMD_sm_sinOG_muts_1.csv')
RMD.data = RMD_WGS[, 3:ncol(RMD_WGS)]
table(is.na(RMD.data))
rownames(RMD.data) = RMD_WGS$sample_id
head(RMD.data)
table(colMeans(RMD.data)==0)

# 3. Compute PCA in R using prcomp()
RMDm = as.matrix(RMD.data)
dim(RMDm)
res=prcomp(RMDm,center=T,scale.=F)
names(res)

fviz_eig(res, addlabels = TRUE, ncp= 100)

# the same but manually
vars = apply(res$x, 2, var)
props <- vars / sum(vars)
sum(props[1:100]) # hasta el PC100 se explica un 50% de la variabilidad 



files = list.files(path = 'results/ssm_remove/RMD/', pattern = '.csv', full.names = TRUE)
files
result_muts = lapply(files, function(x){
  rmd_mut = read.csv(x)
  name = sub('.*(RMD_sm_sinOG_muts_.*.csv)', '\\1', x)
  final = data.frame('sample_id'=rmd_mut$sample_id, 'cancer_type' = rmd_mut$cancer_type)
  
  # projection
  rmd_mut.data = rmd_mut[,3:ncol(rmd_mut)]
  pred <- predict(res, newdata=rmd_mut.data)
  dim(pred)
  proj = cbind(final, pred[,1:100])
  print(proj[1:5,1:5])
  
  write.csv(proj, paste0('results/ssm_remove/RMD_PCA/PCA_', name), row.names = FALSE )

})

