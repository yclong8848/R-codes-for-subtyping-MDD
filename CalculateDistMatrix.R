


########################## Read the file MexSNP ################################# 

data1 = read.csv("MexSNP.csv", header = T, colClasses = "character")
dim(data1) # 399 83900 

L = dim(data1)[1]   # 399 subjects
M = dim(data1)[2]   # 83900 = 2 columns + 83898 common SNPs


########################## Normalized HD distance Matrix #########################

library("e1071")

DistM = matrix(nrow=L, ncol=L) 
subject = list()

for (i in 1:L){
  subject[[i]]=as.character(data1[i, 3:M])
}


for (i in 1:(L-1)){
  for (j in (i+1):L){
    DistM[i, j] =  hamming.distance(subject[[i]], subject[[j]])
  }
}

DistM = DistM/(M-2)

write.csv(DistM, file="DistMatrix83898.csv")

########################## Normalized HD distance Matrix for 19 SNPs ###############

attach(data1)

SNP19 = cbind(exm167893, exm1508600, exm1616604, exm875366, exm283068, exm445797, exm1441979, exm669085, exm75804, exm1355772, exm1044842, exm1435859, exm1325307, exm1505393, exm1369092, exm1293569, exm782507, exm2249659, exm2275308 )

l = dim(SNP19 )[1]   # 399 subjects
m = dim(SNP19 )[2]   # 19 significant SNPs



library("e1071")

Distm = matrix(nrow=l, ncol=l) 
subj = list()

for (i in 1:l){
  subj[[i]]=as.character(SNP19[i, 1:m])
}


for (i in 1:(l-1)){
  for (j in (i+1):l){
    Distm[i, j] =  hamming.distance(subj[[i]], subj[[j]])
  }
}

Distm = Distm/m

write.csv(Distm, file="DistMatrix19.csv")



########################## Normalized HD distance Matrix for MDS and do MDS #########################


library("e1071")

DistM = matrix(nrow=L, ncol=L) 
subject = list()

for (i in 1:L){
  subject[[i]]=as.character(data1[i, 3:M])
}


for (i in 1:L){
  for (j in 1:L){
    DistM[i, j] =  hamming.distance(subject[[i]], subject[[j]])
  }
}

DistM = DistM/(M-2)

write.csv(DistM, file="DistMatrix83898_MDS.csv")

fit <- cmdscale(DistM, eig = TRUE, k = 2)

x <- fit$points[, 1]
y <- fit$points[, 2]

plot(c(x[1:126], x[323:399]), c(y[1:126], y[323:399]), pch = 19, xlab="MDS Coordinate 1", ylab="MDS Coordinate 2", col = 28)
points(x[127:322], y[127:322], pch = 19, xlab="MDS Coordinate 1", ylab="MDS Coordinate 2", col = 34)




########################## Normalized HD distance Matrix for 19 SNPs for MDS and do MDS ###############

attach(data1)

SNP19 = cbind(exm167893, exm1508600, exm1616604, exm875366, exm283068, exm445797, exm1441979, exm669085, exm75804, exm1355772, exm1044842, exm1435859, exm1325307, exm1505393, exm1369092, exm1293569, exm782507, exm2249659, exm2275308 )

l = dim(SNP19 )[1]   # 399 subjects
m = dim(SNP19 )[2]   # 19 significant SNPs



library("e1071")

Distm = matrix(nrow=l, ncol=l) 
subj = list()

for (i in 1:l){
  subj[[i]]=as.character(SNP19[i, 1:m])
}


for (i in 1:l){
  for (j in 1:l){
    Distm[i, j] =  hamming.distance(subj[[i]], subj[[j]])
  }
}

Distm = Distm/m

write.csv(Distm, file="DistMatrix19_MDS.csv")



fit <- cmdscale(Distm, eig = TRUE, k = 2)

x <- fit$points[, 1]
y <- fit$points[, 2]

plot(c(x[1:126], x[323:399]), c(y[1:126], y[323:399]), pch = 19, xlab="MDS Coordinate 1", ylab="MDS Coordinate 2", col = 28)
points(x[127:322], y[127:322], pch = 19, xlab="MDS Coordinate 1", ylab="MDS Coordinate 2", col = 34)

####################################### Mean and SD of distances ########################################


library("e1071")

subject = list()

for (i in 1:L){
  subject[[i]]=as.character(data1[i, 3:M])
}


MDD = list()

for (i in 1:126){
  MDD[[i]]=subject[[i]]  
}
for (i in 127:203){
  MDD[[i]]=subject[[i+196]]
}
# MDD 203
Control = list()

for (i in 1:196){
  Control[[i]]=subject[[i+126]]  
}
# Control 196

DistM1 = matrix(nrow=203, ncol=203) 

for (i in 1:203){
  for (j in 1:203){
    DistM1[i, j] =  hamming.distance(MDD[[i]], MDD[[j]]) 
  }
}

VS1 <- as.vector(DistM1/83898)

mean(VS1)
sd(VS1)


DistM2 = matrix(nrow=196, ncol=196) 

for (i in 1:196){
  for (j in 1:196){
    DistM2[i, j] =  hamming.distance(Control[[i]], Control[[j]]) 
  }
}

VS2 <- as.vector(DistM2/83898)

mean(VS2)
sd(VS2)


DistM3 = matrix(nrow=203, ncol=196) 

for (i in 1:203){
  for (j in 1:196){
    DistM3[i, j] =  hamming.distance(MDD[[i]], Control[[j]]) 
  }
}

VS3 <- as.vector(DistM3/83898)

mean(VS3)
sd(VS3)


####################################### Mean and SD of distances for 19 SNPs ############################


attach(data1)

SNP19 = cbind(exm167893, exm1508600, exm1616604, exm875366, exm283068, exm445797, exm1441979, exm669085, exm75804, exm1355772, exm1044842, exm1435859, exm1325307, exm1505393, exm1369092, exm1293569, exm782507, exm2249659, exm2275308 )

l = dim(SNP19 )[1]   # 399 subjects
m = dim(SNP19 )[2]   # 19 significant SNPs




library("e1071")

subject = list()

for (i in 1:l){
  subject[[i]]=as.character(SNP19[i, 1:m])
}


MDD = list()

for (i in 1:126){
  MDD[[i]]=subject[[i]]  
}
for (i in 127:203){
  MDD[[i]]=subject[[i+196]]
}
# MDD 203
Control = list()

for (i in 1:196){
  Control[[i]]=subject[[i+126]]  
}
# Control 196

DistM1 = matrix(nrow=203, ncol=203) 

for (i in 1:203){
  for (j in 1:203){
    DistM1[i, j] =  hamming.distance(MDD[[i]], MDD[[j]]) 
  }
}

VS1 <- as.vector(DistM1/19)

mean(VS1)
sd(VS1)


DistM2 = matrix(nrow=196, ncol=196) 

for (i in 1:196){
  for (j in 1:196){
    DistM2[i, j] =  hamming.distance(Control[[i]], Control[[j]]) 
  }
}

VS2 <- as.vector(DistM2/19)

mean(VS2)
sd(VS2)


DistM3 = matrix(nrow=203, ncol=196) 

for (i in 1:203){
  for (j in 1:196){
    DistM3[i, j] =  hamming.distance(MDD[[i]], Control[[j]]) 
  }
}

VS3 <- as.vector(DistM3/19)

mean(VS3)
sd(VS3)

#################################################



