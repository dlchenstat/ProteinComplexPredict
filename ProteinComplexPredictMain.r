#R code developed by Donglai Chen for manuscript "A Label Free Mass spectrometry Method to Predict Endogenous Protein Complex Composition". 
# This code returns cluster IDs for split fitted profile clustering for concatenated profiles of two IEX alone, two SEC alone and IEX+SEC.
# Execute file "ProteinComplexPredictFunctions.r" first
# Input files are two IEX datasets, two SEC datasets, list of SEC cytosolic proteins, list of protein contaminents in IEX, (if proteins are already filtered, those two lists are not needed) Gaussian fitting results of four datasets, list of known protein complexes.
# 
# List of input files
# The IEX input file
# IEX_bio1_common_cytosol.csv
# IEX_bio2_common_cytosol.csv
# CytoContaminents.csv (remove low quality IEX proteins)
# 
# The IEX peak detection file (generated from Gaussian fitting Matlab code)
# peakloc-iex-bio1-uma-2015aug.csv
# peakloc-iex-bio2-uma-2015aug.csv
# 
# The SEC input file
# SEC_Bio1_nov.csv
# SEC_Bio2_nov.csv
# SEC_Bio1_Bio2_cytosol_list.csv (select cytosolic proteins in SEC)
# 
# The SEC peak detection file  (generated from Gaussian fitting Matlab code)
# peakloc-sec1-uma-2015nov.csv
# peakloc-sec2-uma-2015nov.csv

conta=read.csv("CytoContaminents.csv")
#read IEX data
iex1 = read.csv("IEX_bio1_common_cytosol.csv", row.names = 1)
iex2 = read.csv("IEX_bio2_common_cytosol.csv", row.names = 1)
iexloc1 = read.csv("peakloc-iex-bio1-uma-2015aug.csv", row.names = 1)
iexloc2 = read.csv("peakloc-iex-bio2-uma-2015aug.csv", row.names = 1)
iex1=iex1[!substr(rownames(iex1),1,9) %in% conta[,1], ]
iex2=iex2[!substr(rownames(iex2),1,9) %in% conta[,1], ]
iexloc1=iexloc1[!substr(rownames(iexloc1),1,9) %in% conta[,1],]
iexloc2=iexloc2[!substr(rownames(iexloc2),1,9) %in% conta[,1],]
iex1 = std(iex1)
iex2 = std(iex2)

iex1=deldup(iex1)
iex2=deldup(iex2)

iexloc1 = dellocname(iexloc1)
iexloc2 = dellocname(iexloc2)

reprotiex = reproprot(iex1, iex2, iexloc1, iexloc2)
#number of reproducible proteins 1569
#delete proteins with duplicate first 9 characters, splicing variants
reiex1 = iex1[reprotiex, ]
reiex2 = iex2[reprotiex, ]

#read SEC data
seclist = read.csv("SEC_Bio1_Bio2_cytosol_list.csv", row.names = 1)
sec1 = read.csv("SEC_Bio1_nov.csv", row.names = 1)
sec2 = read.csv("SEC_Bio2_nov.csv", row.names = 1)
sec1 = sec1[substr(rownames(sec1), 1, 9) %in% rownames(seclist), ]
sec2 = sec2[substr(rownames(sec2), 1, 9) %in% rownames(seclist), ]
sec1=sec1[!substr(rownames(sec1),1,9) %in% conta[,1], ]
sec2=sec2[!substr(rownames(sec2),1,9) %in% conta[,1], ]
sec1 = std(sec1)
sec2 = std(sec2)
secloc1 = read.csv("peakloc-sec1-uma-2015nov.csv", row.names = 1)
secloc2 = read.csv("peakloc-sec2-uma-2015nov.csv", row.names = 1)
secloc1=secloc1[!substr(rownames(secloc1),1,9) %in% conta[,1],]
secloc2=secloc2[!substr(rownames(secloc2),1,9) %in% conta[,1],]
#reproducible by 2 fraction shift
reprotsec = reproprot(sec1, sec2, secloc1, secloc2, 2)

resec1 = sec1[reprotsec, ]
resec2 = sec2[reprotsec, ]

secloc1 = dellocname(secloc1)
secloc2 = dellocname(secloc2)
#use only first 9 character of the protein IDs
rownames(resec1) = substr(rownames(resec1), 1, 9)
rownames(resec2) = substr(rownames(resec2), 1, 9)


known = read.csv("Knowns2.0.csv")
rownames(known) = known[, 1]
knownsec=known[rownames(known) != "AT3G51890_2",]
rownames(knownsec)= substr(knownsec[, 1],1,9)

protname = intersect(rownames(reiex1), rownames(resec1))

#only keep reproducible proteins
reiex1 = reiex1[protname, ]
reiex2 = reiex2[protname, ]
resec1 = resec1[protname, ]
resec2 = resec2[protname, ]


#find reproducible peaks
iexnumpeak = numrepropeak(reiex1, reiex2, iexloc1, iexloc2, 4)
secnumpeak = numrepropeak(resec1, resec2, secloc1, secloc2, 2)

#create split fitted profiles for IEX, SEC
spfitiex1 = fitspprofile(reiex1, iexloc1, iexnumpeak[, 3:7], T)
spfitiex2 = fitspprofile(reiex2, iexloc2, iexnumpeak[, c(3, 8:11)], T)
spfitsec1 = fitspprofile(resec1, secloc1, secnumpeak[, 3:7], F)
spfitsec2 = fitspprofile(resec2, secloc2, secnumpeak[, c(3, 8:11)], F)

#combine split fitted profiles
spfitiex12sec12 = cbind(spfitiex1, spfitiex2 , matrix(0, ncol = 2 * ncol(resec1)))

for (i in rownames(spfitiex12sec12)) {
  if (sum(i %in% rownames(reiex1))) {
    spfitiex12sec12[i, (2 * ncol(reiex1) + 1):ncol(spfitiex12sec12)] = cbind(spfitsec1[i, ], spfitsec2[i, ])
  } else{
    j = substr(i, 1, 9)
    spfitiex12sec12[i, (2 * ncol(reiex1) + 1):ncol(spfitiex12sec12)] = cbind(spfitsec1[j, ], spfitsec2[j, ])
  }
}

#return a table of cluster IDs for different number of clusters
#results of IEX, SEC combined
numclust = seq(20, 600, 10)
clusterIDtablespfitIEX = matrix(nrow = nrow(spfitiex12sec12), ncol = length(numclust))
for (i in numclust) {
  j = which(numclust == i)
  clusterIDtablespfitIEX[, j] = orderlabel(spfitiex12sec12, i)[, 1]
  cat(i, "\r")
}
colnames(clusterIDtablespfitIEX) = paste("X", numclust, sep = "")
rownames(clusterIDtablespfitIEX) = rownames(spfitiex12sec12)


#cluster results of IEX only
spfitiex12=cbind(spfitiex1,spfitiex2)

clusterIDtablespfitIEX12 = matrix(nrow = nrow(spfitiex12), ncol = length(numclust))
for (i in numclust) {
  j = which(numclust == i)
  clusterIDtablespfitIEX12[, j] = orderlabel(spfitiex12, i)[, 1]
  cat(i, "\r")
}
colnames(clusterIDtablespfitIEX12) = paste("X", numclust, sep = "")
rownames(clusterIDtablespfitIEX12) = rownames(spfitiex12)


#cluster results of sec only
spfitsec12=cbind(spfitsec1,spfitsec2)

clusterIDtablespfitSEC12 = matrix(nrow = nrow(spfitsec12), ncol = length(numclust))
for (i in numclust) {
  j = which(numclust == i)
  clusterIDtablespfitSEC12[, j] = orderlabel(spfitsec12, i)[, 1]
  cat(i, "\r")
}
colnames(clusterIDtablespfitSEC12) = paste("X", numclust, sep = "")
rownames(clusterIDtablespfitSEC12) = rownames(spfitsec12)


write.table(clusterIDtablespfitIEX, "clusterID_IEX+SEC.csv")
write.table(clusterIDtablespfitIEX12, "clusterID_IEX.csv")
write.table(clusterIDtablespfitSEC12, "clusterID_SEC.csv")
