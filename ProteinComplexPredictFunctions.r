#R code developed by Donglai Chen for manuscript "A Label Free Mass spectrometry Method to Predict Endogenous Protein Complex Composition". 
#Load libraries and defined function first

#load packages for heat map
library(plotrix)
library(gplots)

#remove duplicate protein IDs, only use first 9 characters of IDs.
deldup = function(iex) {
  rniex1 = substr(row.names(iex), 1, 9)
  crniex1 = table(rniex1)
  dupiex1 = names(crniex1[crniex1 > 1])
  delind = numeric(0)
  for (i in dupiex1) {
    dupmat = iex[rniex1 == i,]
    rsum = apply(dupmat, 1, sum)
    delind = c(delind, which(rniex1 == i)[setdiff(1:dim(dupmat)[1], which.max(rsum))])
  }
  iex = iex[-delind, ]
  rniex1 = substr(row.names(iex), 1, 9)
  rownames(iex) = rniex1
  iex
}

#divide raw profiles by their global max
std = function(prot) {
  sum1 = apply(prot, 1, sum)
  prot = prot[sum1 > 0, ]
  max1 = apply(prot, 1, max)
  prot = prot / max1
}

#Only use first 9 characters of protein IDs
dellocname = function(iexloc1) {
  allname = rownames(iexloc1)
  xx = table(substr(rownames(iexloc1), 1, 9))
  xx = xx[xx > 1]
  for (i in names(xx)) {
    dupname = rownames(iexloc1)[substr(rownames(iexloc1), 1, 9) == i]
    maxname = dupname[which.max(apply(iexloc1[dupname, 15:18], 1, max, na.rm =
                                        T))]
    delname = setdiff(dupname, maxname)
    allname = setdiff(allname, delname)
  }
  iexloc1 = iexloc1[allname, ]
  rownames(iexloc1) = substr(rownames(iexloc1), 1, 9)
  iexloc1
}

#select reproducible protein, proteins within fraction shift of 4
reproprot = function(iex1, iex2, loc1, loc2, fracdiff = 4) {
  protname = intersect(rownames(iex1), rownames(iex2))
  protnew = protname
  for (i in protname) {
    if (sum(iex1[i, ]) == 0 | sum(iex2[i, ]) == 0) {
      protnew = setdiff(protnew, i)
      next
    }
    if (sum(!is.na(loc1[i, 11:14])) > 0) {
      pkloc1 = loc1[i, 11:14]
      pkloc1 = pkloc1[!is.na(pkloc1)]
    } else{
      pkloc1 = which.max(iex1[i, ])
    }
    if (sum(!is.na(loc2[i, 11:14])) > 0) {
      pkloc2 = loc2[i, 11:14]
      pkloc2 = pkloc2[!is.na(pkloc2)]
    } else{
      pkloc2 = which.max(iex2[i, ])
    }
    rppkloc1 = rep(pkloc1, each = length(pkloc2))
    rppkloc2 = rep(pkloc2, length(pkloc1))
    pkdist = abs(rppkloc1 - rppkloc2)
    if (min(pkdist) > fracdiff) {
      protnew = setdiff(protnew, i)
    }
  }
  protnew
}

#fitted deconvolute clustering,  For IEX we will split the peaks into a different group and have multiple entries. For SEC can we use the peak with the largest mass. This will be the peak with the smallest number.
orderlabel = function(commoniex1,
                      numclust,
                      weighto = -apply(commoniex1, 1, which.max)) {
  clustiex1 = hclust(dist(commoniex1), method = "ward.D2")
  labiex1 = cutree(clustiex1, k = numclust)
  reordiex1 = reorder(as.dendrogram(clustiex1), weighto)
  relabiex1 = labiex1[labels(reordiex1)]
  map = numclust:1
  names(map) = unique(relabiex1)
  data.frame(clusterID = map[as.character(labiex1)], row.names = rownames(commoniex1))
  
}

#For IEX, return separate entries of fitted reproducible peaks. For SEC, return the fitted reproducible peak that is on the most left.
fitspprofile = function(reiex1, reiexloc1, iexnumpeak1, ifiex) {
  fitted = c()
  for (i in rownames(reiex1)) {
    if (sum(i %in% rownames(reiexloc1)) > 0 & iexnumpeak1[i, 1] > 0) {
      numfrac = sum(!is.na(reiexloc1[i, 11:14]))
      if (ifiex) {
        repropeaks = (1:numfrac)[reiexloc1[i, 11:14] %in% iexnumpeak1[i, 2:5]]
        #mh = max(reiexloc1[i, 15 + repropeaks - 1])
      } else{#for SEC, choose the most left reproducible peak
        ifrep=reiexloc1[i, 11:14] %in% iexnumpeak1[i, 2:5]
        repropeaks = (1:numfrac)[ifrep][which.min(reiexloc1[i, 11:14][ifrep])]
        #mh = max(reiexloc1[i, 15:18][ifrep])
      }
      for (j in repropeaks) {
        k = which(repropeaks == j)
        tobind = reiex1[i, ] * 0 +  exp(-(1:ncol(reiex1) -reiexloc1[i, 11 + j - 1]) ^ 2 / reiexloc1[i, 19 + j - 1] ^ 2) 
        if (length(repropeaks) > 1) {
          rownames(tobind) = paste(i, k, sep = '_')
        }
        tobind=tobind/max(tobind)
        fitted = rbind(fitted, tobind)
      }
    } else{
      fitted = rbind(fitted, reiex1[i, ])
    }
  }
  fitted
}

#return a table of reproducible peak location. Reproducible peaks are those within 4 fractions in IEX and within 2 fractions in SEC in 2 bio replicates.
numrepropeak = function(iex1, iex2, iexloc1, iexloc2, peakdif) {
  commoniex = intersect(rownames(iex1), rownames(iex2))
  numpeakiex = data.frame(matrix(0, nrow = length(commoniex), ncol = 11), row.names =
                            commoniex)
  colnames(numpeakiex) = c(
    "num peak iex1",
    "num peak iex2",
    "num reproducible peak",
    paste("IEX1 reproducible peak", 1:4),
    paste("IEX2 reproducible peak", 1:4)
  )
  for (i in commoniex) {
    numpeakiex[i, 1] = sum(!is.na(iexloc1[i, 11:14]))
    numpeakiex[i, 2] = sum(!is.na(iexloc2[i, 11:14]))
    if (sum(!is.na(iexloc1[i, 11:14])) > 0 &
        sum(!is.na(iexloc2[i, 11:14])) > 0) {
      peakiex1 = iexloc1[i, 11:14]
      peakiex1 = peakiex1[!is.na(peakiex1)]
      peakiex2 = iexloc2[i, 11:14]
      peakiex2 = peakiex2[!is.na(peakiex2)]
      repropeak = 0
      for (j in peakiex1) {
        for (k in peakiex2) {
          if (abs(j - k) <= peakdif) {
            peakiex1 = setdiff(peakiex1, j)
            peakiex2 = setdiff(peakiex2, k)
            repropeak = repropeak + 1
            numpeakiex[i, 4 + repropeak - 1] = j
            numpeakiex[i, 8 + repropeak - 1] = k
            break
          }
        }
      }
      numpeakiex[i, 3] = repropeak
    }
  }
  numpeakiex
}


#calculate the number of known proteins that are in the same cluster
clustab = function(reiex1,
                   known,
                   numclust,
                   clusterIDtableIEX1) {
  clustab1 = matrix(0, nrow = length(table(known[, 2])), ncol = length(numclust) *
                      2 + 1)
  rownames(clustab1) = names(table(known[, 2]))
  colnames(clustab1) = c(
    "num prot in this complex",
    paste("X", numclust, sep = ''),
    paste("clusterSize", numclust, sep = '')
  )
  for (j in names(table(known[, 2]))) {
    clustab1[j, 1] = sum(rownames(reiex1) %in% as.character(known[known[, 2] %in% j, 1])   )
  }
  for (k in numclust) {
    labiex1 = clusterIDtableIEX1[, paste("X", k, sep = '')]
    names(labiex1) = rownames(clusterIDtableIEX1)
    for (j in names(table(known[, 2]))) {
      numeachclust = table(labiex1[names(labiex1) %in% as.character(known[known[, 2] %in% j, 1])])
      clustab1[j , paste("X", k, sep = '')] = max(numeachclust)
      maxclust1 = names(numeachclust)[which.max(numeachclust)]
      clustab1[j, paste("clusterSize", k, sep = '')] = sum(labiex1 == maxclust1)
      
    }
    cat(k, "\r")
  }
  clustab1
}

#calculate purity and intactness
ratioclust = function(clustab1, useful, numclust) {
  ratio1 = matrix(nrow = length(useful), ncol = length(numclust))
  colnames(ratio1) = paste("X", as.character(numclust), sep = '')
  rownames(ratio1) = useful
  for (i in useful) {
    ratio1[i, paste("X", as.character(numclust), sep = '')] = unlist(clustab1[i, paste("X", as.character(numclust), sep =
                                                                                         '')] / clustab1[i, paste("clusterSize", as.character(numclust), sep = '')])
  }
  ratio1
}
#mean distance to center in known proteins
clusdist = function(reiex1,
                    known,
                    numclust,
                    clusterIDtableIEX1) {
  clustab1 = matrix(0, nrow = length(table(known[, 2])), ncol = length(numclust))
  rownames(clustab1) = names(table(known[, 2]))
  colnames(clustab1) = paste("X", numclust, sep = '')
  for (k in numclust) {
    labiex1 = clusterIDtableIEX1[, paste("X", k, sep = '')]
    names(labiex1) = rownames(clusterIDtableIEX1)
    for (j in names(table(known[, 2]))) {
      numeachclust = table(labiex1[as.character(known[known[, 2] %in% j, 1])])
      subclust = intersect(names(labiex1)[labiex1 == names(numeachclust[which.max(numeachclust)])], as.character(known[known[, 2] %in% j, 1]))
      subiex = reiex1[subclust, ]
      centeriex = apply(subiex, 2, mean)
      clustab1[j , paste("X", k, sep = '')] = mean(apply(t(t(subiex) - centeriex) ^
                                                           2, 1, sum))
    }
    cat(k, "\r")
  }
  clustab1
}

#return a list of distances within clusters for different number of clusters. Then you can show a box plot
distinclust = function(reiex1,
                       numclust,
                       clusterIDtableIEX1,
                       dista = "cmax") {
  i = 1
  distlist = list()
  for (k in numclust) {
    labiex1 = clusterIDtableIEX1[, paste("X", k, sep = '')]
    names(labiex1) = rownames(clusterIDtableIEX1)
    eachclu = numeric(length(unique(labiex1)))
    for (j in sort(unique(labiex1))) {
      subiex = reiex1[names(labiex1[labiex1 == j]), ]
      centeriex = apply(subiex, 2, mean)
      eachclu[j] = switch(
        dista,
        cmax = max(apply(t(
          t(subiex) - centeriex
        ) ^ 2, 1, sum)),
        cmean = mean(apply(t(
          t(subiex) - centeriex
        ) ^ 2, 1, sum)),
        csum = sum(apply(t(
          t(subiex) - centeriex
        ) ^ 2, 1, sum)),
        pmean = mean(dist(subiex)),
        pmax = max(dist(subiex)),
        psum = sum(dist(subiex))
      )
      names(eachclu)[j] = paste("X", j, sep = '')
    }
    distlist[[i]] = eachclu
    names(distlist)[i] = paste("X", k, sep = '')
    i = i + 1
    cat(k, "\r")
  }
  distlist
}

