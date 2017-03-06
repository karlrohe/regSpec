# This code implements regularized spectral clustering for 
# symmetric, directed, and bipartite (i.e. rectangular A).
# It has the following implemented:
# uses irlba package
# checks if symmetric (currently, uses svd either way)
# uses regularization
# project onto sphere
# use kmeans ++
# plots top eigenvalues and eigenvectors
# returns projected singular vectors and kmeans++ clusters.

# TODO:
# 1) For symmetric calculations, it would be faster to find eigenvectors. 
#     With current R libraries, it looks like one would need to use RcppEigen
#     to access the c++ library Eigen. Then, find the relevant function in that library.
# 2) Currently has adjMat as input.  Should also allow igraph object. 
# 3) should we find k+1 eigenvalues to inspect the eigengap?

#  irlba uses default parameters

# Version 0.1.1.  Mar6,2017  karlrohe@stat.wisc.edu

library(irlba)
library(Matrix)




kmpp <- function(X, k) {
  # kmeans ++ 
  # from https://stat.ethz.ch/pipermail/r-help/2012-January/300051.html
  # Hans Werner
  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)
  xx = rowSums(X^2)
  for (i in 2:k) {
    dm <- distmat(X,xx, X[C, ])
    pr <- apply(dm, 1, min); pr[C] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
  }
  
  return(kmeans(X, X[C, ]))
}

distmat = function(X,xx,y){
  #my code
  # used in kmpp to compute distances
  if(length(y) == ncol(X)){
    xy = X%*%y  
    yy = sum(y^2)
    # because of machine error, some zeros are -10^(-15).  Then, sqrt gives some errors +10^(-13) to make positive.
    dm = sqrt(xx + yy - 2*xy+ 10^(-13)) 
    
  }
  
  if(length(y) != ncol(X)){
    xy = X%*%t(y)
    yy = rowSums(y^2)
    dm = sqrt(matrix(xx, ncol = length(yy), nrow = length(xx)) + matrix(yy, byrow = T,ncol = length(yy), nrow = length(xx)) -2* xy  + 10^(-13))
  }
  return(dm)    		
}


kmppi = function(X,k, inits =10){
  # allow for several restarts of kmeans++
  withss = sum(X^2)
  for(i in 1:inits){
    km = kmpp(X,k)
    #     print(c(i,km$tot.withinss))
    if(km$tot.withinss < withss) {
      returnDat = km
      withss = returnDat$tot.withinss
    }
    
  } 
  return(returnDat)
}






regspec = function(A, k, inits = 10, tau = -1, quiet = F, project = T){
  # A is a matrix or a Matrix.
  # k is the number of clusters.  If A is asymmetric, allows a 2 vector.  
  #   First element is number of sending clusters.  Second element is number of receiving clusters.
  # tau is regularization parameter.  Like k, it can be a 2 vectors.  Default is pretty good.  
  # This code currently uses the irlba package for both symmetric and asymmetric calculations.
  #  This finds a partial SVD.  
  # For symmetric calculations, it would be faster to find eigenvectors. 
  #  irlba uses default parameters
  
  #   ensure arguments to function are reasonable
  sym = F
  if(!isSymmetric(object = A)) if(!quiet) print("adjacency matrix is not symmetric. So, this function will perform a directed analysis akin to di-sim.")
  if(isSymmetric(object = A)){
    sym = T
    if(length(k) > 1){
      if(!quiet) print("k is given more than one argument, but A is symmetric.")
      if(!quiet) print(paste("using k = ", k[1]))
      k = k[1]
    }
  }
  
  # compute the (regularized) graph Laplacian.
  nr = nrow(A); nc = ncol(A)
  rs=rowSums(A); cs = colSums(A)
  E = sum(A)
  
  if(length(tau) == 1) if(tau <0) tau = c(E/nr, E/nc)
  
  #     tau = c(1,1)
  #     if(rs > 2* log(nr)) tau[1] = rs/nr
  #     if(cs > 2* log(cs)) tau[2] = cs/nc
  #   }
  #   
  
  if(length(tau)==1) tau = c(tau,tau)
  
  Drow = Diagonal(n = nr, x = 1/sqrt(rs + tau[1]))
  Dcol = Diagonal(n = nc, x = 1/sqrt(cs + tau[2]))
  tmp = Drow%*%A
  L = tmp %*% Dcol
  
  K = min(k)
  
  # find the singular vectors
  s = irlba(L, nu = K, nv = K)
  
  # project singular vectors onto sphere
  # to prevent 0/0, the "projection" adds a bit in the denominator. 
  if(project == T) {
    nu= t(apply(s$u, 1, function(x, nr) return(x/sqrt(sum(x^2)+length(x)/(100*nr))), nr = nr))
    if(!sym) nv= t(apply(s$v, 1, function(x, nr) return(x/sqrt(sum(x^2) + length(x)/(100*nc) )), nr))
  }
  if(project == F){
    nu = s$u
    if(!sym) nv = s$v
  }
  
  
  # run kmeans++
  kmLeft = kmppi(X = nu, k[1],inits)
  #   kmLeft = kmeans(x = nu, centers = k[1],nstart =  100)
  kmLeft$nu = nu
  if(!sym) {
    if(length(k) ==2) kmRight = kmppi(X = nv, k[2], inits)
    if(length(k) ==1) kmRight = kmppi(X = nv, k[1], inits)
    kmRight$nv = nv
  }
  
  # create object that is returned.
  if(!sym) outDat = list(nu = nu, nv=nv, uclst = kmLeft$cluster, vclst = kmRight$cluster, ucent = kmLeft$center, vcent = kmRight$center)
  if(sym) outDat = list(nu = nu, uclst = kmLeft$cluster, ucent = kmLeft$center)
  
  if(!quiet){
    # make some plots.
    plot(-sort(-s$d)/sum(s$d), main = "top singular values")
    print("press 1 for plot of eigenvectors.")
    x = scan()
    
    if(x!=1) return(outDat)
    
    
    # if you want to plot the data, downsample to 1000 data points.
    if(nr > 1000){
      samp = sample(nr,1000)
    }
    if(nr<1001) samp = 1:nr
    if(sym) plot(as.data.frame(nu[samp,]), col = kmLeft$clust[samp], 
                 main ="leading eigenvectors,\nprojected on sphere and colored by cluster")
    if(!sym) plot(as.data.frame(nu[samp,]), col = kmLeft$clust[samp], 
                  main ="leading left singular vectors,\nprojected on sphere and colored by cluster")
  }
  return(outDat)
}