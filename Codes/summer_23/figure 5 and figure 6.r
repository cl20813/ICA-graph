## figure 5 and 6
## run "library and function_6_25.r" file first!

## data generation
X.gen = function(d,n){
  # d=5; n=4000
  Z = matrix(rep(0,d*n),n,d)
  for(j in 1:d){
    # Z[,j] = matrix(rnorm(n*1, 0,1), n, 1)    
    Z[,j] = matrix(rgamma(n*1, shape=j, rate=(3)), n, 1)
    }
  Z = Z - matrix(1, n, 1) %*% apply(Z, 2, mean) # centering Z
  res = eigen(var(Z)) # Whitening Z
  Z = Z %*%  (res$vectors)  %*% diag(res$values^{-1/2}) %*% t(res$vectors)
  trueW = gen_trueW(d) # generate true W
  X = Z %*% solve(t(trueW)) # genereate X

## Whitening X
  X = Auddy_prewhiten(X,trueW)$X
  unwh_trueW = Auddy_prewhiten(X,trueW)$unwh_trueW
  Whitening_mat = Auddy_prewhiten(X,trueW)$Whitening_mat
  unwh_trueW_frob = sqrt(sum(unwh_trueW^2)) #frobenius norm of W after whitening X
  out = list( X=X, unwh_trueW=unwh_trueW, Whitening_mat=Whitening_mat, unwh_trueW_frob=unwh_trueW_frob)
  return(out)
}

## figure 5

mytest= function(X){
    d= dim(X)[2]   
    Mhat = my_m4_hat(X)
    L=10
    top.s.list = rep(0,L)
    top.u.list = list()
    for (l in 1:L){
        G = matrix(rnorm(d^2),d,d)    
        mat = matrix(0,d,d)
        for(i in 1:d){
            for(j in 1:d){
                mat[i,j] = sum(diag(Mhat[,,i,j]@data %*%G)) # matrix(Mhat[,,i,j]@data,1,d^2) %*% matrix(G,d^2,1)
            }
        }
        res = svd(mat)
        top.s.list[l] = sum(res$d)
        top.u.list[[l]] = res$u[,]
    }
    ind = order(top.s.list, decreasing=TRUE)
    W = top.u.list[[ ind[1] ]]
    W = sign_hungarian(unwh_trueW,W)
    return(W)
}

# initial estimate  
n_iter = 50
W_D5 = rep(0,n_iter)
W_D8 = rep(0,n_iter )
W_D10 = rep(0,n_iter )
W_D25 = rep(0,n_iter )

for (k in 1:n_iter){

    X = X.gen(5,2000)$X ; unwh_trueW = X.gen(5,2000)$unwh_trueW
    unwh_trueW1 =unwh_trueW  
    unwh_trueW_frob1 = X.gen(5,2000)$unwh_trueW_frob
  W_d5 = mytest(X)
    X = X.gen(8,2000)$X ; unwh_trueW = X.gen(8,2000)$unwh_trueW
    unwh_trueW2 =unwh_trueW
    unwh_trueW_frob2 = X.gen(5,2000)$unwh_trueW_frob
  W_d8 = mytest(X)
    X = X.gen(10,2000)$X ; unwh_trueW = X.gen(10,2000)$unwh_trueW
    unwh_trueW3 =unwh_trueW
    unwh_trueW_frob3 = X.gen(5,2000)$unwh_trueW_frob
  W_d10 = mytest(X)
    X = X.gen(25,2000)$X ; unwh_trueW = X.gen(25,2000)$unwh_trueW
    unwh_trueW4 =unwh_trueW
    unwh_trueW_frob4 = X.gen(5,2000)$unwh_trueW_frob
  W_d25 = mytest(X)
  
  W_D5[k]  = frob_norm(unwh_trueW1, W_d5)/unwh_trueW_frob1
  W_D8[k]  = frob_norm(unwh_trueW2, W_d8)/unwh_trueW_frob2 
  W_D10[k]  = frob_norm(unwh_trueW3, W_d10)/unwh_trueW_frob3 
  W_D25[k]  = frob_norm(unwh_trueW4, W_d25)/unwh_trueW_frob4 
}

frob= c(W_D5, W_D8,W_D10,W_D25)
tmp = (rep(c('W_D5', 'W_D8', 'W_D10', 'W_D25'), each=n_iter) )
tmp = factor(tmp, levels= unique(tmp)) # preserve the order of factors
data = data.frame(frob,tmp)
ggplot(data = data, mapping= aes(x = tmp, y = frob), ylim=c(0,0.8) ) +
    geom_boxplot(alpha=0.3, width=1/length(unique(data$tmp)), 
    color= c("green","blue","orange","red"), fill=c("green","blue","orange","red"),
    outlier.shape = NA ) +
     labs(title = "",
       y = "" , x="") 


## figure 6.

n_iter = 50
met_1 = rep(0,n_iter)
met_2 = rep(0,n_iter )
met_3 = rep(0,n_iter )
met_4 = rep(0,n_iter )
met_5 = rep(0,n_iter )
L=10

for (k in 1:n_iter){

## data generation
  d=15; n=4000
  Z = matrix(rep(0,d*n),n,d)
  for(j in 1:d){
    # Z[,j] = matrix(rnorm(n*1, 0,1), n, 1)    
    Z[,j] = matrix(rgamma(n*1, shape=j, rate=(3)), n, 1)
    }
  Z = Z - matrix(1, n, 1) %*% apply(Z, 2, mean) # centering Z
  res = eigen(var(Z)) # Whitening Z
  Z = Z %*%  (res$vectors)  %*% diag(res$values^{-1/2}) %*% t(res$vectors)
  trueW = gen_trueW(d) # generate true W
  X = Z %*% solve(t(trueW)) # genereate X

## Whitening X
  X = Auddy_prewhiten(X,trueW)$X
  unwh_trueW = Auddy_prewhiten(X,trueW)$unwh_trueW
  Whitening_mat = Auddy_prewhiten(X,trueW)$Whitening_mat
  unwh_trueW_frob = sqrt(sum(unwh_trueW^2)) #frobenius norm of W after whitening X

  Mhat = my_m4_hat(X)
    top.s.list = rep(0,L)
    top.u.list = list()
    for (l in 1:L){
        G = matrix(rnorm(d^2),d,d)    
        mat = matrix(0,d,d)
        for(i in 1:d){
            for(j in 1:d){
                mat[i,j] = sum(diag(Mhat[,,i,j]@data %*%G)) # matrix(Mhat[,,i,j]@data,1,d^2) %*% matrix(G,d^2,1)
            }
        }
        res = svd(mat)
        top.s.list[l] = sum(res$d)
        top.u.list[[l]] = res$u[,]
    }
    ind = order(top.s.list, decreasing=TRUE)
    W = top.u.list[[ ind[1] ]]
    W2 = sign_hungarian(unwh_trueW,W)
    rg = matrix(rnorm(d^2),d,d)
    JadeW = JADE(X, n.comp = d, eps = 1e-06, maxiter = 100, na.action = na.fail)$W
    # fIC1 = fICA(X, g="pow3",  method="sym")$W
    

    fIC1 = fICA(X, g="pow3", init=JadeW, method="def")$W 
    fIC2 = fICA(X, g="pow3", init=W2, method="def")$W  
    fIC3 = fICA(X, g="pow3", init=rg, method="def")$W 
    fIC4 = fICA(X, g="pow3", init=matrix(rexp(d^2,3),d,d), method="def")$W 
    fIC5 = fICA(X, g="pow3", init= unwh_trueW, method="def")$W 

   # frob_norm(unwh_trueW,W)
  met_1[k] = frob_norm(unwh_trueW, fIC1)/unwh_trueW_frob 
  met_2[k] = frob_norm(unwh_trueW, fIC2)/unwh_trueW_frob 
  met_3[k] = frob_norm(unwh_trueW,  fIC3)/unwh_trueW_frob 
  met_4[k] = frob_norm(unwh_trueW,  fIC4)/unwh_trueW_frob 
  met_5[k] = frob_norm(unwh_trueW,  fIC4)/unwh_trueW_frob 

}
frob= c(met_1, met_2, met_3, met_4,met_5)
tmp = (rep(c('JADE W', 'Gaussian', 'algorithm 1', 'rexp(3)', 'unwh_trueW'), each=n_iter) )
tmp = factor(tmp, levels= unique(tmp)) # preserve the order of factors
data = data.frame(frob,tmp)
ggplot(data = data, mapping= aes(x = tmp, y = frob), ylim=c(0,0.8) ) +
    geom_boxplot(alpha=0.3, width=1/length(unique(data$tmp)), 
    color= c("green","blue","orange","red","black"), 
    fill=c("green","blue","orange","red","black"),
    outlier.shape = NA ) +
     labs(title = "",
       y = "" , x="") 

###################
###################
sum(met_1)
sum(met_2)
sum(met_3)
sum(met_4)
sum(met_5)









