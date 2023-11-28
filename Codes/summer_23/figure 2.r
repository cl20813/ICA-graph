## Investigate relative errors in a function of d.
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
  unwh_trueW_frob = sqrt(sum(unwh_trueW^2)) # frobenius norm of W after whitening X
  out = list( X=X, unwh_trueW=unwh_trueW, Whitening_mat=Whitening_mat, unwh_trueW_frob=unwh_trueW_frob)
  return(out)
}

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

# investigate W_hat in a function of d
n=2000
D = 30
rel.error = rep(0,D)
W.list = list()
unwh_trueW.list = list()
for(i in 4:D){
    set.seed(1)
    d = i
    data = X.gen(d,n)
    X = data$X
    unwh_trueW = data$unwh_trueW
    unwh_trueW_frob = data$unwh_trueW_frob
    set.seed(1)
    W = Auddy_fastICA_mod2(X,unwh_trueW,L=10,T=400)
    rel.error[i] = frob_norm(unwh_trueW, W)/unwh_trueW_frob
 }

rel.error = rel.error[4:D]
d = 4:D
data = data.frame(d=d, rel.error=rel.error)
ggplot(data=data, mapping= aes(x=d , y=rel.error )) +
    geom_line ( color= "green") +
    # geom_point(shape = 21, color = "black", fill= "#69b3a2", size=2)+
    #theme_ipsum() +
    ggtitle("Relative errors d:4~30")



