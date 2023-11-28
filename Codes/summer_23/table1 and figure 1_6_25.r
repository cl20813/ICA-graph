## figure 1 and table 1; test Auddy and Yuan's deflation method
## run "library and function_6_25.r" file first!

n_iter = 50
frob_alg1 = rep(0,n_iter)
frob_alg2 = rep(0,n_iter )
frob_alg3 = rep(0,n_iter )
frob_alg4 = rep(0,n_iter )
frob_alg5 = rep(0,n_iter )
frob_fICA = rep(0,n_iter)

for (k in 1:n_iter){
## data generation
  d=5; n=4000
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
  unwh_trueW_frob = sqrt(sum(unwh_trueW^2)) # frobenius norm of W after whitening X

## Simulation results
 
  fICAW = fICA(X, g="pow3", method="def")$W 
  AuddyW = Auddy_fastICA(X,unwh_trueW,10,400)  # Auddy m4hat Auddy deflation
  myW = Auddy_fastICA_mod2(X,unwh_trueW,10,400) # my_m4hat and my deflation
  myW2 = Auddy_fastICA_mod3(X,unwh_trueW,10,400) # my_m4hat and Auddy deflation
  myW3 = Auddy_fastICA_mod4(X,unwh_trueW,10,400) # Auddy m4hat my deflation
  myW4 = Auddy_fastICA_mod5(X,unwh_trueW,10,400) # (Auddy 4th moment + my m0) and my deflation

  
  frob_alg1[k]  = frob_norm(unwh_trueW, AuddyW)/unwh_trueW_frob 
  frob_alg2[k] = frob_norm(unwh_trueW, myW)/unwh_trueW_frob 
  frob_alg3[k] = frob_norm(unwh_trueW, myW2)/unwh_trueW_frob 
  frob_alg4[k] = frob_norm(unwh_trueW, myW3)/unwh_trueW_frob 
  frob_alg5[k] = frob_norm(unwh_trueW, myW4)/unwh_trueW_frob 
  frob_fICA[k] = frob_norm(unwh_trueW, fICAW)/unwh_trueW_frob
}

## boxplots for relative errors

frob= c(frob_alg1, frob_alg2,frob_alg3,frob_alg4, frob_alg5,frob_fICA)
tmp = (rep(c('frob_alg1', 'frob_alg2', 'frob_alg3', 'frob_alg4', 'frob_alg5', 'frob_fICA'), each=n_iter) )
tmp = factor(tmp, levels= unique(tmp)) # preserve the order of factors
data = data.frame(frob,tmp)
ggplot(data = data, mapping= aes(x = tmp, y = frob), ylim=c(0,0.8) ) +
    geom_boxplot(alpha=0.3, width=1/length(unique(data$tmp)), 
    color= c("green","blue","orange","red","yellow","black"), fill=c("green","blue","orange","red","yellow","black"),
    outlier.shape = NA ) +
     labs(title = "",
       y = "" , x="") 

sum(frob_alg1)