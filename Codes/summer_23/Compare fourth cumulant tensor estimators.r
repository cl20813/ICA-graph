# Compare kurtosis tensor estimator

## run library and function file first !

reps = 50
frob_my_m4hat = rep(0,reps)
frob2_Auddy_pg8_m4 = rep(0,reps)
frob3_Auddy_pg20_m4 = rep(0,reps)
test = rep(0,reps)

## simulation
for (r in 1:reps){

## data generation
##########
    d=4; n=4000
    Z = matrix(rep(0,d*n),n,d)
    for(j in 1:d){
    Z[,j] = matrix(rgamma(n*1, shape=j, rate=(3)), n, 1)
    }
    Z = Z - matrix(1, n, 1) %*% apply(Z, 2, mean)
    res = eigen(var(Z))
    Z = Z %*%  (res$vectors)  %*% diag(res$values^{-1/2}) %*% t(res$vectors)
    round(var(Z),2)
    trueW = gen_trueW(d)
    X = Z %*% solve(t(trueW))
  

## Whitening X
    var.X = var(X)
    res = eigen(var.X)
    Whitening_mat = (res$vectors) %*% diag(res$values^{-1/2}) %*% t((res$vectors))
    X = X %*% Whitening_mat
    round( var(X), 2)
    unwh_trueW =  trueW %*% t( solve(Whitening_mat) )    
##########
##########

    true_kurt_tensor = gamma_kurt_tensor(unwh_trueW)
    true_kurt_tensor = kurt_to_mat(true_kurt_tensor) 
    #kurt_to_mat transform dxdxdx tensor into d^2 x d^2 matrix
        
    My_m4_hat  = my_m4_hat(X)
    My_m4_hat = kurt_to_mat(My_m4_hat) 

    Auddy_pg20 = Auddy_M4hat_pg20(X)
    Auddy_pg20_m4 = kurt_to_mat(Auddy_pg20)
  
    # Auddy3M = fold(Auddy3M, row_idx= c(1,2), col_idx=c(3,4), modes=c(d,d,d,d))

    frob_my_m4hat[r] =  frob_norm(true_kurt_tensor,My_m4_hat)
    frob3_Auddy_pg20_m4[r]  = frob_norm(true_kurt_tensor, Auddy_pg20_m4)
}

## boxplots 

frob= c(frob_my_m4hat, frob3_Auddy_pg20_m4)
tmp = (rep(c('my_m4hat', 'Auddy_pg20_m4'), each=reps))
tmp = factor(tmp, levels= unique(tmp)) # preserve order of factors
data = data.frame(frob,tmp)
ggplot(data = data, mapping= aes(x = tmp, y = frob), ylim=c(0,0.8) ) +
    geom_boxplot(alpha=0.3, width=1/length(unique(data$tmp)), 
    color= c("green","blue"), fill=c("green","blue"),
    outlier.shape = NA ) +
     labs(title = "Frobenius norm, d=4, n=4000",
       y = "" , x="") 
