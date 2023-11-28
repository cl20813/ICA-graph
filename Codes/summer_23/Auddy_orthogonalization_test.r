d=5; n=100
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
round(cov(X),2)
res2 = eigen(cov(X))
X = X %*%  (res2$vectors)  %*% diag(res2$values^{-1/2}) %*% t(res2$vectors)

# check round(cov(X),2) is closed to an identity matrix.
round(cov(X),2)

n = nrow(X)
a_j = solve(t(trueW))[,1]
a_j = a_j/ norm(a_j,"2")
for(i in 1:n){
    X[i,] = X[i,] - as.numeric( t(a_j) %*% X[i,] )*a_j
}

test = rep(0,n)
for(i in 1:n){
    test[i] = X[i,] %*% a_j
}
sum(test) # if zero, X is orthogonalized with respect to a_j. 
round(cov(X),2)

res3 = eigen(cov(X))
res3$values 




