# install.packages("fastmatrix")
# install.packages("rockchalk")

library(rTensor)
library(gnorm)
library(fastmatrix)
library(rockchalk)

demean.fun = function(x)
{
    M= matrix(0, nrow(x),ncol(x))
    for(i in 1: nrow(x))
    {
        M[i,] = colMeans(x)
    }
    x= x-M
    return(x)
}  

##### data generation ##################################################################
d=3
n=2000
alpha = 1
deg = 3
    Z = matrix(rep(0,d*n),n,d)
    dd = ceiling(d/3)
    Z[,1:dd] = matrix(rgnorm(n*dd, alpha=alpha, beta=deg), n, dd)
    Z[,(dd+1):(2*dd)] = matrix(runif(n*dd, min=-1, max=3), n, dd)
    Z[, (2*dd+1):d] = matrix(rgamma(n*(d-2*dd), shape=3, 
    rate=5), n, (d-2*dd))
    Z= demean.fun(Z)
    
    # standardize Z
    var.Z = (1/n)*t(Z) %*% Z 
    sq.inv = eigen(var.Z)$vectors %*%
    diag(eigen(var.Z)$values^{-1/2}) %*% t(eigen(var.Z)$vectors)
    Z = Z %*% sq.inv

    # generate mixing matrix
    trueW = diag(rep(1, d))
    for (j in 1:(d-1)){
        trueW[j, j+1] = -.5
        trueW[j+1, j] = .5
    }

    X = Z %*% solve(t(trueW))
    ## Whitening 
    var.X = var(X)
    sq.inv = eigen(var.X)$vectors %*%
    diag(eigen(var.X)$values^{-1/2}) %*% t(eigen(var.X)$vectors)
    X = X %*% sq.inv
#####################################################################################

#### sanity check for fourth moment tensor: E x o x o x o x
#### Note that X is n by p matrix. 

#### method 1
 M4.hat = matrix(0, d,d^3)
    for(i in 1:nrow(X))
    {
        xi.1 = X[i,]
        xi.2 = kronecker.prod(xi.1,xi.1)
        xi.3 = kronecker.prod(xi.1,xi.2)
        M4.hat = M4.hat + as.vector(xi.1) %*%  t(xi.3)
    }
    M4.hat = 1/nrow(X)* M4.hat

### method 2
M4.hat2 = rand_tensor(modes = c(d, d, d, d))*0
for(s in 1:nrow(X))
{
 X_s.tensor = rand_tensor(modes = c(d, d, d, d))*0
 for (i in 1:d)
 {
    for(j in 1:d)
    {
        for(k in 1:d)
        {
            for(l in 1:d)
            {
                 X_s.tensor[i,j,k,l] = X[s,i]*X[s,j]*X[s,k]*X[s,l]
            }
        }
    }
 }
 M4.hat2 = M4.hat2 +  X_s.tensor
}
M4.hat2 = 1/nrow(X)*M4.hat2

M4.hat2 = unfold(M4.hat2, row_idx = c(1), col_idx = c(2, 3, 4))
M4.hat2 = as.matrix(M4.hat2@data)

sum(M4.hat-M4.hat2) 
##########################################################################################

#### sanity check for M0(Auddy and Yuan): 

### method 1 

 e2 = rand_tensor(modes = c(d, d, d, d))*0
 for (i in 1:d)
 {
    for(j in 1:d)
    {
        for(k in 1:d)
        {
            for(l in 1:d)
            {
                index = c(i,j,k,l)
                uni.index = unique(index)
                if(  length(uni.index) <=2)
                {
                    e2[i,j,k,l] = 1
                }              
            }
        }
    }
 }

M0.Au = unfold(e2, row_idx = c(1), col_idx = c(2, 3, 4))
M0.Au = M0.Au@data

#### method 2
M_1.234_M0 = matrix(0, d,d^3)

    ## d x d^3, returns 1 for [i,   ( (j-1)*d+k-1)*d+l  ] element.

    ## case 1) i=j=k=l,
    for (i in 1:d)
    {
        M_1.234_M0[i,  (((i-1)*d+i)-1)*d+i  ] =1
        # M_1.234_M0[i,  (((i-1)*d+i)-1)*d+i  ] =3
    }

    ## case 2-1) i=j=k,  i \neq l
       for (i in 1:d)
    {
        for (l in 1:d)
        {
            M_1.234_M0[i,  ( (i-1)*d+i-1)*d+l   ] =1
        }      
    }

    ## case 2-2) i=j=l   i \neq k
       for (i in 1:d)
    {
        for (k in 1:d)
        {
            M_1.234_M0[i,  ( (i-1)*d+k-1)*d+i   ] =1
        }      
    }
    ## case 2-3) i=k=l   i \neq j
       for (i in 1:d)
    {
        for (j in 1:d)
        {
            M_1.234_M0[i,  ( (j-1)*d+i-1)*d+i   ] =1
        }      
    }
    ## case 2-4) j=k=l   i \neq j
       for (i in 1:d)
    {
        for (j in 1:d)
        {
            M_1.234_M0[i,  ( (j-1)*d+j-1)*d+j   ] =1
        }      
    }

    ## case 3-1) i=j, k=l  i \neq k
       for (i in 1:d)
    {
        for (k in 1:d)
        {
            M_1.234_M0[i,  (((i-1)*d+k)-1)*d+k  ] =1
        }      
    }
    ## case 3-2) i=k, j=l  i \neq j
        for (i in 1:d)
    {
        for (j in 1:d)
        {
            M_1.234_M0[i,  (((j-1)*d+i)-1)*d+j  ] =1
        }      
    }
    ## case 3-3) i=l, j=k  i \neq j
        for (i in 1:d)
    {
        for (j in 1:d)
        {
            M_1.234_M0[i,  (((j-1)*d+j)-1)*d+i ] = 1
        }      
    }
M0.Au2 = M_1.234_M0

sum(M0.Au - M0.Au2)

########################################################################################

# Sanity check for M0: M_0_{I,j,k,l} =  E[x_ix_j]E[x_kx_l] + E[x_ix_k]E[x_j x_l] + E[x_ix_l]E[x_jx_k]
# method 1
e = rand_tensor(modes = c(d, d, d, d))*0
 for (i in 1:d)
 {
    for(j in 1:d)
    {
        for(k in 1:d)
        {
            for(l in 1:d)
            {
                index = c(i,j,k,l)
                uni.index = unique(index)
                temp = as.matrix(table(index))
                if(  length(uni.index)==1) ## i=j=k=l
                {
                    e[i,j,k,l] = 3
                }
                else if (temp[1]==temp[2] & length(uni.index)==2) ## {i,j,k,l}= {i1,i2}, i1 \neq i2
                {
                    e[i,j,k,l] = 1   
                }             
            }
        }
    }
 }
M0 = unfold(e, row_idx = c(1), col_idx = c(2, 3, 4))
M0 = M0@data

# method 2
M0.n = matrix(0, d,d^3)

    ## d x d^3, returns 1 for [i,   ( (j-1)*d+k-1)*d+l  ] element.
    ## case 1) i=j=k=l,
        for (i in 1:d)
    {
        M0.n[i,  (((i-1)*d+i)-1)*d+i  ] = 3  
    }

    ## case 2-1) i=j, k=l  i \neq k
       for (i in 1:d)
    {
        for (k in 1:d)
        {
            if(i !=k )
            {
            M0.n[i,  (((i-1)*d+k)-1)*d+k  ] = 1
            }    
        }
    }
    ## case 2-2) i=k, j=l  i \neq j
        for (i in 1:d)
    {
        for (j in 1:d)
        {
            if(i !=j )
            {
            M0.n[i,  (((j-1)*d+i)-1)*d+j  ] = 1
            }
        }
    }
    ## case 2-3) i=l, j=k  i \neq j
        for (i in 1:d)
    {
        for (j in 1:d)
        {   
            if(i !=j )
            {
            M0.n[i,  (((j-1)*d+j)-1)*d+i ] = 1
            }
        }
    }
sum(M0 - M0.n)

##############################################################################################
#### compare M4.hat- M0 and  kurtosis excess from package "rockchalk"

v=c(0.5,1,1)
v= v/norm(v,"2")
v= as.matrix(v)
tv= t(v)

M4.hat.Z = rand_tensor(modes = c(d, d, d, d))*0
for(s in 1:nrow(Z))
{
 Z_s.tensor = rand_tensor(modes = c(d, d, d, d))*0
 for (i in 1:d)
 {
    for(j in 1:d)
    {
        for(k in 1:d)
        {
            for(l in 1:d)
            {
                 Z_s.tensor[i,j,k,l] = Z[s,i]*Z[s,j]*Z[s,k]*Z[s,l]
            }
        }
    }
 }
 M4.hat.Z = M4.hat.Z +  Z_s.tensor
}

M4.hat.Z = 1/nrow(Z)*M4.hat.Z
M0.hat = e
tnsr = M4.hat.Z - M0.hat

my.kurtosis = ttl(tnsr= tnsr, list_mat=list(tv,tv,tv,tv), ms=c(1,2,3,4) )

package.kurtosis = kurtosis(Z %*% v, excess=TRUE )
sum(my.kurtosis@data - package.kurtosis)


# fastICA  w.init=Null then a matrix of normal r.v. is used.

