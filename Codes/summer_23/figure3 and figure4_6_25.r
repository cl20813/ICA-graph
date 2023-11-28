## figure 3
## run "library and function_6_25.r" file first!

## data generation
  d=5; n=2000
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

## 
myW = Auddy_fastICA_mod2(X,unwh_trueW,10,400) # my M4hat and my deflation
unwh_myW = myW %*% Whitening_mat

melted_trueW = melt(trueW)
melted_unwh_myW = melt(unwh_myW )
data_trueW = data.frame(melted_trueW, 
category= as.factor(rep("trueW",nrow(melted_trueW))),
group = as.factor(rep("trueW",nrow(melted_trueW))))

data_myW = data.frame(melted_unwh_myW,
 category= as.factor(rep("estimated W",nrow(melted_trueW))),
 group = as.factor(rep("trueW",nrow(melted_trueW))))
data = rbind(data_trueW, data_myW)

ggplot(data, aes(x=Var2, y=Var1, fill = value)) + 
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # lables vertical
        strip.text.y = element_blank()) +  #remove facet bar on y 
  scale_fill_gradient(low = "darkblue", high = "lightblue") +

  ggtitle("d=25, n=4000") +

  facet_grid(rows = vars(data$group), 
             cols = vars(data$category), scales = "free", space="free_y") #facets to add gaps 


## figure 4 

  myW = Auddy_fastICA_mod2(X,unwh_trueW,10,400) # my M4hat and my deflation

  unwh_myW = myW %*% Whitening_mat
  sparse_unwh_myW_01 = hard_thre(unwh_myW ,0.1)
  sparse_unwh_myW_02 = hard_thre(unwh_myW ,0.2)    
  sparse_unwh_myW_03 = hard_thre(unwh_myW ,0.3) 
  sparse_unwh_myW_04 = hard_thre(unwh_myW ,0.4) 
  sparse_unwh_myW_05 = hard_thre(unwh_myW ,0.5) 
 
  g0 = ICA_graph(trueW)
  g1 = ICA_graph(sparse_unwh_myW_01)
  g2 = ICA_graph(sparse_unwh_myW_02)
  g3 = ICA_graph(sparse_unwh_myW_03)
  g4 = ICA_graph(sparse_unwh_myW_04)
  g5 = ICA_graph(sparse_unwh_myW_05)

melted_trueW = melt(g0)
melted_g1 = melt(g1 )
melted_g2 = melt(g2 )
melted_g3 = melt(g3 )
melted_g4 = melt(g4 )
melted_g5 = melt(g5 )

data_trueW = data.frame(melted_trueW, 
category= as.factor(rep("trueW",nrow(melted_trueW))),
group = as.factor(rep(1,nrow(melted_trueW))))

data_myW1 = data.frame(melted_g1,
 category= as.factor(rep("th=0.1",nrow(melted_trueW))),
 group = as.factor(rep(1,nrow(melted_trueW))))

data_myW2 = data.frame(melted_g2,
 category= as.factor(rep("th=0.2",nrow(melted_trueW))),
 group = as.factor(rep(1,nrow(melted_trueW))))

data_myW3 = data.frame(melted_g3,
 category= as.factor(rep("th=0.3",nrow(melted_trueW))),
 group = as.factor(rep(1,nrow(melted_trueW))))

data_myW4 = data.frame(melted_g4,
 category= as.factor(rep("th=0.4",nrow(melted_trueW))),
 group = as.factor(rep(1,nrow(melted_trueW))))

data_myW5 = data.frame(melted_g5,
 category= as.factor(rep("th=0.5",nrow(melted_trueW))),
 group = as.factor(rep(1,nrow(melted_trueW))))

data = rbind(data_trueW, data_myW1, data_myW2,data_myW3,data_myW4,data_myW5)

ggplot(data, aes(x=Var2, y=Var1, fill = value)) + 
  geom_raster() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), # lables vertical
        strip.text.y = element_blank()) +  #remove facet bar on y 
  scale_fill_gradient(low = "darkblue", high = "lightblue") +

  ggtitle("ICA graph, d=25, n=4000") +
  facet_grid(rows = vars(data$group), 
             cols = vars(data$category), scales = "free", space="free_y") #facets to add gaps 
