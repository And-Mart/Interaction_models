#### Script with the two models present in the paper "Models that estimate missing interactions change our understanding of the structure and species roles in the largest seed-dispersal network" ####

#### Packages required ####
library(ape) #phylogeny manipulation
library(igraph) #network and matrix manipulation
library(sparsesvd) #singular value decomposition function
library(Rphylopars) #phylogenetic inference function
library(dplyr) #general data manipulation
library(castor) #used for the function to find the nearest tip in the phylogeny


#### Files required ####
Birds_tree <- read.tree("birds_phylo.txt") # Birds phylogeny with all species we are interested + spp used as input
Plants_tree <- read.tree("plants_phylo.txt") # Plant phylogeny with all species we are interested + spp used as input
In_m <- read.table("input_matrix.txt", row.names = 1, header = T) #Input matrix with selected species
atl_i<- read.table("atlantics_info.txt", header = T) #Trait and interaction frequency for each spp pair obtained from the Atlantics:frugivory database


##############################
#### Trait matching model ####
##############################

## Trait model - function
D.Int<-function(x=c(4,0.01,1), Data=NULL, ret.list = FALSE)
{
  #x=c(4,0.04)
  #free parameters
  alpha <- x[1]
  S <- x[2]
  phi <- x[3] 
  
  A <- Data$bird.gape[1]
  Y <- Data$freq #proportion of studies with documented interactions  
  
  #Y=(Y-min(Y))/(max(Y)-min(Y))+0.01
  B <- Data$seed.size
  
  #trait-matching function
  temp <- (S*(A-B))+exp(-phi)
  den <- sqrt(2*pi)*temp
  P <- (alpha*(exp(-log(temp)^2))/2)/den
  
  J <- data.frame(seed.size = Data$seed.size,Freq = Y, predicted = P)
  Z <- ((J$Freq-J$pred)^2) #computing least squares
  Z.w <- Z*Data$n.studies #weighting distances by the number of studies in which the pair has been recorded
  SSQ <- (sum (Z.w, na.rm=T)) #summing least squares
  
  if (ret.list == TRUE) {list(J,SSQ)
  }else{
    SSQ
  }
}


# Family information used to select which spp to use for parameter optimization
fam_inf<-distinct(atl_i[,c(1,8)]) # bird spp and families in the atlantics dataframe
in_fam<-fam_inf[which(fam_inf$bird_sp %in% colnames(In_m)),]

## Creating a matrix to be filled with the estimated interactions
dat<-matrix(0,length(unique(atl_i$plant_sp)),length(unique(atl_i$bird_sp)))
row.names(dat) <- unique(atl_i$plant_sp)
colnames(dat) <- unique(atl_i$bird_sp)

###### Finding the interaction frequency for each bird spp
for(k in colnames(dat)){
  sub1 <- filter(atl_i, bird_sp == k)#Dataframe to save new freq
  if((k %in% in_fam$bird_sp)==T){ # If spp is in the input, use just its information
    op <- optim(c(3,0.1,1), D.Int, Data = sub1, lower = c(0,0.0001,1), upper = c(4,2,5))
    
    
  }
  if((k %in% in_fam$bird_sp)==F){ #if spp isnt in input, use the nearest spp in the input for optimization
    fam_sp<-fam_inf[which(fam_inf$family == fam_inf[which(fam_inf$bird_sp == k),2]),1]
    cut_t<-keep.tip(Birds_tree, fam_sp)
    near <- find_nearest_tips (cut_t, target_tips = k)
    
    v.near <- order(near$nearest_distance_per_tip)
    ord_n<-cut_t$tip.label[v.near]
    opt_sp<-ord_n[which(ord_n %in% in_fam$bird_sp)]
    op_in<-filter(atl_i, bird_sp==opt_sp[1])
    op <- optim(c(3,0.1,1), D.Int, Data = op_in, lower = c(0,0.0001,1), upper = c(4,2,5))
    
  }
  
  res <- D.Int(x = op$par,Data = sub1, ret.list = TRUE)[[1]] #data frame from function with the parameter set with the best fit
  prob.int <- res[3]
  
  sub3 <- sub1
  sub3$freq<- prob.int$predicted #replacing frequencies from target with that predicted by the trait matching model
  
  dat[,which(colnames(dat)==k)]<-sub3$freq
}

dat[which(is.na(dat))]<-0 # changing NAs due to missing trait data to 0
write.table(dat, "TM_result.txt", col.names = T, row.names = T) #saving the resulting matrix

##############################
#### Latent traits model #####
##############################
#### Preparing the objects we need to run the model ####
In_m <- as.matrix(In_m)

In_b <- colnames(In_m)[which(colnames(In_m) %in% Birds_tree$tip.label)] #Input birds
In_p <- row.names(In_m) #Input plants

Out_b <- Birds_tree$tip.label[-which(Birds_tree$tip.label %in% In_b)] #Birds we want to estimate
Out_p <- Plants_tree$tip.label[-which(Plants_tree$tip.label %in% In_p)] #Plants we want to estimate 

## Input matrix dimensions for future use
mb <- nrow(In_m)
nb <- ncol(In_m)

## Getting the incidence matrix and turning it into a adjacency matrix
Mb.inc <- graph.incidence(In_m, weighted = NULL)
Mb.adj <- as_adjacency_matrix(Mb.inc, sparse=T)

#---------------------------------------------------------
#Latent trait model step 1: obtaining latent traits using t-SVD
#---------------------------------------------------------
# In this step we use singular value decomposition to summarize the information inside the input interaction matrix.

## Testing the function for our matrix
SVD.res.b <- sparsesvd(Mb.adj, rank=0) ## Birds
cumsum(SVD.res.b$d/sum(SVD.res.b$d)) #proportion of structure explained by each rank

# the sum of the matrix products (u %*% d %*% v)  of all ranks recovers M. We need to set the number of ranks we gonna use
# We opted to choose the first rank where the sum becomes equal to 0.75
rank_b <- 27

m1_b <- nrow(SVD.res.b$v)
n1_b <- ncol(SVD.res.b$v)
SVD.res.b$v <- SVD.res.b$v+matrix(rnorm(m1_b*n1_b,0,0.01),m1_b,n1_b)

#-------------------------------------------------------------------------------------------------
#Latent trait model step 2: map latent traits into phylogeny and impute missing traits based on distance
#inferring the values of the left and right subspaces 
#-------------------------------------------------------------------------------------------------

## In this step we use two of the dimensions obtained with the SVD as traits of each species and use phylogenetic imputation
# to estimate values for the species that weren't in our input matrix

# Birds imputation 
traits.new.b <- list() # just creating a list to save the results from the imputation
traits.new.pb <- list()
for(r in 1:rank_b){
  traits_b <- matrix(c(SVD.res.b$u[,r],SVD.res.b$v[,r]), ncol = 2, 
                     nrow = nrow(Mb.adj))
  plant.latent_b <- traits_b[1:length(In_p),] #only plants latent traits
  p.phylopar_b <- data.frame(species=In_p, V1=traits_b[1:length(In_p),1],
                             V2=traits_b[1:length(In_p),2]) #phylopar input format
  birds.latent_b <- traits_b[-(1:length(In_p)),] #only birds latent traits
  f.phylopar_b <- data.frame(species=In_b, V1=traits_b[-(1:length(In_p)),1],
                             V2=traits_b[-(1:length(In_p)),2]) #phylopar input format
  
  f_phy_b <- phylopars(trait_data = f.phylopar_b,
                       tree = Birds_tree, pheno_error = FALSE,
                       phylo_correlated = TRUE,pheno_correlated = TRUE)
  p_phy_b <- phylopars(trait_data = p.phylopar_b,
                       tree = Plants_tree, pheno_error = FALSE,
                       phylo_correlated = TRUE,pheno_correlated = TRUE)
  traits.new.b[[r]] <- f_phy_b$anc_recon[1:length(Birds_tree$tip.label),]
  traits.new.pb[[r]] <- p_phy_b$anc_recon[1:length(Plants_tree$tip.label),]
}

#----------------------------------------------------------------------------------------
#Latent model step 3: Build new matrix with missing interactions inferred from phylogeny
#----------------------------------------------------------------------------------------
# The last step is constructing a interaction matrix with the other species in the phylogenies (the ones that were not in the
# input matrix). After this step is done we can combine both matrix again to have the full frugivores x plants matrix

#using the imputed latent traits to generate a new interaction matrix
new.dim.b <- nrow(Mb.adj)+length(Out_b)+length(Out_p) # matrix dimensions

b_sp <- c(In_p, In_b, Out_p, Out_b) # names for the rows/columns

## constructing the new matrix
Mb.part2 <- array(NA,c(new.dim.b,new.dim.b,rank_b))
for(r in 1:rank_b){
  
  traits.all.b <- rbind(traits.new.b[[r]],traits.new.pb[[r]]) #animal + plant traits
  or.d.b<-match(b_sp, row.names(traits.all.b)) #matching the names order
  traits.all.b <-traits.all.b[or.d.b,] 
  
  temp.b <- traits.all.b[,1] %*% diag(SVD.res.b$d[r],new.dim.b,new.dim.b)
  Mb.part2[,,r] <- t(temp.b) %*% t(traits.all.b[,2]) #
}

Mb.new <- round(rowSums(Mb.part2, dims = 2,na.rm = T),1)
Mb.new[Mb.new < 0] <- 0 #correcting negative values
colnames(Mb.new) <- rownames(traits.all.b)
rownames(Mb.new) <- rownames(traits.all.b)

M.new.b<-Mb.new[match(c(In_p, Out_p), row.names(Mb.new)), 
                match(c(In_b, Out_b), colnames(Mb.new))] # the resulting interaction matrix

write.table(M.new.b, "LT_result.txt", col.names = T, row.names = T) #saving the resulting matrix 