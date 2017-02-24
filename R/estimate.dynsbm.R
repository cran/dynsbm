# This file is part of dynsbm.

# dysbm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# dynsbm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with dynsbm.  If not, see <http://www.gnu.org/licenses/>
		
estimate.dynsbm <- function (Y, present=NULL, Q,
                             edge.type=c("binary","discrete","continuous"), K=-1,
                             directed=FALSE, self.loop=FALSE,
                             init.method=c("kmeans","specc"),
                             nb.cores=1,
                             perturbation.rate=0., iter.max=20) {
    T=dim(Y)[1]
    N=dim(Y)[2]
    if (is.null(present)){
        present = matrix(1L,dim(Y)[2],dim(Y)[1])
    } else{
        if (storage.mode(present)!= "integer") storage.mode(present) <- "integer"
    }
    if (length(init.method)>1) init.method = init.method[1]
    if(init.method=="specc" && any(present==0L)){ ## specc not available if absent nodes
        init.method <- "kmeans"
    }
    ##---- initialization
    if(init.method=="kmeans"){        
        Yaggr=c()
        if(!any(present==0L)){ ## concat all time step if no absent nodes
            for (t in 1:T){
                if(directed){
                    Yaggr = cbind(Yaggr,Y[t,,],t(Y[t,,]))
                } else{
                    Yaggr = cbind(Yaggr,Y[t,,])
                }
            }
        } else{ ## aggregate all time step if absent nodes
            if(directed){
                Yaggr <- cbind(Y[1,,],t(Y[1,,]))
                for (t in 2:T){
                    Yaggr <- Yaggr+cbind(Y[t,,],t(Y[t,,]))
                }
                
            } else{
                Yaggr <- Y[1,,]
                for (t in 2:T){
                    Yaggr <- Yaggr+Y[t,,]
                }
            }
            Yaggr <- Yaggr / rowSums(present)
        }
        if (nb.cores==1){
            km <- kmeans(Yaggr, Q, nstart=100)
            init.cluster <- km$cluster
        } else{
            RNGkind("L'Ecuyer-CMRG")
            pkm <- mclapply(rep(64/nb.cores,nb.cores), function(nstart) kmeans(Yaggr, Q, nstart=nstart), mc.cores=nb.cores)
            i <- sapply(pkm, function(result) result$tot.withinss)
            init.cluster <- pkm[[which.min(i)]]$cluster
        }
    }
    if(init.method=="specc"){        
        eigenv <- matrix(0, N, Q*T)
        for (t in 1:T){
            D=diag(rowSums(Y[t,,]))
            U=D-Y[t,,]
            xi <- eigen(U,symmetric=TRUE)$vectors[,(ncol(Y[t,,])-Q+1):ncol(Y[t,,])] 
            yi <- xi/sqrt(rowSums(xi^2))
            eigenv[,1+((t-1)*Q):(t*Q-1)] = yi
        }        
        init.cluster = kmeans(eigenv, Q, nstart=100)$cluster
    }        
    ##----- perturbation
    v = sample(1:N, perturbation.rate*N)
    vs = sample(v, perturbation.rate*N)
    init.cluster[v] = init.cluster[vs]
    ##---- estimation
    return(dynsbmcore(T, N, Q, as.vector(Y), present, edge.type, K, init.cluster-1, iter.max, nb.cores, directed, self.loop)) # -1 for compatibility with C++    
}
