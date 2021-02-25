library(vegan)
library(parallel)


metaGamma <- function(S, k) {
    mu <- 1
    rate <- k / mu
    p <- seq(1, 1/S, length = S) - 1/(2 * S)
    
    x <- qgamma(p, k, rate)
    
    return(x / sum(x))
}


param <- expand.grid(S = round(exp(seq(log(500), log(5000), length.out = 3))), 
                     k = exp(seq(-2, 5, length.out = 7)), 
                     size = c(1000, 10000))

sim <- mclapply(1:nrow(param), mc.cores = 10, FUN = function(i) {
    # otu table with `nsite` sites; note: sites are really just replicates
    # of estimator becaues estimator is calculated for each site
    nsite <- 1000
    x <- t(rmultinom(nsite, param$size[i], metaGamma(param$S[i], param$k[i])))
    
    # Chao richness estimates
    r <- estimateR(x)[2:3, ]
    
    # Chao lower CI
    rlo <- r[1, ] - 2 * r[2, ]
    
    # Chao upper CI
    rhi <- r[1, ] + 2 * r[2, ]
    
    # output
    o <- quantile(r[1, ], probs = c(0.025, 0.5, 0.975))
    names(o) <- paste0('q', c('lo', 'mid', 'hi'))
    
    # how often does richness estimate CI contain true richness
    o <- c(o, propInclude = mean(param$S[i] >= rlo & param$S[i] <= rhi))
    
    # how many species missed in samples
    o <- c(o, propSamp = mean(rowSums(x > 0) / param$S[i]))
    
    return(o)
})

sim <- as.data.frame(do.call(rbind, sim))

plot(param$k, sim$propInclude, log = 'x')


allS <- sort(unique(param$S))
allSize <- sort(unique(param$size))

layout(matrix(1:(length(allS) * length(allSize)), ncol = length(allSize)))
par(mar = c(1, 1, 0, 0), oma = c(3, 3, 1, 1))

for(n in allSize) {
    for(s in rev(allS)) {
        i <- param$S == s & param$size == n
        plot(param$k[i], sim$propInclude[i], xlim = range(param$k), ylim = 0:1, 
             axes = FALSE, frame.plot = TRUE, log = 'x')
        
        if(n == min(allSize)) axis(2)
        if(s == min(allS)) axis(1)
    }
}

