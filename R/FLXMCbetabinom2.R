# Based on Ivan Kondofersky's code from his bachelor's thesis.

#' @export
FLXMCbetabinom2 = local({

lbeta1 <- function(x, size, a, b) {
    #unique(cbind(x, size-x))
    
    s <- seq(from=0, to=size, by=1)
    uniquelb <- lbeta(a+s, b+size-s)

    res1 <- uniquelb[x+1L]
    #res2 <- lbeta(a+x, b+size-x)

    #if(!all(res1==res2)) browser() else cat("ok\n")

    #sort(unique(res2))
    #sort(unique(uniquelb))

    res1
}

digamma1 <- function(x, a, size) {
    uniquedg <- digamma(seq(from=0, to=size, by=1) + a)
    res1 <- uniquedg[x+1L]
    #res2 <- digamma(x+a)

    #if(!all(res1==res2)) browser()

    # unique(abs(res1-res2))
    #sort(unique(digamma(x+px)))
    res1
}

dbetabinom <- function(x, size, a, b, log=FALSE) {
    #z <- lchoose(size, x) + lbeta(a+x, b+size-x) - lbeta(a,b)
    #z <- lbeta(a+x, b+size-x) - lbeta(a,b)
    z <- lbeta1(x, size, a, b) - lbeta(a,b)
    if(!any(is.finite(z))) {
        cat("z not finite\n")
        return(NA)
    }
    if(log) z else exp(z)
}

BBlogLikGrad <- function(ab, x, size, w=1) {
    c(sum(w*(digamma1(x, ab[1], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[1]) + digamma(ab[1]+ab[2]))),
             
      sum(w*(digamma1(size-x, ab[2], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[2]) + digamma(ab[1]+ab[2]))))
}

BBlogLikGrad2 <- function(ab, x, size, w=rep(1, length(x))) {
    res <- .C("bb_loglik_grad", ab, x, length(x), as.integer(size), w, double(2))
    res[[6]]
}


## the old code
if(FALSE) {
    dbetabinom <- function(x, size, a, b, log=FALSE) {
        z <- lchoose(size, x) + lbeta(a+x, b+size-x) - lbeta(a,b)
        if(log) z else exp(z)
    }

    BBlogLikGrad <- function(ab, x, size, w=1) {
        c(sum(w*(digamma(x+ab[1]) -
                digamma(size+ab[1]+ab[2]) -
                digamma(ab[1]) + digamma(ab[1]+ab[2]))),
                
        sum(w*(digamma(size-x+ab[2]) -
                digamma(size+ab[1]+ab[2]) -
                digamma(ab[2]) + digamma(ab[1]+ab[2]))))
    }
}



BBlogLik <- function(ab, x, size, w=1) {
    sum(w * dbetabinom(x, size, ab[1], ab[2], log=TRUE))
}


BBmle <- function(x, size=NULL, w=1, eps=sqrt(.Machine$double.eps)) {
    if(!all(x>=0)) stop("Negative values in data matrix!")
    N <- ncol(x)

    if(is.null(size)) {
        size <- apply(x, 2, max)
    } else {
        size <- rep(size, length.out=N)
    }

    res <- matrix(NA, nrow=2, ncol=N)
    for(i in 1:N) {
        res[,i] <- optim(c(1,1), fn=BBlogLik, gr=BBlogLikGrad,
                         x=x[,i], size=size[i], w=w, 
                         control=list(fnscale=-1),
                         method="L-BFGS-B", lower=c(eps, eps))$par
    }

    rownames(res) <- c("alpha", "beta")
    colnames(res) <- colnames(x)
    res
}


BBlogLikReg <- function(ab, x, size, w=1, alpha2=0) {
    a2 = alpha2/length(x)
    # bad
    if(FALSE) {
        sum(w * dbetabinom(x, size, ab[1], ab[2], log=TRUE)) +
            sum(a2 * dbetabinom(x, size, ab[1], ab[2], log=TRUE))
    }

    # slightly less bad
    if(FALSE) {
        dens = dbetabinom(x, size, ab[1], ab[2], log=TRUE)
        sum(w * dens) + sum(a2 * dens)
    }

    # bit better
    dens = dbetabinom(x, size, ab[1], ab[2], log=TRUE)
    sum((w+a2)*dens)
}

BBlogLikGradReg <- function(ab, x, size, w=1, alpha2=0) {
    grad1 = (digamma1(x, ab[1], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[1]) + digamma(ab[1]+ab[2]))
    grad2 = (digamma1(size-x, ab[2], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[2]) + digamma(ab[1]+ab[2]))

    a2 = alpha2/length(x)

    if(FALSE) {
        c(sum(w*grad1) + sum(a2*grad1),
          sum(w*grad2) + sum(a2*grad2))
    }

    c(sum((w+a2)*grad1), sum((w+a2)*grad2))
}


BBmle5 <- function(x, size=NULL, w=1, eps=sqrt(.Machine$double.eps)) {
    if(!all(x>=0)) stop("Negative values in data matrix!")
    N <- ncol(x)

    browser()

    if(is.null(size)) {
        size <- apply(x, 2, max)
    } else {
        size <- rep(size, length.out=N)
    }

    res <- matrix(NA, nrow=2, ncol=N)
    for(i in 1:N) {

        abstar = optim(c(1,1), fn=BBlogLik, gr=BBlogLikGrad,
                       x=x[,i], size=size[i], w=1, 
                       control=list(fnscale=-1),
                       method="L-BFGS-B", lower=c(eps, eps))$par

        # which x fits those params best?
        if(FALSE) {
            cand = dbetabinom(seq(0,1,by=0.01), size[i], margin_i[1], margin_i[2])
            which.max(cand)
        }

        #cand = dbetabinom(seq(0,1,by=0.01), size[i], margin_i[1], margin_i[2])

        xs = seq(0,size[i],by=0.05)
        cand = dbetabinom(xs, size[i], margin_i[1], margin_i[2])
        plot(xs, cand, type="o")

        cand[which.max(cand)]

        wold = w
        w = wold

        w = runif(length(w))

        res_reg = optim(c(1,1), fn=BBlogLikReg, gr=BBlogLikGradReg,
                         x=x[,i], size=size[i], w=w, 
                         alpha2=0,
                         control=list(fnscale=-1),
                         method="L-BFGS-B", lower=c(eps, eps))$par


        res[,i] <- optim(c(1,1), fn=BBlogLik, gr=BBlogLikGrad,
                         x=x[,i], size=size[i], w=w, 
                         control=list(fnscale=-1),
                         method="L-BFGS-B", lower=c(eps, eps))$par

        res_reg
        res[,i]

        abstar

        BBlogLik(res_reg, x[,i], size[i], w=w)
        BBlogLik(res[,i], x[,i], size[i], w=w)

        if(FALSE) {

            regs = lapply(0:100, \(a2) {
                optim(c(1,1), fn=BBlogLikReg,
                        x=x[,i], size=size[i], w=w, 
                        alpha2=a2,
                        control=list(fnscale=-1),
                        method="L-BFGS-B", lower=c(eps, eps))$par
            }) %>% do.call(rbind, .)


            paramdf = data.table(a=regs[,1], b=regs[,2], col=1) %>%
                rbind(data.table(a=abstar[1], b=abstar[2], col=2))


            plot(paramdf$a, paramdf$b, type="p", col=paramdf$col)



            loglikdf = data.table(ll = apply(regs, 1, BBlogLik, x[,i], size[i], w=w), col=1)# %>%
                #rbind(data.table(ll=BBlogLik(abstar, x[,i], size[i]), col=2))


            plot(loglikdf$ll, type="l", col=loglikdf$col)


        }

    }

    rownames(res) <- c("alpha", "beta")
    colnames(res) <- colnames(x)
    res
}


BBmle2 <- function(x, size, w=1, eps=sqrt(.Machine$double.eps)) {
    #if(!all(x>=0)) stop("Negative values in data matrix!")
    N <- ncol(x)
    size <- rep(size, length.out=N)

    res <- matrix(.Call("rif_bb_mle", x, size, w), nrow=2)
    rownames(res) <- c("alpha", "beta")
    colnames(res) <- colnames(x)

    res
}


BBlogLik_plaid <- function(ab, xuni, ws, size) {
    ll = dbetabinom(xuni, size[1], ab[1], ab[2], log=TRUE)
    sum(ws * ll)
}

BBlogLikGrad_plaid <- function(ab, xuni, ws, size) {

    if(FALSE) {

        x = x[,i]


        res =
        c(sum(w*(digamma(x+ab[1]) -
                digamma(size+ab[1]+ab[2]) -
                digamma(ab[1]) + digamma(ab[1]+ab[2]))),
        sum(w*(digamma(size-x+ab[2]) -
                digamma(size+ab[1]+ab[2]) -
                digamma(ab[2]) + digamma(ab[1]+ab[2]))))


    }


    const1 = -digamma(size+ab[1]+ab[2]) -
             digamma(ab[1]) + digamma(ab[1]+ab[2])

    const2 = -digamma(size+ab[1]+ab[2]) -
             digamma(ab[2]) + digamma(ab[1]+ab[2])

    s1 = sum(ws * (digamma(xuni+ab[1]) + const1))
    s2 = sum(ws * (digamma(size-xuni+ab[2]) + const2))


    c(s1,s2)
}


BBmle3 <- function(x, size=NULL, w=1, eps=sqrt(.Machine$double.eps)) {
    if(!all(x>=0)) stop("Negative values in data matrix!")
    N <- ncol(x)

    if(is.null(size)) {
        size <- apply(x, 2, max)
    } else {
        size <- rep(size, length.out=N)
    }

    wsum = lapply(seq(ncol(x)), \(col) {
        tapply(w, x[,col], sum)
    })



    ## plaid
    resmat = lapply(1:N, \(i) {
        xuni = unique(x[,i]) %>%
            sort

        res = optim(c(1,1), fn=BBlogLik_plaid, gr=BBlogLikGrad_plaid,
                    xuni = xuni, ws = wsum[[i]], size = size[i],
                    control=list(fnscale=-1),
                    method="L-BFGS-B", lower=c(eps, eps))$par
        res
    }) %>% do.call(cbind, .)

    rownames(resmat) <- c("alpha", "beta")
    colnames(resmat) <- colnames(x)

    resmat
}

BBmle4 <- function(x, xuni, size=NULL, w=1, eps=sqrt(.Machine$double.eps)) {
    #browser()

    #if(!all(x>=0)) stop("Negative values in data matrix!")
    N <- ncol(x)
    size <- rep(size, length.out=N)

    wsfun = \() {
        lapply(seq(ncol(x)), \(col) {
            tapply(w, x[,col], sum)
        })
    }
    #wsum = wsfun()
    wsum=NULL


    optifun = \() .Call("rif_bb_mle_plaid", x, xuni, size, w, wsum)

    tmp = try(optifun())
    #stop("asdf")
    if(is(tmp, "try-error")) browser()

    res <- matrix(tmp, nrow=2)
    rownames(res) <- c("alpha", "beta")
    colnames(res) <- colnames(x)

    return(res)
}


function(formula=.~., size, eps=1e-5, mle="BBmle4") {

    z <- new("FLXMC", weighted=TRUE, formula=formula, dist="mvbetabinom",
             name="model based beta-binomial clustering")

    size <- as.integer(size)
    

    if(is.character(mle)) {
        mlename = mle
        mle = get(mle, mode="function")
    }

    z@preproc.y = \(y) {
        xuni <<- lapply(seq(ncol(y)), \(col) {
            sort(unique(y[,col]))
        })
        y
    }


    z@defineComponent <- expression({
        logLik <- function(x,y) {
            z <- y
            for(k in 1:ncol(y)) {
                z[,k] <- dbetabinom(y[,k], size=size[k], a=ab[1,k], b=ab[2,k], 
                                    log=TRUE)
            }
            rowSums(z, na.rm=TRUE)
        }
        predict <- function(x, ...) {
            matrix(center, nrow=nrow(x), ncol=length(center), byrow=TRUE)
        }

        new("FLXcomponent", 
            parameters=list(a=ab[1,], b=ab[2,], size=size, prob=prob, center=center),
            df=length(ab), logLik=logLik, predict=predict)
    })

    z@fit <- function(x,y,w) {
        para <- list(size=rep(size, length=ncol(y)))
        if(mlename == "BBmle4") {
            para$ab <- mle(y, xuni, size=para$size, w=w, eps=eps)
        } else {
            para$ab <- mle(y, size=para$size, w=w, eps=eps)
        }
        para$prob <- apply(para$ab, 2, function(z) z[1]/sum(z))
        para$center <- para$prob * para$size
        with(para, eval(z@defineComponent))
    }

    z
}

})


