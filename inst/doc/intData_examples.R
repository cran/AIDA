## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----results='hide', message = FALSE, warning=FALSE---------------------------
library(AIDA)

## -----------------------------------------------------------------------------
data(creditcard)
CreditCard_microdata <- creditcard$microdata
CreditCard_min_max <- creditcard$min_max
CreditCard_CR <- creditcard$centers_ranges

## -----------------------------------------------------------------------------
credit_card_int_unif <- intData(CreditCard_min_max, Seq = "LbUb_VarbyVar", 
                                VarNames = colnames(CreditCard_microdata)[3:7])

# Check the parameters of the latent distribution                          
credit_card_int_unif@LatentParam

## ----fig.width=7, fig.height=2------------------------------------------------
credit_agrby <- factor(paste(CreditCard_microdata$Name,CreditCard_microdata$Month, sep = "_"))
credit_card_U <- get_latent_var(CreditCard_microdata[,3:7], CreditCard_min_max, credit_agrby,
                                rownames(CreditCard_min_max), Seq = "LbUb_VarbyVar")

oldpar <- par(no.readonly = TRUE)
par(mfrow=c(1,5), mar=c(2, 2, 2, 1))
for (i in 1:5){
    hist(credit_card_U[,i], xlab = NULL, ylab = NULL, 
            main = colnames(credit_card_U)[i], probability = TRUE)
    lines(density(credit_card_U[,i], na.rm = TRUE), col = '#009de0', lwd = 2)
}
par(oldpar)

## -----------------------------------------------------------------------------
credit_card_int_triang <- intData(CreditCard_min_max, Seq = "LbUb_VarbyVar", LatentDist = "Triang", 
                                    VarNames = colnames(CreditCard_microdata)[3:7])

# Check the parameters of the latent distribution
credit_card_int_triang@LatentParam

head(credit_card_int_triang)

## -----------------------------------------------------------------------------
credit_card_int_triang@Centers[1:5,]
credit_card_int_triang@Ranges[1:5,]
LowerBounds(credit_card_int_triang)[1:5,]
UpperBounds(credit_card_int_triang)[1:5,]

## -----------------------------------------------------------------------------
credit_card_int_KDE <- intData(CreditCard_CR, Seq = "LbUb_VarbyVar", 
                                VarNames = colnames(CreditCard_microdata)[3:7], 
                                LatentCase = "General", LatentDist = "KDE", Umicro = credit_card_U)

# Check the parameters of the latent distribution
credit_card_int_KDE@LatentParam

## -----------------------------------------------------------------------------
credit_card_int_agr <- micro2intData(CreditCard_microdata[,3:7], credit_agrby, LatentCase = "General")

# Check the parameters of the latent distribution
credit_card_int_agr@LatentParam

head(credit_card_int_agr)

## -----------------------------------------------------------------------------
credit_card_cov <- int_cov(credit_card_int_agr)
credit_card_cor <- cov2cor(credit_card_cov)

# Check the covariance matrix
credit_card_cov

## ----fig.width=6, fig.height=6------------------------------------------------
SYMB.pairs.panels(credit_card_int_agr, type = "rectangles", 
                    corr = credit_card_cor, labels = colnames(credit_card_int_agr))

