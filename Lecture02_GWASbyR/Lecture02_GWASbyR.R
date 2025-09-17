# To perform GWAS with correlation
# Input: y-vector of phenotype, X matrix of numeric genotype with rows as individuals and SNPs as column
# Output: vector of probability for SNPs in same order
GWASbyCor <- function(X, y) {
  n <- nrow(X)
  r <- cor(y, X)
  n <- nrow(X)
  t <- r / sqrt((1 - r^2) / (n - 2))
  p <- 2 * (1 - pt(abs(t), n - 2))
  zeros <- p == 0
  p[zeros] <- 1e-10
  return(p)
}


setwd("/Users/Jiabo/Documents/China/SWUN/Conference/2025西安培训/lecture/Code_Materials")
source("gapit_functions.txt")

myGD <- read.table("/Users/sameen/workspace/statistical genomics/data/mdp_numeric.txt", head = T)
myGM <- read.table("/Users/sameen/workspace/statistical genomics/data/mdp_SNP_information.txt", head = T)

source("/Users/sameen/workspace/statistical genomics/function/G2P.R")
source("/Users/sameen/workspace/statistical genomics/function/GWASbyCor.R")
X <- myGD[, -1]
index1to5 <- myGM[, 2] < 6
X1to5 <- X[, index1to5]
set.seed(99164)
mySim <- G2P(X = X1to5, h2 = .75, alpha = 1, NQTN = 10, distribution = "norm")
p <- GWASbyCor(X = X, y = mySim$y)

color.vector <- rep(c("deepskyblue", "orange", "forestgreen", "indianred3"), 10)
m <- nrow(myGM)
plot(t(-log10(p)) ~ seq(1:m), col = color.vector[myGM[, 2]])
abline(v = mySim$QTN.position, lty = 2, lwd = 2, col = "black")

p.obs <- p[!index1to5]
m2 <- length(p.obs)
p.uni <- runif(m2, 0, 1)
order.obs <- order(p.obs)
order.uni <- order(p.uni)

plot(
  -log10(p.uni[order.uni]),
  -log10(p.obs[order.obs])
)
abline(a = 0, b = 1, col = "red")

order.obs <- order(p.obs)
X6to10 <- X[, !index1to5]
Xtop <- X6to10[, order.obs[1]]

boxplot(mySim$y ~ Xtop)

PCA <- prcomp(X)
plot(mySim$y, PCA$x[, 2])
cor(mySim$y, PCA$x[, 2])

set.seed(99164)
s <- sample(length(mySim$y), 10)
plot(mySim$y[s], PCA$x[s, 2])
cor(mySim$y[s], PCA$x[s, 2])

y <- mySim$y
X <- cbind(1, PCA$x[, 2], Xtop)
LHS <- t(X) %*% X
C <- solve(LHS)
RHS <- t(X) %*% y
b <- C %*% RHS
yb <- X %*% b
e <- y - yb
n <- length(y)
ve <- sum(e^2) / (n - 1)
vt <- C * ve
t <- b / sqrt(diag(vt))
p <- 2 * (1 - pt(abs(t), n - 2))

LM <- cbind(b, t, sqrt(diag(vt)), p)
rownames(LM) <- cbind("Mean", "PC2", "Xtop")
colnames(LM) <- cbind("b", "t", "SD", "p")
LM

G <- myGD[, -1]
n <- nrow(G)
m <- ncol(G)
P <- matrix(NA, 1, m)

for (i in 1:m) {
  x <- G[, i]
  if (max(x) == min(x)) {
    p <- 1
  } else {
    X <- cbind(1, PCA$x[, 2], x)
    LHS <- t(X) %*% X
    C <- solve(LHS)
    RHS <- t(X) %*% y
    b <- C %*% RHS
    yb <- X %*% b
    e <- y - yb
    n <- length(y)
    ve <- sum(e^2) / (n - 1)
    vt <- C * ve
    t <- b / sqrt(diag(vt))
    p <- 2 * (1 - pt(abs(t), n - 2))
  } # end of testing variation
  P[i] <- p[length(p)]
} # end of looping for markers

p.obs <- P[!index1to5]
m2 <- length(p.obs)
p.uni <- runif(m2, 0, 1)
order.obs <- order(p.obs)
order.uni <- order(p.uni)

plot(-log10(p.uni[order.uni]),
  -log10(p.obs[order.obs]),
  ylim = c(0, 7)
)
abline(a = 0, b = 1, col = "red")

G <- myGD[, -1]
n <- nrow(G)
m <- ncol(G)
P <- matrix(NA, 1, m)

for (i in 1:m) {
  x <- G[, i]
  if (max(x) == min(x)) {
    p <- 1
  } else {
    X <- cbind(1, PCA$x[, 1:3], x)
    LHS <- t(X) %*% X
    C <- solve(LHS)
    RHS <- t(X) %*% y
    b <- C %*% RHS
    yb <- X %*% b
    e <- y - yb
    n <- length(y)
    ve <- sum(e^2) / (n - 1)
    vt <- C * ve
    t <- b / sqrt(diag(vt))
    p <- 2 * (1 - pt(abs(t), n - 2))
  } # end of testing variation
  P[i] <- p[length(p)]
} # end of looping for markers

p.obs <- P[!index1to5]
m2 <- length(p.obs)
p.uni <- runif(m2, 0, 1)
order.obs <- order(p.obs)
order.uni <- order(p.uni)

plot(-log10(p.uni[order.uni]),
  -log10(p.obs[order.obs]),
  ylim = c(0, 7)
)
abline(a = 0, b = 1, col = "red")

p.obs <- P
m2 <- length(p.obs)
p.uni <- runif(m2, 0, 1)
order.obs <- order(p.obs)
order.uni <- order(p.uni)

plot(
  -log10(p.uni[order.uni]),
  -log10(p.obs[order.obs]),
)
abline(a = 0, b = 1, col = "red")

color.vector <- rep(c("deepskyblue", "orange", "forestgreen", "indianred3"), 10)
m <- nrow(myGM)
plot(t(-log10(P)) ~ seq(1:m), col = color.vector[myGM[, 2]])
abline(v = mySim$QTN.position, lty = 2, lwd = 2, col = "black")

myGD <- read.table(file = "mdp_numeric.txt", head = T)
X <- myGD[, -1]
p <- colMeans(X) / 2
M <- X - 1
Z <- t(M) - 2 * (p - .5)
K <- crossprod((Z), (Z))
adj <- 2 * sum(p * (1 - p))
K <- K / adj

# source("/Users/Jiabo/Dropbox/GAPIT/Functions/GAPIT.library.R")
# source("/Users/Jiabo/Dropbox/GAPIT/Functions/gapit_functions.txt")

myGD <- read.table(file = "mdp_numeric.txt", head = T)
taxa <- myGD[, 1]
favorite <- c("33-16", "38-11", "B73", "B73HTRHM", "CM37", "CML333", "MO17", "YU796NS")
index <- taxa %in% favorite
snps <- myGD[, -1]

# K=GAPIT.kinship.loiselle(t(myGD[,-1]), method="additive", use="all")
K[index, index]

K1 <- GAPIT.kinship.VanRaden(snps)
K1[index, index]

K2 <- GAPIT.kinship.Zhang(snps)
K2[index, index]

heatmap.2(K1, cexRow = .2, cexCol = 0.2, col = rev(heat.colors(256)), scale = "none", symkey = FALSE, trace = "none")
quartz()
heatmap.2(K2, cexRow = .2, cexCol = 0.2, col = rev(heat.colors(256)), scale = "none", symkey = FALSE, trace = "none")

n <- nrow(myGD)
ind.a <- seq(1:(n * n))
i <- 1:n
j <- (i - 1) * n
ind.d <- i + j
par(mfrow = c(1, 3))
plot(K2[ind.a], K1[ind.a], main = "All elements", xlab = "Zhang", ylab = "VanRaden")
lines(K2[ind.d], K1[ind.d], main = "All elements", xlab = "Zhang", ylab = "VanRaden", col = "red", type = "p")
plot(K2[ind.d], K1[ind.d], main = "Diagonals", xlab = "Zhang", ylab = "VanRaden")
plot(K2[-ind.d], K1[-ind.d], main = "Off diag", xlab = "Zhang", ylab = "VanRaden")
