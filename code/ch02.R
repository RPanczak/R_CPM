library("arm")

N <- 1000
sigma <- 10
y <- rnorm(N, 0, sigma)
x1 <- sample(c(-0.5,0.5), N, replace=TRUE)
x2 <- sample(c(-0.5,0.5), N, replace=TRUE)
display(lm(y ~ x1))
display(lm(y ~ x1 + x2 + x1:x2))
# this was with y under the Null

# now with y under the alternative of separate x1 effects for x2 values
# specifically:
# overall effect of x1 = 2.8 * sigma; for x2==-.5: 2.1 * sigma; for x2==0.5: 3.5 * sigma
y[x1==.5 & x2== -.5] <- rnorm(length(y[x1==.5 & x2==-.5]), 2.1*sigma, sigma)
y[x1==.5 & x2== .5] <- rnorm(length(y[x1==.5 & x2== .5]), 3.5*sigma, sigma)
display(lm(y ~ x1)) # SE 0.72
display(lm(y ~ x1 + x2)) # SE 0.69
display(lm(y ~ x1 + x2 + x1:x2)) # SE interaction 1.31; 1.31/0.72 equals 1.82 rather than a factor 2
