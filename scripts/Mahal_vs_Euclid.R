# (2.109) page 53 of linear models book Rencher Schaalje
# function to compute the inverse square root of a matrix
# https://scicomp.stackexchange.com/questions/10375/efficient-computation-of-the-matrix-square-root-inverse
# good response by Borchers
fnMatSqrtInverse = function(mA) {
  ei = eigen(mA)
  d = ei$values
  d = (d+abs(d))/2
  d2 = 1/sqrt(d)
  d2[d == 0] = 0
  return(ei$vectors %*% diag(d2) %*% t(ei$vectors))
}

# rotation parameters
theta <- pi/4
maj.min <- c(1,2)

## A matrix: rot * lambda * rot2
R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),2,2)
A <-  R %*% diag(maj.min) %*% t(R)

## square root of A
C <- eigen(A)$vectors
D <- diag(eigen(A)$values)
all.equal(A,C%*%D%*%t(C))
A_half <- C%*%sqrt(D)%*%t(C)

## A^-1/2 matrix: lambda * rot2
A_neghalf <- diag(maj.min) %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)),2,2)

all.equal(solve(fnMatSqrtInverse(A) %*% fnMatSqrtInverse(A)), A)

A_neghalf <- fnMatSqrtInverse(A)

n <- 100
orig_pts <- cbind(cos((pi/n)*(1:(2*n))),sin((pi/n)*(1:(2*n))))
aniso_pts <- matrix(NA, 2*n, 2)
for(i in 1:(2*n)) aniso_pts[i,] <-  t(A %*% orig_pts[i,])

plot(orig_pts, type="l", xlim = c(-3,3), ylim=c(-3,3), asp=1)
lines(aniso_pts, type="l", xlim = c(-3,3), ylim=c(-3,3), asp=1)

for (i in 1:(2*n)){points((A_neghalf %*% aniso_pts[i,])[1], (A_neghalf %*% aniso_pts[i,])[2], col = "blue")}
for (i in 1:(2*n)){points((A_neghalf %*% A %*% orig_pts[i,])[1], (A_neghalf %*% A %*% orig_pts[i,])[2], col = "purple")}
for (i in 1:(2*n)){points((A %*% orig_pts[i,])[1], (A %*% orig_pts[i,])[2], col = "red")}

# points to eval (compare distances)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    pts_e <- c(i,j)
    
    # calculate orig distance
    sqrt((orig_pts[pts_e[1],1]-orig_pts[pts_e[2],1])^2 + (orig_pts[pts_e[1],2]-orig_pts[pts_e[2],2])^2)
    diff_o <- orig_pts[pts_e[1],] - orig_pts[pts_e[2],]
    sqrt(t(diff_o)%*%diff_o)
    c(dist(orig_pts[pts_e,]))
    
    # sqrt((aniso_pts[1,1]-aniso_pts[2,1])^2 + (aniso_pts[1,2]-aniso_pts[2,2])^2)
    sqrt(mahalanobis(aniso_pts[pts_e[1],], center = aniso_pts[pts_e[2],], cov = A%*%A))
    diff_a <- aniso_pts[pts_e[1],]-aniso_pts[pts_e[2],]
    sqrt(t(diff_a)%*%solve(t(A)%*%A)%*%diff_a)
    
    if(!all.equal(c(dist(orig_pts[pts_e,])),sqrt(mahalanobis(aniso_pts[pts_e[1],], center = aniso_pts[pts_e[2],], cov = A%*%A)))){
      stop("this is trouble; mahal != euclid after adjustment")
    }
    
  }
}
