#eigencalc.R

#Makes a 2x2 matrix with two positive real eigenvalues
eig.make2x2 <- function(){
  repeat{
    c1 <- sample(c(-4:-1,1:4),2,replace = TRUE)
    c2 <- sample(c(-4:-1,1:4),2,replace = TRUE)
    A <- cbind(c1,c2)
    if (det(A)== 0)
      next
    if (det(A) <0)
      A <- cbind(c2,c1)
    trc <- A[1,1]+A[2,2]
    if (trc <0)
      A <- -A
    dt <- det(A)
    if (trc^2>4*dt)
      return(A)
  }
}

#this function will make the inital matrix
eig.make3x3 <- function(){
  repeat{
    c1 <- sample(c(1:4),3,replace = TRUE)
    c2 <- sample(c(1:4),3,replace = TRUE)
    c3 <- sample(c(1:4),3,replace = TRUE)
    A <- cbind(c1,c2,c3)
    
    non_invertible = tryCatch({
      solve(A)
    }, error = function(cond){
      return(T)
    })
    
    eig_A = eigen(A)
    is_imaginary = eig_A[[1]][1] %>% is.complex()
    
    if(is_imaginary | non_invertible){
      next
    } else{
      return(A)
    }
  }
}
#A <- eig.make2x2();A

#Diagonalizes a matrix with real eigenvalues
#The eigenvectors are unit vectors whose components sum to > 0
eig.diagonalize <- function(A) {
  trc <- A[1,1]+A[2,2]
  lam1 <- 0.5*(trc+sqrt(trc^2-4*det(A)))
  lam2 <- 0.5*(trc-sqrt(trc^2-4*det(A)))
  v1 <- c(A[1,1]-lam2,A[2,1])
  v2 <- c(A[1,1]-lam1,A[2,1])
  D <- diag(c(lam1,lam2))
  v1 <- v1/sqrt(sum(v1^2))
  if (sum(v1)<0)
    v1 <- -v1
  v2 <- v2/sqrt(sum(v2^2))
  if (sum(v2)<0)
    v2 <- -v2
  P <- cbind(v1,v2)
  return (list(diag = D, cofb = P))
}
#lst <- eig.diagonalize(A)
#D <- lst$diag;D
#P <- lst$cofb;P
#P%*%D%*%solve(P);A

#Raises a matrix with positive eigenvalues to a power
eig.nthpower <- function(A, power) {
  x <- eig.diagonalize(A)
  D <- x$diag
  P <- x$cofb
  Pinv <- solve(P)
  return(P%*%D^power%*%Pinv)
}
#B <- eig.nthpower(A,0.5);B%*%B;A
#q <- 0.5


#Plots unit square, parallelogram for A, and intermediate power
eig.plot <- function(A,q){
  x <- eig.diagonalize(A)
  P <- x$cofb
  D <- x$diag
  top <- max(A[2,1],A[2,2],A[2,1]+A[2,2],1,D[1,1]*P[2,1])
  bottom <- min(A[2,1],A[2,2],A[2,1]+A[2,2],-1,D[1,1]*P[2,1])
  right <- max(A[1,1],A[1,2],A[1,1]+A[1,2],1,D[1,1]*P[1,1])
  left <- min(A[1,1],A[1,2],A[1,1]+A[1,2],-1,D[1,1]*P[1,1])
  par(mar = c(2,2,1,1))
  plot(NULL, NULL, xlim = c(left, right),  ylim = c(bottom, top), asp = 1)
  #Draw standard basis vectors
  arrows(c(0,0),c(0,0),c(0,1),c(1,0),length = 0.1)
  segments(c(1,1),c(1,1),c(0,1),c(1,0), lty = 3)
  #Draw columns of A
  arrows(c(0,0),c(0,0),A[1,],A[2,],length = 0.1)
  segments(A[1,1]+A[1,2],A[2,1]+A[2,2],A[1,],A[2,], lty = 3)
  #Draw intermediate parallelogram
  G <- eig.nthpower(A,q)
  #Draw columns of G
  arrows(c(0,0),c(0,0),G[1,],G[2,],length = 0.1)
  segments(G[1,1]+G[1,2],G[2,1]+G[2,2],G[1,],G[2,], lty = 3)
  polygon(c(0,G[1,1],G[1,1]+G[1,2],G[1,2]),
          c(0,G[2,1],G[2,1]+G[2,2],G[2,2]),density = NULL,col = "beige")
  #Draw eigenvectors
  
  arrows(0,0,D[1,1]^q*P[1,1],D[1,1]^q*P[2,1],length = 0.1, lwd = 2, col = "blue")
  arrows(0,0,D[2,2]^q*P[1,2],D[2,2]^q*P[2,2],length = 0.1, lwd = 2, col = "blue")
}

#eig.plot(A,0.7)

