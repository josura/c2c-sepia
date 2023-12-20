library(MASS)

matriceTest <- t(matrix(c(0,0.2,0.4,0.5,0,0,0,0,
                        0.2,0,0.4,0.5,0,0,1,0,
                        0.3,0.2,0,0.5,0,0,0,0,
                        0.4,0.2,0.4,0,0,0,0,1,
                        0,1,0,0,0,0,0,0,
                        1,1,1,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0),nrow=8))

normalizedMatrix <- apply(t(matriceTest), 2, function(col) {
  if (sum(col) != 0) {
    col / sum(col)
  } else {
    rep(0, length(col))
  }
})


# Creating the identity matrix with the same dimensions as normalizedMatrix
identityMatrix <- diag(nrow(normalizedMatrix))

# Subtracting the transpose of the normalizedMatrix from the identity matrix
resultMatrix <- identityMatrix - normalizedMatrix

inverted <- ginv(resultMatrix)

input <- c(1,1,0.1,2,0,0,0,0)

perturbation <- inverted %*% input

perturvation.conserverd <- inverted %*% input - (input/2)

perturbation.iterated <- inverted %*% perturbation

for(i in range(0,20)){
  perturbation.iterated <- inverted %*% perturbation.iterated
}
