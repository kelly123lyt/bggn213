
x <- "4AKE"
x <- "1AKE"

pdbplot <- function(x){
  

  s <- read.pdb(x) 

  sc <- trim.pdb(s, chain="A", elety="CA")

  plotb3(sc$atom$b, sse=s1.chainA, typ="l", ylab="Bfactor")
  
  return(sc$atom$b)
}