# saved as bloomfilter.R
# Uso de Librerias
library(digest)
library(magristtr)
library(gmp)
library(Rmpfr)

# funcion de mapeo, usando el algoritmo 1 
vbits <- function(x,m=1000,k=7){
  vec <- rep(FALSE,m)
  for (i in 1:k) {
    for (j in 1:length(x)) {
      hash <- digest(x[j], algo = "murmur32", serialize = FALSE, seed = i) %>%
        Rmpfr::mpfr(base = 16) %% 
        m %>%
        as.integer()
      
      vec[hash + 1] <- TRUE
    }
  }
  attributes(vec) <- list(m = m, k= k)
  
  return(vec)
}

#Generamos el vector de 48 bits con los nombres
#Y con 3 funciones hash murmur 32
vect <- vbits(c("Jesus", "Jhonatan", "Reymundo", "Cristian"), m = 48, k = 3)


#Funcion para testear o buscar un elemento en el vector de bits
# usando el algoritmo 2
is_vect <- function(x, vbits){
  k <- attr(vbits,"k")
  m <- attr(vbits,"m")
  
  for (i in 1:k) {
    hash <- digest(x,algo="murmur32",serialize = FALSE, seed=i) %>%
      Rmpfr::mpfr(base = 16)%%
      m %>%
      as.integer()
    if(!vbits[hash + 1]){
      return(FALSE)
    }
  }
  return(TRUE)
}

# Testeamos una cierta cadena de caracteres 
is_vect("j",vect)
is_vect("jes",vect)
is_vect("Jesus",vect)
