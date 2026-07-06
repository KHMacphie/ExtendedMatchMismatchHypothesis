## Function to organise posteriors

#organise chains into one posterior for each parameter
stanpost <- function(model, parameters){
  
  df <- as.data.frame(extract(model, pars=parameters, permuted=FALSE))
  chns <- model@sim$chains
  npara <- (ncol(df))/chns
  
  # Organise + combine chains of posterior distributions
  remove.start <- function(x, n){  #function to remove first n number of characters in a string
    substr(x, nchar(x)-(nchar(x)-n-1), nchar(x))
  }
  
  fullpost <- data.frame(matrix(NA, nrow = chns*nrow(df), ncol = npara))
  
  for(i in 1:npara){
    fullpost[,i] <- stack(df[,(i*chns-(chns-1)):(i*chns)])[,1]
    colnames(fullpost)[i] <- remove.start(colnames(df)[(i*chns-(chns-1))],8)
  }
  
  return(fullpost)
}  


## Function to approximate probability scale results
diggle_approx <- function(mean,var){
  plogis(mean/sqrt(1+((16*sqrt(3))/(15*pi))^2*var))
}
