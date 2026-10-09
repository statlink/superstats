createalfa <- function(alfax) {
  if ( nrow(alfax) == 1 ) { #if has one row just return the transposed
    Alfa <- t(alfax);
  } else {
    tmp <- createalfa( alfax[-1, , drop = FALSE] );
    Alfa <- NULL
    for ( i in 1:ncol(alfax) )  Alfa <- rbind(Alfa, cbind( rep.int(alfax[1, i], nrow(tmp) ), tmp) );
  }
  Alfa
}
