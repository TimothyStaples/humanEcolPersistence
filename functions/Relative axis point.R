relative.axis.point<-function(proportion, axis){
  ## AUTHOR: Timothy Staples
  ## DATE: 01/07/2016
  
  ## FUNCTION PURPOSE: Give proportion along axis in axis units. 
  ##                   To return meaningful values, this function should be run 
  ##                   after the main plot() command.
  
  ## Arguments: proportion: numeric, generally 0 to 1
  ##            axis: string, define axis to use. Must be "x" or "y"
  
  if(!axis %in% c("x","y")){stop('axis must equal "x" or "y"')}
  if(!is.numeric(proportion)){stop('proportion must be a number')}
  
  if(axis=="x"){
    return(par("usr")[1] + proportion*(par("usr")[2]-par("usr")[1]))
  }
  
  if(axis=="y"){
    return(par("usr")[3] + proportion*(par("usr")[4]-par("usr")[3]))
  }
}
