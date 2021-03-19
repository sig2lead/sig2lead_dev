####################################
# User Input Variable (Gene KD)
####################################
library(httr)

define_knockdown <- function(gene)
{
  usr <- gene
  
  #usr <- "BCL2A1"
  if(input$ILINCS == "Legacy"){
  url1 <- paste("http://ilincs2018.ilincs.org/ilincs/api/SignatureMeta?filter=%7B%22where%22%3A%7B%22treatment%22%3A%22",usr,"%22%2C%20%22libraryid%22%3A%22LIB_6%22%7D%7D", sep="")
  #url1 <- paste("http://ilincs.org/api/SignatureMeta?filter=%7B%22where%22%3A%7B%22treatment%22%3A%22",usr,"%22%2C%20%22libraryid%22%3A%22LIB_6%22%7D%7D", sep="")
  
  raw.result <- GET(url = url1)
  req <- GET(url = url1)
  raw.result$status_code
  sigid <<- list("character")
  
  if(length(content(raw.result)) != 0){
    for (i in 1:length(content(raw.result)))
    {
      sigid[i] <- content(raw.result)[[i]]$signatureid
      
      
      
    }
  } else{
    print("A knockdown signature for youe gene target is not in LINCS.  Would you like to try a different gene target?")
    sigid <- NULL
  }
  
  sigid <<- sigid
  }
  else if(input$ILINCS == "Current"){
    #url1 <- paste("http://www.ilincs.org/api/SignatureMeta?filter=%7B%22where%22%3A%7B%22treatment%22%3A%22",usr,"%22%7D%7D", sep="")
    url1 <- paste("http://www.ilincs.org/api/SignatureMeta?filter=%7B%22where%22%3A%7B%22treatment%22%3A%22",usr,"%22%2C%20%22libraryid%22%3A%22LIB_6%22%7D%7D",sep="")
    req <- GET(url = url1)
    sigid <-fromJSON(httr::content(req,type="text", encoding="UTF-8"))$signatureid
    
    if (length(sigid) == 0){
      sigid <<- NULL
      print("A knockdown signature for youe gene target is not in LINCS.  Would you like to try a different gene target?")
      return()
    }
    else {
      sigid <<- sigid
    }
    
  }
  #return(sigid)
}
