#####################################
# KDnet v1.0                        #
# Last update: Mar 21, 2019         #
# Author: Somchai Chutipongtanate   #
#####################################

# Knock down genes in LINCS; last update Mar 8, 2019 -----------

# KDgenes <- read.csv("~/Desktop/KDgenes.csv", stringsAsFactors = FALSE)
 #KDgenes <- KDgenes$GeneName
 #save(KDgenes, file = 'KDgenes')

#----------------------------------------------------------------


KDnet <- function(x){
  require(visNetwork)
  require(dplyr)
  require(readr) 
  load("./lib/KDgenes")

  #gene_knockdown <- "AAAS" # Aladin; Q9NRG9 (AAAS_HUMAN); Plays a role in the normal development of the peripheral and central nervous system
  #if(!(gene_knockdown %in% KDgenes) ){
  
  print("Not found your knockdowns")
  print("Start protein network analysis")
  api_KDnet <-paste0("http://stitch.embl.de/api/psi-mi-tab/interactions?identifier=", 
                     gene_knockdown, "&required_score=400&species=9606&limit=10")
  print(api_KDnet)
  KDnet <- suppressMessages(read_tsv(api_KDnet, col_names = FALSE)) %>% 
                                filter(X3 == gene_knockdown | X4 == gene_knockdown)
  KDnet <<- KDnet
  
  tmp <- strsplit(KDnet$X15, "|", fixed = TRUE) %>% sapply("[", 1) %>% 
    strsplit(":", fixed = TRUE) %>% sapply("[", 2) %>% 
    as.numeric()
  edge_KDnet <- data.frame("from" = KDnet$X3, "to" = KDnet$X4, "score" = tmp, "width" = (tmp*2)^2, 
                           "title" = paste0("<p>", "Score:", "<br>", tmp, "</p>"), stringsAsFactors = FALSE)
  edge_KDnet <<- edge_KDnet
  
  mol_KDnet <- unique(c(edge_KDnet$from, edge_KDnet$to)) 
  
  group_KDnet <- toupper(mol_KDnet) %>% replace(. %in% toupper(gene_knockdown), "Queried gene") %>%
                  replace(. %in% toupper(KDgenes), "Gene KD present in LINCS") %>%
                  replace(!. %in% c("Queried gene", "Gene KD present in LINCS"), "Not found in LINCS")
  
  link_node_KDnet <- character()
  for (i in seq_along(mol_KDnet)) {
      link_node_KDnet[i] <- paste0('<a target= "_blank"', '', 'href=', '', 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=', 
                                   mol_KDnet[i], '>', mol_KDnet[i], '</a>')}
      
  node_KDnet <- data.frame(label = mol_KDnet, id = mol_KDnet, group = group_KDnet, title = link_node_KDnet)
  node_KDnet <<- node_KDnet
#}
  print("Finish")
  visNetwork(node_KDnet, edge_KDnet) %>%
    visNodes(physics = FALSE) %>%
    visEdges(physics= FALSE, smooth = FALSE, color = list(color = "#8080ff" , highlight = "blue")) %>% #008080
    visGroups(groupname = "Queried gene", shape = "circle",
              font =  list(color = "blue"),
              color = list(background = "fdff00", border = "black",
                           highlight = list(background = "feff8f", 
                                            border = "red")) ) %>%
    visGroups(groupname = "Gene KD present in LINCS", shape = "circle",
              font = list(color = "red", size = 10),
              color = list(background = "caff44" , border = "black", 
                           highlight = list(background = "caff44", 
                                            border = "red"))) %>%
    visGroups(groupname = "Not found in LINCS", shape = "dot", size = 10,
              font = list(color = "grey50", size = 10),
              color = list(background = "ddd" , border = "black", 
                           highlight = list(background = "grey50", 
                                            border = "red")))  %>%
    visIgraphLayout(randomSeed = 9, layout = "layout_nicely", physics = FALSE, smooth = FALSE) 
           
 
}
