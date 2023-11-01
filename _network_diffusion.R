# INITIALIZATION ----------------------------------------------------------
# This is a nice color palette.
cols <- pals::tableau20(20)
# Take a look at it.
scales::show_col(cols)

# We also want to have a certain directory structure to store output and saved
# data objects. For me this is 1. data - 1.1 output; 1.2 RDS Show warnings is 
# disabled so that we are not bothered if the directories already exist.
dir.create("data", showWarnings = FALSE)
dir.create("data/output", showWarnings = FALSE)
dir.create("data/RDS", showWarnings = FALSE)
dir.create("data/RAW", showWarnings = FALSE)

# Create a variable to our working directory for produced outputs
WORKDIR <- getwd()
FILESOURCE <- file.path(WORKDIR, "data")

# LOAD LIBRARIES ----------------------------------------------------------

require(diffuStats)
require(STRINGdb)
require(igraph)

library(GSEABase)
library(msigdb)
library(Category)
library(annotate)
library(org.Hs.eg.db)


# FUNCTIONS ---------------------------------------------------------------
## This section provides you with some helper-functions.

# 'prepareSTRINGdb' makes DB and PPI-graph available to R in a single command
# The function will create 'data/STRINGdb' in the current directory, download
# the database and make both DB and PPI-graph available in 'stringDB' object
prepareSTRINGdb <- function() {
  require(STRINGdb)
  # create a directory to store data
  dir.create("data", showWarnings = FALSE)
  dir.create("data/STRINGdb", showWarnings = FALSE)
  # download database
  string_db <- STRINGdb$new( version="12", species=9606,score_threshold=600, 
                             input_directory=paste(getwd(), "/data/STRINGdb/",sep=""))
  # extract graph
  string.g <- string_db$get_graph()
  node_ensembl <- igraph::V(string.g)$name
  # get protein information
  string_proteins <- string_db$get_proteins()
  rownames(string_proteins) <- string_proteins$protein_external_id	# adjust df rownames
  node_symbol <- string_proteins[node_ensembl, 2]		# map to gene name
  igraph::V(string.g)$node.name <- node_symbol				# replace string_ids with gene-names
  # keep the ensembl ID
  igraph::V(string.g)$node.ensembl <- node_ensembl
  igraph::V(string.g)$name <- node_symbol
  
  sDB.list <- list()
  sDB.list[["string.DB"]] <- string_db
  sDB.list[["string.g"]] <- string.g
  
  assign("stringDB",sDB.list, .GlobalEnv)
}

## return subgraph with specified threshold for degree
#'@param graph, the input-graph to work on
#'@param degree, threshold for degree of connectivity, i.e., max. number of edges
#'@return subgraph, constructed ipgraph::induced.subgraph
getSubDegree <- function( graph, degree ) {
  deg <- sort(igraph::degree(graph), decreasing=T)
  idx <- which(deg < degree)
  idx2 <- which(igraph::V(graph)$name %in% names(idx))
  gi <- igraph::induced.subgraph(graph, idx2)
  
  return(gi)
}


## Function is deprecated and was replaced by 'prepareSTRINGdb'
## extract graph from STRINGdb and transform for further use
#'@param STRINGdb, stringdb object to extract and construct graph from
#'@return string.g, igraph stringDB object
prepareGraph <- function( STRINGdb ) {
  string.g <- STRINGdb$get_graph()
  node_ensembl <- igraph::V(string.g)$name
  # get protein information
  string_proteins <- STRINGdb$get_proteins()
  rownames(string_proteins) <- string_proteins$protein_external_id	# adjust df rownames
  node_symbol <- string_proteins[node_ensembl, 2]		# map to gene name
  igraph::V(string.g)$name <- node_symbol				# replace string_ids with gene-names
  # remove expendable information?!
  string.g <- igraph::remove.edge.attribute(string.g, "neighborhood")
  string.g <- igraph::remove.edge.attribute(string.g, "neighborhood_transferred")
  string.g <- igraph::remove.edge.attribute(string.g, "fusion")
  string.g <- igraph::remove.edge.attribute(string.g, "cooccurence")
  string.g <- igraph::remove.edge.attribute(string.g, "homology")
  string.g <- igraph::remove.edge.attribute(string.g, "coexpression")
  string.g <- igraph::remove.edge.attribute(string.g, "coexpression_transferred")
  string.g <- igraph::remove.edge.attribute(string.g, "experiments")
  string.g <- igraph::remove.edge.attribute(string.g, "experiments_transferred")
  string.g <- igraph::remove.edge.attribute(string.g, "database")
  string.g <- igraph::remove.edge.attribute(string.g, "database_transferred")
  string.g <- igraph::remove.edge.attribute(string.g, "textmining")
  string.g <- igraph::remove.edge.attribute(string.g, "textmining_transferred")
  
  return(string.g)
}


## 'diffuse' performs NWD on given network using specified score-type
#'@param network, reduced STRINGdb network
#'@param K_rl, pre-computed diffusion-kernel
#'@param input, data.frame containing at least two cols 'node_id' and 'score'
#'@param method, character string denoting method of diffusion
#'@return NULL
diffuseWrapper <- function(network, K_rl, input, method) {
  
  x <- input
  # remove duplicates (should be non left actually)
  x <- x[ !duplicated(x$node_id), ]
  # now take remaining Genes and create a binary input_vector for the diffusion
  tmp_input_vec <- rep(0,length(igraph::V(network)$name))
  names(tmp_input_vec) <- igraph::V(network)$name
  idx <- intersect( names(tmp_input_vec), as.character(x$node_id))
  # in case the remaining genes in x cannot be mapped to the network - remove 
  # unmapped genes from x
  x <- x[which(x$node_id %in% idx),]
  tmp_input_vec[match(x$node_id, names(tmp_input_vec))] <- x$node_score
  # diffusion
  if(method == "ber_p") {
    df_z <- diffuStats::diffuse_grid(K = K_rl, scores=tmp_input_vec,
                                   grid_param = expand.grid(method = eval(method)), n.perm = 10)
  } else {
    df_z <- diffuStats::diffuse_grid(K = K_rl, scores=tmp_input_vec,
                                     grid_param = expand.grid(method = eval(method)), n.perm = 1)
    
  }
  result <- list()
  result[["dgrid"]] <- df_z
  result[["x"]] <- x
  
  return(result)
}

# plot all the graphs
#'@param input, output of 'diffuseWrapper' containing input-GOI (x) and diffusion results (df_z)
#'@param network, background-network the kernel is based on
#'@param outdir, path to plot to
#'@param score, String denoting the initial diffusion value
#'@param cutoff, threshold for the the diffusion scores to include in graph (quantile [0,1])
plotDiffNets <- function(input, network, outdir, score, cutoff=0.995) {

  df_z <- input$dgrid
  x <- input$x
  # Create a subnetwork out of the diffused scores
  # The cutoff below is chosen depending on the data, please adjust accordingly
  tmpCutoff <- quantile(abs(df_z$node_score), cutoff)
  idx2 <- which(abs(df_z$node_score) > tmpCutoff)
  diffused_net <- igraph::induced.subgraph(graph=network, igraph::V(network)$name[idx2])
  diffused_net <- igraph::set_vertex_attr(diffused_net, "weight", index = igraph::V(diffused_net), df_z$node_score[idx2])
  
  # use edge_betweenness to cluster in communities and assign different colors
  com <- igraph::cluster_edge_betweenness(diffused_net, modularity=TRUE)
  com$csize <- sapply( 1:length(com), function(x) sum(com$membership==x) )
  vgroups <- com$membership
  colormap <- "blue-darkorange"
  palette.name <- supraHex::visColormap( colormap=colormap )
  mcolors <- palette.name(length( com ))
  vcolors <- mcolors[vgroups]
  com$significance <- dnet::dCommSignif( diffused_net, com )
  # colouring for adults
  mark.groups <- communities(com)
  mark.col <- supraHex::visColoralpha(mcolors, alpha=0.2)
  mark.border <- supraHex::visColoralpha(mcolors, alpha=0.2)
  edge.color <- c("#C0C0C0","#000000")[igraph::crossing(com, diffused_net)+1]
  edge.color <- supraHex::visColoralpha(edge.color, alpha=0.5)
  
  mysizes <- rep(2, length(V(diffused_net)$name))
  names(mysizes) <- V(diffused_net)$name
  
  tmp_shape_vec <- rep("sphere",length(V(diffused_net)$name))
  names(tmp_shape_vec) <- V(diffused_net)$name
  idx <- intersect( names(tmp_shape_vec), x$node_id)
  tmp_shape_vec[idx] <- "csquare"		# scores for GOI set to 1
  
  
  png(file= paste0(outdir, "/net_diff_", score, ".png"), height = 12, width = 12, res=300, units = "in")
  dnet::visNet( diffused_net,
                glayout=layout.fruchterman.reingold,
                newpage=F, vertex.label.cex=0.75,
                vertex.color = diffuStats::scores2colours(df_z$node_score[idx2]),
                vertex.shape = tmp_shape_vec,
                mark.groups=mark.groups, mark.col=mark.col, edge.color=edge.color,
                mark.shape=1, mark.expand=20, vertex.size=2.5,
                main = paste0("Diffusion network - subject: ",score)
  )
  dev.off()

}

# Perform pathway analysis and print result to excel file 
#'@param input, list of subjects contain (thresholded) GOI and list with diffusion scores
#'@param network, background-network the kernel is based on
#'@param score, String denoting the initial diffusion value
#'@param qcutoff, threshold for quantile(node_score, qcutoff)
#'@param pcutoff, threshold for selection of pvalue
#'@param useDB, character vector containing denominator for DBs to use, i.e., HALLMARK,
# REACTOME, MOTIF, GO, IMMU
#'@param path, save directory
printEnrichments <- function(network, input, score, qcutoff, pcutoff, useDB, path) {
  ### DEBUG ###
  # input <- input_list
  # outdir <- FILESOURCE
  # score <- 1
  # qcutoff <- 0.75
  # useDB <- c("HALLMARK", "REACTOME")
  # path <- FILESOURCE
  ### DEBUG ###
  
  print("init DBs")
  # check for default setup
  if( missing(pcutoff) ) {
    pcutoff <- 0.05
  }
  
  fields <- names(input)
  result <- list()
  if( "HALLMARK" %in% useDB ) {
    hm <- MSigDB[["HALLMARK"]]
    idx <- grep("HALLMARK",names(hm))
    hm <- hm[idx]
    gsc_hall <- GeneSetCollection(mapply(function(geneIds,hallname) {
      GeneSet(geneIds,geneIdType=EntrezIdentifier(),
              setName=hallname)
    },hm,names(hm)))
  }
  
  if( "REACTOME" %in% useDB ) {
    c2 <- MSigDB[["C2_CURATED"]] #Reactome
    idx <- grep("REACTOME",names(c2))
    c2 <- c2[idx]
    gsc_reac <- GeneSetCollection(mapply(function(geneIds,reactomename) {
      GeneSet(geneIds,geneIdType=EntrezIdentifier(),
              setName=reactomename)
    },c2,names(c2)))
  }
  
  if( "MOTIF" %in% useDB ) {
    c3 <- MSigDB[["C3_MOTIF"]]
    idx <- grep("",names(c3))
    c3 <- c3[idx]
    gsc_motif <- GeneSetCollection(mapply(function(geneIds,motifname) {
      GeneSet(geneIds,geneIdType=EntrezIdentifier(),
              setName=motifname)
    },c3,names(c3)))
  }
  
  if( "GO" %in% useDB ) {
    c5 <- MSigDB[["C5_GENE_ONTOLOGY"]]
    idx <- grep("GO",names(c5))
    c5 <- c5[idx]
    gsc_go <- GeneSetCollection(mapply(function(geneIds,goname) {
      GeneSet(geneIds,geneIdType=EntrezIdentifier(),
              setName=goname)
    },c5,names(c5)))
  }
  
  if( "IMMU" %in% useDB ) {
    c7 <- MSigDB[["C7_IMMUNOLOGIC_SIGNATURES"]]
    idx <- grep("",names(c7))
    c7 <- c7[idx]
    gsc_immu <- GeneSetCollection(mapply(function(geneIds,immuname) {
      GeneSet(geneIds,geneIdType=EntrezIdentifier(),
              setName=immuname)
    },c7,names(c7)))
  }
  
  print("finished init")
  
  # @param: direction, "over" or "under"
  # @param: tname, name this test
  diffNet <- function(diffused_net, gSC, universe, pcutoff, direction, mincount, tname) {
    # Genes in diffusion subnetwork
    my_symbols <- V(diffused_net)$name
    # prepare object for hyperGTest
    params <- GSEAKEGGHyperGParams(name=tname, geneSetCollection=gSC,
                                   geneIds=my_symbols, universeGeneIds=universe,
                                   pvalueCutoff = pcutoff, testDirection = direction)
    go_over <- hyperGTest(params)
    HyperEnrichment <- subset ( summary(go_over), Count >= mincount)
    
    return(HyperEnrichment)
  }
  
  # @param: direction, "over" or "under"
  # @param: tname, name this test
  diffNetGOI <- function(diffused_net, x, gSC, universe, pcutoff, direction, mincount, tname) {
    my_symbols <- as.character(x$GeneName[ x$GeneName %in% universe ])
    params <- GSEAKEGGHyperGParams(name=tname, geneSetCollection=gSC,
                                   geneIds=my_symbols, universeGeneIds=universe,
                                   pvalueCutoff = pcutoff,testDirection = direction)
    go_over <- hyperGTest(params)
    HyperEnrichment <- subset (summary (go_over), Count >= mincount)
    
    return(HyperEnrichment)
  }
  
  # minimal number of genes
  mincount <- 2
  # All genes
  MyGeneUniverse <- V(network)$name
  for( i in fields ) {
    print(i)
    my_Pathways <- list() # empty list
    df_z <- input[[i]]$dgrid
    x <- input[[i]]$x
    # Create a subnetwork out of the diffused scores
    # The cutoff below is chosed depending on the data, please adjust accordingly
    tmpCutoff <- quantile(abs(df_z$node_score), qcutoff)
    idx2 <- which(df_z$node_score > tmpCutoff)
    
    diffused_net <- induced.subgraph(graph=network, V(network)$name[idx2])
    diffused_net <- igraph::set_vertex_attr(diffused_net, "weight", index = V(diffused_net), df_z$node_score[idx2])
    
    if( "HALLMARK" %in% useDB ) {
      ## 1. HALLMARK
      # Over - Genes in diffusion subnetwork
      my_Pathways[["HALLMARK_Net_over"]] <- 
        diffNet(diffused_net, gsc_hall, MyGeneUniverse, pcutoff, "over", mincount, "My_HALL_1")
      # Over - Genes of Interest in diffusion subnetwork
      my_Pathways[["HALLMARK_GOI_over"]] <- 
        diffNetGOI(diffused_net, x, gsc_hall, MyGeneUniverse, pcutoff, "over", mincount, "My_HALL_2")
      # Under - Genes in diffusion subnetwork
      my_Pathways[["HALLMARK_Net_under"]] <- 
        diffNet(diffused_net, gsc_hall, MyGeneUniverse, pcutoff, "under", mincount, "My_HALL_3")
      # Under - Genes of Interest in diffusion subnetwork
      my_Pathways[["HALLMARK_GOI_under"]] <- 
        diffNetGOI(diffused_net, x, gsc_hall, MyGeneUniverse, pcutoff, "under", mincount, "My_HALL_4")
      print(paste("HALLMARK done for:", i, sep=" "))
    }
    
    if( "REACTOME" %in% useDB ) {
      ## 2. REACTOME
      # Over - Genes in diffusion subnetwork
      my_Pathways[["REACTOME_Net_over"]] <- 
        diffNet(diffused_net, gsc_reac, MyGeneUniverse, pcutoff, "over", mincount, "My_REACTOME_1")
      # Over - Genes of Interest in diffusion subnetwork
      my_Pathways[["REACTOME_GOI_over"]] <- 
        diffNetGOI(diffused_net, x, gsc_reac, MyGeneUniverse, pcutoff, "over", mincount, "My_REACTOME_2")
      # Under - Genes in diffusion subnetwork
      my_Pathways[["REACTOME_Net_under"]] <- 
        diffNet(diffused_net, gsc_reac, MyGeneUniverse, pcutoff, "under", mincount, "My_REACTOME_3")
      # Under - Genes of Interest in diffusion subnetwork
      my_Pathways[["REACTOME_GOI_under"]] <- 
        diffNetGOI(diffused_net, x, gsc_reac, MyGeneUniverse, pcutoff, "under", mincount, "My_REACTOME_3")
      print(paste("REACTOME done for:", i, sep=" "))
    }
    
    if( "MOTIF" %in% useDB ) {
      ## 3. MOTIF
      # Over - Genes in diffusion subnetwork
      my_Pathways[["MOTIF_Net_over"]] <- 
        diffNet(diffused_net, gsc_motif, MyGeneUniverse, pcutoff, "over", mincount, "My_MOTIF_1")
      # Over - Genes of Interest in diffusion subnetwork
      my_Pathways[["MOTIF_GOI_over"]] <- 
        diffNetGOI(diffused_net, x, gsc_motif, MyGeneUniverse, pcutoff, "over", mincount, "My_MOTIF_2")
      # Under - Genes in diffusion subnetwork
      my_Pathways[["MOTIF_Net_under"]] <- 
        diffNet(diffused_net, gsc_motif, MyGeneUniverse, pcutoff, "under", mincount, "My_MOTIF_3")
      # Under - Genes of Interest in diffusion subnetwork
      my_Pathways[["MOTIF_GOI_under"]] <- 
        diffNetGOI(diffused_net, x, gsc_motif, MyGeneUniverse, pcutoff, "under", mincount, "My_MOTIF_4")
      print(paste("MOTIF done for:", i, sep=" "))
    }
    
    if( "GO" %in% useDB ) {
      ## 4. GO
      # Over - Genes in diffusion subnetwork
      my_Pathways[["GO_Net_over"]] <- 
        diffNet(diffused_net, gsc_go, MyGeneUniverse, pcutoff, "over", mincount, "My_GO_1")
      # Over - Genes of Interest in diffusion subnetwork
      my_Pathways[["GO_GOI_over"]] <- 
        diffNetGOI(diffused_net, x, gsc_go, MyGeneUniverse, pcutoff, "over", mincount, "My_GO_2")
      # Under - Genes in diffusion subnetwork
      my_Pathways[["GO_Net_under"]] <- 
        diffNet(diffused_net, gsc_go, MyGeneUniverse, pcutoff, "under", mincount, "My_GO_3")
      # Under - Genes of Interest in diffusion subnetwork
      my_Pathways[["GO_GOI_under"]] <- 
        diffNetGOI(diffused_net, x, gsc_go, MyGeneUniverse, pcutoff, "under", mincount, "My_GO_4")
      print(paste("GO done for:", i, sep=" "))
    }
    
    if( "IMMU" %in% useDB ) {
      ## 5. IMMU
      # Over - Genes in diffusion subnetwork
      my_Pathways[["IMMU_Net_over"]] <- 
        diffNet(diffused_net, gsc_immu, MyGeneUniverse, pcutoff, "over", mincount, "My_IMMU_1")
      # Over - Genes of Interest in diffusion subnetwork
      my_Pathways[["IMMU_GOI_over"]] <- 
        diffNetGOI(diffused_net, x, gsc_immu, MyGeneUniverse, pcutoff, "over", mincount, "My_IMMU_2")
      # Under - Genes in diffusion subnetwork
      my_Pathways[["IMMU_Net_under"]] <- 
        diffNet(diffused_net, gsc_immu, MyGeneUniverse, pcutoff, "under", mincount, "My_IMMU_3")
      # Under - Genes of Interest in diffusion subnetwork
      my_Pathways[["IMMU_GOI_under"]] <- 
        diffNetGOI(diffused_net, x, gsc_immu, MyGeneUniverse, pcutoff, "under", mincount, "My_IMMU_4")
      print(paste("IMMU done for:", i, sep=" "))
    }
    
    openxlsx::write.xlsx(my_Pathways, 
                         file = paste0(path, "/" ,score," ", i,".xlsx"),#rowNames=TRUE,
                         colNames=TRUE, colWidths = c(NA, "auto", "auto"))
    # add it to the output matrix
    result[[i]] <- my_Pathways
    print("saving to file complete")
  }
  return(result)
}




