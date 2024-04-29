library(ggnetwork)
library(GGally)
library(sna)
library(ggplot2)
library(patchwork)

#Used to build the tree branches for visualization
build_tree_branches = function(x, ...){
  
  splitvar = 0L
  p.value = 0L
  id = 0L
  
  ggparty::ggparty(x, terminal_space = 0) +
    ggparty::geom_edge() +
    ggparty::geom_edge_label() +
    ggplot2::theme(legend.position = "none") +
    ggparty::geom_node_label(line_list = list(
      ggplot2::aes(label = splitvar)),
      line_gpar = list(list(size = 10,
                            col = "black",
                            fontface = "bold",
                            alignment = "center")
      ),
      ids = "inner") +
    ggplot2::coord_cartesian(ylim = c(0.1, 1.1))
}

nets_maker <- function(L, tree, graph_method, nodelabel, stages, disc_cont, pal, J=3,
                       sim=100,thresh_edge = 0.85, thresh_dir = 0.5,
                       bl = NULL, alpha = 0.05){
  
  #Set x and y coordinates for nodes in graph so they all have the same structure
  full_mat <- matrix(rep(1,(length(nodelabel)*length(nodelabel))),byrow = TRUE, ncol = length(nodelabel))
  circle_x <- gplot.layout.circle(full_mat,NULL)[,1]
  circle_y <- gplot.layout.circle(full_mat,NULL)[,2]
  
  #Start by making tree and find data behind terminal nodes
  # run this for the condensended rankings using predicted nodes
  nodes_tree <- predict(tree, type = "node")
  
  node_id_tree <- sort(unique(nodes_tree))
  
  nets <- list()
  
  for (i in seq_along(node_id_tree)) {
    
    print(node_id_tree[i])
    
    #Get the data in a specific node
    Li <- L[rep(nodes_tree == node_id_tree[i], each = J),]
    
    if (disc_cont == "discrete") {
      Wi <- X_maker(Li)
      H <- "N"
    } else if (disc_cont == "cont") {
      Wi <- Y_maker(Li)
      H <- "J"
    }
    
    if (graph_method == "pcor" | graph_method == "kcor") {
      #Get GGM
      if (graph_method == "pcor") {
        if (disc_cont == "discrete") {
          mat <- pcor_mat_makerX(Li, alpha = alpha)
        } else if (disc_cont == "cont") {
          mat <- pcor_mat_makerY(Wi, alpha = alpha)
        }
        #Get a matrix based on kendall tau correlations
      } else if (graph_method == "kcor") {
        mat <- kcor_mat_maker(Li)
      }
      
      graph_mat <- mat
      graph_mat[!(graph_mat == 0)] <- 1
      net_i <- network::network(graph_mat, directed = FALSE)
      edge_weights <- rep(0,network.edgecount(net_i))
      edge_color <- rep("grey",network.edgecount(net_i))
      
      k <- 0
      
      for (m in 1:ncol(graph_mat)) {
        for (n in 1:nrow(graph_mat)) {
          if (n > m) {
            if (graph_mat[m,n] != 0) {
              k <- k + 1
              edge_weights[k] <- abs(mat[m,n]*2.5)
              if (mat[m,n] > 0) {
                edge_color[k] <- "blue"
              } else{
                edge_color[k] <- "red"
              }
              
            }
          }
        }
      }
      
      
      net_i %e% "weight" <- edge_weights
      
      #Turn the following two lines off if you don't want a circle
      net_i %v% "x" <- circle_x
      net_i %v% "y" <- circle_y
      
      net_i %v% "Stage" <- stages
      names(pal) <- unique(stages)
      
      nets[[i]] <- ggnet2(net_i, edge.size = "weight", label = nodelabel, label.size = 3,
                          node.color = "Stage", palette = pal, label.color = "white", 
                          edge.color = edge_color, mode = c("x","y"),
                          node.size = 6) +
        labs(title = paste(H, "=", nrow(Wi))) + 
        theme(plot.title = element_text(color="black", size=8, face="bold.italic",hjust = 0.5)) +
        theme(panel.border = element_rect(color = "grey50", fill = NA),
              aspect.ratio = 1)
      
    } else if (graph_method == "BN") {
      
      adj_mat <- matrix(0, nrow = length(nodelabel), ncol = length(nodelabel))
      rownames(adj_mat) = colnames(Wi)[-ncol(Wi)]
      colnames(adj_mat) = colnames(Wi)[-ncol(Wi)]
      graph_mat <- adj_mat
      u <- 0
      
      while (u < sim) {
        
        rd_farm_idx <- sample(seq(1,length(Li$id),J),size = (length(Li$id)/J),replace = TRUE)
        add_idx <- rep(c(0,1,2), times = length(rd_farm_idx))
        rd_farm_idx <- rep(rd_farm_idx,each = J)
        
        bootstrap_L_idx <- rd_farm_idx + add_idx
        
        L_boot <- Li[bootstrap_L_idx,]
        
        if (disc_cont == "discrete") {
          W_boot <- X_maker(L_boot)
        } else if (disc_cont == "cont") {
          W_boot <- NULL
          try(W_boot <- Y_maker(L_boot))
          if (is.null(W_boot)) {
            next            
          }
        }
        
        BN_boot <- NULL
        
        try(BN_boot <- BN_maker(W_boot, bl = bl),TRUE)
        
        if (is.null(BN_boot)) {
          next
        } else{
          u <- u + 1
        }
        
        arcs_boot <- BN_boot$arcs
        n_arcs_boot <- length(arcs_boot)/2
        
        for (j in 1:n_arcs_boot) {
          adj_mat[arcs_boot[j],arcs_boot[j+n_arcs_boot]] <- adj_mat[arcs_boot[j],arcs_boot[j+n_arcs_boot]] + 1
        }
      }
      
      for (l in 1:ncol(adj_mat)) {
        for (j in 1:nrow(adj_mat)) {
          if (((adj_mat[l,j] + adj_mat[j,l])/sim > thresh_edge) & (adj_mat[l,j]/sim > thresh_dir)) {
            graph_mat[l,j] <- 1
          }
        }
      }
      
      print(graph_mat)
      
      net_i <- network::network(graph_mat, directed = TRUE)
      edge_color <- rep("grey",network.edgecount(net_i))
      
      net_i %v% "Stage" <- stages
      
      #Turn the following two lines off if you don't want a circle
      net_i %v% "x" <- circle_x
      net_i %v% "y" <- circle_y
      
      names(pal) <- unique(stages)
      
      if (sum(graph_mat) == 0) {
        nets[[i]] <- ggnet2(net_i, label = nodelabel, label.size = 3,
                            label.color = "white", mode = c("x","y"),
                            node.color = "Stage", palette = pal,
                            node.size = 6) +
          labs(title = paste(H,"=", nrow(Wi))) + 
          theme(plot.title = element_text(color="black", size=8, face="bold.italic",hjust = 0.5)) +
          theme(panel.border = element_rect(color = "grey50", fill = NA),
                aspect.ratio = 1)
      } else{
        nets[[i]] <- ggnet2(net_i, label = nodelabel, label.size = 3,
                            label.color = "white", mode = c("x","y"),
                            node.color = "Stage", palette = pal, edge.color = edge_color,
                            node.size = 6, arrow.size = 6, arrow.gap = 0.04) +
          labs(title = paste(H,"=", nrow(Wi))) + 
          theme(plot.title = element_text(color="black", size=8, face="bold.italic",hjust = 0.5)) +
          theme(panel.border = element_rect(color = "grey50", fill = NA),
                aspect.ratio = 1)
        
      }
      
    }
  }
  
  return(nets)
}

