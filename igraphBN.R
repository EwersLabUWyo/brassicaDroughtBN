# igraphBN.R
# R version 3.3.1 (2016-06-21)
# June 15, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Plotting a BN using igraph. 

#-----------------------------------------------------------------------
library(bnlearn)
library(igraph)
#-----------------------------------------------------------------------

# Load network string. 
dag <- model2network("[SM][M10][M6|M10][gs|SM:M6][Photo|gs:M10]
                      [NSC|Photo][M7|Photo:M6][Starch|NSC][M9|M7:M10]
                     [Fv'/Fm'|Photo:M9][M2|M7:M9][M3|M6:M9]
                     [M1|Fv'/Fm':M2:M10][M5|Fv'/Fm':M3:M9][M4|M1:M9]
                     [M8|M1:M6:M10]")

# Create graph from network arcs. 
g <- graph_from_data_frame(dag$arcs)

# Identify nodes of graph. 
g[1]

# Set node colors.
col = c("lightblue3", "lightcoral", "lightblue3", rep("deepskyblue4", 3), 
        rep("lightblue3", 2), "deepskyblue4", rep("lightblue3", 3), 
        "deepskyblue4", rep("lightblue3", 3))

# Format layout for graph. 
l <- layout_in_circle(g)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

# Plot graph. 
plot(g, vertex.size = 30, 
     vertex.label.cex = 1.3, 
     vertex.color = col, 
     edge.arrow.size = .5,
     vertex.frame.color = "white", 
     vertex.label.color = "white", 
     layout = l * 1.18, rescale = F)




