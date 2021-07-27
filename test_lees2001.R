library(RColorBrewer) #to use brewer.pal
library(fields) #to use designer.colors
source("~/BSC/test_morani.R") ##@## get Hexagon, for debugging

# Function that input test data format and write Visium-style output files
write_visium_data = function(data, visium=visium) {
  for (i in 0:(nrow(data)-1)) { # i:y
    if (i %% 2 == 0) col_ixs = seq(0, ncol(data)-1, 2)
    else col_ixs = seq(1, ncol(data)-1, 2)
    
    for (j in col_ixs) { # j:x
      value = data[i+1, j+1]
      if (is.na(value)) next
      if (value == 0) in_tissue = 0
      else { in_tissue = 1 }
      
      barcode = sprintf("Ai%d_j%d", i, j) # i:y, j:x
      xcoord = 10*(i+1) # i:y
      ycoord = 10*(j+1) # j:x
      
      print_str = sprintf("%s,%d,%d,%d,%d,%d", barcode, in_tissue, 
                          i, j, 10*i, 10*j)
      write(print_str, visium$positions_file, append = T)
      
      # if (value == 0) color = "yellow"
      # else color = colors[value]
      # Hexagon(j+1, nrow(data)-(i+1), col=color) # x, y - draw hexagon
    }
  }
}

# Visium-style file settings
create_visium_list = function(gene) {
  visium = list() # data holder
  visium$positions_path = paste("visium/", gene, ".graph_positions_list.csv", sep='') ##@## replace w/ arg
  visium$barcodes_path = paste("visium/", gene, ".graph_barcodes.tsv", sep='') # AAACAAGTATCTCCCA-1, ...
  visium$features_path = paste("visium/", gene, ".graph_features.tsv", sep='') # ENSMUSG00000051951  Xkr4  Gene  Expression, ...
  visium$matrix_path = paste("visium/", gene, ".graph_matrix.mtx", sep='') # as below
  # %%MatrixMarket matrix coordinate integer general
  # %metadata_json: {"software_version": "spaceranger-1.1.0", "format_version": 2}
  # 32285 2695 15850802 # genes barcodes valuesCnt
  # 9 1 2
  # 11 1 2
  visium$positions_file = file(visium$positions_path, 'w')
  visium$barcodes_file = file(visium$barcodes_path, 'w')
  visium$features_file = file(visium$features_path, 'w')
  visium$matrix_file = file(visium$matrix_path, 'w')
  return(visium)
}

############
##  MAIN  ##
############

# Read test data fig1 graph A/B/C 
setwd('~/BSC')
gene1 = "A"
gene2 = "B"
data1 = read.table(paste('~/local/lees2001/', gene1, '.graph.txt', sep=''), 
                   na.strings=".", sep=" ", row.names=1, header=T)
data2 = read.table(paste('~/local/lees2001/', gene2, '.graph.txt', sep=''),
                   na.strings=".", sep=" ", row.names=1, header=T)
visium1 = create_visium_list(gene='A')
visium2 = create_visium_list(gene='B')

write_visium_data(data1, visium1)
write_visium_data(data2, visium2)

# Read positions_list.csv
positions_colnames = c("in_tissue", "rix", "cix", "prix", "pcix")
con_graph = read.csv(positions_path, header=F, row.names=1)
colnames(con_graph) = positions_colnames
con_graph = subset(con_graph, in_tissue == 1)

