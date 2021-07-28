library(hash)
# library(RColorBrewer) #to use brewer.pal
# library(fields) #to use designer.colors
source("~/BSC/test_morani.R") ##@## get Hexagon, for debugging

# Function that input test data format and write Visium-style positions file
write_positions_data = function(test_position_data, spatial=spatial) {
  data = test_position_data$test_data
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
      write(print_str, spatial$positions_file, append = T)
      
      # if (value == 0) color = "yellow"
      # else color = colors[value]
      # Hexagon(j+1, nrow(data)-(i+1), col=color) # x, y - draw hexagon
    }
  }
}

# Visium-style spatial file settings
Spatial = function(positions_path) {
  spatial = list()
  spatial$positions_path = positions_path
  spatial$positions_file = file(spatial$positions_path, 'w')
  return(spatial)
}

# Visium-style file settings
Visium = function(out_dir="visium/") {
  visium = list() # data holder
  visium$barcodes_path = paste(out_dir, "barcodes.tsv", sep='') # AAACAAGTATCTCCCA-1, ...
  visium$features_path = paste(out_dir, "features.tsv", sep='') # ENSMUSG00000051951  Xkr4  Gene  Expression, ...
  visium$matrix_path = paste(out_dir, "matrix.mtx", sep='') # as below
  # %%MatrixMarket matrix coordinate integer general
  # %metadata_json: {"software_version": "spaceranger-1.1.0", "format_version": 2}
  # 32285 2695 15850802 # genes barcodes valuesCnt
  # 9 1 2
  # 11 1 2
  visium$barcodes_file = file(visium$barcodes_path, 'w')
  visium$features_file = file(visium$features_path, 'w')
  visium$matrix_file = file(visium$matrix_path, 'w')
  return(visium)
}

# Data per gene, given as test data format
GeneTestData = function(gsym, in_dir="~/local/lees2001/") {
  gene = list()
  test_data_path = paste(in_dir, gsym, '.graph.txt', sep='')
  if (!file.exists(test_data_path)) write(paste("[ERROR]", test_data_path, "does not exist."), stderr())
  test_data = read.table(test_data_path, 
                         na.strings=".", sep=" ", row.names=1, header=T)
  return(test_data)
}

# Function that input test data format and write Visium-style output files
#  input: gsyms = ['A', 'B']
#  output: barcodes.tsv, features.tsv, matrix.mtx
write_visium_data = function(gsyms, visium=visium) {
  visium$barcode2bid = hash()
  visium$bid2barcode = hash()
  visium$feature2fid = hash()
  visium$fid2feature = hash()
  barcode_cnt = 0
  feature_cnt = 0
  
  # Set and save features (gene info)
  gene_test_data = list()
  for (gsym in gsyms) {
    gene_test_data[[gsym]] = GeneTestData(gsym)
  
    # ENSMUSG00000051951  Xkr4  Gene  Expression, ...
    gene_info = paste(gsym, gsym, "Gene", "Expression", sep="\t")
    write(gene_info, visium$features_file, append=T)
    feature_cnt = feature_cnt + 1
    visium$feature2fid[gsym] = feature_cnt
    visium$fid2feature[toString(feature_cnt)] = gsym
  }
  nrows = nrow(gene_test_data[[gsym]])
  ncols = ncol(gene_test_data[[gsym]])
  
  for (i in 0:(nrows-1)) { # i:y
    if (i %% 2 == 0) col_ixs = seq(0, ncols-1, 2)
    else col_ixs = seq(1, ncols-1, 2)
  
    for (j in col_ixs) { # j:x
  
      # Set and save barcode
      barcode = sprintf("i%d_j%d", i, j) # i:y, j:x
      # AAACAAGTATCTCCCA-1, ...
      write(barcode, visium$barcodes_file, append=T)
      barcode_cnt = barcode_cnt + 1
      barcode_id = toString(barcode_cnt)
      visium$barcode2bid[barcode] = barcode_cnt
      visium$bid2barcode[barcode_id] = barcode
  
      for (gsym in gsyms) {
        fid = visium$feature2fid[[gsym]] # gene ID
        value = gene_test_data[[gsym]][i+1, j+1]
        # print(value) ##@##
        if (is.na(value)) next
        else if (value == 0) in_tissue = 0
        else in_tissue = 1
  
        # %%MatrixMarket matrix coordinate integer general
        # %metadata_json: {"software_version": "spaceranger-1.1.0", "format_version": 2}
        # 32285 2695 15850802 # genes barcodes valuesCnt
        # 9 1 2
        # 11 1 2
        if (in_tissue == 0) next # discard line if not in tissue
        mtx_line = paste(fid, barcode_id, value, sep=" ")
        write(mtx_line, visium$matrix_file, append=T)
      }
    }
  }
}

read_visium_data = function(data_dir) {
  positions_path = file.path(data_dir, 'positions_list.csv')
  barcodes_path = file.path(data_dir, 'barcodes.tsv')
  features_path = file.path(data_dir, 'features.tsv')
  matrix_path = file.path(data_dir, 'matrix.mtx')
  
  positions = read.table(positions_path, sep=',', row.names=1)
  colnames(positions) = c('in_tissue', 'row', 'col', 'prow', 'pcol')
  
  barcodes = read.table(barcodes_path, sep='\t')[,1] # as vector
  features = read.table(features_path, sep='\t')[,1] # 1st col only
  read_cnts = read.table(matrix_path, sep=' ', col.names=c('feature', 'barcode', 'count'))

  data = list(
    spatial = positions,
    barcodes = barcodes, 
    features = features, 
    counts = list()
  )
  
  for (fix in 1 : (length(features))) {
    feature = features[fix] # gene id -> gene symbol
    data$counts[[feature]] = read_cnts[read_cnts$feature==fix,]$count # save counts for gene
    names(data$counts[[feature]]) = read_cnts[read_cnts$feature==fix,]$barcode # rename rows as barcode ids
    names(data$counts[[feature]]) = barcodes[as.integer(names(data$counts[[feature]]))] # rename rows as barcodes
  }
  
  return(data)
}



############
##  MAIN  ##
############

# Read test data fig1 graph A/B/C 
setwd('~/BSC')
test_position_data = Gene('A')

spatial = Spatial("visium/positions_list.csv")
write_positions_data(test_position_data, spatial)

visium = Visium()
gsyms = c('A', 'B')
write_visium_data(gsyms, visium)

# Read Visium-style data
vdata = read_visium_data(data_dir='./visium')


# # Read positions_list.csv
# positions_colnames = c("in_tissue", "rix", "cix", "prix", "pcix")
# con_graph = read.csv(positions_path, header=F, row.names=1)
# colnames(con_graph) = positions_colnames
# con_graph = subset(con_graph, in_tissue == 1)

