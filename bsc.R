library(Matrix)
library(dplyr)

# Create connectivity matrix -> C
create_connectivity_matrix = function(sdata) { # input Seurat spatial data
  # Seurat - sdata@images$anterior1@coordinates format:
  #                    tissue row col imagerow imagecol
  # AAACAAGTATCTCCCA-1      1  50 102     7474     8500
  # AAACACCAATAACTGC-1      1  59  19     8552     2788
  # AAACAGAGCGACTCCT-1      1  14  94     3163     7950
  # AAACAGCTTTCAGAAG-1      1  43   9     6636     2100
  positions_in_tissue = sdata@images$anterior1@coordinates %>% filter(tissue==1) 
  barcodes_in_tissue = rownames(positions_in_tissue)
  nbarcodes_in_tissue = length(barcodes_in_tissue)
  
  C = Matrix(nrow=nbarcodes_in_tissue, ncol=nbarcodes_in_tissue, data=0, sparse=T) # C: (sparse) connectivity matrix
  rownames(C) = barcodes_in_tissue
  colnames(C) = barcodes_in_tissue
  
  for (barcode in barcodes_in_tissue) {
    conn_mat = list() # data holder for connectivity matrices
    row_i = positions_in_tissue[barcode, 'row']
    col_i = positions_in_tissue[barcode, 'col'] # sdata$counts[['A']][[barcode]] later
    # Set nearby nodes and check connectivity (only in_tissue)
    neighbors = subset(sdata@images$anterior1@coordinates,
                       tissue==1 &
                         (
                           ((row==row_i-1) & (col==col_i-1)) |
                             ((row==row_i-1) & (col==col_i+1)) |
                             ((row==row_i) & (col==col_i-2)) |
                             ((row==row_i) & (col==col_i+2)) |
                             ((row==row_i+1) & (col==col_i-1)) |
                             ((row==row_i+1) & (col==col_i+1))
                         ))
    neighbor_barcodes = rownames(neighbors)
    if (length(neighbor_barcodes) > 0) C[barcode, neighbor_barcodes] = 1
  }
  
  W = C / rowSums(C) # W: weighted connectivity matrix
  W[is.na(W)] = 0
  
  conn_mat[['barcodes_in_tissue']] = barcodes_in_tissue
  conn_mat[['nbarcodes_in_tissue']] = nbarcodes_in_tissue
  conn_mat[['W']] = W
  conn_mat[['C']] = C
  
  return(conn_mat)
}

calculate_bsc = function(feature1, feature2, sdata, conn_mat, assay='SCT') {
  # Set X, Y variables from two features
  # Seurat - sdata@assays$SCT@data: row=feature, col=cell matrix
  # feature1 = 'Gpr88' # debug
  # feature2 = 'Penk' # debug
  mrna = sdata@assays[[assay]]@data # default to assay='SCT'
  X_values = mrna[feature1, conn_mat$barcodes_in_tissue] # 1:X
  Y_values = mrna[feature2, conn_mat$barcodes_in_tissue] # 2:Y
  X = data.frame(value=X_values)
  Xmean = mean(X$value)
  Y = data.frame(value=Y_values)
  Ymean = mean(Y$value)
  
  # Smoothened values
  X['smooth'] = conn_mat$W %*% X[,'value'] # Xsm = W * X
  Y['smooth'] = conn_mat$W %*% Y[,'value'] # Ysm = W * Y
  Xmean_sm = mean(X$smooth) # muX
  Ymean_sm = mean(Y$smooth) # muY
  
  # Calculate Peason's r(X,Y), r(smooth), L_XX, L_YY, L_XY as in Lee S (2001)
  r = (sum((X$value - Xmean) * (Y$value - Ymean))
       / (sqrt(sum((X$value - Xmean)^2)) * sqrt(sum((Y$value - Ymean)^2))))
  r_sm = (sum((X$smooth - Xmean_sm) * (Y$smooth - Ymean_sm))
          / (sqrt(sum((X$smooth - Xmean_sm)^2)) * sqrt(sum((Y$smooth - Ymean_sm)^2))))
  L_XX = (sum((X$smooth - Xmean)^2) / sum((X$value - Xmean)^2))
  L_YY = (sum((Y$smooth - Ymean)^2) / sum((Y$value - Ymean)^2))
  L_XY = sqrt(L_XX) * sqrt(L_YY) * r_sm
  
  # Group into Bivariate Spatial Correlation values list
  bsc = list('r' = r, 'r_sm' = r_sm,
             'L_XX' = L_XX, 'L_YY' = L_YY, 'L_XY' = L_XY)

  return(bsc)
}