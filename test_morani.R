library(RColorBrewer) #to use brewer.pal
library(fields) #to use designer.colors

# Read test data A
setwd('~/local')
data = read.table('./lees2001/A.graph.txt', na.strings = ".", sep = " ",
                  row.names = 1, header=T)

# Write data to Visium tissue_positions_list.csv format

#Initiate the plot window but do show any axes or points on the plot
plot(0, 0, type = "n", axes = FALSE, xlim=c(0, ncol(data)),
     ylim=c(0, nrow(data)), xlab="", ylab= "", asp=1)

for (i in 0:(nrow(data)-1)) {
  for (j in 0:(ncol(data)-1)) {
    value = data[i+1, j+1]
    if (is.na(value)) { next }
    if (value == 0) { in_tissue = 0 }
    else { in_tissue = 1 }
    
    barcode = sprintf("Ai%d_j%d", i, j)
    xcoord = 10*(i+1)
    ycoord = 10*(j+1)
    
    print_str = sprintf("%s,%d,%d,%d,%d,%d", barcode, in_tissue, 
                               i, j, 10*i, 10*j)
    print(print_str)
    
    Hexagon(i+1, j+1, col="black")
  }
}


#Function to create the polygon for each hexagon
Hexagon <- function (x, y, unitcell = 1, col = col) {
  polygon(c(x, x, x + unitcell/2, x + unitcell, x + unitcell,
            x + unitcell/2), c(y + unitcell * 0.125,
                               y + unitcell * 0.875,
                               y + unitcell * 1.125,
                               y + unitcell * 0.875,
                               y + unitcell * 0.125,
                               y - unitcell * 0.125),
          col = col, border=NA)
}#function

#Actual plotting of hexagonal polygons on map
offset <- 0.5 #offset for the hexagons when moving up a row
for (row in 1:SOM_Rows) {
  for (column in 0:(SOM_Columns - 1))
    # Hexagon(column + offset, row - 1, col = ColorCode[row + SOM_Rows * column])
  offset <- ifelse(offset, 0, 0.5)
}