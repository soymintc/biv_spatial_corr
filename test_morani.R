library(RColorBrewer) #to use brewer.pal
library(fields) #to use designer.colors

# Read test data A
setwd('~/local')
data = read.table('./lees2001/A.graph.txt', na.strings = ".", sep = " ",
                  row.names = 1, header=T)

# Color code
colors = c("white", "grey", "black")


#Function to create the polygon for each hexagon
Hexagon <- function (x, y, unitcell = 1, col = col) {
  x = x/2
  polygon(c(x, x, x + unitcell/2, x + unitcell, x + unitcell, x + unitcell/2), 
          c(y + unitcell * 0.125,
            y + unitcell * 0.875,
            y + unitcell * 1.125,
            y + unitcell * 0.875,
            y + unitcell * 0.125,
            y - unitcell * 0.125),
          col = col, border=NA)
} # hexagon fn

# Write data to Visium tissue_positions_list.csv format
plt = plot(0, 0, type = "n", axes = FALSE, xlim=c(0, (ncol(data))), # draw blank plot
           ylim=c(0, (1.2*nrow(data))), xlab="", ylab= "", asp=1)

for (i in 0:(nrow(data)-1)) { # y
  for (j in 0:(ncol(data)-1)) { # x
    value = data[i+1, j+1]
    if (is.na(value)) next
    if (value == 0) in_tissue = 0
    else { in_tissue = 1 }
    
    barcode = sprintf("Ai%d_j%d", i, j) # y, x
    xcoord = 10*(i+1) # y
    ycoord = 10*(j+1) # x
    
    print_str = sprintf("%s,%d,%d,%d,%d,%d", barcode, in_tissue, 
                               i, j, 10*i, 10*j)
    print(print_str)
    
    if (value == 0) color = "yellow"
    else color = colors[value]
    Hexagon(j+1, nrow(data)-(i+1), col=color) # x, y - draw hexagon
  }
}



