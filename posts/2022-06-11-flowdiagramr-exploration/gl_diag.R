## load libraries ----
library(ggplot2)
library(flowdiagramr)


variables <- data.frame(
  id = 1:14,
  name = c("GL", "P", "PTH", "PTC", "PJR", "C", "CTB", "CTS", "CJR", "B", "BTH", "BTS", "BJR", "BE"),
  xmin = c(4, 4, 0, 4, 8, 8, 4, 8, 12, 12, 8, 12, 16, 16),
  xmax = c(7, 7, 3, 7, 11, 11, 7, 11, 15, 15, 11, 15, 19, 19),
  ymin = c(14, 12, 10, 10, 10, 8, 6, 6, 6, 4, 2, 2, 2, 0),
  ymax = c(15, 13, 11, 11, 11, 9, 7, 7, 7, 5, 3, 3, 3, 1),
  xlabel = c(5.5, 5.5, 1.5, 5.5, 9.5, 9.5, 5.5, 9.5, 13.5, 13.5, 9.5, 13.5, 17.5, 17.5),
  ylabel = c(14.5, 12.5, 10.5, 10.5, 10.5, 8.5, 6.5, 6.5, 6.5, 4.5, 2.5, 2.5, 2.5, 0.5),
  outline_color = c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"),
  fill_color = c("#b59dac", "#D9AF6B", "#ACACAC", "#ACACAC", "#ACACAC", "#D9AF6B", "#ACACAC", "#ACACAC", "#ACACAC", "#D9AF6B", "#ACACAC", "#ACACAC", "#ACACAC", "#b59dac"),
  label_text = c("Goldilocks", "Porridge", "Too hot", "Too cold", "Just right", "Chairs", "Too big", "Too small", "Just right", "Beds", "Too hard", "Too soft", "Just right", "Bears!"),
  label_color = c("#585c45", "#585c45", "#585c45", "#585c45", "#585c45", "#585c45", "#585c45", "#585c45", "#585c45", "#585c45", "#585c45", "#585c45", "#585c45", "#585c45"),
  label_size = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4)
)

flows <- data.frame(
  name = c("m_k10B", "m_k11B", "m_k12B", "m_k13BJR", "m_k1GL", "m_k2P", "m_k3P", "m_k4P", "m_k5PJR", "m_k6C", "m_k7C", "m_k8C", "m_k9CJR"),
  type = c("main", "main", "main", "main", "main", "main", "main", "main", "main", "main", "main", "main", "main"),
  id = c(10L, 11L, 12L, 13L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
  from = c("B", "B", "B", "BJR", "GL", "P", "P", "P", "PJR", "C", "C", "C", "CJR"),
  to = c("BTH", "BTS", "BJR", "BE", "P", "PTH", "PTC", "PJR", "C", "CTB", "CTS", "CJR", "B"),
  xstart = c(12, 13.5, 15, 17.5, 5.5, 4, 5.5, 7, 9.5, 8, 9.5, 11, 13.5),
  xend = c(11, 13.5, 16, 17.5, 5.5, 3, 5.5, 8, 9.5, 7, 9.5, 12, 13.5),
  ystart = c(4.5, 4, 4.5, 2, 14, 12.5, 12, 12.5, 10, 8.5, 8, 8.5, 6),
  yend = c(2.5, 3, 2.5, 1, 13, 10.5, 11, 10.5, 9, 6.5, 7, 6.5, 5),
  xlabel = c(11.5, 13.25, 15.5, 17.25, 5.25, 3.5, 5.25, 7.5, 9.25, 7.5, 9.25, 11.5, 13.25),
  ylabel = c(3.6, 3.5, 3.6, 1.5, 13.5, 11.6, 11.5, 11.6, 9.5, 7.6, 7.5, 7.6, 5.5),
  curvature = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  line_color = c("grey25", "grey25", "grey25", "grey25", "grey25", "grey25", "grey25", "grey25", "grey25", "grey25", "grey25", "grey25", "grey25"),
  line_size = c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7),
  line_type = c("solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid"),
  label_text = c("k10*B", "k11*B", "k12*B", "k13*BJR", "k1*GL", "k2*P", "k3*P", "k4*P", "k5*PJR", "k6*C", "k7*C", "k8*C", "k9*CJR"),
  label_color = c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"),
  label_size = c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
  show_label = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  arrow_size = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25),
  show_arrow = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
)


 ## ggplot2 code ----
###
# make the diagram with ggplot2
###
# Start with an empty ggplot2 canvas. The coord_equal function ensures
# that the x and y coordinates are displayed in equal proportions to
# on another (that is, it makes sure that the squares look like squares).
# All layers are added sequentially onto this blank canvas.
diagram_plot <- ggplot() +
  coord_equal(clip = "off")


# LAYER 1: STATE VARIABLES
# plot the states variable nodes as rectangles

# The variables data frame is used to create rectangles, with size determined
# by the xmin, xmax, ymin, and ymax values in the nodes data frame. The
# outline color of the rectangles is defined by var_outline_color; the
# inside color (fill) of the rectangles is defined by var_fill_color.
# The color variables can be a single value or a vector, giving different
# colors to different rectangles/nodes/state variables. If a vector, the
# color and fill vectors must have a length that is equal to the number
# of rows in the nodes data frame (one value for each row).

# create the nodes/boxes/variables
# these are just empty rectangles with no text
for(i in 1:nrow(variables)) {
  diagram_plot <- diagram_plot +  # add new stuff to blank canvas
    geom_rect(
      data = variables[i, ],  # one row of the data frame
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),  # location information
      color = variables[i, "outline_color"],  # border color
      fill = variables[i, "fill_color"]  # internal, fill color
    )
}

# add label text, which goes on top of boxes based on location information
for(i in 1:nrow(variables)) {
  diagram_plot <- diagram_plot +  # add text to boxes
    geom_text(
      data = variables[i, ],
      aes(x = xlabel, y = ylabel, label = label_text),
      size = variables[i, "label_size"],
      color = variables[i, "label_color"]
    )
}

## add in all the flows
# start with the lines/arrows
for(i in 1:nrow(flows)) {
  if(flows[i, "show_arrow"] == TRUE) {
    diagram_plot <- diagram_plot +  # add the lines to the plot with boxes
      geom_curve(  # always use geom_curve, which is straight when cuvature = 1
        data = flows[i, ],
        aes(x = xstart,
            y = ystart,
            xend = xend,
            yend = yend),
        linetype = flows[i, "line_type"],
        arrow = arrow(length = unit(flows[i, "arrow_size"],"cm"), type = "closed"),
        color = flows[i, "line_color"],
        arrow.fill = flows[i, "line_color"],
        lineend = "round",
        size = flows[i, "line_size"],
        curvature = flows[i, "curvature"],
        ncp = 1000  # controls smoothness of curve, larger number = more smooth
      )
  }
}

for(i in 1:nrow(flows)) {
  if(flows[i, "show_label"] == TRUE) {
    diagram_plot <- diagram_plot +  # now add the flow labels to the canvas
      geom_text(
        data = flows[i, ],
        aes(x = xlabel, y = ylabel, label = label_text),
        size = flows[i, "label_size"],
        color = flows[i, "label_color"])
  }
}

# If with_grid == FALSE (default) then void out the theme
# otherwise keep the grey background with grid
# the grid can be useful for updating positions of items
with_grid <- FALSE  # default is false
if(with_grid == FALSE) {
  diagram_plot <- diagram_plot +
    theme_void()  # makes an empty plot theme with no axes, grids, or ticks
} else {
  # The else here may seem silly, but otherwise the returned plot is NULL
  diagram_plot <- diagram_plot  # just returns default ggplot2 theme
}
  


# These lines plot or save the generated diagram. 
# Uncomment them if you want to perform either action. 
# plot(diagram_plot) 
# ggsave('diagram_plot.png',diagram_plot)