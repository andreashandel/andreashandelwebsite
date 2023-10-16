# Code only portion of flowdiagramr exploration post
# https://www.andreashandel.com/posts/flowdiagramr-R-package/


## ---- packages --------
library('ggplot2')
library('flowdiagramr')
library('sysfonts') #for extra fonts
library('showtext') #for extra fonts

## ---- variables --------
# the unique internal labels, just abbreviations of the full names
variables = c("GL","P","PTH","PTC","PJR",
              "C","CTB","CTS","CJR",
              "B","BTH","BTS","BJR",
              "BE")
# the eventual names for all boxes
# we'll use that later to label the boxes
varnames = c("Goldilocks",
             "Porridge",
             "Too hot","Too cold","Just right",
             "Chairs",
             "Too big", "Too small", "Just right",
             "Beds",
             "Too hard","Too soft","Just right",
             "Bears!")
# assigning the varnames the variable labels
# needed later
names(varnames) <- variables
# a matrix specifying the locations on a matrix grid layout
# this mimics the look of the original blog post
varlocations = matrix(data = c("","GL", "", "", "",
                               "","P", "", "", "",
                               "PTH","PTC", "PJR", "", "",
                               "",   "",    "C", "", "",
                               "","CTB",  "CTS", "CJR", "",
                               "",   "",   "",   "B", "",
                               "",   "", "BTH", "BTS", "BJR",
                               "",   "", "",     "", "BE"
                      ), ncol = 5, byrow = TRUE)

## ---- flows --------
# setting up the inflows and outflows (arrows) for each box
flows = list( GL_flows = c("-k1*GL"),
              P_flows = c("k1*GL","-k2*P","-k3*P","-k4*P"),
              PTH_flows = c("k2*P"),
              PTC_flows = c("k3*P"),
              PJR_flows = c("k4*P","-k5*PJR"),
              C_flows = c("k5*PJR","-k6*C","-k7*C","-k8*C"),
              CTB_flows = c("k6*C"),
              CTS_flows = c("k7*C"),
              CJR_flows = c("k8*C","-k9*CJR"),
              B_flows = c("k9*CJR","-k10*B","-k11*B","-k12*B"),
              BTH_flows = c("k10*B"),
              BTS_flows = c("k11*B"),
              BJR_flows = c("k12*B","-k13*BJR"),
              BE_flows = c("k13*BJR")
)


## ---- prepare -------------
# model object
gl_model = list(variables = variables, flows = flows)
# model layout
model_settings = list(varlocations=varlocations,
                      varbox_x_size = 3)
# prepare model
gl_list = flowdiagramr::prepare_diagram(gl_model, model_settings)


## ---- update --------
#set colors that are similar to original blog post
varcolors = c("#b59dac", rep(c("#D9AF6B", "#ACACAC","#ACACAC","#ACACAC"),3), "#b59dac")
# make them a named vector since that's required by update_diagram
names(varcolors) = variables
# list of all style updates we want
diagram_settings = list(var_fill_color = varcolors,
                        var_label_text = varnames,
                        var_label_color = c(all =  "#585c45"),
                        flow_show_label = c(all = FALSE),
                        var_label_size = c(all = 4))

# update the look
gl_list2 <- flowdiagramr::update_diagram(gl_list,diagram_settings)

## ---- makediag --------
# create and plot diagram
gl_diag <- flowdiagramr::make_diagram(gl_list2)
plot(gl_diag)


## ---- updateggplot --------
# get different fonts
sysfonts::font_add_google(name = "Henny Penny", family = "henny")
showtext::showtext_auto()

# update the plot by adding ggplot2 commands
gl_diag2 <- gl_diag  +
       labs(title = "The Goldilocks Decision Tree",
            caption = "Made with flowdiagramr:\n https://andreashandel.github.io/flowdiagramr/") +
       theme_void() +
        theme(plot.margin = unit(c(1, 1, 0.5, 1), "cm"),
        legend.position = "none",
        plot.background = element_rect(colour = "#f2e4c1", fill = "#f2e4c1"),
        panel.background = element_rect(colour = "#f2e4c1", fill = "#f2e4c1"),
        plot.title = element_text(family = "henny", hjust = 0, face = "bold",
                                  size = 45, color = "#585c45",
                                  margin = margin(t = 10, r = 0, b = 10, l = 0)),
        plot.caption = element_text(family = "henny", hjust = 0,
                                    size = 16, color = "#585c45",
                                    margin = margin(t = 10)),
       text = element_text(family = "henny")
       )

## ---- plotupdate --------
plot(gl_diag2)


## ---- writecode -----
write_diagram(gl_list2,filename = "gl_diag.R", always_overwrite = TRUE)


## ---- newdiagcode -----
source("gl_diag_mod.R")
# plot new diagram
plot(diagram_plot)


## ---- flowsnew --------
# need to define again since the above file overwrote it
variables = c("GL","P","PTH","PTC","PJR",
              "C","CTB","CTS","CJR",
              "B","BTH","BTS","BJR",
              "BE")

# more complex flows
flowsnew = list( GL_flows = c("-k1*GL"),
              P_flows = c("k1*GL","-k2*P","-k3*P","-k4*P","-kk1*P*CJR"),
              PTH_flows = c("k2*P","k3a*PTC"),
              PTC_flows = c("k3*P","-k3a*PTC"),
              PJR_flows = c("k4*P","-k5*PJR","kk1*P*CJR"),
              C_flows = c("k5*PJR","-k6*C","-k7*C","-k8*C"),
              CTB_flows = c("k6*C"),
              CTS_flows = c("k7*C"),
              CJR_flows = c("k8*C","-k9*CJR"),
              B_flows = c("k9*CJR","-k10*B","-k11*B","-k12*B"),
              BTH_flows = c("k10*B"),
              BTS_flows = c("k11*B"),
              BJR_flows = c("k12*B","-k13*BJR"),
              BE_flows = c("k13*BJR")
)


## ---- preparemakenew -------------
# model object
gl_model_new = list(variables = variables, flows = flowsnew)
# model layout
model_settings = list(varlocations=varlocations,
                      varbox_x_size = 3)
# prepare model
gl_list_new = flowdiagramr::prepare_diagram(gl_model_new, model_settings)
# update the look
gl_list_new2 <- flowdiagramr::update_diagram(gl_list_new,diagram_settings)
# create and plot diagram
gl_diag_new <- flowdiagramr::make_diagram(gl_list_new2)
plot(gl_diag_new)
ggsave('featured.png',diagram_plot)




