#making a hex sticker for myself

library('hexSticker')
library('ggplot2')
library('here')
library('DSAIDE')
library('tidyr')

#make inside graph for sticker
#use DSAIDE for this
modelsettings =  list(fitmodel = 3, iter = 1000, solvertype = 3, S = 288, b = 0.08, g = 1.6,n=1,t1=10,t2=12)
modelsettings$modeltype = "_fit_"
modelsettings$nplots = 1
modelsettings$simfunction = 'simulate_fit_noro'
result = run_model(modelsettings)
p1 <- generate_ggplot(result) + theme_void() + theme(legend.position="none")
plot(p1)

#building sticker myself using the individual building block functions from hexSticker
mysticker <- ggplot() + geom_hexagon(size = 2, fill = "#0099cc", color = "#0099cc") +
                      geom_subview(subview = p1, x=1, y=1, width=1.3, height=1) +
                      geom_pkgname(package="Andreas\n Handel", x = 1, y = 1, color = "black", family = "sans", size=24) +
                      geom_url(url = 'www.andreashandel.com',family = 'sans', size = 8) +
                      theme_sticker()
plot(mysticker)
ggsave(filename=here('assets/images',"icon.png"), plot = mysticker, family = "sans", width = 200, height = 200, units = "mm")
#ggsave(filename=here('static/img',"icon-32.png"), plot = mysticker, width = 43.9, height = 50.8, bg = 'transparent',  units = "mm", family = "sans")




#ggplot2 example
#p <- ggplot(aes(x = mpg, y = wt), data = mtcars) + geom_point()
#img <- p + theme_void() + theme_transparent()
#sticker(img, package="Andreas \n Handel", p_size=20, s_x=1, s_y=.75, s_width=1.3, s_height=1, filename= here('static/img',"icon.png"))
#img <- here('media',"fig.png")
#sticker(img=NULL, package="Andreas Handel", p_size=20, p_y = 1.3, s_x=1, s_y=.85, s_width=0.8, s_height=0.8, filename=here('media',"icon.png"), )
