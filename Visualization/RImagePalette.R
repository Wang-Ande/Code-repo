# packages
library(pak)
pak("RImagePalette")
library(RImagePalette)
#------------------------ run Example ------------------------------------------
#Load an image
# 使用 JPEG 图像示例
lifeAquatic <- jpeg::readJPEG("your_image.jpg")
# 或 PNG
lifeAquatic <- png::readPNG("your_image.png")
display_image(lifeAquatic)

# Create a palette of 9 colors
palette <- image_palette(lifeAquatic, n = 9)

#------------------------ your data --------------------------------------------
lifeAquatic <- jpeg::readJPEG("火星图标大.jpg")
display_image(lifeAquatic)
palette <- image_palette(lifeAquatic, n = 5)
scales::show_col(palette)

lifeAquaticPalette <- image_palette(lifeAquatic, n=9, choice=median, volume=TRUE)
scales::show_col(lifeAquaticPalette)

library(ggplot2)
#Create plot
p <- ggplot(data = iris, aes(x=Species, y=Sepal.Width, fill=Species)) + geom_bar(stat="identity")
#Apply scale
p + theme_bw() + scale_fill_manual(values=lifeAquaticPalette[c(3,5,2)])
