#### ggplot2 colors : How to change colors automatically and manually?
http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
```
# Use a single color
# box plot
ggplot(ToothGrowth, aes(x=dose, y=len)) +
  geom_boxplot(fill='#A4A4A4', color="darkred")
# scatter plot
ggplot(mtcars, aes(x=wt, y=mpg)) + 
  geom_point(color='darkblue')
```
#### Change colors by groups
```
#Default colors
#The following R code changes the color of the graph by the levels of dose :
# Box plot
bp<-ggplot(ToothGrowth, aes(x=dose, y=len, fill=dose)) +
  geom_boxplot()
bp
# Scatter plot
sp<-ggplot(mtcars, aes(x=wt, y=mpg, color=cyl)) + geom_point()
sp
### The lightness (l) and the chroma (c, intensity of color) of the default (hue) colors can be modified using the functions scale_hue as follow :
# Box plot
bp + scale_fill_hue(l=40, c=35)
# Scatter plot
sp + scale_color_hue(l=40, c=35)
```
####  Change colors manually
A custom color palettes can be specified using the functions :

scale_fill_manual() for box plot, bar plot, violin plot, etc
scale_color_manual() for lines and points
```
# Box plot
bp + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# Scatter plot
sp + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
### Note that, the argument breaks can be used to control the appearance of the legend. This holds true also for the other scale_xx() functions.
# Box plot
bp + scale_fill_manual(breaks = c("2", "1", "0.5"), 
                       values=c("red", "blue", "green"))
# Scatter plot
sp + scale_color_manual(breaks = c("8", "6", "4"),
                        values=c("red", "blue", "green"))
```
#### Use RColorBrewer palettes
The color palettes available in the RColorBrewer package are described here : color in R.
```
# Box plot
bp + scale_fill_brewer(palette="Dark2")
# Scatter plot
sp + scale_color_brewer(palette="Dark2")
```

 
