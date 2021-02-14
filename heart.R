## Import required packages ##
library(ggplot2)
library(dplyr)

## Create data frame with values for both parametric equations ##
heartdf = data_frame(
  t = seq(0, 2*pi, pi/60), 
  x = 16*sin(t)^3, #x values for heart curve
  y = 13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t), #y values for heart curve
  x2 = 3*sqrt(2)*cos(t)/(sin(t)^2 + 1), #x values for leminiscate
  y2 = 20 + 4*sqrt(2)*cos(t)*sin(t)/(sin(t)^2 + 1) #y values for lemniscate
)

## Create plot ##
p = ggplot(data = heartdf, aes(x, y)) + 
  geom_path(aes(group = 1)) +
  geom_path(aes(group = 1, x = x2, y= y2), size = 2, colour = "red") +
  geom_polygon(aes(group = 1), fill = "red") +
  geom_text(aes(x = 0, y = 0, label = "Happy Valentine's Day \n Marioly!!!"), 
            size = 10, colour = "white")

p