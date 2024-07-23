
ggplot(background_df,aes(x1,x2)) +
  geom_tile(aes(),
            alpha=0.65)


ggplot(background_df,aes(x1,x2)) +
  geom_tile(aes(fill=z),
            alpha=0.65)+
  scale_fill_gradient(low="blue", high="orange") +
  geom_contour(aes(z=z),
               bins=5)

# what is the difference between color and fill?
# color: outside is color
# fill: inside is fill

# A key things to learn seems to how to handle 3 variables
#   with x,y, and a z variable
