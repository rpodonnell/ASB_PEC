###########O'Donnell et. al. (2022) Publication Maps######

####load packages####
library(ggmap)
library(ozmaps)
library(sf)
library(viridis)
library(ggspatial)
library(ggrepel)
library(ggplot2)

###Import Data####
path <- ""
setwd(path)
df <- read.csv('prostanthera_study_group_AVH_210922.csv', header = TRUE, stringsAsFactors = TRUE) #AVH occurrences
cities <- read.csv('cities2.csv') #city labels w lat longs
pops <- read.csv('phylicifolia_pops.csv') # populations of P. phylicifolia
morph_pops <- read.csv('morphological.csv', header = TRUE) # populations sampled for morphological analysis
dartpops <- read.csv('dartpops.csv') # populations sampled for morphological analysis

###generate map of Australia with states
oz_states <- ozmaps::ozmap_states
oz_states

##ggmapping
ggplot(oz_states) + 
  geom_sf(fill = "white") + #set state fill 
  theme_bw() + #set theme
  coord_sf(xlim=c(145, 152.3), ylim=c(-39.3, -33.3)) + #set lat long frame
  geom_point(df, mapping = aes(x=Longitude, y=Latitude, 
             group=Scientific.Name, colour=Scientific.Name), size=2) + #plot observation points
  scale_color_brewer(palette="Set2") + #set point colours
  annotation_scale(location = "br", width_hint = 0.3) + #add scale bar
  annotation_north_arrow(location = "br", which_north = "true", #add north arrow 
             pad_x = unit(0.6, "in"), pad_y = unit(0.3, "in"),
             style = north_arrow_fancy_orienteering) + 
  xlab("") + ylab("") + # remove latitude and longitude labels
  geom_point(dartpops, mapping=aes(x=long, y=lat,colour=species), size=8, alpha=0.3) +
  geom_label_repel(pops, mapping=aes(x=Longitude, y=Latitude, label=Population..), size=2.5, fill="grey", fontface = "bold")
  #add dart pops
  geom_text_repel(data = cities, aes(x=long, y=lat, label=city), size=3) +
  geom_point(cities, mapping=aes(x=long, y=lat), size=2) 

+
  theme(legend.position = c(0.81, 0.25), legend.title = element_blank(),legend.key.size = unit(0.15, "cm"), legend.text=element_text(size=7.5)) +
  theme(legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))

theme_

theme_bw
###SHAPES

shapes <- data.frame(
  shape = c(0:19, 22, 21, 24, 23, 20),
  x = 0:24 %/% 5,
  y = -(0:24 %% 5)
)
ggplot(shapes, aes(x, y)) + 
  geom_point(aes(shape = shape), size = 5, fill = "red") +
  geom_text(aes(label = shape), hjust = 0, nudge_x = 0.15) +
  scale_shape_identity() +
  expand_limits(x = 4.1) +
  theme_void()



theme_set(theme_bw())
ggplot(oz_states) + 
  geom_sf(fill = "white") +
  coord_sf(xlim=c(142.44, 154.6), ylim=c(-38.85, -22.5)) +
  geom_point(df, mapping = aes(x=Longitude, y=Latitude, group=Scientific.Name, shape=Scientific.Name, colour = Scientific.Name), size=3) +
  #   scale_shape_manual(values=c(21:24)) +
  #   # scale_color_viridis(discrete = TRUE) +
  # scale_fill_viridis(discrete = TRUE) +
  # scale_color_viridis_b(discrete=TRUE, option="D") +
  scale_color_brewer(palette="Set2") +
  
  annotation_scale(location = "br", width_hint = 0.3) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + 
  xlab("") + ylab("") +
  geom_point(cities, mapping=aes(x=long, y=lat), size=1)+
  geom_label_repel(pops, mapping=aes(x=Longitude, y=Latitude, label=Population..), size=2.5, fill="grey", fontface = "bold") +
  theme(legend.position = c(0.2, 0.089), legend.title = element_blank(),legend.key.size = unit(0.15, "cm"), legend.text=element_text(size=7.5)) +
  theme(legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"))
