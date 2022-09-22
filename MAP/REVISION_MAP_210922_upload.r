###########O'Donnell et. al. (2022) Publication Maps######

####load packages####
library(ggmap)
library(ozmaps)
library(sf)
library(ggspatial)
library(ggrepel)
library(ggplot2)

###Import Data####
path <- ""
setwd(path)
df <- read.csv('prostanthera_study_group_AVH_210922.csv', header = TRUE, stringsAsFactors = TRUE) #AVH occurrences
cities <- read.csv('cap_cities.csv') #city labels w lat longs
pops <- read.csv('phylicifolia_pops.csv') # populations of P. phylicifolia
morph_pops <- read.csv('morphological2.csv', header = TRUE) # populations sampled for morphological analysis
dartpops <- read.csv('dartpops.csv') # populations sampled for morphological analysis
dartpop_labs <- read.csv('dartpop_loc_labels.csv')
cap_cities <- read.csv('cap_cities_only.csv')

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
                         pad_x = unit(0.81, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) + 
  xlab("") + ylab("") + # remove latitude and longitude labels
  geom_point(dartpops, mapping=aes(x=long, y=lat,colour=species), size=8, alpha=0.3) + #add dart pops
  geom_point(morph_pops, mapping=aes(x=Longitude, y=Latitude),shape=3, size=2.5) + #add morph pops
  geom_label_repel(dartpop_labs, mapping=aes(x=long, y=lat, label=city), 
                   size=3, fill=alpha(c("grey"),0.5), fontface = "bold",
                   box.padding = unit(1,"lines"), point.padding = unit(0.3, "lines"),segment.size=0.7) + ## add dartpop labels
  geom_point(cities, mapping=aes(x=long, y=lat),size=1.5) + #add city points
  geom_text_repel(cap_cities, mapping = aes(x=long, y=lat, label=city), size=3) + #add city labels
  theme(legend.position = c(0.81, 0.22), legend.title = element_blank(), ### move legend
        legend.key.size = unit(0.15, "cm"), legend.text=element_text(size=8.5, face="italic")) +
  theme(legend.background = element_rect(fill="white",size=0.5, 
                                         linetype="solid", colour ="black")) ##legend fill and border
