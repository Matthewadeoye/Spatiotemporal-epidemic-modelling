#USA and UK cities Measles biweekly dataset
measles_cases<- read.csv("measlesUKUS.csv")

measles_cases<- subset(measles_cases, select=c(decimalYear,country,cases,loc))
measles_cases<- measles_cases[measles_cases$country !="US",]
measles_cases<- measles_cases[measles_cases$loc ==c("LONDON","BIRMINGHAM","LIVERPOOL","MANCHESTER","LEEDS","SHEFFIELD"),]
table(measles_cases$loc)
names(measles_cases)[names(measles_cases) == 'loc'] <- 'city'


#load required packages quietly
packages<- c("spatPomp", "fingertipsR", "fingertipscharts", "sp", "sf", "ggplot2", "viridis", "rgdal", "tidyverse", "maps", "mapdata", "maptools", "mapview", "ggmap", "rgeos", "broom", "plyr", "dplyr", "viridis")
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))
}))

#UK Measles epidemic
data("measlesUK")
  
measlesUK<- measlesUK[measlesUK$city ==c("LONDON","BIRMINGHAM","LIVERPOOL","MANCHESTER","LEEDS","SHEFFIELD"),]
ggplot(measlesUK, aes(x = year, y = log(cases), colour = city)) +
  geom_line() +
  facet_wrap( ~ city)

measlesUK<- measlesUK[measlesUK$year >=1950,]

indid<- 93119
measles <- fingertips_data(IndicatorID = indid, AreaTypeID = "All")
shapefile<- readOGR("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/CTYUA_MAY_2023_UK_BUC.shp")
shapefile <- fortify(shapefile, region = 'CTYUA23CD')

shp_measles<- sp::merge(shapefile,measles,
                   by.x="id",
                   by.y="AreaCode",
                   all.x=F)

shp_measles <- arrange(shp_measles, order)
measles.map <- ggplot(shp_measles, aes(long, lat, fill = Value, group = group)) +
  geom_polygon(col = "white") +
  scale_fill_viridis(option = "D", direction = -1) +  # Reverse the color scale
  coord_equal() + theme_void() +
  ggtitle('Annual measles incidence from 2012 to 2021') +
  labs(fill = "Measles rate per 100,000")

measles.map

