#Required packages
library(fingertipsR)
library(fingertipscharts)
require(rgdal)
library(sp)
library(sf)
library(tidyverse)
library(ggplot2)
library(maps)
library(mapdata)
library(maptools)
library(rgdal)
library(ggmap)
library(rgeos)
library(broom)
library(plyr)
library(dplyr)
library(viridis)

#Overall Uk TB incidence data
#inds <- select_indicators()
indid <- 91359
indicators()%>% filter(IndicatorID==91359)
TBIncidence <- fingertips_data(IndicatorID = indid, AreaTypeID = "All")
#head(TBIncidence)

ggplot(TBIncidence, aes(x = Timeperiod, y = Value, group = 1)) +
  geom_line() +
  labs(title = "UK Tuberculosis Incidence per 100,000 population",
       x = "Year",
       y = "Value")+
theme(plot.title = element_text(hjust = 0.5))

#UK TB incidence by regions
indid <- 91361
indicators()%>% filter(IndicatorID==91361)
df <- fingertips_data(IndicatorID = indid, AreaTypeID = "All")
#head(df)


#Map
map_data<- df %>%
  group_by(IndicatorID) %>%
  filter(IndicatorID==91361,
                TimeperiodSortable==max(TimeperiodSortable)) %>%
  ungroup() %>%
  mutate(ComparedtoEnglandvalueorpercentiles=factor(ComparedtoEnglandvalueorpercentiles, levels=c("Better","Similar","Worse","Not compared")))

ons_api <- "https://services1.arcgis.com/ESMARspQHYMw9BZ9/arcgis/rest/services/Counties_and_Unitary_Authorities_December_2021_EN_BFC_2022/FeatureServer/0/query?outFields=*&where=1%3D1&f=geojson"
  
TBmap<- fingertipscharts::map(data=map_data,
                          ons_api=ons_api,
                          area_code=AreaCode,
                          fill=ComparedtoEnglandvalueorpercentiles,
                          title="UK TB incidence (three year average)",
                          )
TBmap

TBmapInt<- fingertipscharts::map(data=map_data,
                              ons_api=ons_api,
                              area_code=AreaCode,
                              fill=ComparedtoEnglandvalueorpercentiles,
                              type="interactive",
                              value=Value,
                              name_for_label=AreaName,
                              title="UK TB incidence (three year average)",
)
TBmapInt

######################################################################################################

#UK TB incidence by regions
indid <- 91361
indicators()%>% filter(IndicatorID==91361)
df <- fingertips_data(IndicatorID = indid, AreaTypeID = "All")
#head(df)
head(df[df$AreaCode=="E06000001",])

shapefile<- readOGR("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/CTYUA_MAY_2023_UK_BUC.shp")
shapefile <- fortify(shapefile, region = 'CTYUA23CD')

shp_df<- sp::merge(shapefile,df,
                   by.x="id",
                   by.y="AreaCode",
                   all.x=F)

shp_df <- arrange(shp_df, order)

map <- ggplot(shp_df, aes(long, lat, fill = Value, group = group)) +
  geom_polygon(col = "white") +
  scale_fill_viridis(option = "D", direction = -1) +  # Reverse the color scale
  coord_equal() + theme_void() +
  ggtitle('UK TB incidence (three year average 2000 to 2021)') +
  labs(fill = "Tuberculosis rate per 100,000")

map