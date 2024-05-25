#Mapping Simulated data
set.seed(1223)
source("Simulation.R")

SimulatedData<- SimulationResults[[1]]

#OverallCasesDF<- data.frame(Cases=rowSums(SimulatedData), AreaCode=Westmidlands.shp$LAD23CD, id=c(1:ndept), AreaName=Westmidlands.shp$LAD23NM, months=c(1:ndept))
OverallCasesDF<- data.frame(Cases=rowSums(SimulatedData), nuts3=France.shp$nuts3, id=c(1:ndept), AreaName=France.shp$nom, months=c(1:ndept))
#NWD<- as.vector(t(SimulatedData))
#cumulativeNWD<- t(apply((SimulatedData), MARGIN=1, FUN=cumsum))
#cumulativeNWD<- as.vector(t(cumulativeNWD))
#OverallCasesDF<- data.frame(Cases=cumulativeNWD, AreaCode=rep(Westmidlands.shp$LAD23CD, each=time), month=rep(1:time,ndept), AreaName=rep(Westmidlands.shp$LAD23NM, each=time), id=c(1:(ndept*time)))

#shapefile<- readOGR("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/LAD_MAY_2023_UK_BUC_V2.shp")
#shapefile <- fortify(shapefile, region = 'LAD23CD')
shapefile<- readOGR("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/ModifiedFranceShapefile.shp")
shapefile <- fortify(shapefile, region = 'nuts3')


#shp_OverallCasesDF<- sp::merge(shapefile,OverallCasesDF,
#                               by.x="id",
#                               by.y="AreaCode",
#                               all.x=F)
shp_OverallCasesDF<- sp::merge(shapefile,OverallCasesDF,
                               by.x="id",
                               by.y="nuts3",
                               all.x=F)
shp_OverallCasesDF <- arrange(shp_OverallCasesDF, order)

SimulatedData.map <- ggplot(shp_OverallCasesDF, aes(long, lat, fill = Cases, group = group)) +
  geom_polygon(col = "white") +
  scale_fill_viridis(option = "B", direction = -1) +  # Reverse the color scale
  coord_equal() + theme_void() +
 # ggtitle('Total simulated cases in 60 months') +
  labs(fill = "Prevalence")
 # labs(title = 'Month: {frame_time}', fill = "Cumulative Cases") +
  #transition_time(month) +
  #ease_aes('linear')

#animate(SimulatedData.map, duration = 60)
#anim_save('plot_SimulatedData.map.gif')

SimulatedData.map
