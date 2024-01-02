#Mapping Simulated data
source("Simulation.R")

#OverallCasesDF<- data.frame(OverallCases=rowSums(SimulatedData), AreaCode=Westmidlands.shp$LAD23CD, id=c(1:30), AreaName=Westmidlands.shp$LAD23NM, months=c(1:30))
#NWD<- as.vector(t(SimulatedData))
cumulativeNWD<- t(apply((SimulatedData), MARGIN=1, FUN=cumsum))
cumulativeNWD<- as.vector(t(cumulativeNWD))
OverallCasesDF<- data.frame(Cases=cumulativeNWD, AreaCode=rep(Westmidlands.shp$LAD23CD, each=time), month=rep(1:time,ndept), AreaName=rep(Westmidlands.shp$LAD23NM, each=time), id=c(1:(ndept*time)))

shapefile<- readOGR("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 2/Shapefile/LAD_MAY_2023_UK_BUC_V2.shp")
shapefile <- fortify(shapefile, region = 'LAD23CD')

shp_OverallCasesDF<- sp::merge(shapefile,OverallCasesDF,
                               by.x="id",
                               by.y="AreaCode",
                               all.x=F)
shp_OverallCasesDF <- arrange(shp_OverallCasesDF, order)

SimulatedData.map <- ggplot(shp_OverallCasesDF, aes(long, lat, fill = Cases, group = group)) +
  geom_polygon(col = "white") +
  scale_fill_viridis(option = "B", direction = -1) +  # Reverse the color scale
  coord_equal() + theme_void() +
  ggtitle('Total simulated cases in 24 months') +
  #labs(fill = "Cases")
  labs(title = 'Month: {frame_time}', fill = "Cumulative Cases") +
  transition_time(month) +
  ease_aes('linear')

animate(SimulatedData.map, duration = 60)
#anim_save('plot_SimulatedData.map.gif')

#SimulatedData.map
