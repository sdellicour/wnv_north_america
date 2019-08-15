library(rgdal)
library(rgeos)
library(diagram)
library(adephylo)
library(seraphim)
library(lubridate)

	# 1. Plotting the different environmental rasters
	# 2. Extraction of the spatio-temporal information embedded in 100 trees
	# 3. Generating a dispersal history graph (by time slice or not)
		# 3.1. Building the maximum clade consensus (MCC) with TreeAnnotator
		# 3.2. Extraction of the spatio-temporal information embedded in the MCC tree
		# 3.3. Estimating annual kernel density polygons
		# 3.4. Plotting the global graph of dispersal history
		# 3.5. Plotting dispersal history graph snapshots
	# 4. Estimating and plotting different dispersal statistics
	# 5. Performing RRW simulations along trees (null model)
	# 6. Testing the impact of environmental factors on branch dispersal velocity
	# 7. Testing the impact of migratory flyways on branch dispersal frequency
		# 7.1. Testing the impact of the four US administrative flyways (NAMFs)
		# 7.2. Testing the impact of the three flyways estimated by La Sorte et al. (2014) for terrestrial birds

e_WNV = extent(-128, -63, 19, 55) # extent of study area
localTreesDirectory = "Tree_extractions/WNV4_gamma_100"
nberOfExtractionFiles = 100; mostRecentSamplingDatum = 2016.6475

# 1. Plotting the different environmental rasters

envVariableFiles = c("Elevation","Land_cover_forests","Land_cover_shrublands","Land_cover_savannas","Land_cover_grasslands",
					 "Land_cover_croplands","Land_cover_water","Annual_mean_temperature","Annual_precipitation")
envVariableNames1 = c("","","","","","","Water","Annual","Annual")
envVariableNames2 = c("Elevation","Forests","Shrublands","Savannas","Grasslands","Croplands","areas","mean temp.","precipitation")
rS = list(); cols = list()
for (i in 1:length(envVariableFiles)) rS[[i]] = raster(paste0("Environmental_files/WNV_rasters/",envVariableFiles[i],"_WNV_04.asc"))
rS[[1]][rS[[1]][]<0] = 0; rS[[9]][] = rS[[9]][]/100 # legend: temperature in °C and precipitation in meters
cols[[1]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(100)
colour1 = "white"; colour2 = "olivedrab3"; r = rS[[2]]
cols[[2]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
colour1 = "white"; colour2 = "burlywood3"; r = rS[[3]]
cols[[3]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
colour1 = "white"; colour2 = "brown"; r = rS[[4]]
cols[[4]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
colour1 = "white"; colour2 = "chartreuse4"; r = rS[[5]]
cols[[5]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
colour1 = "white"; colour2 = "navajowhite4"; r = rS[[6]] 
cols[[6]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
colour1 = "white"; colour2 = colorRampPalette(brewer.pal(9,"Blues"))(9)[5]
cols[[7]] = colorRampPalette(c(colour1,colour2),bias=1)((max(r[],na.rm=T)+1)-min(r[],na.rm=T))[1:((max(r[],na.rm=T)+1)-min(r[],na.rm=T))]
cols[[8]] = colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
cols[[9]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
xmin = -125; xmax = -65; ymin = 21; ymax = 53
labsX = c(expression(125*degree*W), expression(65*degree*W))
labsY = c(expression(21*degree*N), expression(53*degree*N))

dev.new(width=7.5, height=5.6); par(mfrow=c(3,3), oma=c(2,2.5,1,0.3), mar=c(0,0,0,0), mgp=c(1,0.2,0), lwd=0.2)
for (i in 1:length(rS))
	{
		plot(rS[[i]], bty="n", box=F, axes=F, legend=F, col=cols[[i]], colNA="#D0D0D0")
		axis(1, c(xmin,xmax), labels=labsX, pos=ymin(rS[[i]]), col="gray30", cex.axis=0.6, col.axis="gray30", lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.02, mgp=c(0,0.10,0))
		axis(2, c(ymin,ymax), labels=labsY, pos=xmin(rS[[i]]), col="gray30", cex.axis=0.6, col.axis="gray30", lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.02, mgp=c(0,0.27,0))
		if (i == 1)
			{
				plot(rS[[i]], legend.only=T, add=T, col=cols[[i]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.78,0.795,0.18,0.45), adj=3,
					 axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0, 
					 mgp=c(0,0.4,0), at=c(0), labels=c(0)), alpha=1, side=3)
				mtext("4651", side=3, adj=1, line=-7.1, at=-66, cex=0.45, font=1, col="gray30")
			}	else	{
				plot(rS[[i]], legend.only=T, add=T, col=cols[[i]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.78,0.795,0.18,0.54), adj=3,
					 axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)
		     }
		if (nchar(envVariableNames1[i] > 0)) mtext(envVariableNames1[i], side=1, adj=0, line=-3.1, at=-126, cex=0.55, font=1, col="gray30")
		if (nchar(envVariableNames2[i] > 0)) mtext(envVariableNames2[i], side=1, adj=0, line=-2.3, at=-126, cex=0.55, font=1, col="gray30")
		rect(xmin(rS[[i]]), ymin(rS[[i]]), xmax(rS[[i]]), ymax(rS[[i]]), lwd=0.2, border="gray30")
	}
dev.copy2pdf(file="Environmental_rasters.pdf")

# 2. Extraction of the spatio-temporal information embedded in 100 trees

allTrees = scan(file="WNV_RRW_100.trees", what="", sep="\n", quiet=T)
burnIn = 0
randomSampling = FALSE
nberOfTreesToSample = 100
coordinateAttributeName = "location"
nberOfCores = 1

treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)

# 3. Generating a dispersal history graph (by time slice or not)

	# 3.1. Building the maximum clade consensus (MCC) with TreeAnnotator
	
system(paste0("BEAST_ver_1.10.4/bin/treeannotator -burninTrees 0 -heights keep WNV_RRW_100.trees WNV_RRW_MCC.tree"))

	# 3.2. Extraction of the spatio-temporal information embedded in the MCC tree

source("Divers_R_functions/mccTreeExtraction.r")
mcc_tre = readAnnotatedNexus("WNV_RRW_MCC.tree")
mcc_tab = mccTreeExtraction(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, "WNV_RRW_MCC.csv", row.names=F, quote=F)

	# 3.3. Estimating annual kernel density polygons

nodes = c()
for (i in 1:nberOfExtractionFiles)
	{
		tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
		startingNodeID = which(!tab[,"node1"]%in%tab[,"node2"])
		startingNode = tab[startingNodeID,c("startYear","startLon","startLat")]; colnames(startingNode) = c("time","lon","lat")
		endingNodes = tab[,c("endYear","endLon","endLat")]; colnames(endingNodes) = c("time","lon","lat")
		nodes = rbind(nodes,startingNode,endingNodes)
	}

endDays = c("31","27","31","30","31","30","31","31","30","31","30","31")
months = c("01","02","03","04","05","06","07","08","09","10","11","12"); c = 0
years = c(1991:2016); timePoints = matrix(nrow=length(years)*length(months), ncol=2)
monthNames = rep(NA, dim(timePoints)[1]); colnames(timePoints) = c("startTime","endTime")
for (i in 1:length(years))
	{
		for (j in 1:length(months))
			{
				c = c+1; monthNames[c] = paste(years[i],months[j],sep="-")
				timePoints[c,"startTime"] = decimal_date(ymd(paste(years[i],months[j],"01",sep="-")))
				timePoints[c,"endTime"] = decimal_date(ymd(paste(years[i],months[j],endDays[j],sep="-")))
			}
	}
row.names(timePoints) = monthNames; timePoints = timePoints[8:308,] # from August 1991 to August 2016
write.csv(timePoints, "WNV_time_points.csv", quote=F)

timePoints = read.csv("WNV_time_points.csv", header=T); percentages = c(0.95)
for (h in 1:length(percentages))
	{
		c = 0; percentage = gsub("\\.","",as.character(percentages[h]))
		for (i in 1:dim(timePoints)[1])
			{
				selectedNodes = nodes[which((nodes[,"time"]>timePoints[i,"startTime"])&(nodes[,"time"]<timePoints[i,"endTime"])),]
				if (dim(selectedNodes)[1] > 0)
					{
						if (dim(unique(selectedNodes))[1] > 2)
							{
								layerName = paste0("WNV_contours_",percentage,"_",row.names(timePoints)[i])
								if (!file.exists(paste0("WNV_contours_",percentage,"/",layerName,".shp")))
									{
										c = c+1; H = Hpi(cbind(selectedNodes[,"lon"],selectedNodes[,"lat"])); print(row.names(timePoints)[i])
										kde = kde(cbind(selectedNodes[,"lon"],selectedNodes[,"lat"]), H=H, compute.cont=T, gridsize=c(1000,1000))
										contourLevel = contourLevels(kde, prob=0.01); polygons = list()
										contourLines = contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
										for (j in 1:length(contourLines)) polygons[[j]] = Polygon(cbind(contourLines[[j]]$x,contourLines[[j]]$y))
										ps = Polygons(polygons,1); contourPolygons = SpatialPolygons(list(ps))
										contourPolygons_df = SpatialPolygonsDataFrame(contourPolygons, data.frame(ID=1:length(contourPolygons)))
										writeOGR(contourPolygons_df, dsn=paste0("./WNV_contours_",percentage), layer=layerName, driver="ESRI Shapefile")
									}
							}	else	{
								write.csv(unique(selectedNodes), paste0("WNV_contours_",percentage,"/WNV_contours_",percentage,"_",row.names(timePoints)[i],".csv"))
							}
					}
			}
	}

	# 3.4. Plotting the global graph of dispersal history

e_Figure = extent(-145,-50,19,54)
elevation = raster("Environmental_files/WNV_rasters/Elevation_WNV_04.asc")
background = raster("Environmental_files/Natural_Earth/Gray_background.tif")
lakes = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="Natural_Earth_lakes"), e_Figure)
borders = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="International_borders"), e_Figure)
background = crop(background, e_Figure); r = background; background[background[]==106] = NA
background_cols = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
flyways = readOGR(dsn="./Environmental_files/WNV_shapefiles/", layer="Migratory_birds_flyways")
mcc = read.csv("WNV_RRW_MCC.csv", header=T); mcc = mcc[order(mcc[,"startYear"]),]
mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]; mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
endYearsM = ((mcc[,"endYear"]-min(mcc[,"startYear"]))/(max(mcc[,"endYear"])-min(mcc[,"startYear"]))*100)+1
startYearM = ((mcc[1,"startYear"]-min(mcc[,"startYear"]))/(max(mcc[,"endYear"])-min(mcc[,"startYear"]))*100)+1
cols = colorRampPalette(brewer.pal(11,'RdYlGn'))(101); cols_mcc = cols[endYearsM]; col_start = cols[startYearM]
rast = elevation; rast[!is.na(rast[])] = NA; rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
files = list.files("WNV_contours_095"); files = files[which((grepl("-",files))&(grepl(".shp",files)))]; pols = list(); cols_pol = list()
for (i in 1:length(files))
	{
		pol = shapefile(paste0("WNV_contours_095/",files[i]))
		pol@proj4string = CRS("+init=epsg:4326"); pols[[i]] = pol
		date = gsub(".shp","",unlist(strsplit(files[i],"_"))[length(unlist(strsplit(files[i],"_")))])
		date = decimal_date(ymd(paste(date,"15",sep="-")))
		yearM = ((as.numeric(date)-min(mcc[,"startYear"]))/(max(mcc[,"endYear"])-min(mcc[,"startYear"]))*100)+1
		cols_pol[[i]] = cols[yearM]
	}
for (i in 1:length(cols_pol)) cols_pol[[i]] = paste0(cols_pol[[i]],"70")
xmin = -123; xmax = -67; ymin = 21; ymax = 52
labsX = c(expression(123*degree*W), expression(67*degree*W)); labsY = c(expression(21*degree*N), expression(52*degree*N))

dev.new(width=10, height=5); par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1,2.8,0.9,0.5), mgp=c(1,0.2,0), lwd=0.3)
plot(background, main="", cex.main=1, bty="n", box=F, axes=F, legend=F, col=background_cols, colNA="#D0D0D0")
plot(borders, add=T, lwd=1, col="white", lty=1); plot(lakes, add=T, lwd=0.7, col="#D0D0D0", border=NA)
plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.911,0.917,0.112,0.892), adj=3,
	 axis.args=list(at=c(1999:2016), cex.axis=0.6, lwd=0, col="gray30", lwd.tick=0.5, col.tick="gray30", tck=-0.5, col.axis="gray30", line=0, mgp=c(0,0.3,0)), alpha=1, side=3)
for (i in dim(mcc)[1]:1)
	{
		curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
	}
for (i in dim(mcc)[1]:1)
	{
		points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=cols_mcc[i], cex=0.7)
		points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", cex=0.7, lwd=0.1)
		if (i == 1)
			{
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=col_start, cex=0.7)
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", cex=0.7, lwd=0.1)
			}
	}
axis(1, c(xmin,xmax), labels=labsX, pos=ymin(r), cex.axis=0.5, col.axis="gray30", lwd=0, lwd.tick=0.3, tck=0.008, col.tick="gray30", mgp=c(0,-1.1,0), col="gray30", asp=2)
axis(2, c(ymin,ymax), labels=labsY, pos=xmax(r), cex.axis=0.5, col.axis="gray30", lwd=0, lwd.tick=0.5, tck=-0.008, col.tick="gray30", mgp=c(0,0.15,0), col="gray30", asp=2)
rect(xmin(background), ymin(background), xmax(background), ymax(background), xpd=T, lwd=0.5, border="gray30")

	# 3.5. Plotting dispersal history graph snapshots

years = c(1999:2016); cols = colorRampPalette(brewer.pal(11,'RdYlGn'))(101)
rast = raster("Environmental_files/WNV_rasters/Elevation_WNV_04.asc"); rast[!is.na(rast)] = 0
lakes = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="Natural_Earth_lakes"), extent(rast))
borders = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="Coastline_borders"), extent(rast))
mcc = read.csv("WNV_RRW_MCC.csv", header=T); mcc = mcc[order(mcc[,"startYear"]),]
mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]; mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
pols_list = list(); cols_pol_list = list(); mccs = list(); cols_mcc_1 = list(); cols_mcc_2 = list()
for (i in 1:(length(years)-1))
	{
		files = list.files("WNV_contours_095"); pols = list(); cols_pol = list()
		files = files[which((grepl(paste0(years[i],"-"),files))&(grepl(".shp",files)))]
		for (j in 1:length(files))
			{
				pol = shapefile(paste0("WNV_contours_095/",files[j]))
				pol@proj4string = CRS("+init=epsg:4326"); pols[[j]] = pol
				date = gsub(".shp","",unlist(strsplit(files[j],"_"))[length(unlist(strsplit(files[j],"_")))])
				date = decimal_date(ymd(paste(date,"15",sep="-")))
				yearM = ((as.numeric(date)-min(mcc[,"startYear"]))/(max(mcc[,"endYear"])-min(mcc[,"startYear"]))*100)+1
				cols_pol[[j]] = cols[yearM]
			}
		for (j in 1:length(cols_pol)) cols_pol[[j]] = paste0(cols_pol[[j]],"70")
		pols_list[[i]] = pols; cols_pol_list[[i]] = cols_pol
		mccs[[i]] = mcc[which(mcc[,"endYear"]<(years[i]+1)),]	
		startYearsM = ((mccs[[i]][,"startYear"]-min(mcc[,"startYear"]))/(max(mcc[,"endYear"])-min(mcc[,"startYear"]))*100)+1
		endYearsM = ((mccs[[i]][,"endYear"]-min(mcc[,"startYear"]))/(max(mcc[,"endYear"])-min(mcc[,"startYear"]))*100)+1
		startYearM = ((mccs[[i]][1,"startYear"]-min(mcc[,"startYear"]))/(max(mcc[,"endYear"])-min(mcc[,"startYear"]))*100)+1
		cols_mcc_1[[i]] = cols[startYearsM]; cols_mcc_2[[i]] = cols[endYearsM]; col_start = cols[startYearM]
	}
xmin = -125; xmax = -65; ymin = 21; ymax = 53
labsX = c(expression(125*degree*W), expression(65*degree*W))
labsY = c(expression(21*degree*N), expression(53*degree*N))

plotBranches = TRUE; cumulative = FALSE
for (h in 1:2)
	{
		pdf(paste0("WNV_figures_&_SI/SI_files/Dispersal_years_",h,"_NEW.pdf"), width=7.5, height=5.6); par(mfrow=c(3,3), oma=c(2,2.5,1,0.3), mar=c(0,0,0,0), mgp=c(1,0.2,0), lwd=0.2)
		# dev.new(width=7.5, height=5.6); par(mfrow=c(3,3), oma=c(2,2.5,1,0.3), mar=c(0,0,0,0), mgp=c(1,0.2,0), lwd=0.2)
		for (i in (((h-1)*9)+(c(1:9))))
			{
				if (file.exists(paste0("WNV_contours_095/WNV_contours_095_",years[i],".shp")))
					{
						plot(rast, bty="n", box=F, axes=F, legend=F, col="white", colNA="#D0D0D0")
						if (cumulative == FALSE)
							{
								for (k in 1:length(pols_list[[i]])) plot(pols_list[[i]][[k]], axes=F, col=cols_pol_list[[i]][[k]], add=T, border=NA)
							}	else	{
								for (j in 1:i)
									{
										for (k in 1:length(pols_list[[j]])) plot(pols_list[[j]][[k]], axes=F, col=cols_pol_list[[j]][[k]], add=T, border=NA)
									}
							}
						plot(rast, bty="n", box=F, axes=F, legend=F, col=rgb(0,0,0,0), colNA="#D0D0D0", add=T)
						plot(lakes, lwd=0.7, col="#D0D0D0", border=NA, add=T)
						if (plotBranches == TRUE)
							{
								if ((cumulative == FALSE)&(i > 1))
									{
										col_mcc_1 = cols_mcc_1[[i]][which(mccs[[i]][,"endYear"]>years[i-1])]
										col_mcc_2 = cols_mcc_2[[i]][which(mccs[[i]][,"endYear"]>years[i-1])]
										mcc = mccs[[i]][which(mccs[[i]][,"endYear"]>=years[i]),]
									}	else	{
										col_mcc_1 = cols_mcc_1[[i]]; col_mcc_2 = cols_mcc_2[[i]]; mcc = mccs[[i]]
									}
								for (j in 1:dim(mcc)[1])
									{
										curvedarrow(cbind(mcc[j,"startLon"],mcc[j,"startLat"]), cbind(mcc[j,"endLon"],mcc[j,"endLat"]), arr.length=0, arr.width=0,
							  						lwd=0.01, lty=1, lcol="gray30", arr.pos=F, curve=0.1, dr=NA, endhead=F, arr.lwd=0.01, arr.col=rgb(0,0,0,0))
							  		}
							  	for (j in 1:dim(mcc)[1])
							  		{
							  			points(mcc[j,"startLon"], mcc[j,"startLat"], pch=16, col=col_mcc_1[j], cex=0.5)
							  			points(mcc[j,"startLon"], mcc[j,"startLat"], pch=1, col="gray30", cex=0.5, lwd=0.01)
									}
							  	for (j in 1:dim(mcc)[1])
							  		{
							  			points(mcc[j,"endLon"], mcc[j,"endLat"], pch=16, col=col_mcc_2[j], cex=0.5)
							  			points(mcc[j,"endLon"], mcc[j,"endLat"], pch=1, col="gray30", cex=0.5, lwd=0.01)
									}
							}
						axis(1, c(xmin,xmax), labels=labsX, pos=ymin(rast), cex.axis=0.6, col.axis="gray30", lwd=0, lwd.tick=0.2, tck=-0.02, col.tick="gray30", mgp=c(0,0.10,0), col="gray30")
						axis(2, c(ymin,ymax), labels=labsY, pos=xmin(rast), cex.axis=0.6, col.axis="gray30", lwd=0, lwd.tick=0.2, tck=-0.02, col.tick="gray30", mgp=c(0,0.27,0), col="gray30")
						mtext(years[i], side=1, adj=0, line=-2.3, at=-126, cex=0.6, font=1, col="gray30") # col=cols_pol_list[[i]][[length(cols_pol_list[[i]])]])
						rect(xmin(rast), ymin(rast), xmax(rast), ymax(rast), lwd=0.2, border="gray30")
					}
			}
		dev.off()
	}

# 4. Estimating and plotting different dispersal statistics

timSlices = 200; onlyTipBranches = FALSE; showingPlots = FALSE; outputName = "WNV"; nberOfCores = 1; slidingWindow = 1
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)

stats = read.table("Dispersal_statistics/WNV_estimated_dispersal_statistics.txt", header=T)
v = stats[,"mean_branch_dispersal_velocity"]; D = stats[,"original_diffusion_coefficient"]; D = D/365
cat("Mean branch dispersal velocity:"); cat("\n"); cat("\t")
cat(paste0(round(mean(v),1)," km/year, 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]"))
cat("Original diffusion coefficient:"); cat("\n"); cat("\t")
cat(paste0(round(mean(D),1)," km2/day, 95% HPD = [",round(quantile(D,0.025),1),",",round(quantile(D,0.975),1),"]"))

dev.new(width=4, height=1.5); par(mgp=c(0,0,0), oma=c(0.8,0.8,0,0), mar=c(1.5,1.5,1,1))
col1 = rgb(100, 100, 100, 255, maxColorValue=255); col2 = rgb(100, 100, 100, 100, maxColorValue=255)
tab1 = read.table("Dispersal_statistics/WNV_median_spatial_wavefront_distance.txt", header=T)
tab2 = read.table("Dispersal_statistics/WNV_95%HPD_spatial_wavefront_distance.txt", header=T)
plot(tab1[,1], tab1[,2], type="l", axes=F, ann=F, ylim=c(0,4300), xlim=c(1998,2010), col=NA)
slicedTimes = tab1[,1]; waveFrontDistances1MedianValue = tab1[,2]; lower_l_1 = tab2[,2]; upper_l_1 = tab2[,3]
xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=col2, border=0)
lines(slicedTimes, waveFrontDistances1MedianValue, lwd=1, col=col1)
axis(side=1, pos=0, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.050, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, pos=1998, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.045, col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab="time (year)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
title(ylab="distance (km)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
dev.copy2pdf(file="Figure1_spatial_wavefronts.pdf")

dev.new(width=4, height=4); par(mgp=c(0,0,0), oma=c(0.8,0.8,0,0), mar=c(1.5,1.5,1,1))
col1 = rgb(100, 100, 100, 255, maxColorValue=255); col2 = rgb(100, 100, 100, 100, maxColorValue=255)
tab1 = read.table("Dispersal_statistics/WNV_median_mean_branch_dispersal_velocity.txt", header=T)
tab2 = read.table("Dispersal_statistics/WNV_95%HPD_mean_branch_dispersal_velocity.txt", header=T)
plot(tab1[,1], tab1[,2], type="l", axes=F, ann=F, ylim=c(0,10000), xlim=c(1998,2010), col=NA)
slicedTimes = tab1[,1]; branchDispersalVelocityMeanValue = tab1[,2]; lower_l_1 = tab2[,2]; upper_l_1 = tab2[,3]
xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=col2, border=0)
lines(slicedTimes, branchDispersalVelocityMeanValue, lwd=1, col=col1)
axis(side=1, pos=0, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.050, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, pos=1998, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.045, col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab="time (year)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
title(ylab="mean branch dispersal velocity", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
ddev.copy2pdf(file="Figure1_branch_dispersal_velocity.pdf")

# 5. Performing RRW simulations along trees (null model)

simulationsDirectory = localTreesDirectory; n1 = 100; n2 = 100; showingPlots = FALSE; newPlot = FALSE; model = "gamma"
trees = read.annotated.nexus("WNV_RRW_100.trees"); log = read.table("WNV_RRW_100_b.log", header=T)
for (i in 1:nberOfExtractionFiles)
	{
		tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
		rast1 = raster("Environmental_files/WNV_rasters/Elevation_WNV_04.asc"); rast1[!is.na(rast1[])] = 0
		 if (i == 1)
		 	{
		 		points1 = tab[,c("startLon","startLat")]; points2 = tab[,c("endLon","endLat")]
		 	}	else		{
		 		points1 = rbind(points1, tab[,c("startLon","startLat")])
		 		points2 = rbind(points2, tab[,c("endLon","endLat")])
		 	}
	}
colnames(points1) = c("lon","lat"); colnames(points2) = c("lon","lat")
points = rbind(points1,points2); hull = chull(points); hull = c(hull,hull[1])
p = Polygon(points[hull,]); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
pointsRaster = rasterize(points, crop(rast1, sps, snap="out"))
pointsRaster[!is.na(pointsRaster[])] = 0; # plot(mask(simRasters[[h]],sps))
hullRaster = crop(rast1, sps, snap="out"); bufferRaster = hullRaster
rast2 = mask(hullRaster, sps, snap="out")
rast2[!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
for (i in 1:nberOfExtractionFiles)
	{
		print(paste0("Simulating tree #",i))
		tree = trees[[i]]; envVariables = list(rast2)
		rates = c(); geoDists = matrix(nrow=dim(tab)[1], ncol=1)
		for (j in 1:length(tree$annotations))
			{
				rates = c(rates, tree$annotations[[j]]$location.rate)
			}
		for (j in 1:dim(tab)[1])
			{
				x1 = cbind(tab[j,"startLon"], tab[j,"startLat"])
				x2 = cbind(tab[j,"endLon"], tab[j,"endLat"])
				geoDists[j,1] = rdist.earth(x1, x2, miles=F, R=NULL)
			}
		ancestID = which(!tab[,"node1"]%in%tab[,"node2"])[1]
		ancestPosition = c(tab[ancestID,"startLon"], tab[ancestID,"startLat"])
		col11 = log[i,"treeLengthPrecision1"]
		col12 = log[i,"treeLengthPrecision3"]
		col22 = log[i,"treeLengthPrecision2"]
		my_prec = c(col11, col12, col12, col22)
		if (model == "cauchy") reciprocalRates = TRUE
		if (model == "gamma") reciprocalRates = TRUE
		if (model == "logN") reciprocalRates = FALSE
		my_var = solve(matrix(my_prec,nrow=2))		
		sigma1 = sqrt(my_var[1,1]); sigma2 = sqrt(my_var[2,2])
		sigmas = c(sigma1, sigma2); # source("simulatorRRW1_mod.r")
		cor = my_var[1,2]/(sqrt(my_var[1,1])*sqrt(my_var[2,2]))
		source("Seraphim_functions/simulatorRRW1.r"); # showingPlots = TRUE; newPlot = TRUE
		output = simulatorRRW1(tree, rates, sigmas, cor, envVariables, mostRecentSamplingDatum,
							   ancestPosition, reciprocalRates, n1, n2, showingPlots, newPlot)
		file = as.character(paste(simulationsDirectory,"/TreeSimulations_",i,".csv",sep=""))
		write.csv(output, file, row.names=F, quote=F)
	}

# 6. Testing the impact of environmental factors on branch dispersal velocity

envVariables = list(); resistances = list(); avgResistances = list(); fourCells = FALSE
nberOfRandomisations = 0; randomProcedure = 3; showingPlots = FALSE; nberOfCores = 10; OS = "Unix"; simulations = FALSE
envVariableNames = c("Elevation","Land_cover_forests","Land_cover_shrublands","Land_cover_savannas","Land_cover_grasslands",
					 "Land_cover_croplands","Land_cover_water","Annual_mean_temperature","Annual_precipitation"); c = 0
for (k in c(10,100,1000))
	{
		for (i in 1:length(envVariableNames))
			{
				c = c+1
				rast = raster(paste("Environmental_files/WNV_rasters/",envVariableNames[i],"_WNV_16.asc",sep=""))
				rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
				names(rast) = paste(envVariableNames[i],"_k",k,sep="")
				envVariables[[c]] = rast; names(envVariables[[c]]) = paste(envVariableNames[i],"_k",k,sep="")
				resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
			}
		for (i in 1:length(envVariableNames))
			{
				c = c+1
				rast = raster(paste("Environmental_files/WNV_rasters/",envVariableNames[i],"_WNV_16.asc",sep=""))
				rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
				names(rast) = paste(envVariableNames[i],"_k",k,sep="")
				envVariables[[c]] = rast; names(envVariables[[c]]) = paste(envVariableNames[i],"_k",k,sep="")
				resistances[[c]] = FALSE; avgResistances[[c]] = FALSE
			}
	}

pathModel = 2; outputName = "WNV_least-cost_extractions"
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations)
pathModel = 3; outputName = "WNV_circuitscape_extractions"
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations)
simulations = TRUE
pathModel = 2; outputName = "WNV_least-cost_simulations"
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations)
pathModel = 3; outputName = "WNV_circuitscape_simulations"
spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations)

extractions = list(); simulations = list(); pathModels = c("Least-cost path model","Circuitscape path model")
extractions[[1]] = read.table("Seraphim_analyses/WNV_least-cost_extractions_16_LR_results.txt", header=T)
extractions[[2]] = read.table("Seraphim_analyses/WNV_circuitscape_extractions_16_LR_results.txt", header=T)
simulations[[1]] = read.table("Seraphim_analyses/WNV_least-cost_simulations_16_LR_results.txt", header=T)
simulations[[2]] = read.table("Seraphim_analyses/WNV_circuitscape_simulations_16_LR_results.txt", header=T)
envVariableNames = c("Elevation","Land_cover_forests","Land_cover_shrublands","Land_cover_savannas","Land_cover_grasslands",
					 "Land_cover_croplands","Land_cover_water","Annual_mean_temperature","Annual_precipitation"); c = 0
allResults = matrix(nrow=length(envVariableNames)*2*2*3, ncol=7); kS = c(10,100,1000); CR = c("C","R"); L = 0
colnames(allResults) = c("Path model","Environmental factor","k","Regression coefficient","Q statistic","p(Q) > 0","BF")
for (i in 1:length(pathModels))
	{
		for (j in 1:length(envVariableNames))
			{
				for (k in 1:length(CR))
					{
						for (l in 1:length(kS))
							{
								L = L+1
								allResults[L,1] = pathModels[i]
								allResults[L,2] = paste0(envVariableNames[j]," (",CR[k],")")
								allResults[L,3] = kS[l]
								index1 = which(grepl("LR_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))&grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
								index2 = which(grepl("delta_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))&grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
								index3 = which(grepl("delta_R2",colnames(simulations[[i]]))&grepl(envVariableNames[j],colnames(simulations[[i]]))&grepl(paste0("k",kS[l],"_",CR[k]),colnames(simulations[[i]])))
								R2 = extractions[[i]][,index1]; Qe = extractions[[i]][,index2]; Qs = simulations[[i]][,index3]; c = 0
								for (m in 1:length(Qe))
									{
										if (Qs[m] < Qe[m]) c = c+1
									}
								p = c/length(Qe); BF = p/(1-p)
								allResults[L,4] = paste0(round(median(R2),3)," [",round(quantile(R2,0.025),3)," - ",round(quantile(R2,0.975),3),"]")
								allResults[L,5] = paste0(round(median(Qe),3)," [",round(quantile(Qe,0.025),3)," - ",round(quantile(Qe,0.975),3),"]")
								allResults[L,6] = sum(Qe>0)/nberOfExtractionFiles
								if (as.numeric(allResults[L,6]) >= 0.9)
									{
										allResults[L,7] = round(BF,1)
									}	else	{
										allResults[L,7] = "-"
									}
							}
					}
			}
	}
write.csv(allResults, "Seraphim_analyses.csv", row.names=F, quote=F)

dev.new(width=6, height=2.5); par(mfrow=c(1,1), mgp=c(1,0.35,0), oma=c(0.5,1,1,1), mar=c(2.5,2,0.5,0), lwd=0.2)
envVariableNames = c("Annual_mean_temperature_C"); envVariableTitle1 = c("Impact of annual mean temperature")
envVariableTitle2 = c("on branch dispersal velocity"); envVariableTitle3 = c("(tested as conductance factor)")
Qe = list(); Qs = list(); cols1 = list(); cols2 = list(); kS = c(100,1000); ltys = c(1,3); xLim=c(-0.17,0.1)
cols1[[1]] = rgb(204,0,0,255,maxColorValue=255); cols2[[1]] = rgb(204,0,0,100,maxColorValue=255) # red
cols1[[2]] = rgb(120,120,120,255,maxColorValue=255); cols2[[2]] = rgb(120,120,120,100,maxColorValue=255)
for (i in 1:length(envVariableNames))
	{
		for (k in 1:length(kS))
			{
				envVariableName = gsub("_C",paste0("_k",kS[k],"_C"),envVariableNames[i])
				Qe[[i]] = extractions_CS[,paste0("Univariate_LR_delta_R2_",envVariableName)]
				Qs[[i]] = simulations_CS[,paste0("Univariate_LR_delta_R2_",envVariableName)]
				if (kS[k] == 100)
					{
						plot(density(Qe[[i]]), lwd=0.7, col=cols1[[1]], lty=ltys[k], xlim=xLim, axes=F, ann=F)
					}	else	{
						lines(density(Qe[[i]]), lwd=0.7, col=cols1[[1]], lty=ltys[k])
					}
				lines(density(Qs[[i]]), lwd=0.7, col=cols1[[2]], lty=ltys[k])
			}
		for (k in 1:1)
			{
				envVariableName = gsub("_C",paste0("_k",kS[k],"_C"),envVariableNames[i])
				Qe[[i]] = extractions_CS[,paste0("Univariate_LR_delta_R2_",envVariableName)]
				Qs[[i]] = simulations_CS[,paste0("Univariate_LR_delta_R2_",envVariableName)]
				polygon(density(Qe[[i]]), col=cols2[[1]], border=NA)
				polygon(density(Qs[[i]]), col=cols2[[2]], border=NA)
			}
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.02,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		title(xlab=expression(italic(Q) == {R^{2}}[env] - {R^{2}}[null]), cex.lab=0.7, mgp=c(1,0,0), col.lab="gray30")
		title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
		legend(x=-0.162, y=35, legend=c("inferred trees","simulated trees",expression(italic(k) == 1000),expression(italic(k) == 100)), 
			   lwd=0.7, cex=0.7, lty=c(1,1,rev(ltys)), col=c(unlist(cols1),"gray30","gray30"), text.col=c(unlist(cols1),"gray30","gray30"),
			   border=NA, x.intersp=0.5, bty="n")
	}
dev.copy2pdf(file="Impact_dispersal_velocity.pdf")

# 7. Testing the impact of migratory flyways on branch dispersal frequency

	# 7.1. Testing the impact of the four US administrative flyways (NAMFs)

flyways = readOGR(dsn="Environmental_files/WNV_shapefiles/", layer="Migratory_birds_flyways")
rast1 = raster("Environmental_files/WNV_rasters/Elevation_WNV_08.asc"); rast1[!is.na(rast1[])] = 0
rast2 = mask(crop(rast1, flyways), flyways); showingPlots = FALSE
ratiosOfChangingFlywayEvents = matrix(nrow=nberOfExtractionFiles, ncol=2)
colnames(ratiosOfChangingFlywayEvents) = c("extractions","simulations")
for (i in 1:nberOfExtractionFiles)
	{
		tabs = list(); print(i)
		tabs[[1]] = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
		tabs[[2]] = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",i,".csv"), header=T)
		for (j in 1:length(tabs))
			{
				if (showingPlots == TRUE)
					{
						dev.new(width=8, height=5.8); par(oma=c(0,0,0,0), mar=c(0,0,0,0))
						plotRaster(rast1, new=F); plot(flyways, add=T, lwd=0.1)
					}
				indices = which((!is.na(extract(rast2,tabs[[j]][,c("startLon","startLat")])))&(!is.na(extract(rast2,tabs[[j]][,c("endLon","endLat")]))))
				tab = tabs[[j]][indices,]; n = 0
				for (k in 1:dim(tab)[1])
					{
						pts = cbind(c(tab[k,"startLon"],tab[k,"endLon"]),c(tab[k,"startLat"],tab[k,"endLat"]))
						line = SpatialLines(list(Lines(Line(pts), ID="a"))); crs(line) = crs(flyways)
						intersections = gIntersects(line, flyways, byid=T)
						if (sum(intersections) > 1)
							{
								n = n+1
								if (showingPlots == TRUE) plot(line, add=T, col="red", lwd=0.3)
							}	else	{
								if (showingPlots == TRUE) plot(line, add=T, col="green3", lwd=0.3)
							}
					}
				ratiosOfChangingFlywayEvents[i,j] = n/dim(tab)[1]
			}
	}
c = 0
for (i in 1:dim(ratiosOfChangingFlywayEvents)[1])
	{
		if (ratiosOfChangingFlywayEvents[i,"extractions"] < ratiosOfChangingFlywayEvents[i,"simulations"]) c = c+1
	}
p = c/dim(ratiosOfChangingFlywayEvents)[1]; BF = p/(1-p); print(BF)

	# 7.2. Testing the impact of the three flyways estimated by La Sorte et al. (2014) for terrestrial birds

rast = raster("Environmental_files/WNV_rasters/Western_flyway_normalised.asc"); rast[!is.na(rast[])] = 0
differenceBetweenFlywayValues = matrix(nrow=nberOfExtractionFiles, ncol=3); flyways = list()
colnames(differenceBetweenFlywayValues) = c("extractions","simulations","randomisations")
flyways[[1]] = raster("Environmental_files/WNV_rasters/Eastern_flyway_normalised.asc")
flyways[[2]] = raster("Environmental_files/WNV_rasters/Central_flyway_normalised.asc")
flyways[[3]] = raster("Environmental_files/WNV_rasters/Western_flyway_normalised.asc")
for (i in 1:nberOfExtractionFiles)
	{
		tabs = list(); print(i)
		tabs[[1]] = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
		tabs[[2]] = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",i,".csv"), header=T)
		for (j in 1:length(tabs))
			{
				indices = which((!is.na(extract(rast,tabs[[j]][,c("startLon","startLat")])))
								&(!is.na(extract(rast,tabs[[j]][,c("endLon","endLat")]))))
				tab = tabs[[j]][indices,]; diffs = c()
				startingValues = matrix(nrow=dim(tab)[1], ncol=length(flyways))
				endingValues = matrix(nrow=dim(tab)[1], ncol=length(flyways))
				for (k in 1:length(flyways))
					{
						startingValues[,k] = extract(flyways[[k]], tab[,c("startLon","startLat")])
						endingValues[,k] = extract(flyways[[k]], tab[,c("endLon","endLat")])
					}	
				for (k in 1:dim(tab)[1])
					{
						startingValue = max(startingValues[k,])
						index = which(startingValues[k,]==startingValue)
						endingValue = endingValues[k,index]
						diffs = c(diffs, endingValue-startingValue)
					}
				differenceBetweenFlywayValues[i,j] = mean(diffs)
			}
	}
c = 0
for (i in 1:dim(differenceBetweenFlywayValues)[1])
	{
		if (differenceBetweenFlywayValues[i,"extractions"] > differenceBetweenFlywayValues[i,"simulations"]) c = c+1
	}
p = c/dim(differenceBetweenFlywayValues)[1]; BF = p/(1-p); print(BF)

