library(adephylo)
library(diagram)
library(maptools)
library(lubridate)
library(rgdal)
library(rgeos)
library(seraphim)

# 1. Extraction of the spatio-temporal information embedded in 100 trees
# 2. Generating and saving annual and monthly progression polygons
# 3. RRW simulations along tree branches (generating the null dispersal model)
# 4. Subdividing lineages by genotypes and using a time cut-off (end of invasion phase)
# 5. Generating a dispersal history graph (by time slice or not)
# 6. Estimating and plotting some dispersal statistics
# 7. Testing the impact of environmental factors on lineage dispersal locations
# 8. Testing the impact of environmental factors on lineage dispersal tendency
# 9. Testing the impact of environmental factors on lineage dispersal velocity
# 10. Testing the impact of migratory flyways on lineage dispersal frequency

e_WNV = extent(-128, -63, 19, 55) # extent of study area
localTreesDirectory = "Tree_extractions/WNV_gamma_all"
nberOfExtractionFiles = 100; mostRecentSamplingDatum = 2016.6475

# 1. Extraction of the spatio-temporal information embedded in 100 trees

allTrees = scan(file="WNV_RRW_100.trees", what="", sep="\n", quiet=T)
burnIn = 0
randomSampling = FALSE
nberOfTreesToSample = 100
coordinateAttributeName = "location"
nberOfCores = 1

treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
   
# 2. Generating and saving annual and monthly progression polygons

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

	# 2.1. Minimum convex hulls as monthly progression polygons

timePoints = read.csv("WNV_time_points.csv", header=T)
for (i in 1:dim(timePoints)[1])
	{
		selectedNodes = nodes[which((nodes[,"time"]<timePoints[i,"endTime"])),]
		if (dim(selectedNodes)[1] > 0)
			{
				if (dim(unique(selectedNodes))[1] > 2)
					{
						layerName = paste0("WNV_convex_hull_",row.names(timePoints)[i])
						if (!file.exists(paste0("WNV_convex_hulls/",layerName,".shp")))
							{
								hull = chull(selectedNodes[,2:3]); hull = c(hull,hull[1])
								p = Polygon(selectedNodes[hull,2:3]); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								contourPolygons_df = SpatialPolygonsDataFrame(sps, data.frame(ID=1:length(sps)))
								writeOGR(contourPolygons_df, dsn=paste0("./WNV_convex_hulls"), layer=layerName, driver="ESRI Shapefile")
							}
					}	else	{
						write.csv(unique(selectedNodes), paste0("WNV_convex_hulls/WNV_convex_hull_",row.names(timePoints)[i],".csv"))
					}
			}
	}

gadm_0_USA = gSimplify(crop(getData("GADM", country="USA", level=0), e_WNV), 0.1)
dev.new(); plot(gadm_0_USA, lwd=0.5, border="blue")
for (i in 1:dim(timePoints)[1])
	{
		layerName = paste0("WNV_convex_hull_",row.names(timePoints)[i])
		if (file.exists(paste0("WNV_convex_hulls/",layerName,".shp")))
			{
				sps = readOGR(dsn=paste0("./WNV_convex_hulls"), layer=layerName); plot(sps, lwd=0.05, add=T)
			}
	}

	# 2.2. Kernel density areas as annual progression polygons

years = c(1991:2016); percentages = c(0.95)
for (h in 1:length(percentages))
	{
		c = 0; percentage = gsub("\\.","",as.character(percentages[h]))
		for (i in 1:length(years))
			{
				startTime = years[i]; endTime = years[i]+1
				selectedNodes = nodes[which((nodes[,"time"]>=startTime)&(nodes[,"time"]<endTime)),]
				if (dim(selectedNodes)[1] > 0)
					{
						if (dim(unique(selectedNodes))[1] > 2)
							{
								layerName = paste0("WNV_contours_",percentage,"_",years[i])
								if (!file.exists(paste0("WNV_contours_",percentage,"/",layerName,".shp")))
									{
										c = c+1; H = Hpi(cbind(selectedNodes[,"lon"],selectedNodes[,"lat"])); print(years[i])
										kde = kde(cbind(selectedNodes[,"lon"],selectedNodes[,"lat"]), H=H, compute.cont=T, gridsize=c(1000,1000))
										contourLevel = contourLevels(kde, prob=0.01); polygons = list()
										contourLines = contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
										for (j in 1:length(contourLines)) polygons[[j]] = Polygon(cbind(contourLines[[j]]$x,contourLines[[j]]$y))
										ps = Polygons(polygons,1); contourPolygons = SpatialPolygons(list(ps)); # plot(contourPolygons, add=T)
										contourPolygons_df = SpatialPolygonsDataFrame(contourPolygons, data.frame(ID=1:length(contourPolygons)))
										writeOGR(contourPolygons_df, dsn=paste0("./WNV_contours_",percentage), layer=layerName, driver="ESRI Shapefile")
									}
							}	else		{
								write.csv(unique(selectedNodes), paste0("WNV_contours_",percentage,"/WNV_contours_",percentage,"_",years[i],".csv"), quote=F, row.names=F)
							}
					}
			}
	}

	# 2.3. Kernel density areas as monthly progression polygons

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

# 3. RRW simulations along tree branches (generating the null dispersal model)

n1 = 100; n2 = 100
simulationsDirectory = localTreesDirectory; showingPlots = FALSE; newPlot = FALSE; model = "gamma"
trees = readAnnotatedNexus("WNV_RRW_100.trees"); log = read.table("WNV_RRW_100_b.log", header=T)
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
pointsRaster[!is.na(pointsRaster[])] = 0
hullRaster = crop(rast1, sps, snap="out"); bufferRaster = hullRaster
rast2 = mask(hullRaster, sps, snap="out")
rast2[!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
registerDoMC(cores=10); buffer = list()
for (i in 1:nberOfExtractionFiles)
	{
		print(paste0("Simulating tree #",i))
		tree = trees[[i]]
		envVariables = list(rast2)
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

# 4. Subdividing lineages by genotypes and using a time cut-off (end of invasion phase)

	# 4.1. Subdividing lineages by genotypes

trees = read.nexus("WNV_RRW_100.trees")
nextStrain_data = read.csv("NextStrain_data.csv", header=T)
tipLabels = trees$tip.label; genotypes = rep(NA, length(tipLabels$tip.label))
for (i in 1:length(tipLabels$tip.label))
	{
		tipLabel = unlist(strsplit(tipLabels$tip.label[i],"_"))[1]
		index = which(nextStrain_data[,"sequenceID"]==tipLabel)
		if (length(index) == 1)
			{
				genotypes[i] = as.character(nextStrain_data[index,"strain"])
			}
	}
genotypes = cbind(tipLabels$tip.label, genotypes); colnames(genotypes) = c("","genotype")
write.table(genotypes, "WNV_genotypes.txt", row.names=F, sep="	", quote=F)

genotypes = cbind(genotypes, rep(NA,dim(genotypes)[1]), rep(NA,dim(genotypes)[1]), rep(NA,dim(genotypes)[1]))
colnames(genotypes) = c("tipLabel","genotype","samplingDate","longitude","latitude")
for (i in 1:dim(genotypes)[1])
	{
		tipLabel = unlist(strsplit(genotypes[i,"tipLabel"],"_"))
		genotypes[i,"samplingDate"] = round(as.numeric(tipLabel[2]),4)
		genotypes[i,"longitude"] = round(as.numeric(tipLabel[6]),4)
		genotypes[i,"latitude"] = round(as.numeric(tipLabel[7]),4)
	}
for (i in 0:nberOfExtractionFiles)
	{
		if (i == 0)
			{
				obs = read.csv(paste0("WNV_RRW_MCC.csv"), header=T)
			}	else		{
				obs = read.csv(paste0("Tree_extractions/WNV_gamma_all/TreeExtractions_",i,".csv"), header=T)
				sim = read.csv(paste0("Tree_extractions/WNV_gamma_all/TreeSimulations_",i,".csv"), header=T)
				ran = read.csv(paste0("Tree_extractions/WNV_gamma_all/TreeRandomisation_",i,".csv"), header=T)
			}
		buffer_obs = matrix(nrow=dim(obs), ncol=2); colnames(buffer_obs) = c("startGenotype","endGenotype")
		if (!"endGenotype"%in%colnames(obs)) obs = cbind(obs, buffer_obs)
		if (i > 0)
			{
				buffer_sim = matrix(nrow=dim(sim), ncol=2); colnames(buffer_sim) = c("startGenotype","endGenotype")
				buffer_ran = matrix(nrow=dim(ran), ncol=2); colnames(buffer_ran) = c("startGenotype","endGenotype")
				if (!"endGenotype"%in%colnames(sim)) sim = cbind(sim, buffer_sim)
				if (!"endGenotype"%in%colnames(ran)) ran = cbind(ran, buffer_ran)
			}
		for (j in 1:dim(obs)[1])
			{
				if (!obs[j,"node2"]%in%obs[,"node1"])
					{
						index = which((as.numeric(genotypes[,"samplingDate"])==round(obs[j,"endYear"],4))
									  &(as.numeric(genotypes[,"longitude"])==round(obs[j,"endLon"],4))
									  &(as.numeric(genotypes[,"latitude"])==round(obs[j,"endLat"],4)))
						if (length(index) == 1)
							{
								obs[j,"endGenotype"] = genotypes[index,"genotype"]
							}
					}
			}
		for (j in 1:dim(obs)[1])
			{
				if ((!obs[j,"node2"]%in%obs[,"node1"])&(is.na(obs[j,"endGenotype"])))
					{
						node1 = obs[j,"node1"]; node2 = obs[j,"node2"]
						index = which((obs[,"node1"]==node1)&(obs[,"node2"]!=node2))
						while (is.na(obs[index,"endGenotype"]))
							{
								temp1 = node1; node1 = obs[which(obs[,"node2"]==node1),"node1"]
								index = which((obs[,"node1"]==node1)&(obs[,"node2"]!=temp1))
								while (obs[index,"node2"]%in%obs[,"node1"])
									{
										index = which(obs[,"node1"]==obs[index,"node2"])[1]
									}
							}
						obs[j,"endGenotype"] = obs[index,"endGenotype"]
					}
			}
		for (j in 1:dim(obs)[1])
			{
				if ((!obs[j,"node2"]%in%obs[,"node1"])&(is.na(obs[j,"endGenotype"]))) print(j)
			}
		tipBranches = which(!obs[,"node2"]%in%obs[,"node1"]); tipBranch_nodes1 = obs[tipBranches,"node1"]
		threeNodeDates = paste0(round(obs[tipBranches,"startYear"],4),"_",round(obs[tipBranches,"endYear"],4))
		for (j in 1:length(threeNodeDates))
			{
				startYear = obs[which(obs[,"node2"]==tipBranch_nodes1[j]),"startYear"]
				threeNodeDates[j] = paste0(round(startYear,4),"_",threeNodeDates[j])
			}
		tabs = list(); tabs[[1]] = obs
		if (i > 0)
			{
				tabs[[2]] = sim
				for (j in 2:length(tabs))
					{
						for (k in 1:dim(tabs[[j]])[1])
							{
								if (!tabs[[j]][k,"node2"]%in%tabs[[j]][,"node1"])
									{
										startYear1 = round(tabs[[j]][k,"startYear"],4); endYear = round(tabs[[j]][k,"endYear"],4)
										startYear2 = tabs[[j]][which(tabs[[j]][,"node2"]==tabs[[j]][k,"node1"]),"startYear"]
										threeNodeDate = paste0(round(startYear2,4),"_",startYear1,"_",endYear)
										indices = which(threeNodeDates==threeNodeDate)
										if (length(indices) > 0)
											{
												if (length(unique(obs[tipBranches[indices],"endGenotype"])) == 1)
													{
														tabs[[j]][k,"endGenotype"] = obs[tipBranches[indices[1]],"endGenotype"]
													}
											}
									}
							}
					}
				for (j in 1:dim(tabs[[2]])[1])
					{
						if ((!tabs[[2]][j,"node2"]%in%tabs[[2]][,"node1"])&(is.na(tabs[[2]][j,"endGenotype"]))) print(j)
					}
			}
		for (j in 1:length(tabs))
			{
				noRemainingGenotypeToIdentify = FALSE
				while (noRemainingGenotypeToIdentify == FALSE)
					{
						for (k in 1:dim(tabs[[j]])[1])
							{
								if (is.na(tabs[[j]][k,"endGenotype"]))
									{
										indices = which(tabs[[j]][,"node1"]==tabs[[j]][k,"node2"])
										if (sum(!is.na(tabs[[j]][indices,"endGenotype"])) == 2)
											{
												if (tabs[[j]][indices[1],"endGenotype"] == tabs[[j]][indices[2],"endGenotype"])
													{
														tabs[[j]][indices[1],"startGenotype"] = tabs[[j]][indices[1],"endGenotype"]
														tabs[[j]][indices[2],"startGenotype"] = tabs[[j]][indices[2],"endGenotype"]
														tabs[[j]][k,"endGenotype"] = tabs[[j]][indices[1],"endGenotype"]
													}	else	{
														genotypeNames = tabs[[j]][indices,"endGenotype"]
														genotypeNames = genotypeNames[order(genotypeNames)]
														if ((genotypeNames[1]=="NY99")&(genotypeNames[2]=="WN02")) genotype = "NY99"
														if ((genotypeNames[1]=="SW03")&(genotypeNames[2]=="WN02")) genotype = "WN02"
														tabs[[j]][indices[1],"startGenotype"] = genotype
														tabs[[j]][indices[2],"startGenotype"] = genotype
														tabs[[j]][k,"endGenotype"] = genotype
													}
											}
									}
							}
						if (sum(is.na(tabs[[j]][,"startGenotype"])) == 2)
							{
								tabs[[j]][which(!tabs[[j]][,"node1"]%in%tabs[[j]][,"node2"]),"startGenotype"] = 
										 tabs[[j]][which(!tabs[[j]][,"node1"]%in%tabs[[j]][,"node2"]),"endGenotype"]
							}				
						if (sum(is.na(tabs[[j]][,"startGenotype"])) == 0) noRemainingGenotypeToIdentify = TRUE
					}
			}
		if (i == 0)
			{
				write.csv(tabs[[1]], "WNV_RRW_MCC.csv", row.names=F, quote=F)
			}	else		{
				write.csv(tabs[[1]], paste0("Tree_extractions/WNV_gamma_all/TreeExtractions_",i,".csv"), row.names=F, quote=F)
				write.csv(tabs[[2]], paste0("Tree_extractions/WNV_gamma_all/TreeSimulations_",i,".csv"), row.names=F, quote=F)
			}
	}

differentGenotypes = unique(genotypes[,"genotype"])
localTreesDirectory1 = "Tree_extractions/WNV_gamma_all"
for (i in 1:nberOfExtractionFiles)
	{
		tab1 = read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",i,".csv"), header=T)
		for (j in 1:length(differentGenotypes))
			{
				tab2 = tab1[which((tab1[,"startGenotype"]==differentGenotypes[j])&(tab1[,"endGenotype"]==differentGenotypes[j])),]
				write.csv(tab2, paste0(gsub("100",differentGenotypes[j],localTreesDirectory1),"/TreeExtractions_",i,".csv"), row.names=F, quote=F)
			}
		tab1 = read.csv(paste0(localTreesDirectory1,"/TreeSimulations_",i,".csv"), header=T)
		for (j in 1:length(differentGenotypes))
			{
				tab2 = tab1[which((tab1[,"startGenotype"]==differentGenotypes[j])&(tab1[,"endGenotype"]==differentGenotypes[j])),]
				write.csv(tab2, paste0(gsub("100",differentGenotypes[j],localTreesDirectory1),"/TreeSimulations_",i,".csv"), row.names=F, quote=F)
			}
	}

localTreesDirectory1 = "Tree_extractions/WNV_gamma_all"
localTreesDirectory2 = "Tree_extractions/WNV_gamma_bf02"
localTreesDirectory3 = "Tree_extractions/WNV_gamma_af02"
for (i in 1:nberOfExtractionFiles)
	{
		tab1 = read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",i,".csv"), header=T)
		tab2 = tab1[which(tab1[,"endYear"]<2002),]; tab3 = tab1[which(tab1[,"startYear"]>=2002),]
		write.csv(tab2, paste0(localTreesDirectory2,"/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		write.csv(tab3, paste0(localTreesDirectory3,"/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		tab1 = read.csv(paste0(localTreesDirectory1,"/TreeSimulations_",i,".csv"), header=T)
		tab2 = tab1[which(tab1[,"endYear"]<2002),]; tab3 = tab1[which(tab1[,"startYear"]>=2002),]
		write.csv(tab2, paste0(localTreesDirectory2,"/TreeSimulations_",i,".csv"), row.names=F, quote=F)
		write.csv(tab3, paste0(localTreesDirectory3,"/TreeSimulations_",i,".csv"), row.names=F, quote=F)
	}

localTreesDirectory1 = "Tree_extractions/WNV_gamma_all"
localTreesDirectory2 = "Tree_extractions/WNV_internal_100"
for (i in 1:nberOfExtractionFiles)
	{
		tab1 = read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",i,".csv"), header=T)
		tab2 = tab1; tab2 = tab2[which(tab2[,"node2"]%in%tab2[,"node1"]),]
		write.csv(tab2, paste0(localTreesDirectory2,"/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		tab1 = read.csv(paste0(localTreesDirectory1,"/TreeSimulations_",i,".csv"), header=T)
		tab2 = tab1; tab2 = tab2[which(tab2[,"node2"]%in%tab2[,"node1"]),]
		write.csv(tab2, paste0(localTreesDirectory2,"/TreeSimulations_",i,".csv"), row.names=F, quote=F)
	}
	
directories = list.files("Tree_extractions")
for (i in 1:length(directories))
	{
		for (j in 1:nberOfExtractionFiles)
			{
				tab1 = read.csv(paste0("Tree_extractions/",directories[i],"/TreeSimulations_",j,".csv"))
				if (!"treeID"%in%colnames(tab1))
					{
						tab2 = cbind(tab1,rep(j,dim(tab1)[1])); colnames(tab2) = c(colnames(tab1),"treeID")
						write.csv(tab2, paste0("Tree_extractions/",directories[i],"/TreeSimulations_",j,".csv"), row.names=F, quote=F)
					}
			}
	}

	# 4.2. Plotting MCC tree with tip nodes coloured by genotype

tree = readAnnotatedNexus("WNV_RRW_MCC.tree"); rootHeight = 18.4631
genotypes = read.table("WNV_genotypes.txt", header=T)
differentGenotypes = c("NY99","WN02","SW03")
cols1 = c("#DE4327","#FAA521","#4676BB"); cols2 = rep(NA, length(tree$tip.label))
for (i in 1:length(tree$tip.label))
	{
		index1 = which(row.names(genotypes)==gsub("'","",tree$tip.label[i]))
		index2 = which(differentGenotypes==genotypes[index1,1]); cols2[i] = cols1[index2]
	}
tree = ape::rotate(tree, node=813)
pdf("WNV_figures_&_SI/SI_files/MCC_tree_genotypes_NEW.pdf", width=8, height=10); par(oma=c(0,0,0,0), mar=c(1,1,0,1), lwd=0.2)
plot(tree, show.tip.label=F, show.node.label=F, edge.width=0.7, cex=0.6, align.tip.label=3, col="gray30", edge.color="gray30")
for (i in 1:dim(tree$edge)[1])
	{
		if (!tree$edge[i,2]%in%tree$edge[,1])
			{
				nodelabels(node=tree$edge[i,2], pch=16, cex=0.8, col=cols2[tree$edge[i,2]])
				nodelabels(node=tree$edge[i,2], pch=1, cex=0.8, col="gray30", lwd=0.5)
			}
	}
axis(at=c(-1998:-2017)+mostRecentSamplingDatum, labels=c(2017:1998), mgp=c(0,0,-0.9), cex.axis=0.6,
	 lwd=0.5, lwd.tick=0.5, tck=-0.006, side=1, line=-1, col.tick="gray30", col.axis="gray30", col="gray30")
dev.off()

# 5. Generating a dispersal history graph (by time slice or not)

	# 5.1. Building the maximum clade consensus (MCC) with TreeAnnotator
	
system(paste0("BEAST_ver_1.10.4/bin/treeannotator -burninTrees 0 -heights keep WNV_RRW_100.trees WNV_RRW_MCC.tree"))

	# 5.2. Extraction of the spatio-temporal information embedded in the MCC tree

source("mccTreeExtraction.r")
mcc_tre = readAnnotatedNexus("WNV_RRW_MCC.tree")
mcc_tab = mccTreeExtraction(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, "WNV_RRW_MCC.csv", row.names=F, quote=F)

	# 5.3. Generating sampling maps (coloured by sampling time and genotypes)

elevation = raster("Environmental_files/WNV_rasters/Elevation_WNV_04.asc"); e_Figure = extent(-128,-63,21,53)
lakes = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="Natural_Earth_lakes"), e_Figure)
borders = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="International_borders"), e_Figure)
background = crop(raster("Environmental_files/Natural_Earth/Gray_background.tif"), e_Figure)
background[background[]==106] = NA; r = background
background_cols = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
mcc = read.csv("WNV_RRW_MCC.csv", header=T); mcc = mcc[sample(1:dim(mcc)[1],dim(mcc)[1],replace=F),]
endYearsM = ((mcc[,"endYear"]-min(mcc[,"startYear"]))/(max(mcc[,"endYear"])-min(mcc[,"startYear"]))*100)+1
startYearM = ((mcc[1,"startYear"]-min(mcc[,"startYear"]))/(max(mcc[,"endYear"])-min(mcc[,"startYear"]))*100)+1
cols1 = colorRampPalette(brewer.pal(11,"RdYlGn"))(101); cols1_mcc = cols1[endYearsM]
differentGenotypes = c("NY99","WN02","SW03"); cols2 = c("#DE4327","#FAA521","#4676BB"); cols2_mcc = rep(NA, dim(mcc)[1])
for (i in 1:dim(mcc)[1]) cols2_mcc[i] = cols2[which(differentGenotypes==mcc[i,"endGenotype"])]
cols_mcc_list = list(); cols_mcc_list[[1]] = cols1_mcc; cols_mcc_list[[2]] = cols2_mcc
rast = elevation; rast[!is.na(rast[])] = NA; rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
xmin = -125; xmax = -70; ymin = 23; ymax = 51
labsX = c(expression(125*degree*W), expression(70*degree*W))
labsY = c(expression(23*degree*N), expression(51*degree*N))

for (h in 1:2)
	{
		dev.new(width=8.2, height=5.1); par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1,2.8,0.9,0.5), mgp=c(1,0.2,0), lwd=0.3)
		plot(background, main="", cex.main=1, bty="n", box=F, axes=F, legend=F, col=background_cols, colNA="#D0D0D0")
		plot(borders, add=T, lwd=1, col="white", lty=1); plot(lakes, add=T, lwd=0.7, col="#D0D0D0", border=NA)
		plot(rast, legend.only=T, add=T, col=cols1, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.891,0.897,0.096,0.908), adj=3,
			 axis.args=list(at=c(1999:2016), cex.axis=0.6, lwd=0, col="gray30", lwd.tick=0.5, col.tick="gray30", tck=-0.7, col.axis="gray30",
			 line=0, mgp=c(0,0.3,0)), alpha=1, side=3)
		for (i in dim(mcc)[1]:1)
			{
				if (!mcc[i,"node2"]%in%mcc[,"node1"])
					{
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=cols_mcc_list[[h]][i], cex=0.8)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", cex=0.8, lwd=0.1)
					}
			}
		axis(1, c(xmin,xmax), labels=labsX, pos=ymin(r), cex.axis=0.5, col.axis="gray30", lwd=0, lwd.tick=0.3, tck=0.008, col.tick="gray30", mgp=c(0,-1.1,0), col="gray30", asp=2)
		axis(2, c(ymin,ymax), labels=labsY, pos=xmin(r), cex.axis=0.5, col.axis="gray30", lwd=0, lwd.tick=0.5, tck=0.008, col.tick="gray30", mgp=c(0,-0.8,0), col="gray30", asp=2)
		rect(xmin(background), ymin(background), xmax(background), ymax(background), xpd=T, lwd=0.5, border="gray30")
	}

	# 5.4. Estimating the HPD region for each time slice

rast = raster("Environmental_files/WNV_rasters/Elevation_WNV_04.asc"); rast[!is.na(rast[])] = 0
precision = 1; nberOfCores = 1; origin = FALSE; startDatum = min(read.csv("WNV_RRW_MCC.csv", header=T)[,"startYear"])

prob = 0.95; spread = spreadGraphic1(localTreesDirectory, nberOfExtractionFiles, rast, prob, startDatum, precision, timeLayers=F, nberOfCores, origin)
spread[is.na(rast)] = NA; writeRaster(spread, paste0("WNV_spread_095.asc"))
prob = 0.95; spreads = spreadGraphic1(localTreesDirectory, nberOfExtractionFiles, rast, prob, startDatum, precision, timeLayers=T, nberOfCores, origin)
for (i in 1:length(spreads))
	{
		r = spreads[[i]]; r[is.na(rast)] = NA; writeRaster(r, paste0("WNV_spreads_095/WNV_spreads_095_",max(r[],na.rm=T),".asc"), overwrite=T)
	}
prob = 0.80; spread = spreadGraphic1(localTreesDirectory, nberOfExtractionFiles, rast, prob, startDatum, precision, timeLayers=F, nberOfCores, origin)
spread[is.na(rast)] = NA; writeRaster(spread, paste0("WNV_spread_080.asc"))
prob = 0.80; spreads = spreadGraphic1(localTreesDirectory, nberOfExtractionFiles, rast, prob, startDatum, precision, timeLayers=T, nberOfCores, origin)
for (i in 1:length(spreads))
	{
		r = spreads[[i]]; r[is.na(rast)] = NA; writeRaster(r, paste0("WNV_spreads_080/WNV_spreads_080_",max(r[],na.rm=T),".asc"), overwrite=T)
	}

	# 5.5. Plotting the global dispersal history graph (without and with HPD polygons)

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
cols = colorRampPalette(brewer.pal(11,"RdYlGn"))(101); cols_mcc = cols[endYearsM]; col_start = cols[startYearM]
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

pdf("WNV_figures_&_SI/SI_files/Dispersal_history_NEW.pdf", width=7.5, height=5); par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1,1.4,0.9,0), mgp=c(1,0.2,0), lwd=0.3)
# dev.new(width=7.5, height=5); par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1,1.4,0.9,0), mgp=c(1,0.2,0), lwd=0.3)
plot(background, main="", cex.main=1, bty="n", box=F, axes=F, legend=F, col="white", colNA="#D0D0D0")
plot(lakes, add=T, lwd=0.7, col="#D0D0D0", border=NA); # plot(borders, add=T, lwd=1, col="white"); # plot(flyways, add=T, lwd=1, lty=1, border="blue")
for (i in length(pols):1) plot(pols[[i]], axes=F, col=cols_pol[[i]], add=T, border=NA)
plot(background, main="", cex.main=1, bty="n", box=F, axes=F, legend=F, col=rgb(0,0,0,0), colNA="#D0D0D0", add=T)
plot(lakes, add=T, lwd=0.7, col="#D0D0D0", border=NA); plot(borders, add=T, lwd=1, col="#D0D0D0", lty=2)
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
axis(1, c(xmin,xmax), labels=labsX, pos=ymin(r), cex.axis=0.5, col.axis="gray30", lwd=0, lwd.tick=0.3, tck=0.01, col.tick="gray30", mgp=c(0,-1.2,0), col="gray30", asp=2)
axis(2, c(ymin,ymax), labels=labsY, pos=xmin(r), cex.axis=0.5, col.axis="gray30", lwd=0, lwd.tick=0.5, tck=0.01, col.tick="gray30", mgp=c(0,-0.9,0), col="gray30", asp=2)
rect(xmin(background), ymin(background), xmax(background), ymax(background), xpd=T, lwd=0.5, border="gray30")
plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.91,0.92,0.060,0.944), adj=3,
	 axis.args=list(at=c(1999:2016), col="gray30", cex.axis=0.6, col.axis="gray30", lwd=0, lwd.tick=0.5, tck=-0.5, col.tick="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)
dev.off()

	# 5.6. Plotting the MCC tree with internal and tip nodes coloured according to their genotype

elevation = raster("Environmental_files/WNV_rasters/Elevation_WNV_04.asc"); e_Figure = extent(-128,-63,21,53)
lakes = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="Natural_Earth_lakes"), e_Figure)
borders = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="International_borders"), e_Figure)
background = crop(raster("Environmental_files/Natural_Earth/Gray_background.tif"), e_Figure)
background[background[]==106] = NA; r = background
background_cols = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
mcc = read.csv("WNV_RRW_MCC.csv", header=T); mcc = mcc[order(mcc[,"startYear"]),]
mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]; mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
differentGenotypes = c("NY99","WN02","SW03"); cols = c("#DE4327","#FAA521","#4676BB"); cols_mcc = rep(NA, dim(mcc)[1])
for (i in 1:dim(mcc)[1]) cols_mcc[i] = cols[which(differentGenotypes==mcc[i,"endGenotype"])]
xmin = -125; xmax = -66; ymin = 23; ymax = 51
labsX = c(expression(125*degree*W), expression(66*degree*W))
labsY = c(expression(23*degree*N), expression(51*degree*N))

dev.new(width=8.2, height=5.1); par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(1,2.8,0.9,0.5), mgp=c(1,0.2,0), lwd=0.3)
plot(background, main="", cex.main=1, bty="n", box=F, axes=F, legend=F, col=background_cols, colNA="#D0D0D0")
plot(borders, add=T, lwd=1, col="white", lty=1); plot(lakes, add=T, lwd=0.7, col="#D0D0D0", border=NA)
for (i in dim(mcc)[1]:1)
	{
		curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
	}
for (i in dim(mcc)[1]:1)
	{
		points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=cols_mcc[i], cex=0.8)
		points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", cex=0.8, lwd=0.1)
	}
axis(1, c(xmin,xmax), labels=labsX, pos=ymin(r), cex.axis=0.5, col.axis="gray30", lwd=0, lwd.tick=0.3, tck=0.008, col.tick="gray30", mgp=c(0,-1.1,0), col="gray30", asp=2)
axis(2, c(ymin,ymax), labels=labsY, pos=xmin(r), cex.axis=0.5, col.axis="gray30", lwd=0, lwd.tick=0.5, tck=0.008, col.tick="gray30", mgp=c(0,-0.8,0), col="gray30", asp=2)
rect(xmin(background), ymin(background), xmax(background), ymax(background), xpd=T, lwd=0.5, border="gray30")

	# 5.7. Plotting dispersal history graph snapshots

years = c(1999:2016); cols = colorRampPalette(brewer.pal(11,"RdYlGn"))(101)
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
		pdf(paste0("WNV_figures_&_SI/SI_files/Dispersal_years_",h,"_NEW.pdf"), width=7.5, height=5.6); par(mfrow=c(3,3),oma=c(2,2.5,1,0.3),mar=c(0,0,0,0),mgp=c(1,0.2,0),lwd=0.2)
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
						axis(1, c(xmin,xmax), labels=labsX, pos=ymin(rast), cex.axis=0.6, col.axis="gray30",
							 lwd=0, lwd.tick=0.2, tck=-0.02, col.tick="gray30", mgp=c(0,0.10,0), col="gray30")
						axis(2, c(ymin,ymax), labels=labsY, pos=xmin(rast), cex.axis=0.6, col.axis="gray30",
							 lwd=0, lwd.tick=0.2, tck=-0.02, col.tick="gray30", mgp=c(0,0.27,0), col="gray30")
						mtext(years[i]+1, side=1, adj=0, line=-2.3, at=-126, cex=0.6, font=1, col="gray30")-
						rect(xmin(rast), ymin(rast), xmax(rast), ymax(rast), lwd=0.2, border="gray30")
					}
			}
		dev.off()
	}

# 6. Estimating and plotting some dispersal statistics

	# 6.1. Estimating dispersal statistics

timeSlices = 200; onlyTipBranches = F; showingPlots = F; nberOfCores = 10; slidingWindow = 1
localTreesDirectory = "Tree_extractions/WNV_gamma_all"; outputName = "WNV_all_obs"; simulations = FALSE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_bf02"; outputName = "WNV_bf02_obs"; simulations = FALSE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_af02"; outputName = "WNV_af02_obs"; simulations = FALSE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_NY99"; outputName = "WNV_NY99_obs"; simulations = FALSE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_WN02"; outputName = "WNV_WN02_obs"; simulations = FALSE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_SW03"; outputName = "WNV_SW03_obs"; simulations = FALSE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_all"; outputName = "WNV_all_sim"; simulations = TRUE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_bf02"; outputName = "WNV_bf02_sim"; simulations = TRUE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_af02"; outputName = "WNV_af02_sim"; simulations = TRUE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_NY99"; outputName = "WNV_NY99_sim"; simulations = TRUE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_WN02"; outputName = "WNV_WN02_sim"; simulations = TRUE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_SW03"; outputName = "WNV_SW03_sim"; simulations = TRUE
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow, simulations)
localTreesDirectory = "Tree_extractions/WNV_gamma_all"

	# 6.2. Plotting dispersal statistics for the inferred trees (for Figure 1)

stats = read.table("Dispersal_statistics/WNV_all_obs_estimated_dispersal_statistics.txt", header=T)
v = stats[,"mean_branch_dispersal_velocity"]; D = stats[,"original_diffusion_coefficient"]; D = D/365
cat("Mean branch dispersal velocity:"); cat("\n"); cat("\t")
cat(paste0(round(mean(v),1)," km/year, 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]"))
cat("Original diffusion coefficient:"); cat("\n"); cat("\t")
cat(paste0(round(mean(D),1)," km2/day, 95% HPD = [",round(quantile(D,0.025),1),",",round(quantile(D,0.975),1),"]"))
v = stats[,"weighted_branch_dispersal_velocity"]; D = stats[,"weighted_diffusion_coefficient"]
cat("Weighted branch dispersal velocity:"); cat("\n"); cat("\t")
cat(paste0(round(mean(v),1)," km/year, 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]"))
cat("Weighted diffusion coefficient:"); cat("\n"); cat("\t")
cat(paste0(round(mean(D))," km2/year, 95% HPD = [",round(quantile(D,0.025)),",",round(quantile(D,0.975)),"]"))

dev.new(width=4, height=1.5); par(mgp=c(0,0,0), oma=c(0.8,0.8,0,0), mar=c(1.5,1.5,1,1))
col1 = rgb(100, 100, 100, 255, maxColorValue=255); col2 = rgb(100, 100, 100, 100, maxColorValue=255)
tab1 = read.table("Dispersal_statistics/WNV_all_obs_median_spatial_wavefront_distance.txt", header=T)
tab2 = read.table("Dispersal_statistics/WNV_all_obs_95%HPD_spatial_wavefront_distance.txt", header=T)
plot(tab1[,1], tab1[,2], type="l", axes=F, ann=F, ylim=c(0,4300), xlim=c(1998,2010), col=NA)
slicedTimes = tab1[,1]; waveFrontDistances1MedianValue = tab1[,2]; lower_l_1 = tab2[,2]; upper_l_1 = tab2[,3]
xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=col2, border=0)
lines(slicedTimes, waveFrontDistances1MedianValue, lwd=1, col=col1)
axis(side=1, pos=0, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.050, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, pos=1998, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.045, col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab="time (year)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
title(ylab="distance (km)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
dev.copy2pdf(file="WNV_figures_&_SI/Figures/Figure1_spatial_wavefronts_TEMP.pdf")

dev.new(width=4, height=4); par(mgp=c(0,0,0), oma=c(0.8,0.8,0,0), mar=c(1.5,1.5,1,1))
col1 = rgb(100, 100, 100, 255, maxColorValue=255); col2 = rgb(100, 100, 100, 100, maxColorValue=255)
tab1 = read.table("Dispersal_statistics/WNV_all_obs_median_mean_branch_dispersal_velocity.txt", header=T)
tab2 = read.table("Dispersal_statistics/WNV_all_obs_95%HPD_mean_branch_dispersal_velocity.txt", header=T)
plot(tab1[,1], tab1[,2], type="l", axes=F, ann=F, ylim=c(0,10000), xlim=c(1998,2010), col=NA)
slicedTimes = tab1[,1]; branchDispersalVelocityMeanValue = tab1[,2]; lower_l_1 = tab2[,2]; upper_l_1 = tab2[,3]
xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=col2, border=0)
lines(slicedTimes, branchDispersalVelocityMeanValue, lwd=1, col=col1)
axis(side=1, pos=0, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.050, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, pos=1998, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.045, col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab="time (year)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
title(ylab="mean branch dispersal velocity", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
dev.copy2pdf(file="WNV_figures_&_SI/Figures/Figure1_dispersal_velocity_TEMP.pdf")

dev.new(width=4, height=3); par(mgp=c(0,0,0), oma=c(0.8,0.8,0,0), mar=c(1.5,1.5,1,1))
col1 = rgb(100, 100, 100, 255, maxColorValue=255); col2 = rgb(100, 100, 100, 100, maxColorValue=255)
tab1 = read.table("Dispersal_statistics/WNV_all_obs_median_weighted_branch_dispersal_velocity.txt", header=T)
tab2 = read.table("Dispersal_statistics/WNV_all_obs_95%HPD_weighted_branch_dispersal_velocity.txt", header=T)
plot(tab1[,1], tab1[,2], type="l", axes=F, ann=F, ylim=c(0,600), xlim=c(1998,2010), col=NA)
slicedTimes = tab1[,1]; branchDispersalVelocityMeanValue = tab1[,2]; lower_l_1 = tab2[,2]; upper_l_1 = tab2[,3]
xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=col2, border=0)
lines(slicedTimes, branchDispersalVelocityMeanValue, lwd=1, col=col1)
axis(side=1, pos=0, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.050, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, pos=1998, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.045, col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab="time (year)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
title(ylab="weighted branch dispersal velocity", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
dev.copy2pdf(file="WNV_figures_&_SI/Figures/Figure1_dispersal_velocity_TEMP.pdf")

	# 6.3. Comparing the evolution of dispersal velocity between inferred and simulated trees

dev.new(width=4, height=4); par(mgp=c(0,0,0), oma=c(0.8,0.8,0,0), mar=c(1.5,1.5,1,1)); cols1 = list(); cols2 = list()
cols1[[1]] = rgb(204,0,0,255,maxColorValue=255); cols2[[1]] = rgb(204,0,0,100,maxColorValue=255) # red
cols1[[2]] = rgb(120,120,120,255,maxColorValue=255); cols2[[2]] = rgb(120,120,120,100,maxColorValue=255)
tab1a = read.table("Dispersal_statistics/WNV_all_obs_median_mean_branch_dispersal_velocity.txt", header=T)
tab2a = read.table("Dispersal_statistics/WNV_all_obs_95%HPD_mean_branch_dispersal_velocity.txt", header=T)
tab1b = read.table("Dispersal_statistics/WNV_all_sim_median_mean_branch_dispersal_velocity.txt", header=T)
tab2b = read.table("Dispersal_statistics/WNV_all_sim_95%HPD_mean_branch_dispersal_velocity.txt", header=T)
plot(tab1a[,1], tab1a[,2], type="l", axes=F, ann=F, ylim=c(0,10000), xlim=c(1998,2010), col=NA)
slicedTimes = tab1a[,1]; branchDispersalVelocityMeanValue = tab1a[,2]; lower_l_1 = tab2a[,2]; upper_l_1 = tab2a[,3]
xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=cols2[[1]], border=0)
lines(slicedTimes, branchDispersalVelocityMeanValue, lwd=1, col=cols1[[1]])
slicedTimes = tab1b[,1]; branchDispersalVelocityMeanValue = tab1b[,2]; lower_l_1 = tab2b[,2]; upper_l_1 = tab2b[,3]
xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=cols2[[2]], border=0)
lines(slicedTimes, branchDispersalVelocityMeanValue, lwd=1, col=cols1[[2]])
axis(side=1, pos=0, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.050, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, pos=1998, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.045, col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab="time (year)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
title(ylab="mean branch dispersal velocity", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")

	# 6.4. Comparing dispersal statistics before and after 2002

subsets = c("bf02","af02")
for (i in 1:length(subsets))
	{
		stats = read.table(paste0("Dispersal_statistics/WNV_",subsets[i],"_obs_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"mean_branch_dispersal_velocity"]; D = stats[,"original_diffusion_coefficient"]; D = D/365
		cat("Mean branch dispersal velocity:"); cat("\n"); cat("\t")
		cat(paste0(round(mean(v),1)," km/year, 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]\n"))
		cat("Original diffusion coefficient:"); cat("\n"); cat("\t")
		cat(paste0(round(mean(D),1)," km2/day, 95% HPD = [",round(quantile(D,0.025),1),",",round(quantile(D,0.975),1),"]\n"))
		v = stats[,"weighted_branch_dispersal_velocity"]; D = stats[,"weighted_diffusion_coefficient"]
		cat("Weighted branch dispersal velocity:"); cat("\n"); cat("\t")
		cat(paste0(round(mean(v),1)," km/year, 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]\n"))
		cat("Weighted diffusion coefficient:"); cat("\n"); cat("\t")
		cat(paste0(round(mean(D))," km2/year, 95% HPD = [",round(quantile(D,0.025)),",",round(quantile(D,0.975)),"]\n\n"))
	}

	# 6.5. Comparing dispersal statistics between the different subsets
	
analyses = c("NY99","WN02","SW03","bf02","af02","all"); cols1 = list(); cols2 = list()
cols1[[1]] = rgb(222,67,39,255,maxColorValue=255); cols2[[1]] = rgb(222,67,39,100,maxColorValue=255) # red
cols1[[2]] = rgb(250,165,33,255,maxColorValue=255); cols2[[2]] = rgb(250,165,33,100,maxColorValue=255) # orange
cols1[[3]] = rgb(70,118,187,255,maxColorValue=255); cols2[[3]] = rgb(70,118,187,100,maxColorValue=255) # blue
cols1[[4]] = rgb(150,150,150,255,maxColorValue=255); cols2[[4]] = rgb(150,150,150,100,maxColorValue=255) # light grey
cols1[[5]] = rgb(76,76,76,255,maxColorValue=255); cols2[[5]] = rgb(60,60,60,100,maxColorValue=255) # dark grey
cols1[[6]] = rgb(77,77,77,255,maxColorValue=255); cols2[[6]] = rgb(0,0,0,0,maxColorValue=255) # black (and transparent)

for (i in length(analyses):1)
	{
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_obs_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"mean_branch_dispersal_velocity"]; D = stats[,"original_diffusion_coefficient"]; D = D/365
		cat(paste0("Mean branch dispersal velocity for genotype ",analyses[i],":\n\t"))
		cat(paste0(round(mean(v),1)," km/year, 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]\n"))
		if (i == length(analyses))
			{
				dev.new(width=4.5, height=2.5); par(mfrow=c(1,1), mgp=c(1,0.35,0), oma=c(0.5,1,1,1), mar=c(2.5,2,0.5,0), lwd=0.2)
				plot(density(v), lwd=0.7, col=NA, axes=F, ann=F, xlim=c(0,7000), ylim=c(0,0.002))
			}
		polygon(density(v), col=cols2[[i]], border=NA)
	}
for (i in length(analyses):1)
	{
		if (analyses[i] != "all") { LTY = 1; LWD = 0.7 }
		if (analyses[i] == "all") { LTY = 3; LWD = 1.2 }
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_obs_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"mean_branch_dispersal_velocity"]; lines(density(v), lwd=LWD, col=cols1[[i]], lty=LTY)
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_sim_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"mean_branch_dispersal_velocity"]; # lines(density(v), lwd=1.2, col=cols1[[i]], lty=3)
	}
axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.02,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30", at=c(0,0.0008,0.0016), labels=c(0,0.0008,0.0016))
title(xlab="mean lineage dispersal velocity (km/year)", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
for (i in length(analyses):1)
	{
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_obs_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"weighted_branch_dispersal_velocity"]; D = stats[,"original_diffusion_coefficient"]; D = D/365
		cat(paste0("Weighted branch dispersal velocity for genotype ",analyses[i],":\n\t"))
		cat(paste0(round(mean(v),1)," km/year, 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]\n"))
		if (i == length(analyses))
			{
				dev.new(width=4.5, height=2.5); par(mfrow=c(1,1), mgp=c(1,0.35,0), oma=c(0.5,1,1,1), mar=c(2.5,2,0.5,0), lwd=0.2)
				plot(density(v), lwd=0.7, col=NA, axes=F, ann=F, xlim=c(120,500), ylim=c(0,0.14))
			}
		polygon(density(v), col=cols2[[i]], border=NA)
	}
for (i in length(analyses):1)
	{
		if (analyses[i] != "all") { LTY = 1; LWD = 0.7 }
		if (analyses[i] == "all") { LTY = 3; LWD = 1.2 }
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_obs_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"weighted_branch_dispersal_velocity"]; lines(density(v), lwd=LWD, col=cols1[[i]], lty=LTY)
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_sim_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"weighted_branch_dispersal_velocity"]; # lines(density(v), lwd=1.2, col=cols1[[i]], lty=3)
	}
axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.02,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab="weighted lineage dispersal velocity (km/year)", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")

for (i in length(analyses):1)
	{
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_obs_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"original_diffusion_coefficient"]; D = stats[,"original_diffusion_coefficient"]; D = D/365
		cat(paste0("Original diffusion coefficient for genotype ",analyses[i],":\n\t"))
		cat(paste0(round(mean(v),1)," km2/year, 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]\n"))
		if (i == length(analyses))
			{
				dev.new(width=5, height=2.5); par(mfrow=c(1,1), mgp=c(1,0.35,0), oma=c(0.5,1,1,1), mar=c(2.5,2,0.5,0), lwd=0.2)
				plot(density(v), lwd=0.7, col=NA, axes=F, ann=F, xlim=c(0,3*10^6), ylim=c(0,0.000005))
			}
		polygon(density(v), col=cols2[[i]], border=NA)
	}
for (i in length(analyses):1)
	{
		if (analyses[i] != "all") { LTY = 1; LWD = 0.7 }
		if (analyses[i] == "all") { LTY = 3; LWD = 1.2 }
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_obs_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"original_diffusion_coefficient"]; lines(density(v), lwd=LWD, col=cols1[[i]], lty=LTY)
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_sim_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"original_diffusion_coefficient"]; # lines(density(v), lwd=1.2, col=cols1[[i]], lty=3)
	}
axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.02,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab="mean lineage dispersal velocity (km/year)", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
for (i in length(analyses):1)
	{
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_obs_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"weighted_diffusion_coefficient"]; D = stats[,"original_diffusion_coefficient"]; D = D/365
		cat(paste0("Weighted diffusion coefficient for genotype ",analyses[i],":\n\t"))
		cat(paste0(round(mean(v),1)," km2/year, 95% HPD = [",round(quantile(v,0.025),1),",",round(quantile(v,0.975),1),"]\n"))
		if (i == length(analyses))
			{
				dev.new(width=6, height=2.5); par(mfrow=c(1,1), mgp=c(1,0.35,0), oma=c(0.5,1,1,1), mar=c(2.5,2,0.5,0), lwd=0.2)
				plot(density(v), lwd=0.7, col=NA, axes=F, ann=F, xlim=c(0,1.5*10^5), ylim=c(0,0.00025))
			}
		polygon(density(v), col=cols2[[i]], border=NA)
	}
for (i in length(analyses):1)
	{
		if (analyses[i] != "all") { LTY = 1; LWD = 0.7 }
		if (analyses[i] == "all") { LTY = 1; LWD = 1.2 }
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_obs_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"weighted_diffusion_coefficient"]; lines(density(v), lwd=LWD, col=cols1[[i]], lty=LTY)
		stats = read.table(paste0("Dispersal_statistics/WNV_",analyses[i],"_sim_estimated_dispersal_statistics.txt"), header=T)
		v = stats[,"weighted_diffusion_coefficient"]; # lines(density(v), lwd=1.2, col=cols1[[i]], lty=3)
	}
axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.02,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab="weighted lineage dispersal velocity (km/year)", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")

# 7. Testing the impact of environmental factors on lineage dispersal locations

startTime = 1999; endTime = 2016+(8/12); months = c("01","02","03","04","05","06","07","08","09","10","11","12")
envVariables = list(); timeSlices = (endTime-startTime)*12; slidingWindow = 1/12; showingPlots = FALSE; nberOfCores = 50
envVariableNames = c("Elevation","Land_cover_forests","Land_cover_shrublands","Land_cover_savannas","Land_cover_grasslands",
					 "Land_cover_croplands","Land_cover_urban_areas","Temperature_interpolated_rasters","Precipitation_interpolated_rasters")
envVariableTitles = c("Elevation","Forests","Shrublands","Savannas","Grasslands","Croplands","Urban areas","Temperature","Precipitation")
for (i in 1:length(envVariableNames))
	{
		if (grepl("interpolated",envVariableNames[i]))
			{
				rasters_list = list(); raster_files = list.files(paste0("Environmental_files/WNV_rasters/",envVariableNames[i],"/"))
				for (j in 1:length(raster_files))
					{
						rast = raster(paste0("Environmental_files/WNV_rasters/",envVariableNames[i],"/",raster_files[j]))
						year1 = unlist(strsplit(raster_files[j],"_"))[1]; month1 = gsub("month","",unlist(strsplit(raster_files[j],"_"))[2])
						startDate = decimal_date(ymd(paste(year1,month1,"01",sep="-")))
						month_index = which(months==month1)
						if (month_index < 12)
							{
								month2 = months[month_index+1]; endDate = decimal_date(ymd(paste(year1,month2,"01",sep="-")))
							}	else	{
								month2 = "01"; year2 = as.numeric(year1)+1; endDate = decimal_date(ymd(paste(year2,month2,"01",sep="-")))
							}
						names(rast) = paste(unlist(strsplit(envVariableNames[[i]],"_"))[1],startDate,endDate,sep="_"); rasters_list[[j]] = rast
					}
				envVariables[[i]] = stack(rasters_list)
			}	else		{
				envVariables[[i]] = raster(paste0("Environmental_files/WNV_rasters/",envVariableNames[[i]],"_WNV_08.asc"))
			}
	}

analyses = c("WNV_gamma_all","WNV_internal_all","WNV_gamma_bf02","WNV_gamma_af02","WNV_gamma_NY99","WNV_gamma_WN02","WNV_gamma_SW03")
for (i in 1:1)
	{
		localTreesDirectory = paste0("Tree_extractions/",analyses[i])
		for (j in 1:nberOfExtractionFiles)
			{
				obs = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), header=T)
				sim = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",j,".csv"), header=T)
				envValues_obs = matrix(nrow=dim(obs)[1], ncol=length(envVariables))
				envValues_sim = matrix(nrow=dim(sim)[1], ncol=length(envVariables))
				colnames(envValues_obs) = envVariableNames
				colnames(envValues_sim) = envVariableNames
				for (k in 1:length(envVariables))
					{
						if (dim(envVariables[[k]])[3] > 1)
							{
								time_intervals = matrix(nrow=length(names(envVariables[[k]])), ncol=2)
								for (l in 1:length(names(envVariables[[k]])))
									{
										time_intervals[l,1] = as.numeric(unlist(strsplit(names(envVariables[[k]])[l],"_"))[2])
										time_intervals[l,2] = as.numeric(unlist(strsplit(names(envVariables[[k]])[l],"_"))[3])
									}
							}
						if (dim(envVariables[[k]])[3] == 1)
							{
								envValues_obs[,k] = raster::extract(envVariables[[k]], SpatialPoints(obs[,c("endLon","endLat")]))
								envValues_sim[,k] = raster::extract(envVariables[[k]], SpatialPoints(sim[,c("endLon","endLat")]))
							}	else	{
								indices_obs = which(obs[,"endYear"]<time_intervals[l,1])
								indices_sim = which(sim[,"endYear"]<time_intervals[l,1])
								if (length(indices_obs) > 0)
									{
										envValues_obs[indices_obs,k] = raster::extract(envVariables[[k]][[l]], SpatialPoints(obs[indices_obs,c("endLon","endLat")]))
									}
								if (length(indices_sim) > 0)
									{
										envValues_sim[indices_sim,k] = raster::extract(envVariables[[k]][[l]], SpatialPoints(sim[indices_sim,c("endLon","endLat")]))					
									}
								for (l in 1:dim(time_intervals)[1])
									{
										indices_obs = which((obs[,"endYear"]>time_intervals[l,1])&(obs[,"endYear"]<=time_intervals[l,2]))
										indices_sim = which((sim[,"endYear"]>time_intervals[l,1])&(sim[,"endYear"]<=time_intervals[l,2]))
										if (length(indices_obs) > 0)
											{
												envValues_obs[indices_obs,k] = raster::extract(envVariables[[k]][[l]], SpatialPoints(obs[indices_obs,c("endLon","endLat")]))
											}
										if (length(indices_sim) > 0)
											{
												envValues_sim[indices_sim,k] = raster::extract(envVariables[[k]][[l]], SpatialPoints(sim[indices_sim,c("endLon","endLat")]))			
											}
									}
							}
					}
				write.csv(envValues_obs, paste0(localTreesDirectory,"/EnvValues_obs_",j,".csv"), row.names=F, quote=F)
				write.csv(envValues_sim, paste0(localTreesDirectory,"/EnvValues_sim_",j,".csv"), row.names=F, quote=F)
			}
	}
for (i in 2:length(analyses))
	{
		localTreesDirectory = paste0("Tree_extractions/",analyses[i])
		for (j in 1:nberOfExtractionFiles)
			{
				envValues_obs1 = read.csv(paste0("Tree_extractions/WNV_gamma_all/EnvValues_obs_",j,".csv"), header=T)
				envValues_sim1 = read.csv(paste0("Tree_extractions/WNV_gamma_all/EnvValues_sim_",j,".csv"), header=T)
				obs1 = read.csv(paste0("Tree_extractions/WNV_gamma_all/TreeExtractions_",j,".csv"), header=T)
				sim1 = read.csv(paste0("Tree_extractions/WNV_gamma_all/TreeSimulations_",j,".csv"), header=T)
				obs2 = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), header=T)
				sim2 = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",j,".csv"), header=T)
				indices1 = c(); indices2 = c()
				for (k in 1:dim(obs1)[1])
					{
						if (sum((obs2[,"node1"]==obs1[k,"node1"])&(obs2[,"node2"]==obs1[k,"node2"])) == 1) indices1 = c(indices1, k)
					}
				for (k in 1:dim(sim1)[1])
					{
						if (sum((sim2[,"node1"]==sim1[k,"node1"])&(sim2[,"node2"]==sim1[k,"node2"])) == 1) indices2 = c(indices2, k)
					}
				envValues_obs2 = envValues_obs1[indices1,]; envValues_sim2 = envValues_sim1[indices2,]
				write.csv(envValues_obs2, paste0(localTreesDirectory,"/EnvValues_obs_",j,".csv"), row.names=F, quote=F)
				write.csv(envValues_sim2, paste0(localTreesDirectory,"/EnvValues_sim_",j,".csv"), row.names=F, quote=F)
			}
	}
BFs_lower = matrix(nrow=length(envVariableNames), ncol=length(analyses))
BFs_higher = matrix(nrow=length(envVariableNames), ncol=length(analyses))
row.names(BFs_lower) = envVariableNames; colnames(BFs_lower) = analyses
row.names(BFs_higher) = envVariableNames; colnames(BFs_higher) = analyses
meanEnvValues_obs_list1 = list(); meanEnvValues_sim_list1 = list()
for (i in 1:length(envVariableNames))
	{
		meanEnvValues_obs_list2 = list(); meanEnvValues_sim_list2 = list()
		for (j in 1:length(analyses))
			{
				meanEnvValues_obs = rep(NA, nberOfExtractionFiles); meanEnvValues_sim = rep(NA, nberOfExtractionFiles)
				localTreesDirectory = paste0("Tree_extractions/",analyses[j]); lowerEnvValues = 0
				for (k in 1:nberOfExtractionFiles)
					{
						meanEnvValues_obs[k] = mean(read.csv(paste0(localTreesDirectory,"/EnvValues_obs_",k,".csv"))[,envVariableNames[i]], na.rm=T)
						meanEnvValues_sim[k] = mean(read.csv(paste0(localTreesDirectory,"/EnvValues_sim_",k,".csv"))[,envVariableNames[i]], na.rm=T)
						if (meanEnvValues_obs[k] < meanEnvValues_sim[k]) lowerEnvValues = lowerEnvValues+1				
					}
				p = lowerEnvValues/nberOfExtractionFiles; BFs_lower[i,j] = round((p/(1-p))/(0.5/(1-0.5)),1)
				p = (1-lowerEnvValues/nberOfExtractionFiles); BFs_higher[i,j] = round((p/(1-p))/(0.5/(1-0.5)),1)
				meanEnvValues_obs_list2[[j]] = meanEnvValues_obs; meanEnvValues_sim_list2[[j]] = meanEnvValues_sim
			}
		meanEnvValues_obs_list1[[i]] = meanEnvValues_obs_list2; meanEnvValues_sim_list1[[i]] = meanEnvValues_sim_list2
	}
write.csv(BFs_lower, "EnvValues_BF_lower.csv", quote=F); write.csv(BFs_higher, "EnvValues_BF_higher.csv", quote=F)

cols1 = list(); cols2 = list(); plottingNullDispersalModelCurve = FALSE
cols1[[1]] = rgb(77,77,77,255,maxColorValue=255); cols2[[1]] = rgb(0,0,0,0,maxColorValue=255) # black (and transparent)
cols1[[3]] = rgb(150,150,150,255,maxColorValue=255); cols2[[3]] = rgb(150,150,150,100,maxColorValue=255) # light grey
cols1[[4]] = rgb(76,76,76,255,maxColorValue=255); cols2[[4]] = rgb(60,60,60,100,maxColorValue=255) # dark grey
cols1[[5]] = rgb(222,67,39,255,maxColorValue=255); cols2[[5]] = rgb(222,67,39,100,maxColorValue=255) # red
cols1[[6]] = rgb(250,165,33,255,maxColorValue=255); cols2[[6]] = rgb(250,165,33,100,maxColorValue=255) # orange
cols1[[7]] = rgb(70,118,187,255,maxColorValue=255); cols2[[7]] = rgb(70,118,187,100,maxColorValue=255) # blue
dev.new(width=7.5, height=4); par(mfrow=c(3,3), oma=c(3,3,1.0,0.3), mar=c(1,1,1,1.5), mgp=c(1,0.2,0), lwd=0.2)
for (i in 1:length(envVariableNames))
	{
		xMin = 9999; xMax = 0; yMin = 9999; yMax = 0
		for (j in 1:length(analyses))
			{
				if (j == 1)
					{
						xMin = min(meanEnvValues_obs_list1[[i]][[j]]); xMax = max(meanEnvValues_obs_list1[[i]][[j]])
						yMin = min(density(meanEnvValues_obs_list1[[i]][[j]])$y); yMax = max(density(meanEnvValues_obs_list1[[i]][[j]])$y)
					}	else	{
						if (xMin > min(meanEnvValues_obs_list1[[i]][[j]])) xMin = min(meanEnvValues_obs_list1[[i]][[j]])
						if (xMax < max(meanEnvValues_obs_list1[[i]][[j]])) xMax = max(meanEnvValues_obs_list1[[i]][[j]])
						if (yMin > min(density(meanEnvValues_obs_list1[[i]][[j]])$y)) yMin = min(density(meanEnvValues_obs_list1[[i]][[j]])$y)
						if (yMax < max(density(meanEnvValues_obs_list1[[i]][[j]])$y))	yMax = max(density(meanEnvValues_obs_list1[[i]][[j]])$y)				
					}	
			}
		if (i == 3) yMax = 2
		for (j in 1:length(analyses))
			{
				if (analyses[j] != "WNV_gamma_all") { LWD = 0.7; LTY=1 }
				if (analyses[j] == "WNV_gamma_all") { LWD = 1.2; LTY=3 }
				if (j == 1) plot(density(meanEnvValues_obs_list1[[i]][[j]]), lwd=0, col=cols2[[j]], xlim=c(xMin,xMax), ylim=c(yMin,yMax), axes=F, ann=F)
				if (j >= 2) polygon(density(meanEnvValues_obs_list1[[i]][[j]]), border=NA, col=cols2[[j]])
				if (plottingNullDispersalModelCurve == TRUE) lines(density(meanEnvValues_sim_list1[[i]][[j]]), lwd=LWD, lty=3, col=cols1[[j]])
				lines(density(meanEnvValues_obs_list1[[i]][[j]]), lwd=LWD, lty=LTY, col=cols1[[j]])
			}
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.02,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
		if (i%in%c(7,8,9)) title(xlab="mean environmental value at tree node positions", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
		if (i%in%c(1,4,7)) title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
		mtext(envVariableTitles[i], line=-1, cex=0.55, font=1, col="gray30")
	}

# 8. Testing the impact of environmental factors on lineage dispersal tendency

pathModel = 0; envVariables = list(); resistances = list(); avgResistances = list(); fourCells = FALSE
nberOfRandomisations = 1; randomProcedure = 2; showingPlots = FALSE; nberOfCores = 1; OS = "Unix"; simulations = FALSE
envVariableNames = c("Elevation","Land_cover_forests","Land_cover_shrublands","Land_cover_savannas","Land_cover_grasslands",
					 "Land_cover_croplands","Annual_mean_temperature","Annual_precipitation"); c = 0
for (i in 1:length(envVariableNames))
	{
		rast = raster(paste("Environmental_files/WNV_rasters/",envVariableNames[i],"_WNV_16.asc",sep=""))
		c = c+1; envVariables[[c]] = rast; names(envVariables[[c]]) = envVariableNames[i]
		resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
	}
for (i in 1:length(envVariableNames))
	{
		rast = raster(paste("Environmental_files/WNV_rasters/",envVariableNames[i],"_WNV_16.asc",sep=""))
		c = c+1; envVariables[[c]] = rast; names(envVariables[[c]]) = envVariableNames[i]
		resistances[[c]] = FALSE; avgResistances[[c]] = FALSE
	}
analyses = c("WNV_gamma_all","WNV_internal_all","WNV_gamma_bf02","WNV_gamma_af02","WNV_gamma_NY99","WNV_gamma_WN02","WNV_gamma_SW03")
for (i in 1:length(analyses))
	{
		localTreesDirectory = paste0("Tree_extractions/",analyses[i]); outputName = paste0(analyses[i],"_dispersal_tendency")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
			 		  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F)
	}

# 9. Testing the impact of environmental factors on lineage dispersal velocity

c = 0; envVariables = list(); resistances = list(); avgResistances = list(); fourCells = FALSE
nberOfRandomisations = 0; randomProcedure = 3; showingPlots = FALSE; nberOfCores = 50; OS = "Unix"; randomisations = FALSE
envVariableNames = c("Elevation","Land_cover_forests","Land_cover_shrublands","Land_cover_savannas","Land_cover_grasslands",
	 				 "Land_cover_croplands","Land_cover_urban_areas","Annual_mean_temperature","Annual_precipitation")
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
analyses = c("WNV_gamma_all","WNV_gamma_bf02","WNV_gamma_af02","WNV_gamma_NY99","WNV_gamma_WN02","WNV_gamma_SW03")
for (i in 1:length(analyses))
	{
		localTreesDirectory = paste0("Tree_extractions/",analyses[i])
		pathModel = 2; simulations = FALSE; outputNames = paste0(analyses[i],"_least-cost_extractions")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
		 nberOfRandomisations,randomProcedure,outputName=analyses[i],showingPlots,nberOfCores,OS,simulations=F,randomisations=F)	
		pathModel = 2; simulations = TRUE; outputNames = paste0(analyses[i],"_circuitscape_simulations")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
		 nberOfRandomisations,randomProcedure,outputName=analyses[i],showingPlots,nberOfCores,OS,simulations=T,randomisations=F)	
		pathModel = 3; simulations = FALSE; outputNames = paste0(analyses[i],"_least-cost_extractions")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
		 nberOfRandomisations,randomProcedure,outputName=analyses[i],showingPlots,nberOfCores,OS,simulations=F,randomisations=F)	
		pathModel = 3; simulations = TRUE; outputNames = paste0(analyses[i],"_circuitscape_simulations")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
		 nberOfRandomisations,randomProcedure,outputName=analyses[i],showingPlots,nberOfCores,OS,simulations=T,randomisations=F)	
	}

envVariableNames = c("Elevation","Land_cover_forests","Land_cover_shrublands","Land_cover_savannas","Land_cover_grasslands",
					 "Land_cover_croplands","Land_cover_urban_areas","Annual_mean_temperature","Annual_precipitation"); c = 0
analyses = c("WNV_gamma_all","WNV_gamma_bf02","WNV_gamma_af02","WNV_gamma_NY99","WNV_gamma_WN02","WNV_gamma_SW03")
for (i in 1:length(analyses))
	{
		extractions_LC = read.table(paste0("Seraphim_analyses/",analyses[i],"_least-cost_extractions_LR_results.txt"), header=T)
		extractions_CS = read.table(paste0("Seraphim_analyses/",analyses[i],"_circuitscape_extractions_LR_results.txt"), header=T)
		simulations_LC = read.table(paste0("Seraphim_analyses/",analyses[i],"_least-cost_simulations_LR_results.txt"), header=T)
		simulations_CS = read.table(paste0("Seraphim_analyses/",analyses[i],"_circuitscape_simulations_LR_results.txt"), header=T)
		results_Q1 = matrix(nrow=length(envVariableNames)*3, ncol=4); colnames(results_Q1) = c("LC-C","LC-R","CS-C","CS-R")
		results_Q2 = matrix(nrow=length(envVariableNames)*3, ncol=4); colnames(results_Q2) = c("LC-C","LC-R","CS-C","CS-R")
		results_BF = matrix(nrow=length(envVariableNames)*3, ncol=4); colnames(results_BF) = c("LC-C","LC-R","CS-C","CS-R")
		for (j in 1:length(envVariableNames))
			{
				for (k in c(10,100,1000))
					{
						if ((j == 1)&(k == 10))
							{
								n = 0; rowNames = rep(NA, length(envVariableNames)*3)
							}
						n = n+1; rowNames[n] = paste0(envVariableNames[j],"_k",k); Qe = list(); Qs = list()
						index1 = which(grepl("delta_R2",colnames(extractions_LC))&grepl(envVariableNames[j],colnames(extractions_LC))
									   &grepl(paste0("k",k,"_C"),colnames(extractions_LC)))
						index2 = which(grepl("delta_R2",colnames(simulations_LC))&grepl(envVariableNames[j],colnames(simulations_LC))
									   &grepl(paste0("k",k,"_C"),colnames(simulations_LC)))
						Qe[[1]] = extractions_LC[,index1]; Qs[[1]] = simulations_LC[,index2]
						index1 = which(grepl("delta_R2",colnames(extractions_LC))&grepl(envVariableNames[j],colnames(extractions_LC))
									   &grepl(paste0("k",k,"_R"),colnames(extractions_LC)))
						index2 = which(grepl("delta_R2",colnames(simulations_LC))&grepl(envVariableNames[j],colnames(simulations_LC))
									   &grepl(paste0("k",k,"_R"),colnames(simulations_LC)))
						Qe[[2]] = extractions_LC[,index1]; Qs[[2]] = simulations_LC[,index2]
						index1 = which(grepl("delta_R2",colnames(extractions_CS))&grepl(envVariableNames[j],colnames(extractions_CS))
									   &grepl(paste0("k",k,"_C"),colnames(extractions_CS)))
						index2 = which(grepl("delta_R2",colnames(simulations_CS))&grepl(envVariableNames[j],colnames(simulations_CS))
									   &grepl(paste0("k",k,"_C"),colnames(simulations_CS)))
						Qe[[3]] = extractions_CS[,index1]; Qs[[3]] = simulations_CS[,index2]
						index1 = which(grepl("delta_R2",colnames(extractions_CS))&grepl(envVariableNames[j],colnames(extractions_CS))
									   &grepl(paste0("k",k,"_R"),colnames(extractions_CS)))
						index2 = which(grepl("delta_R2",colnames(simulations_CS))&grepl(envVariableNames[j],colnames(simulations_CS))
									   &grepl(paste0("k",k,"_R"),colnames(simulations_CS)))
						Qe[[4]] = extractions_CS[,index1]; Qs[[4]] = simulations_CS[,index2]
						for (l in 1:length(Qe))
							{
								results_Q1[n,l] = sum(Qe[[l]]>0); results_Q2[n,l] = mean(Qe[[l]]); c = 0
								for (k in 1:length(Qe[[l]]))
									{
										if (Qs[[l]][k] < Qe[[l]][k]) c = c+1
									}
								p = c/length(Qe[[l]]); results_BF[n,l] = p/(1-p)
							}
					}
			}
		results_Q2 = round(results_Q2,2); results_Q2[which(results_BF<20)] = "-"
		cat(analyses[i],":\n"); results_BF = round(results_BF,1)
		results_BF[which(results_BF<3.0)] = "-"; results_BF[which(results_Q1<0.9)] = "-"
		row.names(results_Q1) = rowNames; row.names(results_Q2) = rowNames; row.names(results_BF) = rowNames
		# cat("\n\t% of Q distributions >0:\n\n"); print(round(results_Q1,0))
		# cat("\n\tAsspciated BF supports:\n\n"); print(results_BF); cat("\n")
		cat("\n\tMean Q values:\n\n"); print(results_Q2); cat("\n")
	}

tabs = list()
tabs[[1]] = read.table("Seraphim_analyses/WNV_gamma_all_least-cost_extractions_LR_results.txt", header=T)
tabs[[2]] = read.table("Seraphim_analyses/WNV_gamma_all_circuitscape_extractions_LR_results.txt", header=T)
tabs[[3]] = read.table("Seraphim_analyses/WNV_gamma_all_least-cost_simulations_16_LR_results.txt", header=T)
tabs[[4]] = read.table("Seraphim_analyses/WNV_gamma_all_circuitscape_simulations_16_LR_results.txt", header=T)
for (i in 1:length(tabs))
	{
		tab = tabs[[i]]; colIndices = 1:dim(tab)[2]
		lastIndex = which(grepl("null_raster.1",colnames(tab)))
		indices1 = which(grepl("water",colnames(tab))); colnames1 = colnames(tab)[indices1]
		indices2 = which(grepl("urban",colnames(tab))); colnames2 = colnames(tab)[indices2]
		for (j in 1:length(indices1))
			{
				colIndices[indices1[j]] = colIndices[indices2[j]]
			}
		tab = tab[,colIndices]; tab = tab[,1:(lastIndex-7)]
		if (i == 1) write.table(tab, "Seraphim_analyses/WNV_gamma_all_least-cost_extractions_LR_results.txt", quote=F, row.names=F)
		if (i == 2) write.table(tab, "Seraphim_analyses/WNV_gamma_all_circuitscape_extractions_LR_results.txt", quote=F, row.names=F)
		if (i == 3) write.table(tab, "Seraphim_analyses/WNV_gamma_all_least-cost_simulations_LR_results.txt", quote=F, row.names=F)
		if (i == 4) write.table(tab, "Seraphim_analyses/WNV_gamma_all_circuitscape_simulations_LR_results.txt", quote=F, row.names=F)
	}

allResults1 = matrix(nrow=length(envVariableNames)*2*2*3, ncol=7); kS = c(10,100,1000); CR = c("C","R"); L = 0
colnames(allResults1) = c("Path model","Environmental factor","k","Regression coefficient","Q statistic","p(Q) > 0","BF")
pathModels = c("Least-cost path model","Circuitscape path model"); extractions = list(); simulations = list()
extractions[[1]] = read.table("Seraphim_analyses/WNV_gamma_all_least-cost_extractions_LR_results.txt", header=T)
extractions[[2]] = read.table("Seraphim_analyses/WNV_gamma_all_circuitscape_extractions_LR_results.txt", header=T)
simulations[[1]] = read.table("Seraphim_analyses/WNV_gamma_all_least-cost_simulations_16_LR_results.txt", header=T)
simulations[[2]] = read.table("Seraphim_analyses/WNV_gamma_all_circuitscape_simulations_16_LR_results.txt", header=T)
for (i in 1:length(pathModels))
	{
		for (j in 1:length(envVariableNames))
			{
				for (k in 1:length(CR))
					{
						for (l in 1:length(kS))
							{
								L = L+1; c = 0
								allResults1[L,1] = pathModels[i]; allResults1[L,3] = kS[l]
								allResults1[L,2] = paste0(envVariableNames[j]," (",CR[k],")")
								index1 = which(grepl("LR_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))
											   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
								index2 = which(grepl("delta_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))
											   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
								index3 = which(grepl("delta_R2",colnames(simulations[[i]]))&grepl(envVariableNames[j],colnames(simulations[[i]]))
											   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(simulations[[i]])))
								R2 = extractions[[i]][,index1]; Qe = extractions[[i]][,index2]; Qs = simulations[[i]][,index3]
								for (m in 1:length(Qe))
									{
										if (Qs[m] < Qe[m]) c = c+1
									}
								p = c/length(Qe); BF = p/(1-p)
								allResults1[L,4] = paste0(round(median(R2),3)," [",round(quantile(R2,0.025),3)," - ",round(quantile(R2,0.975),3),"]")
								allResults1[L,5] = paste0(round(median(Qe),3)," [",round(quantile(Qe,0.025),3)," - ",round(quantile(Qe,0.975),3),"]")
								allResults1[L,6] = sum(Qe>0)/nberOfExtractionFiles
								if (as.numeric(allResults1[L,6]) >= 0.9)
									{
										allResults1[L,7] = round(BF,1)
									}	else	{
										allResults1[L,7] = "-"
									}
							}
					}
			}
	}
write.csv(allResults, "Seraphim_results1_NEW.csv", row.names=F, quote=F)

kS = c(10,100,1000); CR = c("C","R"); colNames = c()
allResults2 = matrix(nrow=length(envVariableNames)*2*2*3, ncol=3+(2*6))
for (i in 1:length(analyses)) colNames = c(colNames, c("Q_mean","BF"))
colnames(allResults2) = c("Path model","Environmental factor","k",colNames)
analyses = c("WNV_gamma_all","WNV_gamma_bf02","WNV_gamma_af02","WNV_gamma_NY99","WNV_gamma_WN02","WNV_gamma_SW03")
pathModels = c("Least-cost path model","Circuitscape path model"); extractions = list(); simulations = list()
for (a in 1:length(analyses))
	{
		suffix = unlist(strsplit(analyses[a],"_"))[length(unlist(strsplit(analyses[a],"_")))]; L = 0; colIndex = 3+((a-1)*2)+1
		extractions[[1]] = read.table(paste0("Seraphim_analyses/WNV_gamma_",suffix,"_least-cost_extractions_LR_results.txt"), header=T)
		extractions[[2]] = read.table(paste0("Seraphim_analyses/WNV_gamma_",suffix,"_circuitscape_extractions_LR_results.txt"), header=T)
		simulations[[1]] = read.table(paste0("Seraphim_analyses/WNV_gamma_",suffix,"_least-cost_simulations_LR_results.txt"), header=T)
		simulations[[2]] = read.table(paste0("Seraphim_analyses/WNV_gamma_",suffix,"_circuitscape_simulations_LR_results.txt"), header=T)
		for (i in 1:length(pathModels))
			{
				for (j in 1:length(envVariableNames))
					{
						for (k in 1:length(CR))
							{
								for (l in 1:length(kS))
									{
										L = L+1; c = 0
										allResults2[L,1] = pathModels[i]; allResults2[L,3] = kS[l]
										allResults2[L,2] = paste0(envVariableNames[j]," (",CR[k],")")
										index1 = which(grepl("LR_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))
													   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
										index2 = which(grepl("delta_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))
													   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
										index3 = which(grepl("delta_R2",colnames(simulations[[i]]))&grepl(envVariableNames[j],colnames(simulations[[i]]))
													   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(simulations[[i]])))
										R2 = extractions[[i]][,index1]; Qe = extractions[[i]][,index2]; Qs = simulations[[i]][,index3]
										for (m in 1:length(Qe))
											{
												if (Qs[m] < Qe[m]) c = c+1
											}
										p = c/length(Qe); BF = p/(1-p)
										allResults2[L,colIndex] = round(mean(Qe),3)
										if (as.numeric(sum(Qe>0)/nberOfExtractionFiles) >= 0.9)
											{
												allResults2[L,colIndex+1] = round(BF,1)
											}	else	{
												allResults2[L,colIndex+1] = "-"
											}
									}
							}
					}
			}
	}
write.csv(allResults2, "Seraphim_results2_NEW.csv", row.names=F, quote=F)

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
dev.copy2pdf(file="WNV_figures_&_SI/Figures/Impact_dispersal_velocity_NEW.pdf")

# 10. Testing the impact of migratory flyways on branch dispersal frequency

	# 10.1. Testing the impact of the four US administrative flyways
		
flyways = readOGR(dsn="Environmental_files/WNV_shapefiles/", layer="Migratory_birds_flyways")
rast1 = raster("Environmental_files/WNV_rasters/Elevation_WNV_08.asc"); rast1[!is.na(rast1[])] = 0
rast2 = mask(crop(rast1, flyways), flyways); showingPlots = FALSE
analyses = c("WNV_gamma_all","WNV_gamma_bf02","WNV_gamma_af02",
			 "WNV_gamma_NY99","WNV_gamma_WN02","WNV_gamma_SW03")
for (i in 1:length(analyses))
	{
		if (analyses[i] == "WNV_gamma_all")
			{
				ratiosOfChangingFlywayEvents = matrix(nrow=nberOfExtractionFiles, ncol=3)
				colnames(ratiosOfChangingFlywayEvents) = c("extractions","simulations","randomisations")
			}	else	{
				ratiosOfChangingFlywayEvents = matrix(nrow=nberOfExtractionFiles, ncol=2)
				colnames(ratiosOfChangingFlywayEvents) = c("extractions","simulations")		
			}
		localTreesDirectory = paste0("Tree_extractions/",analyses[i]); c1 = 0; c2 = 0
		for (j in 1:nberOfExtractionFiles)
			{
				tabs = list(); print(c(i,j))
				tabs[[1]] = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), header=T)
				tabs[[2]] = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",j,".csv"), header=T)
				if (analyses[i] == "WNV_gamma_all")
					{
						tabs[[3]] = read.csv(paste0(localTreesDirectory,"/TreeRandomisation_",j,".csv"), header=T)
					}
				for (k in 1:length(tabs))
					{
						if (showingPlots == TRUE)
							{
								dev.new(width=8, height=5.8); par(oma=c(0,0,0,0), mar=c(0,0,0,0))
								plotRaster(rast1, new=F); plot(flyways, add=T, lwd=0.1)
							}
						indices = which((!is.na(raster::extract(rast2,tabs[[k]][,c("startLon","startLat")])))
										&(!is.na(raster::extract(rast2,tabs[[k]][,c("endLon","endLat")]))))
						tab = tabs[[k]][indices,]; n = 0
						for (l in 1:dim(tab)[1])
							{
								pts = cbind(c(tab[l,"startLon"],tab[l,"endLon"]),c(tab[l,"startLat"],tab[l,"endLat"]))
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
						ratiosOfChangingFlywayEvents[j,k] = n/dim(tab)[1]
					}
			}
		fileName = paste0("Seraphim_analyses/",analyses[i],"_ratios_of_changing_flyway_events.csv")
		write.csv(ratiosOfChangingFlywayEvents, fileName, quote=F, row.names=F)
		for (j in 1:dim(ratiosOfChangingFlywayEvents)[1])
			{
				if (ratiosOfChangingFlywayEvents[j,"extractions"] < ratiosOfChangingFlywayEvents[j,"simulations"]) c1 = c1+1
				if (analyses[i] == "WNV_gamma_all")
					{
						if (ratiosOfChangingFlywayEvents[j,"extractions"] < ratiosOfChangingFlywayEvents[j,"randomisations"]) c2 = c2+1
					}
			}
		p1 = c1/dim(ratiosOfChangingFlywayEvents)[1]; BF1 = p1/(1-p1); cat(analyses[i],": ",round(BF1,3),"\n")
			# BF1 = 0.020 (all), 0.042 (bf02), 0.887 (af02), 0.299 (NY99), 0.099 (WN02), 0.235 (SW03)
		if (analyses[i] == "WNV_gamma_all")
			{
				p2 = c2/dim(ratiosOfChangingFlywayEvents)[1]; BF2 = p2/(1-p2); cat(analyses[i],": ",round(BF2,3),"\n") # BF2 = 0.075 (all)
			}
	}

ratiosOfChangingFlywayEvents = read.csv("Seraphim_analyses/WNV_gamma_all_ratios_of_changing_flyway_events.csv", header=T)
nberOfBranches = dim(read.table("Tree_extractions/WNV_gamma_all/TreeExtractions_1.csv", header=T))[1]
nberOfChangingFlywaysEvents = ratiosOfChangingFlywayEvents*nberOfBranches; cols1 = list(); cols2 = list()
dev.new(width=6, height=2.5); par(mfrow=c(1,1), mgp=c(1,0.35,0), oma=c(0.5,1,1,1), mar=c(2.5,2,0.5,0), lwd=0.2)
cols1[[1]] = rgb(204,0,0,255,maxColorValue=255); cols2[[1]] = rgb(204,0,0,100,maxColorValue=255) # red
cols1[[2]] = rgb(120,120,120,255,maxColorValue=255); cols2[[2]] = rgb(120,120,120,100,maxColorValue=255)
plot(density(nberOfChangingFlywaysEvents[,2]), lwd=0.7, col=cols2[[2]], lty=1, axes=F, ann=F, xlim=c(80,180), ylim=c(0,0.07))
polygon(density(nberOfChangingFlywaysEvents[,2]), col=cols2[[2]], border=NA)
lines(density(nberOfChangingFlywaysEvents[,1]), lwd=0.7, col=cols1[[1]], lty=1)
polygon(density(nberOfChangingFlywaysEvents[,1]), col=cols2[[1]], border=NA)
axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.02,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025, col.tick="gray30", col.axis="gray30", col="gray30")
title(xlab=expression(italic(N)), cex.lab=0.7, mgp=c(1,0,0), col.lab="gray30")
title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")

	# 10.2. Testing the impact of the three flyways estimated by La Sorte et al. (2014) for terrestrial birds
		
		# 10.2.1. Loading the two ".RData" files
		
flyways = c("Eastern_flyway","Central_flyway","Western_flyway")
load("Environmental_files/WNV_rasters/Spring_flyways_LaSorte.RData") # dat.spring
load("Environmental_files/WNV_rasters/Autumn_flyways_LaSorte.RData") # dat.autumn
source("Divers_R_functions/decreaseResolution.r")
template = raster("Environmental_files/WNV_rasters/Elevation_WNV_04.asc")
template = decreaseResolution(template, R=5); template[!is.na(template[])] = 0
for (j in 1:length(flyways))
	{
		rast1 = rasterize(dat.spring[,c("lon","lat")], template, dat.spring[,i+2], fun=mean)
		writeRaster(rast1, paste0("Environmental_files/WNV_rasters/",flyways[i],"_spring_20.asc"), overwrite=T)
		rast2 = rasterize(dat.autumn[,c("lon","lat")], template, dat.autumn[,i+2], fun=mean)
		writeRaster(rast2, paste0("Environmental_files/WNV_rasters/",flyways[i],"_autumn_20.asc"), overwrite=T)
		rast3 = rast1; rast3[] = (rast1[]+rast2[])/2; # plotRaster(rast3)
		writeRaster(rast3, paste0("Environmental_files/WNV_rasters/",flyways[i],"_average_20.asc"), overwrite=T)
	}
rast1 = raster(paste0("Environmental_files/WNV_rasters/",flyways[1],"_average_20.asc"))
rast2 = raster(paste0("Environmental_files/WNV_rasters/",flyways[2],"_average_20.asc"))
rast3 = raster(paste0("Environmental_files/WNV_rasters/",flyways[3],"_average_20.asc"))
rast_tot = rast1; rast_tot[] = rast1[]+rast2[]+rast3[]
rast1[] = rast1[]/rast_tot[]; rast2[] = rast2[]/rast_tot[]; rast3[] = rast3[]/rast_tot[]
writeRaster(rast1, paste0("Environmental_files/WNV_rasters/",flyways[1],"_normalised.asc"), overwrite=T)
writeRaster(rast2, paste0("Environmental_files/WNV_rasters/",flyways[2],"_normalised.asc"), overwrite=T)
writeRaster(rast3, paste0("Environmental_files/WNV_rasters/",flyways[3],"_normalised.asc"), overwrite=T)

		# 10.2.2. Comparing the inferred and simulated numbers of changing flyway events

rast = raster("Environmental_files/WNV_rasters/Western_flyway_normalised.asc")
rast[!is.na(rast[])] = 0; flyways = list()
flyways[[1]] = raster("Environmental_files/WNV_rasters/Eastern_flyway_normalised.asc")
flyways[[2]] = raster("Environmental_files/WNV_rasters/Central_flyway_normalised.asc")
flyways[[3]] = raster("Environmental_files/WNV_rasters/Western_flyway_normalised.asc")
analyses = c("WNV_gamma_all","WNV_gamma_bf02","WNV_gamma_af02",
			 "WNV_gamma_NY99","WNV_gamma_WN02","WNV_gamma_SW03")
for (i in 1:length(analyses))
	{
		localTreesDirectory = paste0("Tree_extractions/",analyses[i]); c1 = 0; c2 = 0
		if (analyses[i] == "WNV_gamma_all")
			{
				differenceBetweenFlywayValues = matrix(nrow=nberOfExtractionFiles, ncol=3)
				colnames(differenceBetweenFlywayValues) = c("extractions","simulations","randomisations")
			}	else	{
				differenceBetweenFlywayValues = matrix(nrow=nberOfExtractionFiles, ncol=2)
				colnames(differenceBetweenFlywayValues) = c("extractions","simulations")		
			}
		for (j in 1:nberOfExtractionFiles)
			{
				tabs = list(); print(c(i,j))
				tabs[[1]] = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), header=T)
				tabs[[2]] = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",j,".csv"), header=T)
				if (analyses[i] == "WNV_gamma_all")
					{
						tabs[[3]] = read.csv(paste0(localTreesDirectory,"/TreeRandomisation_",j,".csv"), header=T)
					}
				for (k in 1:length(tabs))
					{
						indices = which((!is.na(raster::extract(rast,tabs[[k]][,c("startLon","startLat")])))
										&(!is.na(raster::extract(rast,tabs[[k]][,c("endLon","endLat")]))))
						tab = tabs[[k]][indices,]; diffs = c()
						startingValues = matrix(nrow=dim(tab)[1], ncol=length(flyways))
						endingValues = matrix(nrow=dim(tab)[1], ncol=length(flyways))
						for (l in 1:length(flyways))
							{
								startingValues[,l] = raster::extract(flyways[[l]], tab[,c("startLon","startLat")])
								endingValues[,l] = raster::extract(flyways[[l]], tab[,c("endLon","endLat")])
							}	
						for (l in 1:dim(tab)[1])
							{
								startingValue = max(startingValues[l,])
								index = which(startingValues[l,]==startingValue)
								endingValue = endingValues[l,index]
								diffs = c(diffs, endingValue-startingValue)
							}
						differenceBetweenFlywayValues[j,k] = mean(diffs)
					}
			}
		fileName = paste0("Seraphim_analyses/",analyses[i],"_difference_between_flyway_values.csv")
		write.csv(differenceBetweenFlywayValues, fileName, quote=F, row.names=F)
		for (j in 1:dim(differenceBetweenFlywayValues)[1])
			{
				if (differenceBetweenFlywayValues[j,"extractions"] > differenceBetweenFlywayValues[j,"simulations"]) c1 = c1+1
				if (analyses[i] == "WNV_gamma_all")
					{
						if (differenceBetweenFlywayValues[j,"extractions"] > differenceBetweenFlywayValues[j,"randomisations"]) c2 = c2+1
					}
			}
		p1 = c1/dim(differenceBetweenFlywayValues)[1]; BF1 = p1/(1-p1); print(BF1)
			# BF1 = 0.000 (all), 0.333 (bf02), 0.053 (af02), 1.564 (NY99), 0.000 (WN02), 0.563 (SW03)
		if (analyses[i] == "WNV_gamma_all")
			{
				p2 = c2/dim(differenceBetweenFlywayValues)[1]; BF2 = p2/(1-p2); print(BF2) # BF2 = 0.010 (all)
			}
	}

		# 10.2.3. Plotting the six different rasters

rS = list(); cols = list()
rS[[1]] = raster("Environmental_files/WNV_rasters/Western_flyway_spring_20.asc"); cols[[1]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
rS[[2]] = raster("Environmental_files/WNV_rasters/Central_flyway_spring_20.asc"); cols[[2]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(100)
rS[[3]] = raster("Environmental_files/WNV_rasters/Eastern_flyway_spring_20.asc"); cols[[3]] = colorRampPalette(brewer.pal(9,"YlGn"))(100)
rS[[4]] = raster("Environmental_files/WNV_rasters/Western_flyway_autumn_20.asc"); cols[[4]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
rS[[5]] = raster("Environmental_files/WNV_rasters/Central_flyway_autumn_20.asc"); cols[[5]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(100)
rS[[6]] = raster("Environmental_files/WNV_rasters/Eastern_flyway_autumn_20.asc"); cols[[6]] = colorRampPalette(brewer.pal(9,"YlGn"))(100)
rS[[7]] = raster("Environmental_files/WNV_rasters/Western_flyway_average_20.asc"); cols[[7]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
rS[[8]] = raster("Environmental_files/WNV_rasters/Central_flyway_average_20.asc"); cols[[8]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(100)
rS[[9]] = raster("Environmental_files/WNV_rasters/Eastern_flyway_average_20.asc"); cols[[9]] = colorRampPalette(brewer.pal(9,"YlGn"))(100)
rasterNames1 = c("Western flyway","Central flyway","Eastern flyway","Western flyway","Central flyway","Eastern flyway","Western flyway","Central flyway","Eastern flyway")
rasterNames2 = c("(spring)","(spring)","(spring)","(autumn)","(autumn)","(autumn)","(average)","(average)","(average)")
r = rS[[1]]; r[!is.na(r)] = 0; p = rasterToPolygons(r, dissolve=T)
dev.new(width=7.5, height=5.6); par(mfrow=c(3,3), oma=c(2,2.5,1,0.3), mar=c(0,0,0,0), mgp=c(1,0.2,0), lwd=0.2)
for (i in 1:length(rS))
	{
		plot(rS[[i]], bty="n", box=F, axes=F, legend=F, col=cols[[i]], colNA="white") # colNA="#D0D0D0")
		plot(rS[[i]], legend.only=T, add=T, col=cols[[i]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.78,0.795,0.18,0.54), adj=3,
		     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.8, col.tick="gray30", col.axis="gray30", col="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)
		lines(p, col="gray30", lwd=0.2)
		if (nchar(rasterNames1[i]) > 0) mtext(rasterNames1[i], side=1, adj=0, line=-4.0, at=-126, cex=0.55, font=1, col="gray30")
		if (nchar(rasterNames2[i]) > 0) mtext(rasterNames2[i], side=1, adj=0, line=-3.2, at=-122, cex=0.55, font=1, col="gray30")
	}
dev.copy2pdf(file="WNV_figures_&_SI/SI_files/Flyways_LaSorte_2014_NEW1.pdf")

rS = list(); cols = list()
rS[[1]] = raster("Environmental_files/WNV_rasters/Western_flyway_average_20.asc"); cols[[1]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
rS[[2]] = raster("Environmental_files/WNV_rasters/Central_flyway_average_20.asc"); cols[[2]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(100)
rS[[3]] = raster("Environmental_files/WNV_rasters/Eastern_flyway_average_20.asc"); cols[[3]] = colorRampPalette(brewer.pal(9,"YlGn"))(100)
rS[[4]] = raster("Environmental_files/WNV_rasters/Western_flyway_normalised.asc"); cols[[4]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
rS[[5]] = raster("Environmental_files/WNV_rasters/Central_flyway_normalised.asc"); cols[[5]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(100)
rS[[6]] = raster("Environmental_files/WNV_rasters/Eastern_flyway_normalised.asc"); cols[[6]] = colorRampPalette(brewer.pal(9,"YlGn"))(100)
rasterNames1 = c("Western flyway","Central flyway","Eastern flyway","Western flyway","Central flyway","Eastern flyway","Western flyway","Central flyway","Eastern flyway")
rasterNames2 = c("(average)","(average)","(average)","(normalised)","(normalised)","(normalised)")
r = rS[[1]]; r[!is.na(r)] = 0; p = rasterToPolygons(r, dissolve=T)
dev.new(width=7.5, height=5.6); par(mfrow=c(3,3), oma=c(2,2.5,1,0.3), mar=c(0,0,0,0), mgp=c(1,0.2,0), lwd=0.2)
for (i in 1:length(rS))
	{
		plot(rS[[i]], bty="n", box=F, axes=F, legend=F, col=cols[[i]], colNA="white") # colNA="#D0D0D0")
		plot(rS[[i]], legend.only=T, add=T, col=cols[[i]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.78,0.795,0.18,0.54), adj=3,
		     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.8, col.tick="gray30", col.axis="gray30", col="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)
		lines(p, col="gray30", lwd=0.2)
		if (nchar(rasterNames1[i]) > 0) mtext(rasterNames1[i], side=1, adj=0, line=-4.0, at=-126, cex=0.55, font=1, col="gray30")
		if (nchar(rasterNames2[i]) > 0) mtext(rasterNames2[i], side=1, adj=0, line=-3.2, at=-122, cex=0.55, font=1, col="gray30")
	}
dev.copy2pdf(file="WNV_figures_&_SI/SI_files/Flyways_LaSorte_2014_NEW2.pdf")

