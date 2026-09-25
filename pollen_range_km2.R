pollen_range_km2 = function(data_subset,taxon_name="Dummy",npts=100,individual_range_km2=4622)
{ ############## NB! Dummy value for ind range. Fix for trees, see triton_range_km2.R
  
  X <- data_subset
  print(dim(X))
  
  e = order(X$time); X = X[e,];
  
  all_times = unique(X$time); 
  #all_times = all_times[all_times-min(all_times)>1]; 
  all_times = all_times[!is.na(all_times)]; 
  
  sample_times = numeric(0); 
  target_times = seq(min(all_times),max(all_times)-1,length=npts);
  for(j in 1:npts) {
    e = which(abs(all_times-target_times[j])==min(abs(all_times-target_times[j]))); 
    sample_times = c(sample_times,all_times[e])
  }
  
  print(sample_times)
  sample_times = unique(sample_times);
  
  nt = length(sample_times); areas = mean_rel_abun = numeric(nt); 
  
  #if(!exists("species_name")) {species_name="kernel"} 
  #pdf(file=paste0(species_name,"_range_maps.pdf"));
  pdf(file=paste0(taxon_name,"_range_maps.pdf"));
  par(mfrow=c(3,3)); 
  for(j in 1:nt) {
    #Xj = X|>filter(time<=sample_times[j]+1)|>filter(time>=sample_times[j]-1); 
    Xj = X|>filter(time<=sample_times[j]+10000)|>filter(time>=sample_times[j]-10000); 
    mean_rel_abun[j] = mean(Xj$rel.abun,na.rm=TRUE); 
    print(Xj)
    
    #  df = data.frame(longitude=Xj$pal.long,latitude=Xj$pal.lat); 
    df = data.frame(longitude=Xj$long,latitude=Xj$lat); 
    # df$animal_id = species_name; 
    df$animal_id = taxon_name; 
    
    print(length(df$longitude))
    print(length(df$latitude))
    print(summary(df$longitude))
    print(summary(df$latitude))
    
    ## If all observations are at the same location, 
    ## the range area is zero, DONE.  
    vlon = diff(range(na.omit(df$longitude))) 
    vlat = diff(range(na.omit(df$latitude))) 
    if(vlon*vlat==0){ 
      areas[j] == individual_range_km2
      print("all observations are at the same location")
    }
    else{
      
      ## Check if there are only two observation locations. 
      n_long = length(unique(df$longitude))
      n_lat = length(unique(df$latitude)); 
      if( abs( (n_long-2)*(n_lat-2) )==0) {
        areas[j]=2*individual_range_km2
      }
      else{	
        
        ## More than two locations, find the kernel-estimated range 
        ## Rename columns to what ctmm requires
        df_formatted <- data.frame(
          `individual-local-identifier` = df$animal_id,
          longitude = df$longitude,
          latitude  = df$latitude,
          timestamp = seq(from = as.POSIXct("2026-01-01 00:00:00", tz = "UTC"), 
                          length.out = nrow(df), 
                          by = "DSTday"),
          check.names = FALSE)
        
        ## Convert to a ctmm telemetry object
        telemetry_data <- as.telemetry(df_formatted)
        
        ##	Independent and Identical Distribution model 
        ##  (suitable for opportunistic observation times)
        iid_model <- ctmm.fit(telemetry_data, ctmm(model = "IID"))
        
        ## Calculate the spherical Kernel Density Estimate
        ## Be a bit conservative (90%) to reduce spillover
        ## of the kernel onto land masses 
        kde_global <- akde(telemetry_data, iid_model)
        # units=FALSE in the summary() gives area in m^2. 
        out = summary(kde_global, level.UD = 0.90, units = FALSE)
        areas[j] = out$CI[2]/(10^6); 
        
        
        plot(telemetry_data, UD = kde_global,main=round(sample_times[j],digits=2), sub=round(log10(areas[j]),digits=1), level.UD = 0.90)
        
      }}
    
    cat(j,nt,sample_times[j],areas[j],mean_rel_abun[j],"\n"); 
  } 
  dev.off(); 
  return(list(sample_times=sample_times,area_km2 = areas,mean_rel_abun=mean_rel_abun)); 
  
}
