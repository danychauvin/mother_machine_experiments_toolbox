# IMPORT DEEP MOMA DATA INTO R DATAFRAME
# Author: Dany Chauvin

# Set defaults variables
dl <- 0.065 # µm/pixel
vertical_cutoff <- 150 #Necessary to get rid of cells that are touching the top of the growth lane where segmentation is not good anymore.

# Set up a parallel environments for multidyplr

if(!exists("mycluster")){
mycluster <- min(30, parallel::detectCores()-1) %>%  # do not use more than 30 cores
  default_cluster() %>% cluster_library( # load currently loaded packages on each core
    names(sessionInfo()$otherPkgs))
  
cluster_copy(mycluster,c("fit_exp_elongation_predict","fit_exp_elongation_slope","fit_exp_elongation_sd_slope","fit_exp_elongation_intercept"))
}

# Importing conditions from dataList

myconditions <- readr::read_csv(path_to_MM_data_summary,
                               col_types = cols(
                                 date=col_character(),
                                 description=col_character(),
                                 f_start=col_character(),
                                 f_end=col_character(),
                                 condition=col_character(),
                                 t_interval=col_double(),
                                 data_path=col_character(),
                                 folder=col_character(),
                                 promoter=col_character(),
                                 vector=col_character()))

# Generate the list of curated data files to fetch
# Computing t_start and t_end in seconds
myconditions <- myconditions %>% 
  mutate(f_start=as.double(f_start),
         f_end=as.double(f_end)) %>% 
  group_by(date,description,data_path,folder) %>% 
  arrange(f_start) %>%
  mutate(t_end=ifelse(row_number()==1,(f_end)*t_interval*60,NA)) %>% 
  mutate(t_start=ifelse(row_number()==1,(f_start)*t_interval*60,NA)) %>% 
  do((function(.df){
    new_df <- .df
    if(nrow(new_df)>=2){
      for(i in c(2:nrow(new_df))) 
      {
        new_df$t_start[i]=new_df$t_end[i-1]
        new_df$t_end[i]=new_df$t_start[i]+(new_df$f_end[i]-new_df$f_start[i])*new_df$t_interval[i]*60}
    }
    return(new_df)
  })(.)) %>%
  ungroup()

#Find the files
myconditions <- myconditions %>% 
  group_by(date,description,data_path,folder) %>% 
  do((function(.df){
    data_folder <- paste(unique(.df$data_path),unique(.df$folder),sep="/")
    new_df <- find.files(data_folder, "ExportedCellStats_*.csv") %>% 
      data.frame(file_path=., stringsAsFactors=FALSE)
    new_df <- crossing(.df,new_df)
    return(new_df)
  } )(.)) %>% 
  ungroup()

# Fetching information about positions to propagate them in the final data table
mypositions <- myconditions %>% 
  ungroup() %>% 
  distinct(data_path) %>% 
  .$data_path %>% 
  lapply(function(.l) as_tibble(data.table::fread(paste(.l,"/positions.csv",sep=""),skip=0,sep=","))) %>% 
  do.call(rbind, .) %>% 
  mutate_at(vars(date), factor) %>% 
  mutate_at(vars(pos), factor)
  

#Define proper column names
nFiles <- myconditions %>% distinct(file_path) %>% nrow()
properColNames <- c("lane_ID","cell_ID","frame","cell_rank","genealogy","type_of_end","parent_ID","cells_in_lane","bbox_top_px","bbox_bottom_px","touches_detection_roi_top","center_x_px","center_y_px","width_px","length_px","tilt_rad","area_px","bgmask_area_px","phc_total_intensity_au","phc_intensity_coefficient_of_variation","label:dead","label:dying","label:fading","label:shrinking","label:dividing","fluo_cellmask_ch_1","fluo_bgmask_ch_1","fluo_ampl_ch_1","fluo_bg_ch_1","fluo_cellmask_ch_2","fluo_bgmask_ch_2","fluo_ampl_ch_2","fluo_bg_ch_2","oriented_bbox_center_x_px","oriented_bbox_center_y_px","oriented_bbox_width_px","oriented_bbox_length_px","oriented_bbox_orientation_angle_rad","contour_area__px","contour_centroid_x__px","contour_centroid_y__px","contour_variance_x__px2","contour_variance_y__px2","contour_covariance__px2")
properColNames <- c(properColNames,"file_path")

#Reading curated data table content
myframes <- myconditions %>% 
  ungroup() %>% 
  distinct(file_path) %>% 
  .$file_path %>% 
  lapply(function(.l) as_tibble(data.table::fread(.l,skip=0,sep=",")) %>% 
  mutate(file_path=.l) %>% 
  mutate(type_of_end=ifelse(type_of_end=="",NA,type_of_end))) %>% 
  lapply(function(.l) .l %>% select(all_of(properColNames))) %>% 
           do.call(rbind, .)

#Parsing positions and growth lanes number
myframes <- myframes %>%
  left_join(myconditions,by=c('file_path')) %>% 
  group_by(file_path,condition) %>% 
  filter(data.table::between(frame,f_start,f_end-1)) %>% 
  tidyr::extract(file_path, c("date","pos", "gl"), ".*(\\d{8})_.*[Pp]os(\\d+).*_GL(\\d+).*", remove=FALSE, convert=TRUE) %>%
  mutate(pos=as.integer(gl)) %>% 
  mutate(pos=as.integer(pos)) %>% 
  mutate(pos=as.integer(date)) %>% 
  ungroup()

#Renaming columns and computing useful variables
myframes <- myframes %>% 
  rename(id=cell_ID,
         end_type=type_of_end,
         parent_id=parent_ID,
         vertical_bottom='bbox_bottom_px',
         vertical_top='bbox_top_px',
         center_x_pixel='center_x_px',
         center_y_pixel='center_y_px',
         width_pixel='width_px',
         length_pixel='length_px',
         tilt_radian='tilt_rad',
         area_pixel2='area_px',
         background_mask_aread_pixel2='bgmask_area_px',
         fluo_cell_mask_ch1='fluo_cellmask_ch_1',
         fluo_background_mask_ch1='fluo_bgmask_ch_1',
         fluo_amplitude='fluo_ampl_ch_1',
         fluo_background_ch1='fluo_bg_ch_1') %>% 
  #Convert date to factors
  mutate_at(vars(date), factor) %>% 
  mutate_at(vars(pos), factor) %>% 
  mutate_at(vars(gl), factor) %>% 
  #Set cell ref
  mutate(cell=paste(date,pos,gl,id,sep='.'),
         parent=paste(date,pos,gl,parent_id,sep='.'),
         lane_ID=paste(date,pos,gl,sep='.')) %>% 
  #Compute biomass estimations
  mutate(length_um_oriented_bbox=oriented_bbox_length_px*dl,
         width_um_oriented_bbox=oriented_bbox_width_px*dl) %>%
  #Compute vertical center
  mutate(vertical_center=(vertical_bottom + vertical_top)/2) %>%
  # Propagating end_type information
  mutate(gl_id=paste(date, pos, gl, sep='.')) %>% 
  group_by(cell) %>% 
  fill(end_type,.direction="up") %>%
  mutate(discard_top=(vertical_top < vertical_cutoff)) %>%
  mutate(end_type=ifelse(any(discard_top), 'touchtop', end_type)) %>%
  ungroup()

#Computing single-cell growth-rate based on oriented bounding box length, also define cell_rank_beg, the rank of the cell at birth
myframes <- myframes %>% 
  #filter(!(cell %in% wrong_cells)) %>%
  mutate(time_sec=t_start+(frame-f_start)*t_interval*60) %>% 
  group_by(cell) %>% 
  partition(cluster=mycluster) %>%
  mutate(cell_rank_beg=first(cell_rank)) %>% 
  mutate(start_time=first(time_sec), end_time=last(time_sec),
         length_predict_oriented_bbox=fit_exp_elongation_predict(time_sec, length_um_oriented_bbox),
         alpha_oriented_bbox=fit_exp_elongation_slope(time_sec,length_um_oriented_bbox),
         sd_alpha_oriented_bbox=fit_exp_elongation_sd_slope(time_sec,length_um_oriented_bbox),
         ) %>%
  collect() %>% 
  arrange(cell, frame) %>% # sort data after `partition()`
  ungroup()

#Add the mean_position of the cell between two divisions in the growth lane.
myframes <- myframes %>% 
  group_by(cell) %>%
  mutate(mean_position_um=mean(vertical_center)*dl) %>% 
  ungroup()

# Add column switch to cells that are experiencing a switch
## LARGE NUMBER OF CELLS LOST HERE
myframes <- myframes %>%
  group_by(cell) %>% 
  arrange(frame) %>% 
  mutate(switch=ifelse(any((time_sec-t_start)<2*3600),TRUE,FALSE)) %>% 
  ungroup()

# An observation "i" (i.e a row), belonging to a cell cycle "c" is in myframes_to_mycells if:
# - myframes_to_mycells contain all other observations about "c", from birth to division.
# - "c" 

myframes_to_mycells <- myframes %>% 
  group_by(cell) %>% 
  filter(parent_id!=-1) %>% #we should see the birth of the cell
  filter(end_type=='div') %>% #we should see cell division
  filter(switch==FALSE) %>% #if an observation occurs during a switch (2h< after the beginning of the experiment, or after a switch), the full cell cycle is discarded.
  ungroup()

mycells <- myframes_to_mycells %>% 
  group_by(date,pos,gl,cell,description,condition,vector) %>% 
  partition(cluster=mycluster) %>%
  do((function(.df) {
    .mod_ll_t <- fastLmPure( cbind(1, .df$time_sec), log(.df$length_um_oriented_bbox) )
    .mod_l_t <- fastLmPure( cbind(1, .df$time_sec), .df$length_um_oriented_bbox)
    .time_birth <- first(.df$time_sec)
    .time_div <- last(.df$time_sec)
    .cell_rank_beg <- first(.df$cell_rank)
    .alpha_oriented_bbox=first(.df$alpha_oriented_bbox)
    .sd_alpha_oriented_bbox=first(.df$sd_alpha_oriented_bbox)
    data.frame(npoints=.mod_ll_t$df.residual+1,
               time_birth=.time_birth, time_div=.time_div, div_time=.time_div-.time_birth,
               l_birth=first(.df$length_predict_oriented_bbox), l_div=last(.df$length_predict_oriented_bbox),
               w_birth=first(.df$width_um_oriented_bbox),w_div=last(.df$width_um_oriented_bbox),
               logl_time_slope=.mod_ll_t$coefficients[2], logl_time_slopesd=.mod_ll_t$stderr[2],
               logl_time_r2=cor(.df$time_sec, log(.df$length_um_oriented_bbox))^2,
               l_time_slope=.mod_l_t$coefficients[2], l_time_slopesd=.mod_l_t$stderr[2], 
               l_time_r2=cor(.df$time_sec, .df$length_um_oriented_bbox)^2,
               cell_rank_beg=.cell_rank_beg,
               alpha_oriented_bbox=.alpha_oriented_bbox,sd_alpha_oriented_bbox=.sd_alpha_oriented_bbox)
  })(.) ) %>% 
  collect() %>% 
  #arrange(condition, date, pos, gl, id) %>% 
  ungroup()

mycells <- mycells %>% 
  # filter by r2 of exponential fit
  #filter(logl_time_r2>0.95) %>% 
  # create new variables of interest
  mutate(dl=l_div - l_birth,
         alpha_len=log(l_div/l_birth) / div_time)
# ok since l_xxx are fitted
#c_birth=g_birth/l_birth,
#c_div=g_div/l_div,
#dg=g_div - g_birth,
#dcdt=(g_div/l_div - g_birth/l_birth) / div_time,
#g=log(g_div/g_birth) / div_time, # would be cleaner from a fit but difficult to compare
#gamma=dg/dl,
#q=dg/dl*alpha)

myframes %>% distinct(cell) %>% nrow()
myframes_to_mycells %>% distinct(cell) %>% nrow()
mycells %>%  distinct(cell) %>% nrow()
myframes_to_mycells <- myframes_to_mycells %>% 
  semi_join(mycells,by=c("cell"))

myframes <- myframes %>%
  left_join(mypositions %>% 
              mutate(pos=as.integer(pos)) %>% 
              select(date,ch,pos,width) %>%
              mutate_at(vars(date), factor) %>% 
              mutate_at(vars(pos), factor) %>% 
              select(date,ch,pos,width),by=c("date","pos"))
  
myframes_to_mycells <- myframes_to_mycells %>% 
  left_join(mypositions %>% 
              mutate(pos=as.integer(pos)) %>% 
              select(date,ch,pos,width) %>%
              mutate_at(vars(date), factor) %>% 
              mutate_at(vars(pos), factor),
            by=c("date","pos"))

