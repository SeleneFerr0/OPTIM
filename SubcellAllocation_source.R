######################################################################
###                                                                ###    
###                Micro Rep Analysis - Subcell Allocation         ###
###                        Main Codes (Part 1)                     ###
###                Modified on 14 Oct 2016 by Xiaoyu              ###
######################################################################

## source R script that contains all the customized functions to be used
source("SubcellAllocation_source_functions.R")

## Exclude cells in RES (sf0) and Panel (pf0) but not in sample design - cells with no Target.SS available
sf0 <- sf0[cell%in%unique(sd[, cell])]
pf0 <- pf0[cell%in%unique(sd[, cell])]

## Specify "mbd" and "cell" as characters
sf0[, ':='(mbd = as.character(mbd), cell = as.character(cell))] 
pf0[, ':='(mbd = as.character(mbd), cell = as.character(cell))]
sd[, ':='(cell = as.character(cell))]
## Remove comma style in numeric variables
sf0[, pf_univ := as.numeric(gsub(",", "", pf_univ))]
pf0[, pf_panel := as.numeric(gsub(",", "", pf_panel))]
sd[, ':='(zuniv = as.numeric(gsub(",", "", zuniv)), avg_acv = as.numeric(gsub(",", "", avg_acv)), 
          ideal_ss = as.numeric(gsub(",", "", ideal_ss)), target_ss = as.numeric(gsub("," , "", target_ss)))]

## Extrat required columns in sample design file and rename the columns
sd <- sd[, list(cell, zuniv, avg_acv, ideal_ss, target_ss)]
setnames(sd, c("Cell", "zUniv", "Avg.ACV", "Ideal.SS", "Target.SS"))

## Set Avg.ACV = 0.01 if Avg.ACV=0
sd[Avg.ACV==0, Avg.ACV := 0.01]

## Extract the list of all secondary characteristic variables
sec_str_char_list <- as.vector(names(sf0)[grepl("sec_var_", names(sf0))])

## Create dummy data base to store the outputs
subcell_all <- NULL      # results at subcell level with default sample rotation
cell_all <- NULL         # results at cell level with default sample rotation
subcell2_all <- NULL     # results at subcell level with user selected sample rotation
cell2_all <- NULL        # results at cell level with user selected sample rotation
kpi_report_all <- NULL   # kpr report
summary_mbd_all <- NULL  # summary at MBD level

## Extract ANOVA score/variable selection score from codebook
anova_score <- codebook[mbd!="total"&!is.na(value)&value!="", list(anova_score = mean(score)), by=list(mbd, variable)]
## Make sure "mbd" column is "Character"
anova_score[, mbd := as.character(mbd)]
setnames(anova_score, "variable", "sec_str_char")

###############################################################################################
# Perform Chi-squre Test and sample allocation over all secondary characteristics - to meet NSC
for (sec_str_char in sec_str_char_list){
  print(sec_str_char)
  sf <- copy(sf0)  # make a copy of sf0
  pf <- copy(pf0)  # make a copy of pf0
  
  ## Create Subcells - cellcode + secondary characteristics + secondary characeristics values
  ### Sample Frame
  sf[, Subcell := paste0(cell, "-", sec_str_char, sf[[sec_str_char]])]
  
  ### Panel Frame
  pf[, Subcell := paste0(cell, "-", sec_str_char, pf[[sec_str_char]])]
  
  ## Create subsets with missing sec char and non-missing sec char
  sfNA <- sf[is.na(sf[[sec_str_char]])|sf[[sec_str_char]]==""]    # missing sec char in sf
  sf <- sf[!outlet%in%sfNA[, outlet]]                             # non-missing sec char in sf
  pfNA <- pf[is.na(pf[[sec_str_char]])|pf[[sec_str_char]]==""]    # missing sec char in pf
  pf <- pf[!outlet%in%pfNA[, outlet]]                             # non-missing sec char in pf
  
  ## Aggregate the total store counts from universe and panel at cell level
  ### Universe: use the projected value (round to integer)
  sf_count <- sf0[, list(Count.SF_cell = round(sum(pf_univ), 0)), by=cell]
  ### panel: use raw store count
  panel_count <- pf0[, list(Current.SS_cell = .N), by=cell]
  cell <- merge(sf_count, panel_count, by = "cell", all=T)
  setnames(cell, "cell", "Cell")
  
  remove(sf_count, panel_count)
  
  ###########################################################
  ## To create subcell dataset 
  ### Create all possible combinations of cell code with secondary store charactistic values
  ### using "sf" to extract sec char levels to avoid NA
  all_comb <- expand.grid(unique(sf0$cell), unique(sf[[sec_str_char]]))
  all_comb$Subcell <- paste0(all_comb[, 1], "-", sec_str_char, all_comb[, 2])
  all_comb <- as.data.table(all_comb)
  
  ### Aggregate store counts at subcell level
  #### stores with non missing secondary characteristic
  ##### Universe: use the projected value (round to integer)
  sf_count <- sf[, list(Count.SF_nomissing = round(sum(pf_univ, na.rm = T), 0)), by=Subcell]
  ##### panel: use raw store count
  panel_count <- pf[,list(Current.SS_nomissing = .N), by=Subcell]
  subcell <- merge(all_comb[, list(Subcell)], sf_count, by="Subcell", all=TRUE)
  subcell <- merge(subcell, panel_count, by="Subcell", all=TRUE)
  
  remove(all_comb, sf_count, panel_count)
  
  ##### Add Cell code to subcell database
  subcell[, Cell := sapply(sapply(Subcell, strsplit, "-"), '[', 1)]
  ##### merge total store counts (including stores with missing sec char value) at cell level to subcell dataset
  subcell <- merge(subcell, cell, by="Cell", all.x=T)
  
  #### Stores with missing secondary characteristic
  ##### Universe: use the projected value (round to integer)
  sf_countNA <- sfNA[, list(Count.SF_missing_cell = round(sum(pf_univ, na.rm = T), 0)), by=cell]
  ##### panel: use raw store count
  panel_countNA <- pfNA[, list(Current.SS_missing_cell = .N), by=cell]
  cellNA <- merge(sf_countNA, panel_countNA, by="cell", all=TRUE)
  setnames(cellNA, "cell", "Cell")
  subcell <- merge(subcell, cellNA, by="Cell", all.x=TRUE)
  
  remove(sf_countNA, panel_countNA, sfNA, pfNA, cellNA)
  
  ### Replace NA in subcell with 0
  subcell[is.na(subcell)] <- 0
  
  ## Aggregate store counts (nonmissing) at cell level
  subcell[, ':='(Count.SF_nomissing_cell = sum(Count.SF_nomissing, na.rm = T)
                 , Current.SS_nomissing_cell = sum(Current.SS_nomissing, na.rm = T))
          , by = list(Cell)]
  
  ## add mbd from sf0 to subcell 
  sf_mbd <- sf0[, list(dummy = .N), by = list(mbd, cell)][, list(cell, mbd)]
  setnames(sf_mbd, "cell", "Cell")
  subcell <- merge(subcell, sf_mbd, by="Cell", all.x=T)
  
  ## extract value of secondary characterisitcs variable
  subcell[, sec_var := sapply(sapply(Subcell, strsplit, sec_str_char), '[', 2)]
  
  ## Aggregate store counts (nonmissing) at mbd X sec_var, and mbd level
  subcell[, ":="(Count.SF_nomissing_mbd_subcell = sum(Count.SF_nomissing)
                 , Current.SS_nomissing_mbd_subcell = sum(Current.SS_nomissing))
          , by = list(mbd, sec_var)][
            , ":="(Count.SF_nomissing_mbd = sum(Count.SF_nomissing)
                   , Current.SS_nomissing_mbd = sum(Count.SF_nomissing))
            , by = list(mbd)]
  
  ## Proportionally allocate stores with missing sec char to subcells
  ### Universe
  subcell[, per_sf := 0][Count.SF_nomissing_cell!=0
                         , per_sf := Count.SF_nomissing/Count.SF_nomissing_cell][
                           , Count.SF_ms_al0 := Count.SF_missing_cell*per_sf]
  #### there may be cases when all stores in a cell in Universe have missing sec char value (this is not expected though),
  #### proportionally allocate stores based on distribution at mbd level in Universe
  subcell[Count.SF_nomissing_cell==0&Count.SF_missing_cell!=0
          , Count.SF_ms_al0 := Count.SF_missing_cell*Count.SF_nomissing_mbd_subcell/Count.SF_nomissing_mbd]
  
  #### the following codes are used to handle the situation when 
  #### sum of proportionally allocated missing store counts of subcells != actual missing store counts at cell level due to rounding
  subcell[, Count.SF_ms_al := round(Count.SF_ms_al0, 0)][
    , Count.SF_ms_al1 := round((Count.SF_ms_al0 - as.integer(Count.SF_ms_al0))*100, 0)]    #2 digits after decimal, the larger the value is, the closer it is to the next integer
  
  setorder(subcell, Cell, -Count.SF_ms_al1, -Count.SF_nomissing)
  
  allocation_to_adjust <- subcell[, list(Count.SF_ms_al_sum = sum(Count.SF_ms_al, na.rm = T)
                                         , Count.SF_missing_cell_mean = mean(Count.SF_missing_cell, na.rm = T))
                                  , by = Cell]
  
  if(sum(allocation_to_adjust[, Count.SF_ms_al_sum] != allocation_to_adjust[, Count.SF_missing_cell_mean])>0){
    
    cell_to_adjust <- allocation_to_adjust[Count.SF_ms_al_sum != Count.SF_missing_cell_mean, Cell]
    setkey(subcell, Cell)
    
    data.table(iter = cell_to_adjust)[
      , {
        while(sum(subcell[iter, Count.SF_ms_al], na.rm = T) != mean(subcell[iter, Count.SF_missing_cell], na.rm = T)){
          
          if(sum(subcell[iter, Count.SF_ms_al], na.rm = T) < mean(subcell[iter, Count.SF_missing_cell], na.rm = T)){
            ## 1 will be added to subcell with larger Count.SF_ms_al1
            ## in case of tie in Count.SF_ms_al1, 1 will be added to subcell with larger Count.SF_nomissing 
            max_subcell <- subcell[iter, Subcell][which.max(subcell[iter, Count.SF_ms_al1])]
            subcell[Subcell%in%max_subcell, Count.SF_ms_al := Count.SF_ms_al + 1]
          }
          
          if(sum(subcell[iter, Count.SF_ms_al], na.rm = T) > mean(subcell[iter, Count.SF_missing_cell], na.rm = T)){
            ## 1 will be extracted from subcell with smaller Count.SF_ms_al1
            ## in case of tie in Count.SF_ms_al1, 1 will be extracted from subcell with smaller Count.SF_nomissing 
            ## subcells with Count.SF_ms_al = 0 will be excluded from this step
            min_subcell <- subcell[iter][order(-Count.SF_ms_al1, Count.SF_nomissing)][Count.SF_ms_al!=0, list(Subcell, Count.SF_ms_al1)]
            min_subcell <- min_subcell[which.min(min_subcell[, Count.SF_ms_al1]), Subcell]
            subcell[Subcell%in%min_subcell, Count.SF_ms_al:= Count.SF_ms_al - 1]
          }
        }  
      }
      , by = iter
      ]
  }
  
  remove(allocation_to_adjust)
  
  #### calculate the total Count.SF which is the sum of no. of stores with non-missing value and allocated no. of stores with missing value 
  subcell[, ':='(Count.SF = Count.SF_nomissing + Count.SF_ms_al, Count.SF_ms_al0 = NULL, Count.SF_ms_al1 = NULL)]
  
  ### Panel
  subcell[, per_panel := 0][Current.SS_nomissing_cell!=0
                            , per_panel := Current.SS_nomissing/Current.SS_nomissing_cell][
                              , Current.SS_ms_al0 := Current.SS_missing_cell*per_panel]
  
  #### there may be cases when all stores in a cell in panel have missing sec char value (this is not expected though),
  #### proportionally allocate stores based on distribution at mbd level in panel
  subcell[Current.SS_nomissing_cell==0&Current.SS_missing_cell!=0
          , Current.SS_ms_al0 := Current.SS_missing_cell*Current.SS_nomissing_mbd_subcell/Current.SS_nomissing_mbd]
  
  #### the loop below is used to handle the situation when 
  #### sum of proportionally allocated missing store counts of subcells != actual missing store counts at cell level due to rounding
  subcell[, Current.SS_ms_al := round(Current.SS_ms_al0, 0)][
    , Current.SS_ms_al1 := round((Current.SS_ms_al0 - as.integer(Current.SS_ms_al0))*100, 0)]   #2 digits after decimal, the larger the value is, the closer it is to the next integer
  
  setorder(subcell, Cell, -Current.SS_ms_al1, -Current.SS_nomissing)
  
  allocation_to_adjust <- subcell[, list(Current.SS_ms_al_sum = sum(Current.SS_ms_al, na.rm = T)
                                         , Current.SS_missing_cell_mean = mean(Current.SS_missing_cell, na.rm = T))
                                  , by = Cell]
  
  if(sum(allocation_to_adjust[, Current.SS_ms_al_sum] != allocation_to_adjust[, Current.SS_missing_cell_mean])>0){
    
    cell_to_adjust <- allocation_to_adjust[Current.SS_ms_al_sum != Current.SS_missing_cell_mean, Cell]
    setkey(subcell, Cell)
    
    data.table(iter = cell_to_adjust)[
      , {
        while(sum(subcell[iter, Current.SS_ms_al], na.rm = T) != mean(subcell[iter, Current.SS_missing_cell], na.rm = T)){
          
          if(sum(subcell[iter, Current.SS_ms_al], na.rm = T) < mean(subcell[iter,Current.SS_missing_cell], na.rm = T)){
            ## 1 will be added to subcell with larger Current.SS_ms_al1
            ## in case of tie in Current.SS_ms_al1, 1 will be added to subcell with larger Current.SS_nomissing
            max_subcell <- subcell[iter, Subcell][which.max(subcell[iter, Current.SS_ms_al1])]
            subcell[Subcell%in%max_subcell, Current.SS_ms_al:= Current.SS_ms_al + 1]
          }
          
          if(sum(subcell[iter, Current.SS_ms_al], na.rm = T)>mean(subcell[iter,Current.SS_missing_cell], na.rm = T)){
            ## 1 will be extracted from subcell with smaller Current.SS_ms_al1
            ## in case of tie in Current.SS_ms_al1, 1 will be extracted from subcell with smaller Current.SS_nomissing
            ## subcells with Current.SS_ms_al = 0 will be excluded from this step
            min_subcell <- subcell[iter][order(-Current.SS_ms_al1, Current.SS_nomissing)][Current.SS_ms_al!=0, list(Subcell, Current.SS_ms_al1)]
            min_subcell <- min_subcell[which.min(min_subcell[, Current.SS_ms_al1]), Subcell]
            subcell[Subcell%in%min_subcell, Current.SS_ms_al := Current.SS_ms_al - 1]
          }
        }    
      }
      , by = iter
      ]
  }
  
  remove(allocation_to_adjust)
  
  #### calculate the total Current.SS which is the sum of no. of stores with non-missing value and allocated no. of stores with missing value 
  subcell[, ':='(Current.SS = Current.SS_nomissing + Current.SS_ms_al
                 , Current.SS_ms_al0 = NULL, Current.SS_ms_al1 = NULL)]
  
  ### Keep a copy of subcell with the orignal missing store counts information as subcell0
  ### to be used in the MBD importance calculation part
  subcell0 <- copy(subcell)
  ### Exclude unnecessary columns from subcell
  subcell <- subcell[, list(Cell, Subcell, mbd, sec_var, Count.SF, Current.SS, Count.SF_cell, Current.SS_cell)]
  
  ###########################################################
  ## to create cell dataset
  ### Extract corresponding cells from sample design file
  cell <- sd[Cell%in%unique(subcell$Cell)]
  
  ### Set Adjusted TargetSS equal to Target SS to start with
  cell[, Adjusted.Target.SS := Target.SS]
  
  ### Merge cell Target.SS to subcell file
  subcell <- merge(subcell, cell[, list(Cell, Target.SS)], by="Cell", all.x=T)
  setnames(subcell, "Target.SS", "Target.SS_cell")
  
  ### calculate the Target SS at subcell level
  #### if Current SS cell != 0, Target SS at subcell level = Target SS Cell * (Current SS Subcell/Current SS Cell)
  #### if Current SS cell = 0, Target SS at subcell level = Target SS Cell * (Universe SS Subcell/Universe SS Cell)
  subcell[, Target.SS := ifelse(Current.SS_cell==0, round(Count.SF/Count.SF_cell*Target.SS_cell, 0)
                                , round(Current.SS/Current.SS_cell*Target.SS_cell, 0))]
  
  ### Target.SS1 is created to handle the rounding issues in the for loop below
  subcell[, Target.SS1 := ifelse(Current.SS_cell==0, round((Count.SF/Count.SF_cell*Target.SS_cell - as.integer(Count.SF/Count.SF_cell*Target.SS_cell))*100, 0)
                                 , round((Current.SS/Current.SS_cell*Target.SS_cell - as.integer(Current.SS/Current.SS_cell*Target.SS_cell))*100, 0))]
  
  ### Alternative option:
  ### Target SS at subcell level = Target SS Cell * (Universe SS Subcell/Universe SS Cell)
  # subcell[, Target.SS := ifelse(Current.SS_cell==0, 0, round(Count.SF/Count.SF_cell*Target.SS_cell, 0))]
  # subcell[, Target.SS1 := ifelse(Current.SS_cell==0, 0, round((Count.SF/Count.SF_cell*Target.SS_cell - as.integer(Count.SF/Count.SF_cell*Target.SS_cell))*100, 0))]
  
  setorder(subcell, Cell, -Target.SS1, -Target.SS)
  
  allocation_to_adjust <- subcell[, list(Target.SS_sum = sum(Target.SS, na.rm = T)
                                         , Target.SS_cell_mean = mean(Target.SS_cell, na.rm = T))
                                  , by = Cell]
  if(sum(allocation_to_adjust[, Target.SS_sum] != allocation_to_adjust[, Target.SS_cell_mean])>0){
    
    cell_to_adjust <- allocation_to_adjust[Target.SS_sum != Target.SS_cell_mean, Cell]
    setkey(subcell, Cell)
    
    data.table(iter = cell_to_adjust)[
      , {
        while(sum(subcell[iter, Target.SS], na.rm = T) != mean(subcell[iter, Target.SS_cell], na.rm = T)){
          
          if(sum(subcell[iter, Target.SS], na.rm = T) < mean(subcell[iter, Target.SS_cell], na.rm = T)){
            ## 1 will be added to subcell with larger Target.SS1
            ## in case of tie in Target.SS1, 1 will be added to subcell with larger Target.SS
            max_subcell <- subcell[iter, Subcell][which.max(subcell[iter, Target.SS1])]
            subcell[Subcell%in%max_subcell,Target.SS:= Target.SS + 1]
          }
          
          if(sum(subcell[iter, Target.SS], na.rm = T) > mean(subcell[iter, Target.SS_cell], na.rm = T)){
            ## 1 will be extracted from subcell with smaller Target.SS1
            ## in case of tie in Target.SS1, 1 will be extracted from subcell with smaller Target.SS
            ## subcells with Target.SS = 0 will be excluded from this step
            min_subcell <- subcell[iter][order(-Target.SS1, Target.SS)][Target.SS!=0, list(Subcell, Target.SS1)]
            min_subcell <- min_subcell[which.min(min_subcell[, Target.SS1]), Subcell]
            subcell[Subcell%in%min_subcell, Target.SS:= Target.SS - 1]
          }
        }
        
      }
      , by = iter
      ]
  }
  
  remove(allocation_to_adjust)
  
  subcell[,':='(Target.SS_cell=NULL, Target.SS1= NULL)]
  
  ### set the minimum retain sample size equal to min. of current SS and Target SS or 0, based on user selection
  if(min_retain_ss=="Current.SS") {
    subcell[, min_Retain_SS := min(Current.SS, Target.SS), by = list(Subcell)]
  }  else {
    subcell[, min_Retain_SS := 0]
  }
  
  ### number of subcell and cell
  subcell_n <- nrow(subcell)
  cell_n <- nrow(cell)
  
  ## compute the Randomized_Count_SF in Subcell 
  subcell[, Randomized_Count_SF := Count.SF + runif(subcell_n, min=0, max=1)/10 - 0.05]
  
  ## compute xUniv, xUniv% and Cum xUniv% at cell level
  cell[, xUniv := zUniv*Avg.ACV] 
  setorder(cell, -xUniv)
  
  cell[, xUniv_per := xUniv/sum(xUniv)]
  setorder(cell, -xUniv_per)
  cell[, Cum_xUniv_per := sprintf("%1.2f%%", 100*cumsum(xUniv_per))]
  
  ## Aggregate subcell level information at cell level
  Count_of_Subcell <- subcell[, list(Count_of_Subcell = .N), by=Cell]
  cell <- merge(cell, Count_of_Subcell, by = "Cell", all=T)
  
  compare_cell <- subcell[, list(min_Retain_SS_cell = sum(min_Retain_SS, na.rm = T)
                                 , Randomized_Count_SF_cell = sum(Randomized_Count_SF, na.rm = T)
                                 , Current.SS_cell = sum(Current.SS, na.rm = T)
                                 , Count.SF_cell = sum(Count.SF, na.rm = T))
                          , by=list(Cell)]
  cell <- merge(cell, compare_cell, by="Cell", all=T)
  
  names_to_change <- c("zUniv", "Avg.ACV", "Ideal.SS", "Target.SS", "Adjusted.Target.SS", "Count_of_Subcell")
  setnames(cell, names_to_change, paste0(names_to_change, "_cell"))
  
  remove(Count_of_Subcell, compare_cell)
  
  #### Retain Check - suggested by Tamas to turn it off, not need the check as we set min Retain SS = min of target ss and current ss
  # cell[, Retain_Check:= ifelse(Target.SS_cell >= min_Retain_SS_cell, 0, 1)]
  # if(sum(cell[, Retain_Check])>0) stop(paste("For Secondary Store Charastistics-",
  #                                         sec_str_char,", Ideal Sample Size is less than Minumum Retain in cell(s)",
  #                                         list(cell[Target.SS_cell<min_Retain_SS_cell,Cell]),
  #                                         ". Please manually change the Current SS / Ideal SS."))
  
  subcell[, ':='(Current.SS_cell = NULL, Count.SF_cell = NULL)]
  subcell <- merge(subcell, cell, by="Cell", all.x=TRUE)
  
  ## Compute subcell allocation at subcell level
  subcell[, ':='(Allocated_SS_prop = ifelse(Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell<Count.SF
                                            , Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell, Count.SF))]
  subcell[, ':='(Allocated_SS_Retian = ifelse(Target.SS_cell==min_Retain_SS_cell, min_Retain_SS
                                              , ifelse(min_Retain_SS>Allocated_SS_prop, min_Retain_SS, Allocated_SS_prop)))]
  subcell[, ':='(Final_Subcell_Allocation = round(Allocated_SS_Retian, 0))]
  
  ## Aggregate subcell allocation at cell level
  compare_cell_2 <- subcell[, list(Final_Subcell_Allocation_cell = sum(Final_Subcell_Allocation, na.rm = T)), by=Cell]
  cell <- merge(cell, compare_cell_2, by="Cell", all.x=TRUE)
  
  remove(compare_cell_2)
  
  ## Perform subcell allocation - the calculation follows the excel template
  ### Put Adjusted Target SS same as Target SS to start with
  cell[, Adjusted.Target.SS_cell := Target.SS_cell]
  allocation <- allocation_function(subcell, cell)  # allocation_function is defined in the source functions script 
  subcell <- allocation[[1]]
  cell <- allocation[[2]]
  
  ## add sec_str_char to subcell and cell
  subcell[, sec_str_char := sec_str_char]
  cell[, sec_str_char := sec_str_char]
  
  ####################################################################
  ## MBD importance score computation 
  
  ### Take zUniv from sample design, recalculate projection factor based on new sample size
  sd_cell_zUniv <- merge(sd[, list(Cell, zUniv)], cell[, list(Cell, Final_Subcell_Allocation_cell)], by="Cell", all.y=T)
  sd_cell_zUniv[, cell_projection_factor := zUniv/Final_Subcell_Allocation_cell]
  
  ce5 <- subcell[, list(Cell, sec_var, mbd, Subcell, Final_Subcell_Allocation)]
  ce5 <- merge(ce5, sd_cell_zUniv, by = "Cell", all.x=T)
  
  ### Merge with original universe and panel store counts
  ce5 <- merge(ce5, subcell0[, list(Subcell, Count.SF_ms_al, Current.SS_ms_al, Current.SS)], by="Subcell", all.x=T)
  ce5[, ':='(sample_subcell_nonmissing = (Final_Subcell_Allocation - Current.SS_ms_al)
             , zuniv_subcell_nonmissing = (Final_Subcell_Allocation - Current.SS_ms_al)*cell_projection_factor)]
  
  ### SUM OF UNIVERSE AT CELL LEVEL
  ce_cell <- ce5[, list(sample_cell_nonmissing = sum(sample_subcell_nonmissing, na.rm = T)
                        , zuniv_cell_nonmissing = sum(zuniv_subcell_nonmissing, na.rm = T))
                 , by=list(Cell)]
  ce5 <- merge(ce5, ce_cell, by="Cell", all.x=)
  ce5[, ':='(sample_cell = Final_Subcell_Allocation_cell, zuniv_cell = zUniv)]
  
  ### SUM OF UNIVERSE AT MBD X SEC Char level
  ce_mbdsec <- ce5[, list(sample_mbd_secvar = sum(sample_subcell_nonmissing, na.rm = T)
                          , zuniv_mbd_secvar = sum(zuniv_subcell_nonmissing, na.rm = T))
                   , by=list(mbd, sec_var)]
  ce5 <- merge(ce5, ce_mbdsec, by=c("mbd", "sec_var"), all.x=T) 
  
  ### SUM OF UNIVERSE AT MBD level
  ce_mbd <- ce5[, list(sample_mbd_nonmissing = sum(sample_subcell_nonmissing, na.rm = T)
                       , zuniv_mbd_nonmissing = sum(zuniv_subcell_nonmissing, na.rm = T)
                       , sample_mbd = sum(Final_Subcell_Allocation, na.rm = T)
                       , zuniv_mbd = sum(Final_Subcell_Allocation*cell_projection_factor, na.rm = T))
                , by=list(mbd)]
  ce5 <- merge(ce5, ce_mbd, by="mbd", all.x=T) 
  
  remove(sd_cell_zUniv, ce_cell, ce_mbd, subcell0)
  
  ### calulate parameters
  ce5[, ':='(Final_Subcell_Allocation=NULL, cell_projection_factor=NULL)][                              # remove unnecessary columns
    , ':='(p0 = ifelse(sample_cell_nonmissing==0, 0, sample_subcell_nonmissing/sample_cell_nonmissing)
           , p = ifelse(zuniv_mbd_nonmissing==0, 0, zuniv_mbd_secvar/zuniv_mbd_nonmissing)
           , w = pmin(sample_cell_nonmissing/5, 1))][
             , ':='(usedp = w*p0+(1-w)*p)][
               , ':='(v = ifelse(sample_cell_nonmissing==0|zuniv_cell==0, 0
                                 , (zuniv_cell*zuniv_cell/sample_cell_nonmissing*usedp*(1-usedp))*(1-sample_cell_nonmissing/zuniv_cell))
                      , charuniv_sample = ifelse(sample_subcell_nonmissing==0, zuniv_cell*usedp, zuniv_cell*p0))]
  
  ### VARIANCE AGGREGATED AT mbd & sec char Value level
  kpi0 <- ce5[, list(charuniv_sample = sum(charuniv_sample, na.rm = T), v = sum(v, na.rm = T))
              , by=list(mbd, sec_var)][
                , sey := sqrt(v)]                #SAMPLING ERROR CALCULATION
  
  ### UNIVERSE
  univ <- sf[, .SD, .SDcols = c(sec_str_char, "mbd", "pf_univ")]
  setnames(univ, sec_str_char, "sec_var")
  univ <- univ[, sec_var := as.character(sec_var)][!(is.na(sec_var)|sec_var=="")
                                                   , list(resuniv = sum(pf_univ, na.rm = T)), by=list(mbd, sec_var)]
  ### to avoid missing sec char values in RES dropped from the calculation 
  univ <- merge(kpi0[, list(mbd,sec_var)], univ, by=c("mbd", "sec_var"), all=T)
  univ[is.na(resuniv), resuniv := 0]
  
  ### CREATE U0 WITH UNIVERSE AGGREGATED AT mbd LEVEL (INCLUDING MISSING CHAR VAL)
  u0 <- sf[, list(univ_full = sum(pf_univ, na.rm = T)), by=mbd]
  ### CREATE U1 WITH UNIVERSE AGGREGATED AT mbd LEVEL (EXCLUDING MISSING CHAR VAL) 
  u1 <- sf[!is.na(sf[[sec_str_char]]), list(univ_full_nonmissing = sum(pf_univ, na.rm = T)), by=mbd]
  
  ### RATIO BY mbd LEVEL (DESIRED TO BE CLOSE TO 1); WARNING TO BE DISPLAYED IF RATIO > 1.2
  univcorr <- merge(u0, u1, by="mbd", all=T)
  univcorr[, ratio := univ_full/univ_full_nonmissing]
  
  ### UNIVERSE ADJUSTED ACCOUNTING FOR MISSING CHAR VAL
  univ <- merge(univ, univcorr[, list(mbd, ratio)], by="mbd", all.x=T)
  univ[, resuniv_uplifted := resuniv*ratio]
  kpi0 <- merge(kpi0, univ[, list(mbd, sec_var, resuniv_uplifted, ratio)], by=c("mbd", "sec_var"), all=T)
  
  ### Sample Completeness at mbd level
  #### take the zuniv_mbd figure from the smple design
  chlevel0 <- merge(sd[, list(Cell, zUniv)], ce5[, list(Cell, mbd, sec_var, Current.SS, Current.SS_ms_al)], by="Cell", all.y=T)
  
  chlevel1 <- merge(chlevel0
                    , chlevel0[, list(Current.SS_cell = sum(Current.SS, na.rm = T), Current.SS_cell_nonmissing = sum(Current.SS - Current.SS_ms_al)), by=list(Cell)]
                    , by="Cell", all.x=T)
  chlevel1[, projection_facor_current := ifelse(Current.SS_cell==0, 0, zUniv/Current.SS_cell)]
  
  chlevel <- chlevel1[, list(zuniv_mbd = sum(Current.SS*projection_facor_current, na.rm = T)
                             , sample_mbd = sum(Current.SS, na.rm = T)
                             , zuniv_mbd_nonmissing = sum((Current.SS - Current.SS_ms_al)*projection_facor_current, na.rm = T)
                             , sample_mbd_nonmissing = sum(Current.SS - Current.SS_ms_al, na.rm = T))
                      , by = mbd]
  chlevel[, ':='(completeness1 = sample_mbd_nonmissing/sample_mbd
                 , completeness2 = zuniv_mbd_nonmissing/zuniv_mbd)]
  
  ### COMBINING RES AND PANEL SUMMARIES AT mbd LEVEL
  comb <- merge(u0, chlevel[, sample_mbd := NULL], by = "mbd")
  
  comb[, ratio2 := zuniv_mbd/univ_full]
  
  kpi1 <- merge(kpi0, ce_mbdsec[ , list(mbd, sec_var, sample_mbd_secvar)], by=c("mbd", "sec_var"), all.x=T)
  
  setnames(kpi1, c("charuniv_sample", "resuniv_uplifted", "sample_mbd_secvar")
           , c("panel_estimate", "universe0", "panel_size"))
  
  kpi2 <- merge(kpi1, chlevel[, list(mbd,zuniv_mbd, completeness1, completeness2)], by = "mbd", all.x=T)
  kpi2 <- kpi2[is.na(panel_estimate), ':='(panel_estimate=0, panel_size=0)]
  
  kpi1c <- merge(kpi2, comb[, list(mbd, ratio2)], by="mbd")
  
  kpi1c[, ':='(universe = universe0*ratio2)][
    , ':='(diff = panel_estimate-universe)][
      , ':='(z = diff/sey)][
        , ':='(comprat = ifelse(universe==0, 0, panel_estimate/universe), rse = ifelse(sey==0, 0, sey/panel_estimate), univprop=universe/zuniv_mbd)]
  kpi1c[v==0, z := ifelse(univprop>0.05, -9, 0.1)]
  
  kpi1c[ , ':='(t = as.integer(z))][
    , ':='(tt = abs(t))][
      , ':='(invrat = 1/ratio)][
        , ':='(limit = 1-apply(cbind(invrat, completeness2), 1, min))]
  
  remove(ce5, univ, univcorr, u0, u1, chlevel0, chlevel1, chlevel,comb, kpi0, kpi1, kpi2, ce_mbdsec)
  
  ### Concerns with the data
  ####   A will refer to our uncertainty, whether the incomplete RES &panel observations for the specific mbd*secchar combination 
  ####   could be the source of the high structural observed difference; B indicates that panel RES was incomplete, both in number 
  ####   of stores and numerical importance of the missing observations within the mbd; C is less serious condition, as the missing 
  ####   panel observations have less than 5% numerical weight in the respective universe for the secondary characreristic;
  ####   D is an indication that universe instruction process is broken; we shall have less confidence in the outcome;
  ####   E indicates an even bigger discrepancy between the universe figures, trigger for input correction;
  ####   F indicates that the calculation outcome has higher than 10% sampling error, i.e. the conclusions are based on relatively low sample sizes;
  ####   G indicates that the calculation outcome has higher than 15% sampling error, i.e. the action should get low priority;
  ####   H indicates that the given characteristic has relatively low importance within the reported universe (less than 5%),therefore the action should get low priority;
  ####   I indicates that the importance of the group is still relatively low, the "punishment" will be slightly less, but still indicating the lower importance;
  ####   J indicates that the relative difference betweeen the "true" and represented mbd*secondary characreristic value universe is beyond 20%, 
  ####   therefore the action should get higher priority, if all other scores are the same.
  
  kpi1c[, ':='(A = ifelse(limit>=abs(diff)/zuniv_mbd, "A", " "))][
    , ':='(BC = ifelse(is.na(A)&completeness1<.95&completeness2<0.95
                       , "B", ifelse(is.na(A)&completeness1<.95&completeness2>=0.95, "C", " ")))][
                         , ':='(DE = ifelse(ratio2<0.8|ratio2>1.2, "E", ifelse(ratio2<0.9|ratio2>1.1, "D", " ")))][
                           , FG := ifelse(rse>0.15, "G", ifelse(rse>0.1, "F", " "))][
                             , ':='(HI = ifelse(univprop<=0.05, "H", ifelse(univprop<=0.10, "I", " ")))][
                               , ':='(J = ifelse(comprat>1.2|comprat<0.8, "J", " "))][
                                 , ':='(concern = paste0(A,BC, DE, FG, HI, J))]
  
  ### Compute MBD importance score based on the concerns assigned
  kpi1c[, ':='(score = ifelse(tt<2, 0, ifelse(tt>=3, 100, 75)))][
    , ':='(score = ifelse(A=="A", score-100, ifelse(BC=="B", score-50,ifelse(BC=="C", score-5, score))))][
      , ':='(score = ifelse(DE=="D", score-20, ifelse(DE=="E", score-100, score)))][
        , ':='(score = ifelse(FG=="F", score-30, ifelse(FG=="G", score-74, score)))][
          , ':='(score = ifelse(HI=="H", score-75, ifelse(HI=="I", score-50, score)))][
            , ':='(score = ifelse(concern=="     J", score+20, score))]
  
  ### if sey = 0, no concern should be rised, score2=0
  kpi1c[, score2 := 0][sey!=0, ':='(score2 = ifelse(tt<2, 0, ifelse(tt>=3, 70+10*tt, 75+(tt-2)*25)))]  
  kpi1c[sey!=0,':='(score2 = ifelse(A=="A", score2-80-pmin(40,(limit*zuniv_mbd - diff)/sey*40)
                                    , ifelse(BC=="B", score2-pmin(0.1, (1-completeness1))*400 - pmin(0.1, (1-completeness2))*400-10
                                             , ifelse(BC=="C", score2-5-(pmin((1-completeness1),0.5)*100-5)/3, score2))))]
  kpi1c[sey!=0, ':='(score2 = ifelse(DE=="D", score2-20-8*((pmin(abs(ratio2-1),0.2))*100-10)
                                     , ifelse(DE=="E", score2-100, score2)))]
  kpi1c[sey!=0, ':='(score2 = ifelse(FG=="F", score2-30-(rse-0.1)*880, ifelse(FG=="G", score2-74, score2)))]
  kpi1c[sey!=0, ':='(score2 = ifelse(HI=="H", score2-85+univprop*500, ifelse(HI=="I", score2-85+univprop*500, score2)))]
  kpi1c[sey!=0, ':='(score2 = ifelse(concern=="     J", score2+pmin(abs((comprat-1)*100), 40), score2))]
  
  kpi_report <- kpi1c[, list(mbd, sec_var, panel_size, panel_estimate, universe, z
                             , completeness1, completeness2, ratio, ratio2,concern, score, score2)]
  
  ### Get the max score within MBD
  max_score_mbd <- kpi1c[, list(max_score = max(score), max_score2 = max(score2)), by=list(mbd)]
  kpi_report <- merge(kpi_report, max_score_mbd, by="mbd", all.x=T)
  kpi_report[, ':='(sec_str_char = sec_str_char, key = paste(mbd, sec_var))]
  
  ### Merge MBD to subcell file
  subcell <- merge(subcell, max_score_mbd, by="mbd", all.x=T)
  
  setkey(anova_score, sec_str_char)
  subcell <- merge(subcell, anova_score[sec_str_char], by=c("mbd", "sec_str_char"))
  
  ### Sample rotation is recommended if MBD importance score >0, anove score > 1.5 and p-value of Final Chi-sqr test < pvalue
  subcell[, sample_rotation := ifelse((max_score>0)&anova_score>1.5&Final_Chi_Test<pvalue, 1, 0)]
  
  ### Merge mbd names to cell
  cell <- merge(cell, sf_mbd, by="Cell", all.x=T)
  
  #############################################################
  ###### sample Rotation with user's specification
  ### create a copy to cell and subcell to store the results from sample rotation
  cell2 <- copy(cell)
  subcell2 <- copy(subcell)
  
  if(sum(subcell$sample_rotation)>0&is.null(mbd_rotation)){
    ### Automatically perform sample rotation if subcell sample rotation is flagged 
    rotation <- rotation_function(subcell2, cell2, pvalue, mbd_rotation)
    subcell2 <- rotation[[1]]
    cell2 <- rotation[[2]]
  }else{
    subcell2 <- NULL
    cell2 <- NULL
  }
  
  ### Perform subcell allocation with sample rotation  
  if(!is.null(subcell2)&!is.null(cell2)&sum(subcell2$min_Retain_SS_org0>subcell2$min_Retain_SS)>0){
    
    allocation2 <- allocation_function(subcell2, cell2)
    subcell2 <- allocation2[[1]]
    cell2 <- allocation2[[2]]
  }
  
  ### summary at MBD level
  summary_mbd <- subcell[, list(max_score = mean(max_score, na.rm = T)
                                , max_score2 = mean(max_score2, na.rm = T)
                                , min_Final_Chi_Test = min(Final_Chi_Test, na.rm = T))
                         , by=list(mbd)]
  
  ### Consolidate concerns at MBD X Sec Char level to MBD level
  summary_mbd[, concern_mbd0 := ""]
  
  for(i in unique(kpi_report$mbd)){
    
    for(j in unique(kpi_report$sec_var)){
      summary_mbd[mbd==i, concern_mbd0 := paste0(concern_mbd0, kpi_report[mbd==i&sec_var==j, concern])]
    }
    summary_mbd[mbd==i, concern_mbd :=  paste0(union(sort(unlist(strsplit(gsub(" ", "", concern_mbd0), ""))), ""), collapse="")]
  }
  
  summary_mbd[, concern_mbd0 := NULL][
    , sec_str_char := sec_str_char]
  
  ### Aggregate the outputs from each sec char 
  subcell_all <- rbind(subcell_all, subcell, use.names=TRUE, fill=TRUE)
  cell_all <- rbind(cell_all, cell, use.names=TRUE, fill=TRUE)
  # cell_all[, ":="(Curr_Chi_Test = as.numeric(Curr_Chi_Test))]
  subcell2_all <- rbind(subcell2_all, subcell2)
  cell2_all <- rbind(cell2_all, cell2)
  
  kpi_report_all <- rbind(kpi_report_all, kpi_report)
  summary_mbd_all <- rbind(summary_mbd_all, summary_mbd)
  
  remove(sf_mbd, sf, pf, subcell, cell, subcell2, cell2, kpi1c, kpi_report, summary_mbd)
}

###################################################
# Generating and formating output

## MBD summry output
### Specify "mbd" column to "Character", in case the variable types are different in two data sets
summary_mbd_all[, mbd := as.character(mbd)]

output_mbdreport <- merge(summary_mbd_all[, list(mbd, sec_str_char, max_score, max_score2, min_Final_Chi_Test, concern_mbd)]
                          , anova_score, by=c("mbd", "sec_str_char"))
### Sample rotation recommendation - ANOVA Score > 1.5, Max Score > 0  and Minimum p value of Chi Sqr test within MBD < 0.05
output_mbdreport[, sample_rotation := ifelse(anova_score>1.5&(max_score>0)&min_Final_Chi_Test<pvalue, "Yes", "No")]
setorder(output_mbdreport, sec_str_char)
output_mbdreport[, ':='(max_score = round(as.numeric(max_score), 2)
                        , anova_score = round(as.numeric(anova_score), 2)
                        , min_Final_Chi_Test = ifelse(as.numeric(min_Final_Chi_Test)<0.0001, "<0.0001", round(as.numeric(min_Final_Chi_Test), 4)))]
output_mbdreport <- output_mbdreport[, list(sec_str_char, mbd, anova_score, max_score, min_Final_Chi_Test, concern_mbd, sample_rotation)]
setnames(output_mbdreport, c("Secondary Characteristic","MBD","Variable Selection Score", "Max. MBD Importance Score"
                             , "Min. Chi-square p-value", "Concern", "Sample Rotation Recommendation"))

## Subcell Allocation output
output_subcell1_all <- subcell_all[, list(sec_str_char, mbd, Cell, sec_var, Count.SF, Current.SS, Final_Exp_SS, min_Retain_SS
                                          , Final_Chi_Sqr_Build, Final_Subcell_Allocation, Additional_SS, Drop_SS, Priority_Rank)]
setorder(output_subcell1_all,sec_str_char)
setnames(output_subcell1_all, c("Secondary Characteristic", "MBD", "Cell", "Secondary Characteristic Value"
                                , " Count SF", "Current SS", "Expected SS", "Min Retain SS", "Final Chi Sqr Build" 
                                , "Final Sub-cell Allocation", "Additional SS","Drop SS", "Sub Cell Priority Rank"))
## Micro Rep Check Output
output_chisq1_all<- cell_all[, list(sec_str_char, mbd, Cell, Count_of_Subcell_cell, Count.SF_cell, Target.SS_cell, Current.SS_cell, min_Retain_SS_cell
                                    , Final_Subcell_Allocation_cell,  Additional_SS, Drop_SS, Curr_Chi_Test, Final_Chi_Test)]
setorder(output_chisq1_all, sec_str_char)
setnames(output_chisq1_all, c("Secondary Characteristic", "MBD", "Cell", "Number of Sub-cells"
                              , "Count SF", "Target SS", "Current SS", "Min. Retain SS", "Final Subcell Allocation"
                              , "Additional SS", "Drop SS", "P-value (current)", "P-value (after NSC)"))
## KPI report output
kpi_report_all <- kpi_report_all[, ':='(panel_estimate = round(as.numeric(panel_estimate), 1)
                                        , universe = round(as.numeric(universe), 1)
                                        , z = round(as.numeric(z), 4)
                                        , completeness1 = round(as.numeric(completeness1), 4)
                                        , completeness2 = round(as.numeric(completeness2), 4)
                                        , ratio = round(as.numeric(ratio), 4)
                                        , ratio2 = round(as.numeric(ratio2), 4)
                                        , score2 = round(as.numeric(score2), 2))]
output_kpireport_all <- kpi_report_all[, list(sec_str_char, mbd, sec_var, panel_size, panel_estimate, universe, z
                                              , completeness1, completeness2, ratio, ratio2, concern, score)]
setnames(output_kpireport_all, c("Secondary Characteristic", "MBD", "Secondary Characteristic Value"
                                 , "Panel Size", "Panel Estimate", "Universe", "Z Score", "Sample Completeness"
                                 , "Universe Completeness", "Ratio 1", "Ratio 2", "Concern", "MBD Importance Score"))
## Subcell Allocation output after sample rotation on recommended cells
if(!is.null(subcell2_all)&!is.null(cell2_all)){
  
  output_subcell2_all <- subcell2_all[, list(sec_str_char, mbd, Cell, sec_var, Count.SF, Current.SS,Final_Exp_SS,min_Retain_SS
                                             , Final_Chi_Sqr_Build, Final_Subcell_Allocation, Additional_SS, Drop_SS, Priority_Rank)]
  setnames(output_subcell2_all, c("Secondary Characteristic", "MBD", "Cell", "Secondary Characteristic Value"
                                  , "Count SF", "Current SS", "Expected SS", "Min Retain SS", "Final Chi Sqr Build"
                                  , "Final Subcell Allocation", "Additional SS", "Drop SS", "Sub Cell Priority Rank"))
  ## Micro Rep Check Output
  output_chisq2_all <- cell2_all[, list(sec_str_char, mbd, Cell, Count_of_Subcell_cell, Count.SF_cell, Target.SS_cell, min_Retain_SS_cell, Current.SS_cell
                                        , Final_Subcell_Allocation_cell, Additional_SS, Drop_SS, Curr_Chi_Test, Final_Chi_Test)]
  
  setnames(output_chisq2_all, c("Secondary Characteristic", "MBD", "Cell", "Number of Sub-cells", "Count SF", "Target SS"
                                , "Min. Retain SS", "Current SS", "Final Subcell Allocation", "Additional SS"
                                , "Drop SS", "P-value (current)", "P-value (after rotation)"))
  
} else {
  output_subcell2_all <- NULL
  output_chisq2_all <- NULL
}

