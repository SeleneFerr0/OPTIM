

######################################################################
###                                                                ###    
###                Micro Rep Analysis - Subcell Allocation         ###
###                        Main Codes (Part2)                      ###
###                Modified on 14 Oct 2016 by Xiaoyu               ###
###                                                                ###
######################################################################

##############################################
# sample Rotation based on the sec. char. and MBD(s) selected by user

## source R script that contains all the customized functions to be used
source("SubcellAllocation_source_functions.R")

## filter subcell and cell with the sec char selected by user
subcell <- subcell_all[sec_str_char==sec_char]
cell <- cell_all[sec_str_char==sec_char]
output_chisq1 <- output_chisq1_all[`Secondary Characteristic`==sec_char, list(`Secondary Characteristic`, `MBD`, `Cell`, `Number of Sub-cells`
                                                                              , `Count SF`, `Target SS`, `Current SS`, `Min. Retain SS`
                                                                              , `Final Subcell Allocation`, `Additional SS`, `Drop SS`
                                                                              , `P-value (current)`, `P-value (after NSC)`)]

output_subcell1 <- output_subcell1_all[`Secondary Characteristic`==sec_char, list(`Secondary Characteristic`, `MBD`, `Cell`, `Secondary Characteristic Value`
                                                                                  , ` Count SF`, `Current SS`, `Expected SS`, `Min Retain SS`
                                                                                  , `Final Chi Sqr Build`, `Final Sub-cell Allocation` 
                                                                                  , `Additional SS`,`Drop SS`, `Sub Cell Priority Rank`)]

output_kpireport <- output_kpireport_all[`Secondary Characteristic`==sec_char, list(`Secondary Characteristic`, `MBD`, `Secondary Characteristic Value`
                                                                                    , `Panel Size`, `Panel Estimate`, `Universe`, `Z Score`
                                                                                    , `Sample Completeness`, `Universe Completeness`, `Ratio 1`, `Ratio 2`
                                                                                    , `Concern`, `MBD Importance Score`)]

if(sum(subcell$sample_rotation)>0&is.null(mbd_rotation)){
  
  ## Automatically perform sample rotation if subcell sample rotation was flagged 
  subcell2 <- subcell2_all[sec_str_char==sec_char]
  cell2 <- cell2_all[sec_str_char==sec_char]
  
  output_chisq2 <- output_chisq2_all[`Secondary Characteristic`==sec_char, list(`Secondary Characteristic`, `MBD`,`Cell`,`Number of Sub-cells`,`Count SF`
                                                                                ,`Target SS`,`Min. Retain SS`,`Current SS`, `Final Subcell Allocation`,
                                                                                `Additional SS`, `Drop SS`,`P-value (current)`, `P-value (after rotation)`)]
  output_chisq2 <- merge(output_chisq2, output_chisq1[, list(Cell, `P-value (after NSC)`)], by="Cell", all.x=T)
  output_chisq2 <- output_chisq2[,list(`Secondary Characteristic`, `MBD`,`Cell`, `Number of Sub-cells`, `Count SF`, `Target SS`, `Min. Retain SS`
                                       , `Current SS`,`Final Subcell Allocation`, `Additional SS`, `Drop SS`, `P-value (current)`
                                       , `P-value (after NSC)`, `P-value (after rotation)`)]
  
  output_subcell2 <- output_subcell2_all[`Secondary Characteristic`==sec_char, list(`Secondary Characteristic`, `MBD`,`Cell`, `Secondary Characteristic Value`
                                                                                    , `Count SF`, `Current SS`, `Expected SS`, `Min Retain SS`, `Final Chi Sqr Build`
                                                                                    , `Final Subcell Allocation`, `Additional SS`,`Drop SS`, `Sub Cell Priority Rank`)]
  
}else if(!is.null(mbd_rotation)){
  
  if(sum(cell[mbd%in%mbd_rotation, Final_Chi_Test]<pvalue)>0){
    
    cell2 <- copy(cell)
    subcell2 <- copy(subcell)
    
    ## Execute sample roation and subcell allocation in MBD selected by user 
    ## Modification from the current Excel template - run sample rotation after every deduction of min retain SS
    ## - set iterations equal to the max of Target SS at cell level
    x <- 1
    while(sum(cell2[mbd%in%mbd_rotation, Final_Chi_Test]<pvalue)>0&x<=max(cell2$Target.SS_cell)){
      
      ### perform sample rotation on MBD(s) selected if cells within the MBD have 
      ### pvalue of Chi-sqr test < pvalue defined by user (e.g. 0.05)
      rotation2 <- rotation_function(as.data.table(subcell2), as.data.table(cell2), pvalue, mbd_rotation)
      subcell2 <- rotation2[[1]]
      cell2 <- rotation2[[2]]
      
      ### Perform subcell allocation with sample rotation  
      if(sum(subcell2$min_Retain_SS_org0>subcell2$min_Retain_SS)>0){
        
        allocation2 <- allocation_function(as.data.table(subcell2), as.data.table(cell2))
        subcell2 <- allocation2[[1]]
        cell2 <- allocation2[[2]]
      }
      x <- x+1
    }
    
    output_subcell2 <- subcell2[, list(sec_str_char, mbd, Cell, sec_var, Count.SF, Current.SS, Final_Exp_SS, min_Retain_SS
                                       , Final_Chi_Sqr_Build, Final_Subcell_Allocation, Additional_SS, Drop_SS, Priority_Rank)]
    
    setnames(output_subcell2, c("Secondary Characteristic", "MBD","Cell","Secondary Characteristic Value"
                                , "Count SF","Current SS", "Expected SS", "Min Retain SS" , "Final Chi Sqr Build"
                                , "Final Subcell Allocation", "Additional SS", "Drop SS", "Sub Cell Priority Rank"))
    
    ## Micro Rep Check Output
    output_chisq2 <- cell2[, list(sec_str_char, mbd,Cell, Count_of_Subcell_cell, Count.SF_cell
                                  , Target.SS_cell, min_Retain_SS_cell, Current.SS_cell, Final_Subcell_Allocation_cell
                                  , Additional_SS, Drop_SS, Curr_Chi_Test, Final_Chi_Test)]
    output_chisq2[, Curr_Chi_Test := as.numeric(Curr_Chi_Test)]
    output_chisq2 <- merge(output_chisq2, output_chisq1[, list(Cell, `P-value (after NSC)`)], by="Cell", all.x=T)
    setnames(output_chisq2, c("Cell", "Secondary Characteristic", "MBD", "Number of Sub-cells"
                              , "Count SF", "Target SS", "Min. Retain SS", "Current SS", "Final Subcell Allocation"
                              , "Additional SS", "Drop SS","P-value (current)", "P-value (after rotation)"
                              , "P-value (after NSC)"))
    output_chisq2 <- output_chisq2[, list(`Secondary Characteristic`,`MBD`, `Cell`, `Number of Sub-cells`, `Count SF`, `Target SS`
                                          , `Current SS`, `Min. Retain SS`, `Final Subcell Allocation`, `Additional SS`, `Drop SS`
                                          , `P-value (current)`, `P-value (after NSC)`, `P-value (after rotation)`)]
  }else{
    subcell2 <- NULL
    cell2 <- NULL
    output_chisq2 <- NULL
    output_subcell2 <- NULL
  }
  
}else{
  subcell2 <- NULL
  cell2 <- NULL
  output_chisq2 <- NULL
  output_subcell2 <- NULL
}

# Merge the selected MBD for sample roation into one string - to prepare for the log output
if(length(mbd_rotation)<=1){
  mbd_rotation_lst <- mbd_rotation
} else {
  mbd_rotation_lst <- mbd_rotation[1]
  for(i in mbd_rotation[-1]){
    mbd_rotation_lst <- paste(mbd_rotation_lst, i, sep = ", ")
  }
}

# Generate the log output
if(min_retain_ss=="Current.SS"){
  min_retain_ss = "Minimum of Current Sample Size and Target Sample Size"
} 
output_log <- data.table(Input = c("Secondary Characteristic", "MBD Selected for Sample Rotation"
                                   , "Min. Retain Sample Size", "P Value", "Time")
                         , Value = list(sec_char, mbd_rotation_lst, min_retain_ss, pvalue, format(Sys.time(), format = "%Y%m%d-%H%M%S")))