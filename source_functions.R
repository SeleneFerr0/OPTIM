

pass_function <- function(subcell, cell){
  subcell[, Allocated_SS_Retian := ifelse(Target.SS_cell==min_Retain_SS_cell, min_Retain_SS
                                          , ifelse(min_Retain_SS > Allocated_SS_prop, min_Retain_SS, Allocated_SS_prop))][
                                            , Final_Subcell_Allocation := round(Allocated_SS_Retian, 0)]
  
  compare_cell_2 <- subcell[, list(Final_Subcell_Allocation_cell = sum(Final_Subcell_Allocation, na.rm = T)), by = Cell]
  
  cell <- merge(cell[, .SD, .SDcols = -c("Final_Subcell_Allocation_cell")], compare_cell_2, by="Cell")
  remove(compare_cell_2)
  
  cell[, ':='(Diff = Target.SS_cell - Final_Subcell_Allocation_cell
              , ft = ifelse(Target.SS_cell - min_Retain_SS_cell > 0, round(Target.SS_cell/(Target.SS_cell - min_Retain_SS_cell), 0)+1, 1)
              , Check = ifelse(Target.SS_cell > Final_Subcell_Allocation_cell, 1, ifelse(Target.SS_cell < Final_Subcell_Allocation_cell, -1, 0)))]
  
  cell[Target.SS_cell < Current.SS_cell, Check := 0]
  
  return(list(subcell, cell))
}



allocation_function <- function(subcell, cell){
  ## subcell and cell are data.table
  cell[, ':='(Diff = Target.SS_cell - Final_Subcell_Allocation_cell
              , ft = ifelse(Target.SS_cell - min_Retain_SS_cell>0, round(Target.SS_cell/(Target.SS_cell - min_Retain_SS_cell), 0)+1, 1)
              , Check = ifelse(Target.SS_cell > Final_Subcell_Allocation_cell, 1, ifelse(Target.SS_cell < Final_Subcell_Allocation_cell, -1, 0)))]
  
  cell[Target.SS_cell < Current.SS_cell, Check := 0]
  ft_max <- max(cell$ft)
  
  #####################################################################

  cell[, Adjust_TSS_iterator := Adjusted.Target.SS_cell + Diff][
    , Adjusted.Target.SS_cell := Adjust_TSS_iterator]
  
  data.table(i = 1:ft_max)[
    , {
      subcell <<- merge(subcell[, .SD, .SDcols = -c("Adjusted.Target.SS_cell")], cell[, list(Cell, Adjusted.Target.SS_cell)], by="Cell")
      
      subcell[, Allocated_SS_prop := ifelse(Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell < Count.SF
                                            , Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell, Count.SF)]
      
      pass <- pass_function(subcell, cell)
      subcell <<- pass[[1]]
      cell <<- pass[[2]]
      
      cell[, Adjust_TSS_iterator := Adjusted.Target.SS_cell + Diff][
        , Adjusted.Target.SS_cell := Adjust_TSS_iterator]
    }
    , by = i
    ]
  
  # for(i in 1:ft_max){
  #   
  #   subcell <- merge(subcell[, Adjusted.Target.SS_cell := NULL], cell[, list(Cell, Adjusted.Target.SS_cell)], by="Cell")
  #   
  #   subcell[, Allocated_SS_prop := ifelse(Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell < Count.SF
  #                                         , Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell, Count.SF)]
  #   
  #   pass <- pass_function(subcell, cell)
  #   subcell <- pass[[1]]
  #   cell <- pass[[2]]
  #   
  #   cell[, Adjust_TSS_iterator := Adjusted.Target.SS_cell + Diff][
  #     , Adjusted.Target.SS_cell := Adjust_TSS_iterator]
  # }
  
  ### Second Pass
  Diff_max <- max(abs(cell$Diff))
  
  data.table(i = 1:(Diff_max+1))[
    , {
      subcell[, Allocated_SS_prop := ifelse(Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell<Count.SF
                                            , Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell, Count.SF)]
      subcell <<- merge(subcell[, .SD, .SDcols = -c("Adjusted.Target.SS_cell")], cell[, list(Cell, Adjusted.Target.SS_cell)], by="Cell")
      
      pass <- pass_function(subcell, cell)
      subcell <<- pass[[1]]
      cell <<- pass[[2]]
      
      cell[, Adjust_TSS_iterator := Adjusted.Target.SS_cell + Check][
        , Adjusted.Target.SS_cell := Adjust_TSS_iterator]   
    }
    , by = i
    ]
  
  # for(i in 1:(Diff_max+1)){
  #   
  #   subcell[, Allocated_SS_prop := ifelse(Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell<Count.SF
  #                                         , Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell, Count.SF)]
  #   subcell <- merge(subcell[, Adjusted.Target.SS_cell := NULL], cell[, list(Cell, Adjusted.Target.SS_cell)], by="Cell")
  #   
  #   pass <- pass_function(subcell, cell)
  #   subcell <- pass[[1]]
  #   cell <- pass[[2]]
  #   
  #   cell[, Adjust_TSS_iterator := Adjusted.Target.SS_cell + Check][
  #     , Adjusted.Target.SS_cell := Adjust_TSS_iterator]
  # }
  
  ### Third Pass
  n <- 2
  while(sum(cell$Check!=0)>0&n<2^100){
    
    subcell <- merge(subcell[, .SD, .SDcols = -c("Adjusted.Target.SS_cell")], cell[, list(Cell, Adjusted.Target.SS_cell)], by="Cell")
    
    subcell[, Allocated_SS_prop := ifelse(Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell<Count.SF
                                          , Randomized_Count_SF*Adjusted.Target.SS_cell/Randomized_Count_SF_cell, Count.SF)]
    
    pass <- pass_function(subcell, cell)
    subcell <- pass[[1]]
    cell <- pass[[2]]
    
    cell[, Adjust_TSS_iterator := Adjusted.Target.SS_cell + Check/n][
      , Adjusted.Target.SS_cell := Adjust_TSS_iterator]
    
    n <- n*2
  }
  
  ### number of stores to be dropped / added
  subcell[, ':='(Drop_SS = ifelse(Final_Subcell_Allocation < Current.SS, Final_Subcell_Allocation - Current.SS, 0)
                 , Additional_SS = ifelse(Final_Subcell_Allocation > Current.SS, Final_Subcell_Allocation - Current.SS, 0))]
  drop_add_cell <- subcell[, list(Drop_SS = sum(Drop_SS), Additional_SS = sum(Additional_SS)), by = Cell]
  
  if(sum(names(cell)%in%c("Drop_SS", "Additional_SS"))>0){
    cell[, ':='(Drop_SS = NULL, Additional_SS = NULL)]
  }
  cell <- merge(cell, drop_add_cell, by="Cell", all.x=T)
  
  subcell[, ':='(zUniv_prop = Count.SF*zUniv_cell/Count.SF_cell)][
    , ':='(xUniv_prop = zUniv_prop*Avg.ACV_cell)][
      , ':='(Priority_Index = (xUniv_prop/sum(xUniv_prop))*(log((Final_Subcell_Allocation + 1)/(Current.SS + 1))))][
        , Priority_Index := ifelse(Priority_Index==0, NA, Priority_Index)][
          , Priority_Rank := ifelse(!is.na(Priority_Index), rank(Priority_Index, ties.method="max"), 0)]
  
  rank_max <- max(subcell$Priority_Rank)
  
  subcell[, Priority_Rank := ifelse(Priority_Rank==0, "", rank_max - Priority_Rank + 1)][
    , Final_Exp_SS := Ideal.SS_cell*Count.SF/Count.SF_cell][
      , Final_Chi_Sqr_Build := ifelse(Final_Exp_SS>0, (Final_Subcell_Allocation - Final_Exp_SS)^2/Final_Exp_SS, 0)][
        , Curr_Exp_SS := Current.SS_cell*Count.SF/Count.SF_cell][
          , Curr_Chi_Sqr_Build := ifelse(Curr_Exp_SS>0, (Current.SS - Curr_Exp_SS)^2/Curr_Exp_SS, 0)]
  
  Final_Chi_Sqr_Build_cell <- subcell[, list(Final_Chi_Sqr_Build_cell = sum(Final_Chi_Sqr_Build)), by = Cell]
  
  if("Final_Chi_Sqr_Build_cell"%in%names(cell)){
    cell <- cell[, .SD, .SDcols = -c("Final_Chi_Sqr_Build_cell")]
  }
  cell <- merge(cell, Final_Chi_Sqr_Build_cell, by="Cell")
  
  cell[, Final_Chi_Test := ifelse(Count_of_Subcell_cell>1&Final_Subcell_Allocation_cell>0
                                  , 1 - pchisq(Final_Chi_Sqr_Build_cell, df = Count_of_Subcell_cell - 1, lower.tail=T), 0)]
  
  Curr_Chi_Sqr_Build_cell <- subcell[, list(Curr_Chi_Sqr_Build_cell = sum(Curr_Chi_Sqr_Build)), by = Cell]
  
  if("Curr_Chi_Sqr_Build_cell"%in%names(cell)){
    cell <- cell[, .SD, .SDcols = -c("Curr_Chi_Sqr_Build_cell")]
  }
  cell <- merge(cell, Curr_Chi_Sqr_Build_cell, by="Cell")
  
  cell[, Curr_Chi_Test := ifelse(Count_of_Subcell_cell>1&Current.SS_cell>0
                                 , 1 - pchisq(Curr_Chi_Sqr_Build_cell, df = Count_of_Subcell_cell - 1, lower.tail=T), 0)]
  
  ## Update Final_Chi_Test in subcell
  if("Final_Chi_Test"%in%names(subcell)){
    subcell[,  Final_Chi_Test := NULL]
  } 
  
  subcell <- merge(subcell, cell[, list(Cell, Final_Chi_Test)], by = "Cell")
  
  return(list(subcell,cell))
}



rotation_function <- function(subcell, cell, pvalue, mbd_rotation){
  

  subcell[, min_Retain_SS_org0 := min_Retain_SS][
    , min_Retain_SS_org := min_Retain_SS]     # set min Retain SS after sample rotation as min Retain SS to start with
  

  if(is.null(mbd_rotation)){
    
    chisq_ck_cell <- unique(subcell[sample_rotation==1, Cell])
    
    #### subcell with sample rotation flag
    chisq_ck_subcell <- subcell[sample_rotation==1]
  }
  )
  if(!is.null(mbd_rotation)){
    
    chisq_ck_cell <- unique(subcell[(mbd%in%mbd_rotation)&(Final_Chi_Test<pvalue), Cell])

    chisq_ck_subcell <- subcell[(mbd%in%mbd_rotation)&(Final_Chi_Test<pvalue)]
  }
  
  data.table(i = chisq_ck_cell)[
    , {
      chisq_ck_subcell_i <- chisq_ck_subcell[Cell==as.character(i)]
      subcell_i <- chisq_ck_subcell_i[(Final_Chi_Sqr_Build==max(Final_Chi_Sqr_Build))&(min_Retain_SS>Final_Exp_SS), Subcell] 
      subcell[Subcell%in%subcell_i, min_Retain_SS:=min_Retain_SS_org-1]   
    }
    , by = i
    ]
  

  
  if(sum(subcell$min_Retain_SS_org>subcell$min_Retain_SS)>0){

    subcell[,min_Retain_SS_cell := sum(min_Retain_SS),by=Cell][
      , Allocated_SS_Retian := ifelse(Target.SS_cell==min_Retain_SS_cell, min_Retain_SS
                                      , ifelse(min_Retain_SS>Allocated_SS_prop, min_Retain_SS, Allocated_SS_prop))][
                                        , Final_Subcell_Allocation := round(Allocated_SS_Retian, 0)][
                                          , Final_Chi_Sqr_Build := ifelse(Final_Exp_SS>0, (Final_Subcell_Allocation - Final_Exp_SS)^2/Final_Exp_SS, 0)]
    
    
    #### Update Final_Chi_Sqr_Build_cell, Final_Subcell_Allocation_cell in cell
    Final <- subcell[, list(Final_Chi_Sqr_Build_cell = sum(Final_Chi_Sqr_Build, na.rm = T)
                            , Final_Subcell_Allocation_cell = sum(Final_Subcell_Allocation, na.rm = T))
                     , by=Cell]
    cell[, ":="(Final_Chi_Sqr_Build_cell = NULL, Final_Subcell_Allocation_cell = NULL)]
    
    cell <- merge(cell, Final, by = "Cell", all.x = T)
    

    cell[, Final_Chi_Test := NULL][
      , Final_Chi_Test := ifelse (Count_of_Subcell_cell>1&Final_Subcell_Allocation_cell>0
                                  , 1-pchisq(Final_Chi_Sqr_Build_cell, df=Count_of_Subcell_cell-1, lower.tail=T), 0)]

    if("Final_Chi_Test"%in%names(subcell)){
      subcell[,  Final_Chi_Test := NULL]
    } 
    subcell <- merge(subcell, cell[, list(Cell, Final_Chi_Test)], by = "Cell")
  }
  

  cell[, min_Retain_SS_cell := NULL]
  cell <- merge(cell, unique(subcell[, list(Cell, min_Retain_SS_cell)]), by="Cell", all=F)
  subcell <- subcell
  cell <- cell
  return(list(subcell, cell))
}



output_chisq_func <- function(output_chisq, pvalue){
  
  output_chisq <- as.data.frame(output_chisq)
  
  output_chisq[,"P-value (current)"] <- ifelse(output_chisq[, "P-value (current)"]<0.0001, "<0.0001*"
                                               , ifelse(output_chisq[, "P-value (current)"]<pvalue&output_chisq[, "P-value (current)"]>=0.0001
                                                        , paste0(round(as.numeric(output_chisq[, "P-value (current)"]), 4), "*")
                                                        , round(as.numeric(output_chisq[, "P-value (current)"]), 4)))
  if("P-value (after NSC)"%in%names(output_chisq)){
    output_chisq[, "P-value (after NSC)"] <- ifelse(output_chisq[, "P-value (after NSC)"]<0.0001,"<0.0001*"
                                                    , ifelse(output_chisq[, "P-value (after NSC)"]<pvalue&output_chisq[, "P-value (after NSC)"]>=0.0001
                                                             , paste0(round(as.numeric(output_chisq[, "P-value (after NSC)"]), 4), "*")
                                                             , round(as.numeric(output_chisq[, "P-value (after NSC)"]), 4)))
  }
  if("P-value (after rotation)"%in%names(output_chisq)){
    output_chisq[, "P-value (after rotation)"] <- ifelse(output_chisq[, "P-value (after rotation)"]<0.0001, "<0.0001*"
                                                         , ifelse(output_chisq[, "P-value (after rotation)"]<pvalue&output_chisq[, "P-value (after rotation)"]>=0.0001
                                                                  , paste0(round(as.numeric(output_chisq[, "P-value (after rotation)"]), 4), "*")
                                                                  , round(as.numeric(output_chisq[, "P-value (after rotation)"]), 4)))
  }
  return(output_chisq)
}



pvalue_style <- function(wb){
  
  ## if p-value <0.01, highlighted in red; if p-value 0.05-0.1, highlighted in yellow
  pvalue1 <- createCellStyle(wb, name = "pvalue1")
  setFillPattern(pvalue1, fill = XLC$"FILL.SOLID_FOREGROUND")
  setFillForegroundColor(pvalue1, color = XLC$"COLOR.RED")
  ## if p-value 0.01-0.05, highlighted in orange; 
  pvalue2 <- createCellStyle(wb, name = "pvalue2")
  setFillPattern(pvalue2, fill = XLC$"FILL.SOLID_FOREGROUND")
  setFillForegroundColor(pvalue2, color = XLC$"COLOR.ORANGE")
  ## if p-value 0.01-0.05, highlighted in orange; 
  pvalue3 <- createCellStyle(wb, name = "pvalue3")
  setFillPattern(pvalue3, fill = XLC$"FILL.SOLID_FOREGROUND")
  setFillForegroundColor(pvalue3, color = XLC$"COLOR.YELLOW")
  
  return(list(pvalue1,pvalue2,pvalue3))
}



highlight_func <- function(wb, data, colIndex, pvalue, pvaluestyle, sheetname){
  data <- as.data.frame(data)
  rowIndex <- which(data[,colIndex] < pvalue) + 1
  if(length(rowIndex)>0){
    setCellStyle(wb, sheet = sheetname, row = rowIndex, col = colIndex
                 , cellstyle = pvaluestyle)
  }
  return(wb)
}


sample_rotation_note <- function(pvalue, mbdreport, mbd_rotation, output_chisq1, output_chisq2, text2){
  
  if(is.null(mbd_rotation)&is.null(output_chisq2)){
    text1 <- c("No sample rotation is recommended")
  }
  
  if(is.null(mbd_rotation)&!is.null(output_chisq2)){
    text1 <- c("Sample rotation is performed on recommended cells")
  }
  
  if(!is.null(mbd_rotation)&is.null(output_chisq2)){
    mbd_rotation_list <- NULL
    for(i in mbd_rotation){
      mbd_rotation_list <- paste0(mbd_rotation_list, " ", i)
    }
    mbd_rotation_list <- paste0("[",mbd_rotation_list, " ]")
    text1 <- paste("There is no cell in the selected MBD(s):", mbd_rotation_list, ", with p value of Chi Square Test <", pvalue, ", no sample rotation has
                   been performed")
  }
  
  if(!is.null(mbd_rotation)&!is.null(output_chisq2)){
    output_chisq1 <- as.data.table(output_chisq1)
    mbd_rotation_list_yes <- NULL
    mbd_rotation_list_no <- NULL
    for(i in mbd_rotation){
      if(sum(output_chisq1[MBD==i][["P-value (after NSC)"]]<pvalue)>0){
        mbd_rotation_list_yes <- paste0(mbd_rotation_list_yes, " ", i)
      }else{
        mbd_rotation_list_no <- paste0(mbd_rotation_list_no, " ", i)
      }
    }
    
    if(!is.null(mbd_rotation_list_yes)){
      mbd_rotation_list_yes <- paste0("[",mbd_rotation_list_yes, " ]")
      text1 <- paste("Sample rotation is performed on selected MBD(s):", mbd_rotation_list_yes)
    }else{
      text1 <- NULL
    }
    
    if(!is.null(mbd_rotation_list_no)){
      mbd_rotation_list_no <- paste0("[",mbd_rotation_list_no, " ]")
      text2 <-paste("There is no cell in the selected MBD(s):", mbd_rotation_list_no, ", with p value of Chi Square Test <", pvalue, ", no sample rotation has
                    been performed")
    }else{
      text2 <- NULL
    }
  }
  
  return(list(text1, text2))
  }
