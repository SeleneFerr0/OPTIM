library(plyr)

source("lib_mra_anova.R")

# mbds
mbds <- c("total",unique(chars.DT[,mbd]))

# merging sales and chars
setkey(sales.DT,outlet)
setkey(chars.DT,outlet)

hk.DT = chars.DT[sales.DT, nomatch = 0]
hk.DT$share = hk.DT$brand_sales/hk.DT$category_sales

#all_chars <- names(hk.DT)[!(names(hk.DT) %in% c(
#  "outlet","mbd","type","category","brand","category_sales","brand_sales","outlet_sales","share"))]
all_chars <- names(hk.DT)[c(which(grepl("sec_var_",names(hk.DT))))]

n_all_chars <- length(all_chars)

# coping with NAs !!!verify!!!
for (col in all_chars) hk.DT[is.na(get(col)), (col) := 0]

hk = as.data.frame(hk.DT)

# setting the type of individual characteristics to factor
for (col in all_chars) hk[,col] = as.factor(hk[,col])

WGT.POWER = 0.5

for (cur_mbd in mbds){
  for (CRIT in all_chars){
    incProgress(1/(n_all_chars*length(mbds)), detail = paste("Doing part:", cur_mbd, ", ",CRIT))
    
    filter <- cur_mbd
    if (cur_mbd == "total") { filter <- mbds}
    
    zz <- ddply(
      hk[hk$mbd %in% filter,]
      , .(category, brand)
      , .fun = cat.brand.char.procfunction2
      , cname = CRIT
      , misval = "<NA>"
    )
    
    if (is.numeric(zz$pval)){
      score = as.integer(as.character(
        cut(
          zz$pval
          , c(0,0.001,0.01,0.1,1)
          , labels = c(3,2,1,0)
          , include.lowest = TRUE
        )
      ))
      
      avg.score = sum(score * zz$sales^WGT.POWER) / sum(zz$sales^WGT.POWER)
      
      zz[zz$pval == 0,"pval"] = min(zz[zz$pval > 0,"pval"])
      
      #ofile = file.path(path2data,paste0("anovaReport_",cur_mbd,"_",CRIT,".html"))
      #rmarkdown::render('doc_anova_report.Rmd',output_file = ofile)
      
      write.csv(
        zz
        , file.path(path,paste0("pval_",gsub("/","-",cur_mbd),"_",CRIT,".csv"))
        , quote = FALSE
        , row.names = FALSE
      )
      
      codebook.DT[mbd == cur_mbd & variable == CRIT,`:=`(score,avg.score)]
    }
  }
}

# ix=(hk$category_code=="FCMP") & (hk$brand_code == "Anlene")
# wwdf=hk[ix,]
# cat.brand.char.procfunction2(wwdf,CRIT,5)
