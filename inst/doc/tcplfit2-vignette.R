params <-
list(my_css = "css/rmdformats.css")

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, class.source="scroll-100", warning = FALSE, message = FALSE-------
# Primary Packages #
library(tcplfit2)
library(tcpl)
# Data Formatting Packages #
library(data.table)
library(DT)
library(htmlTable)
library(dplyr)
library(stringr)
# Plotting Packages #
library(ggplot2)
library(gridExtra)

## ----example1, warning=FALSE--------------------------------------------------
# tested concentrations
  conc <- list(.03,.1,.3,1,3,10,30,100)
# observed responses at respective concentrations
  resp <- list(0,.2,.1,.4,.7,.9,.6, 1.2)
# row object with relevant parameters
  row = list(conc = conc, resp = resp, bmed = 0, cutoff = 1, onesd = .5,name="some chemical")
# execute concentration-response modeling through potency estimation
  res <- concRespCore(row,
                      fitmodels = c("cnst", "hill", "gnls",
                                    "poly1", "poly2", "pow", "exp2", "exp3",
                                        "exp4", "exp5"),
                      conthits = T)

## ----echo=FALSE---------------------------------------------------------------
htmlTable::htmlTable(head(res),
        align = 'l',
        align.header = 'l',
        rnames = FALSE  ,
        css.cell =  ' padding-bottom: 5px;  vertical-align:top; padding-right: 10px;min-width: 5em ')

## ----example1 plot, fig.height = 4.55, fig.width = 8--------------------------
# plot the winning curve from example 1, add a title
concRespPlot2(res, log_conc = TRUE) + ggtitle("Example 1: Chemical A")

## ----example2_load, warning=FALSE---------------------------------------------
# read in the data
# Loading in the level 3 example data set from invitrodb
  data("mc3")
  head(mc3)

## ----example2, warning=FALSE--------------------------------------------------
# determine the background variation
# chosen as logc <= -2 in this example but will be assay/application specific
  temp <- mc3[mc3$logc<= -2,"resp"]
  bmad <- mad(temp)
  onesd <- sd(temp)
  cutoff <- 3*bmad

# select six chemical samples. Note that there may be more than one sample processed for a given chemical
  spid.list <- unique(mc3$spid)
  spid.list <- spid.list[1:6]
  
# create empty objects to store results and plots
  model_fits <- NULL
  result_table <- NULL
  plt_lst <- NULL

# loop over the samples to perform concentration-response modeling & hitcalling
  for(spid in spid.list) {
    # select the data for just this sample
    temp <- mc3[is.element(mc3$spid,spid),]

    # The data file stores concentrations in log10 units, so back-transform
    conc <- 10**temp$logc
    # Save the response values
    resp <- temp$resp

    # pull out all of the chemical identifiers and the assay name
    dtxsid <- temp[1,"dtxsid"]
    casrn <- temp[1,"casrn"]
    name <- temp[1,"name"]
    assay <- temp[1,"assay"]
    
    # Execute curve fitting
    # Input concentrations, responses, cutoff, a list of models to fit, and other model fitting requirements
    # force.fit is set to true so that all models will be fit regardless of cutoff
    # bidirectional = FALSE indicates only fit models in the positive direction.
    # if using bidirectional = TRUE the coff only needs to be specified in the positive direction.
    model_fits[[spid]] <- tcplfit2_core(conc, resp, cutoff, force.fit = TRUE, 
                                        fitmodels = c("cnst", "hill", "gnls", 
                                                      "poly1", "poly2", "pow", 
                                                      "exp2","exp3", "exp4", "exp5"),
                                        bidirectional = FALSE)
    # Get a plot of all curve fits
    plt_lst[[spid]] <- plot_allcurves(model_fits[[spid]], 
                                      conc = conc, resp = resp, log_conc = TRUE)
    
    # Pass the output from 'tcplfit2_core' to 'tcplhit2_core' along with
    # cutoff, onesd, and any identifiers
    out <- tcplhit2_core(model_fits[[spid]], conc, resp, bmed = 0,
                         cutoff = cutoff, onesd = onesd, 
                         identifiers = c(dtxsid = dtxsid, casrn = casrn, 
                                         name = name, assay = assay))
    # store all results in one table
    result_table <- rbind(result_table,out)
  }


## ----example 2 fit results----------------------------------------------------
# shows the structure of the output object from tcplfit2_core (only top level)
str(model_fits[[1]],max.lev = 1)

## -----------------------------------------------------------------------------
str(model_fits[[1]][["hill"]])

## ----example2 plot1, fig.height = 9, fig.width = 7----------------------------
grid.arrange(grobs=plt_lst,ncol=2)

## -----------------------------------------------------------------------------
htmlTable::htmlTable(head(result_table),
        align = 'l',
        align.header = 'l',
        rnames = FALSE  ,
        css.cell =  ' padding-bottom: 5px;  vertical-align:top; padding-right: 10px;min-width: 5em ')

## ----example2 plot2-----------------------------------------------------------
# plot the first row
concRespPlot2(result_table[1,],log_conc = TRUE) + 
  ggtitle(paste(result_table[1,"dtxsid"], result_table[1,"name"]))

## ----example3_init, fig.height = 6, fig.width = 7, message=FALSE, warning = FALSE----
# Loading in the Level 0 example data set from invitrodb
data("mc0")
data.table::setDTthreads(2)
dat <- mc0

## ----echo=FALSE---------------------------------------------------------------
htmlTable::htmlTable(head(dat[wllt=='t',]),
        align = 'l',
        align.header = 'l',
        rnames = FALSE  ,
        css.cell =  ' padding-bottom: 5px;  vertical-align:top; padding-right: 10px;min-width: 5em ')

## ----example3_cndx, class.source="scroll-100", fig.height = 6, fig.width = 7, warning=FALSE----
# Order by the following columns
setkeyv(dat, c('acid', 'srcf', 'apid', 'coli', 'rowi', 'spid', 'conc'))

# Define a temporary replicate ID (rpid) column for test compound wells
# rpid consists of the sample ID, well type (wllt), source file, assay plate ID, and 
# concentration.
nconc <- dat[wllt == "t" , ## denotes test well as the well type (wllt)
             list(n = lu(conc)), #total number of unique concentrations
             by = list(acid, apid, spid)][ , list(nconc = min(n)), by = acid]
dat[wllt == "t" & acid %in% nconc[nconc > 1, acid],
    rpid := paste(acid, spid, wllt, srcf, apid, "rep1", conc, sep = "_")]
dat[wllt == "t" & acid %in% nconc[nconc == 1, acid],
    rpid := paste(acid, spid, wllt, srcf, "rep1", conc, sep = "_")]

# Define rpid column for non-test compound wells
dat[wllt != "t",
    rpid := paste(acid, spid, wllt, srcf, apid, "rep1", conc, sep = "_")]

# set the replicate index (repi) based on rowid 
# increment repi every time a replicate ID is duplicated
dat[, dat_rpid := rowid(rpid)]
dat[, rpid := sub("_rep[0-9]+.*", "",rpid, useBytes = TRUE)]
dat[, rpid := paste0(rpid,"_rep",dat_rpid)]

# For each replicate, define concentration index
# by ranking the unique concentrations
indexfunc <- function(x) as.integer(rank(unique(x))[match(x, unique(x))])
# the := operator is a data.table function to add/update rows
dat[ , cndx := indexfunc(conc), by = list(rpid)]

## ----echo=FALSE---------------------------------------------------------------
# tcplConf(user="_dataminer", pass="pass", db="invitrodb", drvr="MySQL", host="ccte-mysql-res.epa.gov")

## ----example3_mc2, fig.height = 6, fig.width = 7------------------------------
# If no adjustments are required for the data, the corrected value (cval) should be set as original rval
dat[,cval := rval]

# Poor well quality (wllq) wells should be removed
dat <- dat[!wllq == 0,]

##Fitting generally cannot occur if response values are NA therefore values need to be removed
dat <- dat[!is.na(cval),]

## ----echo=FALSE---------------------------------------------------------------
htmlTable::htmlTable(head(tcpl::tcplMthdList(3)),
        align = 'l',
        align.header = 'l',
        rnames = FALSE  ,
        css.cell =  ' padding-bottom: 5px;  vertical-align:top; padding-right: 10px;min-width: 5em ')

## ----example3 normalize-------------------------------------------------------
# calculate bval of the median of all the wells that have a type of n
dat[, bval := median(cval[wllt == "n"]), by = list(apid)]
# calculate pval based on the wells that have type of m or o excluding any NA wells
dat[, pval := median(cval[wllt %in% c("m","o")], na.rm = TRUE), by = list(apid, wllt, conc)]
# take pval as the minimum per assay plate (apid)
dat[, pval := min(pval, na.rm = TRUE), by = list(apid)]

# Calculate normalized responses
dat[, resp := ((cval - bval)/(pval - bval) * 100)]

## ----echo=FALSE---------------------------------------------------------------
htmlTable::htmlTable(head(tcpl::tcplMthdList(4)),
        align = 'l',
        align.header = 'l',
        rnames = FALSE  ,
        css.cell =  ' padding-bottom: 5px;  vertical-align:top; padding-right: 10px;min-width: 5em ')

## ----example3_get_bmad.and.onesd----------------------------------------------
bmad <- mad(dat[cndx %in% c(1, 2) & wllt == "t", resp])
onesd <- sd(dat[cndx %in% c(1, 2) & wllt == "t", resp])

## ----example3_fitting, fig.height = 6, fig.width = 7--------------------------
#do tcplfit2 fitting
myfun <- function(y) {
  res <- tcplfit2::tcplfit2_core(y$conc,
                          y$resp,
                          cutoff = 3*bmad,
                          bidirectional = TRUE,
                          verbose = FALSE,
                          force.fit = TRUE,
                          fitmodels = c("cnst", "hill", "gnls", "poly1",
                                        "poly2", "pow", "exp2", "exp3",
                                        "exp4", "exp5")
                          )
  list(list(res)) #use list twice because data.table uses list(.) to look for values to assign to columns
}

## ----example3_fitting_full, eval=FALSE, echo = FALSE--------------------------
#  # only want to run tcplfit2 for test wells in this case
#  # this chunk doesn't run, fit the curves on the subset below
#  dat[wllt == 't',params:= myfun(.SD), by = .(spid)]

## ----example3_fitting_subset--------------------------------------------------
# create a subset that contains 6 samples and run curve fitting
subdat <- dat[spid %in% unique(spid)[10:15],]
subdat[wllt == 't',params:= myfun(.SD), by = .(spid)]

## ----example3_hitcalling, fig.height = 6, fig.width = 7-----------------------
myfun2 <- function(y) {
  res <- tcplfit2::tcplhit2_core(params = y$params[[1]],
                                 conc = y$conc,
                                 resp = y$resp,
                                 cutoff = 3*bmad,
                                 onesd = onesd
                                 )
  list(list(res))
}

# continue with hitcalling
res <- subdat[wllt == 't', myfun2(.SD), by = .(spid)]

# pivot wider
res_wide <- rbindlist(Map(cbind, spid = res$spid, res$V1))


## ----echo=FALSE---------------------------------------------------------------
htmlTable::htmlTable(head(res_wide),
        align = 'l',
        align.header = 'l',
        rnames = FALSE  ,
        css.cell =  ' padding-bottom: 5px;  vertical-align:top; padding-right: 10px;min-width: 5em ')

## ----example3_plot, fig.height = 8, fig.width = 7-----------------------------
# allocate a place-holder object
  plt_list <- NULL
# plot results using `concRespPlot`
  for(i in 1:nrow(res_wide)){
    plt_list[[i]] <- concRespPlot2(res_wide[i,])
  }
# compile and display winning model plots for concentration-response series
  grid.arrange(grobs=plt_list,ncol=2)

## ----example 4 lower, warning=FALSE-------------------------------------------
# We'll use data from mc3 in this section
data("mc3")

# determine the background variation
# background is defined per the assay.  In this case we use logc <= -2
# However, background should be defined in a way that makes sense for your application
temp <- mc3[mc3$logc<= -2,"resp"]
bmad <- mad(temp)
onesd <- sd(temp)
cutoff <- 3*bmad

# load example data
spid <- unique(mc3$spid)[94]
ex_df <- mc3[is.element(mc3$spid,spid),]

# The data file has stored concentration in log10 form, fix it 
conc <- 10^ex_df$logc # back-transforming concentrations on log10 scale
resp <- ex_df$resp

# modify the data for demonstration purposes 
conc2 <- conc[conc>0.41]
resp2 <- resp[which(conc>0.41)]

# pull out all of the chemical identifiers and the name of the assay
dtxsid <- ex_df[1,"dtxsid"]
casrn <- ex_df[1,"casrn"]
name <- ex_df[1,"name"]
assay <- ex_df[1,"assay"]

# create the row object
row_low <- list(conc = conc2, resp = resp2, bmed = 0, cutoff = cutoff, onesd = onesd,
            assay=assay, dtxsid=dtxsid,casrn=casrn,name=name)

# run the concentration-response modeling for a single sample
res_low <- concRespCore(row_low,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", 
                                          "pow", "exp2", "exp3", "exp4", "exp5"), 
                        bidirectional=F)

concRespPlot2(res_low, log_conc = T) + 
  geom_rect(aes(xmin = log10(res_low[1, "bmdl"]),
                xmax = log10(res_low[1, "bmdu"]),ymin = 0,ymax = 30),
            alpha = 0.05,fill = "skyblue") + 
  geom_segment(aes(x = log10(res_low[, "bmd"]),
                   xend = log10(res_low[, "bmd"]), y = 0, 
                   yend = 30),col = "blue")

## ----example 4 lower-res------------------------------------------------------
# function results
res_low['Min. Conc.'] <- min(conc2)
res_low['Name'] <- name
res_low[1, c("Min. Conc.", "bmd", "bmdl", "bmdu")] <- round(res_low[1, c("Min. Conc.", "bmd", "bmdl", "bmdu")], 3)

## ----example_4_table, echo=FALSE----------------------------------------------
DT::datatable(res_low[1, c("Name","Min. Conc.", "bmd", "bmdl", "bmdu")],rownames = FALSE)

## ----example 4 lower-demo-----------------------------------------------------
# using the argument to set a lower bound for BMD
res_low2 <- concRespCore(row_low,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", 
                                           "pow", "exp2", "exp3", "exp4", "exp5"), 
                         conthits = T, aicc = F, bidirectional=F, bmd_low_bnd = 0.8)

## ----example 4 new lower-res--------------------------------------------------
# print out the new results
# include previous results side by side for comparison 
res_low2['Min. Conc.'] <- min(conc2)
res_low2['Name'] <- paste(name, "after `bounding`", sep = "-")
res_low['Name'] <- paste(name, "before `bounding`", sep = "-")
res_low2[1, c("Min. Conc.", "bmd", "bmdl", "bmdu")] <- round(res_low2[1, c("Min. Conc.", "bmd", "bmdl", "bmdu")], 3)

output_low <- rbind(res_low[1, c('Name', "Min. Conc.", "bmd", "bmdl", "bmdu")], 
                    res_low2[1, c('Name', "Min. Conc.", "bmd", "bmdl", "bmdu")])

## ----example_4_lower_res_table, echo = FALSE----------------------------------
DT::datatable(output_low,rownames = FALSE)

## ----example 4 lower plot, class.source="scroll-100"--------------------------
# generate some concentration for the fitted curve 
logc_plot <- seq(from=-3,to=2,by=0.05)
conc_plot <- 10**logc_plot

# initiate the plot
plot(conc2,resp2,xlab="conc (uM)",ylab="Response",xlim=c(0.001,100),ylim=c(-5,60),
       log="x",main=paste(name,"\n",assay),cex.main=0.9)

# add vertical lines to mark the minimum concentration in the data and the lower threshold set by bmd_low_bnd
abline(v=min(conc2), lty = 1, col = "brown", lwd = 2)
abline(v=res_low2$bmd, lty = 2, col = "darkviolet", lwd = 2)

# add markers for BMD and its boundaries before `bounding`
lines(c(res_low$bmd,res_low$bmd),c(0,50),col="green",lwd=2)
rect(xleft=res_low$bmdl,ybottom=0,xright=res_low$bmdu,ytop=50,col=rgb(0,1,0, alpha = .5), border = NA)
points(res_low$bmd, 0, pch = "x", col = "green")

# add markers for BMD and its boundaries after `bounding`
lines(c(res_low2$bmd,res_low2$bmd),c(0,50),col="blue",lwd=2)
rect(xleft=res_low2$bmdl,ybottom=0,xright=res_low2$bmdu,ytop=50,col=rgb(0,0,1, alpha = .5), border = NA)
points(res_low2$bmd, 0, pch = "x", col = "blue")

# add the fitted curve
lines(conc_plot, exp4(ps = c(res_low$tp, res_low$ga), conc_plot))
legend(1e-3, 60, legend=c("Lowest Dose Tested", "Boundary", "BMD-before", "BMD-after"),
       col=c("brown", "darkviolet", "green", "blue"), lty=c(1,2,1,1))

## ----example 5 upper----------------------------------------------------------
# load example data
spid <- unique(mc3$spid)[26]
ex_df <- mc3[is.element(mc3$spid,spid),]

# The data file has stored concentration in log10 form, so fix that
conc <- 10**ex_df$logc # back-transforming concentrations on log10 scale
resp <- ex_df$resp

# pull out all of the chemical identifiers and the name of the assay
dtxsid <- ex_df[1,"dtxsid"]
casrn <- ex_df[1,"casrn"]
name <- ex_df[1,"name"]
assay <- ex_df[1,"assay"]

# create the row object
row_up <- list(conc = conc, resp = resp, bmed = 0, cutoff = cutoff, onesd = onesd,assay=assay,
            dtxsid=dtxsid,casrn=casrn,name=name)

# run the concentration-response modeling for a single sample
res_up <- concRespCore(row_up,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", 
                                         "pow", "exp2", "exp3", "exp4", "exp5"), 
                       conthits = T, aicc = F, bidirectional=F)

concRespPlot2(res_up, log_conc = T)

## ----example 5 upper-res------------------------------------------------------
# max conc
res_up['Max Conc.'] <- max(conc)
res_up['Name'] <- name
res_up[1, c("Max Conc.", "bmd", "bmdl", "bmdu")] <- round(res_up[1, c("Max Conc.", "bmd", "bmdl", "bmdu")], 3)
# function results

## ----example_5_table, echo = FALSE--------------------------------------------
DT::datatable(res_up[1, c('Name','Max Conc.', "bmd", "bmdl", "bmdu")],rownames = FALSE)

## ----example upper-demo-------------------------------------------------------
# using bmd_up_bnd = 2
res_up2 <- concRespCore(row_up,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", 
                                          "pow", "exp2", "exp3", "exp4", "exp5"), 
                        conthits = T, aicc = F, bidirectional=F, bmd_up_bnd = 2)

## ----example upper-2----------------------------------------------------------
# print out the new results
# include previous results side by side for comparison 
res_up2['Max Conc.'] <- max(conc)
res_up2['Name'] <- paste(name, "after `bounding`", sep = "-")
res_up['Name'] <- paste(name, "before `bounding`", sep = "-")
res_up2[1, c("Max Conc.", "bmd", "bmdl", "bmdu")] <- round(res_up2[1, c("Max Conc.", "bmd", "bmdl", "bmdu")], 3)

output_up <- rbind(res_up[1, c('Name', "Max Conc.", "bmd", "bmdl", "bmdu")], 
                   res_up2[1, c('Name', "Max Conc.", "bmd", "bmdl", "bmdu")])

## ----example_upper_2_table, echo = FALSE--------------------------------------
DT::datatable(output_up,rownames = FALSE)

## ----example upper plot, class.source="scroll-100"----------------------------
# generate some concentration for the fitting curve 
logc_plot <- seq(from=-3,to=2,by=0.05)
conc_plot <- 10**logc_plot

# initiate plot
plot(conc,resp,xlab="conc (uM)",ylab="Response",xlim=c(0.001,500),ylim=c(-5,40),
       log="x",main=paste(name,"\n",assay),cex.main=0.9)
# add vertical lines to mark the maximum concentration in the data and the upper boundary set by bmd_up_bnd
abline(v=max(conc), lty = 1, col = "brown", lwd=2)
abline(v=160, lty = 2, col = "darkviolet", lwd=2)

# add marker for BMD and its boundaries before `bounding`
lines(c(res_up$bmd,res_up$bmd),c(0,50),col="green",lwd=2)
rect(xleft=res_up$bmdl,ybottom=0,xright=res_up$bmdu,ytop=50,col=rgb(0,1,0, alpha = .5), border = NA)
points(res_up$bmd, 0, pch = "x", col = "green")

# add marker for BMD and its boundaries after `bounding`
lines(c(res_up2$bmd,res_up2$bmd),c(0,50),col="blue",lwd=2)
rect(xleft=res_up2$bmdl,ybottom=0,xright=res_up2$bmdu,ytop=50,col=rgb(0,0,1, alpha = .5), border = NA)
points(res_up2$bmd, 0, pch = "x", col = "blue")

# add the fitting curve
lines(conc_plot, poly1(ps = c(res_up$a), conc_plot))
legend(1e-3, 40, legend=c("Maximum Dose Tested", "Boundary", "BMD-before", "BMD-after"),
       col=c("brown", "darkviolet", "green", "blue"), lty=c(1,2,1,1))

## ----example with hit core----------------------------------------------------
# using the same data, fit curves 
param <- tcplfit2_core(conc2, resp2, cutoff = cutoff)
hit_res <- tcplhit2_core(param, conc2, resp2, cutoff = cutoff, onesd = onesd, 
                         bmd_low_bnd = 0.8)

## ----res-hit core-------------------------------------------------------------
# adding the result from tcplhit2_core to the output table for comparison
hit_res["Name"]<-  paste("Chlorothalonil", "tcplhit2_core", sep = "-")
hit_res['Min. Conc.'] <- min(conc2)
hit_res[1, c("Min. Conc.", "bmd", "bmdl", "bmdu")] <- round(hit_res[1, c("Min. Conc.", "bmd", "bmdl", "bmdu")], 3)

output_low <- rbind(output_low, 
                    hit_res[1, c('Name', "Min. Conc.", "bmd", "bmdl", "bmdu")])

## ----res-hit_table, echo = FALSE----------------------------------------------
DT::datatable(output_low,rownames = FALSE)

## ----example even lower bound-------------------------------------------------
res_low3 <- concRespCore(row_low,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", 
                                           "pow", "exp2", "exp3", "exp4", "exp5"), 
                         conthits = T, aicc = F, bidirectional=F, bmd_low_bnd = 0.4)

## ----example even lower bound-res---------------------------------------------
# print out the new results
# add to previous results for comparison 
res_low3['Min. Conc.'] <- min(conc2)
res_low3['Name'] <- paste("Chlorothalonil", "after `bounding` (two fifths)", sep = "-")
res_low3[1, c("Min. Conc.", "bmd", "bmdl", "bmdu")] <- round(res_low3[1, c("Min. Conc.", "bmd", "bmdl", "bmdu")], 3)

output_low <- rbind(output_low[-3, ], 
                    res_low3[1, c('Name', "Min. Conc.", "bmd", "bmdl", "bmdu")])

## ----lower_bound_res_table, echo = FALSE--------------------------------------
DT::datatable(output_low,rownames = FALSE)

## ----example even lower bound-plot, class.source="scroll-100"-----------------
# initiate the plot
plot(conc2,resp2,xlab="conc (uM)",ylab="Response",xlim=c(0.001,100),ylim=c(-5,60),
       log="x",main=paste(name,"\n",assay),cex.main=0.9)

# add vertical lines to mark the minimum concentration in the data and the lower boundary set by bmd_low_bnd
abline(v=min(conc2), lty = 1, col = "brown", lwd = 2)
abline(v=0.4*min(conc2), lty = 2, col = "darkviolet", lwd = 2)

# add markers for BMD and its boundaries before `bounding`
lines(c(res_low$bmd,res_low$bmd),c(0,50),col="green",lwd=2)
rect(xleft=res_low$bmdl,ybottom=0,xright=res_low$bmdu,ytop=50,col=rgb(0,1,0, alpha = .5), border = NA)
points(res_low$bmd, 0, pch = "x", col = "green")

# add markers for BMD and its boundaries after `bounding`
lines(c(res_low3$bmd,res_low3$bmd),c(0,50),col="blue",lwd=2)
rect(xleft=res_low3$bmdl,ybottom=0,xright=res_low3$bmdu,ytop=50,col=rgb(0,0,1, alpha = .5), border = NA)
points(res_low3$bmd, 0, pch = "x", col = "blue")

# add the fitted curve
lines(conc_plot, exp4(ps = c(res_low$tp, res_low$ga), conc_plot))
legend(1e-3, 60, legend=c("Lowest Dose Tested", "Boundary Dose", "BMD-before", "BMD-after"),
       col=c("brown", "darkviolet", "green", "blue"), lty=c(1,2,1,1))

## ----appendix plt1, fig.height = 6, fig.width = 7, warning = FALSE------------
  # call additional R packages
  library(stringr)  # string management package

  # read in the file
  data("signatures")
  
  # set up a 3 x 2 grid for the plots
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
  par(mfrow=c(3,2),mar=c(4,4,2,2))
    
  # fit 6 observations in signatures
  for(i in 1:nrow(signatures)){
    # set up input data
    row = list(conc=as.numeric(str_split(signatures[i,"conc"],"\\|")[[1]]),
               resp=as.numeric(str_split(signatures[i,"resp"],"\\|")[[1]]),
               bmed=0,
               cutoff=signatures[i,"cutoff"],
               onesd=signatures[i,"onesd"],
               name=signatures[i,"name"],
               assay=signatures[i,"signature"])
    # run concentration-response modeling (1st plotting option)
    out = concRespCore(row,conthits=F,do.plot=T)
    if(i==1){
      res <- out
    }else{
      res <- rbind.data.frame(res,out)
    }
  }

## ----appendix plt2, fig.height = 8, fig.width = 7-----------------------------
  # set up a 3 x 2 grid for the plots
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
  par(mfrow=c(3,2),mar=c(4,4,2,2))
  # plot results using `concRespPlot`
  for(i in 1:nrow(res)){
    concRespPlot(res[i,],ymin=-1,ymax=1)
  }

## -----------------------------------------------------------------------------
# Load the example data set
data("signatures")
# using the first row of signatures data as an example 
signatures[1,]

## -----------------------------------------------------------------------------
# using the first row of signature as an example 
conc <- as.numeric(str_split(signatures[1,"conc"],"\\|")[[1]])
resp <- as.numeric(str_split(signatures[1,"resp"],"\\|")[[1]])
cutoff <- signatures[1,"cutoff"]

# run curve fitting
output <- tcplfit2_core(conc, resp, cutoff)
# show the structure of the output 
summary(output)

## -----------------------------------------------------------------------------
# get plots in normal and in log-10 concentration scale
basic <- plot_allcurves(output, conc, resp)
basic_log <- plot_allcurves(output, conc, resp, log_conc = T)
grid.arrange(basic, basic_log)

## -----------------------------------------------------------------------------
# prepare the 'row' object for concRespCore
row <- list(conc=conc,
           resp=resp,
           bmed=0,
           cutoff=cutoff,
           onesd=signatures[1,"onesd"],
           name=signatures[1,"name"],
           assay=signatures[1,"signature"])

# run concentration-response modeling 
out <-  concRespCore(row,conthits=F)
# show the output
out

## -----------------------------------------------------------------------------
# pass the output to the plotting function
basic_plot <- concRespPlot2(out)
basic_log <- concRespPlot2(out, log_conc = TRUE)
res <- grid.arrange(basic_plot, basic_log)

## -----------------------------------------------------------------------------
# Using the fitted result and plot from the example in the last section
# get the cutoff from the output
cutoff <- out[, "cutoff"]

basic_plot + 
  # Cutoff Band - a transparent rectangle
  geom_rect(aes(xmin = 0,xmax = 30,ymin = -cutoff,ymax = cutoff),
            alpha = 0.1,fill = "skyblue") +
  # Titles
  ggtitle(
    label = paste("Best Model Fit",
                  out[, "name"],
                  sep = "\n"),
    subtitle = paste("Assay Endpoint: ",
                     out[, "assay"])) +
  ## Add BMD and BMR labels
  geom_hline(
    aes(yintercept = out[, "bmr"]),
    col = "blue") +
  geom_segment(
    aes(x = out[, "bmd"], xend = out[, "bmd"], y = -0.5, yend = out[, "bmr"]),
    col = "blue"
  ) + geom_point(aes(x = out[, "bmd"], y = out[, "bmr"], fill = "BMD"), shape = 21, cex = 2.5)

## -----------------------------------------------------------------------------
# Get all potency estimates and the corresponding y value on the curve
estimate_points <- out %>%
  select(bmd, acc, ac50, ac10, ac5) %>%
  tidyr::pivot_longer(everything(), names_to = "Potency Estimates") %>%
  mutate(`Potency Estimates` = toupper(`Potency Estimates`)) 

y <-  c(out[, "bmr"], out[, "cutoff"], rep(out[, "top"], 3))
y <-  y * c(1, 1, .5, .1, .05)
estimate_points <- cbind(estimate_points, y = y)

# add Potency Estimate Points and set colors
basic_plot + geom_point(
  data = estimate_points,
  aes(x = value, y = y, fill = `Potency Estimates`), shape = 21, cex = 2.5
)

## -----------------------------------------------------------------------------
# add Potency Estimate Points and set colors - with plot in log-10 concentration
basic_log + geom_point(
  data = estimate_points,
  aes(x = log10(value), y = y, fill = `Potency Estimates`), shape = 21, cex = 2.5
)

## -----------------------------------------------------------------------------
# maybe want to extract and use the same x's in the base plot 
# to calculate predicted responses 
conc_plot <- basic_plot[["layers"]][[2]][["data"]][["conc_plot"]]

basic_plot +
  # fitted parameter values of another curve you want to add
  geom_line(data=data.frame(x=conc_plot, y=tcplfit2::exp5(c(0.5, 10, 1.2), conc_plot)), aes(x,y,color = "exp5"))+
  # add different colors for comparisons 
  scale_colour_manual(values=c("#CC6666", "#9999CC"),
                      labels = c("Curve 1-exp4", "Curve 2-exp5")) +
  labs(title = "Curve 1 v.s. Curve 2")

## ----example 1----------------------------------------------------------------
# some example data
conc <- list(.03, .1, .3, 1, 3, 10, 30, 100)
resp <- list(0, .2, .1, .4, .7, .9, .6, 1.2)
row <- list(conc = conc,
            resp = resp,
            bmed = 0,
            cutoff = 1,
            onesd = .5)

# AUC is included in the output
concRespCore(row, conthits = TRUE, AUC = TRUE)

## ----example 2, fig.height = 4.55, fig.width = 8------------------------------
# This is taken from the example under tcplfit2_core
conc_ex2 <- c(.03, .1, .3, 1, 3, 10, 30, 100)
resp_ex2 <- c(0, .1, 0, .2, .6, .9, 1.1, 1)

# fit all available models in the package
# show all fitted curves 
output_ex2 <- tcplfit2_core(conc_ex2, resp_ex2, .8)
grid.arrange(plot_allcurves(output_ex2, conc_ex2, resp_ex2),
          plot_allcurves(output_ex2, conc_ex2, resp_ex2, log_conc = TRUE), ncol = 2)

## ----example 2 cont., fig.height = 6, fig.width = 6---------------------------
fit_method <- "hill"
# extract the parameters 
modpars <- output_ex2[[fit_method]][output_ex2[[fit_method]]$pars]

# plug into get_AUC function 
estimated_auc1 <- get_AUC(fit_method, min(conc_ex2), max(conc_ex2), modpars)
estimated_auc1

# extract the predicted responses from the model
pred_resp <- output_ex2[[fit_method]][["modl"]]

# plot to see if the result make sense
# the shaded area is what the function tries to find
plot(conc_ex2, pred_resp)
lines(conc_ex2, pred_resp)
polygon(c(conc_ex2, max(conc_ex2)), c(pred_resp, min(pred_resp)), col=rgb(1, 0, 0,0.5))

## ----example 2 other models---------------------------------------------------
# list of models
fitmodels <- c("gnls", "poly1", "poly2", "pow", "exp2", "exp3", "exp4", "exp5")
mylist <- list()
for (model in fitmodels){

  fit_method <- model
  # extract corresponding model parameters
  modpars <- output_ex2[[fit_method]][output_ex2[[fit_method]]$pars]
  
  # get AUC
  mylist[[fit_method]] <- get_AUC(fit_method, min(conc_ex2), max(conc_ex2), modpars)
  
}
# print AUC's for other models 
data.frame(mylist,row.names = "AUC")

## ----example 3, fig.height = 4.55, fig.width = 8------------------------------
# Taking the code from example 3 in the vignette 
library(stringr)  # string management package
data("signatures")

# use row 5 in the data
conc <- as.numeric(str_split(signatures[5,"conc"],"\\|")[[1]])
resp <- as.numeric(str_split(signatures[5,"resp"],"\\|")[[1]])
cutoff <- signatures[5,"cutoff"]

# plot all models, this is an example of negative curves 
output_negative <- tcplfit2_core(conc, resp, cutoff)
grid.arrange(plot_allcurves(output_negative, conc, resp),
          plot_allcurves(output_negative, conc, resp, log_conc = TRUE), ncol = 2)

## ----example 3 cont., fig.height = 6, fig.width = 6---------------------------
fit_method <- "exp3"

# extract corresponding model parameters and predicted response
modpars <- output_negative[[fit_method]][output_negative[[fit_method]]$pars]
pred_resp <- output_negative[[fit_method]][["modl"]]

estimated_auc2 <- get_AUC(fit_method, min(conc), max(conc), modpars)
estimated_auc2

# plot this curve
pred_resp <- pred_resp[order(conc)]
plot(conc[order(conc)], pred_resp)
lines(conc[order(conc)], pred_resp)
polygon(c(conc[order(conc)], max(conc)), c(pred_resp, max(pred_resp)), col=rgb(1, 0, 0,0.5))

## ----example 3 convert negative AUC-------------------------------------------
get_AUC(fit_method, min(conc), max(conc), modpars, return.abs = TRUE) 

## ----example 4, fig.height = 6, fig.width = 6---------------------------------
# simulate a poly2 curve
conc_sim <- seq(0,3, length.out = 100)
## biphasic poly2 parameters
b1 <- -1.3
b2 <- 0.7
## converted to tcplfit2's poly2 parameters
a <- b1^2/b2
b <- b1/b2

## plot the curve
resp_sim <- poly2(c(a, b, 0.1), conc_sim)
plot(conc_sim, resp_sim, type = "l")
abline(h = 0)

## ----example 4 cont.----------------------------------------------------------
# get AUC for the simulated Polynomial 2 curve 
get_AUC("poly2", min(conc_sim), max(conc_sim), ps = c(a, b))

## ----example 5----------------------------------------------------------------
out <- tcplhit2_core(output_ex2, conc_ex2, resp_ex2, 0.8, onesd = 0.4)
out
post_hit_AUC(out)

## ----setup-2, warning=FALSE---------------------------------------------------
# prepare concentration data for demonstration
ex_conc <- seq(0, 100, length.out = 500)
ex2_conc <- seq(0, 3, length.out = 100)

## ----poly 1, fig.width=5, fig.height=5, warning=FALSE-------------------------
poly1_plot <- ggplot(mapping=aes(ex_conc)) +  
  geom_line(aes(y = 55*ex_conc, color = "a=55")) +
  geom_line(aes(y = 10*ex_conc, color = "a=10")) +
  geom_line(aes(y = 0.05*ex_conc, color = "a=0.05")) +
  geom_line(aes(y = -5*ex_conc, color = "a=(-5)")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.1,0.8)) +
  scale_color_manual(name='a values',
                     breaks=c('a=(-5)', 'a=0.05', 'a=10', 'a=55'),
                     values=c('a=(-5)'='black', 'a=0.05' = 'red', 'a=10'='blue', 'a=55'='darkviolet'))

poly1_plot

## ----poly 2, fig.width=8, fig.height=5, warning=FALSE-------------------------
fits_poly <- data.frame(
  # change a 
  y1 = poly2(ps = c(a = 40, b = 2),x = ex_conc),
  y2 = poly2(ps = c(a = 6, b = 2),x = ex_conc),
  y3 = poly2(ps = c(a = 0.1, b = 2),x = ex_conc),
  y4 = poly2(ps = c(a = -2, b = 2),x = ex_conc),
  y5 = poly2(ps = c(a = -20, b = 2),x = ex_conc),
  
  # change b 
  y6 = poly2(ps = c(a = 4,b = 1.8),x = ex_conc),
  y7 = poly2(ps = c(a = 4,b = 7),x = ex_conc),
  y8 = poly2(ps = c(a = 4,b = 16),x = ex_conc)
)

# shows how changes in parameter 'a' affect the shape of the curve 
poly2_plot1 <- ggplot(fits_poly, aes(ex_conc)) +
  geom_line(aes(y = y1, color = "a=40")) +
  geom_line(aes(y = y2, color = "a=6")) +
  geom_line(aes(y = y3, color = "a=0.1")) +
  geom_line(aes(y = y4, color = "a=(-2)")) +
  geom_line(aes(y = y5, color = "a=(-20)")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.15,0.8)) +
  scale_color_manual(name='a values',
                     breaks=c('a=(-20)', 'a=(-2)', 'a=0.1', 'a=6', 'a=40'),
                     values=c('a=(-20)'='black', 'a=(-2)'='red', 'a=0.1'='blue', 'a=6'='darkviolet', 'a=40'='darkgoldenrod1'))

# shows how changes in parameter 'b' affect the shape of the curve 
poly2_plot2 <- ggplot(fits_poly, aes(ex_conc)) +  
  geom_line(aes(y = y6, color = "b=1.8")) +
  geom_line(aes(y = y7, color = "b=7")) +
  geom_line(aes(y = y8, color = "b=16")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.15,0.8)) +
  scale_color_manual(name='b values',
                     breaks=c('b=1.8', 'b=7', 'b=16'),
                     values=c('b=1.8'='black', 'b=7'='red', 'b=16'='blue'))

grid.arrange(poly2_plot1, poly2_plot2, ncol = 2)

## ----pow, fig.width=8, fig.height=5, warning=FALSE----------------------------
fits_pow <- data.frame(
  # change a
  y1 = pow(ps = c(a = 0.48,p = 1.45),x = ex2_conc),
  y2 = pow(ps = c(a = 7.2,p = 1.45),x = ex2_conc),
  y3 = pow(ps = c(a = -3.2,p = 1.45),x = ex2_conc),
  
  # change p
  y4 = pow(ps = c(a = 1.2,p = 0.3),x = ex2_conc),
  y5 = pow(ps = c(a = 1.2,p = 1.6),x = ex2_conc),
  y6 = pow(ps = c(a = 1.2,p = 3.2),x = ex2_conc)
)

# shows how changes in parameter 'a' affect the shape of the curve
pow_plot1 <- ggplot(fits_pow, aes(ex2_conc)) +  
  geom_line(aes(y = y1, color = "a=0.48")) +
  geom_line(aes(y = y2, color = "a=7.2")) +
  geom_line(aes(y = y3, color = "a=(-3.2)")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.15,0.8)) +
  scale_color_manual(name='a values',
                     breaks=c('a=(-3.2)', 'a=0.48', 'a=7.2'),
                     values=c('a=(-3.2)'='black', 'a=0.48'='red', 'a=7.2'='blue'))

# shows how changes in parameter 'p' affect the shape of the curve
pow_plot2 <- ggplot(fits_pow, aes(ex2_conc)) +  
  geom_line(aes(y = y4, color = "p=0.3")) +
  geom_line(aes(y = y5, color = "p=1.6")) +
  geom_line(aes(y = y6, color = "p=3.2")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.15,0.8)) +
  scale_color_manual(name='p values',
                     breaks=c('p=0.3', 'p=1.6', 'p=3.2'),
                     values=c('p=0.3'='black', 'p=1.6'='red', 'p=3.2'='blue'))



grid.arrange(pow_plot1, pow_plot2, ncol = 2)

## ----Hill, fig.height=5, fig.width=8, warning=FALSE---------------------------
fits_hill <- data.frame(
  # change tp
  y1 = hillfn(ps = c(tp = -200,ga = 5,p = 1.76), x = ex_conc),
  y2 = hillfn(ps = c(tp = 200,ga = 5,p = 1.76), x = ex_conc),
  y3 = hillfn(ps = c(tp = 850,ga = 5,p = 1.76), x = ex_conc),

  # change ga
  y4 = hillfn(ps = c(tp = 120,ga = 4,p = 1.76), x = ex_conc),
  y5 = hillfn(ps = c(tp = 120,ga = 12,p = 1.76), x = ex_conc),
  y6 = hillfn(ps = c(tp = 120,ga = 20,p = 1.76), x = ex_conc),
  
  # change p
  y7 = hillfn(ps = c(tp = 120,ga = 5,p = 0.5), x = ex_conc),
  y8 = hillfn(ps = c(tp = 120,ga = 5,p = 2), x = ex_conc),
  y9 = hillfn(ps = c(tp = 120,ga = 5,p = 5), x = ex_conc)
  
)

# shows how changes in parameter 'tp' affect the shape of the curve
hill_plot1 <- ggplot(fits_hill, aes(log10(ex_conc))) +  
  geom_line(aes(y = y1, color = "tp=(-200)")) +
  geom_line(aes(y = y2, color = "tp=200")) +
  geom_line(aes(y = y3, color = "tp=850")) +
  labs(x = "Concentration in Log-10 Scale", y = "Response") +
  theme(legend.position = c(0.2,0.8),
        legend.key.size = unit(0.5, 'cm')) +
  scale_color_manual(name='tp values',
                     breaks=c('tp=(-200)', 'tp=200', 'tp=850'),
                     values=c('tp=(-200)'='black', 'tp=200'='red', 'tp=850'='blue'))

# shows how changes in parameter 'ga' affect the shape of the curve
hill_plot2 <- ggplot(fits_hill, aes(log10(ex_conc))) + 
  geom_line(aes(y = y4, color = "ga=4")) +
  geom_line(aes(y = y5, color = "ga=12")) +
  geom_line(aes(y = y6, color = "ga=20")) +
  labs(x = "Concentration in Log-10 Scale", y = "Response") +
  theme(legend.position = c(0.8,0.25),
        legend.key.size = unit(0.4, 'cm')) +
  scale_color_manual(name='ga values',
                     breaks=c('ga=4', 'ga=12', 'ga=20'),
                     values=c('ga=4'='black', 'ga=12'='red', 'ga=20'='blue'))

# shows how changes in parameter 'p' affect the shape of the curve
hill_plot3 <- ggplot(fits_hill, aes(log10(ex_conc))) +  
  geom_line(aes(y = y7, color = "p=0.5")) +
  geom_line(aes(y = y8, color = "p=2")) +
  geom_line(aes(y = y9, color = "p=5")) +
  labs(x = "Concentration in Log-10 Scale", y = "Response") +
  theme(legend.position = c(0.8,0.2),
        legend.key.size = unit(0.4, 'cm')) +
  scale_color_manual(name='p values',
                     breaks=c('p=0.5', 'p=2', 'p=5'),
                     values=c('p=0.5'='black', 'p=2'='red', 'p=5'='blue'))


grid.arrange(hill_plot1, hill_plot2, hill_plot3, ncol = 2, nrow = 2)

## ----gnls, fig.width=8, fig.height=5, warning=FALSE---------------------------
fits_gnls <- data.frame(
  # change la
  y1 = gnls(ps = c(tp = 750,ga = 15,p = 1.45,la = 17,q = 1.34), x = ex_conc),
  y2 = gnls(ps = c(tp = 750,ga = 15,p = 1.45,la = 50,q = 1.34), x = ex_conc),
  y3 = gnls(ps = c(tp = 750,ga = 15,p = 1.45,la = 100,q = 1.34), x = ex_conc),

  # change q
  y4 = gnls(ps = c(tp = 750,ga = 15,p = 1.45,la = 20,q = 0.3), x = ex_conc),
  y5 = gnls(ps = c(tp = 750,ga = 15,p = 1.45,la = 20,q = 1.2), x = ex_conc),
  y6 = gnls(ps = c(tp = 750,ga = 15,p = 1.45,la = 20,q = 8), x = ex_conc)
  
)

# shows how changes in parameter 'la' affect the shape of the curve
gnls_plot1 <- ggplot(fits_gnls, aes(log10(ex_conc))) +  
  geom_line(aes(y = y1, color = "la=17")) +
  geom_line(aes(y = y2, color = "la=50")) +
  geom_line(aes(y = y3, color = "la=100")) +
  labs(x = "Concentration in Log-10 Scale", y = "Response") +
  theme(legend.position = c(0.15,0.8)) +
  scale_color_manual(name='la values',
                     breaks=c('la=17', 'la=50', 'la=100'),
                     values=c('la=17'='black', 'la=50'='red', 'la=100'='blue'))

# shows how changes in parameter 'q' affect the shape of the curve
gnls_plot2 <- ggplot(fits_gnls, aes(log10(ex_conc))) +  
  geom_line(aes(y = y4, color = "q=0.3")) +
  geom_line(aes(y = y5, color = "q=1.2")) +
  geom_line(aes(y = y6, color = "q=8")) +
  labs(x = "Concentration in Log-10 Scale", y = "Response") +
  theme(legend.position = c(0.15,0.8)) +
  scale_color_manual(name='q values',
                     breaks=c('q=0.3', 'q=1.2', 'q=8'),
                     values=c('q=0.3'='black', 'q=1.2'='red', 'q=8'='blue'))
  
  
grid.arrange(gnls_plot1, gnls_plot2, ncol = 2)

## ----exp2, fig.width=8, fig.height=5, warning=FALSE---------------------------
fits_exp2 <- data.frame(
  # change a
  y1 = exp2(ps = c(a = 20,b = 12), x = ex2_conc),
  y2 = exp2(ps = c(a = 9,b = 12), x = ex2_conc),
  y3 = exp2(ps = c(a = 0.1,b = 12), x = ex2_conc),
  y4 = exp2(ps = c(a = -3,b = 12), x = ex2_conc),
  
  # change b
  y5 = exp2(ps = c(a = 0.45,b = 4), x = ex2_conc),
  y6 = exp2(ps = c(a = 0.45,b = 9), x = ex2_conc),
  y7 = exp2(ps = c(a = 0.45,b = 20), x = ex2_conc)
  
)

# shows how changes in parameter 'a' affect the shape of the curve 
exp2_plot1 <- ggplot(fits_exp2, aes(ex2_conc)) +  
  geom_line(aes(y = y1, color = "a=20")) +
  geom_line(aes(y = y2, color = "a=9")) +
  geom_line(aes(y = y3, color = "a=0.1")) +
  geom_line(aes(y = y4, color = "a=(-3)")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.15,0.8)) +
  scale_color_manual(name='a values',
                     breaks=c('a=(-3)', 'a=0.1', 'a=9', 'a=20'),
                     values=c('a=(-3)'='black', 'a=0.1'='red', 'a=9'='blue', 'a=20'='darkviolet'))

# shows how changes in parameter 'b' affect the shape of the curve 
exp2_plot2 <- ggplot(fits_exp2, aes(ex2_conc)) +  
  geom_line(aes(y = y5, color = "b=4")) +
  geom_line(aes(y = y6, color = "b=9")) +
  geom_line(aes(y = y7, color = "b=20")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.15,0.8)) +
  scale_color_manual(name='b values',
                     breaks=c('b=4', 'b=9', 'b=20'),
                     values=c('b=4'='black', 'b=9'='red', 'b=20'='blue'))

grid.arrange(exp2_plot1, exp2_plot2, ncol = 2)

## ----exp3, fig.width=5, fig.height=5, warning=FALSE---------------------------
fits_exp3 <- data.frame(
  # change p
  y1 = exp3(ps = c(a = 1.67,b = 12.5,p = 0.3), x = ex2_conc),
  y2 = exp3(ps = c(a = 1.67,b = 12.5,p = 0.9), x = ex2_conc),
  y3 = exp3(ps = c(a = 1.67,b = 12.5,p = 1.2), x = ex2_conc)
  
)

# shows how changes in parameter 'p' affect the shape of the curve 
exp3_plot <- ggplot(fits_exp3, aes(ex2_conc)) +  
  geom_line(aes(y = y1, color = "p=0.3")) +
  geom_line(aes(y = y2, color = "p=0.9")) +
  geom_line(aes(y = y3, color = "p=1.2")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.15,0.8)) +
  scale_color_manual(name='p values',
                     breaks=c('p=0.3', 'p=0.9', 'p=1.2'),
                     values=c('p=0.3'='black', 'p=0.9'='red', 'p=1.2'='blue'))


exp3_plot

## ----exp4, fig.width=8, fig.height=5, warning=FALSE---------------------------
fits_exp4 <- data.frame(
  # change tp  
  y1 = exp4(ps = c(tp = 895,ga = 15),x = ex_conc),
  y2 = exp4(ps = c(tp = 200,ga = 15),x = ex_conc),
  y3 = exp4(ps = c(tp = -500,ga = 15),x = ex_conc),
  
  # change ga
  y4 = exp4(ps = c(tp = 500,ga = 0.4),x = ex_conc),
  y5 = exp4(ps = c(tp = 500,ga = 10),x = ex_conc),
  y6 = exp4(ps = c(tp = 500,ga = 20),x = ex_conc)
  
)

# shows how changes in parameter 'tp' affect the shape of the curve 
exp4_plot1 <- ggplot(fits_exp4, aes(ex_conc)) +  
  geom_line(aes(y = y1, color = "tp=895")) +
  geom_line(aes(y = y2, color = "tp=200")) +
  geom_line(aes(y = y3, color = "tp=(-500)")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.1,0.8)) +
  scale_color_manual(name='tp values',
                     breaks=c('tp=(-500)', 'tp=200', 'tp=895'),
                     values=c('tp=(-500)'='black', 'tp=200'='red', 'tp=895'='blue'))


# shows how changes in parameter 'ga' affect the shape of the curve 
exp4_plot2 <- ggplot(fits_exp4, aes(ex_conc)) +  
  geom_line(aes(y = y4, color = "ga=0.4")) +
  geom_line(aes(y = y5, color = "ga=10")) +
  geom_line(aes(y = y6, color = "ga=20")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.8,0.2)) +
  scale_color_manual(name='ga values',
                     breaks=c('ga=0.4', 'ga=10', 'ga=20'),
                     values=c('ga=0.4'='black', 'ga=10'='red', 'ga=20'='blue'))


grid.arrange(exp4_plot1, exp4_plot2, ncol = 2)

## ----exp5, fig.width=5, fig.height=5, warning=FALSE---------------------------
fits_exp5 <- data.frame(
  # change p
  y1 = exp5(ps = c(tp = 793,ga = 6.25,p = 0.3), x = ex_conc),
  y2 = exp5(ps = c(tp = 793,ga = 6.25,p = 3.4), x = ex_conc),
  y3 = exp5(ps = c(tp = 793,ga = 6.25,p = 8), x = ex_conc)
  
)

# shows how changes in parameter 'p' affect the shape of the curve 
exp5_plot <- ggplot(fits_exp5, aes(ex_conc)) +  
  geom_line(aes(y = y1, color = "p=0.3")) +
  geom_line(aes(y = y2, color = "p=3.4")) +
  geom_line(aes(y = y3, color = "p=8")) +
  labs(x = "Concentration", y = "Response") +
  theme(legend.position = c(0.8,0.2)) +
  scale_color_manual(name='p values',
                     breaks=c('p=0.3', 'p=3.4', 'p=8'),
                     values=c('p=0.3'='black', 'p=3.4'='red', 'p=8'='blue'))


exp5_plot

## ----appendix-table, echo=FALSE, warning=FALSE--------------------------------
# First column - tcplfit2 available models.
Model <- c(
  "Constant", "Linear", "Quadratic","Power", "Hill", "Gain-Loss",
  "Exponential 2", "Exponential 3","Exponential 4", "Exponential 5"
)
# Second column - model abbreviations used in invitrodb & tcplfit2.
Abbreviation <- c(
  "cnst", "poly1", "poly2","pow", "hill", "gnls",
  "exp2", "exp3", "exp4", "exp5"
)
# Third column - model equations.
Equations <- c(
  "$f(x) = 0$", # constant
  "$f(x) = ax$", # linear
  "$f(x) = a(\\frac{x}{b}+(\\frac{x}{b})^{2})$", # quadratic
  "$f(x) = ax^p$", # power
  "$f(x) = \\frac{tp}{1 + (\\frac{ga}{x})^{p}}$", # hill
  "$f(x) = \\frac{tp}{(1 + (\\frac{ga}{x})^{p} )(1 + (\\frac{x}{la})^{q} )}$", # gain-loss
  "$f(x) = a*(exp(\\frac{x}{b}) - 1)$", # exp 2
  "$f(x) = a*(exp((\\frac{x}{b})^{p}) - 1)$", # exp 3
  "$f(x) = tp*(1-2^{\\frac{-x}{ga}})$", # exp 4
  "$f(x) = tp*(1-2^{-(\\frac{x}{ga})^{p}})$" # exp 5
)
# Fourth column - model parameter descriptions.
OutputParameters <- c(
  "", # constant
  "a (y-scale)", # linear,
  "a (y-scale) </br> b (x-scale)", # quadratic
  "a (y-scale) </br> p (power)", # power
  "tp (top) </br> ga (gain AC50) </br> p (gain-power)", # hill
  "tp (top) </br> ga (gain AC50) </br> p (gain power) </br> la (loss AC50) </br> q (loss power)", # gain-loss
  "a (y-scale) </br> b (x-scale)", # exp2
  "a (y-scale) </br> b (x-scale) </br> p (power)", # exp3
  "tp (top) </br> ga (AC50)", # exp4
  "tp (top) </br> ga (AC50) </br> p (power)" # exp5
)
# Fifth column - additional model details.
Details <- c(
  "Parameters always equals 'er'.", # constant
  "", # linear 
  "", # quadratic
  "", # power
  "Concentrations are converted internally to log10 units and optimized with f(x) = tp/(1 + 10^(p*(gax))), then ga and ga_sd are converted back to regular units before returning.", # hill
  "Concentrations are converted internally to log10 units and optimized with f(x) = tp/[(1 + 10^(p*(gax)))(1 + 10^(q*(x-la)))], then ga, la, ga_sd, and la_sd are converted back to regular units before returning." , # gain-loss
  "", # exp2
  "", # exp3
  "", # exp4
  "") # exp5
# Consolidate all columns into a table.
output <- 
  data.frame(Model, Abbreviation, Equations,
             OutputParameters, Details)
# Export/print the table into an html rendered table.
htmlTable(output,
        align = 'l',
        align.header = 'l',
        rnames = FALSE  ,
        css.cell =  ' padding-bottom: 5px;  vertical-align:top; padding-right: 10px;min-width: 5em ',
        caption="*tcplfit2* model details.",
        tfoot = "Model descriptions are pulled from tcplFit2 manual at <https://cran.r-project.org/package=tcplfit2/tcplfit2.pdf>."
)

