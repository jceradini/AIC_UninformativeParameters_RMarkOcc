
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AIC and uninformative parameters functions
# Joe Ceradini, University of Wyoming, Utah Valley University (Capitol Reef Field Station)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# To streamline and improve comparison of nested (and non-nested) models when using AIC.

# WORK IN PROGRESS. 
# Only works with single-season occupancy models from RMark. 
# - Double check output, especially for models with interactions. 
# - Would be great to have a much improved version of this that worked with multiple R model types.

# NESTED: Both models must have identical terms and 
#         one of the models must have one or more extra terms.
# www.statisticshowto.com/nested-model/

# Refs: 
#   Arnold. 2010. Uninformative parameters and model selection using Akaike's information criterion. Journal of Wildlife Management 74:1175-1178.
#   Murtaugh. 2014. In defense of P values. Ecology 95:611-617.
#   Cooch and White. 2016. Ch 4: Building & comparing models. Sidebar from 4-61 - 4-62. Program MARK - A Gentle introduction.
#   Burnham and Anderson. 2002. Model selection and multimodel inference. 2nd edition.

# Model set can be a mix of nested and nonnested models.
# For models that differ by a single predictor, function will determine
#   if the additional predictor in the larger model is informative based
#   on the chosen AICc threshold. The AICc weights will then be recalculated
#   with only the models determined to be informative. A LRT is also calculated. 
# Models that differ by a single predictor can differ by 1 or more
#   parameters. For example, if the larger model contains a factor with >2
#   levels.

# aic.thresh is the chosen threshold that larger model needs to overcome to
#   be considered informative. aic.thresh is multiplied by the difference in 
#   degrees of freedom between the nested models being compared. So, if
#   models differ by 1 DF then aic.thresh stays the same (since it is multiplied by 1),
#   but if the DF difference is >1 then aic.thresh is larger. 
# However, if aic.thresh = 0 then the thresh will always be 0 regardless of 
#   DF diff btw models being compared.  

# Uninformative based on AICc vs LRT
# - Extent to which these agree depends on,
#   1) sample size, since the AICc penalty per additional parameter 
#     increases as sample size decreases.
#   2) how many parameters are being added by the additional
#     predictor (i.e., is it continuous or categorical/factor)
# - See: 
#   Cooch and White. 2016. Ch 4: Building & comparing models. Sidebar from 4-61 - 4-62.
#     Program MARK - A Gentle introduction.
#   Murtaugh. 2014. In defense of P values. Ecology 95:611-617. Especially Fig. 2.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pr(chi-sq_df >= (2 x df)) = the P-value of the likelihood ratio test when 
# the AIC values of the 2 models are equal.
# - Differ by 1 df, so 1*2 = 2
(1 - pchisq(2, df = 1))

# - Differ by 2 df, so 2*2 = 5
(1 - pchisq(4, df = 2))

# - Differ by 3 df, so 3*2 = 6
(1 - pchisq(6, df = 3))

# - Differ by 4 df, so 4*2 = 8
(1 - pchisq(8, df = 4))
(1 - pchisq(5.87361, df = 4))


require(RMark)
#require(lmtest)
require(AICcmodavg)
require(tidyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# aic.nested.occup ####

# - currently only test that has worked is single-season occupancy models
# - function is expecting the list of models that is returned from the RMark 
#   function approach.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Capture history (CH) file
# setwd("C:/Users/jceradin/Dropbox (UW WYNDD)/Proj_WGFD_WTPDogInventory/Data/rmark_input_files")
# ch <- read.table("ch_noNA.txt", sep = ",", header = TRUE)
# str(ch) # 1 row per site = 440 rows
# ch$ch <- as.character(ch$ch)
# unique(nchar(ch$ch)) # good - all have 3 occasions
# 
# # Process data
# wtpd.proc <- process.data(ch, model = "Occupancy", groups = c("stratum"), begin.time=1)
# str(wtpd.proc)
# 
# # Design data
# wtpd.ddl <- make.design.data(wtpd.proc)
# str(wtpd.ddl)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run models via function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setwd("C:/Users/jceradin/Box Sync/Research/Analysis/Mark-Recap/test")
# combos.s1 <- function() # START FUNCTION
# {
#   # Detection: nested
#   p.ppt.slope.tempF <- list(formula = ~ppt.cat + slope.sd + tempF)
#   
#   p.ppt.slope <- list(formula = ~ppt.cat + slope.sd)
#   p.ppt.tempF <- list(formula = ~ppt.cat + tempF)
#   p.slope.tempF <- list(formula = ~slope.sd + tempF)
#   #p.strat.tempF <- list(formula = ~stratum + tempF)
#   
#   p.ppt <- list(formula = ~ppt.cat)
#   p.slope <- list(formula = ~slope.sd)
#   p.tempF <- list(formula = ~tempF)
#   p.1 <- list(formula = ~1)
#   #p.strat <- list(formula = ~stratum)
#   
#   # Detection: non-nested
#   p.grnd.jdate <- list(formula = ~grnd + jdate)
#   
#   # Occupancy
#   Psi.test <- list(formula = ~ stratum + slope.mn + bare.indx + ndvi16.ucl)
#   
#   
#   # Tell MARK the model type and create DF of all combos of formulas 
#   cml <- create.model.list("Occupancy") 
#   
#   # Wrapper for mark function that works with the clm DF to run all combos, name them and save into a list of models / AIC model set
#   return(mark.wrapper(cml, data = wtpd.proc, ddl = wtpd.ddl, adjust = FALSE, output = FALSE)) # adjust...
#   
# } # END FUNCTION
# 
# # Run function to run all combos of p and Psi
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mods.stage1 <- combos.s1()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FUNCTION
#~~~~~~~~~~~~~
# Laake on LRT in RMark: www.phidot.org/forum/viewtopic.php?f=21&t=2457
# - Note that the lnl field is actually -2LnL
# - And for many mods (such as those with indiv covariates or occupancy mods), 
#   deviance is same as lnl
# --> www.phidot.org/forum/viewtopic.php?f=21&t=2084&hilit=lnl+instead+of+deviance
# --> Gentle MARK book, sidebar starting on 5-2


# SPRING 2018
#~~~~~~~~~~~~~
# Function comparing correctly except not comparing interactions of form (x + y + x*y) BUT
#   is comparing quadratics correctly now that added "main" argument.

# - terms within int.term and main arguments need to be in same order
# --> int.term = c(slope.sq, temp.sq), main = c("slope", "temp") *NOT* main = c("temp", "slope")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


aic.nested.occup <- function(model.list, aic.thresh, param, int.term=NULL, main=NULL){
  {
    # model.list <- mods.s1.test2
    # aic.thresh <- 0
    # param <- "p"
    #int.term <- c("slp.mn.sq")
    #main <- c("slope.mn")

    # Separate models and AIC table, which are in same list
    marklist <- model.list # will need original marklist object for remove.mark near end of function
    set <- model.table(model.list, adjust = FALSE, use.lnl = TRUE, model.name = FALSE)
    set <- set[ , c(param, "model", "npar", "AICc", "DeltaAICc", "weight", "Neg2LnL")]
    names(set)[1] <- "Modnames" 
    set$Modnames <- gsub("~", "", set$Modnames)
    set$mod.num <- row.names(set)
    model.list <- model.list[-length(model.list)] # remove model.table from model.list
    
    aic.thresh <- ifelse(aic.thresh > 0, -aic.thresh, aic.thresh) # convert aic.thresh to negative
    
    # Predictor counting
    # - Factors are tricky. Make DF with # of predictors rather
    #   than # parameters (so a factor is counted as 1 rather than nlevels-1)
    pred.df <- data.frame(Modnames = set$Modnames,
                          npred = sapply(strsplit(as.character(set$Modnames), "+", fixed=TRUE), length))
    pred.df <- pred.df[order(pred.df$npred), ]
    pred.df$Modnames <- gsub("~", "", pred.df$Modnames)
    
    # Separate model name into separate columns, so each pred is in a col
    pred.df <- separate(pred.df, Modnames, into = paste("X", 1:max(pred.df$npred), sep=""), 
                        sep = "\\+", remove=FALSE, fill = "right")
    # Trim white space
    pred.df[ ,grep("X", names(pred.df))] <- apply(pred.df[ ,grep("X", names(pred.df))], 2, function(x) trimws(x, "both"))
    df.max <- max(pred.df$npred)
    
    # When int.term isn't NULL: function compares models that differ by 1 npred but a model with and w/o an 
    # interaction differs by 2 preds, so substract 1 from npred of models that contain int.term so function will
    # compare those models to nested models without the interaction.
    # - This is different than how I'm dealing with factors, since it's easy to count factors as 1 pred even when
    #   they're multiple parameters (since it is a single predictor in the model). 
    # - LRT's etc. that rely on true # params are based on param counts, not npred.
     # if(!is.null(int.term)){
     #  for(i in 1:length(int.term)){
     #    term.i <- int.term[i]
     #    pred.df$npred[grep(term.i, pred.df$Modnames)] <- pred.df$npred[grep(term.i, pred.df$Modnames)]-1
     #  }
     # }
    
    # Merge
    set <- merge(set, pred.df, by = "Modnames", all = TRUE, sort = FALSE)
    #set[ ,c("Modnames", "npar", "npred")]
    
    # npreds
    npred.p1 <- sort(unique(set$npred + 1))
    npred.p1 <- intersect(unique(set$npred), npred.p1) # intersection will drop comparisons can't/don't wanna make
    
    # Another post-hoc fix so will compare quadratics to univariate term that's not squared
    # (the "-1" npred step above makes it so "need" to do this fix)
    #base.int.model <- paste(main, int.term, sep = " + ")
    
    # which(set$Modnames %in% base.int.model)
    # which(set$Modnames %in% c("slope.mn + slp.mn.sq","tempF + slope.mn"))
    
    #set[ ,c("Modnames", "npar", "npred")]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Nested loop 1: For models that differ by a single pred (not param), calculate AIC 
    #               differences and LRT for all model comparisons, nested or not.
    #               Remove non-nested comparisons in nested loop 2.
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    compare <- data.frame() # to hold 1st nested loop ouput
    
    # OUTER LOOP: loop npred.p1 to subset to groups of predictors
    for(i in 1:length(npred.p1)){ # START OUTER
      # Smaller mods
      aic.small <- set[set$npred == npred.p1[i]-1, ] # subset to npred i mods
      preds <- aic.small$Modnames # subset to npred i predictors
      
      # Bigger mods
      aic.big <- set[set$npred == npred.p1[i], ] # subset to mods where npred = npred+1
      #aic.big <- set[set$npred == npred.p1[i] | set$Modnames %in% base.int.model, ] # subset to mods where npred = npred+1
      aic.big <- aic.big[ ,!names(aic.big) %in% c("DeltaAICc", "weight")] # drop some extra cols
      
      # INNER LOOP: loop thru predictors extracted in outer loop
      for(j in 1:length(preds)){ # START INNER
        # Of the smaller npred[i] mods, grab ones that contain pred j
        aic.small.j <- aic.small[aic.small$Modnames == preds[j], ]
        aic.big$mod.small <- aic.small.j$Modnames # Add smaller mod name to bigger set
        
        # AIC diff, log-lik (LL) diff , DF diff
        aic.big$AIC.diff <- aic.big$AICc - aic.small.j$AICc
        aic.big$Neg2LL.diff <- aic.big$Neg2LnL - aic.small.j$Neg2LnL # may need to make argument for dev vs LL ? 
        aic.big$npar.small <- aic.small.j$npar
        aic.big$DF.diff <- aic.big$npar - aic.big$npar.small
        
        # LRT: DF diff must be 0 or greater
        not.neg <- aic.big[aic.big$DF.diff >= 0, ]
        not.neg$LRT <- (1 - pchisq((aic.small.j$Neg2LnL - not.neg$Neg2LnL), 
                                   df = not.neg$npar - not.neg$npar.small))
        aic.big$LRT <- not.neg$LRT[match(aic.big$Modnames, not.neg$Modnames)]
        
        # Save as go
        compare <- rbind(compare, aic.big)
        
      } # END INNER
      
    } # END OUTER
    compare$npred <- NULL
    #head(compare[ ,c("Modnames", "npar", "mod.small","npar.small", "DF.diff", "LRT")], 15)
    #head(compare[ ,c("Modnames", "mod.small","DF.diff", "LRT")], 5)
    
    
    # Prep for nested loop 2
    #~~~~~~~~~~~~~~~~~~~~~~~~
    #pred.nest <- ifelse(is.null(int.term), df.max-1, df.max)
    compare$npred.big <- pred.df$npred[match(compare$Modnames, pred.df$Modnames)]
    compare$npred.small <- pred.df$npred[match(compare$mod.small, pred.df$Modnames)]
    
    # Split preds of smaller mods into separate columns
    # compare <- separate(compare, mod.small, into = c(paste("small", 1:pred.nest, sep = ".")), 
    #                     sep = "\\+", fill = "right", remove=FALSE)
    compare <- separate(compare, mod.small, into = c(paste("small", 1:df.max, sep = ".")), 
                        sep = "\\+", fill = "right", remove=FALSE)
    compare[ ,grep("small.", names(compare))] <- apply(compare[ ,grep("small.", names(compare))], 2, function(x) trimws(x, "both"))
    
    # Subset
    compare.sub <- cbind(Modnames = compare$Modnames, mod.small = compare$mod.small, compare[ ,grep("small.", names(compare))],
                         npred.big = compare$npred.big, npred.small = compare$npred.small)
    
    # Vector of smaller mod pred names
    cols <- grep("small.", names(compare.sub), value = TRUE)
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Nested loop 2: Identify nested and non-nested comparisons from
    #                 nested loop 1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    out.outer <- compare.sub[ , c("Modnames", "npred.big", "mod.small", "npred.small")] # hold loop 2 output
    
    for(i in 1:length(cols)){ # START OUTER
      
      temp <- compare.sub[ , c("Modnames", cols[i])]
      names(temp)[2] <- "small.x"
      test.nest <- vector() # create inner loop empty vector each iter of i
      
      for(j in 1:length(compare.sub$Modnames)){ # START INNER
        
        # Are the preds in the smaller model also in the larger model?
        small.x.j <- ifelse(is.na(temp$small.x[j]), temp$small.x[j], 
                            ifelse(temp$small.x[j] == 1, NA, temp$small.x[j]))
        
        nest <- grepl(small.x.j, temp$Modnames[j], fixed = TRUE)# fixed = TRUE
        nest <- ifelse(is.na(nest), FALSE, nest)
        test.nest <- c(test.nest, nest)
        
      } # END INNER
      
      temp$nest <- test.nest
      names(temp)[2:3] <- c(cols[i], paste("nest", i, sep = "."))
      temp$Modnames <- NULL
      out.outer <- cbind(out.outer, temp)
      
    } # END OUTER

    # Subset loop output to only the mods that are nested and for which
    # the # of nested preds is one less than npred of bigger model
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 1) Dynamic count of TRUE's per row: sum the nest.x logicals
    out.outer$nest.TRUE <- apply(out.outer[ ,grep("nest", names(out.outer))], 
                                 1, sum)
    
    # For when int.term is not NULL, otherwise the nest.TRUE count for models with the interaction will be same
    # as npred.big, resulting in nest.diff of 0, which doesn't make sense given the next steps.
    #if(length(int.term)>1){int.term <- paste(int.term, collapse = "|")}
    
    # Subtract 1 from nest.TRUE for all rows with an int.term in mod.small
    # - needed to deal with quadratics in form of: x1 + x1.sq (where square x1 before run model, not in model)
    #if(!is.null(int.term)){
    #  out.outer$nest.TRUE[grep(int.term, out.outer$mod.small)] <- out.outer$nest.TRUE[grep(int.term, out.outer$mod.small)]-1}
    
    # Add 1 back to nest.TRUE for all rows with an "*" indicating are part of an interaction of the form, x1 + x2 + x1*x2
    #if(!is.null(int.term)){
    #  out.outer$nest.TRUE[grep("\\*", out.outer$mod.small)] <- out.outer$nest.TRUE[grep("\\*", out.outer$mod.small)]+1}
    
    # 2) Diff btw npred of big model and count of TRUE for nested test
    out.outer$nest.diff <- out.outer$npred.big - out.outer$nest.TRUE
    
    #out.outer[,c("Modnames", "npred.big", "mod.small", "npred.small", "nest.TRUE", "nest.diff")]

    # 3) Subset to nested only
    nested.only <- out.outer[out.outer$nest.TRUE != 0, ]
    # big.small.same <- which(as.character(nested.only$Modnames) == as.character(nested.only$mod.small))
    # nested.only <- nested.only[-big.small.same, ]
    # nested.only$nest.diff[nested.only$Modnames %in% base.int.model] <- nested.only$nest.diff[nested.only$Modnames %in% base.int.model] + 1
    nested.only <- nested.only[nested.only$nest.diff == 1, ]
    
    nested.only[,c("Modnames", "npred.big", "mod.small", "npred.small", "nest.TRUE", "nest.diff")]
    #out.outer[ ,grep("small.", names(out.outer), invert = TRUE)]
    
    # For nested stats, make table with only the nested mod comparisons
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    compare.nested <- merge(nested.only, compare, 
                            by = c("Modnames", "mod.small", cols), 
                            all.x = TRUE, all.y = FALSE, sort = FALSE)
    
    # Sort based on mod name and K
    compare.nested <- compare.nested[order(compare.nested$npar, compare.nested$Modnames), ]
    
    # Cols to keep
    compare.nested <- compare.nested[ ,names(compare.nested) %in% 
                                        c("Modnames", "model","mod.small", "AIC.diff", "Neg2LL.diff", "DF.diff", "LRT", "mod.num")]
    
    # Rename
    names(compare.nested)[1] <- "mod.big"

    
    # What predictor is being tested? (Big vs small mod)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #test.of <- data.frame()
    test.of <- vector()
    for(i in 1:nrow(compare.nested)){
      small <- as.character(compare.nested$mod.small[i])
      small2 <- gsub(" ", "", small, fixed = TRUE)
      big <- as.character(compare.nested$mod.big[i])
      big2 <- gsub(" ", "", big, fixed = TRUE)
      
      test.of.i <- setdiff(strsplit(big2,"+", fixed = TRUE)[[1]],
                           strsplit(small2,"+", fixed = TRUE)[[1]])
      test.of.i <- ifelse(length(test.of.i)>1, paste(test.of.i, collapse="."), test.of.i)
      
      #temp <- data.frame(big=big, small=small, test.of.i=test.of.i, i=i)
      test.of <- c(test.of, test.of.i)
    }
    
    # Add to compare.nested DF
    compare.nested$pred.tested <- test.of
    
    #compare.nested[compare.nested$DF.diff > 0, c("mod.big", "mod.small", "DF.diff", "pred.tested")]
    compare.nested <- compare.nested[compare.nested$DF.diff > 0, ]

    #compare.nested[, c("mod.big", "mod.small", "DF.diff", "pred.tested")]


    # Finding uninformative mods and adjusting AIC table
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # aic.thresh is the chosen threshold that larger model needs to overcome to
    #   be considered informative. Multiply the threshold by the DF.diff so if
    #   mods diff by 1 param then thresh is the same but if diff by >1 param, 
    #   then thresh is larger. 
    # Then calculate: (AIC big mod - AIC small mod) - (DF diff * threshold diff).
    #   If diff is less than zero, then the bigger model is informative, otherwise
    #   bigger mod is uninformative
    compare.nested$threshXdf.diff <- compare.nested$DF.diff * aic.thresh
    compare.nested$AICdiff.threshDF <- compare.nested$AIC.diff - compare.nested$threshXdf.diff
    uninform <- as.numeric(compare.nested$mod.num[compare.nested$AICdiff.threshDF >= 0])
    uninform <- uninform[!duplicated(uninform)]

    if(length(uninform)==0){stop("\n Not really an error \n All models are informative -- that's just surprising")}

      # Remove uninformative mods from list and recalc AIC weights
      set.inform <- remove.mark(marklist, uninform)
      set.inform <- model.table(set.inform, adjust = FALSE, use.lnl = TRUE, model.name = FALSE)
      mods.uninform <- data.frame(mod.num = uninform)
      mods.uninform$model <- marklist$model.table[[param]][match(mods.uninform$mod.num, row.names(marklist$model.table))]
      
      # AIC table with all the mods but uninform mods have NA for weight
      aic.corrected <- set
      aic.corrected$weight <- NA
      
      # Add weights from set.inform
      aic.corrected$weight <- set.inform$weight[match(aic.corrected$model, set.inform$model)]
      aic.corrected$wt2 <- aic.corrected$weight
      aic.corrected$wt2[is.na(aic.corrected$wt2)] <- 0
      aic.corrected$cum.w <- cumsum(aic.corrected$wt2)
      
      # Clean-up tables
      x.cols <- grep("X", names(set), value = TRUE)
      set <- set[ ,!names(set) %in% c(x.cols, "model", "mod.num", "npred")]
      aic.corrected <- aic.corrected[ ,!names(aic.corrected) %in% c(x.cols, "model", "mod.num", "npred", "wt2")]
      compare.nested <- compare.nested[ , c("mod.big", "mod.small", "pred.tested", "AIC.diff", "DF.diff", "threshXdf.diff", 
                                            "AICdiff.threshDF", "Neg2LL.diff", "LRT")]
      # set[ ,-c(1,2)] <- round(set[ ,-c(1,2)], 2)
      # aic.corrected[ ,-c(1,2)] <- round(aic.corrected[ ,-c(1,2)], 2)
      # compare.nested[ ,-c(1,2)] <- round(compare.nested[ ,-c(1,2)], 2)
      
    }
    
    return(list(AICc.original = set, AICc.adjusted = aic.corrected, 
                compare.nested = compare.nested, mods.uninform = mods.uninform))

  
} # END FUNCTION


# # source("C:/Users/jceradin/Box Sync/Research/Analysis/Mark-Recap/aic_nested_occup.R")
# nest.test <- aic.nested.occup(mods.s3a.rm, aic.thresh = 0,
#                               param = "Psi",
#                               int.term = c("ndviGw15Rt", "anthro1k"))
# write.xlsx(mods.s3a.rm$model.table, "nest_test.xlsx", "aic.rmark", append=TRUE)
# write.xlsx(nest.test$compare.nested, "nest_test.xlsx", "compare.nested", append=TRUE)



# # Single-season WTPD model works
# out <- aic.nested.occup(model.list = mods.stage1, aic.thresh = 2, param = "p")
# out

# aic.nested.occup(model.list = mods.s1b, aic.thresh = 0, 
#                  param = "p", int.term = "slope.sd")
# 
# aic.nested.occup(model.list = mods.s1b, aic.thresh = 0, param = "p")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Could add comparison btw univariate models and intercept only
# - although poorly performing univariate models should be ranked
#   low and have little AIC weight regardless
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(lmtest)
data(mtcars)
head(mtcars)

int <- lm(mpg ~ 1, mtcars)
summary(int)

hp <- lm(mpg ~ hp, mtcars)
summary(hp)

lrtest(int, hp)


