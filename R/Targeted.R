

#' Return dilution curve limit of detection
#' @param dCurve data.frame
#' @param method c("blank","low")
#' @return lod
#' @export
#' @note  No note
#' @import magrittr
#' @references NA 
#' @examples print("No examples")
getLOD = function(dCurve,method="blank"){
  
  ### getLOD 
  maxConc =  max(dCurve$concentration,na.rm=T)
  blankInt = dCurve[dCurve$concentration == 0,]$intensity 
  maxIntMedian =  median( dCurve[dCurve$concentration %in% maxConc,]$intensity ,na.rm=T) 
  
  #### translate intenstiy to concentration
  # based on max con intensity only
  #blankMeasuredConc = (blankInt/maxIntMedian)*maxConc # scale intensities
  # using regression of high conc (above median conc) data
  slope = (lm( concentration ~ intensity  , data = subset(dCurve, concentration > concentration %>% unique %>% median(.,na.rm=T)  )) %>% coef)[2]
  blankMeasuredConc = blankInt*slope
  
  if("blank" %in% method ){
    
    #		# LIMIT OF DETECTION & LIMIT OF QUANTIFICATION
    #		
    #		# Blank Sample Method, Mani et al. http://www.biomedcentral.com/1471-2105/13/S16/S9
    #		# Assuming that random measurement errors are normally distributed, and with 5% 
    #		# risk of incorrectly claiming detection in the absence of analyte (alpha) or missing the detection of analyte (beta),
    #		# LOD = 3.29 sdB and LOQ = 3 x LOD = 10 sigmaB where sigmaB is the standard deviation of the blank sample.
    #		#
    #		# OR
    #		#Currie LA: Limits for qualitative detection and quantitative determination. Application to radiochemistry. Analytical Chemistry 1968, 40(3):586-593.
    #		
    #		# In Mani et al. the LOD is calculated as the sd of "blank (light r.t. region) to heavy peptide ratio". This given the lod as a fraction of the heavy peptide concentration.
    #		# We calculate the LOD as a fraction of the most intense peptide concentration
    
    lod = 3.29 * sd(blankMeasuredConc,na.rm=T)
    return(lod)
  }else{
    
    #		# Blank and low Concentration Sample Method, Mani et al. http://www.biomedcentral.com/1471-2105/13/S16/S9 and http://www.nature.com/nbt/journal/v27/n7/full/nbt.1546.html
    #		# LOD = meanB +t(1-b) (sdB + sdS)/sqrt(n)
    
    lowConc = min(dCurve$concentration[dCurve$concentration > 0],na.rm=T)
    lowInt = dCurve[dCurve$concentration %in% lowConc,]$intensity
    lowMeasuredConc = (lowInt/maxIntMedian)*maxConc
    n = length(blankMeasuredConc)+length(lowMeasuredConc)
    lod = mean(blankMeasuredConc,na.rm=T) +  (qt(0.95,n-1) * (sd(blankMeasuredConc,na.rm=T)+sd(lowMeasuredConc,na.rm=T))/sqrt(n))
    return(lod)
  }
}

#' Return dotProduct of two vectors
#' @param u vector 1
#' @param v vector 2
#' @param normalize dp TRUE/FALSE 
#' @return dp
#' @export
#' @note  No note
#' @references NA 
#' @examples print("No examples")
dotProduct <- function(u,v,norm=F){
  
  dp <- v %*% u
  
  if(norm){
    # dpNorm = (u * v) / (sqrt(v*v) * sqrt(u*u) )  
    dpNorm <- dp / (sqrt(dotProduct(v,v))* sqrt(dotProduct(u,u)))
    # dp norm is the dp of uN * vN, where uN is u scaled to unit length
    # i.e. uN = u / sqrt(u*u)  
    
    # OR
    # cos(alpha), where alpha is the angle between u and v
    return(as.vector(unlist(dpNorm)))
  }
  
  return(as.vector(unlist(dp)))
  
}

#' Return dotProducts to most transition intensities of most intense runs 
#' @param eset ExpressionSet
#' @param u vector 1
#' @param v vector 2
#' @param nbRefRuns (default top  4) 
#' @param dp matrix of dp's per peptide and run 
#' @return dp
#' @export
#' @import magrittr, Biobase
#' @note  No note
#' @references NA 
#' @examples print("No examples")
getAllDotProduct = function(eset, nbRefRuns = 4){
  dp = data.frame()
  allPeptides = unique(fData(eset)$peptide)
  for(peptide in  allPeptides){
    
    intensities = subset(exprs(eset), fData(eset)$peptide == peptide )
    intensities[is.na(intensities)] = 0  
    
    # ref Vec is sum across 4 most intenst runs
    refVec =  intensities[,order(colSums(intensities,na.rm=T),decreasing=T)][,1:min(nbRefRuns,ncol(intensities))] %>% rowSums(.,na.rm=T)
    dp %<>% rbind(apply(intensities,2,function(v){dotProduct(v,refVec, norm=T)}))
  }
  colnames(dp) = colnames(intensities)
  rownames(dp) = allPeptides
  return(dp)
}

#' Plot dilution curve
#' @param dCurve data.frame columns concentration, intensity
#' @param lod limit of detection
#' @return 
#' @export
#' @import magrittr, ggplot2
#' @note  No note
#' @references NA 
#' @examples print("No examples")
ggDilutionCurve = function(dCurve,lod, title=""){
  
  lodNq = data.frame(limit=c(lod,lod*3) %>% round(1) ,type=c("lod","loq"))
  
  p = ggplot(dCurve,
             aes( x =concentration , y=intensity)) 
  
  # add title
  p = p +ggtitle(label =title )
  
  # add points
  p = p + geom_point(size=3)
  
  # log scale
  p = p + scale_x_log10(name="fmol on Column")
  p = p + scale_y_log10(name="ms1 Intensity", limits=range(df$intensity))
  
  # ablines
  if(!is.na(lod)){
    p = p + geom_vline(data=lodNq 
                       ,aes(xintercept = limit, color=type)
                       ,linetype="longdash"
                       ,size=1)
  }
  # legend
  # p = p + scale_color_manual(values=c("red","blue","black")
  #                            , labels=c(paste("LOD:",lodNq$limit[1] %>% signif(.,3)  %>% format(., scientific=T)) ,paste("LOQ:",lodNq$limit[2] %>% signif(.,3)  %>% format(., scientific=T) ), "R2") ) 
  # 
  #if(!is.na(lod)){
    p = p + scale_color_manual(values=c("red","blue")
                               , labels=c(paste("LOD:",lodNq$limit[1] %>% signif(.,3)  %>% format(., scientific=T)) ,paste("LOQ:",lodNq$limit[2] %>% signif(.,3)  %>% format(., scientific=T) )) )
    
  #}
  
  
  # reg curve
  p = p + geom_smooth(data =subset(dCurve, concentration > lodNq$limit[2])
                      , method = "lm"
                      , se = F
                      , fullrange=T
                      , color="grey"
                      ,size=1
  )

  # r2 labels
  p = p + annotate("text", x = min(dCurve$concentration) , y = max(dCurve$intensity)
                   , label = paste("R^2: ", with(log(dCurve), cor(concentration,intensity,use="pairwise.complete.obs"))^2 %>% round(.,2)) 
                   , parse=T, size=8, hjust=0, vjust=1)  
  # theme
  p = p + theme_bw()
  p = p + theme(title=element_text(size=12, face='bold')
                ,axis.title = element_text(size = 20)
                ,axis.text =  element_text(size = 20)
                ,axis.text.y =  element_text(angle=45)
                ,legend.text = element_text(size = 20)
                ,legend.title = element_blank()
                ,legend.justification=c(1,0)
                ,legend.position=c(0.975,0.1)
                ,legend.background = element_rect(fill="transparent")
  ) 
  return(p)  
}


#calibrationCurve <- function(eset,method="blank"){
#	
#	out <- list()
#	class(out) <- "calibrationCurve"
#	
#	### store the curve
#	out$curve <- data.frame(concentration=rep(fData(eset)$concentration[fData(eset)$concentration > 0],ncol(eset))
#			, intensity= as.vector(unlist(exprs(eset)[fData(eset)$concentration > 0,]))	)
#	logCurve <- log10(out$curve)
#	out$fit <- lm(intensity ~ concentration, data=logCurve)
#	
#	out$lod <- NA
#	out$loq <- NA
#	### store ExpressionSet matching assay
#	out$eset <- eset
#	#print(names( fData(out$eset)))
#	#[1] "Peptide.Sequence" "Protein.Name"     "Precursor.Mz"     "Precursor.Charge"
#	#[5] "Retention.Time"   "concentration"    "dilutionCurveId" 
#
#	out$label <- paste(fData(out$eset)$dilutionCurveId[1],fData(out$eset)$Protein.Name[1])
#	
#	if("blank" %in% method ){
#		
#		# LIMIT OF DETECTION & LIMIT OF QUANTIFICATION
#		
#		# Blank Sample Method, Mani et al. http://www.biomedcentral.com/1471-2105/13/S16/S9
#		# Assuming that random measurement errors are normally distributed, and with 5% 
#		# risk of incorrectly claiming detection in the absence of analyte (alpha) or missing the detection of analyte (beta),
#		# LOD = 3.29 sdB and LOQ = 3 x LOD = 10 sigmaB where sigmaB is the standard deviation of the blank sample.
#		#
#		# OR
#		#Currie LA: Limits for qualitative detection and quantitative determination. Application to radiochemistry. Analytical Chemistry 1968, 40(3):586-593.
#		
#		# In Mani et al. the LOD is calculated as the sd of "blank (light r.t. region) to heavy peptide ratio". This given the lod as a fraction of the heavy peptide concentration.
#		# We calculate the LOD as a fraction of the most intense peptide concentration
#		
#		maxConc <-  max(fData(eset)$concentration,na.rm=T)
#		blankInt <- exprs(eset)[fData(eset)$concentration == 0,] 
#		maxIntMedian <-  median( exprs(eset)[fData(eset)$concentration %in% maxConc,] ,na.rm=T) 
#		blankMeasuredConc <- (blankInt/maxIntMedian)*maxConc
#		
#		out$lod <- 3.29 * sd(blankMeasuredConc,na.rm=T)
#		out$loq <- 3 * out$lod
#	}else if("low" %in% method ){
#		
#		# Blank and low Concentration Sample Method, Mani et al. http://www.biomedcentral.com/1471-2105/13/S16/S9 and http://www.nature.com/nbt/journal/v27/n7/full/nbt.1546.html
#		# LOD = meanB +t(1-b) (sdB + sdS)/sqrt(n)
#		
#		maxConc <-  max(fData(eset)$concentration,na.rm=T)
#		blankInt <- exprs(eset)[fData(eset)$concentration == 0,] 
#		lowConc <- min(fData(eset)$concentration[fData(eset)$concentration > 0],na.rm=T)
#		lowInt <- exprs(eset)[fData(eset)$concentration %in% lowConc,]
#		maxIntMedian <-  median( exprs(eset)[fData(eset)$concentration %in% maxConc,] ,na.rm=T) 
#		blankMeasuredConc <- (blankInt/maxIntMedian)*maxConc
#		lowMeasuredConc <- (lowInt/maxIntMedian)*maxConc
#		n <- length(blankMeasuredConc)+length(lowMeasuredConc)
#		out$lod <- mean(blankMeasuredConc,na.rm=T) +  (qt(0.95,n-1) * (sd(blankMeasuredConc,na.rm=T)+sd(lowMeasuredConc,na.rm=T))/sqrt(n))
#		out$loq <- 3 * out$lod
#			
#	}
#	
#	return(out)
#	
#}


#plot.calibrationCurve <- function(x, pch=19, lwd=2,cex.axis=1.5,cex.lab=1.5,cex=1.5,cex.main=1.5,xlab="Concentration",ylab="Area Under Curve",...){
#	
#	### include lod if out of range
#	xlim=c(min(x$lod,min(x$curve[,1],na.rm=T),na.rm=T),max(x$curve[,1],na.rm=T))
#	
#	plot.default(x$curve
#			, log="xy"
#			, main=x$label
# 			, pch=pch
#			, cex.axis=cex.axis 
#			, cex.lab=cex.lab
#			, cex=cex
#			, cex.main=1.5
#			, xlab=xlab
#			, ylab=ylab
#			, xlim=xlim
#			, ...
#		
#	)
#	
#	### linear model for points above loq
#	slope <- NA
#	selData <- log10(x$curve[x$curve[,1] > x$loq  ,])
#	if(length(unique(selData[,1])) > 1){
#		fit <- lm(intensity ~ concentration, data=selData)
#		slope <- coef(fit)[2]
#		abline(fit, col="darkgrey",lwd=lwd)
#	}
#	
#	abline(v=x$lod, col="red",lwd=lwd, lty=2)
#	abline(v=x$loq, col="blue",lwd=lwd, lty=2)
#	
#	### bottom right 
#	digits <- 3
#	legend("bottomright"
#			,legend=c(as.expression(bquote(R^2*"" == .(round(summary(x$fit)$r.squared,digits)))
#					)
#					,paste("K=",round(slope,digits)) 
#					,paste("LOD=",round(x$lod,digits)) 
#					,paste("LOQ=",round(x$loq,digits)) 
#			)
#			,text.col=c("black","darkgrey","red","blue"), box.lwd=NA, cex=1.5 )
#	
#}
