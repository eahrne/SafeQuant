# TODO: Add comment
# 
# Author: ahrnee-adm
###############################################################################


#' get SQ root dir
#' @return dir path String
#' @export
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
getSQRootDir <- function(){
	
	scriptPath <-  dirname(sys.frame(1)$ofile)
	dirs <- c(strsplit(dirname(normalizePath(scriptPath,mustWork = T)),
					"/|\\\\")[[1]], basename(scriptPath))
	#sqRootDirNb <-  max(which(dirs == "SafeQuant"))
	sqRootDirNb <-  max(which(grepl("SafeQuant",dirs))) # in-case check match SafeQunat.Rcheck)
	
	return(paste(paste(dirs[1:sqRootDirNb],collapse="/"),c("/"),collapse="",sep=""))
}

