#' @title Check Executable
#' 
#' Test if a file is executable
#' 
#' @param exe.path Character string with full path to the file to check.
#' @return Integer 0 if exe.path is executable, otherwise 1.
#' @export check.if.executable
check.if.executable <- function(exe.path){
	res <- system(paste("command -v",paste0("'",exe.path,"'"),">/dev/null 2>&1 || { echo >&2 ''; exit 1; }"))
	res
}