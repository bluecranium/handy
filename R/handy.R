# -------------------------------------------------------------------	
# Return Vector of UCSC-style chromosome names, including the 22 human autosomes and the X and Y chromosomes
chrs <- function()
{
	paste("chr",c(1:22,"X","Y"),sep="")
}
# -------------------------------------------------------------------

# ------------------------------------------------------------------- 
# My basic ggplot theme - removes the gridlines for cleaner appearance
ggnice <- function()
{
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x=element_line(color="black"), axis.line.y=element_line(color="black"), axis.text = element_text(color="black"))
}
# -------------------------------------------------------------------	

# -------------------------------------------------------------------
# torso - grab the middle of a data.frame
torso <- function(x, n=6)
{
	mix <- round(nrow(x)/2,digits=0)
	x[mix:(mix+n),]
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# rtorso - grab random sample of rows from a data.frame
rtorso <- function(x, n=6)
{
	x[sample(1:nrow(x),n),]	
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# ctorso - grab random sample of consecutive rows from a data.frame
ctorso <- function(x, n=2, rows=3)
{
	s <- sample(1:nrow(x),n)
	i <- as.vector(sapply(s,FUN=function(x) seq(x,x+rows)))
	x[i,]	
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Improved list of objects, including memory usage and lengths
# Code from: http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(print(object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Automatically set optimal width
# Tip from: http://suppressingfire.livejournal.com/39392.html
.adjustWidth <- function(...){
       options(width=Sys.getenv("COLUMNS"))
       TRUE
}

nicewidth <- function()
{
	.adjustWidthCallBack <<- addTaskCallback(.adjustWidth)
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# morecolors - expand a set of colors to any number of colors
# Input: n=number of colors to output, rand=shuffle color order or not
# Output: Vector of hex color codes
morecolors <- function(n, rand=FALSE)
{
	# Brewer's qualitative palette "Set1" only has 9 values
	# Extrapolate from these to create palettes of any size
	pal <- colorRampPalette(brewer.pal(9,"Set1"))(n)
	if(rand==TRUE){pal <- sample(pal)}
	pal
}
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# gc - Calculate GC content of a DNAString
# Input: DNAStringSet object
# Output: Vector of GC contents for each sequence in the DNAStringSet
getgc <- function(seq)
{
	g <- alphabetFrequency(seq)[,3]
	c <- alphabetFrequency(seq)[,2]
	a <- alphabetFrequency(seq)[,1]
	t <- alphabetFrequency(seq)[,4]
	gc <- (g+c) / (a+t+g+c)
	gc
}
# -------------------------------------------------------------------

# --------------------------------------------------------------------
# print recursive summary of a list with nested lists (list of lists)
nsummary <- function(mylist,level=1)
{
	cla <- sapply(mylist,class)
	for(i in 1:length(mylist))
	{
		print(paste0(paste(rep("-> ",level),collapse=""),"[[",i,"]] '",names(mylist)[i],"' (",cla[i],")"))
		if(cla[i]=="list")
		{
			nsummary(mylist[[i]],level=level+1)
		}
	}
}
# --------------------------------------------------------------------
