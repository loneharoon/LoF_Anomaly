rm(list=ls())
Sys.setenv(TZ="Asia/Kolkata")
library(xts)
library(Rlof) 

normalizedata <- function(x){
  # function used to normalize data
  xx<- vector(mode="numeric",length=length(x))
  minval= min(x)
  maxval= max(x)
  for (i in 1:length(x)){
    xx[i] <- (x[i]-minval)/(maxval-minval)
  }
  return (xx)
}

outlierfactor <- function(daymat){
  daymat <- daymat[complete.cases(daymat),] # removes rows containing NAs
  daymat1 <- as.xts(t(apply(daymat,1, function(x) abs(fft(x-mean(x)))^2/length(x))) ) #apply FFT
  daymat <- daymat1
  dis_mat <- dist(daymat) # compute distance matrix
  fit <- cmdscale(dis_mat, eig = TRUE, k = 2) # apply MDS
  x <- scale(fit$points[, 1])
  y <- scale(fit$points[, 2])
  daymat <- data.frame(x,y)
  
  df.lof2 <- lof(daymat,c(4:8),cores = 2) # apply LOF
  df.lof2 <-apply(df.lof2,2,normalizedata) # nrormaization
  anom_max <-  apply(df.lof2,1,function(x) round(max(x,na.rm = TRUE),2) ) #feature bagging for outlier detection
  return(anom_max)
}

## MAIN CODE STARTS####
fil <- "Path to directory containing dataset"
fname <- "filename within fil directory"
df <- read.csv(file=paste0(fil,fname),head=TRUE,sep=",")
df_xts <- xts(df$power,as.POSIXct(strptime(df$timestamp,"%Y-%m-%d %H:%M:%S"),origin="1970-01-01"))
daydata <-  split.xts(df_xts,f="days",drop=FALSE,K=1)
ful_daydata <-    lapply(daydata,function(x) {
  # this function handles missing data via interpolation
  timebounds <- paste(as.Date(start(x),tz="Asia/Kolkata"),"/",as.Date(end(x),tz="Asia/Kolkata")," 23",sep="")
  xtstime <- timeBasedSeq(timebounds)
  ob <- xts(1:length(xtstime), xtstime) 
  output <- na.approx(cbind(ob,x))# filling NA values
  output <- output[,-1]
  output
})

raw_mat <- matrix(0,nrow=length(ful_daydata),ncol = length(ful_daydata[[1]])) # should be no of observations/days X 24 (hours of day)
raw_mat <- t(sapply(ful_daydata,function(x) coredata(x)))
dates <- lapply(ful_daydata,function(x) as.Date(index(x[1]), tz="Asia/Kolkata",origin="1970-01-01" ))
row_names <- do.call(as.Date,list(unlist(dates)) )  
rownames(raw_mat) <- as.character(row_names)
col_names <- strftime(index(ful_daydata[[1]]),format="%H:%M:%S",tz="Asia/Kolkata" )
colnames(raw_mat) <- col_names

daydata <- raw_mat[1:31,] # get selected data from raw_mat data matrix
outlier_score <- outlierfactor(daydata) # This contains the anomaly score computed with the algorithm

