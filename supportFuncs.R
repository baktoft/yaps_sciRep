getDatTmb <- function(Hx, Hy, toa){
	datTmb <- list(
		H = matrix(c(Hx, Hy), ncol=2),
		nh = length(Hx),
		np = ncol(toa),
		toa=toa
	)
	return(datTmb)
}


getParams <- function(datTmb){
	list(
		XY = matrix(c(mean(datTmb$H[,1]) + rnorm(ncol(datTmb$toa), sd=50), mean(datTmb$H[,2]) + rnorm(ncol(datTmb$toa), sd=50)), ncol=2),	#positions
		top = na.approx(apply(datTmb$toa, 2, function(k) {median(k[k != -9999])}), rule=2),								#time of ping
		ss=rnorm(ncol(datTmb$toa), 1415, 5),																			#speed of sound
		logD_xy = -2,				#diffusivity of transmitter movement (D_xy in ms)
		# logD_b = -15,				#diffusivity of burst interval (D_b in ms)
		logSigma_bi = -5,			#sigma  burst interval (sigma_bi in ms)
		logD_v = -3,				#diffusivity of speed of sound (D_v in ms)
		logSigma_toa = -8,			#sigma for Gaussian 
		logScale = -3,				#scale parameter for t-distribution
		log_t_part = -3				#Mixture ratio between Gaussian and t
	)
}


#Inits should be in a credible range...
getInits <- function() {
	init_logD_xy <- -2
	# init_logD_b <- -15
	init_logSigma_bi <- -5
	init_logD_v <- -10
	init_logSigma_toa <- -5
	init_logScale <- -1
	init_log_t_part <- -2
	# inits <- c(init_logD_xy, init_logD_b, init_logD_v, init_logSigma_toa, init_logScale, init_log_t_part)
	inits <- c(init_logD_xy, init_logSigma_bi, init_logD_v, init_logSigma_toa, init_logScale, init_log_t_part)
	return(inits)
}


runYAPS <- function(datTmb,inits, params, silent){
	#Compile and run TMB-model
	dyn.load(dynlib("yaps"))
	obj <- MakeADFun(datTmb,params,DLL="yaps",random=c("XY","ss","top"),inner.control = list(maxit = 500000), silent=silent)
	opt <- nlminb(inits,obj$fn,obj$gr)

	#Obtain parameter estimates and standard deviations
	obj$fn()
	pl <- obj$env$parList()
	jointrep <- sdreport(obj, getJointPrecision=TRUE)
	param_names <- rownames(summary(jointrep))
	sds <- summary(jointrep)[,2]
	summ <- data.frame(param=param_names, sd=sds)
	plsd <- split(summ[,2], f=summ$param)

	#Extract data 
	sd_xy <- matrix(plsd$XY, ncol=2)
	yapsRes <- data.frame(x=pl$XY[,1], y=pl$XY[,2], top=pl$top+T0, sd_x=sd_xy[,1], sd_y=sd_xy[,2])
	return(yapsRes)

}

plotRes <- function(yapsRes){
	#load data for gps-track
	gps <- read.table('gps.txt', header=TRUE)
	#load data for track estimated using U-MAP version 1.3.3, Lotek Wireless Inc., Newmarket, Ontario, Canada
	umap <- read.table('umap.txt', header=TRUE)


	#Number of hydros detecting each ping
	nObs <- apply(toa, 2, function(k) sum(k!= -9999))

	# plot results
	par(mfrow=c(2,2))
	plot(y~x, data=gps, type="l", asp=1)
	points(Hy[-noHs]~Hx[-noHs], col="blue", pch=20, cex=2)
	lines(y~x, data=umap, col="green", lty=2)
	lines(y~x, data=yapsRes, col="red")

	plot(x ~ t, data=gps, type="l", lwd=2, xlab="Time", ylab=("x-coordinate"))
	lines(x ~ top, data=yapsRes, type="l", col="red")

	plot(y ~ t, data=gps, type="l", lwd=2, xlab="Time", ylab=("y-coordinate"))
	lines(y ~ top, data=yapsRes, type="l", col="red")

	plot(nObs ~ yapsRes$top, type="h", ylab="#detecting hydros", xlab="Time", col=(nObs<3) + 1, xlim=range(gps$t))
	points(nObs[which(nObs == 0)]+0.5 ~ yapsRes$top[which(nObs == 0)], col="blue", type="h") #to include pings detected by zero hydros
}