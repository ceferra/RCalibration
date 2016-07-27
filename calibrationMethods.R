
### Code for the paper
##
#Applied Intelligence
#June 2013, Volume 38, Issue 4, pp 566–585
#On the effect of calibration in classifier combination
#Antonio Bella, Cèsar Ferri, José Hernández-Orallo, María José Ramírez-Quintana

# Function to assign the calibrated probability from the "$interval" and "$calProb"
# INPUTS
# ------
# cal: list with $interval and $calProb
# pTest: estimated probability of 1 class (array)

# OUTPUT
# ------
# pCal: calibrated probability associated with pTest (array)

asignarPCal <- function(cal,pTest){
	interval <- c(0.0,cal[[1]])
	pCal <- pTest
	for(i in 1:(length(interval)-1)){
		pCal[which((pTest>=interval[i]) & (pTest<interval[i+1]))]<-cal[[2]][i]
	}
	return(pCal)
}


# =================== #
# CALIBRATION METHODS #
# =================== #

# INPUTS
# ------
# yVal: real class in the validation set (array)
# pVal: estimated probability in the validation set (matrix, 1 column per class)
# pTest: estimated probability in the test set (matrix, 1 column per class)

# OUTPUT
# ------
# Calibrated probability associated with pTest (matrix, 1 column per class)

# ----------------------------------------------------------------- #
# Isotonic Regression Calibration Method (CORElearn implementation) #
# ----------------------------------------------------------------- #

#PAVCal <- function(yVal,pVal,pTest){
#	pCal <- asignarPCal(calibrate(yVal,pVal[,1],method="isoReg"),pTest[,1])
#	return(cbind(pCal,1-pCal))
#}

#yVal<-c(0,0,0,0,0,1,1,0,1,1,1,1,0,1,1,1,1,1,0,0,0,0,1,0,1,0,0,0,1,0)
#pVal<-matrix(c(0.2,0.8,0.3,0.7,1,0,0.2,0.8,0.2,0.8,0.8,0.2,0.9,0.1,0.9,0.1,0.5,0.5,0.6,0.4,0.9,0.1,0.3,0.7,0.6,0.4,0.8,0.2,0.6,0.4,0.8,0.2,0.1,0.9,0.8,0.2,1,0,0.6,0.4,0.1,0.9,0.8,0.2,0.1,0.9,0.4,0.6,0.1,0.9,0.8,0.2,0.1,0.9,0.9,0.1,0.4,0.6,0.9,0.1),ncol=2,byrow=TRUE)

PAVCal <- function(yVal,pVal,pTest){
	p<-pVal[,1]
	pVal<-p
	yI <- as.integer(yVal)
	yI[yI==2] <- 0
	# Sort by estimated probability decreasingly
	pInd <- order(pVal)
	yOrd <- yI[pInd]
	pOrd <- pVal[pInd]
	
	auxP <- pOrd
	pj <- {}
	nj <- {}
	rj <- {}
	while(length(auxP)>0){
		pj <- c(pj,auxP[1])
		repes <- which(auxP==auxP[1])
		nj <- c(nj,length(repes))
		rj <- c(rj,sum(yOrd[which(pOrd==auxP[1])])/length(repes))
		auxP <- auxP[-repes]	
	}
	cambia<-TRUE
	while(cambia){
		cambia<-FALSE
		pjAux<-njAux<-rjAux<-{}
		if(length(pj)>1){
			for(i in 1:(length(pj)-1)){
				if(rj[i+1]<rj[i]){
					pjAux<-c(pjAux,(pj[i]*nj[i]+pj[i+1]*nj[i+1])/(nj[i]+nj[i+1]))			
					njAux<-c(njAux,nj[i]+nj[i+1])
					rjAux<-c(rjAux,(rj[i]*nj[i]+rj[i+1]*nj[i+1])/(nj[i]+nj[i+1]))
					if((i+2)<=length(pj)){
						for(j in (i+2):length(pj)){
							pjAux<-c(pjAux,pj[j])
							njAux<-c(njAux,nj[j])
							rjAux<-c(rjAux,rj[j])
						}
					}
					cambia<-TRUE
					pj<-pjAux
					nj<-njAux
					rj<-rjAux
					break
				}
				else{
					pjAux<-c(pjAux,pj[i])
					njAux<-c(njAux,nj[i])
					rjAux<-c(rjAux,rj[i])	
				}	
			}
		}
	}
	
	if(length(pj)==1)
		interval<-c(1.0)
	else{
		interval <- 1:(length(pj)-1)
		for(i in 1:(length(pj)-1))
			interval[i] <- mean(pj[i:(i+1)])
	}
	
	pCal <- asignarPCal(list(interval=interval,calProb=rj),pTest[,1])
	return(cbind(pCal,1-pCal))
}


# ----------------------------------------------------- #
# Binning Calibration Method (CORElearn implementation) #
# ----------------------------------------------------- #

binningCal <- function(yVal,pVal,pTest,bins=10){
	pCal <- asignarPCal(calibrate(yVal,pVal[,1],method="binning",noBins=bins),pTest[,1])
	return(cbind(pCal,1-pCal))
}

####################################################################

# -------------------------- #
# Platt's Calibration Method #
# -------------------------- #

plattCal <- function(yVal,pVal,pTest){
	out <- pVal[,1]
	classes <- levels(yVal)
	target <- 	yVal==classes[1]
	len <- length(target)
	prior1 <- length(which(target==TRUE))
	prior0 <- len-prior1
	
	A <- 0
	B <- log((prior0+1)/(prior1+1))
	hiTarget <- (prior1+1)/(prior1+2)
	loTarget <- 1/(prior0+2)
	lambda <- 1e-3
	olderr <- 1e300
	pp <- rep((prior1+1)/(prior0+prior1+2),len)
	count <- 0
	for(it in 1:100){
		a<-b<-c<-d<-e<-0
		for(i in 1:len){
			if(target[i]==TRUE) 
				t <- hiTarget 
			else 
				t <- loTarget
			d1 <- pp[i]-t
			d2 <- pp[i]*(1-pp[i])
			a <- a + out[i]*out[i]*d2
			b <- b + d2
			c <- c + out[i]*d2
			d <- d + out[i]*d1
			e <- e + d1	
		}
		if(abs(d)<1e-9 & abs(e)<1e-9) 
			break
		oldA <- A
		oldB <- B
		err <- 0
		while(TRUE){
			det <- (a+lambda)*(b+lambda)-c*c
			if(det==0) 
				lambda <- lambda*10
			A <- oldA + ((b+lambda)*d-c*e)/det
			B <- oldB + ((a+lambda)*e-c*d)/det
			err <- 0
			for(i in 1:len){
				p <- 1/(1+exp(out[i]*A+B))
				pp[i] <- p
				err <- err - t*log(p) + (1-t)*log(1-p)
			}
			if(err<olderr*(1+1e-7)){
				lambda <- lambda*0.1
				break	
			}
			lambda <- lambda*10
			if(lambda>=1e6) 
				break
		}
		diff <- err-olderr
		scale <- 0.5*(err+olderr+1)
		if(diff>-1e-3*scale & diff<1e-7*scale) 
			count <- count+1 
		else 
			count <- 0
		olderr <- err
		if(count==3) 
			break	
	}
	
	pCal <- 1/(1+exp(A*pTest[,1]+B))
	return(cbind(pCal,1-pCal))
}

# ------------------------------------------- #
# Binning Calibration Method (my improvement) #
# ------------------------------------------- #
#yVal<-c(0,0,0,0,0,1,1,0,1,1,1,1,0,1,1,1,1,1,0,0,0,0,1,0,1,0,0,0,1,0)
#pVal<-matrix(c(0.2,0.8,0.3,0.7,1,0,0.2,0.8,0.2,0.8,0.8,0.2,0.9,0.1,0.9,0.1,0.5,0.5,0.6,0.4,0.9,0.1,0.3,0.7,0.6,0.4,0.8,0.2,0.6,0.4,0.8,0.2,0.1,0.9,0.8,0.2,1,0,0.6,0.4,0.1,0.9,0.8,0.2,0.1,0.9,0.4,0.6,0.1,0.9,0.8,0.2,0.1,0.9,0.9,0.1,0.4,0.6,0.9,0.1),ncol=2,byrow=TRUE)
binningCalIm <- function(yVal,pVal,pTest,bins=10){
	newLen<-as.integer(length(yVal)/bins)*bins
	sobran<-length(yVal)-newLen
	quitar<-sample(1:length(yVal),sobran)
	if(sobran>0)
		yVal<-yVal[-quitar]
	p<-pVal[,1]
	pVal<-p
	if(sobran>0)
		pVal<-pVal[-quitar]
	yI <- as.integer(yVal)
	yI[yI==2] <- 0
	# Sort by estimated probability decreasingly
	pInd <- order(pVal)
	yOrd <- yI[pInd]
	pOrd <- pVal[pInd]
	
	auxP <- pOrd
	pj <- {}
	nj <- {}
	rj <- {}
	while(length(auxP)>0){
		pj <- c(pj,auxP[1])
		repes <- which(auxP==auxP[1])
		nj <- c(nj,length(repes))
		rj <- c(rj,sum(yOrd[which(pOrd==auxP[1])])/length(repes))
		auxP <- auxP[-repes]	
	}
	newY<-{}
	newP<-{}
	for(i in 1:length(pj)){
		newY<-c(newY,rep(pj[i],nj[i]))
		newP<-c(newP,rep(rj[i],nj[i]))
	}
	lenBin<-newLen/bins
	interval <- rep(0.0,bins)
	calProb <- rep(0.0,bins)
	for(i in 0:(bins-1)){
		interval[i+1] <- mean(newY[((lenBin*i)+1):(lenBin*(i+1))])		
		calProb[i+1] <- mean(newP[((lenBin*i)+1):(lenBin*(i+1))])
	}
	intervalNew <- 1:(bins-1)
	for(i in 1:(bins-1))
		intervalNew[i] <- mean(interval[i:(i+1)])
	
	pCal <- asignarPCal(list(interval=intervalNew,calProb=calProb),pTest[,1])
	return(cbind(pCal,1-pCal))
}

# ------------------------------------------- #
# Binning Calibration Method (my improvement + GLOBAL) #
# ------------------------------------------- #
#yVal<-c(0,0,0,0,0,1,1,0,1,1,1,1,0,1,1,1,1,1,0,0,0,0,1,0,1,0,0,0,1,0)
#pVal<-matrix(c(0.2,0.8,0.3,0.7,1,0,0.2,0.8,0.2,0.8,0.8,0.2,0.9,0.1,0.9,0.1,0.5,0.5,0.6,0.4,0.9,0.1,0.3,0.7,0.6,0.4,0.8,0.2,0.6,0.4,0.8,0.2,0.1,0.9,0.8,0.2,1,0,0.6,0.4,0.1,0.9,0.8,0.2,0.1,0.9,0.4,0.6,0.1,0.9,0.8,0.2,0.1,0.9,0.9,0.1,0.4,0.6,0.9,0.1),ncol=2,byrow=TRUE)
binningCalGlobal <- function(yVal,pVal,pTest,bins=10){
	newLen<-as.integer(length(yVal)/bins)*bins
	sobran<-length(yVal)-newLen
	quitar<-sample(1:length(yVal),sobran)
	if(sobran>0)
		yVal<-yVal[-quitar]
	p<-pVal[,1]
	pVal<-p
	if(sobran>0)
		pVal<-pVal[-quitar]
	yI <- as.integer(yVal)
	yI[yI==2] <- 0
	# Sort by estimated probability decreasingly
	pInd <- order(pVal)
	yOrd <- yI[pInd]
	pOrd <- pVal[pInd]
	
	auxP <- pOrd
	pj <- {}
	nj <- {}
	rj <- {}
	while(length(auxP)>0){
		pj <- c(pj,auxP[1])
		repes <- which(auxP==auxP[1])
		nj <- c(nj,length(repes))
		rj <- c(rj,sum(yOrd[which(pOrd==auxP[1])])/length(repes))
		auxP <- auxP[-repes]	
	}
	newY<-{}
	newP<-{}
	for(i in 1:length(pj)){
		newY<-c(newY,rep(pj[i],nj[i]))
		newP<-c(newP,rep(rj[i],nj[i]))
	}
	lenBin<-newLen/bins
	interval <- rep(0.0,bins)
	calProb <- rep(0.0,bins)
	for(i in 0:(bins-1)){
		interval[i+1] <- mean(newY[((lenBin*i)+1):(lenBin*(i+1))])		
		#calProb[i+1] <- mean(newP[((lenBin*i)+1):(lenBin*(i+1))])
		calProb[i+1] <- mean(newP[((lenBin*i)+1):(lenBin*(i+1))])-mean(newY[((lenBin*i)+1):(lenBin*(i+1))])
	}
	intervalNew <- 1:(bins-1)
	for(i in 1:(bins-1))
		intervalNew[i] <- mean(interval[i:(i+1)])
	
	interval <- c(0.0,intervalNew)
	pTest<-pTest[,1]
	pCal <- pTest
	for(i in 1:(length(interval)-1)){
		pCal[which((pTest>=interval[i]) & (pTest<interval[i+1]))]<-pCal[which((pTest>=interval[i]) & (pTest<interval[i+1]))]+calProb[i]
	}
	
	pCal[pCal<0.0]<-0.0
	pCal[pCal>1.0]<-1.0
	#pCal <- asignarPCal(list(interval=intervalNew,calProb=calProb),pTest[,1])
	return(cbind(pCal,1-pCal))
}


# ------------------------------- #
# MULTI-CLASS CALIBRATION METHODS #
# ------------------------------- #

# ----------------------------------------------- #
# Similarity Binning Averaging Calibration Method #
# ----------------------------------------------- #

SBACal <- function(model,val,test,k,w){
	noAtts<-length(val[1,])
	pred<-predictClassProb(model,val)
	pVal<-pred$probabilities
	class <- val[,noAtts]
	newVal <- cbind(val[,1:(noAtts-1)],pVal,class)
	
	pred<-predictClassProb(model,test)
	pTest<-pred$probabilities
	class <- test[,noAtts]
	newTest <- cbind(test[,1:(noAtts-1)],pTest,class)
	
	m<-learnIBk(newVal,k,w)
	pred<-predictClassProb(m,newTest)
	pCal<-pred$probabilities
	
	return(pCal)
}

SBACal2 <- function(valL,yValL,pValL,testL,yTestL,pTestL,k,w){
	noAtts<-length(valL[1,])
	class<-yValL
	newVal <- cbind(valL[,1:(noAtts-1)],pValL,class)
	class<-yTestL
	newTest <- cbind(testL[,1:(noAtts-1)],pTestL,class)
	
	m<-learnIBk(newVal,k,w)
	pred<-predictClassProb(m,newTest)
	pCal<-pred$probabilities
	
	return(pCal)
}


normaliseP <- function(p){
	largo<-length(p[1,])			
	for(i in 1:length(p[,1])){
		sP<-sum(p[i,]) #Peta cuando esta suma es 0
		if(sP==0 || is.nan(sP)){
			for(j in 1:largo){
				p[i,j]<-1/largo
			}
		}
		else{
			for(j in 1:largo){
				p[i,j]<-p[i,j]/sP
			}
		}
	}
	return(p)	
}

sumNto1<-function(n){
	s<-0
	for(i in (n-1):1)	
		s<-s+i
	return(s)
}

PAV1vs1 <- function(yVal,pVal,pTest){
	classes <- levels(yVal)
	noClass <- length(classes)
	lenT<-length(pTest[,1])
	num<-0
	m.pCal<-{}
	for(i in 1:(noClass-1)){
		for(j in (i+1):noClass){
			yV<-yVal[yVal==classes[i] | yVal==classes[j]]
			pV<-pVal[(yVal==classes[i] | yVal==classes[j]),c(i,j)]
			pV<-normaliseP(pV)
			pT<-pTest[,c(i,j)]
			pT<-normaliseP(pT)
			
			yV<-factor(yV,classes[c(i,j)])
			
			pC<-PAVCal(yV,pV,pT)
			num<-num+1
			pCal<-matrix(c(1:lenT,rep(0,noClass*lenT)),ncol=noClass+1)
			pCal[,i+1]<-pC[,1]
			pCal[,j+1]<-pC[,2]
			m.pCal<-rbind(m.pCal,pCal)
		}	
	}
	pCal<-{}
	for(i in 1:lenT){
		aux<-matrix(m.pCal[m.pCal[,1]==i],ncol=noClass+1)
		for(j in 1:noClass){
			pCal<-c(pCal,mean(aux[,j+1]))	
		}
	}
	pCal<-matrix(pCal,byrow=TRUE,ncol=noClass)
	return(pCal)
}

#yVal<-c(3,1,3,2,2)
#yVal<-as.factor(yVal)
#pVal<-matrix(c(0.1,0.3,0.6,0.5,0.2,0.3,0.4,0.1,0.5,0.2,0.6,0.2,0.3,0.4,0.3),ncol=3,byrow=TRUE)
#pTest<-pVal
PAV1vsAll <- function(yVal,pVal,pTest){
	classes <- levels(yVal)
	noClass <- length(classes)
	lenT<-length(pTest[,1])
	pCal<-{}
	for(i in 1:noClass){
		yV<-yVal
		yV<-as.factor(yV!=classes[i])
		pV<-pVal[,i]
		pV<-c(pV,1-pV)
		pV<-matrix(pV,ncol=2)
		pT<-pTest[,i]
		pT<-c(pT,1-pT)
		pT<-matrix(pT,ncol=2)
				
		pC<-PAVCal(yV,pV,pT)
		pCal<-c(pCal,pC[,1])
	}
	pCal<-matrix(pCal,ncol=noClass)
	pCal<-normaliseP(pCal)
	return(pCal)
}


Platt1vs1 <- function(yVal,pVal,pTest){
	classes <- levels(yVal)
	noClass <- length(classes)
	lenT<-length(pTest[,1])
	num<-0
	m.pCal<-{}
	for(i in 1:(noClass-1)){
		for(j in (i+1):noClass){
			yV<-yVal[yVal==classes[i] | yVal==classes[j]]
			pV<-pVal[(yVal==classes[i] | yVal==classes[j]),c(i,j)]
			pV<-normaliseP(pV)
			pT<-pTest[,c(i,j)]
			pT<-normaliseP(pT)
			
			yV<-factor(yV,classes[c(i,j)])
			
			pC<-plattCal(yV,pV,pT)
			num<-num+1
			pCal<-matrix(c(1:lenT,rep(0,noClass*lenT)),ncol=noClass+1)
			pCal[,i+1]<-pC[,1]
			pCal[,j+1]<-pC[,2]
			m.pCal<-rbind(m.pCal,pCal)
		}	
	}
	pCal<-{}
	for(i in 1:lenT){
		aux<-matrix(m.pCal[m.pCal[,1]==i],ncol=noClass+1)
		for(j in 1:noClass){
			pCal<-c(pCal,mean(aux[,j+1]))	
		}
	}
	pCal<-matrix(pCal,byrow=TRUE,ncol=noClass)
	return(pCal)
}

Platt1vsAll <- function(yVal,pVal,pTest){
	classes <- levels(yVal)
	noClass <- length(classes)
	lenT<-length(pTest[,1])
	pCal<-{}
	for(i in 1:noClass){
		yV<-yVal
		yV<-as.factor(yV!=classes[i])
		pV<-pVal[,i]
		pV<-c(pV,1-pV)
		pV<-matrix(pV,ncol=2)
		pT<-pTest[,i]
		pT<-c(pT,1-pT)
		pT<-matrix(pT,ncol=2)
				
		pC<-plattCal(yV,pV,pT)
		pCal<-c(pCal,pC[,1])
	}
	pCal<-matrix(pCal,ncol=noClass)
	pCal<-normaliseP(pCal)
	return(pCal)
}


Binning1vs1 <- function(yVal,pVal,pTest){
	classes <- levels(yVal)
	noClass <- length(classes)
	lenT<-length(pTest[,1])
	num<-0
	m.pCal<-{}
	for(i in 1:(noClass-1)){
		for(j in (i+1):noClass){
			yV<-yVal[yVal==classes[i] | yVal==classes[j]]
			pV<-pVal[(yVal==classes[i] | yVal==classes[j]),c(i,j)]
			pV<-normaliseP(pV)
			pT<-pTest[,c(i,j)]
			pT<-normaliseP(pT)
			
			yV<-factor(yV,classes[c(i,j)])
			
			pC<-binningCalIm(yV,pV,pT)
			num<-num+1
			pCal<-matrix(c(1:lenT,rep(0,noClass*lenT)),ncol=noClass+1)
			pCal[,i+1]<-pC[,1]
			pCal[,j+1]<-pC[,2]
			m.pCal<-rbind(m.pCal,pCal)
		}	
	}
	pCal<-{}
	for(i in 1:lenT){
		aux<-matrix(m.pCal[m.pCal[,1]==i],ncol=noClass+1)
		for(j in 1:noClass){
			pCal<-c(pCal,mean(aux[,j+1]))	
		}
	}
	pCal<-matrix(pCal,byrow=TRUE,ncol=noClass)
	return(pCal)
}

Binning1vsAll <- function(yVal,pVal,pTest){
	classes <- levels(yVal)
	noClass <- length(classes)
	lenT<-length(pTest[,1])
	pCal<-{}
	for(i in 1:noClass){
		yV<-yVal
		yV<-as.factor(yV!=classes[i])
		pV<-pVal[,i]
		pV<-c(pV,1-pV)
		pV<-matrix(pV,ncol=2)
		pT<-pTest[,i]
		pT<-c(pT,1-pT)
		pT<-matrix(pT,ncol=2)
				
		pC<-binningCalIm(yV,pV,pT)
		pCal<-c(pCal,pC[,1])
	}
	pCal<-matrix(pCal,ncol=noClass)
	pCal<-normaliseP(pCal)
	return(pCal)
}

Binning1vsAllGlobal <- function(yVal,pVal,pTest){
	classes <- levels(yVal)
	noClass <- length(classes)
	lenT<-length(pTest[,1])
	pCal<-{}
	for(i in 1:noClass){
		yV<-yVal
		yV<-as.factor(yV!=classes[i])
		pV<-pVal[,i]
		pV<-c(pV,1-pV)
		pV<-matrix(pV,ncol=2)
		pT<-pTest[,i]
		pT<-c(pT,1-pT)
		pT<-matrix(pT,ncol=2)
				
		pC<-binningCalGlobal(yV,pV,pT)
		pCal<-c(pCal,pC[,1])
	}
	pCal<-matrix(pCal,ncol=noClass)
	pCal<-normaliseP(pCal)
	return(pCal)
}
