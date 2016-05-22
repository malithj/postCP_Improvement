GreedySegmente=function(Y,binmax) {
	
	n = length(Y)

	csY = cumsum(Y)

	Vlines=list()
	Vcols=list() # matrice triangulaire supérieure avec en [i,j] la valeur V1(i,j)
					# V1[i,j]=V1[i,j-1]+Delta[i,j-1] // Delta[i,i]=0

#############################################	
##### LOCAL FUNCTIONS #######################	
#############################################	
			
	fillV1Line=function(i,jmax){ ### local function
		##### Gain en C ????????????
		J = (i+1):(jmax-1)
		if (i==1) left=0 else left=csY[i-1]
		Delta=(Y[J]-(csY[J]-left)/(J-i+1))^2+(1-1/(J-i+1))/(J-i+1)*(Y[J]-(csY[J-1]-left)/(J-i))^2
		Vlines[[i]] <<- list()
		Vlines[[i]][c(i+1,J+1)] <<- c(0,cumsum(Delta))
		}

	fillV1Column=function(j,imin){ ### local function
		##### Gain en C ??????????????
		I = (j-2):(imin)
		if (imin==1) right=c(csY[I[1:(length(I)-1)]-1],0) else right=csY[I-1]
		Delta=(Y[I]-(csY[j-1]-right)/(j-I))^2+(1-1/(j-I))/(j-I)*(Y[I]-(csY[j-1]-csY[I])/(j-I-1))^2
		Vcols[[j]]<<-list()
		Vcols[[j]][c(j-1,I)] <<- c(0,cumsum(Delta))
		}
	
	calcCtr=function(l,r) { # calculate contrast where l and r are breaks
		if ((l+1)==r) return(NULL)
		index=(l+1):(r-1)
		return(-Vcols[[r]][[l]]+unlist(Vlines[[l]][index])+unlist(Vcols[[r]][index]))
		}
		
	decrCtr = function(newBreak,left,right) {
		#compute possible decreases of contrast by splitting bins
		
		# fill V1 line newbreak up to column right
		if (newBreak<right-1) fillV1Line(newBreak,right)		# fill V1 column newbreak from line left
		if (newBreak>left+1) fillV1Column(newBreak,left)	
		# new contrast variations	
		dCtr=c(calcCtr(left,newBreak), Inf, calcCtr(newBreak,right)) 
		}
#############################################	

	
	breaks=c(1,n+1) 
	fillV1Line(1,n+1) 		# fill V1 line 1 up to column n+1
	fillV1Column(n+1,1)	# fill V1 column n+1 from line 1

	decrement=c(+Inf,dd<-calcCtr(1,n+1),+Inf)

	while ((length(breaks)<(binmax+1))) {
		# new break = new minimal decrement
		newBreak=which.min(decrement)
		# find future position of new break into list 
		here=sum(breaks<newBreak) 
		# new list of breaks
		breaks=c(breaks[1:here], newBreak, breaks[(here+1):length(breaks)])
		# update decrement list
		left=breaks[here]
		right=breaks[here+2] ####### YVES !!!! 1->2
		decrement=c(decrement[1:left],decrCtr(newBreak,left,right),decrement[right:length(decrement)])
		}
	
	return(breaks)
	
}
