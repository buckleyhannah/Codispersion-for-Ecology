#####################################
### Codispersion analysis functions
#####################################

#### Modified codispersion function (modified from Cuevas et al. 2013)
#### See 'Box 1' in Buckley et al. (2016) for a detailed explanation.
#### Buckley, HL., Case, BS., Ellison, AM. 2016. Using codispersion analysis to characterize spatial patterns in species co-occurrences. Ecology 97(1): 32–39.

Codisp.Kern<-function(X, Y, h, k, gamma = 1)
{
  Kernel<-function(u,gamma)
  {
    v=0
    v=ifelse(abs(u)<=1,(1/beta(0.5,gamma+1))*(1-u^2)^gamma,0)
  }
  ifelse(X$coords==Y$coords,1,
         {
           break
           print("The coordinates of X and Y are different")
         })
  
  n=length(X$data)
  mX <- matrix(X$data,nrow=n,ncol=n,byrow=FALSE)
  mY <- matrix(Y$data,nrow=n,ncol=n,byrow=FALSE)
  MatriXX <- (mX - t(mX))^2
  MatriYY <- (mY - t(mY))^2
  MatriXY <- (mX - t(mX))*(mY - t(mY))
  mX <- matrix(X$coords[,1],nrow=n,ncol=n,byrow=FALSE)
  DesignX <- mX - t(mX)
  mY <- matrix(X$coords[,2],nrow=n,ncol=n,byrow=FALSE)
  DesignY <- mY - t(mY)
  
  KERNMATRIXX=Kernel((h[1]-DesignX)/k[1],gamma)*Kernel((h[2]-DesignY)/k[1],gamma)
  
  if(k[1]==k[2]&k[1]==k[3]){
    KERNMATRIYY=KERNMATRIXX
    KERNMATRIXY=KERNMATRIXX } else{
      KERNMATRIYY=Kernel((h[1]-DesignX)/k[2],gamma)*Kernel((h[2]-DesignY)/k[2],gamma)
      KERNMATRIXY=Kernel((h[1]-DesignX)/k[3],gamma)*Kernel((h[2]-DesignY)/k[3],gamma) 
    }
  
  Numerador=sum(KERNMATRIXY*MatriXY)/(2*sum(KERNMATRIXY))
  Denominador1=sum(KERNMATRIYY*MatriYY)/(2*sum(KERNMATRIYY))
  Denominador2=sum(KERNMATRIXX*MatriXX)/(2*sum(KERNMATRIXX))
  v1=Denominador1
  v2=Denominador2
  v3=Numerador
  v4=Numerador/sqrt(Denominador1*Denominador2)
  print(c(v1,v2,v3,v4))
}  

### Function to run codispersion window analysis (modified from Cuevas et al. 2013)

# geodata1 = first input data object (a geoR geodata object)
# geodata2 = second input object
# k = c(k1, k2, k3) = a vector of three bandwidth values for X, Y and XY
# max.window.size = the maximum lag distance
# lx = is the number of divisions in the lags in x (up to the max.window.size) that the kernal is applied over. Half of these divisions are in the 'left', or positive direction, and half are in the 'right', or negative x direction.
# ly = is the number of divisions in the lags in y (up to the max.window.size) that the kernal is applied over in the 'up' direction of the plot

codisp.fn <- function(geodata1, geodata2, k = k, max.window.size = max.window.size, lx = 20, ly = 10) {
  X = geodata1  # input data process 1
  Y = geodata2  # input data process 2
  k = c(k[1],k[2],k[3]) # Set the bandwith for the kernel
  
  
  h_range <- max.window.size     # set the spatial lags over which to calculate codisp
  h1 = seq(-h_range, h_range, l = lx)  # x-axis values for codispersion graph (lags)
  h2 = seq(min(k), h_range, l = ly)    # y-axis values for codispersion graph (lags)
  
  
  MCodisp = matrix(0, ncol = ly, nrow = lx) # loop through the lags
  for(i in 1:lx)     # 'left-right' lags
  {
    for(j in 1:ly)   # 'up' lags
    {
      MCodisp[i,j] = Codisp.Kern(X, Y, c(h1[i], h2[j]), k)[4]; # calculate codisp
    }
  }
  Codispersion <- as.numeric(MCodisp) # save codisp object as output
  X <- rep(h1, length(h2))            # write out values for x-axis
  Y <- rep(h2, each = length(h1))       # write out values for y-axis
  graphing.data <- data.frame(X, Y, Codispersion) # graphing object
  
  # output the dataframe
  return(graphing.data) 
  
  }  # end function
