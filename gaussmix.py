import math 

class gaussmix:
    """guassian mixture for pacbio dme"""

    ################################
    def __init__( self, mixfrac, mean, cov, baseline):
        self.mixfrac = mixfrac
        self.mean = mean
        self.cov = cov
        self.baseline = baseline

    ################################
    def logcompprob(self, x,y):
        """return of the log probability of the 5 components for the input point x,y (green, red) [-,_A,_C,_G,_T]
        /home/UNIXHOME/mbrown/mbrown/workspace2015Q4/sequel-mixturemodel/README_covariance-ellipse.html
        https://en.wikipedia.org/wiki/Multivariate_normal_distribution
        http://mathworld.wolfram.com/BivariateNormalDistribution.html

        log p(x,y) = -log(2*pi*sx*sy*sqrt(1-rho^2)) + (-0.5* z/(1-rho^2))
        z = (x-mx)^2/sx^2 - 2*rho*(x-mx)*(y-my)/(sx*sy) + (y-my)^2/sy^2
        rho = cov12/(sx*sy)
        
        self.mean[component] = [mx,my]
        self.cov[component] = [sx^2,sy^2,cov12]

        install.packages("mvtnorm",repos="http://cran.us.r-project.org")
        library(mvtnorm)
        dmvnorm(x=c(1,1), mean=c(1,1), sigma= cbind(c(1,0),c(0,1)), log=T)
        [1] -1.837877
        dmvnorm(x=c(40,50), mean=c(199.76504517,96.82637024  ), sigma= cbind(c(955.24627686,261.15966797),c(261.15966797,531.82269287)), log=T)
        [1] -21.70605

        test = gaussmix.gaussmix( [1.0], [[1.0,1.0]], [[1.0,1.0,0.0]])
        test.logcompprob(1,1)
        [-1.8378762217451237]

        test = gaussmix.gaussmix( [1.0], [[199.76504517,   96.82637024]], [[955.24627686,   531.82269287,   261.15966797]])
        test.logcompprob(40,50)
        [-21.706051802815622]

        """
        # subtract off baseline
        x = x - self.baseline[0]
        y = y - self.baseline[1]

        lp = []
        for cc in range(len(self.mixfrac)):
            (mx,my) = self.mean[cc]
            sx = math.sqrt(self.cov[cc][0])
            sy = math.sqrt(self.cov[cc][1])
            cov12 = self.cov[cc][2]
            rho = cov12/(sx*sy)
            z = ((x-mx)**2)/(sx**2) - (2*rho*(x-mx)*(y-my))/(sx*sy) + ((y-my)**2)/(sy**2)
            lp.append( -math.log(2*3.14159*sx*sy*math.sqrt(1-rho**2)) + ((-0.5*z)/(1-rho**2)) + math.log(self.mixfrac[cc]) )
        return(lp)

    ################################
    def qv(self, logprobs ):
        "return 10*log10(p1/p2) where p1 is the most likely and p2 is the 2nd"
        ss = sorted(logprobs, key= lambda x: -x) # decreasing so biggest is at 0
        return(10.0*(ss[0] - ss[1])/math.log(10.0))

        
