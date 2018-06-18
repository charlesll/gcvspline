********************
 Available Functions
********************


-------------------
 Function gcvspline
-------------------

Call it as:

gcvspline(x,y,ese,VAL,**options)
    
Input:
	
    x (1D numpy array) are the independent variables
	
    y (2D array as (len(x),1)) are the observations
    (we assume here that you want to use this 
    spline only on 1 dataset... see gcvspl.f if not)
	
    ese are the errors on y
	
    VAL (double) depends on the mode, see gcvspl.f for details (=ese**2 for mode3)

Options:
	
    NC is the number of output spline coefficients, NC >= len(y), default: NC = len(y)
		
    splorder is the half order of the required B-splines. The values 
    splorder = 1,2,3,4 correspond to linear, cubic, quintic, and 
    heptic splines, respectively.
    default: splorder = 2 (cubic) 
		
    splmode is the Optimization mode switch, default: splmode = 2 (General Cross Validated)
	
	splmode = 1: Prior given value for p in VAL (VAL.ge.ZERO). This is the fastest use of GCVSPL, since no iteration is performed in p.
	splmode = 2: Generalized cross validation.
    splmode = 3: True predicted mean-squared error, with prior given variance in VAL.
    splmode = 4: Prior given number of degrees of freedom in VAL (ZERO.le.VAL.le.N-M).
    splmode  < 0: It is assumed that the contents of X, W, M, N, and WK have not been modified since the previous invocation of GCVSPL. If MD < -1, WK(4) is used as an initial estimate for the smoothing parameter p.  At the first call to GCVSPL, MD must be > 0. Other values for MD, and inappropriate values for VAL will result in an error condition, or cause a default value for VAL to be selected.  After return from MD.ne.1, the same number of degrees of freedom can be obtained, for identical weight factors and knot positions, by selecting MD=1, and by copying the value of p from WK(4) into VAL. In this way, no iterative optimization is required when processing other data in Y.     
    
Output:

    c: the spline array
	
    wk: work vector, see gcvspl.f
	
    ier: error parameter 
    ier = 0: Normal exit 
    ier = 1: M.le.0 .or. N.lt.2*M
    ier = 2: Knot sequence is not strictly increasing, or some weight factor is not positive.
    ier = 3: Wrong mode parameter or value.
	
see gcvspl.f for more informations

-----------------------
 Function splderivative
-----------------------

Call it as:

splderivative(xfull,xparse,cparse,**options):

Inputs:
    
    xfull (1D array) contains the entire x range where the spline has to be evaluated
    xparse (1D array) contains the x values of interpolation regions; WARNING: xparse[0] <= xfull[0] <= xparse[n] 
    cparse (1D array) is the evaluated spline coefficients returned by gcvspl for xparse
    
Options:
    
    splineorder (integer): is the spline order, default: splineorder = 2 (cubic)
		
    L (integer): see gcvspl.f for details, default: L = 1
	
    IDER: the Derivative order required, with 0.le.IDER and IDER.le.2*M. If IDER.eq.0, the function value is returned; otherwise, the IDER-th derivative of the spline is returned.
        
see gcvspl.f for more informations
 