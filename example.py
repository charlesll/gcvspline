import numpy as np
from matplotlib import pyplot as plt
from gcvspline import gcvspline, splderivative

# generating a dumb dataset for using the spline
x = np.linspace(0,100,1000)
y = 0.5 + x*30 - 0.001*x**2 + 4.5e-4 * x**3 + 1.64e-6*x**4 +500*np.cos(x/8)
noise =np.random.normal(scale=300,size=len(y))# adding noise , a lot!
y = y+ noise # adding noise , a lot!

# now using the spline, we are going to store results for different smoothing
# factors in an array y_fit
runs = 6
y_fit = np.zeros((np.size(x,0),runs)) # 10 runs
ese = np.sqrt(np.abs(y)) # we are not supposed to know the noise, so we guess (badly) that it's something like that...
smoothing = 0.1 # this is the interesting parameter that will allow to grow of shrink the error level, such that we have control over the spline smooting habilities
VAL = ese**2 # in mode 3, see the docs and Holtring article for further details.
for i in range(0,runs):
    c, wk, ier = gcvspline(x,y,smoothing*ese,VAL,splmode = 3) # gcvspl with mode 3 and smooth factor applied on data errors
    y_fit[:,i] = splderivative(x,x,c) 
    smoothing = smoothing*3 # we double the smoothing factor at each run

plt.figure()
plt.plot(x,y,'k.')
smoothing = 0.01
for i in range(0,runs): 
    if i > 0: 
        thickness = 3.0 
    else: 
        thickness = 0.5
    plt.plot(x,y_fit[:,i],'--',linewidth=thickness,color=[0+float(i)/float(runs-1.0),0,0],label=str(round(smoothing)))
    smoothing = smoothing*3 # we double the smoothing factor at each run

plt.xlabel('X values')
plt.ylabel('Y values')
plt.title('A simple example')
plt.legend(loc=2)
