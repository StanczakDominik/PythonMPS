import numpy as np
import matplotlib.pyplot as plt

def f(r):
    return np.exp(-r**2)

def AnglesToCartesian(angles, radiuses=1):
    return radiuses*np.cos(angles), radiuses*np.sin(angles)

Circle=2*np.pi
NumberOfRegions=8
AngleRange=Circle/NumberOfRegions

RegionAngles=np.linspace(0,Circle,NumberOfRegions+1)

#plot regions
x,y = AnglesToCartesian(RegionAngles)
xlines = np.vstack((np.zeros_like(x), x))
ylines = np.vstack((np.zeros_like(y), y))
plt.plot(xlines, ylines, "k--", label="Regions")

#plot points
RandomAngles = np.arange(0,Circle,AngleRange)+AngleRange*np.random.random(NumberOfRegions)
RandomRadii=2*np.random.random(NumberOfRegions)+0
xpoint, ypoint = AnglesToCartesian(RandomAngles, RandomRadii)
xpointlines = np.vstack((np.zeros_like(xpoint), xpoint))
ypointlines = np.vstack((np.zeros_like(ypoint), ypoint))
plt.plot(xpointlines, ypointlines, "bo-")

#get function values
FunctionValues=f(RandomRadii)


plt.show()
