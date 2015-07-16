import numpy as np
import matplotlib.pyplot as plt

output=open("./output.dat", 'w+')

print("Funkcja: e^(-x^2-y)")
print("Funkcja: e^(-x^2-y)", file=output)
def f(x,y):
    return np.exp(-x**2-y)
def fprime(x,y):
    return -2*x*f(x,y), -f(x,y)

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
maxradius=1
minradius=0
RandomAngles = np.arange(0,Circle,AngleRange)+AngleRange*np.random.random(NumberOfRegions)
RandomRadii=(maxradius-minradius)*np.random.random(NumberOfRegions)+minradius
xpoint, ypoint = AnglesToCartesian(RandomAngles, RandomRadii)
xpointlines = np.vstack((np.zeros_like(xpoint), xpoint))
ypointlines = np.vstack((np.zeros_like(ypoint), ypoint))
plt.plot(xpointlines, ypointlines, "bo-", label="Wylosowane punkty")

#get function values
R = np.sum(RandomRadii)/NumberOfRegions #scaling radius
alpha_R=R/RandomRadii #scaling

xpointScaled = xpoint*alpha_R
ypointScaled = ypoint*alpha_R

plt.plot(xpointScaled, ypointScaled, "ro", label="Przeskalowane punkty")

FunctionValuesScaled=f(xpointScaled, ypointScaled)
FunctionAtZero=f(0,0)


alpha1=np.sum(xpointScaled**2)
alpha2=np.sum(xpointScaled*ypointScaled)
alpha3=np.sum(ypointScaled**2)
delta1=np.sum((FunctionValuesScaled-FunctionAtZero)*xpointScaled)
delta2=np.sum((FunctionValuesScaled-FunctionAtZero)*ypointScaled)

a = -(delta1*alpha3-delta2*alpha2)/(alpha2**2-alpha1*alpha3) #d/dx
b = -(delta2*alpha1-delta1*alpha2)/(alpha2**2-alpha1*alpha3) #d/dy

print("Wyniki analityczne")
print("%.4f, %.4f" %fprime(0,0))

print("Wartosci pochodnych ze skalowaniem")
print("%.4f, %.4f" %(a,b))
print("Roznice")
print("%.4f, %.4f" %(abs(a-fprime(0,0)[0]), abs(b-fprime(0,0)[1])))

print("Wyniki analityczne", file=output)
print("%.4f, %.4f" %fprime(0,0), file=output)

print("Wartosci pochodnych ze skalowaniem", file=output)
print("%.4f, %.4f" %(a,b), file=output)
print("Roznice", file=output)
print("%.4f, %.4f" %(abs(a-fprime(0,0)[0]), abs(b-fprime(0,0)[1])), file=output)

FunctionValues=f(xpoint, ypoint)

alpha1=np.sum(xpoint**2)
alpha2=np.sum(xpoint*ypoint)
alpha3=np.sum(ypoint**2)
delta1=np.sum((FunctionValues-FunctionAtZero)*xpoint)
delta2=np.sum((FunctionValues-FunctionAtZero)*ypoint)

a = -(delta1*alpha3-delta2*alpha2)/(alpha2**2-alpha1*alpha3) #d/dx
b = -(delta2*alpha1-delta1*alpha2)/(alpha2**2-alpha1*alpha3) #d/dy

print("Wartosci pochodnych bez skalowania")
print("%.4f, %.4f" %(a,b))
print("Roznice")
print("%.4f, %.4f" %(abs(a-fprime(0,0)[0]), abs(b-fprime(0,0)[1])))


print("Wartosci pochodnych bez skalowania", file=output)
print("%.4f, %.4f" %(a,b), file=output)
print("Roznice", file=output)
print("%.4f, %.4f" %(abs(a-fprime(0,0)[0]), abs(b-fprime(0,0)[1])), file=output)

plt.xlim(-maxradius, maxradius)
plt.ylim(-maxradius, maxradius)
plt.grid()
plt.savefig("plot.png")
plt.show()
plt.clf()
