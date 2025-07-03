# A simple energy minimization program that uses steepest descent 
# and a force field to minimize the energy of water in internal coordinates
# Written by Jan H. Jensen, 2013, released in the the GNU GPL license

kOH = 50.0
rOHe = 0.95
kHOH = 50.0
thetaHOHe = 104.5

c = 0.005
n_steps = 20

#starting geometry
rOH = 10.0
thetaHOH = 180.0

def Eandg(rOH,thetaHOH):
    E = 2*kOH*(rOH-rOHe)**2 + kHOH*(thetaHOH-thetaHOHe)**2
    grOH = 2*kOH*(rOH-rOHe)
    gthetaHOH = 2*kHOH*(thetaHOH-thetaHOHe)
    return (E,grOH,gthetaHOH)


for i in range(n_steps):
    (E,grOH,gthetaHOH) = Eandg(rOH,thetaHOH)
    if (abs(grOH) >0.001/c or abs(gthetaHOH) > 0.01/c ):
        rOH = rOH - c*grOH
        thetaHOH = thetaHOH - c*gthetaHOH

if (abs(grOH) >0.001/c or abs(gthetaHOH) > 0.01/c ):
    print "not converged"
else:
    print "converged"

print E,rOH,thetaHOH