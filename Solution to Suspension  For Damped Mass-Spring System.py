#       Code to solve suspension of damped mass spring system using 4th order Runge-Kutta Method for ODE
#       Solved by:
#       Shashwata Sourav Roy Samya - Reg no: 2014132090
#       Md. Sadik Hasan - Reg no: 2015132029
#       Hrithik Barua - Reg no: 2015132048
#       Mohammad Rashik Zaman - Reg no: 2015132097

import matplotlib.pyplot as plt

t = 0                   # Initial value of time
tc = 0
to=0
y1 = [0.1016]           # used for under damping, at t = 0, y(0) = 0.1016
y1c = [0.1016]          # used for critical damping , at t = 0, y(0) = 0.1016
y1o = [0.1016]          # used for over damping , at t = 0, y(0) = 0.1016
y2 = [3.048]		# used for under damping , at t = 0, y'(0) = 10
y2c = [0.1016]          # used for cricital damping , at t = 0, y'(0) = 10
y2o = [0.1016]          # used for over damping , at t = 0, y'(0) = 10
time = [0]
timec = [0]
timeo = [0]
T = 0.01		# Change in time
N=10                   # No. of iterations

#function for f1(y,t)

def f1(y1, y2, t):
	return y2
def f1c(y1c, y2c, tc):
	return y2c
def f1o(y1o, y2o, to):
	return y2o

#function for f2(y,t)

def f2(y1, y2, t):
	return -4*y2-100*y1
def f2c(y1c, y2c, t):
	return -20*y2c-100*y1c
def f2o(y1o, y2o, t):
	return -28*y2o-100*y1o

def allC(g, y1, y2, t, T):                                      #Calculating c1, c2, c3 & c4 coefficient for using in Runge-Kutta Method
	c1 = T * g(y1, y2, t)
	c2 = T * g(y1+c1/2, y2+c1/2, t+T/2)
	c3 = T * g(y1+c2/2, y2+c2/2, t+T/2)
	c4 = T * g(y1+c3, y2+c3, t+T)
	
	return (c1+2*c2+2*c3+c4)/6

def allCc(g, y1c, y2c, tc, T):                                  #Calculating c1, c2, c3 & c4 coefficient for using in Runge-Kutta Method
	c1c = T * g(y1c, y2c, tc)
	c2c = T * g(y1c+c1c/2, y2c+c1c/2, tc+T/2)
	c3c = T * g(y1c+c2c/2, y2c+c2c/2, tc+T/2)
	c4c = T * g(y1c+c3c, y2c+c3c, tc+T)
	
	return (c1c+2*c2c+2*c3c+c4c)/6

def allCo(g, y1o, y2o, to, T):                                  #Calculating c1, c2, c3 & c4 coefficient for using in Runge-Kutta Method
	c1o = T * g(y1o, y2o, to)
	c2o = T * g(y1o+c1o/2, y2o+c1o/2, to+T/2)
	c3o = T * g(y1o+c2o/2, y2o+c2o/2, to+T/2)
	c4o = T * g(y1o+c3o, y2o+c3o, to+T)
	
	return (c1o+2*c2o+2*c3o+c4o)/6



while t<N:                                                              # loop for under damping, y(i+1) = y(i) + dy ; dy= (c1+2*c2+2*c3+c4)/6
	y1.append(y1[-1] + allC(f1, y1[-1], y2[-1], t, T))       
	y2.append(y2[-1] + allC(f2, y1[-1], y2[-1], t, T))
	time.append(t)
	t=t+T

while tc<N:                                                             # loop for critical damping, y(i+1) = y(i) + dy ; dy= (c1+2*c2+2*c3+c4)/6
	y1c.append(y1c[-1] + allCc(f1c, y1c[-1], y2c[-1], tc, T))       
	y2c.append(y2c[-1] + allCc(f2c, y1c[-1], y2c[-1], tc, T))
	timec.append(tc)
	tc=tc+T

while to<N:                                                             #loop for over damping, y(i+1) = y(i) + dy ; dy= (c1+2*c2+2*c3+c4)/6
	y1o.append(y1o[-1] + allCo(f1o, y1o[-1], y2o[-1], to, T))       
	y2o.append(y2o[-1] + allCo(f2o, y1o[-1], y2o[-1], to, T))
	timeo.append(to)
	to=to+T

plt.figure(1)
plt.plot(time,y1,'r', label="Under Damping")
plt.plot(timec,y1c,'b', label="Critical Damping")
plt.plot(timeo,y1o,'k', label="Over Damping")
plt.xlabel('time (sec)')
plt.ylabel('Dispalcement (m) ')
plt.title("Suspension of damped mass spring system")
plt.grid()
plt.legend()

plt.show()
