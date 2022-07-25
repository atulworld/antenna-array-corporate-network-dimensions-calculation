#!/usr/bin/env python
# coding: utf-8

# In[50]:


import numpy
import math
initial_population=10
number_of_parents=(initial_population/2)
number_of_parents=int(number_of_parents)

g=(initial_population,2)
num_vari=2;
er=float(input("input the er: "))
h=float(input("input the h: "))
fr=float(input("input the fr in ghz: "))


c=3*10**11;
frq=fr*10**9
num_parents=number_of_parents;
number_of_veri=2;

X=numpy.random.uniform(5,50,size=g)
print("The initial population is   =", X)

num_generations = 500
for generation in range(num_generations):
    print("Generation : ", generation)
    L = numpy.delete(X, 1, 1)
    #print (a)
    W = numpy.delete(X, 0, 1)
    #print (b)
    kulk=(1+((12*h)/W))
    #print (kulk)
    kulk1= (numpy.sqrt(kulk))
    #print(kulk1)
    ef=((er+1)/2)+(((er-1)/2)*(1/(kulk1)))
    #print(ef)
    e1=ef+0.3
    e2=ef-0.258
    #print (e1)
    #print(e2)
    w1=((W/h)+0.264)
    #print (w1)
    w2=((W/h)+0.813)
    #print(w2)
    dw=(0.412*h)*((e1*w1)/(e2*w2))
    dw=2*dw
    #print(dw)
    m=2*(L+dw)
    #print(m)
    n=(numpy.sqrt(ef))
    #print(n)
    fitness=c/(m*n)
    print ("the fitness wrt pop is =", fitness)
    
    topp=numpy.empty((initial_population,1))
    mpty=topp*0
    #print(mpty)
    topp1=mpty+(frq)
    #print(topp1)
    error=((topp1-fitness)**2)/2*(10**-9)
    print ("the deviation in req freq wrt population  =", error)
    def select_fit_parents(X, error, numparent):
        parents = numpy.empty((num_parents, number_of_veri))
        for i in range(num_parents):
            min_error_idx = numpy.where(error == numpy.min(error))
            min_error_idx = min_error_idx[0]
            parents[i, :] = X[min_error_idx,:]
            error[min_error_idx]=(9*10**22)
        return parents
    parents=select_fit_parents(X,error,num_parents)
    print ("The fit veri values from the population are   =", parents)
    print("modified error values  =", error)
    offspring_size=(g[0]-parents.shape[0],num_vari)
    
    def crossover(parents, offspringsize):
        offspring = numpy.empty(offspring_size)
        crossover_point = numpy.uint8(offspring_size[1]/2)
        for k in range(offspring_size[0]):
            parent1_idx = k%parents.shape[0]
            parent2_idx = (k+1)%parents.shape[0]
            offspring[k, 0:crossover_point] = parents[parent1_idx, 0:crossover_point]
            offspring[k, crossover_point:] = parents[parent2_idx, crossover_point:]
        return(offspring)
    offspring_crossover=crossover(parents,offspring_size)
    #print ("The offsprings after crossover are  =", offspring_crossover)
    
    def mutation(offspring_crossover):
        for idx in range(offspring_crossover.shape[0]):
            random_value = numpy.random.uniform(1.0, 4.0, 1)
            offspring_crossover[idx, 1] = offspring_crossover[idx, 1] + random_value
        return offspring_crossover
    offspring_mutation = mutation(offspring_crossover)
    #print ("The offsprings after mutation are  =",  offspring_mutation )
    
    X[0:parents.shape[0], :] = parents
    X[parents.shape[0]:, :] = offspring_mutation
    print("The new population is   =", X)
    
    L = numpy.delete(X, 1, 1)
    #print (a)
    W = numpy.delete(X, 0, 1)
    #print (b)
    kulk=(1+((12*h)/W)) 
    #print (kulk)
    kulk1= (numpy.sqrt(kulk))
    #print(kulk1)
    ef=((er+1)/2)+(((er-1)/2)*(1/(kulk1)))
    #print(ef)
    e1=ef+0.3
    e2=ef-0.258
    #print (e1)
    #print(e2)
    w1=((W/h)+0.264)
    #print (w1)
    w2=((W/h)+0.813)
    #print(w2)
    dw=(0.412*h)*((e1*w1)/(e2*w2))
    dw=2*dw
    #print(dw)
    m=2*(L+dw)
    #print(m)
    n=(numpy.sqrt(ef))
    #print(n)
    fitness=c/(m*n)
    #print (fitness)
    topp=numpy.empty((initial_population,1))
    mpty=topp*0
    #print(mpty)
    topp1=mpty+(frq)
    #print(topp1)
    error=((topp1-fitness)**2)/2*(10**-9)
    #error=((topp1-fitness)**2)/2
    #print ("the deviation in req freq wrt new population  =", error)
    print("Best result : ", numpy.min(error))

best_match_idx = numpy.where(error == numpy.min(error))
print(best_match_idx)

print("Best solution : ", X[best_match_idx, :])
at=fitness
ul=topp1
ra=(at-ul)/10**9
print("For the above best Best solution, deviation from the req fre will be  : ", ra[best_match_idx])
print (fitness)
print (topp1)
Lpatch=(X[best_match_idx])
list1 = Lpatch.tolist()

#width cal
d1=c/(2*frq)
d2=(numpy.sqrt(2/(er+1)))
wpatch=d1*d2

#effective er cal

kulk=(1+((10*h)/wpatch)) 
#print (kulk)
kulk1= (numpy.sqrt(kulk))
    #print(kulk1)
ef=((er+1)/2)+(((er-1)/2)*(1/(kulk1)))


print("the relative permitivity : ",er)
print("the thickness of substrate : ",h)
print("the frequency : ",fr)
print("the effective permitivity : ",ef)
print("the width of patch : ",wpatch)
print("the length of patch : ",list1)

# cal of l and w of feed, 5oohmline, 100ohmline and qwt
import math
import numpy
import scipy
from scipy.special import*
from scipy.integrate import quad
pi=3.14
la=((3*10**8)/frq)*1000;
k=(2*pi)/la;
print(k)
x=k*(wpatch);
(si,ci) = sici(x)
i1=-2+math.cos(x)+(x*si)+(math.sin(x)/x);
g1=i1/(120*pi*pi);

L=list1[0]

def integrand(th, x, k,L):
    return (((math.sin((x/2)*math.cos(th))/math.cos(th))**2)*(numpy.i0((k*L*math.sin(th))))*(math.sin(th))**3)
a1=quad(integrand, 0, pi, args=(x,k,L))[0]

g12=a1/(120*pi*pi);
r_in=1/(2*(g1+g12))

print ("the input impedance of the patch is: ", r_in)

Lg_min=(6*h)+L;  
print ("for single patch min ground length is: ",Lg_min)

Wg_min=(6*h)+wpatch;
print ("for single patch min ground width is: ",Wg_min)

z=50
yo=(numpy.sqrt(z/r_in))
inset=(L/pi)*(math.acos(yo))

print ("the inset feed length is: ",inset)

xo=numpy.sqrt(er)
B=(377*pi)/(2*z*xo); 
m1=(2*B)-1;
m=math.log(m1);
n1=B-1;
n=math.log(n1);

#based on Microwave Engineering book by David M. Pozar Fourth Edition
Wfeed=(2*h/pi)*(B-1-m+(((er-1)/(2*er))*(n+((0.39*0.61)/er))));
print("the inset feed width is: ",Wfeed)

#based on journal by M A Matin (2010)
g=((3*10**8)*(4.65*10**-9)*10**12)/(numpy.sqrt(2*ef)*frq)*10**-3;
print("The gap of the feed line is: ",g)

#for 100 ohm , qwt
z1=100
z2=70.71
B1=(60*pi*pi)/(z1*xo);
m1=(2*B1)-1;
m=math.log(m1);
n1=B1-1;
n=math.log(n1);
W100ohm=(2*h/pi)*(B1-1-m+(((er-1)/(2*er))*(n+((0.39*0.61)/er))));
print("The width of the 100ohm line is: ",W100ohm)


B2=(60*pi*pi)/(z2*xo);
m1=(2*B2)-1;
m=math.log(m1);
n1=B2-1;
n=math.log(n1);
Wqwt=(2*h/pi)*(B2-1-m+(((er-1)/(2*er))*(n+((0.39*0.61)/er))));
print("The width of the qwt is: ",Wqwt)

lam_g=la/(numpy.sqrt(ef))
lqwt=lam_g/4;
print("The length of the qwt is: ",lqwt)


# In[ ]:




