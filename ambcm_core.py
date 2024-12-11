import numpy as np
import scipy.integrate as integrate




def ambrcm_test(z, lambd): 
    try:
        iterator = iter(z)
    except TypeError:
        # not iterable
            res = np.exp(-z/np.pi*integrate.fixed_quad(lambda u: np.log(1-lambd *np.arctan(u)/u)/(1+(z**2)*(u**2)), 0, 60,n=200)[0])
    else:
        #iterable
        res=[]
        for z_el in z:
           res.append(np.exp(-z_el/np.pi*integrate.fixed_quad(lambda u: np.log(1-lambd *np.arctan(u)/u)/(1+(z_el**2)*(u**2)), 0, 60,n=200)[0]))
    return np.array(res)



def G(u, lamb):
    if u <= 2.5:
        res = (-u**5)*(np.log(1-lamb*np.arctan(u)/u) + lamb*np.arctan(u)/u+((lamb*np.arctan(u)/u)**2)/2
                       + ((lamb*np.arctan(u)/u)**3)/3 + ((lamb*np.arctan(u)/u)**4)/4)
    else:
        res = np.sum(u**5 * ((lamb*np.arctan(u)/u)**(k+5))/(k+5)
                  for k in range(0, 10001))

    return res


def F(u, lamb):

    A1 = lambda u: np.arctan(1/u)
    A3 = lambda u: A1(u)-1/u
    A5 = lambda u: A3(u)+1/(3*u**3)


    res = ((np.pi*lamb**2)/2)*u**3 * A3(u)-((np.pi*lamb**3)/2) * (u*A1(u))**2
    + (((np.pi**3)*(lamb**4))/8)*(u*A1(u))+lamb*(u**4)*A5(u) - \
        (lamb**2-((np.pi**2)*(lamb**3)/4))*u**2 * A3(u) + ((lamb**3)/3)*A1(u)*(u*A1(u))**2 -\
        (3*(np.pi**2)*(lamb**4)/8)*A1(u)*u*A1(u)+(np.pi*(lamb**4)/2)*u*A1(u)*A1(u)**2 - \
        (lamb**2)/2*A3(u)*u**3*A3(u)-(lamb**4)/4*(A1(u)**3)*u*A1(u)

    return res


def F1(u, lamb):

    A1 = lambda u: np.arctan(1/u)
    A3 = lambda u: A1(u)-1/u
    A5 = lambda u: A3(u)+1/(3*u**3)

    res = lamb**2*np.pi/2 * (u**3 * A3(u)-lamb*(u*A1(u))**2 + lamb**2*np.pi**2/4 * (u*A1(u))) +\
        lamb/u * (u**5*A5(u) - lamb * (1 - lamb * np.pi**2/4)*u**3*A3(u)+(lamb**2)/3*(u*A1(u)-9/8*lamb*np.pi**2)
                  * (u*A1(u))**2)+lamb**4*np.pi/(2*u**2)*(u*A1(u))**3-lamb**2/(2*u**3)*((u**3*A3(u))**2+lamb**2/2*(u*A1(u))**4)

    return res

def ambrcm(z, lamb):  # main function
    res = np.e**(-z/np.pi*integrate.quad(lambda u: np.log(1-lamb *
                 np.arctan(u)/u)/(1+(z**2)*(u**2)), 0, np.inf)[0])
    return res


def ambrcm_mid(z, lamb):  # z_1<=z <= 1

    res = 0
    um = 1
    b = [-lamb*np.pi/2, lamb*(1-lamb*np.pi**2/8),
         lamb**2 * np.pi / 2 * (1-lamb*np.pi**2/12), -lamb*(1./3+lamb/2-(lamb**2)*(np.pi**2)/4 + (lamb**3)*(np.pi**4)/64)]

    res += b[0] * 1/2 * np.log((1+z**2*um**2)/(z**2*um**2))  # b1f1
    res += b[1] * (1/um-z*(np.pi/2-np.arctan(z*um)))  # b2f2
    res += b[2] * (1/(2*um**2)-z**2/2 *
                   np.log((1+z**2*um**2)/(z**2*um**2)))  # b3f3
    res += b[3] * (1/(3*um**3)-z**2/um + z**3 *
                   (np.pi/2-np.arctan(z*um)))  # b4f4
    res += 1/(6*um**4)*integrate.quad(lambda t: (F1(um*np.e**(t/6), lamb)-G(um*np.e**(t/6), lamb)) *
                                      np.exp(-t)/(z**2*um**2+np.e**(-t/3)), 0, 700, epsabs=1.49e-08,limit=50)[0]


    if(lamb <=1+1e-6 and lamb >=1-1e-6):
        res+=conserv_flux(z,um)
    else:
        res += integrate.quad(lambda u: np.log(1-lamb*np.arctan(u)/u)/(1+(z**2)*(u**2)), 0, um)[0]


    #res+=integrate.quad(lambda u: (lamb*np.pi/(2*u))/(1+(z**2)*(u**2)), um, np.inf)[0]

    #res-=lamb*(-np.log(z)+1/2*np.log((1+z**2*um**2)/um**2))

    res *= -z/np.pi
    res = np.exp(res)

    return res



def ambrcm_small(z, lamb):  # small z
    if(z>1e-8):
        res = 0
        um = 1
        b = [-lamb*np.pi/2, lamb*(1-lamb*np.pi**2/8),
            lamb**2 * np.pi / 2 * (1-lamb*np.pi**2/12), -lamb*(1./3+lamb/2-(lamb**2)*(np.pi**2)/4 + (lamb**3)*(np.pi**4)/64)]

        res += b[0] * 1/2 * np.log((1+z**2*um**2)/(z**2*um**2))  # b1f1
        res += b[1] * (1/um-z*(np.pi/2-np.arctan(z*um)))  # b2f2
        res += b[2] * (1/(2*um**2)-z**2/2 *
                   np.log((1+z**2*um**2)/(z**2*um**2)))  # b3f3
        res += b[3] * (1/(3*um**3)-z**2/um + z**3 *
                    (np.pi/2-np.arctan(z*um)))  # b4f4
    

        res += z**4/6*integrate.quad(lambda t: (F1(np.e**(t/6)/z, lamb)-G(np.e**(t/6)/z, lamb)) *
                                      np.exp(-t)/(1+np.e**(-t/3)), 0, 700, epsabs=1.49e-08,limit=50)[0]
    
        res+=(1/um-z)*integrate.quad(lambda x: (F1(1/(z+(1/um-z)*x), lamb)-G(1/(z+(1/um-z)*x), lamb))*\
                                 (z+(1/um-z)*x)**3/(1+z**2 * 1/(z+(1/um-z)*x)**2) , 0, 1, epsabs=1.49e-08,limit=50)[0]

        if(lamb <=1+1e-6 and lamb >=1-1e-6):
            res+=conserv_flux(z,um)
        else:
            res += integrate.quad(lambda u: np.log(1-lamb*np.arctan(u)/u)/(1+(z**2)*(u**2)), 0, um)[0]

        res *= -z/np.pi
        res = np.exp(res)
    
    else:
        res=ambrcm(z,lamb)

    return res



def ambrcm_auto(z, lamb):
    if(z<=0.2):
        return ambrcm_small(z,lamb)
    else: 
        return ambrcm_mid(z,lamb)



def conserv_flux(z, um):
    res=0
    res+=np.log(um**2/3)*np.arctan(z*um)/z-2/z*np.sum((-1)**k*(z*um)**(2*k+1)/(2*k+1)**2
                  for k in range(0, 10001))
    res+=integrate.quad(lambda u: np.log(1-(1-3*(1-np.arctan(u)/u)/u**2))/(1+z**2*u**2), 0, um)[0]
    return res