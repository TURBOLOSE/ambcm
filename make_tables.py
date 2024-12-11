import numpy as np
import pandas as pd

from ambcm_core import *


z = np.arange(0, 1.05, 0.05)

phi7 = np.empty(0)  # ambrcm at lambda=0.7
phi8 = np.empty(0)  # ambrcm at lambda=0.8
phi10 = np.empty(0)  # ambrcm at lambda=1.

phi7_p = np.empty(0)  # ambrcm at lambda=0.7
phi8_p = np.empty(0)  # ambrcm at lambda=0.8
phi10_p = np.empty(0)  # ambrcm at lambda=1.


for i in z:
    phi7 = np.append(phi7, ambrcm(i, 0.7))
    phi8 = np.append(phi8, ambrcm(i, 0.8))
    phi10 = np.append(phi10, ambrcm(i, 1.))
    phi7_p = np.append(phi7_p, ambrcm_auto(i, 0.7))
    phi8_p = np.append(phi8_p, ambrcm_auto(i, 0.8))
    phi10_p = np.append(phi10_p, ambrcm_auto(i, 1.))


err7 = np.empty(0)
err8 = np.empty(0)
err10 = np.empty(0)

for i in range(len(z)-1):  # check for errors with equation
    err7 = np.append(err7, phi7[1+i]-(1+0.7/2 * z[1+i]*phi7[1+i] *
                     integrate.quad(lambda y: ambrcm(y, 0.7)/(z[1+i]+y), 0, 1)[0]))
    err8 = np.append(err8, phi8[1+i]-(1+0.8/2 * z[1+i]*phi8[1+i] *
                     integrate.quad(lambda y: ambrcm(y, 0.8)/(z[1+i]+y), 0, 1)[0]))
    err10 = np.append(err10, phi10[1+i]-(1+1./2 * z[1+i]*phi10[1+i]
                      * integrate.quad(lambda y: ambrcm(y, 1.)/(z[1+i]+y), 0, 1)[0]))


# writing files

d = {'(lambda=0.7) z': z, 'H(z)': phi7, 'H(z)_1 ': phi8}
df = pd.DataFrame(data=d)

np.savetxt('zfi0708_my_main.txt', df.values, fmt='%1.12f',
           header='lambda=0.7 and 0.8 \n z           H(z)       H   (z)')


d = {'z': z, 'H(z)': phi10}
df = pd.DataFrame(data=d)

np.savetxt('zfi10_my_main.txt', df.values, fmt='%1.12f',
           header='lambda=0.7 and 0.8 \n z           H(z)')


d = {'err7': err7, 'err8': err8, 'err10': err10}
df = pd.DataFrame(data=d)

np.savetxt('errors.txt', df.values, fmt='%e',
           header='equation errors \n lambda=0.7  lambda=0.8  lambda=1.0')




# writing files
d = {'(lambda=0.7) z': z, 'H(z)': phi7_p, 'H(z)_1 ': phi8_p}
df = pd.DataFrame(data=d)

np.savetxt('zfi0708_my_prec.txt', df.values, fmt='%1.12f',
           header='lambda=0.7 and 0.8 \n z           H(z)       H   (z)')


d = {'z': z, 'H(z)': phi10_p}
df = pd.DataFrame(data=d)

np.savetxt('zfi10_my_prec.txt', df.values, fmt='%1.12f',
           header='lambda=0.7 and 0.8 \n z           H(z)')

