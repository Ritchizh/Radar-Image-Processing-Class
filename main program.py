import numpy as np
import matplotlib.pyplot as plt

import methods #import class with the processing methods
pr = methods.PROCESS()
    
pr.load_data('d_R.npy', 2) #filename, freq number

pr.mean_sub() #substract mean

plt.figure(1)
pr.plot_IQ()

plt.figure(2)
pr.focus(pr.d) #raw hologram at focus depth = d

#Processing:
pr.padding(25, 'constant')
pr.gauss_wind(200, 200)
pr.interpolate(2, 'linear') #interpol factor X = 2
pr.butter_filter('band-pass', D_h=6, D_l=300)


plt.figure(3)
pr.focus(pr.d) #processed hologram at focus depth = d
plt.title('Am processed')
plt.show()

