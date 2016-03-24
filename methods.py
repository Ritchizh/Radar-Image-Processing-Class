import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.signal

class PROCESS:
        def __init__(self):
                self.data       = []
                self.E          = []
                self.f          = []
                self.f1         = []
                self.f2         = []
                self.nf         = []
                self.step_m     = []
                self.d          = []
                self.frq_num    = 0
                self.frq_val    = []

        def load_data(self, file_name, frq_num): #loads data from file into variables
                file            = open(file_name, 'rb')
                self.data       = np.load(file)
                self.f1         = np.load(file)
                self.f2         = np.load(file)
                self.nf         = np.load(file)
                self.step_m     = np.load(file)
                self.d          = np.load(file)
                file.close()
                
                self.frq_num    = frq_num #freq number
                ff              = np.linspace(self.f1, self.f2, self.nf) 
                self.frq_val    = ff[self.frq_num] #freq value
                
                self.E          = self.data[frq_num, :, :]
                self.E          = np.ma.conjugate(self.E)
                print('d = ', self.d)
                
                print('f = ',  self.frq_val)

#~~MEAN SUBSTRACTION~~
        def mean_sub(self):
                Background = np.average(self.E)
                self.E = self.E - Background

        def mean_row_sub(self): #substracts an averaged data row from each row
                av_profile = np.zeros((1,self.E.shape[1]), dtype=complex)
                for j in range (0,self.E.shape[1]):
                        av_profile[0][j] = np.average(self.E[:,j])
                for i in range (0, self.E.shape[0]):
                        self.E[i,:] = self.E[i,:] - av_profile

        def mean_column_sub(self): #substracts an averaged data column from each column
                av_profile = np.zeros((self.E.shape[0],1), dtype=complex)
                for i in range (0,self.E.shape[0]):
                        av_profile[i][0] = np.average(self.E[i,:])
                for j in range (0, self.E.shape[1]):
                        self.E[:,j] = self.E[:,j] - av_profile

#~~PLOTTING~~
# Don't forget writing 'plt.show()' in the end of your script.
        def plot_IQ(self): #plots I and Q components
                plt.subplot(121)
                plt.imshow(self.E.real, cm.gray)
                plt.title('In-phase signal')
                plt.tight_layout()
                plt.axis('off')

                plt.subplot(122)
                plt.imshow(self.E.imag, cm.gray)
                plt.title('Quadrature signal')
                plt.tight_layout()
                plt.axis('off')
                

        def focus(self, d): #focuses data at passed depth = d
                F       = np.fft.fft2(self.E) #FFT
                F       = np.fft.fftshift(F)
        
                fx      = np.linspace(-1/(2*self.step_m[1]), 1/(2*self.step_m[1]), self.E.shape[1]) #Spacial frequences (1/dx)
                fy      = np.linspace(-1/(2*self.step_m[0]), 1/(2*self.step_m[0]), self.E.shape[0])
                Fx, Fy  = np.meshgrid(fx, fy)
                k       = (2*np.pi*self.frq_val)/(3*10**8)
        
                TransMx = np.exp(np.lib.scimath.sqrt(4*k**2 - (2*np.pi*Fx)**2 - (2*np.pi*Fy)**2)*self.d*1j)
                S       = F*TransMx
                S       = np.fft.fftshift(S)
                self.f  = np.fft.ifft2(S)
                plt.imshow(abs(self.f), cm.gray) #plots abs
                plt.title('Am')
                plt.tight_layout()
                plt.axis('off')

#~~RADAR_IMAGE PROCESSING~~
        def padding(self, width_div = 15, mode = 'constant'): #pad width is data width/division coeff
                if mode == 'constant':
                        self.E = np.lib.pad(self.E, (self.E.shape[0]//width_div,), 'constant', constant_values=(0)) #with zeros
                if mode == 'mean':
                        self.E = np.lib.pad(self.E, (self.E.shape[1]//width_div,), 'mean') #with mean value

        def gauss_wind(self, stdx = 600, stdy = 600): #window widths in x, y directions
                gausx = scipy.signal.gaussian(self.E.shape[1], std = stdx)
                gausy = scipy.signal.gaussian(self.E.shape[0], std = stdy)
                for i in range(0, self.E.shape[0]-1):
                        for j in range (0, self.E.shape[1]-1):
                                self.E[i,:] = self.E[i,:]*gausx
                                self.E[:,j] = self.E[:,j]*gausy

        def interpolate(self, X, mode = 'linear'): #interpol factor X
                e0              = (self.E.shape[0]-1)*self.step_m[0]
                e1              = (self.E.shape[1]-1)*self.step_m[1]
                points1         = np.linspace(0, e0, self.E.shape[0])
                points2         = np.linspace(0, e1, self.E.shape[1])
                grid_x, grid_y  = np.mgrid[0:e0:(self.E.shape[0]*X*1j), 0:e1:(self.E.shape[1]*X*1j)]
                self.E               = scipy.interpolate.interpn((points1,points2), self.E, (grid_x, grid_y), method = mode) 
                #mode: linear, nearest, splinef2d
                self.step_m     = self.step_m/X


        def butter_filter(self, mode, n=1, D_h=0, D_l=0): #n is filter order; insert D for a high/low pass filter; D_h and D_l for a band-pass
                fx      = np.linspace(-1/(2*self.step_m[1]), 1/(2*self.step_m[1]), self.E.shape[1])
                fy      = np.linspace(-1/(2*self.step_m[0]), 1/(2*self.step_m[0]), self.E.shape[0])
                Fx, Fy  = np.meshgrid(fx, fy)
                F       = np.fft.fft2(self.E)
                F       = np.fft.fftshift(F)
                if mode == 'high-pass':
                        Butter = 1/(1 + (D_h/((Fx**2 + Fy**2)**0.5))**(2*n))
                if mode == 'low-pass':
                        Butter = 1/(1 + (((Fx**2 + Fy**2)**0.5)/D_l)**(2*n))
                if mode =='band-pass':
                        Butter_h = 1/(1 + (D_h/((Fx**2 + Fy**2)**0.5))**(2*n))
                        Butter_l = 1/(1 + (((Fx**2 + Fy**2)**0.5)/D_l)**(2*n))
                        Butter = Butter_h*Butter_l
                Butter  = np.nan_to_num(Butter)       
                y       = F*Butter
                y1      = np.fft.fftshift(y)
                self.E      = np.fft.ifft2(y1)
