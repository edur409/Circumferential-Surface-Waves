from scipy.optimize import least_squares
from scipy import optimize
from scipy.special import jv, jvp, yv, iv
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed   
import os 
 
def Sato_matrix62(omega, p, cl, ct, R):
    '''
    Equation 2.1 and 2.2 of Odaka and Usami (1978).
    "Some properties of spheroidal Modes of a Homogeneous Elastic sphere with Special 
    Reference to Radial Dependence of Displacement. Journal of Computational Physics, 29, 431-445, 1978
    omega = frequency
    p = angular wavenumber
    cl = longitudinal wave velocity of material
    ct = shear wave velocity of material
    '''
    # re-scaling for numerical convenience
    R = R*1000
    cl = cl/1000
    ct = ct/1000
    omega = omega/1e6
    # Code starts
    kl = omega/cl
    kt = omega/ct
    x = kl*R
    y = kt*R
    # Matrix form
    a = np.sqrt(2)*np.sqrt(np.pi)*(np.sqrt(2)*np.sqrt(np.pi)*(-p*jv(p, y) + y*jv(p - 1, y) -3*jv(p, y)/2)*np.sqrt(1/y)/R**2 + np.sqrt(2)*np.sqrt(np.pi)*(-2*y*jv(p - 1, y) + (2*p + 1)*jv(p, y) - (2*y**2 -(2*p - 1)*(2*p + 1))*jv(p, y)/4)*np.sqrt(1/y)/R**2)
    b = 2*np.pi*(p - 1/2)*(p + 1/2)*(-p*jv(p, x) + x*jv(p - 1, x) - 3*jv(p, x)/2)
    c = (-p*jv(p, y) + y*jv(p - 1, y) - 3*jv(p, y)/2)*np.sqrt(1/x)*np.sqrt(1/y)/R**4
    d = (-2*x*jv(p - 1, x) + (2*p + 1)*jv(p, x) - (2*y**2 - (2*p - 1)*(2*p + 1))*jv(p, x)/4)*np.sqrt(1/x)/R**2    
    M = np.array([[a, b], [c, d]])
    return M

class cmasw():
    def __init__(self, filename, xinital,  lb, ub, R, Print_r=True, Savefile= True):
        self.filename = filename
        self.R = R
        self.lb = lb
        self.ub = ub
        self.init = xinital 
        self.Print_r = Print_r
        self.Savefile = Savefile
   
        self.measure = np.loadtxt(filename) 
        self.omega = self.measure[:,0]*2*np.pi
        self.y_data =self.measure[:,1] 
        self.loss = []
        self.para = []
        self.jac = [] # save Jacovian matrix (added on 21/08/2024
        self.num = (len(lb) + 1) // 4 

     
    def aimfuc(self, x, omega, y_data):
        '''
        Calculate the loss function. Perform the forward model and save 
        outputs of a fitting iteration.
        '''
        vp = np.array(x[0])
        vs = np.array(x[1])
        yt =  self.velocity(omega, vp, vs) 
        loss = np.abs(yt - y_data)    
        
        if self.Print_r: 
            plt.clf()  
            plt.scatter(omega/(2*np.pi), y_data, label= 'measurement')    
            plt.plot(omega/(2*np.pi), yt, ':', label='Prediction')    
            self.loss.append(np.sqrt(np.sum(loss**2))) 
            plt.title('TRF Algorithm, Iteration Number: %s; L2 Norm: %s' %(len(self.loss), int(self.loss[-1])))
            plt.legend() 
            plt.xlabel('frequency (Hz)') 
            plt.ylabel('Phase Velocity (m/s)')
            plt.savefig('dispersion_update.png')  
            #print(np.sqrt(np.sum(loss ** 2)))
            if len(self.para) > 0:
                print('L2 Norm: %6.5f , Iteration Parameters: Vp: %5.2f, Vs: %5.2f' % (np.sqrt(np.sum(loss ** 2)), self.para[-1][0], self.para[-1][1]))                 
        if self.Savefile: 
            np.savetxt('prediction.txt', yt)
            self.para.append(x)
        return loss  


    def inverse(self):
        popt = least_squares(self.aimfuc, self.init, bounds=(self.lb, self.ub), method='trf', args=(self.omega, self.y_data))
        print("inversion is done") 
        r = popt.x
        jac = popt.jac # Added 21/08/2024
        self.jac.append(jac) # Added 21/08/2024
        Vp = np.array(r[0])
        Vs = np.array(r[1])
                
        print('Shear Velocity', Vs)         
        print('P wave velocity:', Vp)        

        if self.Savefile:
            np.savetxt('loss.csv', self.loss)
            np.savetxt('parameters.csv', self.para)
            
    def velocity(self, omega, vp, vs):
        '''
        Function to perform the forward model.
        '''
        def opt(i):
            def fun(k): 
                matrix1 = Sato_matrix62(i, k , vp, vs, self.R) #Viktorov_matrix(i, k , vp, vs, self.R) 
                matrix2 = matrix1 #np.real(matrix1)#/10**11
                sign, logdet = np.linalg.slogdet(matrix2)
                return sign * np.exp(logdet)
            c_test = 0.5*np.min(self.y_data) # 90% of the min. vel. of the dispersion curve is a sensible starting choice
            incre = self.R*i / c_test
            root = 0.00001
            for j in range(10**6):
                past = incre
                val1 = fun(incre)
                incre =  incre - 0.01
                val2 = fun(incre)
                if (np.real(val1) * np.real(val2) <= 0):
                    root =  optimize.brentq(fun,incre, past)      
                    break 
            return (self.R*i/root)     #give one value at a frequency
    
        def final(n):
            y = opt(omega[n])
            return y   
        
        z2 = Parallel(n_jobs=-1)(delayed(final)(i) for i in range(len(omega))) 
        return z2

        
 
 

 
