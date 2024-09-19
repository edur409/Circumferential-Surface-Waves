import numpy as np
from scipy.linalg import eigh, eig
from scipy.special import factorial, factorial2

# def calculate_elastic_tensor_full(Vp, Vs, rho):
#     # Calculate the bulk modulus K and shear modulus mu
#     K = rho * Vp**2
#     mu = rho * Vs**2
    
#     # Calculate Lamé parameters lambda and mu
#     lam = K - (2 / 3) * mu
    
#     # Initialize the 6x6 elastic tensor in Voigt notation
#     C_voigt = np.zeros((6, 6))
    
#     # Fill the elastic tensor components using Voigt notation
#     C_voigt[0, 0] = lam + 2 * mu # C_1111
#     C_voigt[1, 1] = lam + 2 * mu # C_2222
#     C_voigt[2, 2] = lam + 2 * mu # C_3333
#     C_voigt[0, 1] = lam # C_1122
#     C_voigt[0, 2] = lam # C_1133
#     C_voigt[1, 0] = lam # C_2211
#     C_voigt[1, 2] = lam # C_2233
#     C_voigt[2, 0] = lam # C_3311
#     C_voigt[2, 1] = lam # C_3322
#     C_voigt[3, 3] = mu # C_1212
#     C_voigt[4, 4] = mu # C_2323
#     C_voigt[5, 5] = mu # C_3131
    
#     # Create a 4th-order tensor
#     C_full = np.zeros((3, 3, 3, 3))
    
#     # Fill the 4th-order tensor based on the Voigt notation
#     for i in range(3):
#         for j in range(3):
#             for k in range(3):
#                 for l in range(3):
#                     C_full[i, j, k, l] = C_voigt[i, j] if k == l else 0
#                     C_full[i, j, k, l] += C_voigt[k, l] if i == j else 0
#                     C_full[i, j, k, l] += C_voigt[i, k] if j == l else 0
#                     C_full[i, j, k, l] += C_voigt[j, l] if i == k else 0
    
#     return C_full

# def V_l(lam, mu, rho):
#     return np.sqrt((lam + 2*mu)/rho)

# def V_t(mu, rho):
#     return np.sqrt(mu/rho)

# lam = 0.54*1e12
# mu = 0.27*1e12 # erg/cm^3
# rho = 2.7065 # Density in g/cm^3

# # Velocities:
# Vp = V_l(lam, mu, rho) # Compressional velocity in m/s
# Vs = V_t(mu, rho) # Shear velocity in m/s

# Constants
TWOPI = 2 * np.pi
RHO = 1.0

# Elastic tensor Cijkl for the standard isotropic material
data = [
    3., 0., 0., 0., 1., 0., 0., 0., 1., 
    0., 1., 0. ,1., 0., 0., 0., 0., 0.,
    0., 0., 1., 0., 0., 0., 1., 0., 0., 
    0., 1., 0., 1., 0., 0., 0., 0., 0., 
    1., 0., 0., 0., 3., 0., 0., 0., 1.,
    0., 0., 0., 0., 0., 1., 0., 1., 0.,
    0., 0., 1., 0., 0., 0., 1., 0., 0.,
    0., 0., 0., 0., 0., 1., 0., 1., 0.,
    1., 0., 0., 0., 1., 0., 0., 0., 3          
]

# Convert to numpy array and reshape
C = np.array(data).reshape((3, 3, 3, 3), order='F')#*1e12
# Calculate the full 4th-order elastic tensor
# C = calculate_elastic_tensor_full(Vp, Vs, rho)

# Function to compute f(p, q, r) for the corner prism
def F(ip, iq, ir):
    A = B = C = 1.0
    return (A**(ip + 1) * B**(iq + 1) * C**(ir + 1) *
            factorial(ip) * factorial(iq) * factorial(ir) /
            factorial(ip + iq + ir + 3))

#def main():
# Read NN from the user
NN = int(input("PLEASE INPUT NN: "))

# Initialize arrays
MAX_SIZE = int((NN + 1)*(NN + 2)*(NN + 3)/2) # R in the original paper
GAMMA = np.zeros((MAX_SIZE, MAX_SIZE))
E = np.zeros((MAX_SIZE, MAX_SIZE))
W = np.zeros(MAX_SIZE)
FV = np.zeros(MAX_SIZE)
FW = np.zeros(MAX_SIZE)
Z = np.zeros((MAX_SIZE, MAX_SIZE))
IERR = 0

# Initialize indices and basis functions
IC = []
LB = []
MB = []
NB = []

IG = 0
for I in range(1, 4):
    for L in range(1, NN+2):
        for M in range(1, NN+2):
            for N in range(1, NN+2):
                if L + M + N > NN + 3:
                    break 
                else:
                    IG += 1
                    IC.append(I)
                    LB.append(L - 1)
                    MB.append(M - 1)
                    NB.append(N - 1)
            
NR = IG
IC = np.array(IC)
LB = np.array(LB)
MB = np.array(MB)
NB = np.array(NB)

# Compute the elements of the E and GAMMA matrices
for IG in range(NR):
    for JG in range(IG, NR):
        I = IC[IG]
        J = IC[JG]
        LS = LB[IG] + LB[JG]
        MS = MB[IG] + MB[JG]
        NS = NB[IG] + NB[JG]

        GAMMA[IG, JG] = (
            C[I-1, 0, J-1, 0] * LB[IG] * LB[JG] * F(LS-2, MS, NS) +
            C[I-1, 1, J-1, 1] * MB[IG] * MB[JG] * F(LS, MS-2, NS) +
            C[I-1, 2, J-1, 2] * NB[IG] * NB[JG] * F(LS, MS, NS-2) +
            (C[I-1, 0, J-1, 1] * LB[IG] * MB[JG] + C[I-1, 1, J-1, 0] * MB[IG] * LB[JG]) * F(LS-1, MS-1, NS) +
            (C[I-1, 0, J-1, 2] * LB[IG] * NB[JG] + C[I-1, 2, J-1, 0] * NB[IG] * LB[JG]) * F(LS-1, MS, NS-1) +
            (C[I-1, 1, J-1, 2] * MB[IG] * NB[JG] + C[I-1, 2, J-1, 1] * NB[IG] * MB[JG]) * F(LS, MS-1, NS-1)
        )
        GAMMA[JG, IG] = GAMMA[IG, JG]
        if I == J:
            E[IG, JG] = F(LS, MS, NS)
        E[JG, IG] = E[IG, JG]

# Solve the eigenvalue problem
W, Z = eigh(GAMMA, E, driver='gvd')

# Compute frequencies from eigenvalues
W = np.sqrt(np.maximum(0.0, W) / RHO) / TWOPI

# Print the lowest 36 frequencies
print(f"FREQUENCIES FOR CORNER PRISM, NN= {NN}")
for i in range(min(36, len(W))):
    print(f"{W[i]:12.5g}", end=" ")
print()

# if __name__ == "__main__":
#     main()
