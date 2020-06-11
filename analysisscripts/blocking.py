# Common imports
import os
# Where to save the figures and data files

from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, sqrt, concatenate
from numpy.linalg import inv

import pandas as pd
from pandas import DataFrame

def block(x):
    # preliminaries
    n = len(x)
    d = int(log2(n))
    s, gamma = zeros(d), zeros(d)
    mu = mean(x)

    # estimate the auto-covariance and variances 
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])
   
    # generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")

    return mu, s[k]/2**(d-k)


def calculate_mean_squared(x):
    mean_squared = 0.0
    for i in range(len(x)):
        mean_squared += x[i]*x[i]
    mean_squared /= float(len(x))
    return mean_squared

def calculate_uncor_var(x):
    uncor_var = 0.0
    x0 = mean(x)
    for i in range(len(x)):
        uncor_var += (x[i] - x0)**2
    uncor_var /= float(len(x))
    return uncor_var

def exact_energy(alpha, N):
    return (0.5*alpha + 1/(8.0*alpha))*N*3

# exercise c - brute force - no interaction

DATA_ID = "../Output//exercise_d//allEnergies"

Ns = ["2"]#, "50", "100", "500"]#, "100"]#, "500"]

As = ["50", "60", "70", "80", "90", "100", "110", "120", "130", "140"]

methodtype = "analytical"
samplingtype = "brute_force"


# for j in range(len(Ns)):
#     print "N: ", Ns[j]
#     print "$\\alpha$: & $\\left< E_L \\right>$: & SEM: & $\\sigma_B$: & CPU time:"+ "\\" + "\\"
#     cpu_time = loadtxt("../Output//exercise_d//"+samplingtype+".txt")
#     for a in range(len(As)):    
        
#         def data_path(dat_id):
#             return os.path.join(DATA_ID, dat_id)

#         infile = open(data_path("analytical_2d_%sp_alpha_%s_MC_21_%s_energy.txt"%(Ns[j],As[a], samplingtype)),'r')

#         x = loadtxt(infile, skiprows=5)
#         # x = x[:int(2**19)]

#         (mu, variance) = block(x) 
#         std = sqrt(variance)

#         uncor_std = sqrt(calculate_uncor_var(x))/sqrt(len(x)) # sqrt(calculate_mean_squared(x) - mu*mu)
        
#         alpha = float(As[a])/100
#         N = float(Ns[j])

#         print "%.2f & %.5f & %.5f & %.5f & %.5f"%(alpha, mu, uncor_std, std, cpu_time[a] )+ "\\"+ "\\"

#     print "Mean CPU time & & & ", mean(cpu_time)


# -------------------------------------------
# Ground state energies with interaction
# -------------------------------------------

DATA_ID = "../Output//"

def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

# omega = [1.0, 0.5, 0.1, 0.05, 0.01]
# sampling_type = "importance" 


# Brute force values
# alpha = [0.98847, 0.98061, 0.94693, 0.92747, 0.88398]
# beta = [0.39965, 0.31091, 0.17764, 0.13815, 0.07287]
# mean_distance = [1.63619, 2.48111, 6.6947, 10.3885, 29.1774]
# E_kin = [0.894379, 0.448803, 0.100298, 0.0532853, 0.0129223]
# E_pot = [1.29896, 0.705148,  0.176739, 0.0997244, 0.0283832]
# E_int = [0.813457, 0.513494, 0.171606, 0.107577, 0.0363629]

# Importance values

# N = 2

# alpha = [0.98846, 0.98082, 0.94734, 0.92262, 0.88305]
# beta = [0.39954, 0.31068, 0.17810, 0.14090, 0.07366]
# mean_distance = [1.64345, 2.48083, 6.72361, 10.3333, 29.2928 ]
# E_kin = [0.893123, 0.454685, 0.0988959, 0.0495128, 0.0130627]
# E_pot = [1.30517, 0.699729, 0.178708,  0.102437, 0.0283079]
# E_int = [0.808643,  0.51295, 0.170973, 0.109083, 0.0362095]

# MC = 24:
# alpha = [0.98846, 0.98082, 0.94734, 0.92262, 0.88305]
# beta = [0.39954, 0.31068, 0.17810, 0.14090, 0.07366]
# mean_distance = [1.63094, 2.4714,  6.74423, 10.3693,28.9974 ]
# E_kin = [0.900155, 0.45106, 0.098532, 0.0535869, 0.0130799]
# E_pot = [1.29065, 0.700326, 0.178782,  0.100298, 0.0276834]
# E_int = [0.815996,  0.515984, 0.170941, 0.107271,  0.0367146]


# N = 6

# alpha = [0.71567, 0.75823, 0.78852, 0.76518, 0.64907]
# beta = [0.49372, 0.34260, 0.15041, 0.10733, 0.05085]
# mean_distance = [0, 0, 0, 0, 0]
# E_kin = [2.34291, 1.32263,  0.29513, 0.117831, 0.00208272]
# E_pot = [10.7076, 5.80937, 1.7035, 1.0882, 0.380349]
# E_int = [7.39876, 4.85476, 1.65562, 1.01623, 0.336716]

# alpha = [0.71567, 0.75823, 0.78852, 0.76518, 0.64907]
# beta = [0.49372, 0.34260, 0.15041, 0.10733, 0.05085]
# mean_distance = [2.76801, 3.99339,  9.91115, 15.2152,  38.5545]
# E_kin = [ 2.40112,  1.33601, 0.283704, 0.132386, -0.0211001]
# E_pot = [ 10.6658, 5.79185, 1.70936, 1.06189, 0.399848]
# E_int = [ 7.39113, 4.85717,  1.65167,  1.02717,0.329007]

# for o in range(len(omega)):

#     infile = open(data_path("exercise_g//one_body//interaction_ground_state_"+ sampling_type + "_%ip_omega_%i_MC_24_energy.txt"%(N,int(omega[o]*100))),'r')

#     x = loadtxt(infile, skiprows=5)

#     (mu, variance) = block(x) 
#     std = sqrt(variance)

#     uncor_std = sqrt(calculate_uncor_var(x))/sqrt(len(x)) # sqrt(calculate_mean_squared(x) - mu*mu)

#     print "%.2f & %.5f & %.5f & %.4f & %.5f & %.5f & %.3f & %.4f & %.4f & %.4f"%(omega[o], alpha[o], beta[o], mu, uncor_std, std, mean_distance[o], E_kin[o], E_pot[o], E_int[o])+ "\\"+ "\\"


# ---------------------------------------------
# Importance sampling 
# ---------------------------------------------

methodtype="analytical"

samplingtype = "brute_force"

DATA_ID = "../Output//exercise_d//allEnergies"

print "$left< E_L right>$: & SEM: & $sigma_B$: & Acceptance & CPU time:"+ "\\" + "\\"

M = loadtxt("../Output//exercise_d//analytical_2p_2d_"+ samplingtype + ".txt", skiprows=1)

stepsize = M[:,0]
acc = M[:,1]
cpu_time = M[:,-1]

for a in range(len(stepsize)):    
    
    def data_path(dat_id):
        return os.path.join(DATA_ID, dat_id)

    printablestep = int(stepsize[a]*1000)
    printedstep = str(printablestep)

    infile = open(data_path(methodtype+"_2d_2p_stepsize_" + printedstep + "_MC_21_"+samplingtype + "_energy.txt"),'r')

    x = loadtxt(infile, skiprows=5)
    # x = x[:int(2**19)]

    (mu, variance) = block(x) 
    std = sqrt(variance)

    uncor_std = sqrt(calculate_uncor_var(x))/sqrt(len(x)) # sqrt(calculate_mean_squared(x) - mu*mu)
    

    print "%.3f & %.5f & %.5f & %.3f& %.3f"%(mu, uncor_std, std, acc[a], cpu_time[a] )+ "\\"+ "\\"

print "Mean CPU time & & & ", mean(cpu_time)

# ---------------------------------------------
# Parallelized code
# ---------------------------------------------

# print "Parallel:"

# identifyers = [1,2,3,4]

# x = [{} for m in range(len(identifyers))]

# for n in identifyers:
#     infile = open(data_path("exercise_i/parallell_checking_6p_importance_interaction_MC_22_%i_energy.txt"%(n)),'r')

#     x[n-1] = loadtxt(infile, skiprows=5)

# z = concatenate((x[0], x[1], x[2], x[3]))
# (mu, variance) = block(z) 
# std = sqrt(variance)

# uncor_std = sqrt(calculate_uncor_var(z))/sqrt(len(z)) # sqrt(calculate_mean_squared(x) - mu*mu)

# print "%.4f & %.5f & %.5f "%(mu, uncor_std, std)+ "\\"+ "\\"

# # Not in parallel

# print "Not parallel:"

# infile = open(data_path("exercise_i/parallell_checking_6p_importance_interaction_MC_24_0_energy.txt"),'r')

# y = loadtxt(infile, skiprows=5)

# (mu_0, variance_0) = block(y) 
# std_0 = sqrt(variance_0)

# uncor_std_0 = sqrt(calculate_uncor_var(y))/sqrt(len(y)) # sqrt(calculate_mean_squared(x) - mu*mu)

# print "%.4f & %.5f & %.5f "%(mu_0, uncor_std_0, std_0)+ "\\"+ "\\"

