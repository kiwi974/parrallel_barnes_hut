""" Script that plot the graph corresponding to the weak scaling study """



import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

colors = [(223./255, 109./255, 20./255),(0,0,0)]

# Table with the number of cores
cores = [[1,2,4,8,16,32,64,128],[4,8,16,32,64,128]]

" ************************************************************ "
" ****************** STRONG SCALING DATA ********************* "
" ************************************************************ "

# Table with the average time per step for 100 000 particule in total
strongTimes = [[100000,[18.713017,2.905857,1.467321,0.800800,0.617382,0.469426,0.327056,0.308315]],
               [500000,[8.388025,4.276133,2.249238,1.292339,1.907628,1.495633]]]



" ************************************************************ "
" ****************** STRONG SCALING DATA ********************* "
" ************************************************************ "

# Nomber of particules per core
nbPartCore = [5000,10000]

# Table with times for weak scaling : time according to the work load per processor
weakTimes = [[nbPartCore[0],[0.537553,0.199486,0.210799,0.300724,0.581771,0.836134,1.290899,1.084245]],
             [nbPartCore[1],[2.14,0.849542,0.974254,1.144641,1.397159,2.279065,6.471566,9.7]]]


" ***** Computation of scaling features ***** "

# Compute speedup per processor with
#       Sp = t1/tp
speedup = [[nbPartCore[j],[weakTimes[j][1][0]/tp for tp in weakTimes[j][1]]] for j in range(len(nbPartCore))]


# Compute the efficiency per processor with
#       Ep = Sp/p
efficiency = [[nbPartCore[j],[speedup[j][1][i]/cores[0][i] for i in range(len(cores[0]))]] for j in range(len(speedup))]





" ************************************************************ "
" ****************** STRONG SCALING CURVES ******************* "
" ************************************************************ "

fig = plt.figure()
ax = plt.gca()


for k in range(len(strongTimes)):

    " ***** Plot scaling curvres ***** "
    ax.scatter(cores[k], strongTimes[k][1], s=3, c=colors[k], marker='*',
               label=str(strongTimes[k][0])+' particules')


    " ***** Plot linear regression curvres ***** "

    fit = np.polyfit(cores[k], strongTimes[k][1], 1)
    fit_fn = np.poly1d(fit)

    ax.plot(cores[k],fit_fn(cores[k]),c=colors[k],linewidth=0.5,label=str(strongTimes[k][0]) + ' particules ')

    #slope, intercept, r_value, p_value, std_err = stats.linregress(cores[k],strongTimes[k][1])
    #line = slope*cores[k]+intercept
    #ax.plot(cores[k],line)

" ***** Show ***** "

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel("Number of cores (CPU)")
ax.set_ylabel("Time per step (s)")
ax.set_title("Strong scaling N-Body problem")
ax.legend()




" ************************************************************ "
" ******************* WEAK SCALING CURVES ******************** "
" ************************************************************ "

fig2, (ax1, ax2) = plt.subplots(2,sharey = True)


" ***** Plot scaling curvres ***** "

ax1.scatter(cores[0],weakTimes[0][1], s=3, c='k', marker='*')
ax1.scatter(cores[0],weakTimes[1][1], s=3, c='m', marker='*')

ax2.plot(cores[0],efficiency[0][1],'-',linewidth=0.8,label=str(weakTimes[0][0]) + " part. per proc")
ax2.plot(cores[0],efficiency[1][1],'-',linewidth=0.8,label=str(weakTimes[1][0]) + " part. per proc")

" ***** Plot linear regression curvres ***** "

fit1 = np.polyfit(cores[0], weakTimes[0][1], 1)
fit_fn1 = np.poly1d(fit1)
ax1.plot(cores[0],fit_fn1(cores[0]),'-',linewidth=0.8,label=str(weakTimes[0][0]) + " part. per proc")


fit2 = np.polyfit(cores[0], weakTimes[1][1], 1)
fit_fn2 = np.poly1d(fit2)
ax1.plot(cores[0],fit_fn2(cores[0]),'-',linewidth=0.8,label=str(weakTimes[1][0]) + " part. per proc")

" ***** Show ***** "

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel("Number of cores (CPU)")
ax1.set_ylabel("Runtime")
ax1.legend()

ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlabel("Number of cores (CPU)")
ax2.set_ylabel("Efficiency")
ax2.legend()

plt.legend()
plt.show()