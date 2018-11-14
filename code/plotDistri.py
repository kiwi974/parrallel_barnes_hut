###### View the evolution of the distributions over the time ######

import matplotlib.pyplot as plt
import time
from operator import itemgetter


# Path until the file containing the distributions
filePath = "distriTest"


distriFile = open(filePath)

contenu = distriFile.read()
contenu = contenu[:-4]

allDistri = contenu.split("<-->")
timeStamp = len(allDistri)

plt.show()

linex = 0

for k in range(timeStamp):
    allDistri[k] = allDistri[k].split("<>")
    absc = []
    ordo = []
    for i in range(len(allDistri[k])):
        position = allDistri[k][i].split(",")
        try:
            allDistri[k][i] = (float(position[0]), float(position[1]))
            #break
        except ValueError:
            allDistri[k][i] = (float(position[0][len(position[0])-8:]), float(position[1]))
    allDistri[k] = sorted(allDistri[k], key =itemgetter(0))
    for i in range(len(allDistri[k])):
        absc.append(allDistri[k][i][0])
        ordo.append(allDistri[k][i][1])

    n = int(len(absc)/2)
    abs1 = absc[:n]
    abs2 = absc[n+1:]
    ord1 = ordo[:n]
    ord2 = ordo[n+1:]

    if (k==0):
        linex = (abs1[len(abs1)-1]+abs2[len(abs2)-1])/2

    plt.xlim(0,1000);
    plt.ylim(0,1000);
    plt.scatter(abs1,ord1,s=2, c='k', marker='*')
    plt.scatter(abs2,ord2,s=2, c='k', marker='*')
    #plt.plot([linex,linex],[0,1000],'k')
    plt.text(900,900,k)
    plt.draw()
    plt.pause(0.001)
    time.sleep(.2)
    plt.clf()




