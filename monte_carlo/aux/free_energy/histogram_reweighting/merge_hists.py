import numpy as np

ne = 500
nv = 3000

H = np.zeros( nv*ne )

for i in range(105, 113):
    
    if(i==111):
        continue

    T_string = "T0p"+str(i)+"0"
    for j in range(50,101,5):
        if(j == 55 or j==95):
            continue
        elif(j < 100):
            P_string = "P0p00" + str(j)
        else:
            P_string = "P0p0" + str(j)
        tot_string = T_string + "/" + P_string + "/hist.out"

        f = open(tot_string, mode='r')

        for k, line in enumerate(f):
            l = float(line.split()[2])
            H[k] = H[k] + l

        f.close()

e_min = -8000
e_max = -7500

v_min =  2000
v_max =  5000

f_hist = open("merged_histogram.dat",mode='w')

for i, v in enumerate(range(v_min, v_max,1)):
    for j, e in enumerate(range(e_min, e_max,1)):
        if(H[ i*ne + j ] > 0):
            f_hist.write( str(v) + " " + str(e) + " " + str(H[ i*ne + j ]) + '\n' )     
