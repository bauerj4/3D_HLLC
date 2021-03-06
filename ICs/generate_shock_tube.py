import sys
import numpy as np

N_CELLS_PER_SIDE = (2**int(sys.argv[1]))

f = open("parent_mesh.dat",'w')

dt = np.dtype('f64')

rho_l = 1.
ux_l = 0.75
uy_l = 0.
uz_l = 0.
p_l = 1.0

rho_r = 0.125
ux_r = 0.0
uy_r = 0.
uz_r = 0.
p_r = 0.1

BARRIER = N_CELLS_PER_SIDE/4
MESH = np.ndarray(shape=(N_CELLS_PER_SIDE, N_CELLS_PER_SIDE, N_CELLS_PER_SIDE, 5), dtype=float)

for i in range(0, N_CELLS_PER_SIDE):
    for j in range(0, N_CELLS_PER_SIDE):
        for k in range(0, N_CELLS_PER_SIDE):
            if (i > BARRIER):
                MESH[i][j][k][0] = rho_r
                MESH[i][j][k][1] = ux_r
                MESH[i][j][k][2] = uy_r
                MESH[i][j][k][3] = uz_r
                MESH[i][j][k][4] = p_r

                #f.write(rho_r)
                #f.write(ux_r)
                #f.write(uy_r)
                #f.write(uz_r)
                #f.write(p_r)
                #print "RIGHT_SIDE"
            else:
                #print "LEFT SIDE"
                MESH[i][j][k][0] = rho_l
                MESH[i][j][k][1] = ux_l
                MESH[i][j][k][2] = uy_l
                MESH[i][j][k][3] = uz_l
                MESH[i][j][k][4] = p_l

                #f.write(rho_l.tobytes())
                #f.write(ux_l)
                #f.write(uy_l)
                #f.write(uz_l)
                #f.write(p_l)
MESH.tofile(f)
