import scipy.io as sio
import numpy as np
import falco_systemID as falco
import falco_systemID2 as falco2
import sys
if __name__ == "__main__":

	path2data = sys.argv[1]
	lr = np.float(sys.argv[2])
	lr2 = np.float(sys.argv[3])
	epoch = np.float(sys.argv[4])
	print_flag = np.bool(sys.argv[5])
	Q0 = np.float(sys.argv[6])
	Q1 = np.float(sys.argv[7])
	R0 = np.float(sys.argv[8])
	R1 = np.float(sys.argv[9])
	if len(sys.argv) > 10:
		R2 = np.float(sys.argv[10])
		IDnet = falco2.vl_net(Q0=Q0, Q1=Q1, R0=R0, R1=R1, R2=R2, path2data=path2data)
		falco2.linear_vl(IDnet, lr=lr, lr2=lr2, epoch=epoch, print_flag=print_flag)
	else:
		IDnet = falco.vl_net(Q0=Q0, Q1=Q1, R0=R0, R1=R1, path2data=path2data)
		falco.linear_vl(IDnet, lr=lr, lr2=lr2, epoch=epoch, print_flag=print_flag)
	# for k in range(len(sys.argv)):
	# 	print(sys.argv[k])