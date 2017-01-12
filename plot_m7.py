import numpy as np 
import matplotlib.pyplot as plt 
import netCDF4 as nc 


def read_SD(infile):
	NUM=[]
	RW=[]
	indataset=nc.Dataset(infile)
	for var in indataset.variables:
		if 'RWET' in var:
			print var
			RW.append(indataset.variables[var][:])
		if 'N_' in var:
			print var
			NUM.append(indataset.variables[var][:])
	RWa=np.asarray(RW)
	NUMa=np.asarray(NUM)
	return RWa,NUMa
def plot_m7(R,N):
	lnr=np.linspace(-50,10,6000)
	pns=np.zeros(max(size(num)),max(size(lnr)))
	pks=np.zeros(max(size(num)),max(size(lnr)))
	pas=np.zeros(max(size(num)),max(size(lnr)))
	pcs=np.zeros(max(size(num)),max(size(lnr)))
	pki=np.zeros(max(size(num)),max(size(lnr)))
	pai=np.zeros(max(size(num)),max(size(lnr)))
	pci=np.zeros(max(size(num)),max(size(lnr)))
	for i in range(len(num[1,:])):
		pns[i,:]=(num[1,i]/(np.sqrt(2*np.pi)*np.log(1.59)))*np.exp(-((lnr-np.log(rwet[1,i]))**2)/(2*(np.log(1.59))**2))
		pks[i,:]=(num[2,i]/(np.sqrt(2*np.pi)*np.log(1.59)))*np.exp(-((lnr-np.log(rwet[2,i]))**2)/(2*(np.log(1.59))**2))
		pas[i,:]=(num[3,i]/(np.sqrt(2*np.pi)*np.log(1.59)))*np.exp(-((lnr-np.log(rwet[3,i]))**2)/(2*(np.log(1.59))**2))
		pcs[i,:]=(num[4,i]/(np.sqrt(2*np.pi)*np.log(2.00)))*np.exp(-((lnr-np.log(rwet[4,i]))**2)/(2*(np.log(2.00))**2))
		pki[i,:]=(num[5,i]/(np.sqrt(2*np.pi)*np.log(1.59)))*np.exp(-((lnr-np.log(rwet[5,i]))**2)/(2*(np.log(1.59))**2))
		pai[i,:]=(num[6,i]/(np.sqrt(2*np.pi)*np.log(1.59)))*np.exp(-((lnr-np.log(rwet[6,i]))**2)/(2*(np.log(1.59))**2))
		pci[i,:]=(num[7,i]/(np.sqrt(2*np.pi)*np.log(2.00)))*np.exp(-((lnr-np.log(rwet[7,i]))**2)/(2*(np.log(2.00))**2))
	data_s=pns+pks+pas+pcs;
	data_i=pki+pai+pci;
	plt.loglog(2*np.exp(lnr),data_s[1,:]*1e-6,'b',lw=3)
	plt.show()
R,N=read_SD('../output/20_28.hourly.nc')
print np.shape(R),np.shape(N)
plt.figure()
R_hyde=np.squeeze(R[:,0,:,14,15])
N_hyde=np.squeeze(N[:,0,:,14,15])
plt.plot(R,N,'r',lw=2)
plt.show()
