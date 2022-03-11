
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import scipy.io as sio
import networkx as nx
import matplotlib
import matplotlib as mpl
from matplotlib.widgets import Button,TextBox
label_size = 8
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
#---------------------------------------------#
def Hamming_vs_nB_2(D,NCells):
    B_max = NCells
    M_max = NCells

    if NCells == 30:
        I  = np.array([2]*10 + [1]*10 + [0]*10) #ideal flag
    elif NCells == 60:
        I  = np.array([2]*20 + [1]*20 + [0]*20) #ideal flag
    elif NCells == 90:
        I  = np.array([2]*30 + [1]*29 + [0]*31) #ideal flag

    M = NCells - np.sum(D == I,1)
    B = np.count_nonzero(np.diff(D,1),1)

    Hist = np.zeros((B_max+1,M_max+1))
    for i in range(D.shape[0]):
        Hist[B[i],M[i]] += 1
    return(Hist)
#-------------------------------------------------------#
def Fraction_Normal(D):
    B = (np.count_nonzero(np.diff(D,1),1) == 2)
    First = (D[:,0] == 2)
    Last  = (D[:,-1] == 0)
    Blue = np.sum(D == 0,1)
    White = np.sum(D == 1,1)
    Red = np.sum(D == 2,1)
    Z = np.count_nonzero(B*First*Last*Blue*White*Red)/D.shape[0]
    return(Z)
#------------------------------------------------------#
def Summarize(D):
    Frac = np.zeros((30,3))
    for i in range(30):
        Frac[i,0] = np.sum(D[:,i] == 0)
        Frac[i,1] = np.sum(D[:,i] == 1)
        Frac[i,2] = np.sum(D[:,i] == 2)
    Frac = Frac/1e5
    #-----------------------------------#
    Var = np.var(D,0)
    #-----------------------------------#
    H = Hamming_vs_nB_2(D,NCells=30)
    tmp = np.sum(H,1)
    nB = np.array([tmp[0],tmp[1],tmp[2],tmp[3],np.sum(tmp[4:])]) #4 and above as one class
    #-----------------------------------#
    fN = Fraction_Normal(D)
    nB[nB == 0] = 1
    nB = np.log10(nB)
    return(fN,Var,Frac,nB)
#--------------------------------------------------------#
def int_to_trinary(a):
    s = ''
    quotient = a ;
    while quotient > 2:
        rem = quotient%3
        quotient = quotient//3
        s = str(rem) + s
    s = str(quotient) + s
    s = '0'*(6-len(s)) + s
    return(s)
#--------------------------------------------------------#
def trinary_to_int(a):
    val = 0
    for i in range(len(a)):
        val += (int(a[-(i+1)])*(3**i))
    return(val)
#---------------------------------------------------------#

fig3 = plt.figure(constrained_layout=True,figsize=(10,15))
#plt.style.use('ggplot')

gs = fig3.add_gridspec(6,10)
ax1 = fig3.add_subplot(gs[1:4,0:3])
ax1b = fig3.add_subplot(gs[1:4,3:6])
ax2 = fig3.add_subplot(gs[0:4,6:],projection='3d')
ax3 = fig3.add_subplot(gs[4,:])
ax4 = fig3.add_subplot(gs[5,:])
#--------------------------------------------------------#
VarData = np.load('Vis_Data/Ultimate_VarData.npy')
FateData = np.load('Vis_Data/Ultimate_FateData.npy')
FF_Frac = np.load('Vis_Data/Ultimate_FF_Frac.npy')
HistData = np.load('Vis_Data/Ultimate_HistData.npy')
nBData = np.sum(HistData,2)
nBData[:,4] = np.sum(nBData[:,4:],1)
nBData = nBData[:,0:5]
nBData[nBData == 0] = 1
nBData = np.log10(nBData)

HistData[HistData == 0] = 1
HistData = np.log10(HistData)

_x = np.arange(0,31)
_y = np.arange(0,15)
_xx, _yy = np.meshgrid(_x, _y)
x, y = _xx.ravel(), _yy.ravel()
width = depth = 0.8
#--------------------------------------------------------#
def submit(expression):
    ctstring = expression.replace('-1','2')
    print(ctstring)
    CID = trinary_to_int(ctstring)
    #---------------------------------#
    ax1.clear()
    ax1b.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    #=====================================#
    ax1.text(0.0,0.0,'$N^b$',fontsize=10)
    ax1.text(1,1,'$B$',fontsize=10)
    ax1.text(1,0,'$W$',fontsize=10)
    ax1.text(1,-1,'$R$',fontsize=10)
    ax1.text(2,0,'$L$',fontsize=10)

    if ctstring[0] == '1':
        ax1.arrow(0.1,0,dx=0.7,dy=0.8,head_width=0.1,color='g')
    elif ctstring[0] == '2':
        ax1.arrow(0.1,0,dx=0.7,dy=0.8,head_width=0.1,color='r')

    if ctstring[1] == '1':
        ax1.arrow(0.1,0,dx=0.7,dy=0.0,head_width=0.1,color='g')
    elif ctstring[1] == '2':
        ax1.arrow(0.1,0,dx=0.7,dy=0.0,head_width=0.1,color='r')

    if ctstring[2] == '1':
        ax1.arrow(0.1,0,dx=0.7,dy=-0.8,head_width=0.1,color='g')
    elif ctstring[2] == '2':
        ax1.arrow(0.1,0,dx=0.7,dy=-0.8,head_width=0.1,color='r')

    if ctstring[3] == '1':
        ax1.arrow(1.1,1,dx=0.7,dy=-0.8,head_width=0.1,color='g')
    elif ctstring[3] == '2':
        ax1.arrow(1.1,1,dx=0.7,dy=-0.8,head_width=0.1,color='r')

    if ctstring[4] == '1':
        ax1.arrow(1.1,0,dx=0.7,dy=0.0,head_width=0.1,color='g')
    elif ctstring[4] == '2':
        ax1.arrow(1.1,0,dx=0.7,dy=0.0,head_width=0.1,color='r')

    if ctstring[5] == '1':
        ax1.arrow(1.1,-1,dx=0.7,dy=0.8,head_width=0.1,color='g')
    elif ctstring[5] == '2':
        ax1.arrow(1.1,-1,dx=0.7,dy=0.8,head_width=0.1,color='r')


    ax1.set_xlim(-0.2,2.2)
    ax1.set_ylim(-1.5,1.5)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_facecolor('#ffe6cc')
    #==============================================================#
    ax1b.bar([0,1,2,3,4],nBData[CID,:],width=0.8,edgecolor='k')
    ax1b.set_yticks([0,1,2,3,4,5])
    ax1b.set_yticklabels(['$10^0$','$10^1$','$10^2$','$10^3$','$10^4$','$10^5$'])
    ax1b.set_xticks([0,1,2,3,4])
    ax1b.set_xticklabels(['0','1','2','3','>3'])
    ax1b.set_facecolor('#ffe6cc')
    ax1b.set_xlabel('$n_B$',fontsize=14)
    #ax1b.yaxis.tick_right()

    #===============================================================#
    H = HistData[CID,:,:]
    bottom = np.zeros(H.shape)
    ax2.bar3d(x, y, bottom.ravel(), width, depth, H.ravel(),color='#66ff33',shade=False,edgecolor='k',linewidth=0.5)
    ax2.grid(False)
    ax2.invert_xaxis()
    ax2.view_init(25,152)
    ax2.set_xticks([0,10,20,30])
    ax2.set_yticks([0,5,10,15])
    ax2.set_zticks([1,3,5])
    ax2.set_zticklabels(['$10^1$','$10^3$','$10^5$'])
    ax2.set_ylabel('$n_B$')
    ax2.set_xlabel('$d_H$')


    #===============================================================#
    Frac = FateData[CID,:,:]
    ax3.bar(np.arange(1,31),Frac[:,0],color='b',width=0.8,edgecolor='k')
    ax3.bar(np.arange(1,31),Frac[:,1],bottom=Frac[:,0],color='w',width=0.8,edgecolor='k')
    ax3.bar(np.arange(1,31),Frac[:,2],bottom=Frac[:,0]+Frac[:,1],color='r',width=0.8,edgecolor='k')
    ax3.set_ylim(0,1)
    ax3.set_yticks([0,0.5,1])
    ax3.set_xticks([])
    ax3.set_xlim(0.5,30.5)
    ax3.set_facecolor('#ffe6cc')
    ax3.set_ylabel('p',fontsize=14)
    #================================================================#
    Var = VarData[CID,:]
    ax4.bar(np.arange(1,31),Var,width=0.8,color='silver',edgecolor='k')
    ax4.set_xticks([1,5,10,15,20,25,30])
    ax4.set_ylabel('$\sigma^2$',fontsize=14)
    ax4.set_xlim(0.5,30.5)
    ax4.set_facecolor('#ffe6cc')
    ax4.set_xlabel('pos',fontsize=14)
#==================================================================#
axbox = fig3.add_axes([0.2, 0.86, 0.2, 0.045])
text_box = TextBox(axbox, "Coupling Type : ")
text_box.on_submit(submit)
text_box.set_val("111-1-1-1")  # Trigger `submit` with the initial string.

plt.show()
