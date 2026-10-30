import numpy as np
import os as _os
D=_os.path.join(_os.path.dirname(_os.path.abspath(__file__)),"original_specsim2")+_os.sep   # folder with the outputs of SPECSIM2.EXE
def vs(f,step,kmax):
    sx,sy=step; out=[]
    for k in range(1,kmax+1):
        dx,dy=k*sx,k*sy
        a=f[max(0,-dx):f.shape[0]-max(0,dx), max(0,-dy):f.shape[1]-max(0,dy)]
        b=f[max(0,dx):f.shape[0]-max(0,-dx), max(0,dy):f.shape[1]-max(0,-dy)]
        out.append(0.5*np.mean((a-b)**2))
    return np.array(out)
lags=np.array([1,2,5,10,20,40])
model=lambda h: 0.3+0.7*(1-np.exp(-h/10.0))
print("model gamma(h) h=1,2,5,10,20,40:",np.round(model(lags),3))
for name in ["310107a7-specsim2.out","042831ad-specsim_2.ou","056a6798-specsim_3.ou","ac514076-TEST1.DAT"]:
    d=np.loadtxt(D+name); z=d[:,2].reshape(128,128)   # x slow, y fast
    gx=vs(z,(1,0),40); gy=vs(z,(0,1),40)
    print(name,"mean %.3f var %.3f | gx"%(z.mean(),z.var()),np.round(gx[lags-1],3),"gy",np.round(gy[lags-1],3))
for tag in ["fp2d_pad0","fp2d_pad1"]:
    raw=np.loadtxt(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)),f"{tag}.out"),skiprows=22)  # title,nreal,20 names
    allx=[];ally=[];vv=[]
    for r in range(raw.shape[1]):
        z=raw[:,r].reshape(128,128)   # x fastest -> reshape (y,x)
        allx.append(vs(z,(1,0),40)); ally.append(vs(z,(0,1),40)); vv.append(z.var())
    print(tag,"mean var %.3f | gx"%np.mean(vv),np.round(np.mean(allx,0)[lags-1],3),"gy",np.round(np.mean(ally,0)[lags-1],3))
