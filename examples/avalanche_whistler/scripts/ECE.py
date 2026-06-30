'''
ECE (Electron Cyclotron Emission) 计算程序

该程序读取 'data.npy' 文件中的多个时间步的等离子体分布函数数据，
并计算每个时间步对应的 ECE 功率。计算结果存储在 Power 数组中，
其中 Power[时间步, 频率索引] 表示对应时间步和频率下的 ECE 功率。

主要功能：
1. 初始化网格和计算相关系数矩阵（如果未缓存）
2. 从 'data.npy' 文件中读取所有时间步的数据
3. 对每个时间步的数据计算 ECE 功率
4. 输出所有时间步的 ECE 功率结果

data_dir 传参方法：
- 默认情况下，程序使用 config/parameters.py 中定义的默认数据目录
- 可通过命令行参数 --data_dir 指定自定义数据目录，例如：
  python src/core/ECE.py --data_dir /path/to/data/

注：程序计算的是多个时间步的 ECE，而非单个特定时刻。
'''


import numpy as np
from numpy import pi,power,exp,sqrt,abs,log,sin,cos,arccos,floor,ceil
import scipy as sp
import scipy.sparse
import scipy.special
import matplotlib.pyplot as plt
import os.path
import ipdb 
import sys
import argparse

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from src.config.parameters import *
from src.utils import mesh
from src.modules import kinetic, source, whistler, exel


def calculate_ECE(data_dir_param=None):
    global data_dir, density, Te, Zeff, vth, R, nump, numxi
    if data_dir_param is not None:
        # 使用parameters.py中的set_data_dir函数来设置新的数据目录
        from src.config.parameters import set_data_dir
        set_data_dir(data_dir_param)
        # 刷新本模块中 import * 复制的变量
        import src.config.parameters as _p
        data_dir = _p.data_dir
        density = _p.density
        Te = _p.Te
        Zeff = _p.Zeff
        vth = _p.vth
        R = _p.R
        nump = _p.nump
        numxi = _p.numxi

    a=0.67  #小半径
    numomega=2
    omegace=whistler.omegace
    omegace_max=whistler.omegace*R/(R-a)
    omegace_min=whistler.omegace*R/(R+a)
    omega_array=np.linspace(2*omegace,3*omegace,numomega)
    nu_e=2.9e-6*(density/1e6)*13*(Te*1000)**(-3.0/2.0)*Zeff
    numr=80;
    nump2=nump*2-1

    mesh.init_mesh()

    if os.path.isfile(data_dir + 'ECE_KImatrix.npz'):
        r_array=np.load(data_dir + 'ECE_r_array.npy')
        ECKI_sparse=sp.sparse.load_npz(data_dir + 'ECE_KImatrix.npz')
        ECK_sparse=sp.sparse.load_npz(data_dir + 'ECE_Kmatrix.npz')
        KIM_vec=np.load(data_dir + 'ECE_KIM_vec.npy')
        KM_vec=np.load(data_dir + 'ECE_KM_vec.npy')
        ECKIO_sparse=sp.sparse.load_npz(data_dir + 'ECE_KIOmatrix.npz')
        ECKO_sparse=sp.sparse.load_npz(data_dir + 'ECE_KOmatrix.npz')
        KIMO_vec=np.load(data_dir + 'ECE_KIMO_vec.npy')
        KMO_vec=np.load(data_dir + 'ECE_KMO_vec.npy')

    else:
        r_array=np.zeros(numomega*numr)
        for i in range(numomega):
            omega=omega_array[i]
            i1=int(ceil(omega/omegace_max))
            i2=int(ceil(omega/omegace_min)+1)
            exlen=0
            rs=R-a
            step=[rs]
            for j in range(i1,i2):
                rp=omegace*R/(omega/j)
                rp1=rp*0.975
                if rp1>R+a:
                    exlen=exlen+(R+a-rs)
                    break
                if rp>R+a:
                    exlen=exlen+(rp1-rs)+20*(R+a-rp1)
                else:
                    if rp1<R-a:
                        exlen=exlen+20*(rp-(R-a))
                    else:
                        exlen=exlen+(rp1-rs)+20*(rp-rp1)
                rs=rp
                step.append(rp1)
                step.append(rp)
            step.append(R+a)
            interval=exlen/(numr-1)
            rr=[]
            rs=R-a
            for j in range(i1,i2):
                rp=omegace*R/(omega/j)
                rp1=rp*0.975
                if rp>=R+a:
                    if rp1<R+a:
                        rr.extend(np.arange(rs,rp1,interval))
                        rs=rp1
                    rr.extend(np.linspace(rs,R+a,numr-len(rr)))
                    break
                else:
                    if rp1<R-a:
                        rr.extend(np.arange(R-a,rp,interval/20))
                        rs=rp
                    else:
                        rr.extend(np.arange(rs,rp1,interval))
                        if ceil((rp-rp1)/(interval/20))+len(rr)>numr-1:
                            rr.extend(np.linspace(rp1,rp,numr-len(rr)-1))
                            rs=R+a
                        else:
                            rr.extend(np.arange(rp1,rp,interval/20))
                            rs=rp
            r_array[i*numr:i*numr+numr]=rr

        omegacearray=omegace*R/r_array
        rhoarray=abs(r_array-R)/a
        densityarray=density*(1-rhoarray**2)**2
        vtharray=vth*sqrt(1.8*(1-rhoarray**2)**2+0.2)
        omegapearray=sqrt(densityarray*1.6e-19**2*9e9*4*pi/9.1e-31)
        omegaarray=np.array([np.ones(numr)*omega_array[x] for x in range(numomega)]).flatten()

        ECKI_row=np.zeros(numomega*numr*(numxi*2)*6*50)
        ECKI_col=np.zeros(numomega*numr*(numxi*2)*6*50)
        ECKI_data=np.zeros(numomega*numr*(numxi*2)*6*50)
        ECK_data=np.zeros(numomega*numr*(numxi*2)*6*50)
        ECKIO_data=np.zeros(numomega*numr*(numxi*2)*6*50)
        ECKO_data=np.zeros(numomega*numr*(numxi*2)*6*50)
        collision_ki=np.zeros(numr*numomega)
        istep=0
        for i in range(numomega*numr):
            #for i in [3918]:
            print(i)
            omega=omegaarray[i]
            omegace0=omegacearray[i]
            omegape0=omegapearray[i]
            #k=omega/cc
            #k=sqrt((omega**7 - 2*omega**5*omegace0**2 + omega**3*omegace0**4 - 2*omega**5*omegape0**2 + 2*omega**3*omegace0**2*omegape0**2 + omega**3*omegape0**4 - omega*omegace0**2*omegape0**4)/
            #        (cc**2*omega**5 - 2*cc**2*omega**3*omegace0**2 + cc**2*omega*omegace0**4 - cc**2*omega**3*omegape0**2 + cc**2*omega*omegace0**2*omegape0**2))
            k=sqrt(((omega**2-omegape0**2)**2-omega**2*omegace0**2)/(omega**2-omegape0**2-omegace0**2)/cc**2)
            print(k)
            #ipdb.set_trace()
            Ex=omegace0*omegape0**2/(omega*(omega**2-omegace0**2-omegape0**2))
            Ey=1
            den=(-Ey**2)*2*k*cc**2  
            k2=sqrt((omega**2-omegape0**2)/cc**2)
            Ez=1
            den2=(-Ez**2)*2*k2*cc**2  
             
            nor=-((1+Ey**2)*(omegape0**2*(omega**2+omegace0**2)/(omega**2-omegace0**2)**2)-2*Ey*(-2*omegape0**2*omega*omegace0/(omega**2-omegace0**2)**2)+Ez**2*(omegape0**2/omega**2))*nu_e*omega
            collision_ki[i]=nor/den
            
            for n in range(2,51):
                gammax=n*omegace0/omega
                if gammax<1:
                    continue
                px=sqrt(gammax**2-1)
                if px==0:
                    continue
                pnode=np.searchsorted(mesh.nodes[0:nump2,0],px)
                if pnode==nump2:
                    continue
                ptris=pnode-2;
                if ptris<0:
                    ptris=0
                elif ptris%2==0:
                    ptris=ptris-1
                if ptris==0:
                    ptris2=2
                elif ptris<nump*2-5:
                    ptris2=ptris+3
                else:
                    ptris2=ptris+2
                kk_array=[]
                for ki in range(numxi-1):
                    kk_array.extend([ptris+ki*(nump-1)*2, ptris2+ki*(nump-1)*2])
                for kk in kk_array:
                    nd=np.zeros((4,2)) # Three triangle points of the element
                    nd[0:3]=np.array([mesh.nodes[ii] for ii in mesh.trisToNodes[kk][0:3]])
                    nd[3]=nd[0]
                    ximax=max(nd[:,1])
                    ximin=min(nd[:,1])
                    xicross=np.zeros(3) # Crossing point of the source term line and the element boundaries
                    for ii in range(3):
                        if abs(nd[ii,0]-nd[ii+1,0])<1e-20:
                            xicross[ii]=-2
                        elif abs(nd[ii,1]-nd[ii+1,1])<1e-20:
                            xicross[ii]=nd[ii,1]
                        else:
                            a_co=(nd[ii+1,1]-nd[ii,1])/(nd[ii+1,0]-nd[ii,0])
                            b_co=nd[ii,1]-a_co*nd[ii,0]
                            xicross[ii]=a_co*px+b_co
                    xicross=np.array([x for x in xicross if (x<=1)and(x>=-1)])
                    if len(xicross)==2:
                        xilocal_min=min(xicross)
                        xilocal_max=max(xicross)
                                    
                        qx=np.zeros((3,2))
                        # 1D quadrature points
                        qx[:,1]=[(1+sqrt(0.6))/2*xilocal_min+(1-sqrt(0.6))/2*xilocal_max,(xilocal_min+xilocal_max)/2,(1-sqrt(0.6))/2*xilocal_min+(1+sqrt(0.6))/2*xilocal_max];
                        qx[:,0]=px
                        # 1D integral weights
                        qw=np.array([5/18,8/18,5/18])*(xilocal_max-xilocal_min)
                        shapes=mesh.eval_shape_values(kk,qx)
                        grads=mesh.eval_shape_grads(kk,qx)
                        kperprho=k*qx[:,0]*sqrt(1-qx[:,1]**2)*cc/omegace0
                        Q=(Ex*n*omegace0/gammax/k/(px/gammax)/cc*sp.special.jn(n,kperprho)-Ey*sqrt(1-qx[:,1]**2)*(sp.special.jn(n+1,kperprho)-sp.special.jn(n-1,kperprho))/2.0)**2
                        kperprho2=k2*qx[:,0]*sqrt(1-qx[:,1]**2)*cc/omegace0
                        Q2=(Ez*qx[:,1]*sp.special.jn(n,kperprho2))**2
                        for ii in range(6):
                            if mesh.boundaryhigh[mesh.trisToNodes[kk][ii]] :
                                continue
                            ECKI_local=omegape0**2*omega*pi*sum((grads[:,ii,0]*px/gammax-grads[:,ii,1]*(n*omegace0/gammax-omega*(1-qx[:,1]**2))/omega/qx[:,1]/gammax)*Q*px**2/(n*omegace0*px/gammax**3)*qw)/den
                            #ECKI_local=omegape**2*omega*pi/(sqrt(2*pi)*vth**3)*sum((grads[:,ii,0]*px/gammax)*Q*px**2/(n*omegace0*px/gammax**3)*qw)/den
                            ECK_local=-2*omegape0**2*omega*pi*sum(shapes[:,ii]*(px/gammax)**2*Q*px**2/(n*omegace0*px/gammax**3)*qw)/den
                            #if mesh.trisToNodes[kk][ii]==11858:
                            #    ipdb.set_trace()
                            ECKI_row[istep]=i
                            ECKI_col[istep]=mesh.trisToNodes[kk][ii]
                            ECKI_data[istep]=ECKI_local
                            ECK_data[istep]=ECK_local
                            ECKI_local=omegape0**2*omega*pi*sum((grads[:,ii,0]*px/gammax-grads[:,ii,1]*(n*omegace0/gammax-omega*(1-qx[:,1]**2))/omega/qx[:,1]/gammax)*Q2*px**2/(n*omegace0*px/gammax**3)*qw)/den2
                            ECK_local=-2*omegape0**2*omega*pi*sum(shapes[:,ii]*(px/gammax)**2*Q2*px**2/(n*omegace0*px/gammax**3)*qw)/den2
                            ECKIO_data[istep]=ECKI_local
                            ECKO_data[istep]=ECK_local
                             
                            istep=istep+1

        ECKI_sparse = sp.sparse.coo_matrix((ECKI_data[0:istep],(ECKI_row[0:istep],ECKI_col[0:istep])),shape=(numomega*numr,len(mesh.nodes)))
        ECKI_sparse = sp.sparse.csc_matrix(ECKI_sparse)
        ECK_sparse = sp.sparse.coo_matrix((ECK_data[0:istep],(ECKI_row[0:istep],ECKI_col[0:istep])),shape=(numomega*numr,len(mesh.nodes)))
        ECK_sparse = sp.sparse.csc_matrix(ECK_sparse)
        ECKIO_sparse = sp.sparse.coo_matrix((ECKIO_data[0:istep],(ECKI_row[0:istep],ECKI_col[0:istep])),shape=(numomega*numr,len(mesh.nodes)))
        ECKIO_sparse = sp.sparse.csc_matrix(ECKIO_sparse)
        ECKO_sparse = sp.sparse.coo_matrix((ECKO_data[0:istep],(ECKI_row[0:istep],ECKI_col[0:istep])),shape=(numomega*numr,len(mesh.nodes)))
        ECKO_sparse = sp.sparse.csc_matrix(ECKO_sparse)
         
        KIM_vec=np.zeros(numr*numomega)
        KM_vec=np.zeros(numr*numomega)
        for i in range(numr*numomega):
            f0=exp(-(sqrt(1+mesh.nodes[:,0]**2)-1)/vtharray[i]**2)/(vtharray[i]**3*sqrt(2*pi))
            KIM_vec[i]=ECKI_sparse[i,:].dot(f0)
            KM_vec[i]=KIM_vec[i]*vtharray[i]**2*2
        KIMO_vec=np.zeros(numr*numomega)
        KMO_vec=np.zeros(numr*numomega)
        for i in range(numr*numomega):
            f0=exp(-(sqrt(1+mesh.nodes[:,0]**2)-1)/vtharray[i]**2)/(vtharray[i]**3*sqrt(2*pi))
            KIMO_vec[i]=ECKIO_sparse[i,:].dot(f0)
            KMO_vec[i]=KIMO_vec[i]*vtharray[i]**2*2

        np.save(data_dir + 'ECE_r_array.npy',r_array)
        sp.sparse.save_npz(data_dir + 'ECE_KImatrix.npz',ECKI_sparse)
        sp.sparse.save_npz(data_dir + 'ECE_Kmatrix.npz',ECK_sparse)
        np.save(data_dir + 'ECE_KIM_vec.npy',KIM_vec)
        np.save(data_dir + 'ECE_KM_vec.npy',KM_vec)
        sp.sparse.save_npz(data_dir + 'ECE_KIOmatrix.npz',ECKIO_sparse)
        sp.sparse.save_npz(data_dir + 'ECE_KOmatrix.npz',ECKO_sparse)
        np.save(data_dir + 'ECE_KIMO_vec.npy',KIMO_vec)
        np.save(data_dir + 'ECE_KMO_vec.npy',KMO_vec)

    #f0=exp(-(sqrt(1+mesh.nodes[:,0]**2)-1)/vth**2)/(vth**3*sqrt(2*pi))
    #ki_array=ECKI_sparse.dot(f0)
    #K_array=ECK_sparse.dot(f0)
    '''开始读入分布数据'''
    fil=open(data_dir + 'data.npy','rb')   #读取data.npy
    time_array=[]
    f_array=[]
    E_array=[]
    current_array=[]
    if enable_whistler:
        whistler.Am_array=[]
    if enable_exel:
        exel.Am_array=[]
    istep=0
    #读取data.npy
    while True:
        if istep>0:
            fil.seek(-1,1)
        istep=istep+1
        time_array.append(np.load(fil))
        f_array.append(np.load(fil))
        E_array.append(np.load(fil))
        current_array.append(np.load(fil))
        if enable_whistler:
            whistler.Am_array.append(np.load(fil))
        if enable_exel:
            exel.Am_array.append(np.load(fil))
        xx=fil.read(1)
        if xx==b'':
            break
    fil.close()

    Power=np.zeros((len(f_array),numomega))
    ref_fac=0.76
    pol_fac=0.20
    #for istep in [1200]:
    for istep in range(len(f_array)):
        ki_array=ECKI_sparse.dot(f_array[istep])
        K_array=ECK_sparse.dot(f_array[istep])
        kio_array=ECKIO_sparse.dot(f_array[istep])
        Ko_array=ECKO_sparse.dot(f_array[istep])
        no_re=(1.8*(1-((r_array-R)/a)**2)**2+0.2)<1
        #no_re[:]=True
        ki_array[no_re]=KIM_vec[no_re]
        K_array[no_re]=KM_vec[no_re]
        kio_array[no_re]=KIMO_vec[no_re]
        Ko_array[no_re]=KMO_vec[no_re]


        sx=np.zeros(numomega*numr)
        so=np.zeros(numomega*numr)
        E_x=np.zeros(numomega*numr)
        E_o=np.zeros(numomega*numr)
        toleft=True
        PP=0
        for ii in range(1):
            for i in range(numomega):
                if toleft:
                    if ii==0:
                        E_x_new=1
                        E_o_new=0
                    else:
                        E_x_new=sqrt((E_x[i*numr+numr-1]**2*(1-pol_fac)+E_o[i*numr+numr-1]**2*pol_fac)*ref_fac)
                        E_o_new=sqrt((E_o[i*numr+numr-1]**2*(1-pol_fac)+E_x[i*numr+numr-1]**2*pol_fac)*ref_fac)
                    E_x[i*numr+numr-1]=E_x_new
                    E_o[i*numr+numr-1]=E_o_new
                    if E_x[i*numr+numr-1]<exp(-5):
                        sx[i*numr:i*numr+numr]=1000
                    else:
                        sx[i*numr+numr-1]=-log(E_x[i*numr+numr-1])
                        for j in range(numr-2,-1,-1):
                            sx[i*numr+j]=sx[i*numr+j+1]+(ki_array[i*numr+j]+ki_array[i*numr+j+1])/2*(r_array[i*numr+j+1]-r_array[i*numr+j])
                        E_x=exp(-sx)
                    if E_o[i*numr+numr-1]<exp(-5):
                        so[i*numr:i*numr+numr]=1000
                    else:
                        so[i*numr+numr-1]=-log(E_o[i*numr+numr-1])
                        for j in range(numr-2,-1,-1):
                            so[i*numr+j]=so[i*numr+j+1]+(kio_array[i*numr+j]+kio_array[i*numr+j+1])/2*(r_array[i*numr+j+1]-r_array[i*numr+j])
                else:
                    E_x_new=sqrt((E_x[i*numr]**2*(1-pol_fac)+E_o[i*numr]**2*pol_fac)*ref_fac)
                    E_o_new=sqrt((E_o[i*numr]**2*(1-pol_fac)+E_x[i*numr]**2*pol_fac)*ref_fac)
                    E_x[i*numr]=E_x_new
                    E_o[i*numr]=E_o_new
                    if E_x[i*numr]<exp(-5):
                        sx[i*numr:i*numr+numr]=1000
                    else:
                        sx[i*numr]=-log(E_x[i*numr])
                        for j in range(1,numr,1):
                            sx[i*numr+j]=sx[i*numr+j-1]+(ki_array[i*numr+j]+ki_array[i*numr+j-1])/2*(r_array[i*numr+j]-r_array[i*numr+j-1])
                        E_x=exp(-sx)
                    if E_o[i*numr]<exp(-5):
                        so[i*numr:i*numr+numr]=1000
                    else:
                        so[i*numr]=-log(E_o[i*numr])
                        for j in range(1,numr,1):
                            so[i*numr+j]=so[i*numr+j-1]+(kio_array[i*numr+j]+kio_array[i*numr+j-1])/2*(r_array[i*numr+j]-r_array[i*numr+j-1])
            E_x=exp(-sx)
            E_o=exp(-so)
            PP=PP+E_x*E_x*K_array+E_o*E_o*Ko_array
            if toleft:
                toleft=False
            else:
                toleft=True

        
        for i in range(numomega):
            s=np.zeros(numr)
            for j in range(2,numr-1):
                s[j]=(r_array[i*numr+j+1]-r_array[i*numr+j-1])/2
            s[0]=(r_array[i*numr+1]-r_array[i*numr])/2
            s[numr-1]=(r_array[i*numr+numr-1]-r_array[i*numr+numr-2])/2
            #ipdb.set_trace()
            Power[istep,i]=sum(s*PP[i*numr:i*numr+numr])

    Power=Power/(1000*1.6e-19/(9.1e-31*3e8**2))


    # 保存计算得到的功率数据
    np.save(data_dir + 'ECE_Power.npy', Power)

    # 如果需要保存更多相关信息，也可以保存时间数组等
    # np.save(data_dir + 'ECE_time_array.npy', time_array)
    # np.save(data_dir + 'ECE_omega_array.npy', omega_array)
    # np.save(data_dir + 'ECE_r_array.npy', r_array)

    print(f"功率数据已保存到 ECE_Power.npy，形状为: {Power.shape}")
    return Power


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate ECE (Electron Cyclotron Emission) power')
    parser.add_argument('--data_dir', type=str, default=None, help='Directory path for data files')
    args = parser.parse_args()
    calculate_ECE(data_dir_param=args.data_dir)
