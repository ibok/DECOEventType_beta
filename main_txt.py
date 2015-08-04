#!/usr/bin/env python
# WIPAC-DECO particle event analyzer (beta)
# Licensed under the GNU v3.0 (open-source)
YQ=float
YN=len
YS=range
Yd=min
YR=max
Yh=None
Yx=int
Ye=False
YF=True
Yv=map
YP=sorted
Ys=set
YA=abs
Yw=type
YK=str
Ya=SystemExit
YI=open
YW=xrange
YG=sum
Yi=enumerate
import argparse,math,os,sys
zT=argparse.ArgumentParser
Yb=sys.maxint
Yz=math.sqrt
import numpy as np
YO=np.vstack
Yu=np.mean
YD=np.array
YC=np.YA
YH=np.arctan2
Yf=np.YG
Yc=np.sqrt
Yn=np.sign
from PIL import Image
YV=Image.YI
from skimage import measure
Ym=measure.find_contours
class zl:
 def __init__(z,m,p):
  z.m=YQ(m)
  z.p=(YQ(p[0]),YQ(p[1]))
 def za(z,x):
  return z.m*(x-z.p[0])+z.p[1]
 def zI(z,o):
  if o.m==z.m:
   return NaN
  else:
   Y=(o.p[1]-z.p[1]+z.m*z.p[0]+o.m*o.p[0])/(z.m-o.m)
   return(Y,z.za(Y))
class zk:
 def __init__(z,x,y):
  z.x=x
  z.y=y
  z.xc=Yu(x)
  z.yc=Yu(y)
  z.area=0.
  n=YN(x)
  for i in YS(0,n):
   z.area+=0.5*(y[i]+y[i-1])*(x[i]-x[i-1])
  z.perimeter=0
  b=YN(z.x)-1
  for i in YS(0,b):
   a=(z.x[i],z.y[i])
   b=(z.y[i+1],z.y[i+1])
   z.perimeter+=zW(a,b)
  z.maxdist=0
  b=YN(z.x)
  for i in YS(0,b):
   for j in YS(i,b):
    u=(z.x[i],z.y[i])
    c=(z.y[j],z.y[j])
    C=zW(u,c)
    if C>z.maxdist:
     z.maxdist=C
 def zW(z,f):
  return Yc((z.xc-f.xc)**2+(z.yc-f.yc)**2)
def zW(a,b):
 return Yz((a[1]-b[1])**2+(a[0]-b[0])**2)
class zo:
 def __init__(z):
  z.blobs=[]
  z.count=0
  z.area=0.
  z.b_area=0.
  z.perimeter=0.
  z.xmin= 1e10
  z.xmax=-1e10
  z.ymin= 1e10
  z.ymax=-1e10
 def zG(z,f):
  z.blobs.append(f)
  z.xmin=Yd(z.xmin,f.x.Yd())
  z.xmax=YR(z.xmax,f.x.YR())
  z.ymin=Yd(z.ymin,f.y.Yd())
  z.ymax=YR(z.ymax,f.y.YR())
  z.cov =Yh
  z.count+=1
  z.b_area+=f.area
  z.perimeter+=f.perimeter
 def zi(z):
  return(z.xmin,z.xmax,z.ymin,z.ymax)
 def zL(z):
  O,H,n,D=(z.xmin,z.xmax,z.ymin,z.ymax)
  xL=YC(H-O)
  yL=YC(D-n)
  if xL>yL:
   n-=0.5*(xL-yL)
   D+=0.5*(xL-yL)
  else:
   O-=0.5*(yL-xL)
   H+=0.5*(yL-xL)
  return(O,H,n,D)
 def zX(z,Q):
  ny,nx=Q.shape
  x0,x1,y0,y1=bg.getBoundingBox()
  z.area=(x1-x0)*(y1-y0)
  i0,i1=[ny-Yx(t)for t in(y1,y0)]
  j0,j1=[Yx(t)for t in(x0,x1)]
  V=1
  i0=0 if i0-V<0 else i0-V
  i1=ny-1 if i1>ny-1 else i1+V
  j0=0 if j0-V<0 else j0-V
  j1=nx-1 if j1>nx-1 else j1+V
  return Q[i0:i1,j0:j1]
 def zM(z,Q,p,q):
  nx,ny=Q.shape
  m=0.
  if p==0 and q==0:
   m=Yf(Q)
  else:
   for i in YS(0,nx):
    x=0.5+i
    for j in YS(0,ny):
     y=0.5+j
     m+=x**p*y**q*Q[i,j]
  return m
 def zg(z,Q):
  if z.cov==Yh:
   N=z.zX(Q).transpose()
   S=z.zM(N,0,0)
   d=z.zM(N,1,0)
   R=z.zM(N,0,1)
   h=z.zM(N,1,1)
   x=z.zM(N,2,0)
   e=z.zM(N,0,2)
   F=d/S
   v=R/S
   z.cov=YO([[x/S-F*F,h/S-F*v],[h/S-F*v,e/S-v*v]])
  return z.cov
 def zE(z,Q):
  P=z.zg(Q)
  s=P[0,0]
  A=P[0,1]
  w=P[1,1]
  K=0.5*YH(2*A,s-w)
  l1=0.5*(s+w)+0.5*Yc(4*A**2+(s-w)**2)
  l2=0.5*(s+w)-0.5*Yc(4*A**2+(s-w)**2)
  return l1,l2,l1/l2,K
def zr(Q,threshold,minArea=2.):
 a=[]
 ny,nx=Q.shape
 I=Ym(Q,threshold)
 for W in I:
  x=W[:,1]
  y=ny-W[:,0]
  f=zk(x,y)
  if f.area>=minArea:
   a.append(f)
 return a
def zJ(a,maxDist):
 n=YN(a)
 G=[]
 if n>=1:
  bg=zo()
  bg.addBlob(a[0])
  G.append(bg)
  for i in YS(1,n):
   bi=a[i]
   i=Ye
   for L in G:
    for bj in L.blobs:
     if bi.distance(bj)<maxDist:
      L.zG(bi)
      i=YF
      break
   if not i:
    bg=zo()
    bg.addBlob(bi)
    G.append(bg)
 return G
def zq(L):
 X=[]
 M=[0]
 g=[0,0,0]
 E=100.
 r=YF
 for b in L.blobs:
  for i in YS(0,YN(b.x)):
   X.append([YQ(b.x[i]),YQ(b.y[i])])
 for i in YS(0,YN(X)):
  for j in YS(i,YN(X)):
   if(zW(X[i],X[j])>M[0]):
    M=[zW(X[i],X[j]),X[i],X[j]]
 if M[1][0]<M[2][0]:
  J=M[1][0]
  y1=M[1][1]
  q=M[2][0]
 else:
  J=M[2][0]
  y1=M[2][1]
  q=M[1][0]
  r=Ye
 U=M[1][1]-M[2][1]
 B=M[1][0]-M[2][0]
 if r:
  return[Yv(Yx,M[1]),Yv(Yx,M[2])]
 else:
  return[Yv(Yx,M[2]),Yv(Yx,M[1])]
def zU(L):
 X=[]
 for b in L.blobs:
  for i in YS(0,YN(b.x)):
   X.append((YQ(b.x[i]),YQ(b.y[i])))
 return YP(Ys(X))
def zB(U,B):
 if B==0 and U==0:
  return NaN
 elif B==0:
  return Yn(U)*Yc(Yb)
 else:
  return YQ(U)/B
def zt(o,point):
 return YA(-o.m*point[0]+point[1]-(o.p[1]-o.m*o.p[0]))/YQ(Yc((-o.m)**2+1))
def zp(zC,Yw):
 print>>f,YK(zC)+','+zj(Yw,0)
def zy():
 raise Ya
p=zT(description="Histogram luminance of a JPG image")
p.add_argument("text_file",nargs=1,help="Image filepath names")
p.add_argument("-x","--threshold",dest="thresh",Yw=YQ,default=Yh,help="Threshold the pixel map")
y=p.add_argument_group("Blob Identification")
y.add_argument("-c","--contours",dest="contours",Yw=YQ,default=Yh,help="Identify blob contours above some threshold")
y.add_argument("-d","--distance",dest="distance",Yw=YQ,default=75.,help="Group blobs within some distance of each other")
y.add_argument("-a","--min-area",dest="area",Yw=YQ,default=2.,help="Remove blobs below some minimum area")
j=p.parse_args()
l=YI(j.text_file[0],'r+')
k=[]
for o in l:
 k.append(o.strip('\n'))
T=[]
zY=[]
if YN(k)==0:
 print('No image files were found.')and zy()
f=YI('classifications.out','w')
def zj(x,opt):
 if opt==0:
  return{1:'spot',2:'worm',3:'track',4:'ambig',5:'big_spot',6:'track_lowconf',7:'noise'}[x]
 else:
  return{1:'_s',2:'_w',3:'_t',4:'_a',5:'_b',6:'_l',7:'_n'}[x]
for zb in k:
 zu=zb[1].replace('//','/')
 zc=zu.split('/')[-1]
 zC=zc.split('.')[0]
 zf=zc.split('.')[-1]
 zO=YV(zu).convert("L")
 Q=[]
 zH=zO.load()
 nx=zO.size[0]
 ny=zO.size[1]
 x0,y0,x1,y1=(0,0,nx,ny)
 for y in YW(ny):
  Q.append([zH[x,y]for x in YW(nx)])
 Q=YD(Q,dtype=YQ)
 zn=YG(YG(Q))/(YN(Q)*YN(Q[0]))
 if zn>4:
  zp(zC,7)
  continue
 a=zr(Q,threshold=args.contours,minArea=args.area)
 G=zJ(a,maxDist=args.distance)
 for g in G:
  T.append(zq(g))
  zY.append(zU(g))
 if j.thresh!=Yh:
  Q[Q<j.thresh]=0.
 if j.contours!=Yh:
  for i,bg in Yi(G):
   X0,X1,Y0,Y1=bg.getSquareBoundingBox()
   l1,l2,r,K=bg.getPrincipalMoments(Q)
   zD=(Yc(l1**2-l2**2)/l1)
   zV=T[i][0]
   zm=T[i][1]
   _i=zl(zB((YQ(T[i][0][1])-T[i][1][1]),(YQ(T[i][0][0])-T[i][1][0])),zV)
   zQ=0.
   zN=zW(zV,zm)
   for pt in zY[i]:
    zQ+=zt(_i,pt)
   zS=zQ/zN
   zd=0.
   zR=[0,0]
   for b in bg.blobs:
    zd+=b.perimeter
    zR[0]+=b.xc*b.perimeter
    zR[1]+=b.yc*b.perimeter
   zh=zR[0]/YQ(bg.count*zd)
   zx=zR[1]/YQ(bg.count*zd)
   ze=zh-bg.xmin
   zF=bg.xmax-zh
   zv=zx-bg.ymin
   zP=bg.ymax-zx
   zs=Yd([ze,zF])/YR([ze,zF])
   zA=Yd([zv,zP])/YR([zv,zP])
   zw=zs*zA
   zK=YA(zs-zA)
   if j.contours==40:
    if zD>.99993 and l1>700:
     Yw=4
    elif bg.b_area<4 or zN<6 or zN<13 and(r>=.2 and zD<.945 or bg.b_area<7)or zD<.7:
     if bg.b_area>50:
      if bg.b_area>62 and zN>10:
       Yw=2
      else:
       Yw=5
     else:
      Yw=1
    else:
     if zK>.55:
      Yw=2
     elif bg.b_area>100 and(l1>100 and l1/10>l2)and zN>30:
      if(zS>9)or(l1/5>zN and zS>3.9)or(80>zN>40 and bg.b_area>100 and zS>5):
       Yw=2
      else:
       Yw=3
     elif zD>.9998 and 130>zN>90:
      Yw=3
     elif zD>.9995 and zN>40 and zw>.8:
      Yw=4
     else:
      if(zw>.978 and zK<.01 or zw>.96 and zK<.0047 or zw>.9 and r<.02 and zK<.055 and zS<5)and zD>.96:
       if bg.b_area>33:
        Yw=3
       else:
        Yw=2
      elif(r<.22 and(zD<.985 and(zK<.015 or zw>.88)and 18>zN>9 or .97<zD<.985 and 18>zN>9 and zw>.83 and bg.b_area<30 or .99>zD>.975 and zS<3.7 and 18>zN>7.6 and zK<.023)or zD>.99 and l1<15 and zw>.86 and zK<.1 and 35>bg.b_area>28):
       Yw=3
      else:
       if(zS>4.6 and bg.b_area<100)or bg.b_area<24 or zD<.979 or(r>.2 and zS>4):
        Yw=2
       else:
        if zK<.7 and zS>6:
         Yw=2
        elif(l1>100 and l2<12 and bg.b_area>40 and bg.b_area<60 and zS<3.8 or(YA(zs)>.99 or YA(zA)>.99)and YA(zs)>.93 and YA(zA)>.93 and zD>.99 and zS<2.95 and zS>2. and zK>.05 or(zw>.9 and zK<.02 and zS<3.1 and zD>.993 and zN>12)and bg.b_area<82 or((zw>.6 and zD>.9923 and zS<3.1 or zw>.88)and(zK<.03 or YA(zs)>.996 or YA(zA)>.996)and bg.b_area<100)):
         Yw=3
        else:
         if zD>.999 and zS<3.14 and bg.b_area>58 and l1/25>l2:
          Yw=3
         elif zD>.999 and .86<zw<.92 and zN>23:
          Yw=2
         elif(zD>.999 and((l1>90 and l2<10)and(zS>2.9 or zS<1.1)or zD>.992 and bg.b_area<50 and .98>YA(zs)>.96 and YA(zA)>.96)):
          Yw=4
         elif zw>.75 and zK<.182 and((bg.b_area>28)or(bg.b_area<28 and zN>17)):
          if((zD>.9996 or r<.028)and zw<.9 and zN<30 and bg.b_area<62 or zD>.99 and l1>400 and l1<600 and l2>60 and zS>3.4 or zD<.99 and zD>.975 and zN<17 and l1<16 and l2>2 and r>.2 or zD>.993 and(zS<3 and 28<zN<40 and .94>zw>.9 and bg.b_area<50 or 3.5<zS<4 and 17<zN<25 and r<.12)):
           Yw=2
          elif(zS<3.76 and zD>.99 and zK<.06 and r<.13 and(bg.b_area>60 or zN>10)and YR(YA(zs),YA(zA))>.935 and Yd(YA(zs),YA(zA))>.875 or zS<4.1 and bg.b_area>30 and zK<0.059 and zN<16 or(zS<4.16 and zw>.74 and zK<.012 and bg.b_area<50 and zN<20 and 12<l1<23 and l2<3)):
           Yw=3
          else:
           Yw=2
         elif zw>.75 and zK<.05 and bg.b_area>30:
          Yw=6
         elif .45<zw<.6 and .2<zK<.5 and .999>zD>.92:
          Yw=3
         else:
          Yw=2
    zp(zC,Yw)
   else:
    print('As of now only image analysis at contour level 40 is supported, sorry.')and zy()
# Created by pyminifier (https://github.com/liftoff/pyminifier)
