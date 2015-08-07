#!/usr/bin/env python
DY=float;DX=len;DV=range;Do=min;DO=max;Dl=None;Ds=int;Di=False;Dz=True;Da=map
DB=abs;Dy=type;Du=str;DM=SystemExit;Df=open;DW=xrange;Dd=sum
DA=enumerate
import argparse,math,os,sys
Dt=argparse.ArgumentParser
DP=os.walk
DC=sys.maxint
Dw=math.sqrt
import numpy as np
DS=np.vstack
DN=np.mean
Dn=np.array
Dr=np.DB
DK=np.arctan2
Dv=np.Dd
DE=np.sqrt
Dq=np.sign
from PIL import Image
Dg=Image.Df
from skimage import measure
DR=measure.find_contours
class tU:
 def __init__(t,m,p):
  t.m=DY(m)
  t.p=(DY(p[0]),DY(p[1]))
 def tu(t,x):
  return t.m*(x-t.p[0])+t.p[1]
 def tM(t,line):
  if line.m==t.m:
   return NaN
  else:
   D=(line.p[1]-t.p[1]+t.m*t.p[0]+line.m*line.p[0])/(t.m-line.m)
   return(D,t.tu(D))
class th:
 def __init__(t,x,y):
  t.x=x
  t.y=y
  t.xc=DN(x)
  t.yc=DN(y)
  t.area=0.
  n=DX(x)
  for i in DV(0,n):
   t.area+=0.5*(y[i]+y[i-1])*(x[i]-x[i-1])
  t.perimeter=0
  w=DX(t.x)-1
  for i in DV(0,w):
   a=(t.x[i],t.y[i])
   b=(t.y[i+1],t.y[i+1])
   t.perimeter+=tf(a,b)
  t.maxdist=0
  w=DX(t.x)
  for i in DV(0,w):
   for j in DV(i,w):
    P=(t.x[i],t.y[i])
    C=(t.y[j],t.y[j])
    N=tf(P,C)
    if N>t.maxdist:
     t.maxdist=N
 def tf(t,E):
  return DE((t.xc-E.xc)**2+(t.yc-E.yc)**2)
def tf(a,b):
 return Dw((a[1]-b[1])**2+(a[0]-b[0])**2)
class tI:
 def __init__(t):
  t.blobs=[]
  t.count=0
  t.area=0.
  t.b_area=0.
  t.perimeter=0.
  t.xmin= 1e10
  t.xmax=-1e10
  t.ymin= 1e10
  t.ymax=-1e10
 def tW(t,E):
  t.blobs.append(E)
  t.xmin=Do(t.xmin,E.x.Do())
  t.xmax=DO(t.xmax,E.x.DO())
  t.ymin=Do(t.ymin,E.y.Do())
  t.ymax=DO(t.ymax,E.y.DO())
  t.cov =Dl
  t.count+=1
  t.b_area+=E.area
  t.perimeter+=E.perimeter
 def td(t):
  return(t.xmin,t.xmax,t.ymin,t.ymax)
 def tA(t):
  r,v,S,K=(t.xmin,t.xmax,t.ymin,t.ymax)
  xL=Dr(v-r)
  yL=Dr(K-S)
  if xL>yL:
   S-=0.5*(xL-yL)
   K+=0.5*(xL-yL)
  else:
   r-=0.5*(yL-xL)
   v+=0.5*(yL-xL)
  return(r,v,S,K)
 def tF(t,g):
  ny,nx=g.shape
  x0,x1,y0,y1=bg.getBoundingBox()
  t.area=(x1-x0)*(y1-y0)
  i0,i1=[ny-Ds(t)for t in(y1,y0)]
  j0,j1=[Ds(t)for t in(x0,x1)]
  q=1
  i0=0 if i0-q<0 else i0-q
  i1=ny-1 if i1>ny-1 else i1+q
  j0=0 if j0-q<0 else j0-q
  j1=nx-1 if j1>nx-1 else j1+q
  return g[i0:i1,j0:j1]
 def tx(t,g,p,q):
  nx,ny=g.shape
  n=0.
  if p==0 and q==0:
   n=Dv(g)
  else:
   for i in DV(0,nx):
    x=0.5+i
    for j in DV(0,ny):
     y=0.5+j
     n+=x**p*y**q*g[i,j]
  return n
 def tm(t,g):
  if t.cov==Dl:
   R=t.tF(g).transpose()
   Y=t.tx(R,0,0)
   X=t.tx(R,1,0)
   V=t.tx(R,0,1)
   o=t.tx(R,1,1)
   O=t.tx(R,2,0)
   l=t.tx(R,0,2)
   s=X/Y
   i=V/Y
   t.cov=DS([[O/Y-s*s,o/Y-s*i],[o/Y-s*i,l/Y-i*i]])
   z=s+t.xmin
   a=i+t.ymin
   print>>xy,tn,z,a
  return t.cov
 def tj(t,g):
  B=t.tm(g)
  y=B[0,0]
  u=B[0,1]
  M=B[1,1]
  f=0.5*DK(2*u,y-M)
  l1=0.5*(y+M)+0.5*DE(4*u**2+(y-M)**2)
  l2=0.5*(y+M)-0.5*DE(4*u**2+(y-M)**2)
  return l1,l2,l1/l2,f
def te(g,threshold,minArea=2.):
 W=[]
 ny,nx=g.shape
 d=DR(g,threshold)
 for A in d:
  x=A[:,1]
  y=ny-A[:,0]
  E=th(x,y)
  if E.area>=minArea:
   W.append(E)
 return W
def tH(W,maxDist):
 n=DX(W)
 F=[]
 if n>=1:
  bg=tI()
  bg.addBlob(W[0])
  F.append(bg)
  for i in DV(1,n):
   bi=W[i]
   x=Di
   for m in F:
    for bj in m.blobs:
     if bi.distance(bj)<maxDist:
      m.tW(bi)
      x=Dz
      break
   if not x:
    bg=tI()
    bg.addBlob(bi)
    F.append(bg)
 return F
def tc(m):
 j=[]
 e=[0]
 for b in m.blobs:
  for i in DV(0,DX(b.x)):
   j.append((b.x[i],b.y[i]))
 for i in DV(0,DX(j)):
  for j in DV(i,DX(j)):
   if(tf(j[i],j[j])>e[0]):
    e=[tf(j[i],j[j]),j[i],j[j]]
 H=e[1]if e[1][0]<e[2][0]else e[2]
 c=e[1]if not e[1]==H else e[2]
 return(Da(Ds,H),Da(Ds,c))
def tT(m):
 j=[]
 for b in m.blobs:
  for i in DV(0,DX(b.x)-1):
   j.append((DY(b.x[i]),DY(b.y[i])))
 return j
def tG(num,dem):
 if dem==0 and num==0:
  return NaN
 elif dem==0:
  return Dq(num)*DE(DC)
 else:
  return DY(num)/dem
def tk(bg):
 T=0.
 G=[0,0]
 for b in bg.blobs:
  T+=b.perimeter
  G[0]+=b.xc*b.perimeter
  G[1]+=b.yc*b.perimeter
 k=G[0]/(bg.count*T)
 L=G[1]/(bg.count*T)
 Q=k-bg.xmin
 b=bg.xmax-k
 p=L-bg.ymin
 J=bg.ymax-L
 U=Do(Q,b)/DO(Q,b)
 h=Do(p,J)/DO(p,J)
 return U*h,DB(U-h)
def tL(line,point):
 return DB(-line.m*point[0]+point[1]-(line.p[1]-line.m*line.p[0]))/DY(DE((-line.m)**2+1))
def tQ(tn,Dy):
 print>>f,Du(tn)+','+tp(Dy)
def tb():
 raise DM
p=Dt(description="Histogram luminance of a JPG image")
p.add_argument("folder",nargs=1,help="JPG folder name")
p.add_argument("-x","--threshold",dest="thresh",Dy=DY,default=Dl,help="Threshold the pixel map")
tw=p.add_argument_group("Blob Identification")
tw.add_argument("-c","--contours",dest="contours",Dy=DY,default=Dl,help="Identify blob contours above some threshold")
tw.add_argument("-d","--distance",dest="distance",Dy=DY,default=100.,help="Group blobs within some distance of each other")
tw.add_argument("-a","--min-area",dest="area",Dy=DY,default=2.,help="Remove blobs below some minimum area")
tP=p.parse_args()
if not tP.contours==40:
 print('Only image analysis at contour level 40 is supported.')and tb()
tC=DP(tP.folder[0])
tN=[]
tE=('jpg jpeg bmp png eps gif im j2k j2p jpx msp pcx png ppm pgm pbm'+'spi tiff webp xbm xv').split(' ')
for tr,afile,i in tC:
 for tv in i:
  if tv.split('.')[-1]in tE:
   tN.append((tr,tr+'/'+tv))
if DX(tN)==0:
 print('No image files were found.')and tb()
f=Df('classifications.out','w')
xy=Df('xandyCent.out','w')
def tp(x):
 return{1:'spot',2:'worm',3:'ambig',4:'track',5:'noise'}[x]
def tJ(old,new):
 return new if new>old else old
for tS in tN:
 tK=tS[1].replace('//','/')
 tq=tK.split('/')[-1]
 tn=tq.split('.')[0]
 tg=tq.split('.')[-1]
 tR=Dg(tK).convert("L")
 g=[]
 tY=tR.load()
 nx=tR.size[0]
 ny=tR.size[1]
 x0,y0,x1,y1=(0,0,nx,ny)
 for y in DW(ny):
  g.append([tY[x,y]for x in DW(nx)])
 g=Dn(g,dtype=DY)
 tX=Dd(Dd(g))/(DX(g)*DX(g[0]))
 if tX>4:
  tQ(tn,5)
  continue
 W=te(g,threshold=args.contours,minArea=args.area)
 F=tH(W,maxDist=args.distance)
 tV=[]
 to=[]
 for g in F:
  tV.append(tc(g))
  to.append(tT(g))
 if tP.thresh!=Dl:
  g[g<tP.thresh]=0.
 Dy=0
 for i,bg in DA(F):
  l1,l2,r,f=bg.getPrincipalMoments(g)
  tO=(DE(l1**2-l2**2)/l1)
  tl=tV[i][0]
  ts=tV[i][1]
  _i=tU(tG((DY(tV[i][0][1])-tV[i][1][1]),(DY(tV[i][0][0])-tV[i][1][0])),tl)
  ti=0.
  tz=tf(tl,ts)
  for pt in to[i]:
   ti+=tL(_i,pt)
  ta=ti/tz
  tB,ty=tk(bg)
  if tO>.99993 and l1>700:
   Dy=tJ(Dy,3)
  elif bg.b_area<4 or tz<6 or tz<13 and(r>=.2 and tO<.945 or bg.b_area<7)or tO<.7:
   if bg.b_area>50:
    if bg.b_area>62 and tz>10:
     Dy=tJ(Dy,2)
    else:
     Dy=tJ(Dy,1)
   else:
    Dy=tJ(Dy,1)
  else:
   if ty>.55:
    Dy=tJ(Dy,2)
   elif bg.b_area>100 and(l1>100 and l1/10>l2)and tz>30:
    if(ta>9)or(l1/5>tz and ta>3.9)or(80>tz>40 and bg.b_area>100 and ta>5):
     Dy=tJ(Dy,2)
    else:
     Dy=4
   elif tO>.9998 and 130>tz>90:
    Dy=4
   elif tO>.9995 and tz>40 and tB>.8:
    Dy=tJ(Dy,3)
   else:
    if(tB>.978 and ty<.01 or tB>.96 and ty<.0047 or tB>.9 and r<.02 and ty<.055 and ta<5)and tO>.96:
     if bg.b_area>33:
      Dy=4
     else:
      Dy=tJ(Dy,2)
    elif(r<.22 and(tO<.985 and(ty<.015 or tB>.88)and 18>tz>9 or .97<tO<.985 and 18>tz>9 and tB>.83 and bg.b_area<30 or .99>tO>.975 and ta<3.7 and 18>tz>7.6 and ty<.023)or tO>.99 and l1<15 and tB>.86 and ty<.1 and 35>bg.b_area>28):
     Dy=4
    else:
     if(ta>4.6 and bg.b_area<100)or bg.b_area<24 or tO<.979 or(r>.2 and ta>4):
      Dy=tJ(Dy,2)
     else:
      if ty<.7 and ta>6:
       Dy=tJ(Dy,2)
      elif(l1>100 and l2<12 and bg.b_area>40 and bg.b_area<60 and ta<3.8 or(DB(U)>.99 or DB(h)>.99)and DB(U)>.93 and DB(h)>.93 and tO>.99 and ta<2.95 and ta>2. and ty>.05 or(tB>.9 and ty<.02 and ta<3.1 and tO>.993 and tz>12)and bg.b_area<82 or((tB>.6 and tO>.9923 and ta<3.1 or tB>.88)and(ty<.03 or DB(U)>.996 or DB(h)>.996)and bg.b_area<100)):
       Dy=4
      else:
       if tO>.999 and ta<3.14 and bg.b_area>58 and l1/25>l2:
        Dy=4
       elif tO>.999 and .86<tB<.92 and tz>23:
        Dy=tJ(Dy,2)
       elif(tO>.999 and((l1>90 and l2<10)and(ta>2.9 or ta<1.1)or tO>.992 and bg.b_area<50 and .98>DB(U)>.96 and DB(h)>.96)):
        Dy=tJ(Dy,3)
       elif tB>.75 and ty<.182 and((bg.b_area>28)or(bg.b_area<28 and tz>17)):
        if((tO>.9996 or r<.028)and tB<.9 and tz<30 and bg.b_area<62 or tO>.99 and l1>400 and l1<600 and l2>60 and ta>3.4 or tO<.99 and tO>.975 and tz<17 and l1<16 and l2>2 and r>.2 or tO>.993 and(ta<3 and 28<tz<40 and .94>tB>.9 and bg.b_area<50 or 3.5<ta<4 and 17<tz<25 and r<.12)):
         Dy=tJ(Dy,2)
        elif(ta<3.76 and tO>.99 and ty<.06 and r<.13 and(bg.b_area>60 or tz>10)and DO(DB(U),DB(h))>.935 and Do(DB(U),DB(h))>.875 or ta<4.1 and bg.b_area>30 and ty<0.059 and tz<16 or(ta<4.16 and tB>.74 and ty<.012 and bg.b_area<50 and tz<20 and 12<l1<23 and l2<3)):
         Dy=4
        else:
         Dy=tJ(Dy,2)
       elif tB>.75 and ty<.05 and bg.b_area>30:
        Dy=4
       elif .45<tB<.6 and .2<ty<.5 and .999>tO>.92:
        Dy=4
       else:
        Dy=tJ(Dy,2)
  if Dy==4:
   break
 tQ(tn,Dy)
# Created by pyminifier (https://github.com/liftoff/pyminifier)
