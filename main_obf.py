#!/usr/bin/env python
# WIPAC DECO event labeling script. GNU 3.0 Open-sourced
ws=float;wQ=len;wj=range;wS=min;wU=max;wT=None;wd=int
wK=False;wB=True;wE=map;wR=sorted;wl=set;wb=abs;wg=type
wo=str;wM=SystemExit;wk=open;wx=xrange;wD=sum;wY=enumerate
import argparse,math,os,sys
wH=argparse.ArgumentParser
wJ=sys.maxint
wa=math.sqrt
import numpy as np
wi=np.vstack
wP=np.mean
wm=np.array
wz=np.wb
wu=np.arctan2
wy=np.wD
wn=np.sqrt
wF=np.sign
from PIL import Image
wV=Image.wk
from skimage import measure
we=measure.find_contours
class tC:
 def __init__(t,m,p):
  t.m=ws(m)
  t.p=(ws(p[0]),ws(p[1]))
 def to(t,x):
  return t.m*(x-t.p[0])+t.p[1]
 def tM(t,tH):
  if tH.m==t.m:
   return NaN
  else:
   w=(tH.p[1]-t.p[1]+t.m*t.p[0]+tH.m*tH.p[0])/(t.m-tH.m)
   return(w,t.to(w))
class tX:
 def __init__(t,x,y):
  t.x=x
  t.y=y
  t.xc=wP(x)
  t.yc=wP(y)
  t.area=0.
  n=wQ(x)
  for i in wj(0,n):
   t.area+=0.5*(y[i]+y[i-1])*(x[i]-x[i-1])
  t.perimeter=0
  H=wQ(t.x)-1
  for i in wj(0,H):
   a=(t.x[i],t.y[i])
   b=(t.y[i+1],t.y[i+1])
   t.perimeter+=tk(a,b)
  t.maxdist=0
  H=wQ(t.x)
  for i in wj(0,H):
   for j in wj(i,H):
    a=(t.x[i],t.y[i])
    J=(t.y[j],t.y[j])
    P=tk(a,J)
    if P>t.maxdist:
     t.maxdist=P
 def tk(t,n):
  return wn((t.xc-n.xc)**2+(t.yc-n.yc)**2)
def tk(a,b):
 return wa((a[1]-b[1])**2+(a[0]-b[0])**2)
class wt:
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
 def tx(t,n):
  t.blobs.append(n)
  t.xmin=wS(t.xmin,n.x.wS())
  t.xmax=wU(t.xmax,n.x.wU())
  t.ymin=wS(t.ymin,n.y.wS())
  t.ymax=wU(t.ymax,n.y.wU())
  t.cov =wT
  t.count+=1
  t.b_area+=n.area
  t.perimeter+=n.perimeter
 def tD(t):
  return(t.xmin,t.xmax,t.ymin,t.ymax)
 def tY(t):
  z,y,i,u=(t.xmin,t.xmax,t.ymin,t.ymax)
  xL=wz(y-z)
  yL=wz(u-i)
  if xL>yL:
   i-=0.5*(xL-yL)
   u+=0.5*(xL-yL)
  else:
   z-=0.5*(yL-xL)
   y+=0.5*(yL-xL)
  return(z,y,i,u)
 def tI(t,V):
  ny,nx=V.shape
  x0,x1,y0,y1=bg.getBoundingBox()
  t.area=(x1-x0)*(y1-y0)
  i0,i1=[ny-wd(t)for t in(y1,y0)]
  j0,j1=[wd(t)for t in(x0,x1)]
  F=1
  i0=0 if i0-F<0 else i0-F
  i1=ny-1 if i1>ny-1 else i1+F
  j0=0 if j0-F<0 else j0-F
  j1=nx-1 if j1>nx-1 else j1+F
  return V[i0:i1,j0:j1]
 def th(t,V,p,q):
  nx,ny=V.shape
  m=0.
  if p==0 and q==0:
   m=wy(V)
  else:
   for i in wj(0,nx):
    x=0.5+i
    for j in wj(0,ny):
     y=0.5+j
     m+=x**p*y**q*V[i,j]
  return m
 def tq(t,V):
  if t.cov==wT:
   e=t.tI(V).transpose()
   s=t.th(e,0,0)
   Q=t.th(e,1,0)
   j=t.th(e,0,1)
   S=t.th(e,1,1)
   U=t.th(e,2,0)
   T=t.th(e,0,2)
   d=Q/s
   K=j/s
   t.cov=wi([[U/s-d*d,S/s-d*K],[S/s-d*K,T/s-K*K]])
   B=d+t.xmin
   E=K+t.ymin
   print>>xy,tz,B,E
  return t.cov
 def tc(t,V):
  R=t.tq(V)
  l=R[0,0]
  b=R[0,1]
  g=R[1,1]
  o=0.5*wu(2*b,l-g)
  l1=0.5*(l+g)+0.5*wn(4*b**2+(l-g)**2)
  l2=0.5*(l+g)-0.5*wn(4*b**2+(l-g)**2)
  return l1,l2,l1/l2,o
def tv(V,threshold,minArea=2.):
 M=[]
 ny,nx=V.shape
 k=we(V,threshold)
 for x in k:
  x=x[:,1]
  y=ny-x[:,0]
  n=tX(x,y)
  if n.area>=minArea:
   M.append(n)
 return M
def tO(M,maxDist):
 n=wQ(M)
 D=[]
 if n>=1:
  bg=wt()
  bg.addBlob(M[0])
  D.append(bg)
  for i in wj(1,n):
   bi=M[i]
   Y=wK
   for I in D:
    for bj in I.blobs:
     if bi.distance(bj)<maxDist:
      I.tx(bi)
      Y=wB
      break
   if not Y:
    bg=wt()
    bg.addBlob(bi)
    D.append(bg)
 return D
def tW(I):
 h=[]
 q=[0]
 c=[0,0,0]
 v=100.
 O=wB
 for b in I.blobs:
  for i in wj(0,wQ(b.x)):
   h.append([ws(b.x[i]),ws(b.y[i])])
 for i in wj(0,wQ(h)):
  for j in wj(i,wQ(h)):
   if(tk(h[i],h[j])>q[0]):
    q=[tk(h[i],h[j]),h[i],h[j]]
 if q[1][0]<q[2][0]:
  W=q[1][0]
  y1=q[1][1]
  G=q[2][0]
 else:
  W=q[2][0]
  y1=q[2][1]
  G=q[1][0]
  O=wK
 L=q[1][1]-q[2][1]
 f=q[1][0]-q[2][0]
 if O:
  return[wE(wd,q[1]),wE(wd,q[2])]
 else:
  return[wE(wd,q[2]),wE(wd,q[1])]
def tG(I):
 h=[]
 for b in I.blobs:
  for i in wj(0,wQ(b.x)):
   h.append((ws(b.x[i]),ws(b.y[i])))
 return wR(wl(h))
def tL(L,f):
 if f==0 and L==0:
  return NaN
 elif f==0:
  return wF(L)*wn(wJ)
 else:
  return ws(L)/f
def tf(tH,point):
 return wb(-tH.m*point[0]+point[1]-(tH.p[1]-tH.m*tH.p[0]))/ws(wn((-tH.m)**2+1))
def tN(tz,wg):
 print>>f,wo(tz)+','+tr(wg)
def tA():
 raise wM
p=wH(description="Histogram luminance of a JPG image")
p.add_argument("text_file",nargs=1,help="Image path list file")
p.add_argument("-x","--threshold",dest="thresh",wg=ws,default=wT,help="Threshold the pixel map")
r=p.add_argument_group("Blob Identification")
r.add_argument("-c","--contours",dest="contours",wg=ws,default=wT,help="Identify blob contours above some threshold")
r.add_argument("-d","--distance",dest="distance",wg=ws,default=75.,help="Group blobs within some distance of each other")
r.add_argument("-a","--min-area",dest="area",wg=ws,default=2.,help="Remove blobs below some minimum area")
p=p.parse_args()
C=('jpg jpeg bmp png eps gif im j2k j2p jpx msp pcx png ppm pgm pbm'+'spi tiff webp xbm xv').split(' ')
X=wk(p.text_file[0],'r+')
tw=[]
for tH in X:
 ta=tH.strip()
 if ta.split('.')[-1]in C:
  tw.append(ta)
tw=wR(wl(tw))
tJ=[]
tP=[]
if wQ(tw)==0:
 print('No image files were found.')and tA()
f=wk('classifications.out','w')
xy=wk('xandyCent.out','w')
def tr(x):
 return{1:'spot',2:'worm',3:'ambig',4:'track'}[x]
def tp(old,new):
 return new if new>old else old
for tn in tw:
 tz=tn.split('/')[-1].split('.')[0]
 ty=wV(tn).convert("L")
 V=[]
 ti=ty.load()
 nx=ty.size[0]
 ny=ty.size[1]
 x0,y0,x1,y1=(0,0,nx,ny)
 for y in wx(ny):
  V.append([ti[x,y]for x in wx(nx)])
 V=wm(V,dtype=ws)
 tu=wD(wD(V))/(wQ(V)*wQ(V[0]))
 if tu>4:
  tN(tz,7)
  continue
 M=tv(V,threshold=args.contours,minArea=args.area)
 D=tO(M,maxDist=args.distance)
 for g in D:
  tJ.append(tW(g))
  tP.append(tG(g))
 if p.thresh!=wT:
  V[V<p.thresh]=0.
 if p.contours==40:
  wg=0
  for i,bg in wY(D):
   X0,X1,Y0,Y1=bg.getSquareBoundingBox()
   l1,l2,r,o=bg.getPrincipalMoments(V)
   tF=(wn(l1**2-l2**2)/l1)
   tm=tJ[i][0]
   tV=tJ[i][1]
   _i=tC(tL((ws(tJ[i][0][1])-tJ[i][1][1]),(ws(tJ[i][0][0])-tJ[i][1][0])),tm)
   te=0.
   ts=tk(tm,tV)
   for pt in tP[i]:
    te+=tf(_i,pt)
   tQ=te/ts
   tj=0.
   tS=[0,0]
   for b in bg.blobs:
    tj+=b.perimeter
    tS[0]+=b.xc*b.perimeter
    tS[1]+=b.yc*b.perimeter
   tU=tS[0]/ws(bg.count*tj)
   tT=tS[1]/ws(bg.count*tj)
   td=tU-bg.xmin
   tK=bg.xmax-tU
   tB=tT-bg.ymin
   tE=bg.ymax-tT
   tR=wS([td,tK])/wU([td,tK])
   tl=wS([tB,tE])/wU([tB,tE])
   tb=tR*tl
   tg=wb(tR-tl)
   if tF>.99993 and l1>700:
    wg=tp(wg,3)
   elif bg.b_area<4 or ts<6 or ts<13 and(r>=.2 and tF<.945 or bg.b_area<7)or tF<.7:
    if bg.b_area>50:
     if bg.b_area>62 and ts>10:
      wg=tp(wg,2)
     else:
      wg=tp(wg,1)
    else:
     wg=tp(wg,1)
   else:
    if tg>.55:
     wg=tp(wg,2)
    elif bg.b_area>100 and(l1>100 and l1/10>l2)and ts>30:
     if(tQ>9)or(l1/5>ts and tQ>3.9)or(80>ts>40 and bg.b_area>100 and tQ>5):
      wg=tp(wg,2)
     else:
      wg=4
    elif tF>.9998 and 130>ts>90:
     wg=4
    elif tF>.9995 and ts>40 and tb>.8:
     wg=tp(wg,3)
    else:
     if(tb>.978 and tg<.01 or tb>.96 and tg<.0047 or tb>.9 and r<.02 and tg<.055 and tQ<5)and tF>.96:
      if bg.b_area>33:
       wg=4
      else:
       wg=tp(wg,2)
     elif(r<.22 and(tF<.985 and(tg<.015 or tb>.88)and 18>ts>9 or .97<tF<.985 and 18>ts>9 and tb>.83 and bg.b_area<30 or .99>tF>.975 and tQ<3.7 and 18>ts>7.6 and tg<.023)or tF>.99 and l1<15 and tb>.86 and tg<.1 and 35>bg.b_area>28):
      wg=4
     else:
      if(tQ>4.6 and bg.b_area<100)or bg.b_area<24 or tF<.979 or(r>.2 and tQ>4):
       wg=tp(wg,2)
      else:
       if tg<.7 and tQ>6:
        wg=tp(wg,2)
       elif(l1>100 and l2<12 and bg.b_area>40 and bg.b_area<60 and tQ<3.8 or(wb(tR)>.99 or wb(tl)>.99)and wb(tR)>.93 and wb(tl)>.93 and tF>.99 and tQ<2.95 and tQ>2. and tg>.05 or(tb>.9 and tg<.02 and tQ<3.1 and tF>.993 and ts>12)and bg.b_area<82 or((tb>.6 and tF>.9923 and tQ<3.1 or tb>.88)and(tg<.03 or wb(tR)>.996 or wb(tl)>.996)and bg.b_area<100)):
        wg=4
       else:
        if tF>.999 and tQ<3.14 and bg.b_area>58 and l1/25>l2:
         wg=4
        elif tF>.999 and .86<tb<.92 and ts>23:
         wg=tp(wg,2)
        elif(tF>.999 and((l1>90 and l2<10)and(tQ>2.9 or tQ<1.1)or tF>.992 and bg.b_area<50 and .98>wb(tR)>.96 and wb(tl)>.96)):
         wg=tp(wg,3)
        elif tb>.75 and tg<.182 and((bg.b_area>28)or(bg.b_area<28 and ts>17)):
         if((tF>.9996 or r<.028)and tb<.9 and ts<30 and bg.b_area<62 or tF>.99 and l1>400 and l1<600 and l2>60 and tQ>3.4 or tF<.99 and tF>.975 and ts<17 and l1<16 and l2>2 and r>.2 or tF>.993 and(tQ<3 and 28<ts<40 and .94>tb>.9 and bg.b_area<50 or 3.5<tQ<4 and 17<ts<25 and r<.12)):
          wg=tp(wg,2)
         elif(tQ<3.76 and tF>.99 and tg<.06 and r<.13 and(bg.b_area>60 or ts>10)and wU(wb(tR),wb(tl))>.935 and wS(wb(tR),wb(tl))>.875 or tQ<4.1 and bg.b_area>30 and tg<0.059 and ts<16 or(tQ<4.16 and tb>.74 and tg<.012 and bg.b_area<50 and ts<20 and 12<l1<23 and l2<3)):
          wg=4
         else:
          wg=tp(wg,2)
        elif tb>.75 and tg<.05 and bg.b_area>30:
         wg=4
        elif .45<tb<.6 and .2<tg<.5 and .999>tF>.92:
         wg=4
        else:
         wg=tp(wg,2)
   if wg==4:
    break
  tN(tz,wg)
 else:
  print('Only image analysis at contour level 40 is supported, sorry.')and tA()
# Created by pyminifier (https://github.com/liftoff/pyminifier)
