#!/usr/bin/env python
# WIPAC DECO event labeling script. GNU 3.0 Open-sourced
CG=float
Cm=len
CI=range
Cz=min
Cl=max
CU=None
Cw=int
Cy=False
Cu=True
Ca=map
CB=sorted
Ce=set
CK=abs
Cg=type
Ci=str
Ck=SystemExit
CL=open
Cq=xrange
Ct=sum
Cp=enumerate
import argparse,math,os,sys
fT=argparse.ArgumentParser
fM=sys.maxint
fj=math.sqrt
import numpy as np
Cr=np.vstack
fO=np.mean
CE=np.array
CD=np.CK
CV=np.arctan2
Cs=np.Ct
Cf=np.sqrt
CW=np.sign
from PIL import Image
CX=Image.CL
from skimage import measure
Cb=measure.find_contours
class fo:
 def __init__(f,m,p):
  f.m=CG(m)
  f.p=(CG(p[0]),CG(p[1]))
 def fk(f,x):
  return f.m*(x-f.p[0])+f.p[1]
 def fL(f,M):
  if M.m==f.m:
   return NaN
  else:
   C=(M.p[1]-f.p[1]+f.m*f.p[0]+M.m*M.p[0])/(f.m-M.m)
   return(C,f.fk(C))
class fF:
 def __init__(f,x,y):
  f.x=x
  f.y=y
  f.xc=fO(x)
  f.yc=fO(y)
  f.area=0.
  n=Cm(x)
  for i in CI(0,n):
   f.area+=0.5*(y[i]+y[i-1])*(x[i]-x[i-1])
  f.perimeter=0
  D=Cm(f.x)-1
  for i in CI(0,D):
   a=(f.x[i],f.y[i])
   b=(f.y[i+1],f.y[i+1])
   f.perimeter+=fq(a,b)
  f.maxdist=0
  D=Cm(f.x)
  for i in CI(0,D):
   for j in CI(i,D):
    s=(f.x[i],f.y[i])
    r=(f.y[j],f.y[j])
    V=fq(s,r)
    if V>f.maxdist:
     f.maxdist=V
 def fq(f,W):
  return Cf((f.xc-W.xc)**2+(f.yc-W.yc)**2)
def fq(a,b):
 return fj((a[1]-b[1])**2+(a[0]-b[0])**2)
class fd:
 def __init__(f):
  f.blobs=[]
  f.count=0
  f.area=0.
  f.b_area=0.
  f.perimeter=0.
  f.xmin= 1e10
  f.xmax=-1e10
  f.ymin= 1e10
  f.ymax=-1e10
 def ft(f,W):
  f.blobs.append(W)
  f.xmin=Cz(f.xmin,W.x.Cz())
  f.xmax=Cl(f.xmax,W.x.Cl())
  f.ymin=Cz(f.ymin,W.y.Cz())
  f.ymax=Cl(f.ymax,W.y.Cl())
  f.cov =CU
  f.count+=1
  f.b_area+=W.area
  f.perimeter+=W.perimeter
 def fp(f):
  return(f.xmin,f.xmax,f.ymin,f.ymax)
 def fh(f):
  E,X,b,G=(f.xmin,f.xmax,f.ymin,f.ymax)
  xL=CD(X-E)
  yL=CD(G-b)
  if xL>yL:
   b-=0.5*(xL-yL)
   G+=0.5*(xL-yL)
  else:
   E-=0.5*(yL-xL)
   X+=0.5*(yL-xL)
  return(E,X,b,G)
 def fQ(f,z):
  ny,nx=z.shape
  x0,x1,y0,y1=bg.getBoundingBox()
  f.area=(x1-x0)*(y1-y0)
  i0,i1=[ny-Cw(t)for t in(y1,y0)]
  j0,j1=[Cw(t)for t in(x0,x1)]
  m=1
  i0=0 if i0-m<0 else i0-m
  i1=ny-1 if i1>ny-1 else i1+m
  j0=0 if j0-m<0 else j0-m
  j1=nx-1 if j1>nx-1 else j1+m
  return z[i0:i1,j0:j1]
 def fH(f,z,p,q):
  nx,ny=z.shape
  I=0.
  if p==0 and q==0:
   I=Cs(z)
  else:
   for i in CI(0,nx):
    x=0.5+i
    for j in CI(0,ny):
     y=0.5+j
     I+=x**p*y**q*z[i,j]
  return I
 def fR(f,z):
  if f.cov==CU:
   l=f.fQ(z).transpose()
   U=f.fH(l,0,0)
   w=f.fH(l,1,0)
   y=f.fH(l,0,1)
   u=f.fH(l,1,1)
   a=f.fH(l,2,0)
   B=f.fH(l,0,2)
   e=w/U
   K=y/U
   f.cov=Cr([[a/U-e*e,u/U-e*K],[u/U-e*K,B/U-K*K]])
  return f.cov
 def fJ(f,z):
  g=f.fR(z)
  i=g[0,0]
  k=g[0,1]
  L=g[1,1]
  q=0.5*CV(2*k,i-L)
  l1=0.5*(i+L)+0.5*Cf(4*k**2+(i-L)**2)
  l2=0.5*(i+L)-0.5*Cf(4*k**2+(i-L)**2)
  return l1,l2,l1/l2,q
def fY(z,threshold,minArea=2.):
 t=[]
 ny,nx=z.shape
 p=Cb(z,threshold)
 for h in p:
  x=h[:,1]
  y=ny-h[:,0]
  W=fF(x,y)
  if W.area>=minArea:
   t.append(W)
 return t
def fS(t,maxDist):
 n=Cm(t)
 Q=[]
 if n>=1:
  bg=fd()
  bg.addBlob(t[0])
  Q.append(bg)
  for i in CI(1,n):
   bi=t[i]
   H=Cy
   for R in Q:
    for bj in R.blobs:
     if bi.distance(bj)<maxDist:
      R.ft(bi)
      H=Cu
      break
   if not H:
    bg=fd()
    bg.addBlob(bi)
    Q.append(bg)
 return Q
def fn(R):
 J=[]
 Y=[0]
 S=[0,0,0]
 n=100.
 N=Cu
 for b in R.blobs:
  for i in CI(0,Cm(b.x)):
   J.append([CG(b.x[i]),CG(b.y[i])])
 for i in CI(0,Cm(J)):
  for j in CI(i,Cm(J)):
   if(fq(J[i],J[j])>Y[0]):
    Y=[fq(J[i],J[j]),J[i],J[j]]
 if Y[1][0]<Y[2][0]:
  A=Y[1][0]
  y1=Y[1][1]
  c=Y[2][0]
 else:
  A=Y[2][0]
  y1=Y[2][1]
  c=Y[1][0]
  N=Cy
 x=Y[1][1]-Y[2][1]
 v=Y[1][0]-Y[2][0]
 if N:
  return[Ca(Cw,Y[1]),Ca(Cw,Y[2])]
 else:
  return[Ca(Cw,Y[2]),Ca(Cw,Y[1])]
def fN(R):
 J=[]
 for b in R.blobs:
  for i in CI(0,Cm(b.x)):
   J.append((CG(b.x[i]),CG(b.y[i])))
 return CB(Ce(J))
def fA(x,v):
 if v==0 and x==0:
  return NaN
 elif v==0:
  return CW(x)*Cf(fM)
 else:
  return CG(x)/v
def fc(M,point):
 return CK(-M.m*point[0]+point[1]-(M.p[1]-M.m*M.p[0]))/CG(Cf((-M.m)**2+1))
def fx(fr,Cg):
 print>>f,Ci(fr)+','+fP(Cg,0)
def fv():
 raise Ck
p=fT(description="Histogram luminance of a JPG image")
p.add_argument("text_file",nargs=1,help="Image filepath names")
p.add_argument("-x","--threshold",dest="thresh",Cg=CG,default=CU,help="Threshold the pixel map")
F=p.add_argument_group("Blob Identification")
F.add_argument("-c","--contours",dest="contours",Cg=CG,default=CU,help="Identify blob contours above some threshold")
F.add_argument("-d","--distance",dest="distance",Cg=CG,default=75.,help="Group blobs within some distance of each other")
F.add_argument("-a","--min-area",dest="area",Cg=CG,default=2.,help="Remove blobs below some minimum area")
d=p.parse_args()
T=CL(d.text_file[0],'r+')
j=[]
for M in T:
 j.append(M.strip('\n'))
O=[]
fC=[]
if Cm(j)==0:
 print('No image files were found.')and fv()
f=CL('classifications.out','w')
def fP(x,opt):
 if opt==0:
  return{1:'spot',2:'worm',3:'track',4:'ambig',5:'big_spot',6:'track_lowconf',7:'noise'}[x]
 else:
  return{1:'_s',2:'_w',3:'_t',4:'_a',5:'_b',6:'_l',7:'_n'}[x]
for fD in j:
 fs=CX(fD).convert("L")
 fr=fD.split('/')[-1].split('.')[0]
 z=[]
 fV=fs.load()
 nx=fs.size[0]
 ny=fs.size[1]
 x0,y0,x1,y1=(0,0,nx,ny)
 for y in Cq(ny):
  z.append([fV[x,y]for x in Cq(nx)])
 z=CE(z,dtype=CG)
 fW=Ct(Ct(z))/(Cm(z)*Cm(z[0]))
 if fW>4:
  fx(fr,7)
  continue
 t=fY(z,threshold=args.contours,minArea=args.area)
 Q=fS(t,maxDist=args.distance)
 for g in Q:
  O.append(fn(g))
  fC.append(fN(g))
 if d.thresh!=CU:
  z[z<d.thresh]=0.
 if d.contours!=CU:
  for i,bg in Cp(Q):
   X0,X1,Y0,Y1=bg.getSquareBoundingBox()
   l1,l2,r,q=bg.getPrincipalMoments(z)
   fE=(Cf(l1**2-l2**2)/l1)
   fX=O[i][0]
   fb=O[i][1]
   _i=fo(fA((CG(O[i][0][1])-O[i][1][1]),(CG(O[i][0][0])-O[i][1][0])),fX)
   fG=0.
   fm=fq(fX,fb)
   for pt in fC[i]:
    fG+=fc(_i,pt)
   fI=fG/fm 
   fz=0.
   fl=[0,0]
   for b in bg.blobs:
    fz+=b.perimeter
    fl[0]+=b.xc*b.perimeter
    fl[1]+=b.yc*b.perimeter
   fU=fl[0]/CG(bg.count*fz)
   fw=fl[1]/CG(bg.count*fz)
   fy=fU-bg.xmin
   fu=bg.xmax-fU
   fa=fw-bg.ymin
   fB=bg.ymax-fw
   fe=Cz([fy,fu])/Cl([fy,fu])
   fK=Cz([fa,fB])/Cl([fa,fB])
   fg=fe*fK
   fi=CK(fe-fK)
   if d.contours==40 or d.contours==20:
    if fE>.99993 and l1>700:
     Cg=4
    elif bg.b_area<4 or fm<6 or fm<13 and(r>=.2 and fE<.945 or bg.b_area<7)or fE<.7:
     if bg.b_area>50:
      if bg.b_area>62 and fm>10:
       Cg=2
      else:
       Cg=5
     else:
      Cg=1
    else:
     if fi>.55:
      Cg=2
     elif bg.b_area>100 and(l1>100 and l1/10>l2)and fm>30:
      if(fI>9)or(l1/5>fm and fI>3.9)or(80>fm>40 and bg.b_area>100 and fI>5):
       Cg=2
      else:
       Cg=3
     elif fE>.9998 and 130>fm>90:
      Cg=3
     elif fE>.9995 and fm>40 and fg>.8:
      Cg=4
     else:
      if(fg>.978 and fi<.01 or fg>.96 and fi<.0047 or fg>.9 and r<.02 and fi<.055 and fI<5)and fE>.96:
       if bg.b_area>33:
        Cg=3
       else:
        Cg=2
      elif(r<.22 and(fE<.985 and(fi<.015 or fg>.88)and 18>fm>9 or .97<fE<.985 and 18>fm>9 and fg>.83 and bg.b_area<30 or .99>fE>.975 and fI<3.7 and 18>fm>7.6 and fi<.023)or fE>.99 and l1<15 and fg>.86 and fi<.1 and 35>bg.b_area>28):
       Cg=3
      else:
       if(fI>4.6 and bg.b_area<100)or bg.b_area<24 or fE<.979 or(r>.2 and fI>4):
        Cg=2
       else:
        if fi<.7 and fI>6:
         Cg=2
        elif(l1>100 and l2<12 and bg.b_area>40 and bg.b_area<60 and fI<3.8 or(CK(fe)>.99 or CK(fK)>.99)and CK(fe)>.93 and CK(fK)>.93 and fE>.99 and fI<2.95 and fI>2. and fi>.05 or(fg>.9 and fi<.02 and fI<3.1 and fE>.993 and fm>12)and bg.b_area<82 or((fg>.6 and fE>.9923 and fI<3.1 or fg>.88)and(fi<.03 or CK(fe)>.996 or CK(fK)>.996)and bg.b_area<100)):
         Cg=3
        else:
         if fE>.999 and fI<3.14 and bg.b_area>58 and l1/25>l2:
          Cg=3
         elif fE>.999 and .86<fg<.92 and fm>23:
          Cg=2
         elif(fE>.999 and((l1>90 and l2<10)and(fI>2.9 or fI<1.1)or fE>.992 and bg.b_area<50 and .98>CK(fe)>.96 and CK(fK)>.96)):
          Cg=4
         elif fg>.75 and fi<.182 and((bg.b_area>28)or(bg.b_area<28 and fm>17)):
          if((fE>.9996 or r<.028)and fg<.9 and fm<30 and bg.b_area<62 or fE>.99 and l1>400 and l1<600 and l2>60 and fI>3.4 or fE<.99 and fE>.975 and fm<17 and l1<16 and l2>2 and r>.2 or fE>.993 and(fI<3 and 28<fm<40 and .94>fg>.9 and bg.b_area<50 or 3.5<fI<4 and 17<fm<25 and r<.12)):
           Cg=2
          elif(fI<3.76 and fE>.99 and fi<.06 and r<.13 and(bg.b_area>60 or fm>10)and Cl(CK(fe),CK(fK))>.935 and Cz(CK(fe),CK(fK))>.875 or fI<4.1 and bg.b_area>30 and fi<0.059 and fm<16 or(fI<4.16 and fg>.74 and fi<.012 and bg.b_area<50 and fm<20 and 12<l1<23 and l2<3)):
           Cg=3
          else:
           Cg=2
         elif fg>.75 and fi<.05 and bg.b_area>30:
          Cg=6
         elif .45<fg<.6 and .2<fi<.5 and .999>fE>.92:
          Cg=3
         else:
          Cg=2
    fx(fr,Cg)
   else:
    print('As of now only image analysis at contour level 40 is supported, sorry.')and fv()
# Created by pyminifier (https://github.com/liftoff/pyminifier)

