
#include <cstdlib>
void gauss(Mat& B,Mark& r,Mat& V){V={1,0,0,1};Z od=az(det(B));for(int k=0;k<200;k++){
 Z a=Z(B[0])*B[0]+Z(B[1])*B[1],b=Z(B[2])*B[2]+Z(B[3])*B[3];
 if(a>b){swap(B[0],B[2]);swap(B[1],B[3]);swap(r[0],r[1]);swap(V[0],V[2]);swap(V[1],V[3]);a=b;}
 need(a>0,"nonzero basis row");I t=safe(fl(2*(Z(B[0])*B[2]+Z(B[1])*B[3])+a,2*a));
 if(!t){need(az(det(B))==od&&az(det(V))==1,"Gauss determinant");return;}
 B[2]=safe(Z(B[2])-Z(t)*B[0]);B[3]=safe(Z(B[3])-Z(t)*B[1]);r[1]=mod(Z(r[1])-Z(t)*r[0],Q);V[2]=safe(Z(V[2])-Z(t)*V[0]);V[3]=safe(Z(V[3])-Z(t)*V[1]);}
 throw runtime_error("Gauss termination");}
W witness(int c,int h){Mat Y=mm(br[h].H,cells[c].C),V;Mark m=br[h].r;gauss(Y,m,V);Z best=Z(1)<<120;I ii=0,jj=0;
 auto test=[&](I i,I j){if(!units[mod(Z(i)*m[0]+Z(j)*m[1],Q)])return;Z a=Z(i)*Y[0]+Z(j)*Y[2],b=Z(i)*Y[1]+Z(j)*Y[3];need(az(a)<(Z(1)<<59)&&az(b)<(Z(1)<<59),"witness size");Z val=a*a+b*b;if(val<best){best=val;ii=i;jj=j;}};
 test(1,0);Z nu=Z(Y[0])*Y[0]+Z(Y[1])*Y[1],dot=Z(Y[0])*Y[2]+Z(Y[1])*Y[3];
 for(I j=1;j<=6;j++){I k=safe(fl(-j*dot,nu));for(I i=k-12;i<=k+12;i++)test(i,j);}
 need(best<(Z(1)<<120),"admissible witness found");I ci=safe(Z(ii)*V[0]+Z(jj)*V[2]),cj=safe(Z(ii)*V[1]+Z(jj)*V[3]);need(units[mod(Z(ci)*br[h].r[0]+Z(cj)*br[h].r[1],Q)],"original child witness mark");
 I u=safe(Z(ci)*br[h].H[0]+Z(cj)*br[h].H[2]),v=safe(Z(ci)*br[h].H[1]+Z(cj)*br[h].H[3]);need(az(u)<1000000&&az(v)<1000000,"parent coefficients within arithmetic domain");return {u,v,ci,cj};}
