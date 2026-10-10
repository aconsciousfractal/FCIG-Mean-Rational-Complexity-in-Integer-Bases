#include "kernel.hpp"
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
int main(int argc,char**argv){try{need(argc==7,"r s p exponent data output-prefix");init(stoi(argv[1]),stoi(argv[2]),stoi(argv[3]),stoi(argv[4]),string(argv[5])+"/log_table.txt");load(argv[5]);string pre=argv[6];double st=omp_get_wtime();const char* t=getenv("CSTAR_THREADS");int threads=t?stoi(t):2;need(threads>=1&&threads<=8,"thread count");omp_set_num_threads(threads);
 vector<I>ubs(cells.size());vector<int>labs(cells.size());int failed=0;string failure;ofstream w(pre+".witness",ios::binary);need(bool(w),"witness output path");
 #pragma omp parallel for schedule(dynamic,1)
 for(int c=0;c<(int)cells.size();c++){try{vector<I>b(NL);vector<unsigned char> records(size_t(NB)*8);for(int h=0;h<NB;h++){auto ww=witness(c,h);need(ww.ci>=INT32_MIN&&ww.ci<=INT32_MAX&&ww.cj>=INT32_MIN&&ww.cj<=INT32_MAX,"coefficient storage");for(int j=0;j<4;j++){records[8*h+j]=(uint32_t(ww.ci)>>(8*j))&255;records[8*h+4+j]=(uint32_t(ww.cj)>>(8*j))&255;}b[br[h].kernel]+=bound(cells[c],ww);}
 #pragma omp critical
 {w.seekp((streamoff)c*NB*8);w.write((char*)records.data(),NB*8);need(bool(w),"complete witness chunk");}tensor(b);auto it=max_element(b.begin(),b.end());ubs[c]=safe(ce(*it,K));labs[c]=it-b.begin();}catch(exception&e){
 #pragma omp critical
 {failed=1;failure=string("cell ")+to_string(c)+": "+e.what();}}}
 if(failed)throw runtime_error(failure);
 ofstream f(pre+".cells");need(bool(f),"profile output path");
 for(int c=0;c<(int)cells.size();c++)f<<c<<" "<<cells[c].lo<<" "<<cells[c].hi<<" "<<ubs[c]<<" "<<labs[c]<<"\n";
 w.flush();need(bool(f)&&bool(w),"complete outputs");
 cerr<<"PASS Q="<<Q<<" K="<<K<<" NL="<<NL<<" NB="<<NB<<" cells="<<cells.size()<<" witnesses="<<size_t(cells.size())*NB<<" seconds="<<omp_get_wtime()-st<<"\n";
 }catch(exception&e){cerr<<"FAIL "<<e.what()<<"\n";return 2;}}
