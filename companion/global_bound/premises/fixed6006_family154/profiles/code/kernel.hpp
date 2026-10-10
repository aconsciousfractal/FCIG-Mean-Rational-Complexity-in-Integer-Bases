// One-step full-mark upper certificate. All acceptance arithmetic is integral.
// Euclidean witness norm is an upper bound for the original infinity norm.
#include <array>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <omp.h>
using namespace std;
using I=long long;using Z=__int128_t;using UZ=__uint128_t;
using Mat=array<I,4>;using Mark=array<int,2>;
constexpr I S=1000000000, VS=100000000000000LL;
int Q,NL,K,NB;vector<int>ps;vector<int>units;vector<I>logs;
void need(bool b,const char*m){if(!b)throw runtime_error(m);}
I safe(Z z){need(z>=INT64_MIN&&z<=INT64_MAX,"64-bit range");return I(z);}
Z az(Z z){return z<0?-z:z;}
I mod(Z a,int p){I r=I(a%p);return r<0?r+p:r;}
Z fl(Z a,Z b){need(b>0,"positive divisor");Z q=a/b,r=a%b;return q-(r<0);}
Z ce(Z a,Z b){return -fl(-a,b);}
int bits(Z z){need(z>0,"positive bits");auto hi=(UZ)z>>64;return hi?128-__builtin_clzll((unsigned long long)hi):64-__builtin_clzll((unsigned long long)z);}
int inv(int a,int p){for(int i=1;i<p;i++)if(a*i%p==1)return i;throw runtime_error("nonunit inverse");}
int lab(Mark r){int n=0;for(int p:ps){int a=mod(r[0],p),b=mod(r[1],p);need(a||b,"primitive mark");n=n*(p+1)+(a?b*inv(a,p)%p:p);}return n;}
int crt(vector<int>x){Z z=0;for(int i=0;i<int(ps.size());i++)z+=Z(x[i])*(Q/ps[i])*inv((Q/ps[i])%ps[i],ps[i]);return mod(z,Q);}
Z det(Mat A){return Z(A[0])*A[3]-Z(A[1])*A[2];}
Mat mm(Mat A,Mat B){return {safe(Z(A[0])*B[0]+Z(A[1])*B[2]),safe(Z(A[0])*B[1]+Z(A[1])*B[3]),safe(Z(A[2])*B[0]+Z(A[3])*B[2]),safe(Z(A[2])*B[1]+Z(A[3])*B[3])};}
I lup(Z n){need(n>0&&bits(n)+bits(S)<126,"positive log and 128-bit product range");int k=bits(n)-1,sh=max(k-10,0);Z v=Z(1)<<sh;I i=safe(k>=10?n>>sh:n<<(10-k));Z rem=k>=10?n-Z(i)*v:0,base=k>=10?Z(i)*v:1; // rem/base < 1/1024
 return safe(2*Z(logs[i])+2*Z(k)*logs[2048]+ce(rem*S,base));}
I llo(Z n){need(n>0&&bits(n)+bits(S)<126,"positive log and 128-bit product range");int k=bits(n)-1,sh=max(k-10,0);Z v=Z(1)<<sh;I i=safe(k>=10?n>>sh:n<<(10-k));Z rem=k>=10?n-Z(i)*v:0,base=k>=10?Z(i)*v:1;
 I xu=safe(ce(rem*S,base)),xl=safe(fl(rem*S,base));Z tab=i==1024?0:2*Z(logs[i]-1);
 return safe(tab+2*Z(k)*(logs[2048]-1)+xl-ce(Z(xu)*xu,2*S));}
// floor/ceiling multiplication by VS without multiplying an unbounded numerator.
Z ratio(Z n,Z d){need(n>=0&&d>0&&bits(d)+bits(VS)<126,"ratio denominator bound");Z a=n/d,r=n%d;need(a==0||bits(a)+bits(VS)<126,"ratio quotient bound");return a*VS+ce(r*VS,d);}
struct Branch{Mat H;Mark r;int kernel;};vector<Branch> br;
void init(int p,int unused1,int unused2,int multiplier,string path){
 need(p==1,"fixed6006 mark selector");ps={2,3,7,11,13};need(multiplier==1||multiplier==2,"marked6006 index multiplier");
 Q=1;NL=1;for(int v:ps){for(int d=2;d*d<=v;d++)need(v%d,"prime input");Q*=v;NL*=v+1;}
 K=Q*multiplier;units.resize(Q);for(int i=0;i<Q;i++)units[i]=gcd(i,Q)==1;
 ifstream f(path);I v;while(f>>v)logs.push_back(v);need(logs.size()==2049,"log table size");
 vector<int>seen(NL);
 for(int a=1;a<=K;a++)if(K%a==0)for(int b=0;b<a;b++){int d=K/a;vector<int>ia(ps.size()),ib(ps.size()),ka(ps.size()),kb(ps.size());bool primitive=true;
  for(int i=0;i<int(ps.size());i++){int t=ps[i];if(a%t||b%t){ia[i]=a%t;ib[i]=b%t;}else{ia[i]=0;ib[i]=d%t;}
   if(!ia[i]&&!ib[i]){primitive=false;break;}
   if(a%t){ka[i]=0;kb[i]=1;}else if(b%t){ka[i]=mod(-d,t);kb[i]=b%t;}else{ka[i]=1;kb[i]=0;}}
  if(!primitive)continue;
  int ix=lab({crt(ka),crt(kb)});seen[ix]++;br.push_back({{a,0,b,d},{crt(ia),crt(ib)},ix});}
 for(int n:seen)need(n==K/Q,"prime-power kernel multiplicity");
 sort(br.begin(),br.end(),[](const Branch&a,const Branch&b){return a.kernel!=b.kernel?a.kernel<b.kernel:a.H<b.H;});NB=br.size();need(NB==(K/Q)*NL,"primitive branch count");}
void tensor(vector<I>&x){int stride=1;for(int k=int(ps.size())-1;k>=0;k--){int n=ps[k]+1;for(int st=0;st<NL;st+=stride*n)for(int off=0;off<stride;off++){Z sum=0;for(int t=0;t<n;t++)sum+=x[st+off+t*stride];for(int t=0;t<n;t++)x[st+off+t*stride]=safe(sum-x[st+off+t*stride]);}stride*=n;}}
struct Cell{Mat C;I lo,hi;int tl,th,k,t;};vector<Cell>cells;I DEN;
void load(string root){ifstream f(root+"/grid.tsv");f>>DEN;need(bool(f)&&DEN>0&&DEN<(I(1)<<54),"grid scale");int ix,tl,th;I lo,hi,a,b,d;
 while(f>>ix>>lo>>hi>>tl>>th>>a>>b>>d){need(ix==(int)cells.size()&&lo>0&&hi>lo&&hi<2*DEN&&tl>=-20&&th<=20&&tl<=th,"grid row contract");need(a>0&&d>0&&a<1000000000&&az(b)<1000000000&&d<1000000000,"centre arithmetic range");cells.push_back({{a,0,b,d},lo,hi,tl,th,ix/21-125,ix%21-10});}
 need(f.eof()&&cells.size()==2688,"complete grid");}
struct W{I u,v,ci,cj;};
I bound(Cell c,W w){I ans=INT64_MIN;for(I x:{c.lo,c.hi})for(int tt:{c.tl,c.th}){Z lin=Z(40)*w.u+Z(tt)*w.v;Z sq=ratio(Z(x)*lin*lin,Z(1600)*DEN*K)+ratio(Z(DEN)*w.v*w.v,Z(x)*K);need(sq>0,"positive squared norm");I b=safe(ce(Z(lup(sq))-llo(VS),2));ans=max(ans,b);}return ans;}
