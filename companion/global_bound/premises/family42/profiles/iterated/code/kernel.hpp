/* Integer kernel for selected mark-30 and mark-30p directional envelopes.
 * Row-basis and mark changes are simultaneous; ambient orientation is retained.
 * No numerical radical-210 result is imported. See proofs/DOSSIER.md.
 */
#include <array>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <omp.h>
#include <cstdlib>
using namespace std;
using I=long long;
using Z=__int128_t;
using UZ=__uint128_t;
using Mat=array<I,4>;
using Mark=array<int,2>;
#ifndef PRIME
#define PRIME 5
#endif
#ifndef FOUR_PRIMES
#define FOUR_PRIMES 1
#endif
constexpr int P=PRIME;
#if FOUR_PRIMES
constexpr int Q=42*P, NL=96*(P+1), NPR=4;
#else
constexpr int Q=6*P, NL=12*(P+1), NPR=3;
#endif
constexpr bool prime_test(int n){for(int d=2;d*d<=n;d++)if(n%d==0)return false;return n>=2;}
static_assert((FOUR_PRIMES ? (P>=11 && P<=31) : (P==5)) && prime_test(P),"declared prime domain");
const I SCALE=1000000000;
int ps[]={
#if FOUR_PRIMES
  2,3,7,P
#else
  2,3,P
#endif
}
;
int units[Q];
vector<I> logs;
I PDEN;
array<I,164> radius_bounds;
array<I,21> angle_bounds;
void need(bool b,const string&m){
  if(!b)throw runtime_error(m);
}
I safe(Z x){
  if(x<INT64_MIN||x>INT64_MAX)throw runtime_error("int64 range");
  return (I)x;
}
string text(Z x){
  if(!x)return "0";
  bool neg=x<0;
  if(neg)x=-x;
  string s;
  while(x){
    s.push_back('0'+x%10);
    x/=10;
  }
  if(neg)s+='-';
  reverse(s.begin(),s.end());
  return s;
}
int mod(I x,int p){
  x%=p;
  return x<0?x+p:x;
}
I ab(I x){
  return x<0?-x:x;
}
Z abZ(Z x){
  return x<0?-x:x;
}
I idivfloor(Z a,Z b){
  need(b>0,"divisor sign");
  Z q=a/b,r=a%b;
  if(r<0)q--;
  need(q>=INT64_MIN&&q<=INT64_MAX,"quotient overflow");
  return (I)q;
}
I ceildivZ(Z a,Z b){
  return -idivfloor(-a,b);
}
int invp(int x,int p){
  for(int i=1;i<p;i++)if(x*i%p==1)return i;
  throw runtime_error("nonunit inverse");
}
Z det(Mat B){
  return Z(B[0])*B[3]-Z(B[1])*B[2];
}
int label(Mark r){
  int ret=0;
  for(int p:ps){
    int a=mod(r[0],p),b=mod(r[1],p);
    need(a||b,"nonprimitive mark");
    ret=ret*(p+1)+(a?b*invp(a,p)%p:p);
  }
  return ret;
}
Mark crtmark(array<int,NPR> aa,array<int,NPR> bb){
  Mark r={
    0,0
  }
  ;
  for(int k=0;k<NPR;k++){
    int p=ps[k],w=Q/p*invp((Q/p)%p,p);
    r[0]+=aa[k]*w;
    r[1]+=bb[k]*w;
  }
  r[0]%=Q;
  r[1]%=Q;
  return r;
}
struct Branch{
  Mat H;
  Mark r;
  int kernel;
}
;
array<Branch,NL> branches;
array<array<int,Q>,NL> admitted;
Mat mm(Mat H,Mat G){
  Mat out={
    safe(Z(H[0])*G[0]+Z(H[1])*G[2]),safe(Z(H[0])*G[1]+Z(H[1])*G[3]),safe(Z(H[2])*G[0]+Z(H[3])*G[2]),safe(Z(H[2])*G[1]+Z(H[3])*G[3])
  }
  ;
  return out;
}
void init(const string&f){
  ifstream in(f);
  need(bool(in),"log table missing");
  I v;
  while(in>>v)logs.push_back(v);
  need(logs.size()==12545,"log table size");
  for(int i=0;i<Q;i++)units[i]=gcd(i,Q)==1;
  int cnt=0;
  array<int,NL> used{
  }
  ;
  for(int a=1;a<=Q;a++)if(Q%a==0)for(int b=0;b<a;b++){
    int d=Q/a;
    array<int,NPR> aa,bb,ka,kb;
    for(int k=0;k<NPR;k++){
      int p=ps[k];
      if(a%p ||b%p){
        aa[k]=a%p;
        bb[k]=b%p;
      }
      else {
        aa[k]=0;
        bb[k]=d%p;
      }
      if(a%p){
        ka[k]=0;
        kb[k]=1;
      }
      else if(b%p){
        ka[k]=mod(-d,p);
        kb[k]=b%p;
      }
      else {
        ka[k]=1;
        kb[k]=0;
      }
    }
    int k=label(crtmark(ka,kb));
    branches[k]={
      {
        a,0,b,d
      }
      ,crtmark(aa,bb),k
    }
    ;
    used[k]++;
    cnt++;
  }
  need(cnt==NL,"branch count");
  for(auto u:used)need(u==1,"kernel bijection");
  for(int i=0;i<NL;i++){
    array<int,NPR> a;
    int z=i;
    for(int k=NPR-1;k>=0;k--){
      a[k]=z%(ps[k]+1);
      z/=ps[k]+1;
    }
    int n=0;
    for(int j=0;j<NL;j++){
      int y=j;
      bool ok=true;
      for(int k=NPR-1;k>=0;k--){
        int b=y%(ps[k]+1);
        y/=ps[k]+1;
        if(a[k]==b)ok=false;
      }
      if(ok)admitted[i][n++]=j;
    }
    need(n==Q,"admitted count");
  }
}
void gauss(Mat &B,Mark&r,Mat*V=nullptr){
  Z old=abZ(det(B));
  if(V)*V={
    1,0,0,1
  }
  ;
  for(int it=0;it<200;it++){
    Z nu=Z(B[0])*B[0]+Z(B[1])*B[1],nv=Z(B[2])*B[2]+Z(B[3])*B[3];
    if(nu>nv){
      swap(B[0],B[2]);
      swap(B[1],B[3]);
      swap(r[0],r[1]);
      if(V){
        swap((*V)[0],(*V)[2]);
        swap((*V)[1],(*V)[3]);
      }
      nu=nv;
    }
    need(nu>0,"zero shortest");
    I t=idivfloor(2*(Z(B[0])*B[2]+Z(B[1])*B[3])+nu,2*nu);
    if(!t){
      need(abZ(det(B))==old,"Gauss det");
      return;
    }
    B[2]=safe(Z(B[2])-Z(t)*B[0]);
    B[3]=safe(Z(B[3])-Z(t)*B[1]);
    r[1]=mod(r[1]-t*r[0],Q);
    if(V){
      (*V)[2]=safe(Z((*V)[2])-Z(t)*(*V)[0]);
      (*V)[3]=safe(Z((*V)[3])-Z(t)*(*V)[1]);
    }
  }
  throw runtime_error("Gauss did not terminate");
}
int bitsZ(Z n){
  need(n>0,"bits positive");
  unsigned long long hi=(UZ)n>>64;
  if(hi)return 128-__builtin_clzll(hi);
  return 64-__builtin_clzll((unsigned long long)n);
}
void tensor(array<I,NL>&x){
  int stride=1;
  for(int kk=NPR-1;kk>=0;kk--){
    int n=ps[kk]+1;
    for(int start=0;start<NL;start+=stride*n)for(int off=0;off<stride;off++){
      I sum=0;
      for(int t=0;t<n;t++)sum+=x[start+off+t*stride];
      for(int t=0;t<n;t++)x[start+off+t*stride]=sum-x[start+off+t*stride];
    }
    stride*=n;
  }
}
