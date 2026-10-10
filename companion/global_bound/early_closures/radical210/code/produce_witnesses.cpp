#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using I=long long;
using W=__int128_t;
struct V { I a,d; };
W norm(V x) {return W(x.a)*x.a+W(x.d)*x.d;}
W dot(V x,V y) {return W(x.a)*y.a+W(x.d)*y.d;}
I fl(W x,I y) { if(y<=0) throw std::runtime_error("nonpositive divisor"); W z=x/y; if(x<0 && x%y) --z; return I(z); }
I ce(W x,I y) {return -fl(-x,y);}
I absI(I x) {return x<0?-x:x;}
struct Stats {uint64_t rows=0,probes=0,gauss_steps=0;};

V lattice_witness(I M,I s,Stats& stats) {
  I B=std::max<I>(1,std::min(s,M-s));
  V chosen{s<=M/2?s:s-M,1};
  if(B<=1) return chosen;
  V u{M,0},v{s,1};
  for(;;) {
    if(norm(v)<norm(u)) std::swap(u,v);
    W n=norm(u),z=dot(u,v), az=z<0?-z:z;
    if(2*az<=n) break;
    W qr=(az+n/2)/n; if(z<0) qr=-qr;
    v={I(W(v.a)-qr*u.a),I(W(v.d)-qr*u.d)};
    ++stats.gauss_steps;
  }
  if(u.d<0) {u.a=-u.a;u.d=-u.d;}
  if(W(u.a)*v.d-W(u.d)*v.a<0) {v.a=-v.a;v.d=-v.d;}
  if(W(u.a)*v.d-W(u.d)*v.a!=M || u.d<=0) throw std::runtime_error("basis invariant");
  I g=std::gcd(u.d,M);
  auto row=[&](I k) {
    ++stats.rows;
    if(std::gcd(absI(k),g)!=1) return;
    W ad=W(v.a)*k, dd=W(v.d)*k;
    I tmin=ce(1-dd,u.d);
    I sg=(u.a>0)-(u.a<0);
    I t0=std::max(fl(1-dd,u.d),fl(-W(sg)*ad-dd,absI(u.a)+u.d));
    // On either side of this continuous minimizer the objective is monotone.
    for(int dir:{-1,1}) {
      I t=dir<0?t0:std::max(tmin,t0+1);
      for(;t>=tmin;t+=dir) {
        W aw=W(u.a)*t+ad,dw=W(u.d)*t+dd;
        if(dw<=0) throw std::runtime_error("positive denominator invariant");
        W ab=aw<0?-aw:aw, h=std::max(ab,dw);
        if(h>=B) break;
        ++stats.probes;
        if(std::gcd(I(dw),M)==1) {B=I(h); chosen={I(aw),I(dw)}; break;}
      }
    }
  };
  row(0);row(-1);row(1);
  I lo=ce(W(std::min<I>(u.a,0)-u.d)*B,M);
  I hi=fl(W(std::max<I>(u.a,0)+u.d)*B,M);
  for(I k=lo;k<=hi;++k) if(k < -1 || k > 1) row(k);
  if(std::max(absI(chosen.a),chosen.d)!=B || chosen.d<=0 || std::gcd(chosen.d,M)!=1 || (W(chosen.a)-W(s)*chosen.d)%M) throw std::runtime_error("invalid generated witness");
  return chosen;
}


void word(std::ostream&f,uint32_t x){for(int j=0;j<4;j++){char c=char((x>>(8*j))&255);f.write(&c,1);}}
int main(int argc,char**argv){try{
 if(argc!=3)throw std::runtime_error("args: moduli.txt output-directory");
 std::ifstream in(argv[1]);if(!in)throw std::runtime_error("missing domain");
 I M;uint64_t total=0;Stats st;
 while(in>>M){if(M<210||M>=150000||M%210)throw std::runtime_error("modulus domain");
  std::ofstream out(std::string(argv[2])+"/w"+std::to_string(M)+".bin",std::ios::binary);
  if(!out)throw std::runtime_error("witness output");
  for(I s=0;s<=M/2;s++){V w=lattice_witness(M,s,st);word(out,uint32_t(w.a));word(out,uint32_t(w.d));}
  if(!out)throw std::runtime_error("witness output incomplete");total+=M;
 }
 std::cout<<"PRODUCED full-state mass "<<total<<" rows "<<st.rows<<" probes "<<st.probes<<" gauss "<<st.gauss_steps<<"\n";
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 2;}}
