/* Witness-only producer/consumer. The consumer has no lattice reduction.
   Every file has precisely floor(M/2)+1 pairs (a,d), in increasing state order.
   C++17, little-endian signed 64-bit integers, bound M<2000000. */
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <sstream>
using I=std::int64_t;using Z=__int128_t;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
I fl(Z a,Z b){need(b>0,"positive divisor");Z q=a/b;if(a%b<0)--q;need(q>=INT64_MIN&&q<=INT64_MAX,"quotient bound");return I(q);}
I absI(I a){need(a!=INT64_MIN,"absolute bound");return a<0?-a:a;}
I mod(I x,I n){x%=n;return x<0?x+n:x;}
void writeI(std::ofstream&f,I x){std::uint64_t u=std::uint64_t(x);for(int j=0;j<8;j++)f.put(char(u>>(8*j)));need(bool(f),"write failure");}
I readI(std::ifstream&f){std::uint64_t u=0;for(int j=0;j<8;j++){int c=f.get();need(c!=EOF,"truncated witness");u|=std::uint64_t((unsigned char)c)<<(8*j);}if(u<=INT64_MAX)return I(u);std::uint64_t v=~u;need(v<=INT64_MAX,"negative range");return -1-I(v);}
int main(int argc,char**argv){try{
 need(argc==4,"arguments produce/consume domain-file work-directory");std::string mode=argv[1],dir=argv[3];need(mode=="produce"||mode=="consume","mode");
 std::ifstream dom(argv[2]);need(bool(dom),"domain");I M,q;std::uint64_t records=0,mass=0;std::map<I,bool> seen;
 std::string line;
 while(std::getline(dom,line)){
  std::istringstream row(line);need(bool(row>>M>>q),"domain parse");std::string extra;need(!(row>>extra),"extra domain fields");
  need(M>=2&&M<2000000,"modulus domain");need(q>=2&&q<=M,"radical domain");need(!seen[M],"duplicate modulus");seen[M]=true;
  I m=M,rad=1;for(I p=2;p*p<=m;p++)if(m%p==0){rad*=p;do{m/=p;}while(m%p==0);}if(m>1)rad*=m;
  need(rad==q,"exact radical support");
  std::string path=dir+"/M"+std::to_string(M)+".bin";
  if(mode=="produce"){
   std::ofstream f(path,std::ios::binary);need(bool(f),"output witness file");std::vector<bool> units(q);for(I j=0;j<q;j++)units[j]=std::gcd(j,q)==1;
   for(I s=0;s<=M/2;s++){
    I ux=M,uy=0,vx=s,vy=1;bool term=false;
    for(int it=0;it<200;it++){
     Z nu=Z(ux)*ux+Z(uy)*uy,nv=Z(vx)*vx+Z(vy)*vy;
     if(nu>nv){std::swap(ux,vx);std::swap(uy,vy);nu=nv;}
     need(nu>0,"short vector zero");I t=fl(2*(Z(ux)*vx+Z(uy)*vy)+nu,2*nu);
     if(!t){term=true;break;}vx-=t*ux;vy-=t*uy;
    }
    need(term,"reduction termination");need((Z(ux)*vy-Z(uy)*vx==M)||(Z(ux)*vy-Z(uy)*vx==-M),"determinant");
    I aa=s,dd=1,best=std::max(absI(aa),dd);
    auto test=[&](I i,I j){I d=i*uy+j*vy;if(!units[mod(d,q)])return;I a=i*ux+j*vx;if(d<0){d=-d;a=-a;}need(d>0,"nonzero unit denominator");I h=std::max(absI(a),d);if(h<best){best=h;aa=a;dd=d;}};
    test(1,0);Z nu=Z(ux)*ux+Z(uy)*uy,dot=Z(ux)*vx+Z(uy)*vy;
    for(I j=1;j<=6;j++){I k=fl(-j*dot,nu);for(I i=k-6;i<=k+6;i++)test(i,j);}
    writeI(f,aa);writeI(f,dd);
   }
  }else{
   std::ifstream f(path,std::ios::binary);need(bool(f),"missing witness");std::map<I,I> hist;I hm=0;
   for(I s=0;s<=M/2;s++){
    I a=readI(f),d=readI(f);need(d>0&&std::gcd(d,M)==1,"full-modulus denominator unit");need((Z(a)-Z(s)*d)%M==0,"state congruence");I h=std::max(absI(a),d);need(h>=1&&h<=M/2,"height bound");I wt=(s==0||2*s==M)?1:2;hist[h]+=wt;hm+=wt;++records;
   }
   need(f.get()==EOF,"extra witness records");need(hm==M,"full mass");mass+=hm;
   std::ofstream out(dir+"/M"+std::to_string(M)+".hist");for(auto kv:hist)out<<kv.first<<" "<<kv.second<<"\n";need(bool(out),"histogram write");
  }
  std::cout<<"PASS "<<mode<<" M="<<M<<" representatives="<<M/2+1<<"\n";
 }
 need(dom.eof(),"domain read");need(!seen.empty(),"empty domain");std::cout<<"COMPLETE "<<mode<<" moduli="<<seen.size()<<" records="<<records<<" mass="<<mass<<"\n";
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}}
