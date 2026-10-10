// Consumer independent of Gauss, oracle optimisation and direction models.
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <string>
int64_t readword(std::istream& f){uint64_t x=0;for(int j=0;j<4;j++){int c=f.get();if(c==EOF)throw std::runtime_error("truncated witness");x|=uint64_t(c)<<(8*j);}return x<(1ULL<<31)?int64_t(x):int64_t(x)-(1LL<<32);}
int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("args domain witness-dir histogram.csv");
 std::ifstream domain(argv[1]);std::ofstream out(argv[3]);if(!domain||!out)throw std::runtime_error("missing file");
 out<<"M,height,count\n";int64_t M,previous=0;uint64_t mass=0,stored=0;int moduli=0;
 while(domain>>M){if(M<=previous||M%210||M>=150000)throw std::runtime_error("domain ordering/range");previous=M;
  std::ifstream f(std::string(argv[2])+"/w"+std::to_string(M)+".bin",std::ios::binary);if(!f)throw std::runtime_error("missing witness file");
  std::vector<int64_t> hist(M/2+1);int64_t sum=0;
  for(int64_t s=0;s<=M/2;s++){
   int64_t a=readword(f),d=readword(f);if(d<1||d>M/2||a< -M/2||a>M/2||std::gcd(d,M)!=1||(a-s*d)%M)throw std::runtime_error("invalid witness: congruence/unit/range");
   int64_t h=std::max(a<0?-a:a,d),w=(s==0||2*s==M)?1:2;if(h<1||h>M/2)throw std::runtime_error("height range");
   hist[h]+=w;sum+=w;stored++;
  }
  if(f.get()!=EOF)throw std::runtime_error("trailing witness bytes");if(sum!=M)throw std::runtime_error("mass mismatch");
  for(int h=1;h<=M/2;h++)if(hist[h])out<<M<<","<<h<<","<<hist[h]<<"\n";
  mass+=sum;moduli++;
 }
 if(!out)throw std::runtime_error("histogram write failure");
 std::cout<<"PASS consumer moduli "<<moduli<<" stored witnesses "<<stored<<" full mass "<<mass<<"\n";
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 2;}}
