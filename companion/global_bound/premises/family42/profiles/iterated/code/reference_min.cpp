// Independent definition-based minima: no lattice reduction, profiles, or geometric search.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("M q output required");
 const int M=std::stoi(argv[1]),q=std::stoi(argv[2]);
 if(M<2||M>1000000||q<2||q>M)throw std::runtime_error("reference bank domain");
 std::vector<unsigned char>unit(q);for(int i=0;i<q;i++)unit[i]=std::gcd(i,q)==1;
 std::vector<std::uint32_t>R(M);std::uint64_t tested=0;
 for(int s=0;s<=M/2;s++){
  int best=std::max(1,s);
  for(int d=2;d<best;d++)if(unit[d%q]){
   int a=(std::int64_t(s)*d)%M;a=std::min(a,M-a);
   best=std::min(best,std::max(a,d));++tested;
  }
  R[s]=best;if(s)R[M-s]=best;
 }
 std::ofstream f(argv[3],std::ios::binary);if(!f)throw std::runtime_error("output");
 for(auto x:R)for(int j=0;j<4;j++)f.put(char(x>>(8*j)));
 if(!f)throw std::runtime_error("write");
 std::cout<<"PASS reference M="<<M<<" q="<<q<<" states="<<M<<" representatives="<<M/2+1<<" unit_denominators="<<tested<<"\n";
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 2;}}
