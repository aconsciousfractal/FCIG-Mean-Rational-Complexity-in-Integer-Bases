// Exhaustive consumer. No lattice reduction, witness search, or minimality claim.
#include "kernel.hpp"
#include <cstdlib>
I signed32(const unsigned char*p){uint32_t z=0;for(int j=0;j<4;j++)z|=uint32_t(p[j])<<(8*j);return z<=INT32_MAX?I(z):-1-I(~z);}
int main(int argc,char**argv){try{
 need(argc==7,"r s p exponent data prefix");init(stoi(argv[1]),stoi(argv[2]),stoi(argv[3]),stoi(argv[4]),string(argv[5])+"/log_table.txt");load(argv[5]);string pre=argv[6];
 const char* t=getenv("CSTAR_THREADS");int nt=t?stoi(t):2;need(nt>=1&&nt<=8,"thread count");omp_set_num_threads(nt);
 streamoff size=streamoff(cells.size())*NB*8;
 {ifstream in(pre+".witness",ios::binary|ios::ate);need(bool(in),"witness file exists");need(in.tellg()>=size,"truncated witness records");need(in.tellg()==size,"extra witness records");}
 vector<I>ubs(cells.size());vector<int>labs(cells.size());int failed=0;string failure;
 #pragma omp parallel
 {ifstream f(pre+".witness",ios::binary);vector<unsigned char>raw(size_t(NB)*8);
  #pragma omp for schedule(dynamic,1)
  for(int c=0;c<(int)cells.size();c++){try{
   f.seekg((streamoff)c*NB*8);f.read((char*)raw.data(),raw.size());need(f.gcount()==(streamsize)raw.size(),"truncated witness records");vector<I>x(NL);
   for(int h=0;h<NB;h++){
    I a=signed32(&raw[8*h]),b=signed32(&raw[8*h+4]);need(units[mod(Z(a)*br[h].r[0]+Z(b)*br[h].r[1],Q)],"nonunit full mark");
    I u=safe(Z(a)*br[h].H[0]+Z(b)*br[h].H[2]),v=safe(Z(a)*br[h].H[1]+Z(b)*br[h].H[3]);need(az(u)<1000000&&az(v)<1000000,"parent coefficient range");
    x[br[h].kernel]=safe(Z(x[br[h].kernel])+bound(cells[c],{u,v,a,b}));
   }
   tensor(x);auto it=max_element(x.begin(),x.end());ubs[c]=safe(ce(*it,K));labs[c]=it-x.begin();
  }catch(exception&e){
   #pragma omp critical
   {failed=1;failure=string("cell ")+to_string(c)+": "+e.what();}
  }}
 }
 if(failed)throw runtime_error(failure);
 ofstream out(pre+".checked.cells");need(bool(out),"profile output");
 for(int c=0;c<(int)cells.size();c++)out<<c<<" "<<cells[c].lo<<" "<<cells[c].hi<<" "<<ubs[c]<<" "<<labs[c]<<"\n";
 need(bool(out),"profile written");unsigned long long records=(unsigned long long)cells.size()*NB;
 cout<<"PASS consumer Q="<<Q<<" K="<<K<<" records="<<records<<" corners="<<4*records<<"\n";
}catch(exception&e){cerr<<"FAIL "<<e.what()<<"\n";return 2;}}
