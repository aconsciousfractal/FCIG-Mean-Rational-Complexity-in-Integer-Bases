#include "kernel.hpp"
I get32(const unsigned char*p){uint32_t z=0;for(int j=0;j<4;j++)z|=uint32_t(p[j])<<(8*j);return z<=INT32_MAX?I(z):-1-I(~z);}
int main(int argc,char**argv){try{need(argc==7,"r s p exponent data prefix");init(stoi(argv[1]),stoi(argv[2]),stoi(argv[3]),stoi(argv[4]),string(argv[5])+"/log_table.txt");load(argv[5]);string pre=argv[6];ifstream f(pre+".witness",ios::binary);ofstream o(pre+".audit");need(bool(f)&&bool(o),"audit files");
 for(int c:{0,1,20,21,137,512,1024,1536,2048,2667,2686,2687})for(int h:{0,1,NB/2,NB-2,NB-1}){
  unsigned char rec[8];f.seekg(((streamoff)c*NB+h)*8);f.read((char*)rec,8);need(f.gcount()==8,"audit record exists");I a=get32(rec),b=get32(rec+4);I u=safe(Z(a)*br[h].H[0]+Z(b)*br[h].H[2]),v=safe(Z(a)*br[h].H[1]+Z(b)*br[h].H[3]);I z=bound(cells[c],{u,v,a,b});o<<c<<" "<<h<<" "<<a<<" "<<b<<" "<<u<<" "<<v<<" "<<z<<"\n";
 }need(bool(o),"audit written");cout<<"PASS audit Q="<<Q<<" K="<<K<<" rows=60\n";
}catch(exception&e){cerr<<"FAIL "<<e.what()<<"\n";return 2;}}
