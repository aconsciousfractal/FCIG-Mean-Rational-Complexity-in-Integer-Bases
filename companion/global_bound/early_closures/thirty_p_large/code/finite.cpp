// Every witness is checked after it has been written to disk.
// The consumer uses no lattice reduction, search, or minimality assumption.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <chrono>
#include <cstdio>
#include <omp.h>
using I=int64_t;using Z=__int128_t;
constexpr I S=1000000000,ALT=100000000000000000LL;constexpr int LIMIT=4000000,BLOCK=65536;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
I absI(I x){need(x!=INT64_MIN,"absolute range");return x<0?-x:x;}
I mod(I x,I q){x%=q;return x<0?x+q:x;}
I flo(I a,I b){need(b>0,"positive divisor");I q=a/b,r=a%b;return q-(r<0);}
I ceilZ(Z a,Z b){need(b>0,"positive divisor");Z q=a/b,r=a%b;if(r>0)q++;need(q>=INT64_MIN&&q<=INT64_MAX,"quotient range");return I(q);}
I floorZ(Z a,Z b){return -ceilZ(-a,b);}
int bits(I n){need(n>0,"positive integer logarithm");return 64-__builtin_clzll(n);}
std::vector<I>tab,LO,HI;
I logup(I n){int k=bits(n)-1,sh=std::max(k-10,0);I v=I(1)<<sh,i=k>=10?n>>sh:n<<(10-k);I rem=k>=10?n-i*v:0,base=k>=10?i*v:1;return 2*tab[i]+2*I(k)*tab[2048]+ceilZ(Z(rem)*S,base);}
I loglo(I n){int k=bits(n)-1,sh=std::max(k-10,0);I v=I(1)<<sh,i=k>=10?n>>sh:n<<(10-k);I rem=k>=10?n-i*v:0,base=k>=10?i*v:1;I xu=ceilZ(Z(rem)*S,base),xl=floorZ(Z(rem)*S,base);return (i==1024?0:2*(tab[i]-1))+2*I(k)*(tab[2048]-1)+xl-ceilZ(Z(xu)*xu,2*S);}
void initlogs(std::string data){std::ifstream f(data+"/log_table.txt");I v;while(f>>v)tab.push_back(v);need(tab.size()==2049,"log table size");LO.resize(LIMIT+1);HI.resize(LIMIT+1);
 for(int n=1;n<=LIMIT;n++){LO[n]=loglo(n);HI[n]=logup(n);need(LO[n]<=HI[n],"ordered log interval");}
 // Independent logarithm method: rational seeds, then alternating log(1+1/n).
 std::ifstream seeds(data+"/log_seeds.tsv");need(bool(seeds),"secondary seed file");I n=0,lo=0,hi=0,at,l,u;while(seeds>>at>>l>>u){n++;need(at==n&&n<=1024&&l<=u,"ordered secondary seeds");need(Z(LO[n])*(ALT/S)<=l&&u<=Z(HI[n])*(ALT/S),"seed enclosure in primary log");lo=l;hi=u;}need(n==1024&&seeds.eof(),"all secondary seeds");
 for(I j=1024;j<LIMIT;j++){
  Z power=j;I lin=0,uin=0;
  for(I k=1;k<=5;k++){
   Z den=k*power;I dn=floorZ(ALT,den),up=ceilZ(ALT,den);
   if(k<=4)lin+=(k%2?dn:-up);
   uin+=(k%2?up:-dn);
   if(k<5)power*=j;
  }
  lo+=lin;hi+=uin;need(lo<=hi,"ordered independent log recurrence");
  need(Z(LO[j+1])*(ALT/S)<=lo&&hi<=Z(HI[j+1])*(ALT/S),"dense independent logarithm enclosure");
 }
 std::cerr<<"PASS dense logarithms arguments="<<LIMIT<<" independent_seeds=1024 recurrence_steps="<<LIMIT-1024<<"\n";
}
struct Rec{int32_t a,d;};static_assert(sizeof(Rec)==8,"record width");
struct Row{I M,p;};
struct Result{I sumhi=0,num=0,maxh=0,reps=0,firstnum=0;int stencil=0;};
void domain(Row r){need(6330<=r.M&&r.M<LIMIT,"modulus domain");need(211<=r.p&&r.p<20000,"prime domain");for(I d=2;d*d<=r.p;d++)need(r.p%d,"prime input");I m=r.M;for(I p:{I(2),I(3),I(5),r.p}){need(m%p==0,"all four prime factors present");do{m/=p;}while(m%p==0);}need(m==1,"exact radical 30p");}
void produce(Row row,std::string pre,int stencil,const std::vector<uint8_t>&unit){I M=row.M,q=30*row.p,reps=M/2+1;std::ofstream out(pre+".bin",std::ios::binary);need(bool(out),"witness output");std::vector<Rec>buf(BLOCK);
 for(I off=0;off<reps;off+=BLOCK){int len=std::min<I>(BLOCK,reps-off);
  for(int k=0;k<len;k++){I s=off+k,ux=M,uy=0,vx=s,vy=1;bool term=false;
   for(int it=0;it<100;it++){
    I nu=ux*ux+uy*uy,nv=vx*vx+vy*vy;if(nu>nv){std::swap(ux,vx);std::swap(uy,vy);nu=nv;}need(nu>0,"zero reduced row");I dot=ux*vx+uy*vy,z=flo(2*dot+nu,2*nu);if(!z){term=true;break;}vx-=z*ux;vy-=z*uy;need(absI(vx)<=M&&absI(vy)<=M,"reduced row range");
   }
   need(term,"Gauss termination");need(absI(ux*vy-uy*vx)==M,"Gauss determinant");I aa=s,dd=1,best=std::max<I>(s,1);
   auto test=[&](I i,I j){I d=i*uy+j*vy;if(!unit[mod(d,q)])return;I a=i*ux+j*vx;if(d<0){d=-d;a=-a;}I h=std::max(absI(a),d);if(h<best){best=h;aa=a;dd=d;}};
   test(1,0);
   if(stencil){I nu=ux*ux+uy*uy,dot=ux*vx+uy*vy;for(I j=1;j<=6;j++){I a=flo(-j*dot,nu);for(I i=a-6;i<=a+6;i++)test(i,j);}}
   else if(!unit[mod(uy,q)]){for(I i=-3;i<=3;i++)test(i,1);for(I i=-1;i<=1;i++)test(i,2);test(-1,3);test(1,3);test(-1,4);test(1,4);}
   need(dd>0&&best<=M/2,"witness range");buf[k]={int32_t(aa),int32_t(dd)};
  }
  out.write((char*)buf.data(),std::streamsize(len)*8);need(bool(out),"witness write");
 }
}
Result consume(Row row,std::string pre){I M=row.M,reps=M/2+1;std::ifstream in(pre+".bin",std::ios::binary|std::ios::ate);need(bool(in),"missing witnesses");need(in.tellg()>=std::streamoff(reps)*8,"truncated witnesses");need(in.tellg()==std::streamoff(reps)*8,"extra witnesses");in.seekg(0);std::vector<Rec>buf(BLOCK);Result r;I mass=0;
 for(I off=0;off<reps;off+=BLOCK){int len=std::min<I>(BLOCK,reps-off);in.read((char*)buf.data(),std::streamsize(len)*8);need(in.gcount()==std::streamsize(len)*8,"truncated witness records");for(int k=0;k<len;k++){I s=off+k,a=buf[k].a,d=buf[k].d;need(d>0&&d<=M/2,"denominator range");need(std::gcd(d,M)==1,"full denominator unit");need((a-s*d)%M==0,"state congruence");I h=std::max(absI(a),d);need(h>=1&&h<=M/2,"height range");I w=(s==0||2*s==M)?1:2;r.sumhi+=w*HI[h];r.maxh=std::max(r.maxh,h);mass+=w;r.reps++;}}
 need(mass==M&&r.reps==reps,"complete state mass");need(in.peek()==EOF,"extra witness tail");r.num=2*r.sumhi-M*LO[M];return r;
}
std::vector<Row> load(std::string path){std::ifstream f(path);need(bool(f),"domain file");std::vector<Row>rr;std::set<I>seen;std::string line;while(std::getline(f,line)){std::istringstream st(line);Row r;std::string x;need(bool(st>>r.M>>r.p)&&!(st>>x),"domain fields");domain(r);need(seen.insert(r.M).second,"duplicate modulus");rr.push_back(r);}need(f.eof()&&!rr.empty(),"complete nonempty domain");return rr;}
int main(int argc,char**argv){try{
 uint16_t e=1;need(*(uint8_t*)&e==1,"little endian");need(argc>=2,"mode");std::string mode=argv[1];
 if(mode=="logs"){need(argc==3,"logs data");initlogs(argv[2]);std::cout<<"PASS logarithms\n";return 0;}
 if(mode=="produce"||mode=="consume"){need(argc==7,"mode M p data prefix stencil");Row r{std::stoll(argv[2]),std::stoll(argv[3])};domain(r);int stencil=std::stoi(argv[6]);need(stencil==0||stencil==1,"stencil");if(mode=="produce"){std::vector<uint8_t>unit(30*r.p);for(I i=0;i<30*r.p;i++)unit[i]=(i%2&&i%3&&i%5&&i%r.p);produce(r,argv[5],stencil,unit);std::cout<<"PASS produced\n";}
 else{initlogs(argv[4]);auto z=consume(r,argv[5]);std::cout<<"PASS consumed "<<z.num<<" "<<z.reps<<"\n";}return 0;}
 need(mode=="run"&&argc==7,"run domain data temporary-directory output threads");int nt=std::stoi(argv[6]);need(1<=nt&&nt<=4,"thread range");omp_set_num_threads(nt);auto rr=load(argv[2]);initlogs(argv[3]);std::string tmp=argv[4],out=argv[5];std::map<I,std::vector<int>>by;for(int i=0;i<(int)rr.size();i++)by[rr[i].p].push_back(i);std::vector<std::pair<I,std::vector<int>>> groups(by.begin(),by.end());std::vector<Result>results(rr.size());int failed=0;std::string failure;auto start=std::chrono::steady_clock::now();
 #pragma omp parallel for schedule(dynamic,1)
 for(int g=0;g<(int)groups.size();g++){try{I p=groups[g].first;std::vector<uint8_t>unit(30*p);for(I i=0;i<30*p;i++)unit[i]=(i%2&&i%3&&i%5&&i%p);
  for(int i:groups[g].second){auto r=rr[i];std::string pre=tmp+"/M"+std::to_string(r.M);produce(r,pre,0,unit);auto z=consume(r,pre);z.firstnum=z.num;
   if(!(Z(85)*std::max<I>(0,z.num)<Z(4)*r.M*LO[r.M])){I first=z.num;
    #pragma omp critical
    {std::cerr<<"INSUFFICIENT stencil0 M="<<r.M<<" numerator="<<first<<"\n";}
    produce(r,pre,1,unit);z=consume(r,pre);z.firstnum=first;z.stencil=1;
   }
   need(Z(85)*std::max<I>(0,z.num)<Z(4)*r.M*LO[r.M],"finite normalized gate");results[i]=z;need(std::remove((pre+".bin").c_str())==0,"delete consumed temporary stream");
  }
 }catch(std::exception&e){
  #pragma omp critical
  {failed=1;failure=std::string("p=")+std::to_string(groups[g].first)+": "+e.what();}}}
 if(failed)throw std::runtime_error(failure);std::ofstream file(out);need(bool(file),"results output");file<<"M,p,states,representatives,stencil,max_height,sum_log_upper_ticks,log_M_lower_ticks,mean_numerator,first_numerator\n";I states=0,reps=0;int refined=0;
 for(size_t i=0;i<rr.size();i++){auto r=rr[i];auto z=results[i];file<<r.M<<","<<r.p<<","<<r.M<<","<<z.reps<<","<<z.stencil<<","<<z.maxh<<","<<z.sumhi<<","<<LO[r.M]<<","<<z.num<<","<<z.firstnum<<"\n";states+=r.M;reps+=z.reps;refined+=z.stencil;}
 need(bool(file),"results complete");std::cout<<"PASS complete moduli="<<rr.size()<<" states="<<states<<" representatives="<<reps<<" refinements="<<refined<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 }catch(std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}}
