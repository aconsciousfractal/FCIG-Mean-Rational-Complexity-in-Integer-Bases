// Direct one-step full-mark certificate, with angular coverage kept explicit.
#include "kernel.hpp"
#ifndef BUILD_CONSUMER
#include "search.hpp"
#endif
#include <cstdlib>
#ifndef ANGULAR
#define ANGULAR 4
#endif
#ifndef DEGREE_EXP
#define DEGREE_EXP 1
#endif
constexpr int AG=ANGULAR, AN=2*AG+1;
struct Costs {array<array<I,4>,AN> denlog;};
vector<Costs> costs;
void prepare(){costs.resize(cells.size());for(size_t c=0;c<cells.size();c++)for(int n=-AG;n<=AG;n++){
 int j=0;for(I x:{cells[c].lo,cells[c].hi})for(int tt:{cells[c].tl,cells[c].th}){
  costs[c].denlog[n+AG][j++]=llo(TD*TD)+llo(DEN)+llo(AG*AG+n*n)+llo(K)+llo(x);
 }}}
I angled(int c,W w,int n){I ans=INT64_MIN;int j=0;auto ce=cells[c];for(I x:{ce.lo,ce.hi})for(int tt:{ce.tl,ce.th}){
 Z lin=TD*Z(w.u)+Z(tt)*w.v;
 Z aa=AG*Z(x)*lin+TD*Z(DEN)*n*w.v, bb=-n*Z(x)*lin+TD*Z(DEN)*AG*w.v;
 Z R=max(az(aa),az(bb));
 I val=safe(::ce(2*Z(lup(R))-costs[c].denlog[n+AG][j++],2));ans=max(ans,val);
 }return ans;}
void write32(unsigned char*p,I a){need(a>=INT32_MIN&&a<=INT32_MAX,"coefficient encoding");for(int i=0;i<4;i++)p[i]=(uint32_t(a)>>(8*i))&255;}
I signed32(const unsigned char*p){uint32_t z=0;for(int j=0;j<4;j++)z|=uint32_t(p[j])<<(8*j);return z<=INT32_MAX?I(z):-1-I(~z);}
int main(int argc,char**argv){try{need(argc==5,"produce|consume p data output-prefix");string mode=argv[1];need(mode=="produce"||mode=="consume","mode");bool prod=mode=="produce";
#ifdef BUILD_CONSUMER
 need(!prod,"consumer cannot produce");
#else
 need(prod,"producer cannot consume");
#endif
int p=stoi(argv[2]);init(p,0,0,DEGREE_EXP,string(argv[3])+"/log_table.txt");load(argv[3]);prepare();string pre=argv[4];double st=omp_get_wtime();const char* tc=getenv("CSTAR_THREADS");int nt=tc?stoi(tc):2;need(nt>=1&&nt<=4,"threads");omp_set_num_threads(nt);
 if(!prod){ifstream in(pre+".witness",ios::binary|ios::ate);need(bool(in),"witness file exists");need(in.tellg()>=streamoff(cells.size())*NB*8,"truncated witness file");need(in.tellg()==streamoff(cells.size())*NB*8,"extra witness data");}
 ofstream out;if(prod){out.open(pre+".witness",ios::binary);need(bool(out),"output stream");}
 vector<I> ub(cells.size());vector<int> labs(cells.size()),angs(cells.size());int failed=0;string failure;
 #pragma omp parallel
 {ifstream in;if(!prod)in.open(pre+".witness",ios::binary);vector<unsigned char> records(size_t(NB)*8);vector<W> ws(NB);vector<I> b(NL);
 #pragma omp for schedule(dynamic,1)
 for(int c=0;c<(int)cells.size();c++){try{
  if(!prod){in.seekg(streamoff(c)*NB*8);in.read((char*)records.data(),records.size());need(in.gcount()==streamsize(records.size()),"truncated witness records");}
  for(int h=0;h<NB;h++){
   if(prod){
#ifndef BUILD_CONSUMER
 ws[h]=witness(c,h);write32(&records[8*h],ws[h].ci);write32(&records[8*h+4],ws[h].cj);
#endif
 }
   else{I a=signed32(&records[8*h]),b=signed32(&records[8*h+4]);need(units[mod(Z(a)*br[h].r[0]+Z(b)*br[h].r[1],Q)],"full mark unit");I u=safe(Z(a)*br[h].H[0]+Z(b)*br[h].H[2]),v=safe(Z(a)*br[h].H[1]+Z(b)*br[h].H[3]);need(az(u)<1000000&&az(v)<1000000,"coefficient bounds");ws[h]={u,v,a,b};}}
  if(prod){
   #pragma omp critical
   {out.seekp(streamoff(c)*NB*8);out.write((char*)records.data(),records.size());need(bool(out),"complete witness block");}}
  I best=INT64_MIN;int labmax=-1,angmax=-1;
  for(int n=-AG;n<=AG;n++){fill(b.begin(),b.end(),0);for(int h=0;h<NB;h++)b[br[h].kernel]=safe(Z(b[br[h].kernel])+angled(c,ws[h],n));tensor(b);auto it=max_element(b.begin(),b.end());I val=safe(::ce(*it,K));if(val>best){best=val;labmax=it-b.begin();angmax=n;}}
  ub[c]=best+safe(::ce(S,Z(4)*AG*AG));labs[c]=labmax;angs[c]=angmax;
 }catch(exception&e){
 #pragma omp critical
 {failed=1;failure=string("cell ")+to_string(c)+": "+e.what();}}}}
 if(failed)throw runtime_error(failure);
 ofstream f(pre+(prod?".cells":".checked.cells"));need(bool(f),"profile output");for(int c=0;c<(int)cells.size();c++)f<<c<<" "<<cells[c].lo<<" "<<cells[c].hi<<" "<<ub[c]<<" "<<labs[c]<<" "<<angs[c]<<"\n";
 need(bool(f),"profile complete");
 if(prod)out.flush();
 ofstream audit(pre+(prod?".audit":".checked.audit"));ifstream readw(pre+".witness",ios::binary);need(bool(audit)&&bool(readw),"audit streams");
 for(int z=0;z<96;z++){
  size_t idx=(z==95?size_t(cells.size())*NB-1: (size_t(z)*104729+size_t(z)*z*15485863)%(size_t(cells.size())*NB));
  int c=idx/NB,h=idx%NB;unsigned char a[8];readw.seekg(streamoff(idx)*8);readw.read((char*)a,8);need(readw.gcount()==8,"audit record exists");I ci=signed32(a),cj=signed32(a+4);I u=safe(Z(ci)*br[h].H[0]+Z(cj)*br[h].H[2]),v=safe(Z(ci)*br[h].H[1]+Z(cj)*br[h].H[3]);
  audit<<c<<" "<<h<<" "<<ci<<" "<<cj<<" "<<u<<" "<<v;for(int n=-AG;n<=AG;n++)audit<<" "<<angled(c,{u,v,ci,cj},n);audit<<"\n";
 }
 ofstream algebra(pre+(prod?".branches":".checked.branches"));for(auto h:br)algebra<<h.kernel<<" "<<h.H[0]<<" "<<h.H[2]<<" "<<h.H[3]<<" "<<h.r[0]<<" "<<h.r[1]<<"\n";
 need(bool(audit)&&bool(algebra),"audit outputs complete");cerr<<"PASS "<<mode<<" q="<<Q<<" p="<<p<<" cells="<<cells.size()<<" branches="<<NB<<" angular_nodes="<<AN<<" witnesses="<<size_t(cells.size())*NB<<" seconds="<<omp_get_wtime()-st<<"\n";
 }catch(exception&e){cerr<<"FAIL "<<e.what()<<"\n";return 2;}}
