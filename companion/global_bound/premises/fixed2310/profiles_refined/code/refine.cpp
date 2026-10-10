// Mark 2310, index 13860. Fixed old-centre witnesses; four radial subrectangles.
#define INDEX_MULTIPLIER 6
#include "kernel.hpp"
#ifndef BUILD_CONSUMER
#include "search.hpp"
#endif
#include <cstdlib>
constexpr int AG=16,AN=33,FIRST=90,LAST=126,PARTS=4,OLD_COUNT=777,NEW_COUNT=3108;
constexpr I CROP=5200308914369309LL;
using Costs=array<array<I,4>,AN>;
I signed32(const unsigned char*p){uint32_t z=0;for(int j=0;j<4;j++)z|=uint32_t(p[j])<<(8*j);return z<=INT32_MAX?I(z):-1-I(~z);}
void write32(unsigned char*p,I a){need(a>=INT32_MIN&&a<=INT32_MAX,"coefficient encoding");for(int i=0;i<4;i++)p[i]=(uint32_t(a)>>(8*i))&255;}
Cell subcell(int c,int part){Cell x=cells[c];I top=min(x.hi,CROP),width=top-x.lo;need(width>=PARTS,"nonempty cropped cell");x.hi=x.lo+width*(part+1)/PARTS;x.lo+=width*part/PARTS;return x;}
Costs prepare(Cell x){Costs ret;for(int n=-AG;n<=AG;n++){int j=0;for(I rho:{x.lo,x.hi})for(int t:{x.tl,x.th})ret[n+AG][j++]=llo(1600)+llo(DEN)+llo(AG*AG+n*n)+llo(K)+llo(rho);}return ret;}
I angled(Cell x,const Costs&cost,W w,int n){I ans=INT64_MIN;int j=0;for(I rho:{x.lo,x.hi})for(int t:{x.tl,x.th}){Z lin=40*Z(w.u)+Z(t)*w.v;Z a=AG*Z(rho)*lin+40*Z(DEN)*n*w.v,b=-n*Z(rho)*lin+40*Z(DEN)*AG*w.v;Z R=max(az(a),az(b));I v=safe(::ce(2*Z(lup(R))-cost[n+AG][j++],2));ans=max(ans,v);}return ans;}
struct Result{I cap,cap16,cap8;int lab16,ang16,lab8,ang8;};
Result evaluate(Cell x,const vector<W>&w){Costs cost=prepare(x);vector<I>b(NL);I best=INT64_MIN,even=INT64_MIN;int l16=-1,n16=-1,l8=-1,n8=-1;for(int n=-AG;n<=AG;n++){fill(b.begin(),b.end(),0);for(int h=0;h<NB;h++)b[br[h].kernel]=safe(Z(b[br[h].kernel])+angled(x,cost,w[h],n));tensor(b);auto it=max_element(b.begin(),b.end());I z=safe(::ce(*it,K));if(z>best){best=z;l16=it-b.begin();n16=n;}if(n%2==0&&z>even){even=z;l8=it-b.begin();n8=n;}}
 I c16=best+safe(::ce(S,1024)),c8=even+safe(::ce(S,256));return {min(c16,c8),c16,c8,l16,n16,l8,n8};}
int main(int argc,char**argv){try{need(argc==4,"produce|consume data prefix");string mode=argv[1],pre=argv[3];bool prod=mode=="produce";need(prod||mode=="consume","mode");
#ifdef BUILD_CONSUMER
 need(!prod,"consumer cannot produce");
#else
 need(prod,"producer cannot consume");
#endif
 init(11,0,0,1,string(argv[2])+"/log_table.txt");load(argv[2]);need(Q==2310&&K==13860&&NL==6912&&NB==41472,"fixed algebra contract");need(3*Z(CROP)*CROP>=4*Z(DEN)*DEN&&3*Z(CROP-1)*(CROP-1)<4*Z(DEN)*DEN,"exact Gauss ceiling");need(cells[LAST*21].lo<CROP&&CROP<cells[LAST*21].hi&&cells[(LAST+1)*21].lo>CROP,"last cropped band");
 int nt=getenv("CSTAR_THREADS")?stoi(getenv("CSTAR_THREADS")):4;need(nt>=1&&nt<=4,"threads");omp_set_num_threads(nt);double st=omp_get_wtime();ofstream out;if(prod){out.open(pre+".witness",ios::binary);need(bool(out),"producer output");}else{ifstream in(pre+".witness",ios::binary|ios::ate);need(bool(in),"witness exists");need(in.tellg()>=streamoff(OLD_COUNT)*NB*8,"truncated witness file");need(in.tellg()==streamoff(OLD_COUNT)*NB*8,"extra witness data");}
 vector<Result>results(NEW_COUNT);int failed=0;string failure;
#pragma omp parallel
 {ifstream in;if(!prod)in.open(pre+".witness",ios::binary);vector<unsigned char>raw(size_t(NB)*8);vector<W>w(NB);
#pragma omp for schedule(dynamic,1)
 for(int j=0;j<OLD_COUNT;j++){try{int c=FIRST*21+j;if(!prod){in.seekg(streamoff(j)*NB*8);in.read((char*)raw.data(),raw.size());need(in.gcount()==streamsize(raw.size()),"complete witness block");}
 for(int h=0;h<NB;h++){if(prod){
#ifndef BUILD_CONSUMER
 W z=witness(c,h);write32(&raw[8*h],z.ci);write32(&raw[8*h+4],z.cj);
#endif
 }else{I a=signed32(&raw[8*h]),b=signed32(&raw[8*h+4]);need(units[mod(Z(a)*br[h].r[0]+Z(b)*br[h].r[1],Q)],"full mark unit");I u=safe(Z(a)*br[h].H[0]+Z(b)*br[h].H[2]),v=safe(Z(a)*br[h].H[1]+Z(b)*br[h].H[3]);need(az(u)<1000000&&az(v)<1000000,"coefficient bounds");w[h]={u,v,a,b};}}
 if(prod){
#pragma omp critical
 {out.seekp(streamoff(j)*NB*8);out.write((char*)raw.data(),raw.size());need(bool(out),"complete producer block");}}
 else for(int part=0;part<PARTS;part++){int ix=(j/21*PARTS+part)*21+j%21;results[ix]=evaluate(subcell(c,part),w);}
 }catch(exception&e){
#pragma omp critical
 {failed=1;failure=string("cell ")+to_string(j+FIRST*21)+": "+e.what();}}}}
 if(failed)throw runtime_error(failure);
 if(prod){out.flush();need(bool(out),"producer flush");}else{
 ofstream f(pre+".cells");for(int i=0;i<NEW_COUNT;i++){int c=(FIRST+i/(21*PARTS))*21+i%21,part=i/21%PARTS;Cell x=subcell(c,part);auto z=results[i];f<<i<<" "<<c<<" "<<part<<" "<<x.lo<<" "<<x.hi<<" "<<x.tl<<" "<<x.th<<" "<<z.cap<<" "<<z.cap16<<" "<<z.cap8<<" "<<z.lab16<<" "<<z.ang16<<" "<<z.lab8<<" "<<z.ang8<<"\n";}need(bool(f),"complete cells");
 ofstream a(pre+".audit");ifstream raw(pre+".witness",ios::binary);for(int z=0;z<96;z++){size_t ix=z==95?size_t(NEW_COUNT)*NB-1:(size_t(z)*104729+size_t(z)*z*15485863)%(size_t(NEW_COUNT)*NB);int i=ix/NB,h=ix%NB,c=(FIRST+i/(21*PARTS))*21+i%21,part=i/21%PARTS;unsigned char b[8];raw.seekg((streamoff(c-FIRST*21)*NB+h)*8);raw.read((char*)b,8);need(raw.gcount()==8,"audit record");I ci=signed32(b),cj=signed32(b+4),u=safe(Z(ci)*br[h].H[0]+Z(cj)*br[h].H[2]),v=safe(Z(cj)*br[h].H[3]);Cell x=subcell(c,part);Costs cost=prepare(x);a<<i<<" "<<h<<" "<<ci<<" "<<cj<<" "<<u<<" "<<v;for(int n=-AG;n<=AG;n++)a<<" "<<angled(x,cost,{u,v,ci,cj},n);a<<"\n";}need(bool(a),"complete audit");
 ofstream b(pre+".branches");for(auto h:br)b<<h.kernel<<" "<<h.H[0]<<" "<<h.H[2]<<" "<<h.H[3]<<" "<<h.r[0]<<" "<<h.r[1]<<"\n";need(bool(b),"complete branches");}
 cerr<<"PASS "<<mode<<" index="<<K<<" old_cells="<<OLD_COUNT<<" subcells="<<NEW_COUNT<<" witnesses="<<size_t(OLD_COUNT)*NB<<" seconds="<<omp_get_wtime()-st<<"\n";return 0;
 }catch(exception&e){cerr<<"FAIL "<<e.what()<<"\n";return 2;}}
