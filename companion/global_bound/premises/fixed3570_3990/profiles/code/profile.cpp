// Complete 210p marked profiles. Producer emits only coefficients.
// Consumer checks all full-mark coefficients and all continuous rectangles.
#include "kernel.hpp"
#ifndef BUILD_CONSUMER
#include "search.hpp"
#endif
#include <cstdlib>
#ifndef FINE_ANGULAR
#define FINE_ANGULAR 8
#endif
constexpr int FIRST=90,LAST=126,PARTS=4,OLD_COUNT=2667;
constexpr I CROP=5200308914369309LL;
using Costs=vector<array<I,4>>;
I signed32(const unsigned char*p){uint32_t z=0;for(int j=0;j<4;j++)z|=uint32_t(p[j])<<(8*j);return z<=INT32_MAX?I(z):-1-I(~z);}
void write32(unsigned char*p,I a){need(a>=INT32_MIN&&a<=INT32_MAX,"coefficient encoding");for(int i=0;i<4;i++)p[i]=(uint32_t(a)>>(8*i))&255;}
Cell rectangle(int c,int part,int count){Cell x=cells[c];I top=min(x.hi,CROP),width=top-x.lo;need(width>=count&&part>=0&&part<count,"nonempty cropped rectangle");x.hi=x.lo+width*(part+1)/count;x.lo+=width*part/count;return x;}
Costs prepare(Cell x,int ag){Costs ret(2*ag+1);for(int n=-ag;n<=ag;n++){int j=0;for(I rho:{x.lo,x.hi})for(int t:{x.tl,x.th})ret[n+ag][j++]=llo(1600)+llo(DEN)+llo(ag*ag+n*n)+llo(K)+llo(rho);}return ret;}
I angled(Cell x,const Costs&cost,W w,int n,int ag){I ans=INT64_MIN;int j=0;for(I rho:{x.lo,x.hi})for(int t:{x.tl,x.th}){Z lin=40*Z(w.u)+Z(t)*w.v;Z a=ag*Z(rho)*lin+40*Z(DEN)*n*w.v,b=-n*Z(rho)*lin+40*Z(DEN)*ag*w.v;Z R=max(az(a),az(b));I v=safe(::ce(2*Z(lup(R))-cost[n+ag][j++],2));ans=max(ans,v);}return ans;}
struct Result{I cap,fine,embedded;int lab,ang,lab_embedded,ang_embedded;};
Result evaluate(Cell x,const vector<W>&w,int ag,bool embed){Costs cost=prepare(x,ag);vector<I>b(NL);I best=INT64_MIN,even=INT64_MIN;int lab=-1,ang=-1,le=-1,ne=-1;for(int n=-ag;n<=ag;n++){fill(b.begin(),b.end(),0);for(int h=0;h<NB;h++)b[br[h].kernel]=safe(Z(b[br[h].kernel])+angled(x,cost,w[h],n,ag));tensor(b);auto it=max_element(b.begin(),b.end());I z=safe(::ce(*it,K));if(z>best){best=z;lab=it-b.begin();ang=n;}if(embed&&n%2==0&&z>even){even=z;le=it-b.begin();ne=n;}}
 I cf=best+safe(::ce(S,Z(4)*ag*ag));if(!embed)return {cf,cf,cf,lab,ang,lab,ang};I ec=even+safe(::ce(S,Z(ag)*ag));return {min(cf,ec),cf,ec,lab,ang,le,ne};}
int outindex(int c,int part,bool hybrid){if(!hybrid||c<FIRST*21)return c;int j=c-FIRST*21;return FIRST*21+(j/21*PARTS+part)*21+j%21;}
int parentcell(int i,bool hybrid){if(!hybrid||i<FIRST*21)return i;int j=i-FIRST*21;return (FIRST+j/(21*PARTS))*21+j%21;}
int subpart(int i,bool hybrid){return !hybrid||i<FIRST*21?0:(i-FIRST*21)/21%PARTS;}
int main(int argc,char**argv){try{need(argc==6,"produce|consume p coarse|hybrid data prefix");string mode=argv[1],style=argv[3],pre=argv[5];bool prod=mode=="produce",hybrid=style=="hybrid",angular=style=="angular_tail";need(prod||mode=="consume","mode");need(hybrid||angular||style=="coarse","grid style");
#ifdef BUILD_CONSUMER
 need(!prod,"consumer cannot produce");
#else
 need(prod,"producer cannot consume");
#endif
 int p=stoi(argv[2]);need(p==13||p==17||p==19,"declared new-prime domain");init(p,0,0,1,string(argv[4])+"/log_table.txt");load(argv[4]);need((hybrid&&INDEX_MULTIPLIER==6)||(angular&&INDEX_MULTIPLIER==2)||(!hybrid&&!angular&&(INDEX_MULTIPLIER==1||INDEX_MULTIPLIER==2)),"style and index contract");need(3*Z(CROP)*CROP>=4*Z(DEN)*DEN&&3*Z(CROP-1)*(CROP-1)<4*Z(DEN)*DEN,"exact Gauss ceiling");need(cells[LAST*21].lo<CROP&&CROP<cells[LAST*21].hi&&cells[(LAST+1)*21].lo>CROP,"last cropped band");
 int nt=getenv("CSTAR_THREADS")?stoi(getenv("CSTAR_THREADS")):2;need(nt>=1&&nt<=4,"threads");omp_set_num_threads(nt);double st=omp_get_wtime();ofstream out;if(prod){out.open(pre+".witness",ios::binary);need(bool(out),"producer output");}else{ifstream in(pre+".witness",ios::binary|ios::ate);need(bool(in),"witness exists");need(in.tellg()>=streamoff(OLD_COUNT)*NB*8,"truncated witness file");need(in.tellg()==streamoff(OLD_COUNT)*NB*8,"extra witness data");}
 int nc=hybrid?FIRST*21+(LAST-FIRST+1)*PARTS*21:OLD_COUNT;vector<Result>results(nc);int failed=0;string failure;
#pragma omp parallel
 {ifstream in;if(!prod)in.open(pre+".witness",ios::binary);vector<unsigned char>raw(size_t(NB)*8);vector<W>w(NB);
#pragma omp for schedule(dynamic,1)
 for(int c=0;c<OLD_COUNT;c++){try{if(!prod){in.seekg(streamoff(c)*NB*8);in.read((char*)raw.data(),raw.size());need(in.gcount()==streamsize(raw.size()),"complete witness block");}
 for(int h=0;h<NB;h++){if(prod){
#ifndef BUILD_CONSUMER
 W z=witness(c,h);write32(&raw[8*h],z.ci);write32(&raw[8*h+4],z.cj);
#endif
 }else{I a=signed32(&raw[8*h]),b=signed32(&raw[8*h+4]);need(units[mod(Z(a)*br[h].r[0]+Z(b)*br[h].r[1],Q)],"full mark unit");I u=safe(Z(a)*br[h].H[0]+Z(b)*br[h].H[2]),v=safe(Z(a)*br[h].H[1]+Z(b)*br[h].H[3]);need(az(u)<1000000&&az(v)<1000000,"coefficient bounds");w[h]={u,v,a,b};}}
 if(prod){
#pragma omp critical
 {out.seekp(streamoff(c)*NB*8);out.write((char*)raw.data(),raw.size());need(bool(out),"complete producer block");}}
 else{bool refine=(hybrid||angular)&&c>=FIRST*21;int np=refine&&hybrid?PARTS:1,ag=refine?FINE_ANGULAR:4;for(int part=0;part<np;part++)results[outindex(c,part,hybrid)]=evaluate(rectangle(c,part,np),w,ag,refine);}
 }catch(exception&e){
#pragma omp critical
 {failed=1;failure=string("cell ")+to_string(c)+": "+e.what();}}}}
 if(failed)throw runtime_error(failure);
 if(prod){out.flush();need(bool(out),"producer flush");}else{
 ofstream f(pre+".cells");for(int i=0;i<nc;i++){int c=parentcell(i,hybrid),part=subpart(i,hybrid);bool refine=(hybrid||angular)&&c>=FIRST*21;Cell x=rectangle(c,part,refine&&hybrid?PARTS:1);auto z=results[i];int ag=refine?FINE_ANGULAR:4;f<<i<<" "<<c<<" "<<part<<" "<<x.lo<<" "<<x.hi<<" "<<x.tl<<" "<<x.th<<" "<<ag<<" "<<z.cap<<" "<<z.fine<<" "<<z.embedded<<" "<<z.lab<<" "<<z.ang<<" "<<z.lab_embedded<<" "<<z.ang_embedded<<"\n";}need(bool(f),"complete cells");
 ofstream a(pre+".audit");ifstream raw(pre+".witness",ios::binary);for(int z=0;z<96;z++){size_t ix=z==95?size_t(nc)*NB-1:(size_t(z)*104729+size_t(z)*z*15485863)%(size_t(nc)*NB);int i=ix/NB,h=ix%NB,c=parentcell(i,hybrid),part=subpart(i,hybrid);bool refine=(hybrid||angular)&&c>=FIRST*21;int ag=refine?FINE_ANGULAR:4;unsigned char b[8];raw.seekg((streamoff(c)*NB+h)*8);raw.read((char*)b,8);need(raw.gcount()==8,"audit record");I ci=signed32(b),cj=signed32(b+4),u=safe(Z(ci)*br[h].H[0]+Z(cj)*br[h].H[2]),v=safe(Z(cj)*br[h].H[3]);Cell x=rectangle(c,part,refine&&hybrid?PARTS:1);Costs cost=prepare(x,ag);a<<i<<" "<<h<<" "<<ci<<" "<<cj<<" "<<u<<" "<<v<<" "<<ag;for(int n=-ag;n<=ag;n++)a<<" "<<angled(x,cost,{u,v,ci,cj},n,ag);a<<"\n";}need(bool(a),"complete audit");
 ofstream b(pre+".branches");for(auto h:br)b<<h.kernel<<" "<<h.H[0]<<" "<<h.H[2]<<" "<<h.H[3]<<" "<<h.r[0]<<" "<<h.r[1]<<"\n";need(bool(b),"complete branches");}
 cerr<<"PASS "<<mode<<" q="<<Q<<" index="<<K<<" source_cells="<<OLD_COUNT<<" output_cells="<<nc<<" witnesses="<<size_t(OLD_COUNT)*NB<<" seconds="<<omp_get_wtime()-st<<"\n";return 0;
 }catch(exception&e){cerr<<"FAIL "<<e.what()<<"\n";return 2;}}
