/* Modulus-6 marked mean: complete cell bounds, five horizons.
   See proofs/UNIFORM_MOD6.md. All acceptances use exact integer arithmetic. */
/* Finite-horizon directional witness envelopes. See proofs/UNIFORM_MOD6.md.
   All decisions and all bounds are integer/rational. Geometric centres are input data.
   No floating-point result is used for acceptance. */
#include "kernel.hpp"
#include <cstring>
#include <memory>
#include <sstream>
constexpr int T=24, J=4*T, MINRAD=-125, MAXRAD=2, NS=(MAXRAD-MINRAD+1)*21;
constexpr I MASS=I(1)<<28, LQ=I(1)<<15, FQ=I(1)<<16;
array<Mark,12> reps;
array<array<int,2>,J> directions;
array<I,65536> bend;
vector<Mat> shapes;
I avgQuant, avgRound;
string audit_prefix;
I ceilrat(Z n,Z d){
  return ceildivZ(n,d);
}
I logup(Z n){
  need(n>0,"logup domain");
  int k=bitsZ(n)-1;
  int sh=max(k-10,0);
  Z baseunit=Z(1)<<sh;
  I i=(k>=10)?safe(n>>sh):safe(n<<(10-k));
  Z rem=(k>=10)?n-Z(i)*baseunit:0;
  Z base=(k>=10)?Z(i)*baseunit:1;
  return 2*logs[i]+2*I(k)*logs[2048]+ceilrat(rem*SCALE,base);
}
I loglo(Z n){
  need(n>0,"loglo domain");
  int k=bitsZ(n)-1;
  int sh=max(k-10,0);
  Z baseunit=Z(1)<<sh;
  I i=(k>=10)?safe(n>>sh):safe(n<<(10-k));
  Z rem=(k>=10)?n-Z(i)*baseunit:0;
  Z base=(k>=10)?Z(i)*baseunit:1;
  I xu=ceilrat(rem*SCALE,base),xl=idivfloor(rem*SCALE,base);
  I tab=i==1024?0:2*(logs[i]-1);
  return tab+2*I(k)*(logs[2048]-1)+xl-ceilrat(Z(xu)*xu,2*SCALE);
}
I halve_up(I a){
  return ceilrat(a,2);
}
Z sqrtfloor(Z n){
  need(n>=0,"sqrt nonnegative");
  if(n==0)return 0;
  Z x=Z(1)<<((bitsZ(n)+1)/2);
  while(true){
    Z y=(x+n/x)/2;
    if(y>=x)return x;
    x=y;
  }
}
Z sqrtceil(Z n){
  Z s=sqrtfloor(n);
  return s+(s*s<n);
}
I scaled_ratio_up(Z n,Z d,I scale){
  need(n>=0&&d>0,"ratio positive");
  if(n==0)return 0;
  int shift=max(0,bitsZ(n)-85);
  if(shift){
    Z unit=Z(1)<<shift;
    n=(n+unit-1)/unit;
    d>>=shift;
    need(d>0,"ratio scaled denominator");
  }
  return ceilrat(n*scale,d);
}
Mat relative_numerator(Mat C,Mat Y){
  return {
    safe(Z(C[3])*Y[0]-Z(C[1])*Y[2]),safe(Z(C[3])*Y[1]-Z(C[1])*Y[3]),safe(-Z(C[2])*Y[0]+Z(C[0])*Y[2]),safe(-Z(C[2])*Y[1]+Z(C[0])*Y[3])
  }
  ;
}
I cond_upper(Mat N){
  Z tr=0;
  for(I a:N){
    need(ab(a)<1000000000000000000LL,"relative matrix bound");
    tr+=Z(a)*a;
  }
  Z d=det(N);
  need(d>0,"relative determinant positive");
  I x=scaled_ratio_up(tr,d,SCALE);
  need(x>=2*SCALE&&x<10*SCALE,"shape condition range");
  Z disc=Z(x)*x-4*Z(SCALE)*SCALE;
  I kappa=ceilrat(Z(x)+sqrtceil(disc),2);
  I c=ceilrat(Z(kappa)*kappa,SCALE);
  need(c>=SCALE&&c<4*SCALE,"condition square bound");
  return c;
}
struct Round {
  uint8_t index;
  uint16_t fraction;
  Z norm;
}
;
Round rounddir(Z a,Z b){
  need(a||b,"zero virtual vector");
  Z R=max(abZ(a),abZ(b)),n;
  int idx;
  if(abZ(a)>=abZ(b)){
    if(a<0){
      a=-a;
      b=-b;
    }
    n=T*(a+b);
    idx=int(n/a);
    n%=a;
  }
  else{
    if(b<0){
      a=-a;
      b=-b;
    }
    n=T*(b-a);
    idx=2*T+int(n/b);
    n%=b;
  }
  idx%=J;
  I f=idivfloor(n*FQ,R);
  need(idx>=0&&idx<J&&f>=0&&f<FQ,"direction bin");
  return {
    (uint8_t)idx,(uint16_t)f,R
  }
  ;
}
int shape_of(Mat B){
  Z D=det(B),nu=Z(B[0])*B[0]+Z(B[1])*B[1],dot=Z(B[0])*B[2]+Z(B[1])*B[3];
  need(D>0&&nu>0,"shape domain");
  int lo=0,hi=164;
  while(lo<hi){
    int m=(lo+hi)/2;
    if(nu*PDEN<=Z(radius_bounds[m])*D)hi=m;
    else lo=m+1;
  }
  need(lo<164,"shape upper coverage");
  int k=lo-160,t=int(idivfloor(40*dot+nu,2*nu));
  need(t>=-10&&t<=10,"shape shear coverage");
  if(k<MINRAD)return -1;
  need(k<=MAXRAD,"shape radial upper coverage");
  return (k-MINRAD)*21+t+10;
}
void reduce(Mat&B,Mark&r){
  gauss(B,r);
  if(det(B)<0){
    B[2]=-B[2];
    B[3]=-B[3];
    r[1]=mod(-r[1],6);
  }
  need(det(B)>0,"reduced determinant");
}
I fallback_base(Mat Y){
  Z D=det(Y),nu=Z(Y[0])*Y[0]+Z(Y[1])*Y[1];
  I K=max(ab(Y[0]),ab(Y[1]));
  return logup(Q-1)+halve_up(logup(D))+logup(K)-loglo(nu);
}
I fallback_error(Mat Y){
  Z nu=Z(Y[0])*Y[0]+Z(Y[1])*Y[1];
  return scaled_ratio_up(4*nu,det(Y),SCALE);
}
#pragma pack(push,1)
struct Route {
  uint16_t fraction;
  uint8_t index;
  int16_t correction;
}
;
#pragma pack(pop)
static_assert(sizeof(Route)==5,"packed route contract");
struct Edge {
  int32_t dest;
  uint32_t cond;
  I base;
  uint32_t error;
  uint8_t index;
  uint16_t fraction;
}
;
vector<Edge> edges;
vector<Route> routes;
void load_shapes(string path,string part){
  ifstream f(part);
  need(bool(f),"partition file");
  f>>PDEN;
  for(auto&x:radius_bounds)f>>x;
  for(auto&x:angle_bounds)f>>x;
  need(bool(f),"partition complete");
  for(int i=1;i<164;i++)need(radius_bounds[i]>radius_bounds[i-1],"partition ordered");
  shapes.resize(NS);
  vector<int> seen(NS);
  ifstream c(path);
  int k,t;
  Mat C;
  while(c>>k>>t>>C[0]>>C[1]>>C[2]>>C[3]){
    need(k>=MINRAD&&k<=MAXRAD&&t>=-10&&t<=10,"centre index range");
    int s=(k-MINRAD)*21+t+10;
    need(!seen[s]++,"duplicate centre");
    need(C[0]>0&&C[1]==0&&C[3]>0,"triangular centre");
    for(auto a:C)need(ab(a)<1000000000,"centre size");
    shapes[s]=C;
  }
  for(int n:seen)need(n==1,"centre domain complete");
  for(int l=0;l<12;l++){
    int z=l;
    array<int,2>a,b;
    for(int j=1;j>=0;j--){
      int v=z%(ps[j]+1);
      z/=ps[j]+1;
      a[j]=(v<ps[j]);
      b[j]=(v<ps[j]?v:1);
    }
    reps[l]=crtmark(a,b);
    need(label(reps[l])==l,"mark inverse");
  }
  for(int j=0;j<J;j++)directions[j]= j<2*T?array<int,2>{
    T,j-T
  }
  :array<int,2>{
    3*T-j,T
  }
  ;
  for(I f0=0;f0<FQ;f0++){
    I v=f0<FQ/2?f0+1:f0;
    bend[f0]=ceilrat(Z(SCALE)*v*(FQ-v),Z(FQ)*FQ*T*T);
  }
  avgQuant=ceilrat(Z(SCALE)*(MASS+Z(J)*FQ),Z(MASS)*FQ);
  avgRound=ceilrat(Z(SCALE)*J,MASS);
}
int16_t logcorrection(Mat N,int j,Z D,Z E){
  auto u=directions[j];
  Z a=Z(u[0])*N[0]+Z(u[1])*N[2],b=Z(u[0])*N[1]+Z(u[1])*N[3],R=max(abZ(a),abZ(b));
  I v=logup(R)-loglo(T)+ceilrat(-Z(loglo(D))-loglo(E),2);
  I l=ceilrat(Z(v)*LQ,SCALE);
  need(l>INT16_MIN&&l<INT16_MAX,"log correction storage");
  return (int16_t)l;
}
void make_edges(){
  double st=omp_get_wtime();
  edges.resize(NS*12);
  routes.resize(size_t(NS)*12*J);
#pragma omp parallel for schedule(dynamic,1)
  for(int c=0;c<NS;c++)for(int h=0;h<12;h++){
    size_t e=size_t(c)*12+h;
    auto br=branches[h];
    Mat Y=mm(br.H,shapes[c]);
    Mark r=br.r;
    reduce(Y,r);
    int j=shape_of(Y);
    if(j<0){
      auto rd=rounddir(-Z(Y[1]),Z(Y[0]));
      edges[e]={
        -1,0,fallback_base(Y),(uint32_t)(fallback_error(Y)+bend[rd.fraction]+ceilrat(SCALE,FQ)+ceilrat(SCALE,MASS)),rd.index,rd.fraction
      }
      ;
      continue;
    }
    Mat C=shapes[j];
    Mat N=relative_numerator(C,Y);
    Z D=det(Y),E=det(C);
    edges[e]={
      j*12+label(r),(uint32_t)cond_upper(N),0,0,0,0
    }
    ;
    for(int z=0;z<J;z++){
      auto u=directions[z];
      Z a=Z(u[0])*N[0]+Z(u[1])*N[2],b=Z(u[0])*N[1]+Z(u[1])*N[3];
      auto rd=rounddir(a,b);
      routes[e*J+z]={
        rd.fraction,rd.index,logcorrection(N,z,D,E)
      }
      ;
    }
  }
  I nf=0,maxc=0;
  for(auto&e:edges){
    nf+=(e.dest<0);
    maxc=max(maxc,I(e.cond));
  }
  cerr<<"EDGES_DONE "<<edges.size()<<" FALLBACK "<<nf<<" COND_MAX "<<maxc<<" sec "<<omp_get_wtime()-st<<endl;
}
struct Model{
  vector<I>b,e;
  vector<uint32_t>w;
  Model():b(NS*12),e(NS*12),w(size_t(NS)*12*J) {
  }
}
;
void init_model(Model&m){
  double st=omp_get_wtime();
  vector<array<int16_t,2>> selected(NS*12);
#pragma omp parallel for schedule(dynamic,2)
  for(int c=0;c<NS;c++)for(int l=0;l<12;l++){
    Mat C=shapes[c];
    Mark r=reps[l];
    Z best=Z(1)<<120;
    I vi=0,vj=0,va=0,vb=0;
    auto test=[&](I i,I j){
      if(!units[mod(i*r[0]+j*r[1],6)])return;
      I a=safe(Z(i)*C[0]+Z(j)*C[2]),b=safe(Z(i)*C[1]+Z(j)*C[3]);
      Z n=Z(a)*a+Z(b)*b;
      if(n<best){
        best=n;
        vi=i;
        vj=j;
        va=a;
        vb=b;
      }
    }
    ;
    test(1,0);
    for(int j=1;j<=6;j++){
      I k=idivfloor(-Z(j)*C[2],C[0]);
      for(I i=k-12;i<=k+12;i++)test(i,j);
    }
    need(best<(Z(1)<<120),"unit witness exists");
    need(ab(vi)<32768&&ab(vj)<32768,"initial coefficient storage");
    selected[c*12+l]={
      (int16_t)vi,(int16_t)vj
    }
    ;
    need(units[mod(vi*r[0]+vj*r[1],6)],"initial witness unit");
    auto rd=rounddir(va,vb);
    size_t row=c*12+l;
    I wt=I(MASS)*rd.fraction/FQ;
    m.w[row*J+rd.index]=MASS-wt;
    m.w[row*J+(rd.index+1)%J]+=wt;
    m.b[row]=logup(rd.norm)+ceilrat(-Z(loglo(det(C))),2);
    m.e[row]=bend[rd.fraction]+ceilrat(SCALE,FQ)+ceilrat(SCALE,MASS);
  }
  ofstream wf(audit_prefix+".witnesses",ios::binary);
  wf.write((char*)selected.data(),selected.size()*sizeof(selected[0]));
  need(bool(wf),"witness dump write");
  cerr<<"INIT_DONE sec "<<omp_get_wtime()-st<<endl;
}
void tensorW(vector<I>&w){
  int stride=1;
  for(int z=1;z>=0;z--){
    int n=ps[z]+1;
    for(int start=0;start<12;start+=stride*n)for(int off=0;off<stride;off++){
      array<I,J> sum{
      }
      ;
      for(int t=0;t<n;t++){
        I*row=&w[size_t(start+off+t*stride)*J];
        for(int j=0;j<J;j++)sum[j]+=row[j];
      }
      for(int t=0;t<n;t++){
        I*row=&w[size_t(start+off+t*stride)*J];
        for(int j=0;j<J;j++)row[j]=sum[j]-row[j];
      }
    }
    stride*=n;
  }
}
void advance(const Model&a,Model&b,int k){
  double st=omp_get_wtime();
  I hk=ceilrat(Z(k)*logup(Q),2);
#pragma omp parallel
  {
    vector<I> w(size_t(12)*J);
    array<I,12>bs,es;
#pragma omp for schedule(dynamic,1)
    for(int c=0;c<NS;c++){
      fill(w.begin(),w.end(),0);
      for(int h=0;h<12;h++){
        size_t edge=size_t(c)*12+h;
        const Edge&ed=edges[edge];
        I*wr=&w[size_t(h)*J];
        if(ed.dest<0){
          I wt=MASS*ed.fraction/FQ;
          wr[ed.index]=MASS-wt;
          wr[(ed.index+1)%J]+=wt;
          bs[h]=ed.base+hk;
          es[h]=ed.error;
          continue;
        }
        size_t row=ed.dest;
        I ls=0,ce=0;
        const uint32_t*mw=&a.w[row*J];
        const Route*rr=&routes[edge*J];
        for(int j=0;j<J;j++){
          I mass=mw[j];
          if(!mass)continue;
          auto route=rr[j];
          I w1=mass*route.fraction/FQ;
          wr[route.index]+=mass-w1;
          wr[(route.index+1)%J]+=w1;
          ls+=mass*route.correction;
          ce+=mass*bend[route.fraction];
        }
        bs[h]=a.b[row]+ceilrat(Z(ls)*SCALE,Z(MASS)*LQ);
        es[h]=ceilrat(Z(a.e[row])*ed.cond,SCALE)+ceilrat(ce,MASS)+avgQuant;
      }
      tensor(bs);
      tensor(es);
      tensorW(w);
      for(int l=0;l<12;l++){
        size_t row=c*12+l;
        I sum=0;
        for(int j=0;j<J;j++){
          I wt=w[size_t(l)*J+j]/Q;
          need(wt>=0&&wt<=MASS,"direction mass range");
          b.w[row*J+j]=(uint32_t)wt;
          sum+=wt;
        }
        need(sum<=MASS&&MASS-sum<J,"direction mass closure");
        b.w[row*J]+=uint32_t(MASS-sum);
        b.b[row]=ceilrat(bs[l],Q);
        b.e[row]=ceilrat(es[l],Q)+avgRound;
        need(b.e[row]>=0&&b.e[row]<1000000000000000LL,"error coefficient range");
        need(ab(b.b[row])<1000000000000LL,"baseline range");
      }
    }
  }
  cerr<<"STEP "<<k+1<<" sec "<<omp_get_wtime()-st<<endl;
}
I point_upper(const Model&m,Mat Y,Mark r,int k){
  reduce(Y,r);
  int j=shape_of(Y);
  if(j<0)return fallback_base(Y)+fallback_error(Y)+ceilrat(Z(k)*logup(Q),2);
  Mat C=shapes[j],N=relative_numerator(C,Y);
  size_t row=j*12+label(r);
  I ls=0;
  for(int z=0;z<J;z++)ls+=I(m.w[row*J+z])*logcorrection(N,z,det(Y),det(C));
  return m.b[row]+ceilrat(Z(ls)*SCALE,Z(MASS)*LQ)+ceilrat(Z(m.e[row])*cond_upper(N),SCALE);
}
void write_audit(){
  ofstream f(audit_prefix+".audit");
  for(int c=0;c<NS;c++)for(int h=0;h<12;h++){
    size_t e=size_t(c)*12+h;
    // Every geometric edge is emitted for the independent Python check.
    auto br=branches[h];
    Mat Y=mm(br.H,shapes[c]);
    Mark r=br.r;
    reduce(Y,r);
    int d=shape_of(Y);
    f<<e<<" "<<d<<" "<<label(r);
    for(auto v:Y)f<<" "<<v;
    auto ed=edges[e];
    f<<" "<<ed.dest<<" "<<ed.cond<<" "<<ed.base<<" "<<ed.error<<" "<<int(ed.index)<<" "<<ed.fraction;
    if(d>=0){
      for(int j:{
        0,T-1,2*T-1,J-1
      }
      ){
        auto rt=routes[e*J+j];
        f<<" "<<j<<" "<<int(rt.index)<<" "<<rt.fraction<<" "<<rt.correction;
      }
    }
    f<<"\n";
  }
  need(bool(f),"audit write");
}
constexpr int AG=16, AN=2*AG+1;
struct Cell { I lo,hi;int tl,th; I K; };
vector<Cell> cells;
vector<I> cellCosts;
void prepare_cells(string path){
  ifstream f(path);int k,t;Cell c;while(f>>k>>t>>c.lo>>c.hi>>c.tl>>c.th>>c.K){
    need(k==MINRAD+int(cells.size()/21)&&t==-10+int(cells.size()%21),"cell ordering");
    need(c.lo==radius_bounds[k+159]&&c.hi==radius_bounds[k+160],"radial endpoints");
    need(c.tl==max(-20,2*t-1)&&c.th==min(20,2*t+1)&&c.K>=SCALE&&c.K<4*SCALE,"cell shear and condition range");
    cells.push_back(c);
  }
  need(cells.size()==NS,"complete cells");
  cellCosts.resize(size_t(NS)*AN*J);
  double st=omp_get_wtime();
#pragma omp parallel for schedule(dynamic,1)
  for(int c=0;c<NS;c++){
    auto C=shapes[c];auto cell=cells[c];
    for(int ii=0;ii<AN;ii++)for(int z=0;z<J;z++){
      int n=ii-AG;I best=INT64_MIN;auto u=directions[z];
      for(I x:{cell.lo,cell.hi})for(int tt:{cell.tl,cell.th}){
        Z A=Z(40)*C[3]*x,B=(Z(C[0])*tt-Z(40)*C[2])*x,D=Z(40)*C[0]*PDEN;
        Z ua=Z(u[0])*A+Z(u[1])*B,ub=Z(u[1])*D;
        Z va=ua*AG+ub*n,vb=-ua*n+ub*AG;
        Z R=max(abZ(va),abZ(vb));need(R>0,"cell norm positive");
        I cost=logup(R)-loglo(T)+ceilrat(-Z(loglo(A))-loglo(D)-loglo(AG*AG+n*n),2);
        best=max(best,cost);
      }
      need(ab(best)<2*SCALE,"cell log-cost bound");
      cellCosts[(size_t(c)*AN+ii)*J+z]=best;
    }
  }
  ofstream audit(audit_prefix+".cellcosts");
  for(size_t id=0;id<cellCosts.size();id++)if(id%997==0||id+1==cellCosts.size()){
    size_t cc=id/(AN*J),ii=(id/J)%AN,zz=id%J;
    audit<<cc<<" "<<ii<<" "<<zz<<" "<<cellCosts[id]<<"\n";
  }
  need(bool(audit),"cell cost audit write");
  cerr<<"CELL_COSTS_DONE sec "<<omp_get_wtime()-st<<endl;
}
void certify_cells(const Model&m,int k){
 double st=omp_get_wtime();vector<I> ub(NS,INT64_MIN);vector<int> wi(NS),wl(NS);
#pragma omp parallel for schedule(dynamic,1)
 for(int c=0;c<NS;c++){
  for(int l=0;l<12;l++){
   size_t row=c*12+l;
   I base=m.b[row]+ceilrat(Z(m.e[row])*cells[c].K,SCALE);
   for(int ii=0;ii<AN;ii++){
    I val=0;auto cost=&cellCosts[(size_t(c)*AN+ii)*J];auto w=&m.w[row*J];
    for(int j=0;j<J;j++)val+=I(w[j])*cost[j];
    I b=base+ceilrat(val,MASS)+ceilrat(SCALE,4*AG*AG);
    if(b>ub[c]){ub[c]=b;wi[c]=ii;wl[c]=l;}
   }
  }
 }
 ofstream f(audit_prefix+".cells"+to_string(k));
 for(int c=0;c<NS;c++)f<<c<<" "<<ub[c]<<" "<<wl[c]<<" "<<wi[c]<<"\n";
 ofstream r(audit_prefix+".radial"+to_string(k));
 for(int j=0;j<128;j++){I v=INT64_MIN;int at=-1;for(int t=0;t<21;t++)if(ub[j*21+t]>v){v=ub[j*21+t];at=j*21+t;}
   r<<j-125<<" "<<cells[j*21].lo<<" "<<cells[j*21].hi<<" "<<v<<" "<<at<<"\n";}
 cerr<<"CELL_CERT_DONE k "<<k<<" sec "<<omp_get_wtime()-st<<endl;
}

int main(int argc,char**argv){
  try{
    need(argc==7,"arguments log-table partition shapes output max-k cells");
    init(argv[1]);
    load_shapes(argv[3],argv[2]);
    int maxk=stoi(argv[5]);
    need(maxk==5,"exactly five certificate horizons");
    double st=omp_get_wtime();
    audit_prefix=argv[4];
    make_edges();
    write_audit();
    prepare_cells(argv[6]);
    Model a,b;
    init_model(a);
    ofstream bank(audit_prefix+".bank");
    ofstream points(audit_prefix+".points");
    ofstream out(argv[4]);
    out<<"N,baseline_num,angular_num,error_num,total_num,den,diagnostic_below_32_over_255\n";
    Mat I0={
      1,0,0,1
    }
    ;
    Mark r0={
      0,1
    }
    ;
    int c0=shape_of(I0),l0=label(r0);
    Mat C0=shapes[c0],N0=relative_numerator(C0,I0);
    I K0=cond_upper(N0);
    for(int k=0;k<=maxk;k++){
      I root=point_upper(a,I0,r0,k);
      if(k>0)certify_cells(a,k);
      points<<k<<" "<<root<<"\n";
      if(k<=1)for(I ss=0;ss<Q;ss++)bank<<k<<" "<<ss<<" "<<point_upper(a,{
        Q,0,ss,1
      }
      ,{
        0,1
      }
      ,k)<<"\n";
      size_t row=c0*12+l0;
      I ls=0;
      for(int z=0;z<J;z++)ls+=I(a.w[row*J+z])*logcorrection(N0,z,1,det(C0));
      I base=a.b[row],angular=ceilrat(Z(ls)*SCALE,Z(MASS)*LQ),error=ceilrat(Z(a.e[row])*K0,SCALE);
      need(root==base+angular+error,"root budget identity");
      bool pass=Z(255)*root<32*Z(SCALE);
      out<<k<<","<<base<<","<<angular<<","<<error<<","<<root<<","<<SCALE<<","<<pass<<"\n";
      out.flush();
      cerr<<setprecision(14)<<"ROOT_BOUND N="<<k<<" upper="<<(long double)root/SCALE<<" base="<<(long double)base/SCALE<<" error="<<(long double)error/SCALE<<" pass="<<pass<<endl;
      if(k<maxk){
        advance(a,b,k);
        swap(a,b);
      }
    }
    cerr<<"COMPLETE sec "<<omp_get_wtime()-st<<endl;
  }
  catch(exception&e){
    cerr<<"ERROR "<<e.what()<<endl;
    return 2;
  }
}
