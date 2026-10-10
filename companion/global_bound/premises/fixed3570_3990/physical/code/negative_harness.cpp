#define main p72_original_main
#include "finite.cpp"
#undef main
#include <functional>
#include <filesystem>
int rejected=0,positive=0;
void reject(const std::string& name,const char* why,std::function<void()> f){
 try{f();}catch(const std::exception&e){need(std::string(e.what())==why,"unexpected rejection reason");rejected++;std::cout<<"REJECT "<<name<<" | "<<e.what()<<"\n";return;}
 throw std::runtime_error(name+": accepted");
}
int main(int argc,char**argv){try{
 need(argc==3,"data and work directory");std::string work=argv[2];std::filesystem::create_directories(work);
 reject("unsupported_mark","supported radical",[]{domain({2730,2730,2730});});
 initlogs(argv[1]);
 for(I q:{I(3570),I(3990)}){
  std::string tag=std::to_string(q)+":";
  domain({q,q,q});domain({2*q,2*q,q});domain({q*q,q,q});
  reject(tag+"not_base_power","modulus is a power of base",[&]{domain({2*q,q,q});});
  I other=q==3570?3990:7140;
  reject(tag+"wrong_base_mark","all five prime factors present",[&]{domain({q*q,other,q});});
  reject(tag+"extra_prime","exact declared radical",[&]{domain({q*11,q*11,q});});
  reject(tag+"outside_modulus_range","modulus domain",[&]{domain({300000000,q,q});});
  auto df=work+"/domain.txt";
  {std::ofstream f(df);f<<q<<" "<<q<<" "<<q<<"\n"<<q<<" "<<q<<" "<<q<<"\n";}
  reject(tag+"duplicate_modulus","duplicate modulus",[&]{load(df);});
  {std::ofstream f(df);f<<q<<" "<<q<<" "<<q<<" extra\n";}
  reject(tag+"extra_domain_field","domain fields",[&]{load(df);});
  {std::ofstream f(df);f<<q<<" "<<q<<"\n";}
  reject(tag+"missing_mark_field","domain fields",[&]{load(df);});
  Row row{q,q,q};I reps=q/2+1;std::vector<Rec>a(reps);
  for(I s=0;s<reps;s++)a[s]={int32_t(s),1};
  auto pre=work+"/fixture";auto write=[&](const std::vector<Rec>& v){std::ofstream f(pre+".bin",std::ios::binary);f.write((char*)v.data(),v.size()*sizeof(Rec));};
  write(a);auto r=consume(row,pre);I expected=HI[q/2]+HI[1];
  for(I s=1;s<q/2;s++)expected+=2*HI[s];
  need(r.reps==reps&&r.sumhi==expected,"exact endpoint weights");
  need(r.num==2*expected-q*loglo(q),"mean uses logM");positive++;
  std::cout<<"PASS full fixture q="<<q<<" reps="<<reps<<" zero_weight=1 midpoint_weight=1\n";
  auto mutate=[&](I s,Rec v){auto b=a;b[s]=v;write(b);};
  mutate(0,{0,11});auto zero=consume(row,pre);need(zero.sumhi==expected+HI[11]-HI[1],"zero exactly weight one");positive++;
  for(I p:{I(2),I(3),I(5),I(7),q/210}){mutate(1,{int32_t(p),int32_t(p)});reject(tag+"denominator_multiple_"+std::to_string(p),"full denominator unit",[&]{consume(row,pre);});}
  mutate(1,{0,0});reject(tag+"zero_denominator","denominator range",[&]{consume(row,pre);});
  mutate(1,{-1,-1});reject(tag+"negative_denominator","denominator range",[&]{consume(row,pre);});
  mutate(1,{2,1});reject(tag+"wrong_congruence","state congruence",[&]{consume(row,pre);});
  mutate(1,{int32_t(q+1),1});reject(tag+"oversized_height","height range",[&]{consume(row,pre);});
  write(a);std::filesystem::resize_file(pre+".bin",a.size()*sizeof(Rec)-1);reject(tag+"truncated_stream","truncated witnesses",[&]{consume(row,pre);});
  write(a);{std::ofstream f(pre+".bin",std::ios::binary|std::ios::app);f.put('x');}reject(tag+"extra_stream_byte","extra witnesses",[&]{consume(row,pre);});
 }
 need(rejected==37&&positive==4,"fixture cardinalities");
 std::cout<<"PASS all37 bounded negative checks and4 positive consumer fixtures\n";return 0;
}catch(const std::exception&e){std::cerr<<"HARNESS_FAIL "<<e.what()<<"\n";return 1;}}
