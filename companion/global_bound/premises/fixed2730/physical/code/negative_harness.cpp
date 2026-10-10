#define main p72_original_main
#include "finite.cpp"
#undef main
#include <functional>
#include <filesystem>
void reject(const char* name,const char* why,std::function<void()> f){
 try{f();}catch(const std::exception&e){need(std::string(e.what())==why,"unexpected rejection reason");std::cout<<"REJECT "<<name<<" | "<<e.what()<<"\n";return;}
 throw std::runtime_error(std::string(name)+": accepted");
}
int main(int argc,char**argv){try{
 need(argc==3,"data and work directory");std::string work=argv[2];std::filesystem::create_directories(work);
 domain({2730,2730});domain({5460,5460});domain({7452900,2730});
 std::cout<<"PASS domains N1 and a valid power N2\n";
 reject("not_a_base_power","modulus is a power of base",[]{domain({5460,2730});});
 reject("wrong_base_radical","all five prime factors present",[]{domain({7452900,4620});});
 reject("outside_modulus_range","modulus domain",[]{domain({300000000,2730});});
 auto df=work+"/domain.txt";
 {std::ofstream f(df);f<<"2730 2730\n2730 2730\n";}
 reject("duplicate_modulus","duplicate modulus",[&]{load(df);});
 {std::ofstream f(df);f<<"2730 2730 extra\n";}
 reject("extra_domain_field","domain fields",[&]{load(df);});
 initlogs(argv[1]);Row row{2730,2730};I reps=row.M/2+1;std::vector<Rec>a(reps);
 for(I s=0;s<reps;s++)a[s]={int32_t(s),1};
 auto pre=work+"/fixture";auto write=[&](const std::vector<Rec>& v){std::ofstream f(pre+".bin",std::ios::binary);f.write((char*)v.data(),v.size()*sizeof(Rec));};
 write(a);auto r=consume(row,pre);I expected=HI[row.M/2]+HI[1];
 for(I s=1;s<row.M/2;s++)expected+=2*HI[s];
 need(r.reps==1366&&r.sumhi==expected,"even compression exact endpoint weights");
 need(r.num==2*expected-row.M*loglo(row.M),"mean uses logM");
 std::cout<<"PASS full fixture | reps=1366 mass=2730 zero_weight=1 midpoint_weight=1\n";
 auto mutate=[&](I s,Rec v){auto b=a;b[s]=v;write(b);};
 mutate(0,{0,11});auto zero=consume(row,pre);need(zero.sumhi==expected+HI[11]-HI[1],"zero has exactly weight one");
 std::cout<<"PASS alternative zero witness | weight=1 midpoint retained\n";
 for(int p:{2,3,5,7,13}){mutate(1,{p,p});std::string name="denominator_multiple_"+std::to_string(p);reject(name.c_str(),"full denominator unit",[&]{consume(row,pre);});
 }
 mutate(1,{0,0});reject("zero_denominator","denominator range",[&]{consume(row,pre);});
 mutate(1,{2,1});reject("wrong_congruence","state congruence",[&]{consume(row,pre);});
 mutate(1,{2731,1});reject("oversized_height","height range",[&]{consume(row,pre);});
 write(a);std::filesystem::resize_file(pre+".bin",a.size()*sizeof(Rec)-1);reject("truncated_stream","truncated witnesses",[&]{consume(row,pre);});
 write(a);{std::ofstream f(pre+".bin",std::ios::binary|std::ios::app);f.put('x');}reject("extra_stream_byte","extra witnesses",[&]{consume(row,pre);});
 std::cout<<"PASS all15 bounded negative checks and2 positive consumer fixtures\n";return 0;
 }catch(const std::exception&e){std::cerr<<"HARNESS_FAIL "<<e.what()<<"\n";return 1;}}
