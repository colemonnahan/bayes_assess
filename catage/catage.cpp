#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <catage.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nyrs.allocate("nyrs");
  nages.allocate("nages");
  obs_catch_at_age.allocate(1,nyrs,1,nages,"obs_catch_at_age");
  effort.allocate(1,nyrs,"effort");
  M.allocate("M");
  relwt.allocate(2,nages);
}

void model_parameters::initializationfunction(void)
{
  log_q.set_initial_value(-6.7);
  log_popscale.set_initial_value(5.1);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_q.allocate(-1,"log_q");
  log_popscale.allocate(1,"log_popscale");
  log_sel_coff.allocate(1,nages-1,-15.,15.,2,"log_sel_coff");
  log_relpop.allocate(1,nyrs+nages-1,-15.,15.,2,"log_relpop");
  effort_devs.allocate(1,nyrs,-5.,5.,3,"effort_devs");
  log_sel.allocate(1,nages,"log_sel");
  #ifndef NO_AD_INITIALIZE
    log_sel.initialize();
  #endif
  log_initpop.allocate(1,nyrs+nages-1,"log_initpop");
  #ifndef NO_AD_INITIALIZE
    log_initpop.initialize();
  #endif
  F.allocate(1,nyrs,1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(1,nyrs,1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,nyrs,1,nages,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  N.allocate(1,nyrs,1,nages,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  C.allocate(1,nyrs,1,nages,"C");
  #ifndef NO_AD_INITIALIZE
    C.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  recsum.allocate("recsum");
  #ifndef NO_AD_INITIALIZE
  recsum.initialize();
  #endif
  initsum.allocate("initsum");
  #ifndef NO_AD_INITIALIZE
  initsum.initialize();
  #endif
  avg_F.allocate("avg_F");
  predicted_N.allocate(2,nages,"predicted_N");
  ratio_N.allocate(2,nages,"ratio_N");
  pred_B.allocate("pred_B");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  // this is just to ``invent'' some relative average
  // weight at age numbers
  relwt.fill_seqadd(1.,1.);
  relwt=pow(relwt,.5);
  relwt/=max(relwt);
}

void model_parameters::userfunction(void)
{
  f =0.0;
 f=0;
  // example of using FUNCTION to structure the procedure section
  get_mortality_and_survivial_rates();
  get_numbers_at_age();
  get_catch_at_age();
  evaluate_the_objective_function();
  // arbitrary prior on this since so correlated. just for testing purposes
  f+= pow((log_popscale-10)/1, 2);
}

void model_parameters::get_mortality_and_survivial_rates(void)
{
  int i, j;
  // calculate the selectivity from the sel_coffs
  for (j=1;j<nages;j++)
  {
    log_sel(j)=log_sel_coff(j);
  }
  // the selectivity is the same for the last two age classes
  log_sel(nages)=log_sel_coff(nages-1);
  // This is the same as F(i,j)=exp(q)*effert(i)*exp(log_sel(j));
  F=outer_prod(mfexp(log_q)*effort,mfexp(log_sel));
  if (active(effort_devs))
  {
    for (i=1;i<=nyrs;i++)
    {
      F(i)=F(i)*exp(effort_devs(i));
    }
  }
  // get the total mortality
  Z=F+M;
  // get the survival rate
  S=mfexp(-1.0*Z);
}

void model_parameters::get_numbers_at_age(void)
{
  int i, j;
  log_initpop=log_relpop+log_popscale;
  for (i=1;i<=nyrs;i++)
  {
    N(i,1)=mfexp(log_initpop(i));
  }
  for (j=2;j<=nages;j++)
  {
    N(1,j)=mfexp(log_initpop(nyrs+j-1));
  }
  for (i=1;i<nyrs;i++)
  {
    for (j=1;j<nages;j++)
    {
      N(i+1,j+1)=N(i,j)*S(i,j);
    }
  }
  // calculated predicted numbers at age for next year
  for (j=1;j<nages;j++)
  {
    predicted_N(j+1)=N(nyrs,j)*S(nyrs,j);
    ratio_N(j+1)=predicted_N(j+1)/N(1,j+1);
  }
  // calculated predicted Biomass for next year for
  // adjusted profile likelihood
  pred_B=(predicted_N * relwt);
}

void model_parameters::get_catch_at_age(void)
{
  C=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
}

void model_parameters::evaluate_the_objective_function(void)
{
  // penalty functions to ``regularize '' the solution
  f+=.01*norm2(log_relpop);
  avg_F=sum(F)/double(size_count(F));
  if (last_phase())
  {
    // a very small penalty on the average fishing mortality
    f+= .001*square(log(avg_F/.2));
  }
  else
  {
    f+= 1000.*square(log(avg_F/.2));
  }
  f+=0.5*double(size_count(C)+size_count(effort_devs))
    * log( sum(elem_div(square(C-obs_catch_at_age),.01+C))
    + 0.1*norm2(effort_devs));
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << "Estimated numbers of fish " << endl;
  report << N << endl;
  report << "Estimated numbers in catch " << endl;
  report << C << endl;
  report << "Observed numbers in catch " << endl;
  report << obs_catch_at_age << endl;
  report << "Estimated fishing mortality " << endl;
  report << F << endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}