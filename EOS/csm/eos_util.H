#ifndef _eos_util_
#define _eos_util_

#include <fundamental_constants.H>
#include <vector>


// Ionization energies (in eV) taken from NIST(2020)
//extern std::vector<std::vector<Real>> chi_eV = {
constexpr Real chi_eV[4][8] = {
	{13.598434599702,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, // Hydrogen
	{24.587389011, 54.417765486,0.0,0.0,0.0,0.0,0.0,0.0}, // Helium
	{11.2602880, 24.383154, // Carbon
	 47.88778, 64.49352,
	 392.090518, 489.993198,0.0,0.0},
	{13.618055, 35.12112, // Oxygen
	 54.93554, 77.41350,
	 113.8990, 138.1189,
	 739.32683, 871.409883}};


// these are approximate since they depend
// on the partition function of the i and i+1th states
// which in turn depend on temperature due to exp(-dE/kT)
// g is only of order unity so shouldnt make too big of a deal
//extern std::vector<std::vector<Real>> gfac = {
constexpr Real gfac[4][8] = {
	{0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
	{2.0, 0.5,0.0,0.0,0.0,0.0,0.0,0.0},
	{0.6,0.2,2.0,0.5,2.0,0.5,0.0,0.0},
	{0.5,2.0,0.6,0.2,2.0,0.5,2.0,0.5}};


constexpr int Z[4] = {1,2,6,8}; // atomic numbwers
constexpr int A[4] = {1,4,12,16}; // atomic weights
constexpr int num_I[4] = {2,3,7,9}; // number of ion species to track
	

constexpr int n_elems = 4; // currently we use [H, He, C, O]


struct brent_args
{
  Real* f;//holds boltzmann terms
  int k;  // Z = k+1
  Real n; //Number density of species k
  Real T; // temperature
  std::vector<Real> ion_frac; // this is a pointer to ion_fracs inside actual_eos

};

//boltzmann function f(T) for jth ion of species k
// this can get quite small depending on T
// so might want to floor it at like 1e-200
// to prevent underflow
inline
Real f_kj(int k, int j, Real T)
{
  Real kT = C::k_B*T;
  Real lf = 2.0_rt*3.14159265359_rt*C::m_e*kT;
  lf /= (C::hplanck*C::hplanck);
  Real lam_fac = pow(lf,1.5_rt);
  Real exp_fac = exp(-chi_eV[k][j]*C::ev2erg/kT);
  return 2.0_rt*gfac[k][j]*lam_fac*exp_fac;
}

// partition ion frac term
inline
Real b_ki(int k, int i, Real Zbar, Real n, Real* f)
{
  Real prod_sum = 1.0;

  for(int j = 0; j <= i; ++j) prod_sum *= f[j];

  return prod_sum /= pow(Zbar*n,i);
}


// transcendental equation for Zbar for each species
// // see Morozova 2015 (SNEC Code)
inline
Real theta(int k, Real n, Real* f, Real Zbar, Real& ion_0)
{
  ion_0 = 0.0;
  Real U = 1.0;


  for(int i = 1; i<num_I[k]; ++i)
  {
    Real beta = b_ki(k,i,Zbar,n,f);

    U += beta;
    ion_0 += i*beta;
  }

  ion_0 = 1.0_rt/ion_0;
  ion_0 *= Zbar;


  return 1.0_rt - ion_0*U;

}

// Brent-Dekker root finding algorithm for average charge Zbar=z_s
inline
int find_root(Real z_a, Real z_b, brent_args* args, Real& z_s)
{

  // fill boltzmann distribution terms
  for(int j = 1; j < num_I[args->k]; ++j) args->f[j-1] = f_kj(args->k,j-1,args->T);

  Real this_z_a = z_a;
  Real this_z_b = z_b;


  int mflag = 1;

  //
  Real this_ion_frac = args->ion_frac[0];

//  std::cout << "k = " << args->k << ", n = " << args->n << ", f0 = " << args->f[0] << ", ifrac0 = " << this_ion_frac <<"\n";

  Real f_a = theta(args->k,args->n,args->f,this_z_a,this_ion_frac);
  Real f_b =  theta(args->k,args->n,args->f,this_z_b,this_ion_frac);

 // std::cout << "f_a = " << f_a << "\n";
  //std::cout << "f_b = " << f_b << "\n";

  if (f_a*f_b >= 0 ) return -1;

  if(std::abs(f_a)<std::abs(f_b)){
    this_z_a = z_b;
    this_z_b = z_a;
	Real f_tmp = f_a;
	f_a = f_b;
	f_b = f_tmp;
  }

  Real z_c = this_z_a;
  Real f_c = f_a;
  Real z_d;
  Real f_s;

  Real delta = 1.e-6;//not sure how delta is different from epsilon..

  int MAXITER = 1000;
  int n_iter = 0;

  // main convergence loop until relative Zbar error is below threshold
  while((std::abs(this_z_b-this_z_a)/this_z_b > delta) && (n_iter < MAXITER))
  {

    if((f_a != f_c) && (f_b != f_c)){ // inverse quadratic interpolation

      z_s = this_z_a*f_b*f_c/(f_a-f_b)/(f_a-f_c);
      z_s += this_z_b*f_a*f_c/(f_b-f_a)/(f_b-f_c);
      z_s += z_c*f_a*f_b/(f_c-f_a)/(f_c-f_b);
    }
    else z_s = this_z_b-f_b*(this_z_b-this_z_a)/(f_b-f_a); // secant

    if (((z_s < (3.*this_z_a+this_z_b)/4.) || z_s > this_z_b)
       || (mflag > 0  && std::abs(z_s-this_z_b)>=0.5_rt*std::abs(this_z_b-z_c))
       || (mflag <= 0  && std::abs(z_s-this_z_b) >= std::abs(z_c-z_d)/2.0)
       || (mflag > 0 && std::abs(this_z_b-z_c)<delta)
       || (mflag <= 0 && std::abs(z_c-z_d)<delta)){
         z_s = 0.5_rt*(this_z_a+this_z_b);
       }
       else mflag = -1;

       f_s = theta(args->k,args->n,args->f,z_s,this_ion_frac);

       z_d = z_c;
       z_c = this_z_b;

       if(f_a*f_s < 0) this_z_b = z_s;
       else this_z_a = z_s;
       if(std::abs(f_a)<std::abs(f_b)){
         Real tmp = this_z_a;
         this_z_a = this_z_b;
         this_z_b = tmp;
       }

       n_iter++;

  }

  if(n_iter >= MAXITER) return 0;

  args->ion_frac[0] = this_ion_frac;
	//technically these are not normalized (i.e. species do not sum to 1)
	// but trying to get them to be normalized means you are adding
	// ion fracs close to 1 with ion fracs clsoe to 1e-200
	// so it's better to just use the ratios
  for(int i = 1; i <= num_I[args->k]; ++i) args->ion_frac[i] = args->ion_frac[i-1]/(z_s*args->n)*args->f[i-1];



  

  return 1;

}

#endif
