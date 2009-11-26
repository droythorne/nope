//#include <globes/globes.h>
#include "./prob_engine.h"

template<typename T> Initialiser<T>::Initialiser(const unsigned capacity=0) 
{
	v_.reserve(capacity);
}

template<typename T> Initialiser<T>& Initialiser<T>::Add(const T& t)
{
	v_.push_back(t);
	return *this;
}







int ProbEngine::my_free_probability_engine()
{
	
	#ifdef VERBOSE 
	printf("my_free_probability_engine\n");
	#endif /*VERBOSE*/

	if (U!=NULL) {gsl_matrix_complex_free(U); U = NULL;}
	if (H0_template!=NULL) { gsl_matrix_complex_free(H0_template);  H0_template = NULL; }
	if (H!=NULL)      { gsl_matrix_complex_free(H);   H = NULL; }
	if (Q!=NULL)      { gsl_matrix_complex_free(Q);   Q = NULL; }
	if (S!=NULL)      { gsl_matrix_complex_free(S);   S = NULL; }
	if (S1!=NULL)      { gsl_matrix_complex_free(S1);   S1 = NULL; }
	if (T0!=NULL)      { gsl_matrix_complex_free(T0);   T0 = NULL; }
	if (NSI_mass_temp!=NULL) { gsl_matrix_complex_free(NSI_mass_temp);   NSI_mass_temp = NULL; }
	if (NSI_flav_temp!=NULL) { gsl_matrix_complex_free(NSI_flav_temp);   NSI_flav_temp = NULL; }
	for(i=0;i<8;i++) {
		if(LAMBDA[i]!=NULL) { gsl_matrix_complex_free(LAMBDA[i]); LAMBDA[i] = NULL;}
	}
	if (lambda!=NULL) {gsl_vector_free(lambda); lambda = NULL;}	
	return 0;
}

int ProbEngine::set_GellMann(void)
{

/*Gell-Mann matrices*/

	LAMBDA[0] = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS,GLB_NU_FLAVOURS);
	gsl_matrix_complex_set(LAMBDA[0],0,1,gsl_complex_rect(1,0));
	gsl_matrix_complex_set(LAMBDA[0],1,0,gsl_complex_rect(1,0));

	LAMBDA[1] = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS,GLB_NU_FLAVOURS);
	gsl_matrix_complex_set(LAMBDA[1],0,1,gsl_complex_rect(0,-1));
	gsl_matrix_complex_set(LAMBDA[1],1,0,gsl_complex_rect(0,1));

	LAMBDA[2] = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS,GLB_NU_FLAVOURS);
	gsl_matrix_complex_set(LAMBDA[2],0,0,gsl_complex_rect(1,0));
	gsl_matrix_complex_set(LAMBDA[2],1,1,gsl_complex_rect(-1,0));

	LAMBDA[3] = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS,GLB_NU_FLAVOURS);
	gsl_matrix_complex_set(LAMBDA[3],0,2,gsl_complex_rect(1,0));
	gsl_matrix_complex_set(LAMBDA[3],2,0,gsl_complex_rect(1,0));

	LAMBDA[4] = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS,GLB_NU_FLAVOURS);
	gsl_matrix_complex_set(LAMBDA[4],0,2,gsl_complex_rect(0,-1));
	gsl_matrix_complex_set(LAMBDA[4],2,0,gsl_complex_rect(0,1));

	LAMBDA[5] = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS,GLB_NU_FLAVOURS);
	gsl_matrix_complex_set(LAMBDA[5],1,2,gsl_complex_rect(1,0));
	gsl_matrix_complex_set(LAMBDA[5],2,1,gsl_complex_rect(1,0));

	LAMBDA[6] = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS,GLB_NU_FLAVOURS);
	gsl_matrix_complex_set(LAMBDA[6],1,2,gsl_complex_rect(0,-1));
	gsl_matrix_complex_set(LAMBDA[6],2,1,gsl_complex_rect(0,1));

	LAMBDA[7] = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS,GLB_NU_FLAVOURS);
	

	gsl_matrix_complex_set(LAMBDA[7],0,0,gsl_complex_rect(1,0));
	gsl_matrix_complex_set(LAMBDA[7],1,1,gsl_complex_rect(1,0));
	gsl_matrix_complex_set(LAMBDA[7],2,2,gsl_complex_rect(-2,0));

	gsl_matrix_complex_scale(LAMBDA[7],gsl_complex_rect(1./M_SQRT3,0));	

	return 0;
}

int ProbEngine::free_GellMann(void)
{
	
	//Is this strictly necessary?
	return 0;
}

int ProbEngine::my_init_probability_engine()
{

	#ifdef VERBOSE
	printf(" my_init_probability_engine()\n");
	#endif /*VERBOSE*/

	//my_free_probability_engine();
	U = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS,GLB_NU_FLAVOURS);
	H = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
	H0_template = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
	Q = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
	S = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
	S1 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
	T0 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
	NSI_mass_temp = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
	NSI_flav_temp = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);

	set_GellMann();
	compute_NSI_matrix(NSI_flav_temp, gmann_flav);
	compute_NSI_matrix(NSI_mass_temp, gmann_mass);

	lambda = gsl_vector_alloc(GLB_NU_FLAVOURS);

	return 0;
}


/***************************************************************************
 * Function compute_NSI_hamiltonian                                        *
 ***************************************************************************
 * Calculates the NSI interaction part of Hamiltonian from Gell-Mann       *
 * matrices, and stores the result in NSI_temp.                            *
 **************************************************************************/
int ProbEngine::compute_NSI_matrix(gsl_matrix_complex *NSI_temp, vector<double> gmann)
{
	gsl_matrix_complex_scale(LAMBDA[0],gsl_complex_rect(gmann[0],0.));
	gsl_matrix_complex_add(NSI_temp, LAMBDA[0]);
	gsl_matrix_complex_scale(LAMBDA[1],gsl_complex_rect(gmann[1],0.));
	gsl_matrix_complex_add(NSI_temp, LAMBDA[1]);
		
	gsl_matrix_complex_scale(LAMBDA[2],gsl_complex_rect(gmann[2],0.));
	gsl_matrix_complex_add(NSI_temp, LAMBDA[2]);

	gsl_matrix_complex_scale(LAMBDA[3],gsl_complex_rect(gmann[3],0.));
	gsl_matrix_complex_add(NSI_temp, LAMBDA[3]);

	gsl_matrix_complex_scale(LAMBDA[4],gsl_complex_rect(gmann[4],0.));
	gsl_matrix_complex_add(NSI_temp, LAMBDA[4]);
	gsl_matrix_complex_scale(LAMBDA[5],gsl_complex_rect(gmann[5],0.));
	gsl_matrix_complex_add(NSI_temp, LAMBDA[5]);
	gsl_matrix_complex_scale(LAMBDA[6],gsl_complex_rect(gmann[6],0.));
	gsl_matrix_complex_add(NSI_temp, LAMBDA[6]);
	gsl_matrix_complex_scale(LAMBDA[7],gsl_complex_rect(gmann[7],0.));
	gsl_matrix_complex_add(NSI_temp, LAMBDA[7]);
	//gsl_matrix_complex_scale(NSI_temp, gsl_complex_rect(GLB_V_FACTOR,0.0));
	return 0;
}
/***************************************************************************
 * Function modified_hamiltonian_cd                                        *
 ***************************************************************************
 * Calculates the Hamiltonian for neutrinos (cp_sign=1) or antineutrinos   *
 * (cp_sign=-1) with energy E, propagating in matter of density V          *
 * (> 0 even for antineutrinos), additional interactions in NON_STD[]      *
 * and stores the result in H.                				   *
 ***************************************************************************/
int ProbEngine::modified_hamiltonian_cd(double E, double V, int cp_sign)
{
  
  #ifdef VERBOSE
  printf(" modified_hamiltonian_cd(double E, double V, int cp_sign)\n");
  #endif /*VERBOSE*/

  double inv_E = 1.0 / E;
  complex<double> (*_H)[GLB_NU_FLAVOURS]
    = (complex<double> (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(H, 0, 0);
  complex<double> (*_H0_template)[GLB_NU_FLAVOURS]
    = (complex<double> (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(H0_template, 0, 0);
  int i, j;

  if (cp_sign > 0)
  {
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        _H[i][j] = _H0_template[i][j] * inv_E;
  }
  else
  {
    for( i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        _H[i][j] = conj(_H0_template[i][j] * inv_E); /* delta_CP -> -delta_CP */
  }
/*Add the standard matter potential to H*/ 
	 _H[0][0] = _H[0][0] + cp_sign * V;

/*Add non_standard interactions to H (proportional to std. matter effect)*/
	gsl_matrix_complex *NSI_flav_norm_temp = gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
	gsl_matrix_complex_memcpy(NSI_flav_norm_temp, NSI_flav_temp);
	gsl_matrix_complex_scale(NSI_flav_norm_temp,gsl_complex_rect(V,0.));
  	gsl_matrix_complex_add(H, NSI_flav_norm_temp);
	gsl_matrix_complex_free(NSI_flav_norm_temp);

	return 0;
}

/***************************************************************************
 * Function my_damped_P_matrix_cd                                          *
 ***************************************************************************
 * Calculates the P matrix for neutrino oscillations in matter of constant *
 * density using a fast eigenvalue solver optimized to 3x3 matrices.       *
 * Stores the result in T0. 						   *
 ***************************************************************************
 * Parameters:                                                             *
 *   E: Neutrino energy                                                    *
 *   L: Baseline                                                           *
 *   V: Matter potential (must be > 0 even for antineutrinos)              *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 ***************************************************************************/

int ProbEngine::my_damped_P_matrix_cd(double E, double L, double V, int cp_sign)
{
  
  #ifdef VERBOSE
  printf(" my_damped_P_matrix_cd(double E, double L, double V, int cp_sign)\n");
  #endif /*VERBOSE*/
  /* Introduce some abbreviations */
  complex<double> (*_S)[3]  = (complex<double> (*)[3]) gsl_matrix_complex_ptr(S,0,0);
  complex<double> (*_Q)[3]  = (complex<double> (*)[3]) gsl_matrix_complex_ptr(Q,0,0);
  complex<double> (*_T0)[3] = (complex<double> (*)[3]) gsl_matrix_complex_ptr(T0,0,0);
  double *_lambda = gsl_vector_ptr(lambda,0);
  int status;
  int i, j, k, l;
  gsl_complex plaq_temp, effective_exp;
  double damping_exponent ; 
  if (V < V_THRESHOLD)                             /* Vacuum */
  {
    /* Use vacuum mixing angles and masses */
    #ifdef VERBOSE
    printf("(V < V_THRESHOLD)\n ");
    #endif /*VERBOSE*/
 
    double inv_E = 0.5/E;
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      _lambda[i] = mq[i] * inv_E;

    if (cp_sign > 0)
    {
	#ifdef VERBOSE
      	printf("(cp_sign > 0)\n");
      	#endif /*VERBOSE*/

      	gsl_matrix_complex_memcpy(Q, U);
    }
    else
    {
      
      #ifdef VERBOSE
      printf("(cp_sign < 0)\n");
      #endif /*VERBOSE*/


      complex<double> (*_U)[3]  = (complex<double> (*)[3]) gsl_matrix_complex_ptr(U,0,0);
      for (i=0; i < GLB_NU_FLAVOURS; i++)
        for (j=0; j < GLB_NU_FLAVOURS; j++)
          _Q[i][j] = conj(_U[i][j]);
    }
  }
  else                                             /* Matter */
  {
    #ifdef VERBOSE
    printf("(V > V_THRESHOLD) - MATTER EFFECTS\n ");
    #endif /*VERBOSE*/
 

    complex<double> (*_H)[3] = (complex<double> (*)[3]) gsl_matrix_complex_ptr(H,0,0);
    
    /* Calculate neutrino Hamiltonian */
    if ((status=modified_hamiltonian_cd(E, V, cp_sign)) != 0)
      return status;
    #ifdef VERBOSE
    printf("Calculate eigenvalues/vectors\n");
    #endif /*VERBOSE*/
    /* Calculate eigenvalues of Hamiltonian */
    if ((status=zheevh3(_H, _Q, _lambda)) != 0)
      return status;
  }

  /*Calculate S-Matrix in mass basis in matter ... */
  double phase;
  gsl_matrix_complex_set_zero(S);
  for (i=0; i < GLB_NU_FLAVOURS; i++)
  {
    phase    = -L * _lambda[i];
    _S[i][i] = complex<double>(cos(phase),sin(phase));
#ifdef VERBOSE
printf("_S[i][i]  %e , %e \n", _S[i][i].real(), _S[i][i].imag());
#endif //VERBOSE
  } 
  
  /* calculate the plaquette factored terms of oscillation probability and store in T0 */
  /* NB. Using gsl_matrix_complex _S as container for effective masses is overkill,    */
  /* but presumably no real savings by using a complex double array or vector type.    */


  double mass_ind_damp_exp =  -1.0 * alpha_damp * pow(E, xi_damp - gamma_damp)* pow(L, beta_damp);
#ifdef VERBOSE
printf("%e \n", mass_ind_damp_exp);
#endif //VERBOSE
  gsl_matrix_complex_set_zero(T0);
  complex<double> *p = &_T0[0][0];
  for (i=0; i < GLB_NU_FLAVOURS; i++)              
  {  
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      for (k=0; k < GLB_NU_FLAVOURS; k++)
      {
        for(l=0; l < GLB_NU_FLAVOURS; l++)
	{  
	plaq_temp = gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(Q,i,l), gsl_complex_conjugate(gsl_matrix_complex_get(Q,j,l))),
gsl_complex_mul(gsl_matrix_complex_get(Q,j,k), gsl_complex_conjugate( gsl_matrix_complex_get(Q,i,k)))); 

	effective_exp = gsl_complex_mul(gsl_complex_rect((_S[k][k]).real(),(_S[k][k]).imag()),gsl_complex_rect((_S[l][l]).real(),-(_S[l][l]).imag()));
	printf(" %e, %e \n", GSL_REAL(plaq_temp), GSL_IMAG(plaq_temp));
	printf(" %e \n", effective_exp);
	plaq_temp = gsl_complex_mul(plaq_temp, effective_exp);
	printf(" %e \n", plaq_temp);
	damping_exponent =  mass_ind_damp_exp * pow(2.0 * fabs(gsl_vector_get(lambda, k) - gsl_vector_get(lambda,l)), xi_damp); 
	printf("damping exponent = %e, k = %d, l= %d \n", damping_exponent, k , l);	

	*p += complex<double> (GSL_REAL(plaq_temp), (GSL_IMAG(plaq_temp))) * exp(damping_exponent);

	}
      }
	printf("%d, %d, %e , %e \n",i,j, p->real(), p->imag());
      p++;
    }
  }
    return 0;
}



/***************************************************************************
 * Function my_S_matrix_cd                                                *
 ***************************************************************************
 * Calculates the S matrix for neutrino oscillations in matter of constant *
 * density using a fat seigenvalue solver optimized to 3x3 matrices.       
 ***************************************************************************
 * Parameters:                                                             *
	itd::cout << *p << std::endl;
 *   E: Neutrino energy                                                    *
 *   L: Baseline                                                           *
 *   V: Matter potential (must be > 0 even for antineutrinos)              *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 ***************************************************************************/

int ProbEngine::my_S_matrix_cd(double E, double L, double V, int cp_sign)
{
  
  #ifdef VERBOSE
  printf(" my_S_matrix_cd(double E, double L, double V, int cp_sign)\n");
  #endif /*VERBOSE*/
  /* Introduce some abbreviations */
  complex<double> (*_S)[3]  = (complex<double> (*)[3]) gsl_matrix_complex_ptr(S,0,0);
  complex<double> (*_Q)[3]  = (complex<double> (*)[3]) gsl_matrix_complex_ptr(Q,0,0);
  complex<double> (*_T0)[3] = (complex<double> (*)[3]) gsl_matrix_complex_ptr(T0,0,0);
  double *_lambda = gsl_vector_ptr(lambda,0);
  int status;
  int i, j, k;
  
  if (V < V_THRESHOLD)                             /* Vacuum */
  {
    /* Use vacuum mixing angles and masses */
    #ifdef VERBOSE
    printf("(V < V_THRESHOLD)\n ");
    #endif /*VERBOSE*/
 
    double inv_E = 0.5/E;
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      _lambda[i] = mq[i] * inv_E;

    if (cp_sign > 0)
    {
	#ifdef VERBOSE
      	printf("(cp_sign > 0)\n");
      	#endif /*VERBOSE*/

      	gsl_matrix_complex_memcpy(Q, U);
    }
    else
    {
      
      #ifdef VERBOSE
      printf("(cp_sign < 0)\n");
      #endif /*VERBOSE*/


      complex<double> (*_U)[3]  = (complex<double> (*)[3]) gsl_matrix_complex_ptr(U,0,0);
      for (i=0; i < GLB_NU_FLAVOURS; i++)
        for (j=0; j < GLB_NU_FLAVOURS; j++)
          _Q[i][j] = conj(_U[i][j]);
    }
  }
  else                                             /* Matter */
  {
    #ifdef VERBOSE
    printf("(V > V_THRESHOLD) - MATTER EFFECTS\n ");
    #endif /*VERBOSE*/
 

    complex<double> (*_H)[3] = (complex<double> (*)[3]) gsl_matrix_complex_ptr(H,0,0);
    
    /* Calculate neutrino Hamiltonian */
    if ((status=modified_hamiltonian_cd(E, V, cp_sign)) != 0)
      return status;
    #ifdef VERBOSE
    printf("Calculate eigenvalues/vectors\n");
    #endif /*VERBOSE*/
    /* Calculate eigenvalues of Hamiltonian */
    if ((status=zheevh3(_H, _Q, _lambda)) != 0)
      return status;
  }

  /* Calculate S-Matrix in mass basis in matter ... */
  double phase;
  gsl_matrix_complex_set_zero(S);
  for (i=0; i < GLB_NU_FLAVOURS; i++)
  {
    phase    = -L * _lambda[i];
    _S[i][i] = complex<double>(cos(phase),sin(phase));
  } 
  
  /* ... and transform it to the flavour basis */
  gsl_matrix_complex_set_zero(T0);
  complex<double> *p = &_T0[0][0];
  for (i=0; i < GLB_NU_FLAVOURS; i++)              /* T0 = S.Q^\dagger */
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      for (k=0; k < GLB_NU_FLAVOURS; k++)
      {
        *p += complex<double>( (_S[i][k]).real()*(_Q[j][k]).real() + (_S[i][k]).imag()*(_Q[j][k]).imag() 
                ,  (_S[i][k]).imag()*(_Q[j][k]).real() - (_S[i][k]).real()*(_Q[j][k]).imag() );
      }
      p++;
    }
  gsl_matrix_complex_set_zero(S);
  p = &_S[0][0];
  for (i=0; i < GLB_NU_FLAVOURS; i++)              /* S = Q.T0 */
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      for (k=0; k < GLB_NU_FLAVOURS; k++)
      {
	  *p += complex<double>( (_Q[i][k]).real()*(_T0[j][k]).real() - (_Q[i][k]).imag()*(_T0[j][k]).imag()
                ,  (_Q[i][k]).imag()*(_T0[j][k]).real() + (_Q[i][k]).real()*(_T0[j][k]).imag() );
      }
      p++;
    }

  return 0;
}

ProbEngine::ProbEngine(vector<double>& p, Baseline &b_line,  void *user_data) :

	baseline(b_line),
	th12(p[ANA_THETA_12]),
	th13(p[ANA_THETA_13]),
	th23(p[ANA_THETA_23]),
	deltacp(p[ANA_DELTA_CP]),

/* Assigns masses upto a constant taking into account hierarchy */

	sdm(p[ANA_DM_21]), //* 1.0e-18; /* Convert to GeV^2 */;
	ldm(p[ANA_DM_31]), //* 1.0e-18 Convert to GeV^2 */
	mq(Initialiser<double>(3)
	.Add(fabs(ldm))
	.Add(fabs(ldm) + sdm)
	.Add(fabs(ldm) + ldm)),

	gmann_flav(Initialiser<double>(8)
	.Add(p[GLB_GELLMANN_FLAV_1])
	.Add(p[GLB_GELLMANN_FLAV_2])
	.Add(p[GLB_GELLMANN_FLAV_4])
	.Add(p[GLB_GELLMANN_FLAV_5])
	.Add(p[GLB_GELLMANN_FLAV_6])
	.Add(p[GLB_GELLMANN_FLAV_7])
	.Add(p[GLB_GELLMANN_FLAV_8])),
	
	gmann_mass(Initialiser<double>(8)
	.Add(p[GLB_GELLMANN_MASS_1])
	.Add(p[GLB_GELLMANN_MASS_2])
	.Add(p[GLB_GELLMANN_MASS_4])
	.Add(p[GLB_GELLMANN_MASS_5])
	.Add(p[GLB_GELLMANN_MASS_6])
	.Add(p[GLB_GELLMANN_MASS_7])
	.Add(p[GLB_GELLMANN_MASS_8])),

	alpha_damp(p[GLB_DAMPING_EXPONENT]) {

	int check; //Simple check that initialisation routines are operating correctly
	
	printf("C1 \n");	

	check = my_init_probability_engine();
	if(check!=0) printf("Error in initialisation of ProbEngine\n");
	
	complex<double> (*_U)[GLB_NU_FLAVOURS]
		= (complex<double> (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(U,0,0);	
/* Compute the vacuum mixing matrix */
	
	_U[0][0] = cos(th12)*cos(th13);
  	_U[0][1] = sin(th12)*cos(th13);
	_U[0][2] = std::polar(sin(th13), -deltacp);
	_U[1][0] = -sin(th12)*cos(th23) - std::polar(cos(th12)*sin(th23)*sin(th13), deltacp);
	_U[1][1] =  cos(th12)*cos(th23) - std::polar(sin(th12)*sin(th23)*sin(th13), deltacp);
	_U[1][2] =  sin(th23)*cos(th13);
	_U[2][0] =  sin(th12)*sin(th23) - std::polar(cos(th12)*cos(th23)*sin(th13), deltacp);
	_U[2][1] = -cos(th12)*sin(th23) - std::polar(sin(th12)*cos(th23)*sin(th13), deltacp);
	_U[2][2] =  cos(th23)*cos(th13);
 
  /* Calculate energy independent matrix H0 * E (ie. the Hamiltonian (pre matter effects) in flavour basis) */

	gsl_matrix_complex_set_zero(H0_template);
	gsl_matrix_complex_set_zero(H); /*Don't think this is necessary here, but will leave for now*/
	for (i=0; i < GLB_NU_FLAVOURS; i++)
	gsl_matrix_complex_set(H0_template, i, i, gsl_complex_rect(0.5*mq[i], 0.0));
	
	gsl_matrix_complex *NSI_mass_norm_temp = gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
	gsl_matrix_complex_memcpy(NSI_mass_norm_temp, NSI_mass_temp);
	gsl_matrix_complex_scale(NSI_mass_norm_temp,gsl_complex_rect(ldm,0.0));
	gsl_matrix_complex_add(H0_template, NSI_mass_norm_temp); 
	gsl_matrix_complex_free(NSI_mass_norm_temp);

	gsl_matrix_complex *T = gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
	gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, H0_template, U, /* T=H0.U^\dagger */
                 GSL_COMPLEX_ZERO, T);
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, U, T,             /* H0=U.T */
                 GSL_COMPLEX_ZERO, H0_template);
	gsl_matrix_complex_free(T);
	

}
/***************************************************************************
 * Function my_probability_matrix                                         *
 ***************************************************************************
 * Calculates the neutrino oscillation probability matrix.                 *
 ***************************************************************************
 * Parameters:                                                             *
 *   P:       Buffer for the storage of the matrix                         *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 *   E:       Neutrino energy (in GeV)                                     *
 *   psteps:  Number of layers in the matter density profile               *
 *   length:  Lengths of the layers in the matter density profile in km    *
 *   density: The matter densities in g/cm^3                               *
 *   filter_sigma: Width of low-pass filter or <0 for no filter  DISABLED  *
 *   user_data: Unused here, should be NULL                                *
 ***************************************************************************/
/*Must incorporate the low pass filter at some point*/

int ProbEngine::my_probability_matrix(double P[3][3], int cp_sign, double E, int psteps ,
				const double *length, const double *density, 
				double filter_sigma, void *user_data)
{

	#ifdef VERBOSE
	printf("Calculate the probability matrix\n");
	#endif /*VERBOSE*/
	int status;
	int i,j;


	/* Convert energy to eV */
	E *= 1.0e9;

	if (psteps > 1 )
	{
		#ifdef VERBOSE
		printf("Multiple density steps\n");
		#endif /*VERBOSE*/
/* In the case of multiple density steps, currently (8/12/07) ignoring damping type effects as they are not straightforward to incorporateat the S-matrix level (see Phys Lett B 274, 87-94 (Giunti, Kim, Lee) for a treatment of intrinsic wavepack decoherence for general MSW type Hamiltonians */

		gsl_matrix_complex_set_identity(S1);
		for (i=0; i<psteps;i++)
		{
			status = my_S_matrix_cd(E,GLB_KM_TO_EV(length[i]), density[i]*GLB_V_FACTOR*GLB_Ne_MANTLE, cp_sign);
			if(status != 0)
				return status;
			gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, S, S1,
					GSL_COMPLEX_ZERO, T0); 
			gsl_matrix_complex_memcpy(S1, T0);
		}
		gsl_matrix_complex_memcpy(S, S1);


/* This just computes the squared modulus of the P matrix elements -> enters into rates calculation */
        #ifdef VERBOSE
        printf("Now compute the square moduli of the P matrix via S-matrix\n");
        #endif /*VERBOSE*/
        complex<double> (*_S)[3] = (complex<double> (*)[3]) gsl_matrix_complex_ptr(S,0,0);
        for (i=0; i < GLB_NU_FLAVOURS; i++)
        for (j=0; j < GLB_NU_FLAVOURS; j++)
        P[j][i] = SQR_ABS(_S[i][j]);

	}
/* In the event of constant matter density, compute the damping effects too */

	else
	{
		#ifdef VERBOSE
		printf("Only 1 density steps\n");
		#endif /*VERBOSE*/


		status = my_damped_P_matrix_cd(E, GLB_KM_TO_EV(length[0]), density[0]*GLB_V_FACTOR*GLB_Ne_MANTLE, cp_sign);
      		if (status != 0)
        	return status;

/* This just computes the squared modulus of the P matrix elements -> enters into rates calculation */
 	#ifdef VERBOSE
	printf("Now compute the square moduli of the P matrix by plaquette method\n");
	#endif /*VERBOSE*/
	complex<double> (*_T0)[3] = (complex<double> (*)[3]) gsl_matrix_complex_ptr(T0,0,0);
  	for (i=0; i < GLB_NU_FLAVOURS; i++)
      	for (j=0; j < GLB_NU_FLAVOURS; j++)
        P[i][j] = (_T0[i][j]).real();
	}
	return 0;
}

int ProbEngine::get_probability(double P[3][3], int cp_sign, double E)
{
	int i = my_probability_matrix(P, cp_sign, E, this->baseline.psteps, this->baseline.steplengths, this->baseline.densities, filter_sigma, 0);
	return i;
}
