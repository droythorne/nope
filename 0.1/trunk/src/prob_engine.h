
#ifndef PROB_ENGINE_H
#define PROB_ENGINE_H
//#include <globes/globes.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <vector>
#include <complex> 
#include <cstdlib>
#include "./eigenfunctions.h"

/* Numerical Constants */
//#define GLB_V_FACTOR        7.56e-14   /* Conversion factor for matter potentials */
#define GLB_V_FACTOR        7.5e-14   /* Conversion factor for matter potentials */
//#define GLB_Ne_MANTLE       0.497      /* Effective electron numbers for calculation */
#define GLB_Ne_MANTLE       0.5        /* Effective electron numbers for calculation */
#define GLB_Ne_CORE         0.468      /*   of MSW potentials                        */
#define V_THRESHOLD 0.001*GLB_V_FACTOR*GLB_Ne_MANTLE  /* The minimum matter potential below */

#define GLB_NU_FLAVOURS	3
#define M_SQRT3  1.73205080756887729352744634151      /* sqrt(3) */
//#define M_PI  3.1415926536
//#define NULL 0
#define GLB_EARTH_RADIUS 6371.0 /* km */

/* Unit conversion */

#define GLB_EV_TO_KM_FACTOR  1.9747235e-10
#define GLB_KM_TO_EV(x)      ((x) / GLB_EV_TO_KM_FACTOR)
#define GLB_EV_TO_KM(x)      ((x) * GLB_EV_TO_KM_FACTOR)

/*Macros*/

#define SQR(x)		((x)*(x))
#define SQR_ABS(x)	(SQR((x).real()) + SQR((x).imag()))
using std::vector;

enum glb_enum_STD_oscillation_parameters
       { ANA_THETA_12 = 0, ANA_THETA_13 = 1, ANA_THETA_23 = 2,
         ANA_DELTA_CP = 3, ANA_DM_21 = 4, ANA_DM_31 = 5 };

enum glb_enum_NSI_oscillation_parameters {
	 GLB_GELLMANN_FLAV_1 = 6, /* Set the index of the new params. */
	 GLB_GELLMANN_FLAV_2 = 7,
	 GLB_GELLMANN_FLAV_3 = 8,
	 GLB_GELLMANN_FLAV_4 = 9,
	 GLB_GELLMANN_FLAV_5 = 10,
	 GLB_GELLMANN_FLAV_6 = 11,
	 GLB_GELLMANN_FLAV_7 = 12,
	 GLB_GELLMANN_FLAV_8 = 13,

	 GLB_GELLMANN_MASS_1 = 14, /* Set the index of the new params. */
	 GLB_GELLMANN_MASS_2 = 15,
	 GLB_GELLMANN_MASS_3 = 16,
	 GLB_GELLMANN_MASS_4 = 17,
	 GLB_GELLMANN_MASS_5 = 18,
	 GLB_GELLMANN_MASS_6 = 19,
	 GLB_GELLMANN_MASS_7 = 20,
	 GLB_GELLMANN_MASS_8 = 21,

	 GLB_DAMPING_EXPONENT = 22 };



template<typename T> 
class Initialiser
{
vector<T> v_; 
public:
Initialiser( const unsigned );
Initialiser& Add( const T& t );
operator vector<T>() const {return v_;}

};
class Baseline {
public:
	int psteps;
	const double *densities;
	const double *steplengths;

	Baseline(int psteps, const double *densities, const double *steplengths):
		psteps(psteps),
		densities(densities),
		steplengths(steplengths) 
	{}
};
class ProbEngine {

private:

        int i; /* Counter for assignment loop */
	static const int filter_sigma = 1;

	Baseline &baseline;
        gsl_matrix_complex *U; /* Mass <-> Flavour basis vacuum mixing matrix (no NSIs) */
        gsl_matrix_complex *H; /* Neutrino Hamiltonian                               */
        gsl_matrix_complex *H0_template;  /* Used in the construction of the vac. Hamiltonian */
        gsl_matrix_complex *Q; /* Eigenvectors of Hamiltonian (= eff. mixing matrix) */
        gsl_matrix_complex *S; /* The neutrino S-matrix                              */

        gsl_matrix_complex *S1, *T0; /* Temporary matrix storage               */
        gsl_matrix_complex *NSI_flav_temp; /*Matrix containing NSI additions to matter potential*/
        gsl_matrix_complex *NSI_mass_temp; /*Matrix containing NSI additions to matter potential*/
        /*Gell-Mann matrices*/
        gsl_matrix_complex *LAMBDA[8];


        gsl_vector *lambda;    /* Eigenvalues of Hamiltonian                         */
        /* Set up doubles to hold the parameters  */

        double th12;
        double th13;
        double th23;
        double deltacp;
        double sdm;
        double ldm;
        vector<double> mq;
        vector<double> gmann_flav;
        vector<double> gmann_mass;
        double alpha_damp; /* Sets the scale of the damping effects */

        /* Parameters which fix the damping model under consideration */
        /* Default to wavepacket decoherence type damping */
        /* Exponent of baseline dependence (for single baseline, not "observable" as can absorb baseline dependence into alpha) */
        static const int beta_damp = 2;
        static const int gamma_damp = 4; /* Exponent of E dependence */
        static const int  xi_damp = 2; /* Exponent of mass squared difference dependence (>0 implies only oscillating terms are damped */

        int my_free_probability_engine();
        int set_GellMann(void);
        int free_GellMann(void);
        int my_init_probability_engine();
        int compute_NSI_matrix(gsl_matrix_complex *NSI_temp, vector<double> gmann);
	int modified_hamiltonian_cd(double E, double V, int cp_sign);

        int my_damped_P_matrix_cd(double E, double L, double V, int cp_sign);
        int my_S_matrix_cd(double E, double L, double V, int cp_sign);
        int my_probability_matrix(double P[3][3], int cp_sign, double E, int psteps ,
                                        const double *length, const double *density,
                                        double filter_sigma, void *user_data);
public:

        ProbEngine(vector<double>& p, Baseline &baseline,  void *user_data);
	int get_probability(double P[3][3], int cp_sign, double E);
};

#endif // PROB_ENGINE_H


