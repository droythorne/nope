#include "./prob_engine.h"
#include <iostream>
#include <fstream>
#include <globes/globes.h>
double true_theta12;
double true_theta13;
double true_theta23;
double true_deltacp;
double true_sdm;
double true_ldm; /* Note: sdm/ldm measured in ev^2 */
double true_damping;

void init()
{
        //Initialize variables
        true_theta12 = asin(sqrt(0.83))/2;
        true_theta13 = asin(sqrt(0.0))/2;
        true_theta23 = M_PI/4;
        true_deltacp = M_PI/4  ;
        true_sdm = 8.2e-5;
        true_ldm = 2.5e-3; /* Note: sdm/ldm measured in ev^2 */
       true_damping = 0.0;
}
int main(int argc, char *argv[]) {
	init();
	vector<double> params(23,0);
	params[GLB_THETA_12] = true_theta12;
	params[GLB_THETA_13] = true_theta13;
	params[GLB_THETA_23] = true_theta23;
	params[GLB_DELTA_CP] = true_deltacp;
	params[GLB_DM_21] = true_sdm;
	params[GLB_DM_31] = true_ldm;
	params[GLB_DAMPING_EXPONENT] = true_damping;

  std::ofstream outp("transition_prob.dat", std::ios_base::app);//Will append to existing file if already present.

                if (! outp){ //opened failed
                        std::cerr << "cannot open c_JIFS_lep.dat for output\n";
                        exit(-1);
                }

	glbInit(argv[0]);
	glbInitExperiment("./NFstandard.glb", &glb_experiment_list[0], &glb_num_of_exps);
	glb_params true_values = glbAllocParams();
	glbDefineParams(true_values, true_theta12, true_theta13, true_theta23, true_deltacp, true_sdm, true_ldm);
	glbSetDensityParams(true_values, 1.0,GLB_ALL);
	glbSetOscillationParameters(true_values);
	glbSetRates();

	double P_matrix[3][3];
	double my_step[10];
        double my_density[10];
        my_step[0] = 100;
        my_step[1] = 500;
        my_density[0]=0.;
        my_density[1]=0.;
	Baseline  my_baseline(1, my_density, my_step);

	ProbEngine prob_gen1(params, my_baseline, 0);
	for(int i = 0; i < 200 ; i++) {
		my_step[0] = 100 + i*50;

		int check = prob_gen1.get_probability(P_matrix, 1, 23.);

		outp << my_step[0] << " " <<  P_matrix[0][1] << " "  <<  glbVacuumProbability(1,2,1,23.,100+i*50) << std::endl;
	}

/*
	 std::ofstream outp2("transition_prob_glb.dat", std::ios_base::app);//Will append to existing file if already present.

                if (! outp2){ //opened failed
                        std::cerr << "cannot open c_JIFS_lep.dat for output\n";
                        exit(-1);
                }

	for(int i=0; i<200; i++) {

		outp2 << glbVacuumProbability(1,2,1,23.,100+i*50) << std::endl;
	}
*/
	return 0;
}
