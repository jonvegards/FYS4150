/* Program for simulating a 2D lattice
 * with the Ising model by using the
 * Metropolis algorithm
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include "lib.h"
#include "time.h"

using namespace std;

/* 1. Create a state, compute the difference
 *      - a matrix rep. one state
 * 2. Flip one spin, compute the energy
 *      - choose a random site in the lattice by using Monte Carlo
 * 3. Compute the energy difference
 * 4. If E_delta <= 0, go to 7.
 * 5. Calc. w = exp(-beta*E_delta)
 * 6. Compare w and random number r, if r<=w, keep the new state
 * 7. Update various expectation values
 * 8. Repeat 2-7 a given number of times
 *
 * b)   initial_temperature = 1.0
 *      L = 2 (i.e. four spins)
 */

// ensuring boundary conditions by
// "coupling togheter" the elements
// on the edges
int periodic(int i, int dimension, int add) {
    return (i+dimension+add) % (dimension);
    // The add-variable is either 1 or -1 depending on which side
    // of the element we are on, e.g.
    // state_matrix[MC_row][periodic(MC_column, L_spins, -1)] is
    // the element below the one we're looking at
}

void metropolis_algo(int **state_matrix, double L_spins, double &energy, double &magnetization, double *w, long &idum);
void SaveResultsFromSeveralMCCycles(double *, double *, double *, double *, int no_MC, double *);
void CreateStartingState(int **, double, double , double , long &, bool);


int main(){
    int**state_matrix;
    double L_spins;
    long idum;
    double w[17]; // Array containing all possible energies for the system kind of sth what
    double energy=0, magnetization=0,initial_temperature=0;
    double values[5]; // Array for keeping our results

    // Initial values
    //double MC_cycles           = 10000;
    int no_MC = 9;
    double MC_cycles[no_MC];
    double e_variance_aka_heat_capacity[no_MC];
    double m_variance_aka_susceptibility[no_MC];
    double mean_magnetization[no_MC];
    MC_cycles[0] = 10;
    for(int i=1; i<no_MC;i++) MC_cycles[i] = MC_cycles[i-1]*5;

    L_spins             = 20.;
    initial_temperature = 1.0;
    idum                = -1;

    // Initialising w
    for(int delta_E=-8;delta_E<=8;delta_E++) w[delta_E+8] = exp(-delta_E/initial_temperature);

    // Allocating memory for state matrix
    state_matrix = (int **) matrix(L_spins, L_spins, sizeof(L_spins));
    // Filling state matrix with either a random state or with every spin up
    bool random = true;
    CreateStartingState(state_matrix, L_spins, energy, magnetization, idum, random);

    clock_t start, finish;
    double time_used[no_MC];

    for(int q=0;q<no_MC;q++){
        start = clock();
        // Setting values to zero
        for(int i=0;i<5;i++) values[i] = 0;
        // Using metropolis algorithm in order to reach an equilibrium
        for(int n=0;n<MC_cycles[q]; n++){
            metropolis_algo(state_matrix, L_spins, energy, magnetization, w, idum);
            values[0] += energy;
            values[1] += energy*energy;
            values[2] += magnetization;
            values[3] += magnetization*magnetization;
            values[4] += fabs(magnetization);
        }
        finish = clock();

        time_used[q] = ((double) (finish-start)/CLOCKS_PER_SEC);

        // Print results oui oui
        double normalization = MC_cycles[q];
        e_variance_aka_heat_capacity[q] = (values[1]/normalization - values[0]*values[0]/(normalization*normalization))/(L_spins*L_spins);
        m_variance_aka_susceptibility[q] = (values[3]/normalization - values[2]*values[2]/(normalization*normalization))/(L_spins*L_spins);
        mean_magnetization[q] = values[4]/(normalization*L_spins*L_spins);
    }
    SaveResultsFromSeveralMCCycles(e_variance_aka_heat_capacity,m_variance_aka_susceptibility, mean_magnetization ,MC_cycles, no_MC, time_used);
    free_matrix(((void **) state_matrix)); // free memory
    return 0;
}

void metropolis_algo(int **state_matrix, double L_spins, double &energy, double &magnetization, double *w, long &idum){
    for(int i=0;i<L_spins*L_spins;i++){
        // Picking random site in lattice
        int MC_row =    (int) (ran1(&idum)*L_spins);
        int MC_column = (int) (ran1(&idum)*L_spins);
        // calculating energy difference if one spin is flipped
        int delta_E = 2*state_matrix[MC_row][MC_column]*
                (state_matrix[MC_row][periodic(MC_column, L_spins, -1)] +
                state_matrix[periodic(MC_row, L_spins, -1)][MC_column] +
                state_matrix[MC_row][periodic(MC_column, L_spins, 1)] +
                state_matrix[periodic(MC_row, L_spins, 1)][MC_column]);
        // comparing random number with energy difference
        if( delta_E <= 0 or ran1(&idum) <= w[delta_E + 8]){
            state_matrix[MC_row][MC_column] *= -1;
            energy += (double) delta_E;
            magnetization += (double) 2*state_matrix[MC_row][MC_column];
        }
    }
}


void CreateStartingState(int **state_matrix, double L_spins, double energy, double magnetization, long &idum, bool random){
    if(random){
        cout << "random" << endl;
        // Filling state matrix with random spins
        for(int row=0;row<L_spins;row++){
            for(int column=0;column<L_spins;column++){
                if(ran1(&idum) <= 0.5) state_matrix[row][column] = 1.;

            }
        }
    }

    else{
        // Filling state matrix with spin one states
        for(int row=0;row<L_spins;row++){
            for(int column=0;column<L_spins;column++){
                state_matrix[row][column] = 1.;

            }
        }
    }
    // Calculating energy and magnetization of initial state
    for(int row=0;row<L_spins;row++){
        for(int column=0;column<L_spins;column++){
            energy -=  (double) state_matrix[row][column]*(state_matrix[periodic(row,L_spins,-1)][column] +
                    state_matrix[row][periodic(column, L_spins, -1)]);
            magnetization += (double) state_matrix[row][column];
        }
    }

}

void SaveResultsFromSeveralMCCycles(double *e_var, double *m_var, double *mean_mag, double *MC, int no_MC, double *time_used){
    ofstream myfile;
    myfile.open("RunningWithSeveralMC_random_start.m");
    myfile << "MC " << "= [";
    for (int i=0; i<no_MC; i++){
        myfile << MC[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "e_var " << "= [";
    for (int i=0; i<no_MC; i++){
        myfile << e_var[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "m_var " << "= [";
    for (int i=0; i<no_MC; i++){
        myfile << m_var[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "mean_mag " << "= [";
    for (int i=0; i<no_MC; i++){
        myfile << mean_mag[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "time_used " << "= [";
    for (int i=0; i<no_MC; i++){
        myfile << time_used[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "semilogx(MC, e_var)" << endl;

    myfile.close();
}
