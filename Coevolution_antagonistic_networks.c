/*
 gcc -o Coevolution_antagonistic_WoL_networks Coevolution_antagonistic_WoL_networks.c -O3 -lm -Wall -Wextra
 ./Coevolution_antagonistic_WoL_networks net_id flag_temporal
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
// #include <unistd.h>

// #include <process.h> //change getpid --> _getpid

double PI = 3.14159265;

/*######################################################

RAN3 ALGORITHM: RANDOM NUMBER GENERATOR BETWEEN (0,1)

########################################################*/

#define MBIG 1000000000
#define MZ 0
#define FAC (1.0 / MBIG)
int MSEED;
long idum = -1;
// Any large MBIG, and any smaller (but still large) MSEED can be substituted for the above values.

void ini_ran()
{
    // MSEED = (int)(time(NULL));
    MSEED = 1234567890;
}

double Random()
{
    // Set idum to any negative value to initialize or reinitialize the sequence.
    static int inext, inextp;
    static long ma[56]; // The value 56 (range ma[1..55]) is special and should not be modified; see Knuth.
    static int iff = 0;
    long mj, mk;
    int i, ii, k;

    if (idum < 0 || iff == 0)
    { // Initialization
        iff = 1;
        mj = labs(MSEED - labs(idum)); // Initialize ma[55] using the seed idum and the large number MSEED.
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;

        for (i = 1; i <= 54; i++)
        { // Now initialize the rest of the table in a slightly random order with numbers that are not especially random.
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ)
                mk += MBIG;
            mj = ma[ii];
        }
        for (k = 1; k <= 4; k++) // We randomize them by “warming up the generator”
            for (i = 1; i <= 55; i++)
            {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ)
                    ma[i] += MBIG;
            }

        inext = 0;   // Prepare indices for our first generated number.
        inextp = 31; // The constant 31 is special; see Knuth.
        idum = 1;
    }

    // Here is where we start, except on initialization.
    if (++inext == 56)
        inext = 1; // Increment inext and inextp, wrapping around 56 to 1.
    if (++inextp == 56)
        inextp = 1;
    mj = ma[inext] - ma[inextp]; // Generate a new random number subtractively.
    if (mj < MZ)
        mj += MBIG; // Be sure that it is in range.
    ma[inext] = mj; // Store it and output the derived uniform deviate.
    return mj * FAC;
}

/*######################################################

BOX-MULLER ALGORITHM: RANDOM NUMBER GENERATOR WITH NORMAL DISTRIBUTION

########################################################*/

void Box_Muller(double *g1, double *g2)
{
    double d1, d2;
    d1 = sqrt(-2.0 * log(Random()));
    d2 = 2.0 * PI * Random();
    *g1 = -d1 * cos(d2);
    *g2 = -d1 * sin(d2);
}
/*################################################

SPECIES STRUCTS DEFINITIONS, WITH ALL VARIABLES USED

##################################################*/
//

typedef struct
{
    int *w_E;             // information of the neighbours of the HP_networks, size = k_E
    int k_victim;         // number of connections in the network
    double z_trait;       // value of trait
    double theta;         // stabilizing selection, theta = z_trait(0)
    double xi_s;          // intensity of environmental selection
    double xi_d;          // intensity of interaction patterns
    double S_i;           // partial selection by environment
    double *M_ij;         // partial selection by interaction patterns, size = k_E
    double fitness_V;     // fitness of the species
    double fitness_amb_V; // ambiental term of fitness
    double fitness_int_V; // interaction term of fitness
} Victim;                 //

typedef struct
{
    int *w_V;             // information of the neighbours of the HP_networks, size = k_V
    int k_exploiter;      // number of connections in the network
    double z_trait;       // value of trait
    double theta;         // stabilizing selection, theta = z_trait(0)
    double xi_s;          // intensity of environmental selection
    double xi_d;          // intensity of interaction patterns
    double S_i;           // environmental selection
    double *M_ij;         // partial selection by interaction patterns, size = k_V
    double fitness_E;     // fitness of the species
    double fitness_amb_E; // ambiental term of fitness
    double fitness_int_E; // interaction term of fitness
} Exploiter;              //

/**####################

FUNCTION DECLARATION

#######################*/

int ini_species();              // species (both V and E) initialization
int ini_species_distribution(); // initial distribution of species' trait each simulation
void Runge_Kutta(double time);
double interaction_function_exploiters(Exploiter *exploiters_net, Victim *victims_net, int id_exploiter, double time); // interaction of the agents
double interaction_function_victims(Exploiter *exploiters_net, Victim *victims_net, int id_victim, double time);       // interaction of the agents
void measure_fitness();

void write_header_global_temporal(char name_1[100], char name_2[100]); // writes a header with important data in the output files
void write_header_global_final(char name_1[100]);                      // writes a header with important data in the output files
void write_data_time(int simul, double t_iter);                        // writes data about temporal information
void write_data_final(int simul);                                      // writes data about temporal information

/**###########################

GLOBAL VARIABLES DEFINITION

##############################*/

int net_id, flag_temporal; // input variables
double xi_d_V, xi_d_E;     // input variables
double xi_d_V_min, xi_d_V_max, xi_d_V_delta;
double xi_d_E_min, xi_d_E_max, xi_d_E_delta;

int N_Sim;      // number of simulation
int N_Max_Iter; // max iteration of each simulation

int M_Victims;    // number of victims in the HP network
int M_Exploiters; // number of exploiters in the HP network

double h, TMAX; // time step - TMAX
double epsilon; // evolutionary response on victim species
double alpha;   // exploiter preference
double delta;   // integration constant to avoid dividing by 0

FILE *f_out_victims, *f_out_exploiters, *f_out_final, *net_file;

Victim *Victim_Species_List;       // List of species with all the information
Exploiter *Exploiter_Species_List; // List of species with all the information

/**############

MAIN PROGRAM

###############*/

int main(int argc, char **argv)
{

    clock_t begin = clock();
    int i;

    ini_ran();
    Random(); // initialization of the RNG

    /**----------------------------------------------------------

                    CHANGE WHEN NECESSARY

    ----------------------------------------------------------*/

    N_Sim = 100;

    alpha = 0.1;
    TMAX = 300;
    delta = 0.000001;
    h = 0.01;
    N_Max_Iter = (int)TMAX / h;

    xi_d_V_min = 0.1;
    xi_d_V_max = 0.9;
    xi_d_V_delta = 0.1;

    xi_d_E_min = 0.1;
    xi_d_E_max = 0.9;
    xi_d_E_delta = 0.1;

    /**----------------------------------------------------------

    Read the variables used in the simulation:
        - xi_V (selection pressure on victims)
        - xi_E (selection pressure on exploiters)

    ----------------------------------------------------------*/

    double variables[2];
    for (i = 1; i < argc; i++)
    {
        sscanf(argv[i], " %lf", &variables[i - 1]);
    }

    net_id = (int)variables[0];
    flag_temporal = (int)variables[1];

    /**
     * Read the values of M_Victims and M_Exploiters from the network with  net_id
     */

    char net_filename[50];
    // sprintf(net_filename, "HP_networks/H_species_0%d.txt", net_id);
    if (net_id < 10)
    {
        sprintf(net_filename, "VE_networks/A_HP_00%d_V_species.csv", net_id);
    }
    else
    {
        sprintf(net_filename, "VE_networks/A_HP_0%d_V_species.csv", net_id);
    }

    net_file = fopen(net_filename, "r");

    if (net_file != NULL) // reads the positions from a file
    {
        fscanf(net_file, "#M_V=%d", &M_Victims);
    }
    fclose(net_file);

    // sprintf(net_filename, "HP_networks/P_species_0%d.txt", net_id);
    if (net_id < 10)
    {
        sprintf(net_filename, "VE_networks/A_HP_00%d_E_species.csv", net_id);
    }
    else
    {
        sprintf(net_filename, "VE_networks/A_HP_0%d_E_species.csv", net_id);
    }
    net_file = fopen(net_filename, "r");

    if (net_file != NULL) // reads the positions from a file
    {
        fscanf(net_file, "#M_E=%d", &M_Exploiters);
    }
    fclose(net_file);

    /**----------------------------------------------------------

    SPECIES NETWORK INITIALIZATION

    ----------------------------------------------------------*/

    ini_species();

    /**----------------------------------------------------------

    SIMULATION

    ----------------------------------------------------------*/

    if (flag_temporal == 1)
    {
        char filename_victims[200], filename_exploiters[200];
        sprintf(filename_victims, "Files/Networks/Fitness/Temporal/Antagonistic_sameIC_temporal_V_alpha=%.1lf_xi_V=%.3lf_xi_E=%.3lf_0%d.csv", alpha, xi_d_V, xi_d_E, net_id);
        sprintf(filename_exploiters, "Files/Networks/Fitness/Temporal/Antagonistic_sameIC_temporal_E_alpha=%.1lf_xi_V=%.3lf_xi_E=%.3lf_0%d.csv", alpha, xi_d_V, xi_d_E, net_id);
        write_header_global_temporal(filename_victims, filename_exploiters);
    }
    else
    {
        char filename_final[100];
        sprintf(filename_final, "Files/A_HP_%d_final_alpha=%.1lf_alt.csv", net_id, alpha);
        write_header_global_final(filename_final);
    }

    for (xi_d_V = xi_d_V_min; xi_d_V < xi_d_V_max + xi_d_V_delta / 10.; xi_d_V += xi_d_V_delta)
    {
        for (xi_d_E = xi_d_E_min; xi_d_E < xi_d_E_max + xi_d_E_delta / 10.; xi_d_E += xi_d_E_delta)
        {
            int sim, iter;

            for (sim = 0; sim < N_Sim; sim++)
            {

                ini_species_distribution();
                measure_fitness();

                iter = 0;
                if (flag_temporal == 1)
                {
                    write_data_time(sim, iter * h);
                }

                do
                {
                    Runge_Kutta(iter * h);
                    measure_fitness();

                    iter++;
                    if (flag_temporal == 1 && iter % 1 == 0)
                    {
                        write_data_time(sim, iter * h);
                    }

                } while (iter < N_Max_Iter);

                if (flag_temporal == 0)
                {
                    write_data_final(sim);
                }

                if (flag_temporal == 1)
                {
                    fprintf(f_out_victims, "\n");
                    fprintf(f_out_exploiters, "\n");
                }

            }
        }

        fprintf(f_out_final, "\n");
    }

    if (flag_temporal == 1)
    {
        fclose(f_out_victims);
        fclose(f_out_exploiters);
    }

    if (flag_temporal == 0)
    {
        fclose(f_out_final);
    }

    /* ------------------------------------- */

    int i_victims, i_exploiters;

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        free(Victim_Species_List[i_victims].w_E);
        free(Victim_Species_List[i_victims].M_ij);
    }
    free(Victim_Species_List);

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        free(Exploiter_Species_List[i_exploiters].w_V);
        free(Exploiter_Species_List[i_exploiters].M_ij);
    }
    free(Exploiter_Species_List);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("Time of simulation: %.2lf\n", time_spent);

    return 0;
}

/**########################

FUNCTION DEFINITIONS

#########################*/

int ini_species()
{
    /**
     * @brief Initializes the bipartite network structure and allocates memory for species.
     *
     * This function dynamically allocates memory for the global arrays of Victim and 
     * Exploiter species. It reads the network topology (adjacency lists) from external 
     * CSV files corresponding to the specific `net_id`. For each species, it stores 
     * its degree (number of connections), populates its list of neighbors (`w_E` for 
     * victims, `w_V` for exploiters), and allocates memory for the interaction 
     * selection arrays (`M_ij`).
     *
     * @return int Returns 0 upon successful initialization and file reading.
     */

    Victim_Species_List = (Victim *)malloc(M_Victims * sizeof(Victim));
    Exploiter_Species_List = (Exploiter *)malloc(M_Exploiters * sizeof(Exploiter));

    int i_victims, i_exploiters;
    int idx_victims, idx_exploiters;
    int aux;

    char net_filename[50];
    if (net_id < 10)
    {
        sprintf(net_filename, "VE_networks/A_HP_00%d_V_species.csv", net_id);
    }
    else
    {
        sprintf(net_filename, "VE_networks/A_HP_0%d_V_species.csv", net_id);
    }

    net_file = fopen(net_filename, "r");

    if (net_file != NULL) // reads the positions from a file
    {
        char header_string[400];
        fgets(header_string, 400, net_file); // first line
        for (i_victims = 0; i_victims < M_Victims; i_victims++)
        {
            fscanf(net_file, "%d,%d,", &aux, &Victim_Species_List[i_victims].k_victim);
            // Victim_Species_List[i_victims].theta = -1. + 2. * Random();

            Victim_Species_List[i_victims].w_E = (int *)calloc(Victim_Species_List[i_victims].k_victim, sizeof(int));
            for (idx_exploiters = 0; idx_exploiters < Victim_Species_List[i_victims].k_victim; idx_exploiters++)
            {
                fscanf(net_file, "%d,", &Victim_Species_List[i_victims].w_E[idx_exploiters]);
            }

            Victim_Species_List[i_victims].M_ij = (double *)calloc(Victim_Species_List[i_victims].k_victim, sizeof(double));
        }
    }
    fclose(net_file);

    if (net_id < 10)
    {
        sprintf(net_filename, "VE_networks/A_HP_00%d_E_species.csv", net_id);
    }
    else
    {
        sprintf(net_filename, "VE_networks/A_HP_0%d_E_species.csv", net_id);
    }

    net_file = fopen(net_filename, "r");

    if (net_file != NULL) // reads the positions from a file
    {
        char header_string[400];
        fgets(header_string, 400, net_file); // first line
        for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
        {
            fscanf(net_file, "%d,%d,", &aux, &Exploiter_Species_List[i_exploiters].k_exploiter);

            // Exploiter_Species_List[i_exploiters].theta = -1. + 2. * Random();

            Exploiter_Species_List[i_exploiters].w_V = (int *)calloc(Exploiter_Species_List[i_exploiters].k_exploiter, sizeof(int));
            for (idx_victims = 0; idx_victims < Exploiter_Species_List[i_exploiters].k_exploiter; idx_victims++)
            {
                fscanf(net_file, "%d,", &Exploiter_Species_List[i_exploiters].w_V[idx_victims]);
            }

            Exploiter_Species_List[i_exploiters].M_ij = (double *)calloc(Exploiter_Species_List[i_exploiters].k_exploiter, sizeof(double));
        }
    }

    fclose(net_file);

    return 0;
}

int ini_species_distribution()
{
    /**
     * @brief Initializes the trait distributions and selection parameters for all species.
     *
     * This function assigns the initial traits (z_trait) for both Victim 
     * and Exploiter species using a normal (Gaussian) distribution generated via the 
     * Box-Muller transform. It sets the environmental optimum (theta) equal to the 
     * initial trait value for each species. Additionally, it initializes the selection 
     * intensities (xi_d for interaction and xi_s for environment) based on global 
     * variables, and resets the interaction partial selection arrays (M_ij) to zero 
     * before the simulation begins.
     *
     * @return int Returns 0 upon successful initialization.
     */

    double gauss1, gauss2;
    gauss1 = 0.;
    gauss2 = 0.;

    int i_victims, i_exploiters;
    int idx_victims, idx_exploiters;

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        if (i_victims % 2 != 0)
        {
            Victim_Species_List[i_victims].z_trait = 0.5 * gauss1; // + Victim_Species_List[i_victims].theta;
        }
        else
        {
            Box_Muller(&gauss1, &gauss2);

            Victim_Species_List[i_victims].z_trait = 0.5 * gauss2; // + Victim_Species_List[i_victims].theta;
        }

        Victim_Species_List[i_victims].theta = Victim_Species_List[i_victims].z_trait;
        Victim_Species_List[i_victims].S_i = 0.;
        Victim_Species_List[i_victims].xi_d = xi_d_V;
        Victim_Species_List[i_victims].xi_s = 1. - Victim_Species_List[i_victims].xi_d;

        for (idx_exploiters = 0; idx_exploiters < Victim_Species_List[i_victims].k_victim; idx_exploiters++)
        {
            Victim_Species_List[i_victims].M_ij[idx_exploiters] = 0.;
        }
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        if (i_exploiters % 2 != 0)
        {
            Exploiter_Species_List[i_exploiters].z_trait = 0.5 * gauss1; // + Exploiter_Species_List[i_exploiters].theta;
        }
        else
        {
            Box_Muller(&gauss1, &gauss2);

            Exploiter_Species_List[i_exploiters].z_trait = 0.5 * gauss2; // + Exploiter_Species_List[i_exploiters].theta;
        }

        Exploiter_Species_List[i_exploiters].theta = Exploiter_Species_List[i_exploiters].z_trait;
        Exploiter_Species_List[i_exploiters].S_i = 0.;
        Exploiter_Species_List[i_exploiters].xi_d = xi_d_E;
        Exploiter_Species_List[i_exploiters].xi_s = 1. - Exploiter_Species_List[i_exploiters].xi_d;

        for (idx_victims = 0; idx_victims < Exploiter_Species_List[i_exploiters].k_exploiter; idx_victims++)
        {
            Exploiter_Species_List[i_exploiters].M_ij[idx_victims] = 0.;
        }
    }

    return 0;
}

void Runge_Kutta(double time)
{
    /**
     * @brief Executes a single step of the 4th-order Runge-Kutta (RK4) numerical integration.
     *
     * This function calculates the temporal evolution of the traits 
     * (z_trait) for all Victim and Exploiter species. It temporarily allocates 
     * auxiliary networks to compute the four intermediate slopes (k1, k2, k3, k4) 
     * by evaluating the differential equations defined in the interaction functions. 
     * Finally, it updates the global trait values for the next time step using the 
     * weighted average of these slopes and frees the allocated memory.
     *
     * @param time The current time step of the simulation.
     */

    Exploiter *aux_exploiters_net;
    Victim *aux_victims_net;

    aux_exploiters_net = (Exploiter *)malloc(M_Exploiters * sizeof(Exploiter));
    aux_victims_net = (Victim *)malloc(M_Victims * sizeof(Victim));

    int i_exploiters, i_victims;

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        aux_exploiters_net[i_exploiters].z_trait = Exploiter_Species_List[i_exploiters].z_trait;
        aux_exploiters_net[i_exploiters].k_exploiter = Exploiter_Species_List[i_exploiters].k_exploiter;
        aux_exploiters_net[i_exploiters].theta = Exploiter_Species_List[i_exploiters].theta;
        aux_exploiters_net[i_exploiters].xi_d = Exploiter_Species_List[i_exploiters].xi_d;
        aux_exploiters_net[i_exploiters].xi_s = Exploiter_Species_List[i_exploiters].xi_s;
        aux_exploiters_net[i_exploiters].w_V = (int *)malloc(aux_exploiters_net[i_exploiters].k_exploiter * sizeof(int));
        aux_exploiters_net[i_exploiters].M_ij = (double *)malloc(aux_exploiters_net[i_exploiters].k_exploiter * sizeof(double));
        int idx_victims;
        for (idx_victims = 0; idx_victims < aux_exploiters_net[i_exploiters].k_exploiter; idx_victims++)
        {
            aux_exploiters_net[i_exploiters].w_V[idx_victims] = Exploiter_Species_List[i_exploiters].w_V[idx_victims];
        }
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        aux_victims_net[i_victims].z_trait = Victim_Species_List[i_victims].z_trait;
        aux_victims_net[i_victims].k_victim = Victim_Species_List[i_victims].k_victim;
        aux_victims_net[i_victims].theta = Victim_Species_List[i_victims].theta;
        aux_victims_net[i_victims].xi_d = Victim_Species_List[i_victims].xi_d;
        aux_victims_net[i_victims].xi_s = Victim_Species_List[i_victims].xi_s;
        aux_victims_net[i_victims].w_E = (int *)malloc(aux_victims_net[i_victims].k_victim * sizeof(int));
        aux_victims_net[i_victims].M_ij = (double *)malloc(aux_victims_net[i_victims].k_victim * sizeof(double));
        int idx_exploiters;
        for (idx_exploiters = 0; idx_exploiters < aux_victims_net[i_victims].k_victim; idx_exploiters++)
        {
            aux_victims_net[i_victims].w_E[idx_exploiters] = Victim_Species_List[i_victims].w_E[idx_exploiters];
        }
    }

    double *k1_exploiters, *k1_victims;
    double *k2_exploiters, *k2_victims;
    double *k3_exploiters, *k3_victims;
    double *k4_exploiters, *k4_victims;

    k1_exploiters = (double *)malloc(M_Exploiters * sizeof(double));
    k2_exploiters = (double *)malloc(M_Exploiters * sizeof(double));
    k3_exploiters = (double *)malloc(M_Exploiters * sizeof(double));
    k4_exploiters = (double *)malloc(M_Exploiters * sizeof(double));
    k1_victims = (double *)malloc(M_Victims * sizeof(double));
    k2_victims = (double *)malloc(M_Victims * sizeof(double));
    k3_victims = (double *)malloc(M_Victims * sizeof(double));
    k4_victims = (double *)malloc(M_Victims * sizeof(double));

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        k1_exploiters[i_exploiters] = h * (interaction_function_exploiters(aux_exploiters_net, aux_victims_net, i_exploiters, time));
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        k1_victims[i_victims] = h * (interaction_function_victims(aux_exploiters_net, aux_victims_net, i_victims, time));
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        aux_exploiters_net[i_exploiters].z_trait = Exploiter_Species_List[i_exploiters].z_trait + 0.5 * k1_exploiters[i_exploiters];
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        aux_victims_net[i_victims].z_trait = Victim_Species_List[i_victims].z_trait + 0.5 * k1_victims[i_victims];
    }

    // obtain k2 for the new networks

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        k2_exploiters[i_exploiters] = h * (interaction_function_exploiters(aux_exploiters_net, aux_victims_net, i_exploiters, time));
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        k2_victims[i_victims] = h * (interaction_function_victims(aux_exploiters_net, aux_victims_net, i_victims, time));
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        aux_exploiters_net[i_exploiters].z_trait = Exploiter_Species_List[i_exploiters].z_trait + 0.5 * k2_exploiters[i_exploiters];
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        aux_victims_net[i_victims].z_trait = Victim_Species_List[i_victims].z_trait + 0.5 * k2_victims[i_victims];
    }

    // obtain k3 for the new networks

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        k3_exploiters[i_exploiters] = h * (interaction_function_exploiters(aux_exploiters_net, aux_victims_net, i_exploiters, time));
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        k3_victims[i_victims] = h * (interaction_function_victims(aux_exploiters_net, aux_victims_net, i_victims, time));
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        aux_exploiters_net[i_exploiters].z_trait = Exploiter_Species_List[i_exploiters].z_trait + k3_exploiters[i_exploiters];
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        aux_victims_net[i_victims].z_trait = Victim_Species_List[i_victims].z_trait + k3_victims[i_victims];
    }

    // obtain k4 for the new networks

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        k4_exploiters[i_exploiters] = h * (interaction_function_exploiters(aux_exploiters_net, aux_victims_net, i_exploiters, time));
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        k4_victims[i_victims] = h * (interaction_function_victims(aux_exploiters_net, aux_victims_net, i_victims, time));
    }

    // update trait values for h + 1

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        Exploiter_Species_List[i_exploiters].z_trait = Exploiter_Species_List[i_exploiters].z_trait + 1. / 6 * (k1_exploiters[i_exploiters] + 2. * k2_exploiters[i_exploiters] + 2. * k3_exploiters[i_exploiters] + k4_exploiters[i_exploiters]);
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        Victim_Species_List[i_victims].z_trait = Victim_Species_List[i_victims].z_trait + 1. / 6 * (k1_victims[i_victims] + 2. * k2_victims[i_victims] + 2. * k3_victims[i_victims] + k4_victims[i_victims]);
    }

    free(k1_exploiters);
    free(k2_exploiters);
    free(k3_exploiters);
    free(k4_exploiters);
    free(k1_victims);
    free(k2_victims);
    free(k3_victims);
    free(k4_victims);

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        free(aux_exploiters_net[i_exploiters].w_V);
        free(aux_exploiters_net[i_exploiters].M_ij);
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        free(aux_victims_net[i_victims].w_E);
        free(aux_victims_net[i_victims].M_ij);
    }

    free(aux_exploiters_net);
    free(aux_victims_net);
}

double interaction_function_exploiters(Exploiter *exploiters_net, Victim *victims_net, int id_exploiter, double time)
{
    /**
     * @brief Computes the rate of trait change for a specific Exploiter species.
     *
     * This function calculates the derivative (slope) of the trait for an 
     * Exploiter species at a given time step. It evaluates the selective pressures 
     * acting on the species, which consist of two main components:
     * 1. Environmental selection: The pull towards the species' optimal environmental trait (theta).
     * 2. Interaction selection: The pressure exerted by its antagonistic interactions 
     * with connected Victim species. For exploiters, this typically drives their traits 
     * to match those of their victims to maximize exploitation efficiency.
     *
     * The result is used by the Runge-Kutta algorithm to update the species' trait.
     *
     * @param id_exploiter The index of the Exploiter species being evaluated.
     * @param time The current time step of the simulation.
     * @return double The calculated rate of change (derivative) of the Exploiter's trait.
     */

    int idx_victims;

    exploiters_net[id_exploiter].S_i = exploiters_net[id_exploiter].xi_s * (exploiters_net[id_exploiter].theta - exploiters_net[id_exploiter].z_trait);
    double aux_denominator_Ei = 0.;

    for (idx_victims = 0; idx_victims < exploiters_net[id_exploiter].k_exploiter; idx_victims++) // loop for the denominator
    {
        int id_victim;
        id_victim = exploiters_net[id_exploiter].w_V[idx_victims];

        double difference = victims_net[id_victim].z_trait - exploiters_net[id_exploiter].z_trait;

        aux_denominator_Ei += exp(-alpha * difference * difference);
    }

    double aux_selection_sum = 0.;
    for (idx_victims = 0; idx_victims < exploiters_net[id_exploiter].k_exploiter; idx_victims++)
    {
        double aux_p_ij = 0.;
        int id_victim;
        id_victim = exploiters_net[id_exploiter].w_V[idx_victims];

        double difference = victims_net[id_victim].z_trait - exploiters_net[id_exploiter].z_trait + delta;
        double abs_difference = fabs(difference);

        aux_p_ij = exp(-alpha * difference * difference);
        exploiters_net[id_exploiter].M_ij[idx_victims] = aux_p_ij;

        aux_selection_sum += exploiters_net[id_exploiter].M_ij[idx_victims] * difference / abs_difference / (aux_denominator_Ei + delta);

        // if (id_exploiter == 0 && time > 8.72)
        // {
        //     printf("\tt=%.2lf - z[%d]=%.5lf - x[%d]=%.5lf - pij=%.5lf - Mij=%.5lf - diff=%.5lf - |diff|=%.5lf - denom = %.5lf\n", time, id_victim, victims_net[id_victim].z_trait, id_exploiter, exploiters_net[id_exploiter].z_trait, aux_p_ij, exploiters_net[id_exploiter].M_ij[idx_victims], difference, abs_difference, aux_denominator_Ei);
        //     getchar();
        // }
    }

    return exploiters_net[id_exploiter].S_i + exploiters_net[id_exploiter].xi_d * aux_selection_sum; // equation for exploiters
}

double interaction_function_victims(Exploiter *exploiters_net, Victim *victims_net, int id_victim, double time)
{
    /**
     * @brief Computes the rate of trait change for a specific Victims species.
     *
     * This function calculates the derivative (slope) of the trait for an
     * Victim species at a given time step. It evaluates the selective pressures 
     * acting on the species, which consist of two main components:
     * 1. Environmental selection: The pull towards the species' optimal environmental trait (theta).
     * 2. Interaction selection: The pressure exerted by its antagonistic interactions 
     * with connected Exploiters species. Unlike exploiters, victims evolve to evade 
     * their antagonists, so this pressure drives their traits away from the traits 
     * of their exploiters to minimize the interaction rate.
     *
     * The result is used by the Runge-Kutta algorithm to update the species' trait.
     *
     * @param id_victim The index of the Victim species being evaluated.
     * @param time The current time step of the simulation.
     * @return double The calculated rate of change (derivative) of the Exploiter's trait.
     */

    int idx_exploiters;

    victims_net[id_victim].S_i = victims_net[id_victim].xi_s * (victims_net[id_victim].theta - victims_net[id_victim].z_trait);

    double aux_denominator_Vi = 0.;

    for (idx_exploiters = 0; idx_exploiters < victims_net[id_victim].k_victim; idx_exploiters++) // loop for the denominator
    {
        int id_exploiter;
        id_exploiter = victims_net[id_victim].w_E[idx_exploiters];

        double difference = victims_net[id_victim].z_trait - exploiters_net[id_exploiter].z_trait;

        aux_denominator_Vi += exp(-alpha * difference * difference);
    }

    double aux_selection_sum = 0.;

    for (idx_exploiters = 0; idx_exploiters < victims_net[id_victim].k_victim; idx_exploiters++)
    {
        int id_exploiter;
        id_exploiter = victims_net[id_victim].w_E[idx_exploiters];
        double difference = victims_net[id_victim].z_trait - exploiters_net[id_exploiter].z_trait + delta;
        double abs_difference = fabs(difference);

        victims_net[id_victim].M_ij[idx_exploiters] = exp(-alpha * difference * difference); // interaction_selection_E_on_V(id_victim, id_exploiters, exploiters_net, victims_net);
        aux_selection_sum += victims_net[id_victim].M_ij[idx_exploiters] * difference / abs_difference / (aux_denominator_Vi + delta);
    }

    return victims_net[id_victim].S_i + victims_net[id_victim].xi_d * aux_selection_sum;
}

void measure_fitness()
{

    /**
     * @brief Calculates the fitness of all Victim and Exploiter species.
     *
     * This function evaluates the current fitness of each species in the network 
     * based on its phenotypic traits. The total fitness is composed of two main 
     * components:
     * 1. Environmental fitness (fitness_amb): Represents the adaptation of the 
     * species to its optimum environmental trait (theta). It is typically modeled 
     * using a decay function based on the difference between the current trait 
     * and the environmental optimum.
     * 2. Interaction fitness (fitness_int): Represents the evolutionary success or 
     * penalty resulting from antagonistic interactions with connected species 
     * (matching or escaping traits).
     *
     * The computed values are stored in the respective structs (`Victim` and 
     * `Exploiter`) to track the adaptiveness of the species over time.
     */
    
    int i_exploiters, i_victims, idx_victims, idx_exploiters;
    double fitness_E, fitness_amb_E, fitness_int_E, diff_amb_E, diff_sel_E, diff_amb_EV;
    double fitness_V, fitness_amb_V, fitness_int_V, diff_amb_V, diff_sel_V;
    double max_xi_E = xi_d_E / (1. - xi_d_E);
    double max_xi_V = xi_d_V / (1. - xi_d_V);

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        fitness_E = 0.;
        diff_amb_E = fabs(Exploiter_Species_List[i_exploiters].z_trait - Exploiter_Species_List[i_exploiters].theta);

        fitness_amb_E = (1. - diff_amb_E / max_xi_E);
        int k_E = Exploiter_Species_List[i_exploiters].k_exploiter;
        fitness_int_E = 0.;
        for (idx_victims = 0; idx_victims < k_E; idx_victims++)
        {
            int id_victim;
            id_victim = Exploiter_Species_List[i_exploiters].w_V[idx_victims];
            diff_amb_EV = fabs(Victim_Species_List[id_victim].theta - Exploiter_Species_List[i_exploiters].theta);
            diff_sel_E = fabs(Victim_Species_List[id_victim].z_trait - Exploiter_Species_List[i_exploiters].z_trait);
            fitness_int_E += 1. / k_E * (1. - diff_sel_E / (max_xi_E + max_xi_V + diff_amb_EV));
        }
        Exploiter_Species_List[i_exploiters].fitness_amb_E = fitness_amb_E;
        Exploiter_Species_List[i_exploiters].fitness_int_E = fitness_int_E;
        Exploiter_Species_List[i_exploiters].fitness_E = (1. - xi_d_E) * fitness_amb_E + xi_d_E * fitness_int_E;
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        fitness_V = 0.;
        diff_amb_V = fabs(Victim_Species_List[i_victims].z_trait - Victim_Species_List[i_victims].theta);

        fitness_amb_V = (1. - diff_amb_V / max_xi_V);
        int k_V = Victim_Species_List[i_victims].k_victim;
        fitness_int_V = 0.;
        for (idx_exploiters = 0; idx_exploiters < k_V; idx_exploiters++)
        {
            int id_exploiter;
            id_exploiter = Victim_Species_List[i_victims].w_E[idx_exploiters];
            diff_amb_EV = fabs(Victim_Species_List[i_victims].theta - Exploiter_Species_List[id_exploiter].theta);
            diff_sel_V = fabs(Victim_Species_List[i_victims].z_trait - Exploiter_Species_List[id_exploiter].z_trait);
            fitness_int_V += 1. / k_V * (diff_sel_V / (max_xi_E + max_xi_V + diff_amb_EV));
        }
        Victim_Species_List[i_victims].fitness_amb_V = fitness_amb_V;
        Victim_Species_List[i_victims].fitness_int_V = fitness_int_V;
        Victim_Species_List[i_victims].fitness_V = (1. - xi_d_V) * fitness_amb_V + xi_d_V * fitness_int_V;
    }
}

/**########################

TEXTFILE FUNCTIONS

#########################*/

void write_header_global_temporal(char name_1[100], char name_2[100])
{
    /*
        File with crossovers times for each simulation
    */

    int i_exploiters, i_victims;

    time_t curtime;
    time(&curtime);

    f_out_victims = fopen(name_1, "w");

    fprintf(f_out_victims, "#Date of simulation: %s", ctime(&curtime));
    fprintf(f_out_victims, "#Seed: %d\n#\n", MSEED);
    fprintf(f_out_victims, "#Number of simulations: %d\n#\n", N_Sim);
    fprintf(f_out_victims, "#----------------------------------------------------------\n");
    fprintf(f_out_victims, "#net_id,xi_V,xi_E,sim,time");
    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        fprintf(f_out_victims, ",theta_V_%d,z_V_%d,fitness_V_%d", i_victims, i_victims, i_victims);
    }
    fprintf(f_out_victims, "\n");

    f_out_exploiters = fopen(name_2, "w");

    fprintf(f_out_exploiters, "#Date of simulation: %s", ctime(&curtime));
    fprintf(f_out_exploiters, "#Seed: %d\n#\n", MSEED);
    fprintf(f_out_exploiters, "#Number of simulations: %d\n#\n", N_Sim);
    fprintf(f_out_exploiters, "#----------------------------------------------------------\n");
    fprintf(f_out_exploiters, "#net_id,xi_V,xi_E,sim,time");
    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        fprintf(f_out_exploiters, ",theta_E_%d,z_E_%d,fitness_E_%d", i_exploiters, i_exploiters, i_exploiters);
    }
    fprintf(f_out_exploiters, "\n");
}

void write_data_time(int simul, double t_iter)
{
    int i_exploiters, i_victims;

    fprintf(f_out_victims, "0%d,%.3lf,%.3lf,%d,%.5lf", net_id, xi_d_E, xi_d_V, simul, t_iter);
    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        fprintf(f_out_victims, ",%.5lf,%.5lf,%.5lf", Victim_Species_List[i_victims].theta, Victim_Species_List[i_victims].z_trait, Victim_Species_List[i_victims].fitness_V);
    }
    fprintf(f_out_victims, "\n");

    fprintf(f_out_exploiters, "0%d,%.3lf,%.3lf,%d,%.5lf", net_id, xi_d_E, xi_d_V, simul, t_iter);
    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        fprintf(f_out_exploiters, ",%.5lf,%.5lf,%.5lf", Exploiter_Species_List[i_exploiters].theta, Exploiter_Species_List[i_exploiters].z_trait, Exploiter_Species_List[i_exploiters].fitness_E);
    }
    fprintf(f_out_exploiters, "\n");
}

void write_header_global_final(char name_1[100])
{
    /*
        File with crossovers times for each simulation
    */

    int i_exploiters, i_victims;

    time_t curtime;
    time(&curtime);

    f_out_final = fopen(name_1, "w");

    fprintf(f_out_final, "#Date of simulation: %s", ctime(&curtime));
    fprintf(f_out_final, "#Seed: %d\n#\n", MSEED);
    fprintf(f_out_final, "#Number of simulations: %d\n#\n", N_Sim);
    fprintf(f_out_final, "#----------------------------------------------------------\n");
    fprintf(f_out_final, "#net_id,xi_E,xi_V,sim");

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        fprintf(f_out_final, ",theta_V_%d,z_V_%d_tinf,fitness_V_%d,fitness_amb_V_%d,fitness_int_V_%d", i_victims, i_victims, i_victims, i_victims, i_victims);
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        fprintf(f_out_final, ",theta_E_%d,z_E_%d_tinf,fitness_E_%d,fitness_amb_E_%d,fitness_int_E_%d", i_exploiters, i_exploiters, i_exploiters, i_exploiters, i_exploiters);
    }

    fprintf(f_out_final, "\n");
}

void write_data_final(int simul)
{
    int i_exploiters, i_victims;
    fprintf(f_out_final, "0%d,%.3lf,%.3lf,%d", net_id, xi_d_E, xi_d_V, simul);
    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        fprintf(f_out_final, ",%.5lf,%.5lf,%.5lf,%.5lf,%.5lf", Victim_Species_List[i_victims].theta, Victim_Species_List[i_victims].z_trait, 
            Victim_Species_List[i_victims].fitness_V, Victim_Species_List[i_victims].fitness_amb_V, Victim_Species_List[i_victims].fitness_int_V);
    }
    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        fprintf(f_out_final, ",%.5lf,%.5lf,%.5lf,%.5lf,%.5lf", Exploiter_Species_List[i_exploiters].theta, Exploiter_Species_List[i_exploiters].z_trait,
            Exploiter_Species_List[i_exploiters].fitness_E, Exploiter_Species_List[i_exploiters].fitness_amb_E, Exploiter_Species_List[i_exploiters].fitness_int_E);
    }
    fprintf(f_out_final, "\n");
}