
//  STDP_functions.hpp
//  Created by Jen Crodelle on 6/10/2019
// holds functions to be called in the main file
//
//

#ifndef STDP_functions_hpp
#define STDP_functions_hpp

#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <sstream>
#include <random>
#include <iomanip>
#include <new>
#include <array>

using namespace std;

double t = 0.0;
// define struct to hold information for each neuron
struct Network
{
    double v, avgVolt, LTD, LTD_cort, backgroundF;
    double o1,r1; // post-synaptic tracer
    double o2, o2_past;
    double iSTDP_trace; //tracer for iSTDP
    int type; // 1 - inhib, 0 - exc
    bool Espike, Ispike; //T/F for neuron spiking at timestep
    
    vector <double> tspN; //spike times of this neuron
    array <double, 1000> synapStrength; //strength of synapse from an input synapse to the cortical cell
    vector < double > tsp_input[1000];  // spike times of the input neuron
    vector < int > synapticConn; //will hold list of LGN SYNAPSES that the neuron receives spikes from

    vector <int> GJConn; //holds GJ connections
    
    double g_excite0, g_inhib0; // conductance
    double g_excite1, g_inhib1; // conductance
    int sisterID;   // sister group number

    vector <int> corticalConn; // hold list of CORTICAL cells that this neuron SENDS spikes to
    array <double, 400> corticalStrength;   // holds the strength with which this neuron SENDS spikes
    
    double countInputSpikes; //count for the number of spikes received by this neuron from all LGN synapses
    double countExternalSpikes; // count for number of spikes received by external drive
    double countExcCorticalSpikes, countInhibCorticalSpikes; //count for number of spikes received by this neuron from other cortical cells
};

// global parameters
double oldv, gE0, gE1,gI0, gI1, t_interval;
double a0,a1,b0,b1,k1,k2,vtilda,k1tilda,k2tilda,tspike;

double vT = -45; //threshold voltage
double vR = -60; //rest voltage
double vReset = -60.0; //reset voltage

double vE =  0.0; //exc reversal potential
double vI = -80.0; //inhib reversal potential
double gL = 1.0;   // leak conductance
double overC = 0.05; // membrane time constant

double backgroundF = 0.02; //background strength
double R0 = 5; //external rate, R0
int N_input = 1000; // number of incoming neurons from LGN

// params for connectivity
double sigma_E = 11.0; //ms
double sigma_I = 15.0; //ms
double oversigmaE = 1.0/sigma_E;
double oversigmaI = 1.0/sigma_I;

// params for learning
double input_A_tripLTP = 0.005; // triplet LTP
double A_trip_LTP_cort = 0.015; // triplet LTP
double A_iSTDP = 0.008; // iSTDP weight
double gmax = 0.02;  // max weight strength
double targetRate = 8.0; // target rate for changing LTD
double tau_LTD = 33.7;   //ms
double tau_LTP = 16.8;  //ms
double tau_inhib = 20.0; //ms
double trip_tau_LTD = 114.0; // ms
double backgroundNu = 0.5;     // in Hz
double probSynapseConn = 0.25;  //
double tau_v = 1000.0; // ms
double overTau_v = 1.0/tau_v;

// for stimulating LGN:
double r = rand();
double u = r/(double)RAND_MAX;
double newXa = (999.0*u) + 1; // initialize this global variable
int numInhib;
// tracer vars
vector <double> r1(1000,0.0);
vector <double> r2_past(1000,0.0);
vector <double> r2(1000,0.0);

// files for data to be exported:
ofstream finalWeights("LGNweights_final.csv"); //final weights from each LGN synapse to the cortical neurons
ofstream W("weightMatrix_synapses.csv"); // file of LGN weights over time
ofstream W_cortical("weightMatrix_cortex.csv"); // file for Cortical weights over time
ofstream spTimes("neuronSpTimes.csv"); //file for spike times of network
ofstream parameterFile("parameters.txt"); // parameters used in these simulations
ofstream electConn_file("electricConnections.csv"); // electric connections

// function prototypes
// set up the neuron structs
void create_neurons(Network* Neuron, int N, int T);
// connect GJs with sister cells
void connectGJs_sisters(Network* Neuron, int N, double probGJs);
//connect neurons with synapses using a probability (no structure) and set all initial connections to 0
void setUpZeroSynapses_allToAll(Network* Neuron, int N,  double probCortConnect);
// connect neurons synaptically within some radius (perdiodic in up/down and left/right)
void setUpZeroSynapses_radius(Network* Neuron, int N, int radiusConn[2][2]);
// give synaptic connections a value.
void connectSynapses(Network* Neuron, int N, double strengthConnect[2][2], double gmax_cortical, int flagForRadius);
// calculate the input from LGN cells and the resulting plasticity (triplet)
void inputFromLGN(Network* Neuron, int N, double t, double dt, double flagNewXa);
// calculate the background drive
void poissonBackground(Network* Neuron, int N, double t, double dt, double oversigmaE, double backgroundNu);
// calculate input from other cotical cells and update plasticity (triplet and iSTDP)
void spikeUpdate(Network* Neuron, int N, double t, double dt, double A_trip_LTP, bool cortexLearnFlag, double gmax_cortical, double spikeletSize);

void create_neurons(Network* Neuron, int N, int T)
{
    double pr1,pr2, startValue, countInhib;
    countInhib = 0;
    for (int i = 0; i < N; i++)
    {
        Network n;
        pr2 = (double)rand()/RAND_MAX;
        if (pr2 < 0.2)
        {
            n.type = 1;
            n.backgroundF = backgroundF;//2.0*backgroundF;
            countInhib = countInhib+1;
            for (int L = 0; L<N_input; L++)
            {
                startValue = (0.0 + 0.14*((double)rand()/RAND_MAX))*gmax; //random start value for synaptic connection
                n.synapStrength[L] = startValue;
            }
            n.sisterID = 0;
        }
        else
        {
            n.type = 0;
            n.backgroundF = backgroundF;
            n.sisterID = rand()%6+1;
            for (int L = 0; L<N_input; L++)
            {
                pr1 = (double)rand()/RAND_MAX;
                if (pr1 < probSynapseConn)
                {
                    startValue = (0.3 + 0.2*((double)rand()/RAND_MAX))*gmax; //random start value for synaptic connection
                    n.synapticConn.push_back(L);
                    n.synapStrength[L] = startValue;
                }
                else
                {  n.synapStrength[L] = 0.0; }
            }
        }
        
        //initialize parameters
        n.v = vR;
        n.o1 = 0.0;
        n.o2 = 0.0;
        n.r1 = 0.0;
        n.o2_past = 0.0;
        n.iSTDP_trace = 0.0;
        n.g_excite0 = 0.0;
        n.g_excite1 = 0.0;
        n.g_inhib0 = 0.0;
        n.g_inhib1 = 0.0;
        n.avgVolt = 0.0;
        n.LTD = 0.0;
        n.LTD_cort = 0.0;
        n.Espike = false;
        n.Ispike = false;
        n.countInputSpikes = 0.0;
        n.countExternalSpikes = 0.0;
        n.countExcCorticalSpikes = 0.0;
        n.countInhibCorticalSpikes = 0.0;
        Neuron[i] = n;
    }
    cout << "there are " << countInhib << " inhibitory and " << N-countInhib << " excitatory neurons " << endl;
    numInhib = countInhib;
}


void connectGJs_sisters(Network* Neuron, int N, double probGJs)
{
    cout << "connect sister cells by GJs" << endl;
    double pr;
    int countForGJ = 0;
    for (int i =0; i<N; i++)
    {   if (Neuron[i].type == 0)
        {
            for (int j =0; j<N; j++)
            { if((Neuron[i].sisterID == Neuron[j].sisterID) && (i != j))
                { pr = (double)rand()/RAND_MAX;
                if (pr < probGJs)
                {
                    Neuron[i].GJConn.push_back(j);
                    Neuron[j].GJConn.push_back(i);
                    countForGJ = countForGJ+1;
                }}
            }
        }
    }
    //gets rid of duplicate connections
    for(int j=0; j<N; j++)
    {
            sort(Neuron[j].GJConn.begin(), Neuron[j].GJConn.end());
            Neuron[j].GJConn.erase(unique(Neuron[j].GJConn.begin(), Neuron[j].GJConn.end()), Neuron[j].GJConn.end());
    }
    cout << "There are " << countForGJ << " many GJs " << endl;
}

void setUpZeroSynapses_radius(Network* Neuron, int N, int radiusConn[2][2])
{
    int X1, Y1, X2, Y2, radiusSynConn;
    int M = sqrt(N);
    for (int i=0; i<N; i++)
    {
        X1 = i%M + 1;
        Y1 = floor(i/M) + 1;
        for (int j = 0; j<N; j++)
        {
        if (i != j) // don't connect to yourself and only to those within radius
        {
            X2 = j%M + 1;
            Y2 = floor(j/M) + 1;
            Neuron[i].corticalStrength[j] = 0.0;
            radiusSynConn = radiusConn[Neuron[i].type][Neuron[j].type];
            if ((abs(X1-X2) < radiusSynConn) && (abs(Y1-Y2) < radiusSynConn))
            { Neuron[i].corticalConn.push_back(j); }
            else if ((abs(X1-X2) < radiusSynConn)  && (abs(Y1-Y2) > M-radiusSynConn))
            {Neuron[i].corticalConn.push_back(j);  }
            else if ((abs(Y1-Y2) < radiusSynConn) && (abs(X1-X2) > M-radiusSynConn ))
            {  Neuron[i].corticalConn.push_back(j); }
            else if ((abs(X1-X2) > M-radiusSynConn) && (abs(Y1-Y2) > M-radiusSynConn ))
            {  Neuron[i].corticalConn.push_back(j); }
        }
        else
        {Neuron[i].corticalStrength[j] = 0.0;}
        }
    }
}


void setUpZeroSynapses_allToAll(Network* Neuron, int N,  double probCortConnect)
{
    for (int i=0; i<N; i++)
    {
        for (int j = 0; j<N; j++)
        {
            double pr2 = (double)rand()/RAND_MAX;
            if (i != j)
            {
                if (pr2 < probCortConnect)
                { Neuron[i].corticalConn.push_back(j); }
            }
            Neuron[i].corticalStrength[j] = 0.0;
        }
    }
}


void connectSynapses(Network* Neuron, int N, double strengthConnect[2][2], double gmax_cortical,int flagForRadius)
{
    int numConn, n1;
    double r1, val;
    for (int i=0; i<N; i++)
    {
        r1 = (double)rand()/RAND_MAX;
        numConn = Neuron[i].corticalConn.size();
        for (int j = 0; j<numConn; j++)
        {
            n1 = Neuron[i].corticalConn[j];
            val = strengthConnect[Neuron[i].type][Neuron[n1].type];
            if ((Neuron[i].type == 0) && (Neuron[n1].type == 0) && (flagForRadius==1))
            {
                val = (0.25 + 0.1*((double)rand()/RAND_MAX))*gmax_cortical;
                Neuron[i].corticalStrength[n1] = val;
            }
            else{Neuron[i].corticalStrength[n1] = val;}
        }
    }
}


void inputFromLGN(Network* Neuron ,int N, double t, double dt, double flagNewXa)
{
    double r, u, newu, newr2, tsp, FR, synStrength,synStrength1, updateStrength;
    double R1 = 20.0;
    double sigmaRate = 80.0;
    int lengthTsp;
    
    if (flagNewXa == 1)
    {
        r = (double) rand();
        u = r/(double)RAND_MAX;
        ::newXa = (999.0*u) + 1;
        newr2 = (double)rand();
        newu = newr2/(double)RAND_MAX;
        ::t_interval = t-log(newu)*20.0;
        flagNewXa = 0;
    }
    
    for(int j=0; j<N_input; j++)
    {
        tsp = t;
        while (tsp <= t+dt)
        {
            r = rand();
            u = r/RAND_MAX;
            while (u == 0) //make sure it's not 0, bc then log(0) is inifinite
            {   r = rand();
                u = r/RAND_MAX;
            }
            FR = R0 + R1*(exp(-((::newXa-j)*(::newXa-j))/(2*sigmaRate*sigmaRate)) + exp(-((::newXa+N_input-j)*(::newXa+N_input-j))/(2*sigmaRate*sigmaRate))+exp(-((::newXa-N_input-j)*(::newXa-N_input-j))/(2*sigmaRate*sigmaRate)));
            FR = FR/1000.0; // convert to per ms
            tsp = -log(u)/FR + tsp;
            if (tsp <= t+dt)
            {
                r1[j] = r1[j] + 1.0;
                r2[j] = r2[j] + 1.0;
                for (int KK = 0; KK<N; KK++)
                {
                    // update conductance of inhibitory neurons:
                    if (Neuron[KK].type==1)
                    { synStrength1 = Neuron[KK].synapStrength[j];
                    Neuron[KK].g_excite1 = Neuron[KK].g_excite1 + synStrength1*exp(-(t - tsp)*oversigmaE);}
                    // if neuron KK receives input fomr LGN cell j
                    if (find(Neuron[KK].synapticConn.begin(), Neuron[KK].synapticConn.end(), j) != Neuron[KK].synapticConn.end())
                    {
                        Neuron[KK].tsp_input[j].push_back(tsp);
                        Neuron[KK].countInputSpikes = Neuron[KK].countInputSpikes + 1;
                        synStrength = Neuron[KK].synapStrength[j];
                        Neuron[KK].g_excite1 = Neuron[KK].g_excite1 + synStrength*exp(-(t - tsp)*oversigmaE);
                        
                        //if this list gets too big, take one out (for memory issues)
                        if (Neuron[KK].tsp_input[j].size() > 10.0)
                        { Neuron[KK].tsp_input[j].erase(Neuron[KK].tsp_input[j].begin());}
            
                        lengthTsp = Neuron[KK].tspN.size();
                        if (lengthTsp > 0)
                        {   //LTD since neuron spiked before LGN
                            updateStrength = -Neuron[KK].o1*Neuron[KK].LTD;
                            Neuron[KK].synapStrength[j] = Neuron[KK].synapStrength[j] + updateStrength*gmax;
                        }
                        if (Neuron[KK].synapStrength[j] < 0.0)
                        {Neuron[KK].synapStrength[j] = 0.0;}
                    }
                }
            }
        }
    }
}

void spikeUpdate(Network* Neuron, int N, double t, double dt, double A_trip_LTP, bool cortexLearnFlag, double gmax_cortical, double spikeletSize)
{
    int numGJ, synConN, synConNI, GJcell, foo_synConN, indexL,numLGNconnect,numSynConn, numSynConn_elect;
    double lengthTsp, updateStrength,updateStrength2,updateStrength_iSTDP, Tsp;
    for (int l = 0; l<N; l++)
    {
        // if neuron l spiked in this time step
        if (Neuron[l].Espike)
        {
            // update this neuron's tracer:
            Neuron[l].o1 = Neuron[l].o1 + 1.0;
            Neuron[l].o2 = Neuron[l].o2 + 1.0;
            Neuron[l].r1 =  Neuron[l].r1 + 1.0;
            Neuron[l].iSTDP_trace =  Neuron[l].iSTDP_trace + 1.0;
            
            // send spike to post-synaptic cells:
            Tsp = Neuron[l].tspN.back();
            numSynConn = Neuron[l].corticalConn.size();

            if (Neuron[l].type == 0)
            {
                Neuron[l].countExcCorticalSpikes = Neuron[l].countExcCorticalSpikes + 1.0;
                for (int ii=0; ii<numSynConn; ii++)
                {
                    synConN = Neuron[l].corticalConn[ii];
                    Neuron[synConN].g_excite1 = Neuron[synConN].g_excite1 + Neuron[synConN].corticalStrength[l]*exp(-(t - Tsp)*oversigmaE);
                    if (cortexLearnFlag)
                    {
                        if (Neuron[synConN].type == 0)
                        {
                            // LTD
                            updateStrength = -Neuron[synConN].o1*Neuron[synConN].LTD_cort;
                            Neuron[synConN].corticalStrength[l] = Neuron[synConN].corticalStrength[l] + updateStrength*gmax_cortical;
                            
                            // LTP
                            updateStrength = Neuron[l].o2_past*Neuron[synConN].r1*A_trip_LTP_cort;
                            Neuron[l].corticalStrength[synConN] = Neuron[l].corticalStrength[synConN] + updateStrength*gmax_cortical;
                            if (Neuron[synConN].corticalStrength[l] < 0.0)
                            {Neuron[synConN].corticalStrength[l] = 0.0;}
                            if (Neuron[l].corticalStrength[synConN] > gmax_cortical)
                            {Neuron[l].corticalStrength[synConN] = gmax_cortical;}
                        }
                        else
                        {
                            updateStrength_iSTDP = A_iSTDP*Neuron[synConN].iSTDP_trace;
                            Neuron[l].corticalStrength[synConN] = Neuron[l].corticalStrength[synConN] + updateStrength_iSTDP*gmax_cortical;
                            
                            if (Neuron[l].corticalStrength[synConN] > 2.0*gmax_cortical)
                            {  Neuron[l].corticalStrength[synConN] = 2.0*gmax_cortical;}
                        }
                    }
                }
            }
            // if neuron that spiked is inhib
            else
            {
                Neuron[l].countInhibCorticalSpikes = Neuron[l].countInhibCorticalSpikes + 1.0;
                for (int newI=0; newI<numSynConn; newI++)
                {
                    synConNI = Neuron[l].corticalConn[newI];
                    Neuron[synConNI].g_inhib1 = Neuron[synConNI].g_inhib1 + Neuron[synConNI].corticalStrength[l]*exp(-(t - Tsp)*oversigmaI);
                    if (cortexLearnFlag)
                    {
                        if (Neuron[synConNI].type == 0)
                        {
                            // iSTDP
                            updateStrength = A_iSTDP*(Neuron[synConNI].iSTDP_trace - 1.5*targetRate*(tau_inhib/1000.0)); //convert tau from ms to sec
                            Neuron[synConNI].corticalStrength[l] = Neuron[synConNI].corticalStrength[l] + updateStrength*gmax_cortical;
                            if (Neuron[synConNI].corticalStrength[l] > 2.0*gmax_cortical)
                            {Neuron[synConNI].corticalStrength[l] = 2.0*gmax_cortical;}
                            if (Neuron[synConNI].corticalStrength[l]<0.0)
                            {Neuron[synConNI].corticalStrength[l] = 0.0;}
                        }
                    }
                }
            }
            
            // GJ-coupled cell updates, if spiked:
            numGJ = Neuron[l].GJConn.size();
            if (numGJ > 0)
            {
                GJcell = Neuron[l].GJConn[0];
                if (Neuron[GJcell].v > vReset)
                    {Neuron[GJcell].v = Neuron[GJcell].v + spikeletSize;}
                    if (Neuron[GJcell].v > vT)
                    {
                        Neuron[GJcell].o1 = Neuron[GJcell].o1 + 1.0;
                        Neuron[GJcell].o2 = Neuron[GJcell].o2 + 1.0;
                        Neuron[GJcell].r1 =  Neuron[GJcell].r1 + 1.0;
                        Neuron[GJcell].iSTDP_trace =  Neuron[GJcell].iSTDP_trace + 1.0;

                        Neuron[GJcell].v = vReset;
                        Neuron[GJcell].tspN.push_back(t+dt);
                        numLGNconnect = Neuron[GJcell].synapticConn.size();
                        for (int L = 0; L<numLGNconnect; L++)
                        {
                            indexL = Neuron[GJcell].synapticConn[L];
                            lengthTsp = Neuron[GJcell].tsp_input[indexL].size();
                            if (lengthTsp > 0)
                            {       // LTP from LGN to GJ cell that was pushed over threshold
                                    updateStrength = r1[indexL]*A_trip_LTP*Neuron[GJcell].o2_past;
                                    Neuron[GJcell].synapStrength[indexL] = Neuron[GJcell].synapStrength[indexL] + updateStrength*gmax;
                            }
                            if (Neuron[GJcell].synapStrength[indexL] > gmax)
                            { Neuron[GJcell].synapStrength[indexL] = gmax;}
                        }
                        numSynConn_elect = Neuron[GJcell].corticalConn.size();
                        for (int ii=0; ii<numSynConn_elect; ii++)
                        {
                            foo_synConN = Neuron[GJcell].corticalConn[ii];
                            Neuron[foo_synConN].g_excite1 = Neuron[foo_synConN].g_excite1 + Neuron[foo_synConN].corticalStrength[GJcell]*exp(-(t - Tsp)*oversigmaE);
                            if (cortexLearnFlag)
                            {
                                if (Neuron[foo_synConN].type == 0)
                                { // LTD
                                    updateStrength = -Neuron[foo_synConN].o1*Neuron[foo_synConN].LTD_cort;
                                    Neuron[foo_synConN].corticalStrength[GJcell] = Neuron[foo_synConN].corticalStrength[GJcell] + updateStrength*gmax_cortical;
                                  // LTP
                                    updateStrength2 = Neuron[GJcell].o2_past*Neuron[foo_synConN].r1*A_trip_LTP_cort;
                                    Neuron[GJcell].corticalStrength[foo_synConN] = Neuron[GJcell].corticalStrength[foo_synConN] + updateStrength2*gmax_cortical;
                                   
                                    if (Neuron[foo_synConN].corticalStrength[GJcell] < 0.0)
                                    {Neuron[foo_synConN].corticalStrength[GJcell] = 0.0;}
                                    if (Neuron[GJcell].corticalStrength[foo_synConN] > gmax_cortical)
                                    {Neuron[GJcell].corticalStrength[foo_synConN] = gmax_cortical;}
                                }
                                else
                                {  // iSTDP
                                    updateStrength_iSTDP = A_iSTDP*Neuron[foo_synConN].iSTDP_trace;
                                    Neuron[GJcell].corticalStrength[foo_synConN] = Neuron[GJcell].corticalStrength[foo_synConN] + updateStrength_iSTDP*gmax_cortical;
                                    if (Neuron[GJcell].corticalStrength[foo_synConN] > 2.0*gmax_cortical)
                                    {  Neuron[GJcell].corticalStrength[foo_synConN] = 2.0*gmax_cortical;}
                                }
                            }
                        }
                    }
            }
            Neuron[l].Espike = false;
        }
    }
}

void poissonBackground(Network* Neuron, int N, double t, double dt)
{
    double t_SP = t;
    double r, u;
    int n;
    while (t_SP <= t+dt)
    {
        r = rand();
        u = r/RAND_MAX;
        while (u == 0)
        {
            r = rand();
            u = r/RAND_MAX;
        }
        t_SP = -log(u)/(N*backgroundNu) + t_SP;
        if (t_SP <= t+dt)
        {
            n  = rand()%N;
            Neuron[n].g_excite1 = Neuron[n].g_excite1 + Neuron[n].backgroundF*exp(-(t - t_SP)*oversigmaE);
            Neuron[n].countExternalSpikes = Neuron[n].countExternalSpikes + 1;
        }
    }
}

#endif
