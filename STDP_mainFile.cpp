//
// STDP main file
//
// triplet rule with homeostatic drive
// includes inhibitory neurons and iSTDP from I --> E
//
//
//  Created by Jen Crodelle on 6/10/2019
//
//
#include "STDP_functions.hpp"
#include <iostream>

int main(int argc, char * argv[]){
    // input params
    int N = atof(argv[1]); // number of neurons
    double Tfin = atof(argv[2]);  // final time (ms)
    double timeForSynapses=atof(argv[3]); //time to start recurrent synapses
    double seed = atof(argv[4]);    //random seed
    double gmax_cortical = atof(argv[5]); // max cortical weight
    double probGJs = atof(argv[6]); // prob of connecting with GJ
    int flagForRadius = atof(argv[7]); // flag to restrict to radius (1 --> rad, 0 --> all-to-all)
   
    if (argc != 8)
    { cout << " there are " << argc << " input arguments, but there should be 6!" << endl;
        return 0; }
// setting up simulation params
    srand(seed);
    int T = (int)20*Tfin;
    double dt = Tfin/T; //step size
    cout << "final time  = " << Tfin/1000.0 << " seconds" << endl;
    
    int radiusConn[2][2] = {{4,10},{10,10}};
    double probCortConnect = 1.0; //prob of cortical cells connecting
    double strengthConnect[2][2] = {{0.0*gmax_cortical,0.1*gmax_cortical},{0.3*gmax_cortical,0.1*gmax_cortical}};
    Network* Neuron = new Network[N];

    // create the struct neurons
    create_neurons(Neuron, N, T);
    cout << "There are " << N << " neurons with " << N_input << " synapses" << endl;
   
    if (flagForRadius == 1)
    {
    // set up synaptic structure with radius
    setUpZeroSynapses_radius(Neuron,N,radiusConn);
    cout << "Cortical synapses in radius" << endl;
    parameterFile << "radius synapses" << endl;
    }
    else
    {
       // set up the syapses as all-to-all
      setUpZeroSynapses_allToAll(Neuron, N,  probCortConnect);
      cout << "Cortical synapses all-to-all" << endl;
      parameterFile << "all-to-all synapses" << endl;
    }
    
   // connect sister cells with GJs
    connectGJs_sisters(Neuron, N, probGJs);
    parameterFile << "Connect sisters with GJs" << endl;
    
    // set initial interval length for stimulating LGN:
    double  r = rand();
    double  u = r/RAND_MAX;
    bool errorflag = false;
    ::t_interval = -log(u)*(20.0); //initial interval..
    
    bool cortexLearnFlag;  // flag for when to turn cortical learning on
    double updateStrength, numGJconn, sumGJ, gc, spikeletSize, A_trip_LTP;
    int indexL, numLGNconnect, lengthTsp, foo;
    
    int onceFlag = 1;
    // time loop:
    for (int i = 0; i < T; i++)
    {
        if (i%200000 == 0) //only sometimes record the weights
        {
            cout << "time: " << t/1000.0 << "sec"<< endl;
            for (int KK = 0; KK < N; KK++)
            {
                for (int L = 0; L < N_input; L++) //LGN weights
                {
                    if (L==N_input-1)
                    {W << Neuron[KK].synapStrength[L] << endl;}
                    else
                    {W << Neuron[KK].synapStrength[L] << ",";}
                }
                for (int J = 0; J < N; J++) //recurrent cortical weights
                {  if (J==N-1)
                    { W_cortical << Neuron[KK].corticalStrength[J]  << endl;}
                    else { W_cortical << Neuron[KK].corticalStrength[J] << ",";}
                }
            }
        }
        // only allow neurons to learn after 3 seconds of background:
        if(t < 3000.0)
        {
            A_trip_LTP = 0.0;
            cortexLearnFlag = false;
            gc = 0.06;
            spikeletSize = 1.0;
        }
        // feedforward synapses and GJs
        else if (t < timeForSynapses)
        {
            A_trip_LTP = input_A_tripLTP;
            cortexLearnFlag = false;
            gc = 0.06;
            spikeletSize = 1.0;
        }
        // now feedforward AND recurrent (turn off GJs)
        else
        {
            A_trip_LTP = input_A_tripLTP;
            cortexLearnFlag = true;
            gc = 0.0;
            spikeletSize = 0.0;
            // once we allow recurrent connections set synaptic connections
            if (onceFlag == 1)
            {
                connectSynapses(Neuron, N, strengthConnect,gmax_cortical,flagForRadius);
                cout << "Connected synapses at time t = " << t/1000.0 << " sec" << endl;
                onceFlag = 0;
            }
        }
        // calculate tracer var for LGN
        for(int L = 0; L<N_input; L++)
        {
            r1[L] = r1[L]*exp(-dt/tau_LTP);
        }
        // update conductance, voltage, etc for V1 neurons
        for (int KK = 0; KK<N; KK++)
        {
            // update conductance:
            Neuron[KK].g_excite1 = Neuron[KK].g_excite0*exp(-dt*oversigmaE);
            Neuron[KK].g_inhib1 = Neuron[KK].g_inhib0*exp(-dt*oversigmaI);
            Neuron[KK].avgVolt = Neuron[KK].avgVolt*exp(-dt*overTau_v);
            // update tracers:
            Neuron[KK].o1 = Neuron[KK].o1*exp(-dt/tau_LTD);
            Neuron[KK].o2 = Neuron[KK].o2*exp(-dt/trip_tau_LTD);
            Neuron[KK].r1 = Neuron[KK].r1*exp(-dt/tau_LTP);
            Neuron[KK].iSTDP_trace = Neuron[KK].iSTDP_trace*exp(-dt/tau_inhib);
            
            // calculate GJ stuff:
            numGJconn = Neuron[KK].GJConn.size();
            if (numGJconn > 0)
            { // loop over GJ connections and calculate sum
          //      vT_local = -45.0;
                sumGJ = 0.0;
                for (int J = 0; J<numGJconn; J++)
                {
                    foo = Neuron[KK].GJConn[J];
                    sumGJ = sumGJ + Neuron[foo].v;
                }
            }
            else
            {
                gc = 0.0;
              //  vT_local = -45.0;
                sumGJ = 0.0;
                numGJconn = 0.0;
            }
            //calculate variables: a0,a1,b0,b1,k1,k2 for RK2 Method
            oldv = Neuron[KK].v;
            gE0 = Neuron[KK].g_excite0;
            gI0 = Neuron[KK].g_inhib0;
            
            a0 = (gL+gE0+gI0+gc*numGJconn)*overC;
            b0 = (gL*vR+gE0*vE+gI0*vI+gc*sumGJ)*overC;
            
            gE1 = Neuron[KK].g_excite1;
            gI1 = Neuron[KK].g_inhib1;
            
            a1 = (gL+gE1+gI1+gc*numGJconn)*overC;
            b1 = (gL*vR+gE1*vE+gI1*vI+gc*sumGJ)*overC;
            
            k1 = -a0*oldv+b0;
            k2 = -a1*(oldv+dt*k1)+b1;
            
            Neuron[KK].v = oldv+(dt/2)*(k1+k2); //RK2 step
            
            // check if voltage went below rest...
            if (Neuron[KK].v < vR)
            {   cout << " We have a problem, v = " << Neuron[KK].v << endl;
                errorflag = true;
                break;
            }
            // if reached threshold, neuron fired, update conductance and synaptic weights:
            if (Neuron[KK].v > vT)
            {
                // calculate spike time
                tspike = t+dt*((vT-oldv)/(Neuron[KK].v-oldv));
                vtilda = (vReset - 0.5*(tspike-t)*(b0+b1-dt*a1*b0))/(1.0+0.5*(tspike-t)*(-a0-a1+dt*a0*a1));
                k1tilda = -a0*vtilda+b0;
                k2tilda = -a1*(vtilda+k1tilda*dt)+b1;
                // new voltage at end of time step
                Neuron[KK].v = vtilda+(dt/2)*(k1tilda+k2tilda);
                
                // if this voltage is still over threshold, error
                if (Neuron[KK].v > vT)
                {   cout << "Two spikes occured in one time step. Vtilda was above threshold." << endl;
                    errorflag = true;
                    break;
                }
                
                //collect spike time for this neuron
                Neuron[KK].tspN.push_back(tspike);
                // calc average voltage (note overTau_v is in ms)
                Neuron[KK].avgVolt = Neuron[KK].avgVolt + 1000.*overTau_v;
                Neuron[KK].Espike = true;
                
                numLGNconnect = Neuron[KK].synapticConn.size();
                // update weight from LGN
                for (int L = 0; L<numLGNconnect; L++)
                {
                    indexL = Neuron[KK].synapticConn[L];
                    lengthTsp = Neuron[KK].tsp_input[indexL].size();
                    // if LGN spiked previously
                    if (lengthTsp > 0)
                    {
                        updateStrength = r1[indexL]*A_trip_LTP*Neuron[KK].o2_past;
                        Neuron[KK].synapStrength[indexL] = Neuron[KK].synapStrength[indexL] + updateStrength*gmax;
                    }
                    // if weight goes above max, set to max:
                    if (Neuron[KK].synapStrength[indexL] > gmax)
                    { Neuron[KK].synapStrength[indexL] = gmax;}
                }
            }
        }
        
        // find incoming spikes from LGN and update synapses onto cortical cells
        if(t < ::t_interval)
        { inputFromLGN(Neuron, N, t, dt, 0);}
        else
        { inputFromLGN(Neuron,N, t, dt, 1);}
        // background drive to V1 neurons
        poissonBackground(Neuron, N, t, dt);
        // update conductances and weights for V1 neurons
        spikeUpdate(Neuron, N, t, dt, A_trip_LTP, cortexLearnFlag, gmax_cortical, spikeletSize);

        //shift all current conductance values to be old g calculate new g at the next timestep
        for (int c = 0; c < N; c++)
        {
            Neuron[c].g_excite0 = Neuron[c].g_excite1;
            Neuron[c].g_inhib0 = Neuron[c].g_inhib1;
            Neuron[c].o2_past = Neuron[c].o2;
            // update the LTD for each postsynaptic neuron (note targetRate is in 1/s)
            Neuron[c].LTD = Neuron[c].avgVolt*Neuron[c].avgVolt*(A_trip_LTP*tau_LTP*trip_tau_LTD)/(targetRate*1000.0*tau_LTD);
            Neuron[c].LTD_cort =  Neuron[c].avgVolt*Neuron[c].avgVolt*(A_trip_LTP_cort*tau_LTP*trip_tau_LTD)/(targetRate*1000.0*tau_LTD);
        }
        if (errorflag == true)
        {break;}
        t = t+dt;
    }
    
    // read in data and calculate some information
    double avgInputSpikes = 0.0;
    double avgExternalSpikes = 0.0;
    double avgExcCorticalSpikes = 0.0;
    double avgInhibCorticalSpikes = 0.0;
    for (int k = 0; k<N; k++)
    {
        avgInputSpikes = avgInputSpikes + Neuron[k].countInputSpikes;
        avgExternalSpikes = avgExternalSpikes + Neuron[k].countExternalSpikes;
        double lengthSp = Neuron[k].tspN.size();
        if (Neuron[k].type == 0)
        {avgExcCorticalSpikes = avgExcCorticalSpikes + lengthSp;}
        else{avgInhibCorticalSpikes = avgInhibCorticalSpikes + lengthSp; }
        spTimes << Neuron[k].type << "," <<  Neuron[k].sisterID << "," << lengthSp << ",";
        if (lengthSp == 0)
        { spTimes << endl;
        }
        for (int j = 0; j<lengthSp; j++)
        {
            if (j == lengthSp - 1)
            { spTimes << Neuron[k].tspN[j] << endl;
            }
            else
            { spTimes << Neuron[k].tspN[j] << ","; }
        }
        int numElect = Neuron[k].GJConn.size();
        electConn_file << numElect << ",";
        if (numElect > 0)
        {
        for (int c=0; c<numElect; c++)
        {
            if (c == numElect-1)
            { electConn_file << Neuron[k].GJConn[c] << endl;}
            else
            { electConn_file << Neuron[k].GJConn[c] << ",";}
        }
        }
        else
        {electConn_file << endl;}
        for (int L = 0; L<N_input; L++)
        {
            if(L == N_input-1)
            {   W << Neuron[k].synapStrength[L] << endl;
                finalWeights << Neuron[k].synapStrength[L] << endl;}
            else
            {    W << Neuron[k].synapStrength[L] << ",";
                finalWeights << Neuron[k].synapStrength[L] << ",";}
        }
    }
    
    avgInputSpikes = (1000.0*avgInputSpikes)/((double)N_input*probSynapseConn*N*Tfin);
    avgExcCorticalSpikes =(1000.0*avgExcCorticalSpikes)/((double)(N-numInhib)*Tfin);
    avgInhibCorticalSpikes =(1000.0*avgInhibCorticalSpikes)/((double)numInhib*Tfin);
    cout << " Avg incoming FR from LGN cells = " << avgInputSpikes << endl;
    cout << " Avg exc cortical FR  = " << avgExcCorticalSpikes << endl;
    cout << " Avg inhib cortical FR  = " << avgInhibCorticalSpikes << endl;

    // read parameters into file
    parameterFile << "N = " << N << endl;
    parameterFile<< "Final time = " << Tfin << endl;
    parameterFile<< "dt = " << dt << endl;
    parameterFile << " rand seed = " << seed << endl;
    
   // neuron parameters
    parameterFile << "NEURON PARAMETERS: " << endl;
    parameterFile << "gL = " << gL << endl;
    parameterFile << "Vrest = "  << vR << endl;
    parameterFile << "gc = 0.06" << endl;
    parameterFile<< "Vthreshold = " << vT << endl;
    parameterFile << "Sigma_E =  " << sigma_E << endl;
    parameterFile << "Sigma_E=I =  " << sigma_I << endl;
    parameterFile << "background F = " << backgroundF << endl;
    parameterFile << "background Nu = " << backgroundNu << endl;
    parameterFile << "external rate R0 = " << R0 << endl;
    // connectivity parameters
     parameterFile << "CONNECTIVITY PARAMETERS: " << endl;
    parameterFile << "radius E-E" <<radiusConn[0][0] << endl;
    parameterFile << "radius I-E" <<radiusConn[0][1] << endl;
    parameterFile << "radius E-I" <<radiusConn[1][0] << endl;
    parameterFile << "radius I-I" <<radiusConn[1][1] << endl;
    parameterFile << "number of LGN synapses " << N_input << endl;
    parameterFile << "prob of GJ = " << probGJs << endl;

    // learning params
    parameterFile << "LEARNING PARAMETERS: " << endl;
    parameterFile<<" A3_LTP = " << A_trip_LTP << endl;
    parameterFile << "A3 LTP cort = " << A_trip_LTP_cort << endl;
    parameterFile<< "tau_LTD = " << tau_LTD << endl;
    parameterFile << "tau_LTP = " << tau_LTP << endl;
    parameterFile << "trip_tau_LTD " << trip_tau_LTD << endl;
    parameterFile<< "gmax = " << gmax << endl;
    parameterFile<< "gmax_cortical = " << gmax_cortical << endl;
    parameterFile << "target rate = " << targetRate << endl;
    // output information
    parameterFile << "SIMULATION OUTPUT INFORMATION: " << endl;
    parameterFile<< "avg input from LGN neuron = " << avgInputSpikes << endl;
    parameterFile << " Avg exc cortical FR  = " << avgExcCorticalSpikes << endl;
    parameterFile << " Avg inhib cortical FR  = " << avgInhibCorticalSpikes << endl;
    parameterFile << "real run time for simulation = " << clock()/CLOCKS_PER_SEC << endl;
    // close csv files
    W.close();
    finalWeights.close();
    spTimes.close();
    parameterFile.close();
    electConn_file.close();
    delete[] Neuron;
    
    cout << clock()/CLOCKS_PER_SEC << " seconds " << endl;
}

