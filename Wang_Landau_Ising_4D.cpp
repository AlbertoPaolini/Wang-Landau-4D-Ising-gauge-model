#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>     
#include <time.h>
#include<vector>
#include <fstream>
#include <filesystem>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <thread>

#define Accuracy 0.9
#define Ns 8
#define Nz 8
#define Nt 8

using namespace std;

double Random_01() {
    thread_local static mt19937_64 rng(random_device{}());  
    thread_local static uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

//*****************************************************************************************************

int Random_int(int a, int b) {
    thread_local static mt19937_64 rng(random_device{}());
    uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}

//*****************************************************************************************************

int Random_int2(int a, int b) {

    thread_local static mt19937_64 rng(random_device{}());
    thread_local static uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}

//*****************************************************************************************************

int Indexposition(int x, int y, int z, int t){
  return ( (x+Ns) %Ns ) * Ns * Nz * Nt + ( (y+Ns) %Ns ) * Nz * Nt + ( (z+Nz) %Nz ) * Nt + ( (t+Nt) %Nt );
}

//*****************************************************************************************************

int sum_plaq(vector< int > &Plaquette){
    int sum = 0;
    int plaquette_size = Plaquette.size();

    for(int i=0; i<plaquette_size; i++){
        sum += Plaquette[i];
    }

    return sum;
}


//*****************************************************************************************************

bool Flat(vector<int> H){
    double Average =0;
    int min = H[0];
    int bin = 0, H_size = H.size();
    int num_zero = 0;
    
    Average += H[0] + H[12] + H[H_size-1] + H[H_size-1-12];
    bin = 4;

    for(int i=20; i<H.size()-1-20; i = i+4){
        Average += H[i];
        bin++;

        if(H[i]< min) min = H[i];
        if(H[i] == 0) num_zero++;
    }

    Average = Average / (double)bin;
    
    //cout<< (double)min/Average <<endl;
    //cout<< num_zero <<endl;


    return ( min > Accuracy * Average);
}

//*****************************************************************************************************

void lattice_creation_minus_frozen(vector< vector<int> > &lattice, vector< vector<int> > &map_nn){ 
  
  int nn;

  for(int i=0; i<Ns; i++){
    for(int j=0; j<Ns; j++){
        for(int k=0; k<Nz; k++){
            for(int l=0; l<Nt; l++){
              vector<int> element;
              element.push_back(-1); // This is the link along the direction 0 from the site of coordinate (i,j,k,l)
              element.push_back(-1); // This is the link along the direction 1 from the site of coordinate (i,j,k,l)
              element.push_back(-1); // This is the link along the direction 2 from the site of coordinate (i,j,k,l)
              element.push_back(-1); // This is the link along the direction 3 from the site of coordinate (i,j,k,l)
              
            //   element.push_back(1); // This is the link along the direction 0 from the site of coordinate (i,j,k,l)
            //   element.push_back(1); // This is the link along the direction 1 from the site of coordinate (i,j,k,l)
            //   element.push_back(1); // This is the link along the direction 2 from the site of coordinate (i,j,k,l)
            //   element.push_back(1); // This is the link along the direction 3 from the site of coordinate (i,j,k,l)

              lattice.push_back(element);

              // Here we initialize a matrix of vector that contain the position of the near-neighbours of each site
              vector<int> nn_position;
              nn = Indexposition( i +1, j, k, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j +1, k, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j, k +1, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j, k, l +1 ); 
              nn_position.push_back(nn);

              nn = Indexposition( i -1, j, k, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j -1, k, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j, k -1, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j, k, l -1 );
              nn_position.push_back(nn);

              map_nn.push_back(nn_position);

            }
        }
    }
  }


}

//*****************************************************************************************************

void lattice_creation_frozen(vector< vector<int> > &lattice, vector< vector<int> > &map_nn){ 
  
  int nn;

  for(int i=0; i<Ns; i++){
    for(int j=0; j<Ns; j++){
        for(int k=0; k<Nz; k++){
            for(int l=0; l<Nt; l++){
              vector<int> element;
              element.push_back(1); // This is the link along the direction 0 from the site of coordinate (i,j,k,l)
              element.push_back(1); // This is the link along the direction 1 from the site of coordinate (i,j,k,l)
              element.push_back(1); // This is the link along the direction 2 from the site of coordinate (i,j,k,l)
              element.push_back(1); // This is the link along the direction 3 from the site of coordinate (i,j,k,l)
              
            //   element.push_back(1); // This is the link along the direction 0 from the site of coordinate (i,j,k,l)
            //   element.push_back(1); // This is the link along the direction 1 from the site of coordinate (i,j,k,l)
            //   element.push_back(1); // This is the link along the direction 2 from the site of coordinate (i,j,k,l)
            //   element.push_back(1); // This is the link along the direction 3 from the site of coordinate (i,j,k,l)

              lattice.push_back(element);

              // Here we initialize a matrix of vector that contain the position of the near-neighbours of each site
              vector<int> nn_position;
              nn = Indexposition( i +1, j, k, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j +1, k, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j, k +1, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j, k, l +1 ); 
              nn_position.push_back(nn);

              nn = Indexposition( i -1, j, k, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j -1, k, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j, k -1, l );
              nn_position.push_back(nn);
              nn = Indexposition( i, j, k, l -1 );
              nn_position.push_back(nn);

              map_nn.push_back(nn_position);

            }
        }
    }
  }


}

//*****************************************************************************************************

void latticeRandom_creation(vector< vector<int> > &lattice, vector< vector<int> > &map_nn){
    
  int nn;
  int spin[2]={1,-1}, random;

  for(int i=0; i<Ns; i++){
    for(int j=0; j<Ns; j++){
        for(int k=0; k<Nz; k++){
            for(int l=0; l<Nt; l++){
              vector<int> element; 

                random=Random_int(0,1);
                element.push_back(spin[random]); // This is the link along the direction 0 from the site of coordinate (i,j,k,l)
                random=Random_int(0,1);
                element.push_back(spin[random]); // This is the link along the direction 1 from the site of coordinate (i,j,k,l)
                random=Random_int(0,1);
                element.push_back(spin[random]); // This is the link along the direction 2 from the site of coordinate (i,j,k,l)
                random=Random_int(0,1);
                element.push_back(spin[random]); // This is the link along the direction 3 from the site of coordinate (i,j,k,l)

                element.push_back(1);
                element.push_back(1);
                element.push_back(1);
                element.push_back(1);

                lattice.push_back(element);

                // Here we initialize a matrix of vector that contain the position of the near-neighbours of each site
                vector<int> nn_position;
                nn = Indexposition( i +1, j, k, l );
                nn_position.push_back(nn);
                nn = Indexposition( i, j +1, k, l );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k +1, l );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k, l +1 );
                nn_position.push_back(nn);

                nn = Indexposition( i -1, j, k, l );
                nn_position.push_back(nn);
                nn = Indexposition( i, j -1, k, l );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k -1, l );
                nn_position.push_back(nn);
                nn = Indexposition( i, j, k, l -1 );
                nn_position.push_back(nn);

                map_nn.push_back(nn_position);

            }
        }
    }
  }
  
//   for(int i=0; i<lattice.size(); i++){
//     for(int direction = 0; direction < 4; direction ++){
//       lattice[map_nn[i][direction]][4 + direction] = lattice[i][direction];
//     }
//   }

}

//*****************************************************************************************************

void plaquette_creation(vector< vector<int> > &lattice, vector< vector<int> > &map_nn, vector< int > &Plaquette, 
    vector< vector< vector<int> > > &link_to_plaquette){

    int lattice_size = lattice.size();
    int position = 0;
    Plaquette.clear();
    Plaquette.resize(6 * lattice_size);
    link_to_plaquette.clear();
    link_to_plaquette.resize(lattice_size);
    
    for (int i = 0; i < lattice_size; ++i) {
      link_to_plaquette[i].resize(4);
    }

    for(int i = 0; i<lattice_size; i++){
        int element, nn1, nn2;
        
        nn1 = map_nn[i][0];
        nn2 = map_nn[i][1];
        element = lattice[i][0]*lattice[i][1] * lattice[nn1][1]*lattice[nn2][0];
        link_to_plaquette[i][0].push_back(position);
        link_to_plaquette[i][1].push_back(position);
        link_to_plaquette[nn2][0].push_back(position);
        link_to_plaquette[nn1][1].push_back(position);
        Plaquette[position] = element;
        position ++;

        nn1 = map_nn[i][0];
        nn2 = map_nn[i][2];
        element = lattice[i][0]*lattice[i][2] * lattice[nn1][2]*lattice[nn2][0];
        link_to_plaquette[i][0].push_back(position);
        link_to_plaquette[i][2].push_back(position);
        link_to_plaquette[nn2][0].push_back(position);
        link_to_plaquette[nn1][2].push_back(position);
        Plaquette[position] = element;
        position ++;

        nn1 = map_nn[i][0];
        nn2 = map_nn[i][3];
        element = lattice[i][0]*lattice[i][3] * lattice[nn1][3]*lattice[nn2][0];
        link_to_plaquette[i][0].push_back(position);
        link_to_plaquette[i][3].push_back(position);
        link_to_plaquette[nn2][0].push_back(position);
        link_to_plaquette[nn1][3].push_back(position);
        Plaquette[position] = element;
        position ++;

        nn1 = map_nn[i][1];
        nn2 = map_nn[i][2];
        element = lattice[i][1]*lattice[i][2] * lattice[nn1][2]*lattice[nn2][1];
        link_to_plaquette[i][1].push_back(position);
        link_to_plaquette[i][2].push_back(position);
        link_to_plaquette[nn2][1].push_back(position);
        link_to_plaquette[nn1][2].push_back(position);
        Plaquette[position] = element;
        position ++;

        nn1 = map_nn[i][1];
        nn2 = map_nn[i][3];
        element = lattice[i][1]*lattice[i][3] * lattice[nn1][3]*lattice[nn2][1];
        link_to_plaquette[i][1].push_back(position);
        link_to_plaquette[i][3].push_back(position);
        link_to_plaquette[nn2][1].push_back(position);
        link_to_plaquette[nn1][3].push_back(position);
        Plaquette[position] = element;
        position ++;

        nn1 = map_nn[i][3];
        nn2 = map_nn[i][2];
        element = lattice[i][2]*lattice[i][3] * lattice[nn1][2]*lattice[nn2][3];
        link_to_plaquette[i][2].push_back(position);
        link_to_plaquette[i][3].push_back(position);
        link_to_plaquette[nn2][3].push_back(position);
        link_to_plaquette[nn1][2].push_back(position);
        Plaquette[position] = element;
        position ++;

    }

}

//*****************************************************************************************************

void Wang_Landau_step(vector<double> &log_g, vector< vector<int> > &Lattice, vector< int > &Plaquette, int Num_Energy_levels, 
    vector< vector< vector< int > > > &link_to_plaquette, vector< int > &H, long double f){
    
    int current_energy = sum_plaq(Plaquette);
    int final_index;
    int lattice_size = Lattice.size();
    
    for(int iteration = 0; iteration < 200000000; iteration++){
                
        int i = Random_int(0, lattice_size-1);
                
        for(int mu = 0; mu < 4; mu++){
                    
            double dE = 0.0;
            for(int j = 0; j< link_to_plaquette[i][mu].size(); j++){
                dE += (double)Plaquette[link_to_plaquette[i][mu][j]];
            }
            dE *= -2.0;
                    
            int current_index = current_energy + 6 * lattice_size;
            int new_energy = current_energy + dE;
            int new_index = new_energy + 6 * lattice_size;
            double dlog = log_g[current_index] - log_g[new_index];

            double prob;
            if( exp(dlog) < 1) prob = exp(dlog);
            else prob = 1.0;

            if(Random_01() < prob){
                current_energy = new_energy;
                final_index = new_index;
                for(int j = 0; j<link_to_plaquette[i][mu].size(); j++){
                    Plaquette[link_to_plaquette[i][mu][j]] *= -1;
                }
            }
            else{
                final_index = current_index;
            }

            log_g[final_index] += log(f);
            H[final_index] += 1;
        }

    }
}
//*****************************************************************************************************

void Wang_Landau(vector<double> &log_g, vector< vector<int> > &Lattice, vector< vector<int> > &Lattice2, vector< vector<int> > &Lattice3, vector< vector<int> > &Lattice4, 
    vector< int > &Plaquette, vector< int > &Plaquette2, vector< int > &Plaquette3, vector< int > &Plaquette4, int Num_Energy_levels, 
    vector< vector< vector< int > > > &link_to_plaquette){

    long double Threshold = 1.000000001;
    double f = exp(1.0);
    bool H_Flat = false;
    int lattice_size = Lattice.size();
    int current_energy = sum_plaq(Plaquette);
    int final_index;
    vector<double> log_g1(Num_Energy_levels, 0.0), log_g2(Num_Energy_levels, 0.0), log_g3(Num_Energy_levels, 0.0), log_g4(Num_Energy_levels, 0.0);

    while(f > Threshold){
        cout<<f<<endl;
        vector<int> H(Num_Energy_levels, 0);
        vector<int> H1(Num_Energy_levels, 0), H2(Num_Energy_levels, 0), H3(Num_Energy_levels, 0), H4(Num_Energy_levels, 0);
        H_Flat = false;
        
        while(H_Flat == false){
            thread t1(Wang_Landau_step, ref(log_g1), ref(Lattice), ref(Plaquette), Num_Energy_levels, ref(link_to_plaquette), ref(H1), f);
            thread t2(Wang_Landau_step, ref(log_g2), ref(Lattice2), ref(Plaquette2), Num_Energy_levels, ref(link_to_plaquette), ref(H2), f);
            thread t3(Wang_Landau_step, ref(log_g3), ref(Lattice3), ref(Plaquette3), Num_Energy_levels, ref(link_to_plaquette), ref(H3), f);
            thread t4(Wang_Landau_step, ref(log_g4), ref(Lattice4), ref(Plaquette4), Num_Energy_levels, ref(link_to_plaquette), ref(H4), f);
            t1.join();
            t2.join();
            t3.join();
            t4.join();

            for(int h = 0; h< H1.size(); h++){
                H[h] = H1[h] + H2[h] + H3[h] + H[4];
            }
            H_Flat = Flat(H);
        }

        f = sqrt(f);
            
    }

    for(int i = 0; i<log_g1.size(); i++){
        log_g[i] = (log_g1[i] + log_g2[i] + log_g3[i] + log_g[4])/4.0; 
    }

}

//*****************************************************************************************************

long double Energy(const vector<double> &log_g, double beta, int shift) {
    int N = (int)log_g.size();
    vector<long double> exps(N);
    long double max_val = -INFINITY;

    for (int i = 0; i < N; i++) {
        int energy = i - shift;
        long double val = (long double)log_g[i] + (long double)beta * (long double)energy;
        exps[i] = val;
        if (val > max_val) max_val = val;
    }

    long double Z = 0.0;
    long double Eavg = 0.0;

    for (int i = 0; i < N; i++) {
        int energy = i - shift;
        long double w = expl(exps[i] - max_val); // usa expl per long double
        Z += w;
        Eavg += (long double)energy * w;
    }

    return Eavg / Z;
}

//*****************************************************************************************************

long double EnergyUp2(const vector<double> &log_g, double beta, int shift) {
    int N = (int)log_g.size();
    vector<long double> exps(N);
    long double max_val = -INFINITY;

    // costruisco i valori dell'esponente (log_g + beta*energy) e prendo il massimo
    for (int i = 0; i < N; i++) {
        int energy = i - shift;
        long double val = (long double)log_g[i] + (long double)beta * (long double)energy;
        exps[i] = val;
        if (val > max_val) max_val = val;
    }

    long double Z = 0.0;
    long double Eavg = 0.0;

    // sommo pesi stabili: w = exp( exps[i] - max_val )
    for (int i = 0; i < N; i++) {
        int energy = i - shift;
        long double w = expl(exps[i] - max_val); // usa expl per long double
        Z += w;
        Eavg += (long double)energy * energy * w;
    }

    return Eavg / Z;
}

//*****************************************************************************************************

int main(){

    ofstream energy_file;
    energy_file.open("energy.txt");

    vector< vector<int> > Lattice, Lattice2, Lattice3, Lattice4;
    vector< vector<int> > NN_Map;

    lattice_creation_frozen(Lattice, NN_Map);
    lattice_creation_minus_frozen(Lattice2, NN_Map);
    latticeRandom_creation(Lattice3, NN_Map);
    latticeRandom_creation(Lattice4, NN_Map);
    vector< int > Plaquette, Plaquette2, Plaquette3, Plaquette4;
    vector< vector< vector< int > > > link_to_plaquette;
    plaquette_creation(Lattice, NN_Map, Plaquette, link_to_plaquette);
    plaquette_creation(Lattice2, NN_Map, Plaquette2, link_to_plaquette);
    plaquette_creation(Lattice3, NN_Map, Plaquette3, link_to_plaquette);
    plaquette_creation(Lattice4, NN_Map, Plaquette4, link_to_plaquette);

    int Num_Energy_levels = 2 * 6 * Lattice.size() + 1;
    vector<double> log_g(Num_Energy_levels, 0.0); 
    Wang_Landau(log_g, Lattice, Lattice2, Lattice3, Lattice4, Plaquette, Plaquette2, Plaquette3, Plaquette4, Num_Energy_levels, link_to_plaquette);
    
    double beta = 0.05;
    for( int measurements = 0; measurements < 100 ; measurements++){
        long double E_avg = Energy(log_g, beta, 6 * Lattice.size());
        long double E2_avg = EnergyUp2(log_g, beta, 6 * Lattice.size());
        cout<<beta<<" Average plaquette value: "<<E_avg/(6.0 * Lattice.size())<<" Susceptibility: "<<(E2_avg - E_avg * E_avg)<<endl;
        beta += 0.01;
        
    }

}