#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_specfun.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

void simuleMB(PnlVect* vect_sim, double T, int N, PnlRng *rng);
double estimeEsp(double T, int N, double a, PnlRng *rng, int M);
double getValExact(double a, double T);
void logBiais(double T, double a, PnlRng *rng, int M, const std::string& chemin_fichier);




double getValExact(double a, double T)
{

    double gamma_plus = pnl_sf_gamma_inc(1/2, a*a/(2*T));
    double gamma_minus = pnl_sf_gamma_inc(-1/2, a*a/(2*T));
    double inv_racine_pi = M_2_SQRTPI /2;

    return T * (1- gamma_plus*inv_racine_pi ) + gamma_minus*a*a*inv_racine_pi/2 ;
            
}


void simuleMB(PnlVect* vect_sim, double T, int N, PnlRng *rng)
{
    double step = T/N ;
    double X = 0.0;
    pnl_vect_set(vect_sim, 0,X);

    double Gi; 
    for(int i = 1; i<N+1; i++)
    {
        Gi = pnl_rng_normal(rng);
        X += sqrt(step) * Gi;
        pnl_vect_set(vect_sim, i,X);

    }

}


double estimeEsp(double T, int N, double a, PnlRng *rng, int M)
{
    double step = T/N ;
    PnlVect * vect_time = pnl_vect_create(M);


    for(int i = 0; i<M; i++)
    {

        PnlVect * vect_sim = pnl_vect_create(N +1); 
        simuleMB(vect_sim, T,N, rng);
        for(int j = 0; j<N+1; j++)
        {
            if(a<= GET(vect_sim,j) || j==N)
            {
                pnl_vect_set(vect_time,i,j*step);
                break;
            }
        }
        pnl_vect_free(&vect_sim);

    }
    return pnl_vect_sum(vect_time)/M ;

    //free
    pnl_vect_free(&vect_time);

}


void logBiais(double T, double a, PnlRng *rng, int M, const std::string& chemin_fichier)
{

    std::ofstream fichier(chemin_fichier);
    if (!fichier.is_open()) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier." << std::endl;
        return;
    }


    // Fermer le fichier
    double esp_estim ;
    double esp_exact ;
    double biais ;
    for(int i = 1; i< 200; i+=20)
    {
        esp_estim = estimeEsp(T,i/10,a,rng,M); 
        // esp_exact = getValExact(a,T);
        esp_exact = 1.915006;
        biais = esp_estim - esp_exact;
        fichier << log(i) << ","<< log(biais) << "\n";
        std::cout << " je suis = " << i << std::endl;

    }

    fichier.close();

}

// void simuleCondMB(double a, double epsilon, PnlRng *rng, double s, double t, double u, PnlVect *res, PnlVect *sim, int M)
// {
//     double Gi = pnl_rng_normal(rng);
//     double zar6a = (t-u)*(a-epsilon)/(t-s)  + (u-s)*a/(t-s)  + sqrt((t-u)*(u-s)/(t-s))*Gi;

//     for(int i=0; i<sim->size; i++)
//     {
//         double wi = pnl_vect_get(sim,i);
//         if(abs(wi - a)<epsilon)
//         {
//             pnl_vect_set(res,0,wi);
//             for(int j = 1; j<M; j++)
//             {
//                 double Bu = pnl_vect_get(j-1) * ()
//                 pnl_vect_set(res,Bu,j);
//             }
//         }
//     }
// }

int main()
{
    PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));

    // PnlVect *vect_sim = pnl_vect_create()
    double a =1;
    double T= 3;


    double N = 20;
    double biais;
    int M = 500000;
    double esp_estim = 0 ;
    // int count =0;
    // for(int N = 5; N<50; N+=5)
    // {
    //     // void simuleMB(PnlVect* vect_sim, double T, int N, PnlRng *rng)
    //     esp_estim += estimeEsp(T,N,a,rng,M); 
    //     count ++;
        
    // }

    // esp_estim = esp_estim/count ;
    // // double esp_exact = getValExact(a,T);
    // double esp_exact = 1.915006;
    // biais = esp_estim - esp_exact;
    // std::cout << " esp_estim = " <<esp_estim << std::endl;

    logBiais(T,a,rng,M,"../../tests/test_mmv_br.txt");

    pnl_rng_free(&rng);
    return 0;
}
// logBiais(double T, double a, PnlRng *rng, const std::string& chemin_fichier)