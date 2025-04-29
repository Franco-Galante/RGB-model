// to compile using gsl library: gcc gptsim.c -lm -lgsl -lgslcblas -o gptsim -Wall -O3

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>


// RANDOM NUMBER GENERATOR
#define SEED 		13L   // default random number generator seed
gsl_rng* rng_;

double uniform01() { return(gsl_rng_uniform_pos(rng_)); }
double uniform(double a, double b) { return(a + gsl_rng_uniform_pos(rng_)*(b-a)); }
unsigned long uniformint(unsigned long a, unsigned long b) { return(a + gsl_rng_uniform_int(rng_,b-a+1)); }
int geom(double p) { return(gsl_ran_geometric(rng_,p)); }	
double exponential(double mu) { return(gsl_ran_exponential(rng_,mu)); } 
double stdgauss(double x) { return(gsl_ran_ugaussian_pdf(x)); } 
double poissonpdf(unsigned int k, double mu) { return(gsl_ran_poisson_pdf (k,mu)); }
unsigned int poisson(double mu) { return gsl_ran_poisson(rng_,mu); }
double gaussian(double m, double sigma) { return(m + gsl_ran_gaussian(rng_,sigma)); }
double pareto(double a, double b) { return gsl_ran_pareto(rng_,a,b); }
long zipf(gsl_ran_discrete_t * zipf) { return((long)gsl_ran_discrete(rng_,zipf)); }
void computeconf(int n, double*b, double *min, double *mean, double *max);

typedef enum {SCEPREGPT,SCEGPT,SCEPOSTGPT,SCETRY} scenario; 
scenario sce = SCEPOSTGPT;
char line[200];

#define FILEOUT1 "gptsim.tra"   // output filename for traces
FILE* fo1;
#define FILEOUT2 "gptsim.met"   // output filename for metrics
FILE* fo2;
#define DEBUG 1	// DEBUG MODE (higher values mean more verbose)

#define MAXT 10000
#define DELTAT 1 // for printing
#define NRUN 400
#define MAXP 10002 // must be at least equal to MAXT/DELTAT + 2

//#define MAXT 200
//#define DELTAT 1 // for printing
//#define NRUN 30000
//#define MAXP 202 // must be at least equal to MAXT/DELTAT + 2

typedef enum {RED,GREEN,BLUE,BLACK} color;
#define anc 4  // number of all colors
#define nc 3   // number of non-black colors
double q_[nc];  // intrinsic (unknown to system) qualities
double p_[nc];  // intrinsic (unknown to system) normalized qualities for relevant answers
char* colornames[anc] = {"RED","GREEN","BLUE","BLACK"};
#define INITBALLS 10	// initial number of non-black balls of each color in each box other than generator
#define INITBLACKS 1 // initial number of black balls in SE box

#define size 5 // number of compartments
typedef enum {GE,TS,WWW,GPT,SE} place;
double mu_[size]; // death rates at places
double Mu_[size]; // aggregate death rates at places
char* placenames[size] = {"GE","TS","WWW","LLM","SE"};

#define nf 8 // flows
typedef enum {H,P,T,I,S,A,U,F} flow;
char flownames[nf] = {'H','P','T','I','S','A','U','F'};
double lam_[nf]; // per-ball rate
double Lam_[nf]; // aggregate rate
double C_[nf];     // biases towards quality
double r_[nc][nf]; // primary replication probs
double alfa;  // lam_[U] = alfa * lam_[A]
int wts; // amplifies lam_[T]
double beta; // frac human feedback: lam_[F] = beta * lam_[A]
double gamma_; // split of Lam_[H] into H' and H''

typedef struct ball {  
    int flow;  // originating flow, -1 if none (ball has been initially added)
    place place;    // current place 
    double gpt_[nc];  // percentage of components originated from GAI system
    double comp_[nc]; // RGB components
} Ball, *Ballptr;
#define MAXBALLSPERPLACE 10000
Ball b_[size][MAXBALLSPERPLACE];
#define NUMBALLSTOMIX 2
#define MIXSE 0  // 0-1: whether to mix balls from SE
#define MIXGAI 0 // 0-1: whether to mix balls from GAI
Ballptr bptrarr_[NUMBALLSTOMIX];
double xibias = 0.5; // xi parameter of model 

// system state: number of balls at places
unsigned long nb_[size];
unsigned long nbl; // number of black balls in SE

// record results
double rnbc_[size][nc][MAXP][NRUN];
unsigned long rnb_[size][MAXP][NRUN];
double rnbl_[MAXP][NRUN];
double acc_[size][MAXP][NRUN];
double div_[size][MAXP][NRUN];


// metrics 
double rfwua_[MAXP][NRUN]; // fraction of used answers in WWW (wrt to used + generated)
double rfqan_[MAXP][NRUN]; // fraction of questions A with no answer
double rfqsn_[MAXP][NRUN]; // fraction of questions S with no answer
double theta_;  // irrelevant threshold
double rfiua_[MAXP][NRUN]; // fraction of irrelevant used answers in WWW
double rairi_[MAXP][NRUN]; // GAI responsability index
double rgpta_[MAXP][NRUN]; // GPT autophagy index
double redrun_[size][MAXP]; // fraction of runs in which majority of answers are irrelevant
double redtot_[size][MAXP]; // denom for redrun


void rplot(void) {
    double C,r;
    FILE* fo3 = fopen("rplot.dat","w");
    fprintf(fo3,"### data for r plot \n");
    fprintf(fo3,"### C\t r_RED \t r_GREEN \t r_BLUE\n");
    for (C=-0.4; C<=10; C+=0.01) {
        fprintf(fo3,"%le",C);
        double denom = 0;
        for(color c=0; c<nc; c++) {
            if(q_[c] + C > 0) denom += (q_[c] + C);
        }
        for(color c=0; c<nc; c++) {
            if(q_[c] + C > 0) r = (q_[c] + C)/denom;
            else r = 0.0;
            fprintf(fo3,"\t%le",r);
        }
        fprintf(fo3,"\n");
    }
    fclose(fo3);
}	


void init(scenario sce) {

    // init random seed
    rng_ = gsl_rng_alloc(gsl_rng_mt19937);	 /* MT19937 Matsumoto, Nishimura */
    gsl_rng_set(rng_,SEED);

    fo1 = fopen(FILEOUT1,"w");
    fo2 = fopen(FILEOUT2,"w");

    // SCENARIO INDEPENDENT PARAMETERS

    /* intrinsic qualities */ 
    q_[RED] = 0.1;
    q_[GREEN] = 0.4;
    q_[BLUE] = 0.5;
    theta_ = 0.2;
    p_[GREEN] = q_[GREEN]/(q_[GREEN]+q_[BLUE]);
    p_[BLUE] = q_[BLUE]/(q_[GREEN]+q_[BLUE]);

    /* amplifies lam_T */
    wts = 3;

    /* biases */
    C_[H] = -0.08;
    C_[P] = 1;
    C_[T] = 0; // not used
    C_[I] = 0;
    C_[S] = 0.1;

    /* death rates */
    mu_[GE] = 0.0;
    mu_[TS] = 0.001;
    mu_[WWW] = 0.01;
    mu_[GPT] = 0.1;
    mu_[SE] = 0.1;

    /* aggregate rates */
    Lam_[P] = 1;
    gamma_ = 0.5;
    
    /* per-ball rates */
    lam_[T] = 0.1;
    lam_[I] = 0.1;

    
    // SCENARIO DEPENDENT PARAMETERS
    switch (sce) {

        case SCEPREGPT:
            alfa = 0;
            beta = 1;
            C_[U] = 0.1;
            C_[F] = -0.08;
            Lam_[H] = 0.1;
            Lam_[S] = 100;
            Lam_[A] = 0;
        break;

        case SCEGPT:
            alfa = 0.4;
            beta = 0.04;
            C_[U] = 0.1;
            C_[F] = 0.1;
            Lam_[H] = 0.1;
            Lam_[S] = 100;
            Lam_[A] = 25;
            break;

        case SCEPOSTGPT:
            alfa = 0.8;
            beta = 0.04;
            C_[U] = 1;
            C_[F] = 0.1;
            Lam_[H] = 1;
            Lam_[S] = 10;
            Lam_[A] = 125;
            break;

        case SCETRY:
            alfa = 0.8;
            beta = 0.04;
            C_[U] = 1;
            C_[F] = 0.1;
            Lam_[H] = 1;
            Lam_[S] = 10;
            Lam_[A] = 125;
            break;    

        default:
            printf("Unrecognized scenario\n");
            exit(1);
    }

    /* derived rates */
    Lam_[U] = alfa * Lam_[A];
    Lam_[F] = beta * Lam_[A];
	
    /* compute r's */
    for(flow f=0; f<nf; f++) {
        double denom = 0;
        for(color c=0; c<nc; c++) {
            if(q_[c] + C_[f] > 0) denom += q_[c] + C_[f];
        }
        for(color c=0; c<nc; c++) {
            if(q_[c] + C_[f] > 0) r_[c][f] = (q_[c] + C_[f])/denom;
            else r_[c][f] = 0.0;
            printf("r_[%d][%d] = %f\n",c,f,r_[c][f]);
        }
    }

    // prints parameters for debug
    if(DEBUG) {
        for(flow f=0; f<nf; f++) {
            printf("Flow %c: lambda: %lf Lambda: %lf C: %lf\n",flownames[f],lam_[f],Lam_[f],C_[f]);
        }
        for(place p=0; p<size; p++) {
            printf("Place %s: mu: %lf\n",placenames[p],mu_[p]);
        }
        printf("alfa: %lf, beta: %lf, wts: %d, gamma: %lf\n",alfa,beta,wts,gamma_); 
        printf("xibias: %lf\n",xibias);  
    }
    
}

void traces(void) {

    double min, mean, max;
    place p;
    color c;
    unsigned long pt; 
    double t;  

    fprintf(fo1,"##### Traces of number of balls\n");

    for(p=1; p<size; p++) {
        fprintf(fo1,"### Place %s\n",placenames[p]);
        fprintf(fo1,"### T\t\t\t\t conf RED \t\t\t\t\t conf GREEN \t\t\t\t\t conf BLUE");
        if(p == SE) fprintf(fo1,"\t\t\t\t\t conf BLACK"); 
        fprintf(fo1,"\t\t\t\t\t conf acc \t\t\t\t\t conf div \t\t\t\t redrun\n");
        for(pt = 1, t = DELTAT; t < MAXT && pt < MAXP; pt++, t+=DELTAT) {
            fprintf(fo1,"%le\t",t);
            for(c=0; c<nc; c++) {
                computeconf(NRUN, rnbc_[p][c][pt], &min, &mean, &max);
                fprintf(fo1,"%le\t%le\t%le\t",min,mean,max);
            }     
            if(p == SE) {
                computeconf(NRUN, rnbl_[pt], &min, &mean, &max);
                fprintf(fo1,"%le\t%le\t%le\t",min,mean,max);
            }        
            computeconf(NRUN, acc_[p][pt], &min, &mean, &max);
            fprintf(fo1,"%le\t%le\t%le\t",min,mean,max);
            computeconf(NRUN, div_[p][pt], &min, &mean, &max);
            fprintf(fo1,"%le\t%le\t%le\t",min,mean,max);
            if(redtot_[p][pt] > 0) fprintf(fo1,"%le\t%le\t",redrun_[p][pt]*1.0/redtot_[p][pt],redtot_[p][pt]); 
            else fprintf(fo1,"%le\t%le\t",-1.0,redtot_[p][pt]);
            fprintf(fo1,"\n");
        }
        fprintf(fo1,"\n\n");
    }    
}

void metrics(void) {

    double min, mean, max;
    unsigned long pt;  
    double t;  

    fprintf(fo2,"### interesting metrics\n");
    fprintf(fo2,"### T\t\t\t\t conf FQAN \t\t\t\t\t conf FQSN \t\t\t\t\t");
    fprintf(fo2," conf FIUA \t\t\t\t\t conf AIRI \t\t\t\t\t conf GPTA \t\t\t\t\t conf FWUA\n");
    for(pt = 1, t = DELTAT; t < MAXT && pt < MAXP; pt++, t+=DELTAT) {
        fprintf(fo2,"%le\t",t);
        computeconf(NRUN, rfqan_[pt], &min, &mean, &max);
        fprintf(fo2,"%le\t%le\t%le\t",min,mean,max);
        computeconf(NRUN, rfqsn_[pt], &min, &mean, &max);
        fprintf(fo2,"%le\t%le\t%le\t",min,mean,max);
        computeconf(NRUN, rfiua_[pt], &min, &mean, &max);
        fprintf(fo2,"%le\t%le\t%le\t",min,mean,max);
        computeconf(NRUN, rairi_[pt], &min, &mean, &max);
        fprintf(fo2,"%le\t%le\t%le\t",min,mean,max);
        computeconf(NRUN, rgpta_[pt], &min, &mean, &max);
        fprintf(fo2,"%le\t%le\t%le\t",min,mean,max);
        computeconf(NRUN, rfwua_[pt], &min, &mean, &max);
        fprintf(fo2,"%le\t%le\t%le\t",min,mean,max);
        fprintf(fo2,"\n");
    }    
}	

char* printball(Ballptr bptr) {
    color j;
    char* s = line;
    sprintf(s, "Ball: ");
    if(bptr->flow >= 0) sprintf(s+strlen(s),"flow %c, place %s, comp: [",flownames[bptr->flow],placenames[bptr->place]);
    else sprintf(s+strlen(s),"initially generated, comp: [");
    for(j=0; j<nc-1; j++) sprintf(s+strlen(s),"%.3f,",bptr->comp_[j]);
    sprintf(s+strlen(s), "%.3f], gpt: [",bptr->comp_[nc-1]);
    for(j=0; j<nc-1; j++) sprintf(s+strlen(s),"%.3f,",bptr->gpt_[j]);
    sprintf(s+strlen(s), "%.3f]",bptr->gpt_[nc-1]);
    return s;
}

void genball(place p, color c) {
    unsigned long cur = nb_[p];
    Ballptr bptr;
    color j;
    if(cur < MAXBALLSPERPLACE) {
        bptr = &b_[p][cur];
        bptr->flow = -1;
        bptr->place = p;
        for(j=0; j<nc; j++) {
            if(j == c) bptr->comp_[j] = 1.0;
            else bptr->comp_[j] = 0.0; 
            bptr->gpt_[j] = 0.0;  
        }
        nb_[p]++;
        if(DEBUG > 2) printf("Generated %s\n",printball(bptr));
    } else {
        printf("Exhausted b vector\n");
        exit(-1);
    }    
}

void copyball(flow f, place p, Ballptr bptr) {
    unsigned long cur = nb_[p];
    double u3,rep;
    color j;
    Ballptr copy;
    // compute average replication probability
    if(f != T) {
        rep = 0.0;
        for(j = 0; j < nc; j++) rep += r_[j][f] * bptr->comp_[j]; 
        u3 = uniform01(); 
        if(u3 < rep) {
            if(cur < MAXBALLSPERPLACE) {
                copy = &b_[p][cur];
                *copy = *bptr;
                copy->flow = f;
                copy->place = p;
                nb_[p]++;
                if(DEBUG>2) printf("Flow %c: copied ball into %s: %s\n",flownames[f],placenames[p],printball(copy));
            } else {
                printf("Exhausted b vector\n");
                exit(-1);
            }    
        }
    } else {
        for(int i = 1; i <= wts; i++) {
            if(cur < MAXBALLSPERPLACE) {
                copy = &b_[p][cur];
                *copy = *bptr;
                copy->flow = f;
                copy->place = p;
                nb_[p]++;
                cur = nb_[p];
                if(DEBUG>2) printf("Flow %c: copied ball into %s: %s\n",flownames[f],placenames[p],printball(copy));
            } else {
                printf("Exhausted b vector\n");
                exit(-1);
            }
        }        
    }   
}

int cmpfunc (const void * a, const void * b) {
    if(*(double*)a > *(double*)b) return 1;
    if(*(double*)a < *(double*)b) return -1;
    return 0;
}

Ball mixballs(int nmix) {
    int i;
    color j;
    Ball mix;
    double alfas2_[nmix];
        
    if(nmix > 2) {
        double alfas1_[nmix+1];
        for(i = 1; i < nmix; i++) alfas1_[i] = uniform01();
        qsort((void*)(&alfas1_[1]),nmix-1,sizeof(double),cmpfunc);
        alfas1_[0] = 0.0;
        alfas1_[nmix] = 1.0;
        for(i = 0; i < nmix; i++) alfas2_[i] = alfas1_[i+1]-alfas1_[i];
    } else {
        alfas2_[0] = uniform01();
        alfas2_[1] = 1.0-alfas2_[0];
    }     
    
    // compute mixture of balls
    for(j = 0; j < nc; j++) {
        mix.comp_[j] = 0.0;
        mix.gpt_[j] = 0.0;
        for(i = 0; i < nmix; i++) {
            if(j == 0 && DEBUG>2) printf("Mixing %s\n",printball(bptrarr_[i]));
            mix.comp_[j] += bptrarr_[i]->comp_[j] * alfas2_[i];
            mix.gpt_[j] += bptrarr_[i]->gpt_[j] * alfas2_[i];
        } 
    }
    mix.flow = -1;
    mix.place = 0;
    if(DEBUG>2) printf("mixballs: mixed ball is %s\n",printball(&mix));
    if(DEBUG > 8) {
        puts("Continue?");
        getchar();
    }
    return mix;
}

Ball demixballs(int nmix, double xibias) {
    color j;
    Ball mix;
    int bestid = -1;
    double bestq = -1;

    // compute mixture of balls: first pass
    for(j = 0; j < nc; j++) {
        mix.comp_[j] = 0.0;
        mix.gpt_[j] = 0.0;
    }
    // compute best ball
    for(int i = 0; i < nmix; i++) {
        double q = 0.0;
        if(DEBUG>2) printf("Demix: considering ball %s\n",printball(bptrarr_[i]));
        for(j = 0; j < nc; j++) q += bptrarr_[i]->comp_[j] * q_[j];
        if(q > bestq) {
            bestq = q;
            bestid = i;
        }
    }    
    if(DEBUG>2) printf("Demix: best ball is %s\n",printball(bptrarr_[bestid]));
    // biased mix
    for(int i = 0; i < nmix; i++) if (i != bestid) {
        for(j = 0; j < nc; j++) {
            mix.comp_[j] += bptrarr_[i]->comp_[j];
            mix.gpt_[j] += bptrarr_[i]->gpt_[j];
        }     
    } 
    for(j = 0; j < nc; j++) {
        mix.comp_[j] *= (1-xibias)*1.0/(nmix-1);
        mix.gpt_[j] *= (1-xibias)*1.0/(nmix-1);
        mix.comp_[j] += (xibias)*bptrarr_[bestid]->comp_[j];
        mix.gpt_[j] += (xibias)*bptrarr_[bestid]->gpt_[j];
    }
    if(DEBUG>2) printf("demixballs: mix is %s\n",printball(&mix));
    if(DEBUG > 8) {
        puts("Continue?");
        getchar();
    }
    return mix;
}

void printstate(double t) {
    place p; 
    printf("%lf: system state: ",t);
    for(p=0; p<size; p++) {
        if(DEBUG > 10) {
            printf("Place %s:\n",placenames[p]);
            for(unsigned int i = 0; i < nb_[p]; i++) printf("\t%s\n",printball(&b_[p][i]));
        } else printf("%s: %ld - ",placenames[p],nb_[p]);
    }
    printf("\n");    
}
        
void run(int r) {

    double t = 0.0;
    double dt,u1,u3;
    unsigned long u2;
    unsigned long pt = 1; // next print point 
    double npt = DELTAT; // next print time 
    double Lambda, LamA, LamS, Mu, Mubl, LamAN, LamSN;
    color c;
    place p,sp;
    double th,lamth;
    Ballptr bptr;
    Ball b;
    unsigned long totb; // for debug 
    unsigned long noanswersA = 0L;
    unsigned long noanswersS = 0L;
    unsigned long questionsA = 0L;
    unsigned long questionsS = 0L;
    double num,numi,numia,den,tmp;
    double ps_[nc];

    printf("**** RUN %d ****\n",r);

    // initial condition
    for(p=0; p<size; p++) nb_[p] = 0L;
    for(c=0; c<nc; c++) {
        genball(GE,c);
        for(p=1; p<size; p++) {
            for(int i=0; i<INITBALLS; i++) genball(p,c); 
        }    
    }    
    nbl = INITBLACKS; 
    
    // main loop    
    while (t < MAXT) {

        /* compute total number of balls */
        if(DEBUG>1) {
            totb = 0;
            for(p=0; p<size; p++) totb += nb_[p];   
            totb += nbl;
            printf("\n%lf: total number of balls %ld\n",t,totb);
        }    

        /* compute variable lambdas */
        Lam_[I] = nb_[WWW] * lam_[I];
        Lam_[T] = nb_[TS] * lam_[T];
        if(nb_[GPT] + nb_[SE] + nbl > 0) { LamA = Lam_[A]; LamAN = 0.0; } else { LamA = 0.0; LamAN = Lam_[A]; }
        if(nb_[SE] + nbl > 0) { LamS = Lam_[S]; LamSN = 0.0; } else { LamS = 0.0; LamSN = Lam_[S]; }

        /* compute aggregate death rates */
        Mu = 0.0;
        for(p=0; p<size; p++) {
            Mu_[p] = mu_[p] * nb_[p];
            Mu += Mu_[p];    
        }
        
        /* death rate of black balls */
        Mubl = nbl * mu_[SE];

        if(DEBUG>1) {
            printf("Lam_[H]: %lf\n",Lam_[H]);
            printf("Lam_[P]: %lf\n",Lam_[P]);
            printf("Lam_[I]: %lf\n",Lam_[I]);
            printf("Lam_[T]: %lf\n",Lam_[T]);
            printf("LamA: %lf\n",LamA);
            printf("LamS: %lf\n",LamS);
            printf("LamAN: %lf\n",LamAN);
            printf("LamSN: %lf\n",LamSN);
            printf("Mu: %lf\n",Mu);
            printf("Mubl: %lf\n",Mubl);
        }

        Lambda = Lam_[H] + Lam_[P] + Lam_[I] + Lam_[T] + LamA + LamS + Mu + Mubl + LamAN + LamSN;
        dt = exponential(1.0/Lambda);
        t += dt;
        if(DEBUG>2) printf("dt = %lf, newt %lf, npt %lf, pt %ld, Lambda = %lf\n",dt,t,npt,pt,Lambda);

        
        while(t > npt && pt < MAXP) {

            // record results 
            for(p=1; p<size; p++) {
                rnb_[p][pt][r] = nb_[p]; 
                for(c = 0; c < nc; c++) {
                    rnbc_[p][c][pt][r] = 0.0;
                    for(unsigned long i = 0; i < nb_[p]; i++) rnbc_[p][c][pt][r] += b_[p][i].comp_[c];
                }
            }    
            rnbl_[pt][r] = nbl;

            // accuracy
            for(p=1; p<size; p++) {
                num = 0;
                den = 0;
                for(unsigned long i = 0; i < nb_[p]; i++) {
                    b = b_[p][i];
                    for(c = 0; c < nc; c++) {
                        tmp = b.comp_[c];
                        if(tmp > 0) {
                            if(q_[c] < theta_) {
                                num += tmp;
                            }    
                            den += tmp;
                        }
                    }
                }
                if(den > 0) {
                    u1 = num/den;
                    if(p != WWW) {
                        if(u1 > 0.5) redrun_[p][pt]++;
                        redtot_[p][pt]++;
                    }    
                    acc_[p][pt][r] = u1;
                } else acc_[p][pt][r] = -1; 
            }
            // diversity
            for(p=1; p<size; p++) {
                for(c = 0; c < nc; c++) if(q_[c] > theta_) ps_[c] = 0;
                den = 0;
                for(unsigned long i = 0; i < nb_[p]; i++) {
                    b = b_[p][i];
                    for(c = 0; c < nc; c++) if(q_[c] > theta_) {
                        tmp = b.comp_[c];
                        ps_[c] += tmp;
                        den += tmp;
                    }
                }
                tmp = 0.0;
                for(c = 0; c < nc; c++) if(q_[c] > theta_) {
                    tmp += fabs(p_[c]-ps_[c]/den);
                }
                div_[p][pt][r] = 1.0-0.5*tmp;            
            }

            // record metrics

            // rfwua
            num = 0.0;
            den = nb_[WWW];
            for(unsigned long i = 0; i < nb_[WWW]; i++) {
                b = b_[WWW][i];
                if(b.flow == U || b.flow == S) num++;
            }    
            if(den > 0) rfwua_[pt][r] = num/den; else rfwua_[pt][r] = -1;

            // noanswers
            if (questionsA > 0) rfqan_[pt][r] = noanswersA*1.0/questionsA; else rfqan_[pt][r] = -1;
            if (questionsS > 0) rfqsn_[pt][r] = noanswersS*1.0/questionsS; else rfqsn_[pt][r] = -1;
            
            // fraction of irrelevant used answers and AI responsability index
            numi = numia = den = 0.0;
            for(unsigned long i = 0; i < nb_[WWW]; i++) {
                b = b_[WWW][i];
                if(b.flow == U || b.flow == S) { 
                    for(c = 0; c < nc; c++) {
                        tmp = b.comp_[c];
                        if(tmp > 0) {
                            if(q_[c] < theta_) {
                                if (b.flow == U) numia += tmp;
                                numi += tmp;
                            }    
                            den += tmp;
                        }
                    }
                }    
            }
            if(den > 0) {
                rfiua_[pt][r] = numi/den;
                if(numi/den > 0.5) redrun_[WWW][pt]++;
                redtot_[WWW][pt]++;
            } else rfiua_[pt][r] = -1;
            if(numi > 0) rairi_[pt][r] = numia/numi; else rairi_[pt][r] = -1;
            // GPT autophagy index
            num = den = 0.0;
            for(unsigned long i = 0; i < nb_[GPT]; i++) {
                for(c = 0; c < nc; c++) {
                    b = b_[GPT][i];
                    num += b.gpt_[c];   
                    den += b.comp_[c];
                }    
            }
            if(den > 0) rgpta_[pt][r] = num/den; else rgpta_[pt][r] = -1;

            npt += DELTAT;
            pt++;
        }

        // consider possible events
        u1 = uniform(0,Lambda);  
        lamth = 0.0;

        // addition to training set
        if(u1 < (lamth += Lam_[H])) {
            u3 = uniform01();
            if(u3 < gamma_) {
                u2 = uniformint(0,nb_[GE]-1);
                bptr = &b_[GE][u2];
                copyball(H,TS,bptr);
            } else {
                if(nb_[SE] > 0) {
                    // printf("trying to call uniformint %ld, %ld\n",0L,nb_[SE] + nbl - 1);
                    u2 = uniformint(0,nb_[SE] + nbl - 1);
                    bptr = NULL;
                    if(u2 < nb_[SE]) bptr = &b_[SE][u2]; 
                    if(bptr != NULL) {
                        copyball(H,TS,bptr);
                    }    
                }       
            }    
        } else if(u1 < (lamth += Lam_[P])) {
        // addition to WWW
            u2 = uniformint(0,nb_[GE]-1);
            bptr = &b_[GE][u2];
            copyball(P,WWW,bptr);
        } else if(u1 < (lamth += Lam_[I])) {
        // indexing event
            u2 = uniformint(0,nb_[WWW]-1);
            bptr = &b_[WWW][u2];
            copyball(I,SE,bptr);  
        } else if(u1 < (lamth += Lam_[T])) { 
        // training event
            u2 = uniformint(0,nb_[TS]-1);
            bptr = &b_[TS][u2];
            copyball(T,GPT,bptr);         
        } else if(u1 < (lamth += LamA)) {
        
            // Question event A
            questionsA++;
            _Bool selectedblack = false; 
            
            // first ball, determining if answer is actually produced
            u2 = uniformint(0,nb_[GPT] + nb_[SE] + nbl - 1);
            if(u2 < nb_[GPT]) bptrarr_[0] = &b_[GPT][u2];
            else if(u2 < nb_[GPT] + nb_[SE]) bptrarr_[0] = &b_[SE][u2-nb_[GPT]];
            else selectedblack = true;
            
            if (selectedblack == true) {
                if(DEBUG>2) printf("Black answer!\n");  
                noanswersA++;  
            } else {
                
                if(MIXGAI && NUMBALLSTOMIX > 1) {
                    for(int i = 1; i < NUMBALLSTOMIX; i++) {
                        u2 = uniformint(0,nb_[GPT] + nb_[SE] - 1);  
                        if(u2 < nb_[GPT]) bptrarr_[i] = &b_[GPT][u2];
                        else bptrarr_[i] = &b_[SE][u2-nb_[GPT]];
                    }
                    b = mixballs(NUMBALLSTOMIX);
                    bptr = &b;
                } else bptr = bptrarr_[0];
                
                // flow F
                u3 = uniform01();
                if(u3 < beta) {
                    copyball(F,GPT,bptr);
                }
                
                // flow U
                u3 = uniform01();
                if(u3 < alfa) {
                    // flow U is autophagous
                    for(color j = 0; j < nc; j++) bptr->gpt_[j] = bptr->comp_[j];
                    copyball(U,WWW,bptr);
                }                          
            }  
            
        } else if(u1 < (lamth += LamS)) {
            
            // Question event S  
            questionsS++;
            u2 = uniformint(0,nb_[SE] + nbl - 1);
            if(u2 < nb_[SE]) {
                bptrarr_[0] = &b_[SE][u2];
                if(MIXSE && NUMBALLSTOMIX > 1) {
                    for(int i = 1; i < NUMBALLSTOMIX; i++) {
                        u2 = uniformint(0,nb_[SE] - 1);  
                        bptrarr_[i] = &b_[SE][u2];
                    }
                    b = demixballs(NUMBALLSTOMIX,xibias);
                    bptr = &b;
                } else bptr = bptrarr_[0];    
                         
                // flow S
                copyball(S,WWW,bptr); 

            } else {
                if(DEBUG>2) printf("Black answer!\n");
                noanswersS++;
            }       

        } else if(u1 < (lamth += Mu)) {
        
        // death event of colored ball
            
            // determine place
            u3 = uniform(0,Mu); 
            sp = -1;
            th = Mu_[0];
            for(p=1; p<size; p++) {
                if(u3 > th && u3 < th + Mu_[p]) {
                    sp = p;
                    break;
                } else th += Mu_[p];
            }
            assert(sp >= 1);

            // determine ball
            u2 = uniformint(0,nb_[sp]-1);
            bptr = &b_[sp][u2];
            if(DEBUG>2) printf("Death of ball at place %s: %s\n",placenames[sp],printball(bptr));
            if(nb_[sp] > 1) {
                *bptr = b_[sp][nb_[sp]-1]; // replace death ball with last ball, if any!
            }
            nb_[sp]--;
        
        } else if (u1 < (lamth += Mubl)) {
            // death event of black ball
            assert(nbl > 0);
            nbl--;
            if(DEBUG>2) printf("%lf: Death of black ball at place %s\n",t,placenames[SE]);     

        } else if(u1 < (lamth += LamAN)) {
            questionsA++;
            noanswersA++;
        } else {
            questionsS++;
            noanswersS++;    
        }

        if(DEBUG > 1) printstate(t);
        if(DEBUG > 10) {
            puts("Continue?");
            getchar();
        }   
    }    
}

int main(int argc, char **argv) {

    init(sce); 
    rplot();

    /* init redrun stuff */
    for(place p=0; p<size; p++) for(unsigned long pt = 0; pt < MAXP; pt++) {
        redrun_[p][pt] = 0.0; // fraction of runs in which majority of answers are irrelevant
        redtot_[p][pt] = 0.0; // denom for redrun
    }   
	
    if(1) for(int r=1; r <= NRUN; r++) {
        run(r);
    }

    traces();
    fclose(fo1);

    metrics();
    fclose(fo2);

    return 0;
}

void computeconf(int n, double *b, double *min, double *mean, double *max)
{
    double tmp, quantile;
    int i;
    double var = 0.0;
    double ave = 0.0;
    int numvalid = 0;
        
    ave = 0.;
    for (i = 1; i<=n; i++) {
        if(b[i] >= 0) { 
            ave += b[i];
            numvalid++;
        }
    }
    if(numvalid > 1) {        
        ave /= numvalid;
        for (i = 1; i<=n; i++) {
            if(b[i] >= 0) var += ((b[i]-ave)*(b[i]-ave));
        }    
        var = sqrt(var/(numvalid-1));
        quantile = gsl_cdf_tdist_Qinv(0.025,numvalid);
        if(DEBUG > 10) printf("Quantile for n = %d, 95 level: %lf\n",numvalid,quantile);
        tmp = quantile/sqrt(numvalid);
    } else {
        ave = var = tmp = 0.0;       
    }    
    (*min) = ave - var*tmp;
    (*max) = ave + var*tmp;
    (*mean) = ave;
}

