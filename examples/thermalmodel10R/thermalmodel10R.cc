#include <array>
#include <iostream>
#include<cmath>
#include "cuddObj.hh"

#include "../../src/TicToc.hh"
#include "../../src/fixedPointMode.hh"
#include "../../src/getAbstraction.hh"
#include "../../src/SymbolicSetSpace.hh"
#include "../../src/abstractionMode.hh"

using namespace std;

const int sDIM = 10; /* System dimension */
const int iDIM = 2; /* Input dimension */
const double T = 25; /* Sampling time */
size_t N = 12; /* Temporal Horizon */

typedef std::array<double,sDIM> state_type;

/* forward declaration of the ode solver */
template<class F>
void ode_solver(F rhs, state_type &x, size_t nint, double h);

/* ODE */
auto  system_post = [](state_type &x, double* u) -> void {
    /* ode describing thermal model*/
    auto rhs=[u](state_type &xx, const state_type &x) -> void {
        const double a=0.05;
        const double ae2=0.005;
        const double ae5=0.005;
        const double ae=0.0033;
        const double ah=0.0036;
        const double te=12;
        const double th=100;
    xx[0]=(-a-ae)*x[0]+a*x[1]+ae*te;
    xx[1]=(-4*a-ae2-ah*u[0])*x[1]+a*x[0]+a*x[6]+a*x[8]+a*x[2]+ae2*te+ah*th*u[0];
    xx[2]=(-2*a-ae)*x[2]+a*x[1]+a*x[3]+ae*te;
    xx[3]=(-2*a-ae)*x[3]+a*x[2]+a*x[4]+ae*te;
    xx[4]=(-4*a-ae5-ah*u[1])*x[4]+a*x[3]+a*x[7]+a*x[5]+a*x[9]+ae5*te+ah*th*u[1];
    xx[5]=(-a-ae)*x[5]+a*x[4]+ae*te;
    xx[6]=(-a-ae)*x[6]+a*x[1]+ae*te;
    xx[7]=(-a-ae)*x[7]+a*x[4]+ae*te;
    xx[8]=(-a-ae)*x[8]+a*x[1]+ae*te;
    xx[9]=(-a-ae)*x[9]+a*x[4]+ae*te;
    };
    size_t nint = 5; /* no. of time step for ode solving */
    double h=T/nint; /* time step for ode solving (T is the sampling time) */
    ode_solver(rhs,x,nint,h); /* Runga Kutte solver */
};

/* defining safe set for the controller
	ul[i] : upper bound on the temperature in ith room
	ll[i] : lower bound on the temperature in ith room */
auto setBounds = [](state_type y) -> bool {
    double ul=21.8;
    double ll=18;
    bool s = true;
        for(int j = 0; j < sDIM; j++)
        {
            if( y[j] >= ul || y[j] <= ll )
            {
                s = false;
                break;
            }
        }
    return s;
};


int  main(){
/*to measure time */
   TicToc tt;

/* CUDD Manager */
   Cudd ddmgr;

/* source state of the system */
   state_type xs;
   xs[0]=17;xs[1]=17;xs[2]=17;xs[3]=17;xs[4]=17;xs[5]=17;xs[6]=17;xs[7]=17;xs[8]=17;xs[9]=17;
    
/************************************************************
     input information
*************************************************************/
    const size_t P = 3;  /* Number of elements in input set*/
    double ud[P][iDIM]={{0,0},{0,1},{1,0}};
    double *UD[P];
    for(size_t i=0;i<P;i++){
        UD[i]=ud[i];
    }
    
/******************************************************
     Symbolic model construction
******************************************************/
    /* defining Symbolic Set for the transition relation */
    SymbolicSetSpace ss(ddmgr,P,N);  
    tt.tic();
    ss.addAbsStates();
    
    /* defining abstraction class */
    getAbstraction<state_type> ab(&ss);
    tt.toc();
    std::cout<<"Number of Elements in the Transition Relation : "<<ab.getSize()<<std::endl;
    
    /* Computing the constrain set */
    std::cout<<"Computing the constrain set ... ";
    tt.tic();
    BDD set = ab.getAbstractSet(system_post,setBounds,xs,iDIM,P,UD);
    tt.toc();
    if(set == ddmgr.bddZero()){
        std::cout << "Set is empty !!";
        return 0;
    }

/******************************************************
     synthesizing the controller
******************************************************/
    std::cout<<"Computing the controller ... ";
    fixedPointMode fp(&ab);
    BDD C;
    tt.tic();
    
    /* controller for safety specification */
    C = fp.safe(set);
    tt.toc();

    if(C == ddmgr.bddZero()){
        std::cout << "No controller is found !!";
        return 0;
    }
    
/******************************************************
     closed-loop simulation
******************************************************/
    /* initial state of the system */
    state_type x0;
    x0[0]=19;x0[1]=19;x0[2]=19;x0[3]=19.8;x0[4]=20.8;x0[5]=19.8;x0[6]=19;x0[7]=19.8;x0[8]=19;x0[9]=19.8;
    
    /* approximation on abstraction */
    double epsilon = 0.1;
    
    /* number of time samples the controller should run */
    size_t t = 20;
    /* find out mode sequence which is very close to the output */
    std::cout<<"Computed initial state is...\n";
    
    /* finding initial state in abstraction for safety specification */
    BDD w0 = ab.getAbsState(system_post,set,x0,xs,epsilon,sDIM);
    
    /* open-loop simulation */
    ab.closedLoopSim(C,w0,system_post,x0,sDIM,t);
    
    return 0;

}

template<class F>
void ode_solver(F rhs, state_type &x, size_t nint, double h) {
    /* runge kutte order 4 */
    state_type k[4];
    state_type tmp;
    
    for(size_t t=0; t<nint; t++) {
        rhs(k[0],x);
        for(size_t i=0;i<sDIM;i++)
            tmp[i]=x[i]+h/2*k[0][i];
        
        rhs(k[1],tmp);
        for(size_t i=0;i<sDIM;i++)
            tmp[i]=x[i]+h/2*k[1][i];
        
        rhs(k[2],tmp);
        for(size_t i=0;i<sDIM;i++)
            tmp[i]=x[i]+h*k[2][i];
        
        rhs(k[3],tmp);
        for(size_t i=0; i<sDIM; i++)
            x[i] = x[i] + (h/6)*(k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i]);
    }
}
