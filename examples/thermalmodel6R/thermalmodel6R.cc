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

const size_t sDIM = 6; /* System dimension */
const size_t iDIM = 2; /* Input dimension */
const double T = 25; /* Sampling time */
const size_t N = 6; /* Temporal Horizon */

typedef std::array<double,sDIM> state_type;

/* forward declaration of the ode solver */
template<class F>
void ode_solver(F rhs, state_type &x, size_t nint, double h);

/* ODE */
auto  system_post = [](state_type &x, double* u) -> void {

    /* ode describing six-room thermal model*/
    auto rhs=[u](state_type &xx, const state_type &x) -> void {
        const double a=0.05;
        const double ae1=0.005;
        const double ae4=0.005;
        const double ae=0.0033;
        const double ah=0.0036;
        const double te=12;
        const double th=100;
    xx[0] = (-3*a-ae1-ah*u[0])*x[0]+a*x[1]+a*x[2]+a*x[4]+ae1*te+ah*th*u[0];
    xx[1] = (-2*a-ae)*x[1]+a*x[0]+a*x[3]+ae*te;
    xx[2] = (-2*a-ae)*x[2]+a*x[0]+a*x[3]+ae*te;
    xx[3] = (-3*a-ae4-ah*u[1])*x[3]+a*x[1]+a*x[2]+a*x[5]+ae4*te+ah*th*u[1];
    xx[4] = (-a-ae)*x[4]+a*x[0]+ae*te;
    xx[5] = (-a-ae)*x[5]+a*x[3]+ae*te;
    };
    size_t nint = 5; /* no. of time step for ode_solving */
    double h=T/nint; /* time step for ode solving (T is an sampling time) */
    ode_solver(rhs,x,nint,h); /* Runga Kutte solver */
};

/* defining target set for the controller
   ul[i] : upper bound on target temperature in ith room
   ll[i] : upper bound on target temperature in ith room */
auto setBounds = [](state_type y) -> bool {
    double ul=21.0;
    double ll=17.5;
    bool s = true;
        for(size_t j = 0; j < sDIM; j++)
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
   xs[0]=17;xs[1]=17;xs[2]=17;xs[3]=17;xs[4]=17;xs[5]=17;

/************************************************************
    input information
***********************************************************/
    /* lower bounds on inputs */
    double lb[sDIM]={0,0};
    /* upper bounds on inputs */
    double ub[sDIM]={1,1};
    /* quantization parameter */
    double eta[sDIM]={.5,1};

/******************************************************
     Symbolic model construction
******************************************************/
    /* defining Symbolic Set for the transition relation */
    SymbolicSetSpace ss(ddmgr,iDIM,lb,ub,eta,N);
    tt.tic();
    ss.addAbsStates();

    /* defining abstraction class */
    getAbstraction<state_type> ab(&ss);
    tt.toc();
    std::cout<<"Number of Elements in the Transition Relation : "<<ab.getSize()<<std::endl;
    
    /* Computing the constrain set */
    std::cout<<"Computing the constrain set ... ";
    tt.tic();
    BDD set = ab.getAbstractSet(system_post,setBounds,xs,iDIM,lb,ub,eta);
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
    
    /* controller for reach and stay specification */
    C = fp.reachStay(set);
    
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
    x0[0]=16;x0[1]=16;x0[2]=16;x0[3]=16;x0[4]=16;x0[5]=16;

    /* approximation on abstraction */
    double epsilon = 0.8;

    /* number of time samples the controller should run */
    size_t t = 40;

    /* find out mode sequence which is very close to the output */
    std::cout<<"Computed the initial state is...\n";
    
    /* finding initial state in abstraction for safety specification */
    BDD w0 = ab.getAbsState(system_post,x0,xs,epsilon,sDIM);
    
    /* closed-loop simulation */
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
