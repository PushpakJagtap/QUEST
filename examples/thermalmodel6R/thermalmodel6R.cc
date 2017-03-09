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

const int sDIM = 6; /* System dimension */
const int iDIM = 2; /* Input dimension */
const double T = 25; /* Sampling time */
size_t P = 3; /* Cardinality of (quantized) input set */
size_t N = 12; /* Temporal Horizon */

typedef std::array<double,sDIM> state_type;

/* forward declaration of the ode solver */
template<class F>
void ode_solver(F rhs, state_type &x, size_t nint, double h);

/* ODE */
auto  system_post = [](state_type &x, int u) -> void {
    
    /* Defining system inputs to system*/
    double ub[iDIM]={0};
    if(u==0)
    {ub[0]=0;ub[1]=0;}
    if(u==1)
    {ub[0]=0;ub[1]=1;}
    if(u==2)
    {ub[0]=1;ub[1]=0;}
    
    /* ode describing thermal model*/
    auto rhs=[ub](state_type &xx, const state_type &x) -> void {
        const double a=0.05;
        const double ae1=0.005;
        const double ae4=0.005;
        const double ae=0.0033;
        const double ah=0.0036;
        const double te=12;
        const double th=100;
    xx[0] = (-3*a-ae1-ah*ub[0])*x[0]+a*x[1]+a*x[2]+a*x[4]+ae1*te+ah*th*ub[0];
    xx[1] = (-2*a-ae)*x[1]+a*x[0]+a*x[3]+ae*te;
    xx[2] = (-2*a-ae)*x[2]+a*x[0]+a*x[3]+ae*te;
    xx[3] = (-3*a-ae4-ah*ub[1])*x[3]+a*x[1]+a*x[2]+a*x[5]+ae4*te+ah*th*ub[1];
    xx[4] = (-a-ae)*x[4]+a*x[0]+ae*te;
    xx[5] = (-a-ae)*x[5]+a*x[3]+ae*te;
    };
    size_t nint = 5; /* no. of time step for ode_solving */
    double h=T/nint;
    ode_solver(rhs,x,nint,h);
};

/* defining constrains for the controller 
   ul[i] : upper limit of the ith dimension 
   ll[i] : lower limit of the ith dimension */
auto constrain = [](state_type y) -> bool {
    double ul=21.0;
    double ll=17.5;
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
   xs[0]=17;xs[1]=17;xs[2]=17;xs[3]=17;xs[4]=17;xs[5]=17;
    
/* initial state of the system */
    state_type x0;
    x0[0]=16;x0[1]=16;x0[2]=16;x0[3]=16;x0[4]=16;x0[5]=16;

    /* approximation on abstraction */
    double epsilon = 0.2;
    
    /* number of time samples the controller should run */
    size_t t = 40;
    
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
    BDD set = ab.getConstrainSet(system_post,constrain,xs);
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
    
    /* find out mode sequence which is very close to the output */
    std::cout<<"Computed the initial state is...\n";
    
    /* finding inital stae in abstraction for safety specification */
    BDD w0 = ab.getAbsState(system_post,x0,xs,epsilon,sDIM);
    
    /* open-loop simulation */
    ab.openLoopSim(C,w0,system_post,x0,sDIM,t);
    
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
