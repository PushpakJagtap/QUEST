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

const int sDIM = 5; /* System dimension */
const int iDIM = 5; /* Input dimension */
size_t P = 3; /* Cardinality of (quantized) input set */
size_t N = 12; /* Temporal Horizon */

typedef std::array<double,sDIM> state_type;

/* ODE */
auto  system_post = [](state_type &x, double* u) -> void {
   state_type z = x;
  /* the ode describing the system */
    const double v = 70 ; 
    const double l = 0.25 ;
    const double T = 0.002777;
    const double q = 0.25;
    x[0] = (1-T*v/l)*z[0] + u[0];
    x[1] = T*v/l*z[0] + (1-T*v/l-q) *z[1] + u[1];
    x[2] = T*v/l*z[1] + (1-T*v/l) *z[2] + u[2];
    x[3] = T*v/l*z[2] + (1-T*v/l) *z[3] + u[3];
    x[4] = T*v/l*z[3] + (1-T*v/l) *z[4] + u[4];
  
};

/* defining safety bounds for the controller
   ul[i] : upper limit of the ith dimension 
   ll[i] : lower limit of the ith dimension */
auto setBounds = [](state_type y) -> bool {
    state_type ul;
    state_type ll;
    for(int i = 0; i < sDIM; i++)
    {
       ul[i] = 15.00;
       ll[i] = 0.00;
    }
    bool s = true;
    for(int i = 0; i < sDIM; i++)
    {
	if( y[i] >= ul[i] || y[i] <= ll[i] )
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
   xs[0] = 3.8570;
   xs[1] = 3.3750;
   xs[2] = 3.3750;
   xs[3] = 8.5177;
   xs[4] = 8.5177;

 /************************************************************
        input information
 *************************************************************/
       const size_t P = 3;  /* Number of elements in input set*/
       double ud[P][iDIM]={{6,0,8,0,0},{6,0,0,0,0},{0,0,8,0,0}};
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
    C = fp.safe(set);
    tt.toc();

    if(C == ddmgr.bddZero()){
        std::cout << "No controller is found !!";
        return 0;
    }
    
/******************************************************
 *  closed-loop simulation
******************************************************/
    /* initial state of the system */
       state_type x0;
       x0[0] = 1.41783;
       x0[1] = 4.92788;
       x0[2] = 10.6711;
       x0[3] = 9.53216;
       x0[4] = 14.5308;
    /* approximation on abstraction */
       double epsilon = 0.00085;
    /* number of times the controller should run */
       size_t t = 14;
    std::cout<<"Computing the initial state ... ";
    
    BDD w0 = ab.getAbsState(system_post,set,x0,xs,epsilon,sDIM);

   /* to get output */
   ab.closedLoopSim(C,w0,system_post,xs,sDIM,t);

   return 0;

}

