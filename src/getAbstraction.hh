#ifndef GETABSTRACTION_HH_
#define GETABSTRACTION_HH_

#include <iostream>
#include <stdexcept>
#include <fstream>

#include "TicToc.hh"

#include "cuddObj.hh"
#include "cudd.h"
#include "dddmp.h"


#include "SymbolicSetSpace.hh"
#include "abstractionMode.hh"

template<class state_type> class getAbstraction: public abstractionMode {

	public:
	using abstractionMode::abstractionMode; 
	template<class F1,class F2>
	BDD getConstrainSet(F1 &system_post,F2 &constrain,state_type x){
	BDD set = mgr_->bddZero();
	const int* minterm;
	size_t t = (nvars/shift-1)/2;
	size_t np = t * shift;
	int* phase = new int [np];
	for(size_t i = 0; i < np; i++)
		phase[i] = 0;
	DdManager *mgr = mgr_->getManager();
	DdNode**  dvars = new DdNode*[np];
	int* element = new int [t];
        for(size_t i = 0; i < np; i++)
      		dvars[i]=Cudd_bddIthVar(mgr,vars_[i]);
	for(begin(transitionRelation_,nvars); !done(); next())
	{
		//begin(transitionRelation_);
		progress();
		minterm = currentMinterm();
		for(size_t k = 0,i =0; i < t; i++) 
		{
		      	size_t idx=0;
		      	for(size_t c = 1, j = 0; j < shift; c*=2,k++, j++)
				idx+=minterm[vars_[k]]*c;
		      	element[i]=idx;
    		}
		state_type y = getOutput(system_post,element,x);
		if(constrain(y))
		{ 
			for(size_t i = 0; i < np; i++)
         			phase[i]=minterm[vars_[i]];			
			BDD current(*mgr_,Cudd_bddComputeCube(mgr,dvars,phase,np));
			set +=current;
		} 
	}
	std::cout<<std::endl;
	std::cout<<" Number of elements in the constrain set : "; 
	std::cout<<set.CountMinterm(np)<<std::endl;
	std::cout<<std::endl;
	return set;
} 



	template<class F1>	
	state_type getOutput(F1 &system_post, int* element, state_type x ){
	state_type y = x;
	size_t t = (nvars/shift-1)/2;
	for(size_t i = 0; i < t; i++)
	{
		system_post(y,element[i]);
	}
	return y;
}
	template<class F1>
	BDD getAbsState(F1 &system_post, state_type x,state_type xs,double eta, int sDIM){
	BDD absState = mgr_->bddZero();
	const int* minterm;
	size_t t = (nvars/shift-1)/2;
	size_t np = t * shift;
	int* phase = new int [np];
	for(size_t i = 0; i < np; i++)
		phase[i] = 0;
	DdManager *mgr = mgr_->getManager();
	DdNode**  dvars = new DdNode*[np];
	int* element = new int [t];
        for(size_t i = 0; i < np; i++)
      		dvars[i]=Cudd_bddIthVar(mgr,vars_[i]);
	for(begin(mgr_->bddOne(),np); !done(); next())
	{
		//begin(transitionRelation_);
		//progress();
		minterm = currentMinterm();
		for(size_t k = 0,i =0; i < t; i++) 
		{
		      	size_t idx=0;
		      	for(size_t c = 1, j = 0; j < shift; c*=2,k++, j++)
				idx+=minterm[vars_[k]]*c;
		      	element[i]=idx;
    		}
		state_type y = getOutput(system_post,element,xs);
		bool s = true;
		for(int i = 0; i < sDIM; i++)
		{
            if(std::abs(y[i]-x[i]) > eta)
			{
				s = false;
				break;
			}
		}
		if(s)
		{ 
			for(size_t i = 0; i < np; i++)
         			phase[i]=minterm[vars_[i]];			
			BDD current(*mgr_,Cudd_bddComputeCube(mgr,dvars,phase,np));
			absState +=current;
			break;
		} 
	}
	return absState;	
}
    template<class F1>
    BDD getAbsState(F1 &system_post, BDD set, state_type x,state_type xs,double eta, int sDIM){
        BDD absState = mgr_->bddZero();
        const int* minterm;
        size_t t = (nvars/shift-1)/2;
        size_t np = t * shift;
        int* phase = new int [np];
        for(size_t i = 0; i < np; i++)
            phase[i] = 0;
        DdManager *mgr = mgr_->getManager();
        DdNode**  dvars = new DdNode*[np];
        int* element = new int [t];
        for(size_t i = 0; i < np; i++)
            dvars[i]=Cudd_bddIthVar(mgr,vars_[i]);
        for(begin(set,np); !done(); next())
        {
            //begin(transitionRelation_);
            //progress();
            minterm = currentMinterm();
            for(size_t k = 0,i =0; i < t; i++)
            {
                size_t idx=0;
                for(size_t c = 1, j = 0; j < shift; c*=2,k++, j++)
                    idx+=minterm[vars_[k]]*c;
                element[i]=idx;
            }
            state_type y = getOutput(system_post,element,xs);
            bool s = true;
            for(int i = 0; i < sDIM; i++)
            {
                if(abs(y[i]-x[i]) > eta)
                {
                    s = false;
                    break;
                }
            }
            if(s)
            { 
                for(size_t i = 0; i < np; i++)
                    phase[i]=minterm[vars_[i]];			
                BDD current(*mgr_,Cudd_bddComputeCube(mgr,dvars,phase,np));
                absState +=current;
                break;
            } 
        }
        return absState;	
    }
    
    template<class F1>
    void openLoopSim(BDD C, BDD D, F1 &system_post, state_type y, int sDIM, int t){
	std::ofstream outfile;
	outfile.open("controller.txt");
	std::vector<int> ms;	
	for(int counter = 0; counter < t; counter++)
	{	
		BDD temp = C & D;
		//temp.print(2,2);
		const int* minterm;
		size_t t = (nvars/shift-1)/2;
		size_t np = t * shift;
		int nextmode = 0;
		int* element = new int [t];
		begin(temp,np+shift);
		minterm = currentMinterm();
		/*for(size_t i = 0; i < np+shift; i++)
		std::cout<<minterm[i]<<" ";
		std::cout<<std::endl;*/
		int* phase = new int [np];
		for(size_t i = 0; i < np; i++)
			phase[i] = 0;
		DdManager *mgr = mgr_->getManager();
		DdNode**  dvars = new DdNode*[np];
		for(size_t i = 0; i < np; i++)
	   		dvars[i]=Cudd_bddIthVar(mgr,vars_[i]);
		size_t k = 0;
		for(size_t i =0; i < t; i++) 
		{
		      	size_t idx=0;
		      	for(size_t c = 1, j = 0; j < shift; c*=2,k++, j++)
				idx+=minterm[vars_[k]]*c;
		      	element[i]=idx;
	    	}
		std::cout<<"Present input sequence : "<<std::endl;
		for(size_t i = 0; i < t; i++)
		std::cout<<element[i]<<" ";
		std::cout<<std::endl;
		for(size_t c = 1, j = 0; j < shift; c*=2,k++, j++)
			nextmode+=minterm[vars_[k]]*c;
		ms.push_back(nextmode);
		std::cout<<"Next input = "<< nextmode<<std::endl;
		y = getOutput(system_post,element,y);
		std::cout<<"Present states : "<<std::endl;
		for(int i = 0 ; i < sDIM; i++)
		std::cout<<y[i]<<" ";
		std::cout<<std::endl;
		for(size_t i = 0; i < np; i++)
			phase[i] = minterm[vars_[i+shift]];
		BDD current(*mgr_,Cudd_bddComputeCube(mgr,dvars,phase,np));
		D = current;
	}
	int l = ms.size();
	std::cout<<"Input sequences : "<<std::endl;
	for(int i = 0; i < l; i++){
	std::cout<<ms[i]<<" ";
	outfile<<ms[i]<<" ";
	}
	outfile.close();
	std::cout<<std::endl;
}


	std::vector<int> getModes(BDD C,BDD w0){
	std::vector<int> modes;
	BDD temp = C & w0;
	size_t t = (nvars/shift-1)/2;
	size_t np = t * shift;
	const int* minterm;
	int element;
	for(begin(temp,np+shift);!done();next())
	{
		minterm = currentMinterm();
		size_t k = np;
		element = 0;
	     	for(size_t c = 1, j = 0; j < shift; c*=2,k++, j++)
			element+=minterm[vars_[j+np]]*c;
		modes.push_back(element);
	}
	int len = modes.size();
	std::cout<<"Possible mode sequences : "<<std::endl;
	for(int i = 0; i < len; i++)
	std::cout<<modes[i]<<" ";
	std::cout<<std::endl;
	return modes;
}



	BDD getNextAbsState(BDD w0, int p){
	size_t t = (nvars/shift-1)/2;
	size_t np = t * shift;	
	int* mint = new int [np];
	begin(w0,np);
	for(size_t i = 0; i < np; i++)
		mint[i] = 0;
	DdManager *mgr = mgr_->getManager();
	BDD set = mgr_->bddZero();
	DdNode**  dvars = new DdNode*[np];
        for(size_t i = 0; i < np; i++)
   		dvars[i]=Cudd_bddIthVar(mgr,vars_[i]);
	int* phase = new int[shift];
	const int* minterm = currentMinterm();
	for(size_t i = 0; i < shift; i++)
		phase[i]=0;
	int *ph = phase;
	int x = p;
	for (; x; x/=2) *(ph++) = 0 + x%2;
	for(size_t i = 0; i < np-shift; i++)
        	mint[i]=minterm[vars_[i+shift]];
	for(size_t i = 0; i < shift; i++)
		mint[i+np-shift] = phase[i];			
	BDD current(*mgr_,Cudd_bddComputeCube(mgr,dvars,mint,np));
	set +=current;
	return set;
}
	/*template<class F1>
	void printOutput(F1 &system_post,BDD C, BDD w0, size_t ilt, int sDIM, state_type xs){
	BDD w = w0;
	for(size_t i = 0; i < ilt; i++)
	{
		const int* minterm;
		size_t t = (nvars/shift-1)/2;
		size_t np = t * shift;
		int nextmode = 0;
		int* element = new int [t];
		begin(w,np);
		minterm = currentMinterm();
		for(size_t i = 0; i < np+shift; i++)
		std::cout<<minterm[i]<<" ";
		std::cout<<std::endl;
		int* phase = new int [np];
		for(size_t i = 0; i < np; i++)
			phase[i] = 0;
		DdManager *mgr = mgr_->getManager();
		DdNode**  dvars = new DdNode*[np];
		for(size_t i = 0; i < np; i++)
	   		dvars[i]=Cudd_bddIthVar(mgr,vars_[i]);
		size_t k = 0;
		for(size_t i =0; i < t; i++) 
		{
		      	size_t idx=0;
		      	for(size_t c = 1, j = 0; j < shift; c*=2,k++, j++)
				idx+=minterm[vars_[k]]*c;
		      	element[i]=idx;
	    	}
		std::cout<<"Present mode sequence : "<<std::endl;
		for(size_t i = 0; i < t; i++)
		std::cout<<element[i]<<" ";
		std::cout<<std::endl;
		state_type y = getOutput(system_post,element,xs);
		std::cout<<"Present densities : "<<std::endl;
		for(int i = 0 ; i < sDIM; i++)
		std::cout<<y[i]<<" ";
		std::cout<<std::endl;
		std::vector<int> modes;
		modes = getModes(C,w);
		std::cout<<"Enter the mode : ";
		std::cin>>nextmode; 
		w = getNextAbsState(w,nextmode);
	}
}*/


};

#endif
