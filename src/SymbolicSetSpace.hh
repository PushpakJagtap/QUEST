#ifndef SYMBOLICSETSPACE_HH_
#define SYMBOLICSETSPACE_HH_

#include <vector>
#include <algorithm>
#include <sstream>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>


#include "dddmp.h"
#include "cuddObj.hh"
#include "CuddMintermIterator.hh"
#include "SymbolicSetSpace.hh"


class SymbolicSetSpace{

private :

/* friend class if any */
friend class abstractionMode;
/*variables */

Cudd* mgr_;
size_t P_;
size_t N_;
size_t *nofBddVars_;
size_t **indBddVars_;
size_t *nofAbsStates_;
size_t nvars_;
BDD symbolicSetSpace_;



public :

/* constructor */ 

	SymbolicSetSpace( Cudd &mgr,const size_t P,const size_t N){
		if ( P == 0 )
		{
			std::ostringstream os;
			os << "Error: SymbolicSetSpace:  each interval must contain at least one mode.";
			throw std::invalid_argument(os.str().c_str());
		}
		mgr_ = &mgr;
		P_ = P;
		N_ = 2*N + 1;
		nvars_ = 0;
		nofBddVars_ = new size_t [N_];
		indBddVars_ = new size_t* [N_];
		nofAbsStates_ = new size_t [N_];
		symbolicSetSpace_ = mgr.bddZero();
		for(size_t j = 0; j < N_; j++)
		{
			nofAbsStates_[j] = P;
			if( P == 1 )
				nofBddVars_[j] = 1;
			else
				nofBddVars_[j] = std::ceil(log2(P));
			indBddVars_[j] = new size_t [nofBddVars_[j]];
			for(size_t i = 0; i < nofBddVars_[j]; i++)
			{
				BDD var = mgr.bddVar();
				indBddVars_[j][i] = var.NodeReadIndex();
			}
		}
		for(size_t i = 0; i < N_; i++)
		{
			for(size_t j = 0; j < nofBddVars_[i]; j++)
			{
				nvars_++;
			}
		} 
		
}
	
	~SymbolicSetSpace(){
	        delete[] nofAbsStates_;
		delete[] nofBddVars_;
		for(size_t i=0; i<N_; i++)
			delete[] indBddVars_[i];
		delete[] indBddVars_;
}


	BDD computeAbsStates(void){
		BDD AbsStates = mgr_->bddOne();
		size_t N = (N_+1)/2;
		BDD** bddVars = new BDD* [N_];
		BDD* symset = new BDD [N_];
		for(size_t i=0; i<N_; i++) 
		{
		      symset[i]=mgr_->bddZero();
		      bddVars[i]= new BDD[nofBddVars_[i]];
		      for(size_t j=0; j<nofBddVars_[i]; j++) 
			bddVars[i][j]=mgr_->bddVar(indBddVars_[i][j]); 
		}
		for(size_t i=0; i<N; i++)
		{
		      int *phase = new int[nofBddVars_[i]];
		      for(size_t j=0;j<nofBddVars_[i];j++)
				phase[j]=0;
		      for(size_t j=0;j<nofAbsStates_[i];j++)
	 		{
				int *p=phase;
				int x=j;
				for (; x; x/=2) *(p++)=0+x%2;
				if( i ==  0 )
				symset[i]+= mgr_->bddComputeCube(bddVars[i],(int*)phase,nofBddVars_[i]);
				else if ( i > 0 )
				{
					BDD current = mgr_->bddComputeCube(bddVars[i],(int*)phase,nofBddVars_[i]);
					BDD post = mgr_->bddComputeCube(bddVars[i+N-1],(int*)phase,nofBddVars_[i+N-1]);
					symset[i]+=current & post;
				}
		      	}
		      delete[] phase;
		      AbsStates&=symset[i];
		}
		for(size_t i=0; i<N_; i++) 
		      delete[] bddVars[i];
		delete[] bddVars;
		delete[] symset;
		return AbsStates;
}

	inline size_t getSize(void) {
    		return symbolicSetSpace_.CountMinterm(nvars_);
}

	void addAbsStates(void) {
	        symbolicSetSpace_|=SymbolicSetSpace::computeAbsStates();
}
	inline void clear() {
    		symbolicSetSpace_ = mgr_->bddZero();
}

	inline void printinfo(void){
		symbolicSetSpace_.print(2,2);
}

	


};

#endif /* SYMBOLICSETSPACE_HH_ */
