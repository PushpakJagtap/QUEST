#ifndef SYMBOLICSETMODE_HH_
#define SYMBOLICSETMODE_HH_



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

#include "SymbolicSetMode.hh"


class SymbolicSetMode{

private :
/* friend class if any */
friend class abstraction;
/*variables */

Cudd* mgr_;
size_t P_;
size_t N_;
size_t *nofBddVars_;
size_t **indBddVars_;
size_t *nofAbsStates_;
size_t nvars_;
BDD symbolicSetMode_;
CuddMintermIterator *iterator_;


public :

/* constructor */ 

	SymbolicSetMode( Cudd &mgr,const size_t P,const size_t N){
		if ( P == 0 )
		{
			std::ostringstream os;
			os << "Error: SymbolicSetMode:  each interval must contain at least one mode.";
			throw std::invalid_argument(os.str().c_str());
		}
		mgr_ = &mgr;
		P_ = P;
		N_ = N;
		nvars_ = 0;
		nofBddVars_ = new size_t [N];
		indBddVars_ = new size_t* [N];
		nofAbsStates_ = new size_t [N];
		symbolicSetMode_ = mgr.bddZero();
    		iterator_=NULL;
		for(size_t j = 0; j < N; j++)
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
		for(size_t i = 0; i < N; i++)
		{
			for(size_t j = 0; j < nofBddVars_[i]; j++)
			{
				nvars_++;
			}
		} 
		
}

	SymbolicSetMode(const SymbolicSetMode& ss,const SymbolicSetMode& im){
		if(ss.mgr_==im.mgr_)
      			mgr_=ss.mgr_;
    		else 
		{
      			std::ostringstream os;
      			os << "Error: SymbolicSetMode ss and im do not have the same dd manager.";
      			throw std::runtime_error(os.str().c_str());
    		}
		if(ss.P_==im.P_)
      			P_=ss.P_;
    		else 
		{
      			std::ostringstream os;
      			os << "Error: SymbolicSetMode ss and im do not have the same modes.";
      			throw std::runtime_error(os.str().c_str());
    		}
		N_ = ss.N_ + im.N_;
		nvars_ = 0;
		nofBddVars_ = new size_t [N_];
		indBddVars_ = new size_t* [N_];
		nofAbsStates_ = new size_t [N_];
		symbolicSetMode_ = ss.symbolicSetMode_ * im.symbolicSetMode_;
    		iterator_=NULL; 
		for (size_t i = 0; i < ss.N_; i++)  
		{
      
      			nofAbsStates_[i] = ss.nofAbsStates_[i];
      			nofBddVars_[i] = ss.nofBddVars_[i];
    		}
    		for (size_t i = 0; i < im.N_; i++)
	  	{
      			nofAbsStates_[ss.N_+i] = im.nofAbsStates_[i];
      			nofBddVars_[ss.N_+i] = im.nofBddVars_[i];
    		}
    		for (size_t i = 0; i < ss.N_; i++)
	 	{
      			indBddVars_[i]=new size_t[ss.nofBddVars_[i]];
      			for(size_t j = 0; j < nofBddVars_[i]; j++) 
        			indBddVars_[i][j]=ss.indBddVars_[i][j];
    		}
    		for (size_t i = 0; i < im.N_; i++)
		{
      			indBddVars_[ss.N_+i]=new size_t [im.nofBddVars_[i]];
      			for(size_t j = 0; j < nofBddVars_[ss.N_+i]; j++) 
        		indBddVars_[ss.N_+i][j]=im.indBddVars_[i][j];
    		}
    		for(size_t i = 0; i < N_; i++)
		{
      			for(size_t j = 0; j < nofBddVars_[i]; j++)
        			nvars_++;
  		}
}
	~SymbolicSetMode(){
	        delete[] nofAbsStates_;
		delete[] nofBddVars_;
		for(size_t i=0; i<N_; i++)
			delete[] indBddVars_[i];
		delete[] indBddVars_;
}

	BDD computeAbsStates(void){
		BDD AbsStates = mgr_->bddOne();
		BDD** bddVars = new BDD* [N_];
		BDD* symset = new BDD [N_];
		for(size_t i=0; i<N_; i++) 
		{
		      symset[i]=mgr_->bddZero();
		      bddVars[i]= new BDD[nofBddVars_[i]];
		      for(size_t j=0; j<nofBddVars_[i]; j++) 
			bddVars[i][j]=mgr_->bddVar(indBddVars_[i][j]); 
		}
		for(size_t i=0; i<N_; i++)
		{
		      int *phase = new int[nofBddVars_[i]];
		      for(size_t j=0;j<nofBddVars_[i];j++)
				phase[j]=0;
		      for(size_t j=0;j<nofAbsStates_[i];j++)
	 		{
				int *p=phase;
				int x=j;
				for (; x; x/=2) *(p++)=0+x%2;
				symset[i]+= mgr_->bddComputeCube(bddVars[i],(int*)phase,nofBddVars_[i]);
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

	void addAbsStates(void) {
	        symbolicSetMode_|=SymbolicSetMode::computeAbsStates();
}

	inline void clear() {
    		symbolicSetMode_ = mgr_->bddZero();
}

	inline void mintermToElement(const int* minterm, size_t* element) const {
		for(size_t i=0; i<N_; i++)
		{
			size_t idx=0;
		      	for(size_t c=1, j=0; j<nofBddVars_[i]; c*=2, j++)
				idx+=minterm[indBddVars_[i][j]]*c;
		      	element[i]=idx;
		}
}

	inline void mintermToElement(const int* minterm, size_t element) const {
		if( N_ != 1 )
		{
			std::ostringstream os;
			os<<"ERROR: this function cannot be used try using other function. ";
			throw std::invalid_argument(os.str().c_str());
		}
		for(size_t i=0; i<N_; i++)
		{
			size_t idx=0;
		      	for(size_t c=1, j=0; j<nofBddVars_[i]; c*=2, j++)
				idx+=minterm[indBddVars_[i][j]]*c;
		      	element=idx;
		}
}

	SymbolicSetMode& operator=(const SymbolicSetMode& ss){
		mgr_ = ss.mgr_;
		N_= ss.N_;
		P_ = ss.P_;
		nvars_ = ss.nvars_;
		nofBddVars_ = new size_t [N_];
		indBddVars_ = new size_t* [N_];
		nofAbsStates_ = new size_t [N_];
		symbolicSetMode_ = mgr_->bddZero();
    		iterator_=NULL;
		size_t nofbits = std::ceil(log2(P_));
		for(size_t j = 0; j < N_; j++)
		{
			nofBddVars_[j] = nofbits;
			nofAbsStates_[j] = P_;
		}
		for(size_t i = 0; i < N_; i++)
		{
			indBddVars_[i] = new size_t [ss.nofBddVars_[i]];
			for(size_t j = 0; j < ss.nofBddVars_[i]; j++)
			{
				indBddVars_[i][j] = ss.indBddVars_[i][j];
			}
		}
		return *this;
}

	SymbolicSetMode(const SymbolicSetMode& ss, size_t flag){
		*this = ss;
		 this->clear();
		if(flag)
		{
		 	for(size_t i = 0; i < N_; i++)
			{
				for(size_t j = 0; j < nofBddVars_[i]; j++)
				{
					BDD var = mgr_->bddVar();
					indBddVars_[i][j] = var.NodeReadIndex();
				}
			}
		}
}


	inline BDD getSymbolicSetMode(void) const {
		return symbolicSetMode_;
}

	inline size_t* getnofBddVars(void) const {
    		return nofBddVars_;
}
  
	inline size_t** getindBddVars(void) const {
    		return indBddVars_;
}
  
  
	inline size_t getSize(void) {
    		return symbolicSetMode_.CountMinterm(nvars_);
}

	inline size_t* getnofAbsStates(void) const {
    		return nofAbsStates_;
}
 
	inline size_t getN(void) const {
    		return N_;
}					

	inline size_t getP(void) const {
    		return P_;
}


	inline void begin(void) {
    		std::vector<size_t> ivars_;
    		ivars_.reserve(nvars_);
    		for(size_t i = 0; i < N_; i++) 
      			for(size_t j = 0; j < nofBddVars_[i]; j++) 
        			ivars_.push_back(indBddVars_[i][j]);
    		iterator_ = new CuddMintermIterator(symbolicSetMode_,ivars_,nvars_);
}
   	inline void next(void) {
    		++(*iterator_);
}
  /*	inline void progress(void) {
    		iterator_->printProgress();
} */
	inline int done(void) {
    		if(iterator_->done())
		{
      			delete iterator_;
      			iterator_=NULL;
      			return 1;
    		}
		else 
    		return 0;
}
	inline const int* currentMinterm(void) {
    		if (iterator_)
      			return iterator_->currentMinterm();
    		else
      			return NULL;
}

};			

#endif /* SYMBOLICSETMODE_HH_ */

		
		




























