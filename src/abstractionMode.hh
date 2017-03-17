#ifndef ABSTRACTION_HH_
#define ABSTRACTION_HH_

#include <iostream>
#include <stdexcept>

#include "cuddObj.hh"
#include "dddmp.h"
#include "TicToc.hh"
#include "cudd.h"
#include "abstractionMode.hh"
#include "SymbolicSetSpace.hh"

class abstractionMode{

protected :
Cudd* mgr_;
SymbolicSetSpace* ss_;
size_t *vars_;
size_t nvars;
size_t shift;
CuddMintermIterator *iterator_;
BDD transitionRelation_;
BDD it;
friend class SymbolicSetSpace; 
friend class fixedPointMode;
public :
		
	abstractionMode( SymbolicSetSpace *ss){
	ss_ = ss;
	mgr_ = ss_->mgr_;
	shift = ss_->nofBddVars_[0];
	nvars = 0;
	for(size_t i = 0; i < ss_->N_; i++)
		for(size_t j = 0; j < shift; j++)
			nvars++;
	if(nvars != ss_->nvars_)
	{
		std::ostringstream os;
		os<<"ERROR: nssvars are not equal. ";
		throw std::invalid_argument(os.str().c_str());
	}
	vars_ = new size_t [nvars];
	size_t k= 0;
	for(size_t i = 0; i < ss_->N_; i++)
		for(size_t j = 0; j < shift; k++,j++)
		vars_[k] = ss_->indBddVars_[i][j];
	transitionRelation_ = ss_->symbolicSetSpace_;
	iterator_ = NULL;
}



inline void begin(BDD transitionRelation_,size_t n) {
        std::vector<size_t> ivars_;
        ivars_.reserve(n);
        for(size_t i=0; i<n; i++) 
        	ivars_.push_back(vars_[i]);
        iterator_ = new CuddMintermIterator(transitionRelation_,ivars_,n);
}

inline void next(void) {
        ++(*iterator_);
}

inline void progress(void) const {
        iterator_->printProgress();
}  

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
  
inline const int* currentMinterm(void) const {
        if (iterator_)
        	return iterator_->currentMinterm();
        else
                return NULL;
}
	
inline void printinfo(void){
		transitionRelation_.print(2,2);
}

	inline double getSize(void) {
    		return transitionRelation_.CountMinterm(nvars);
}

	

};

#endif
