/*
 * fixedPointMode.hh
 *
 *  created on: 09.10.2015
 *      author: rungger
 */

#ifndef FIXEDPOINTMODE_HH_
#define FIXEDPOINTMODE_HH_

#include <iostream>
#include <stdexcept>

#include "cuddObj.hh"
#include "abstractionMode.hh"
class fixedPointMode {
protected:
  /* var: ddmgr_ */
  Cudd *ddmgr_;
  /* var: ab_ 
   * stores the transition relation */
  abstractionMode *ab_;
  /* var: permute 
   * stores the permutation array used to swap pre with post variables */
  int* permute_;
  size_t nssvars;
  size_t nisvars;
  size_t postvars;
  /* helper BDDs */
  /* transition relation */
  BDD R_;  
  /* transition relation with cubePost_ abstracted */
  BDD RR_;  

   /* cubes with input and post variables; used in the existential abstraction  */
  BDD cubePost_;
  BDD cubeInput_;
  friend class abstractionMode;  
public:

  /* function: fixedPointMode 
   *
   * initialize the fixedPointMode object with a <SymbolicModel> containing the
   * transition relation
   */
  fixedPointMode(abstractionMode *ab) {
     ab_=ab;
    ddmgr_=ab_->mgr_;
    nssvars = (ab_->nvars-ab_->shift)/2;
    nisvars = ab_->shift;
    postvars = nssvars+nisvars;
 /* the permutation array */
    size_t n=ddmgr_->ReadSize();
    permute_ = new int[n];
    for(size_t i=0; i<n; i++)
      permute_[i]=i;
    for(size_t i=0; i<nssvars; i++)
      permute_[ab_->vars_[i]]=ab_->vars_[i+postvars];
    /* create a cube with the input Vars */
    BDD* vars = new BDD[nisvars];
    for (size_t i=0; i<nisvars; i++)
      vars[i]=ddmgr_->bddVar(ab_->vars_[i+nssvars]);
    cubeInput_ = ddmgr_->bddComputeCube(vars,NULL,nisvars);
    delete[] vars;
    /* create a cube with the post Vars */
    vars = new BDD[nssvars];
    for (size_t i=0; i<nssvars; i++)
      vars[i]=ddmgr_->bddVar(ab_->vars_[i+postvars]);   
    cubePost_ = ddmgr_->bddComputeCube(vars,NULL,nssvars);
    delete[] vars;

    /* copy the transition relation */
    R_=ab_->transitionRelation_;
    RR_=R_.ExistAbstract(cubePost_);
  }
  ~fixedPointMode() {
    delete[] permute_;
  }

  /* function: pre 
   *
   * computes the enforcable predecessor 
   *  
   * pre(Z) = { (x,u) | exists x': (x,u,x') in transitionRelation 
   *                    and (x,u,x') in transitionRelation  => x' in Z } 
   *
   */
  BDD pre(BDD Z)  {
    /* project onto state alphabet */
    Z=Z.ExistAbstract(cubePost_*cubeInput_);
    /* swap variables */
    Z=Z.Permute(permute_);
    /* find the (state, inputs) pairs with a post outside the safe set */
    BDD nZ = !Z;
    BDD F = R_.AndAbstract(nZ,cubePost_); 
    /* the remaining (state, input) pairs make up the pre */
    BDD nF = !F;
    BDD preZ= RR_.AndAbstract(nF,cubePost_);
    return preZ;
  }

  /* function: safe 
   *  
   * computation of the maximal fixed point mu Z.pre(Z) & S
   *
   */
  BDD safe(BDD S, int verbose=0)  {
    size_t t = (ab_->nvars/ab_->shift-1)/2;
    size_t np = t * ab_->shift;
    if( S.CountMinterm(np) == 0 ){
    	  std::ostringstream os;
          os << "Error: safe: safe set empty.";
          throw std::invalid_argument(os.str().c_str());
      }
    /* check if safe is a subset of the state space */
    std::vector<unsigned int> sup = S.SupportIndices();
    for(size_t i=0; i<sup.size();i++) {
      int marker=0;
      for(size_t j=0; j<nssvars; j++) {
        if (sup[i]==ab_->vars_[j])
          marker=1;
      }
      if(!marker) {
          std::ostringstream os;
          os << "Error: safe: the inital set depends on variables  outside of the state space.";
          throw std::invalid_argument(os.str().c_str());
      }
    }
    if(verbose) 
      std::cout << "Iterations: ";

    BDD Z = ddmgr_->bddZero();
    BDD ZZ = ddmgr_->bddOne();
    /* as long as not converged */
    size_t i;
    for(i=1; ZZ != Z; i++ ) {
      Z=ZZ;
      ZZ=fixedPointMode::pre(Z) & S;

      /* print progress */
      if(verbose) {
        std::cout << ".";
        std::flush(std::cout);
        if(!(i%80))
          std::cout << std::endl;
      }
    }
    if(verbose) 
      std::cout << " number: " << i << std::endl;
    return Z;
  } 
  
  /* function: reach 
   *  
   * computation of the minimal fixed point mu Z.pre(Z) | T
   *
   */
  BDD reach(const BDD &T, int verbose=0)  {
    /* check if target is a subset of the state space */
    std::vector<unsigned int> sup = T.SupportIndices();
    for(size_t i=0; i<sup.size();i++) {
      int marker=0;
      for(size_t j=0; j<nssvars; j++) {
        if (sup[i]==ab_->vars_[j])
          marker=1;
      }
      if(!marker) {
        std::ostringstream os;
        os << "Error: reach: the target set depends on variables outside of the state space.";
        throw std::invalid_argument(os.str().c_str());
      }
    }
    if(verbose) 
      std::cout << "Iterations: ";

    BDD Z = ddmgr_->bddOne();
    BDD ZZ = ddmgr_->bddZero();
    /* the controller */
    BDD C = ddmgr_->bddZero();
    /* as long as not converged */
    size_t i;
    for(i=1; ZZ != Z; i++ ) {
      Z=ZZ;
      ZZ=fixedPointMode::pre(Z) | T;
      /* new (state/input) pairs */
      BDD N = ZZ & (!(C.ExistAbstract(cubeInput_)));
      /* add new (state/input) pairs to the controller */
      C=C | N;
      /* print progress */
      if(verbose) {
        std::cout << ".";
        std::flush(std::cout);
        if(!(i%80))
          std::cout << std::endl;
      }
    }
    if(verbose) 
      std::cout << " number: " << i << std::endl;
    return C;
  }

  /* function: reachAvoid 
   *  
   * controller synthesis to enforce reach avoid specification
   *
   * computation of the minimal fixed point mu Z.(pre(Z) | T) & !A
   *
   */
  BDD reachAvoid(const BDD &T, const BDD &A, int verbose=0)  {
    /* check if target is a subset of the state space */
    std::vector<unsigned int> sup = T.SupportIndices();
    for(size_t i=0; i<sup.size();i++) {
      int marker=0;
      for(size_t j=0; j<nssvars; j++) {
        if (sup[i]==ab_->vars_[j])
          marker=1;
      }
      if(!marker) {
          std::ostringstream os;
          os << "Error: reachAvoid: the inital set depends on variables  outside of the state space.";
          throw std::invalid_argument(os.str().c_str());
      }
    }
    if(verbose) 
      std::cout << "Iterations: ";

    BDD RR=RR_;
    /* remove avoid (state/input) pairs from transition relation */
    RR_= RR & (!A);
    BDD TT= T & (!A);

    BDD Z = ddmgr_->bddOne();
    BDD ZZ = ddmgr_->bddZero();
    /* the controller */
    BDD C = ddmgr_->bddZero();
    /* as long as not converged */
    size_t i;
    for(i=1; ZZ != Z; i++ ) {
      Z=ZZ;
      ZZ=fixedPointMode::pre(Z) | TT;
      /* new (state/input) pairs */
      BDD N = ZZ & (!(C.ExistAbstract(cubeInput_)));
      /* add new (state/input) pairs to the controller */
      C=C | N;
      /* print progress */
      if(verbose) {
        std::cout << ".";
        std::flush(std::cout);
        if(!(i%80))
          std::cout << std::endl;
      }
    }
    if(verbose) 
      std::cout << " number: " << i << std::endl;
    /* restor transition relation */
    RR_=RR;
    return C;
  }
    
    /* func: reachStay */
    BDD reachStay(const BDD &T)  {
        /* check if target is a subset of the state space */
        std::vector<unsigned int> sup = T.SupportIndices();
        for(size_t i=0; i<sup.size();i++) {
            int marker=0;
            for(size_t j=0; j<nssvars; j++) {
                if (sup[i]==ab_->vars_[j])
                    marker=1;
            }
            if(!marker) {
                std::ostringstream os;
                os << "Error: reachAvoid: the inital set depends on variables  outside of the state space.";
                throw std::invalid_argument(os.str().c_str());
            }
        }
        /* we implement the nested fixed point algorithm
         *
         * mu X. nu Y. ( pre(Y) & T ) | pre(X)
         *
         */
        size_t i,j;
        /* outer fp*/
        BDD X=ddmgr_->bddOne();
        BDD XX=ddmgr_->bddZero();
        /* inner fp*/
        BDD Y=ddmgr_->bddZero();
        BDD YY=ddmgr_->bddOne();
        /* the controller */
        BDD C=ddmgr_->bddZero();
        
        /* as long as not converged */
        for(i=1; XX != X; i++) {
            X=XX;
            BDD preX=pre(X);
            /* init inner fp */
            YY = ddmgr_->bddOne();
            for(j=1; YY != Y; j++) {
                Y=YY;
                YY= ( pre(Y) & T ) | preX;
            }
            XX=YY;
            std::cout << "Iterations inner: " << j << std::endl;
            /* remove all (state/input) pairs that have been added
             * to the controller already in the previous iteration * */
            BDD N = XX & (!(C.ExistAbstract(cubeInput_)));
            /* add the remaining pairs to the controller */
            C=C | N;
            //std::cout << C.CountMinterm(17) << std::endl;
        }
        std::cout << "Iterations outer: " << i << std::endl;
        return C;
    }
    
    
    
}; /* close class def */

#endif /* FIXEDPOINTMODE_HH_ */
