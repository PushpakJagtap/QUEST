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
	size_t dim_;
	size_t* xx_;
	double* lb_;
	double* ub_;
	double* eta_;
	double** ud_;
	size_t Z_;
	bool flag_;

	/*getAbstraction(size_t Dim_):dim_(Dim_){ }*/
/******************************************************************
                      For discrete input set
*******************************************************************/
/*************************getAbstractSet********************************/

	template<class F1,class F2>
	BDD getAbstractSet(F1 &system_post,F2 &setBounds,state_type x, const size_t dim, const size_t Z, double **ud){
		flag_=0;
		BDD set = mgr_->bddZero();
		const int* minterm;
		size_t t = (nvars/shift-1)/2;
		size_t np = t * shift;
		dim_=dim;
		Z_=Z;
		ud_ =new double* [Z_];
		for(size_t i=0;i<Z_;i++){
		    ud_[i]=new double[dim];
		}

		for(size_t b=0; b<Z;b++){
			for(size_t a=0; a<dim;a++){
				ud_[b][a]=ud[b][a];
			}
		}
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
        	if(setBounds(y))
        	{
			for(size_t i = 0; i < np; i++)
         			phase[i]=minterm[vars_[i]];			
			BDD current(*mgr_,Cudd_bddComputeCube(mgr,dvars,phase,np));
			set +=current;
        	}
        }
        std::cout<<std::endl;
        std::cout<<" Number of elements in the safe/target set : ";
        std::cout<<set.CountMinterm(np)<<std::endl;
        std::cout<<std::endl;
        return set;
	}

/******************************************************************
	                    For continuous input set
*******************************************************************/

/****************************getAbstractSet********************************/
	template<class F1,class F2>
	BDD getAbstractSet(F1 &system_post,F2 &setBounds,state_type x, const size_t dim, const double* lb, const double* ub, const double* eta){
			lb_ =new double[dim];
			ub_ =new double[dim];
			eta_ =new double[dim];
			xx_ =new size_t[dim-1];
			dim_=dim;
			double Nl, Nu;
		    size_t nInput[dim];
			size_t temp1=1;;
			for(size_t i=0;i<dim;i++){
				lb_[i]=lb[i];
				ub_[i]=ub[i];
				Nl=std::ceil(lb[i]/eta[i]);
				Nu=std::floor(ub[i]/eta[i]);
				nInput[i]=(size_t)std::abs(Nu-Nl)+1;
			}
			for(int j=dim-1;j>=0;j--){
				temp1=temp1*nInput[j];
				xx_[j-1]=temp1;
			}
			for(size_t k=0;k<dim;k++){
				eta_[k]=eta[k];
			}
			flag_=1;
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
				if(setBounds(y))
				{
					for(size_t i = 0; i < np; i++)
						phase[i]=minterm[vars_[i]];
					BDD current(*mgr_,Cudd_bddComputeCube(mgr,dvars,phase,np));
					set +=current;
				}
			}
			std::cout<<std::endl;
			std::cout<<" Number of elements in the safe/target set : ";
			std::cout<<set.CountMinterm(np)<<std::endl;
			std::cout<<std::endl;
			return set;
	}

/***********************************getOutput**************************************/

		template<class F1>
		state_type getOutput(F1 &system_post, int* element, state_type x){
			double u[dim_];
			state_type y = x;
			size_t t = (nvars/shift-1)/2;
			for(size_t j = 0; j < t; j++)
			{
				if (flag_==1){
					size_t temp2,temp3,temp4;
					for(size_t i=0; i<dim_; i++){

						if(i<dim_-1){
							temp2=element[j]/xx_[i];
							temp3=temp2%xx_[i];
							u[i]=lb_[i]+temp3*eta_[i];
						}
						else{
							temp4=element[j]%xx_[i-1];
							u[i]=lb_[i]+temp4*eta_[i];
						}
					}
				}
				else{
					size_t temp5=element[j];
					for(size_t i=0; i<dim_; i++){
						u[i]=ud_[temp5][i];
					}
				}
				system_post(y,u);
			}
			return y;
		}

/***********************************get initial state(reachability)**************************************/

	template<class F1>
	BDD getAbsState(F1 &system_post, state_type x,state_type xs,double epsilon, int sDIM){
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
            if(std::abs(y[i]-x[i]) > epsilon)
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

/***********************************get initial state(safety)**************************************/

   template<class F1>
    BDD getAbsState(F1 &system_post, BDD set, state_type x,state_type xs,double epsilon, int sDIM){
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
                if(abs(y[i]-x[i]) > epsilon)
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

/***********************************closedLoopSim**************************************/

    template<class F1>
    void closedLoopSim(BDD C, BDD D, F1 &system_post, state_type y, int sDIM, int t){
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
		for(size_t k=0;k<dim_;k++)
			std::cout<<"u("<<k<<")"<<'\t';
		std::cout<<std::endl;

		for(size_t j = 0; j < t; j++){
			mapNprint(element[j]);
			std::cout<<std::endl;
		}

		for(size_t c = 1, j = 0; j < shift; c*=2,k++, j++)
			nextmode+=minterm[vars_[k]]*c;
		ms.push_back(nextmode);
		std::cout<<"Next input = ";
		mapNprint(nextmode);
		std::cout<<std::endl;
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
		mapNprint(ms[i]);
		std::cout<<std::endl;
		outfile<<ms[i]<<" ";
	}
	outfile.close();
	std::cout<<std::endl;
}

    void mapNprint(int element){
    	double u[dim_];
    	size_t temp2,temp3,temp4;
    		for(size_t k=0; k<dim_; k++){
    			if(flag_==1){
    				if(k<dim_-1){
    					temp2=element/xx_[k];
    					temp3=temp2%xx_[k];
    					u[k]=lb_[k]+temp3*eta_[k];
    				}
    				else{
    					temp4=element%xx_[k-1];
   						u[k]=lb_[k]+temp4*eta_[k];
   					}

   				}
    			else{
    				u[k]=ud_[element][k];
    			}
    			std::cout<<u[k]<<'\t';
    		}
    }

};

#endif
