/*
    This file is part of dynsbm.

    dysbm is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    dynsbm is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with dynsbm.  If not, see <http://www.gnu.org/licenses/>
 */
#ifndef DYNSBM_DYNSBMDISCRETE_H
#define DYNSBM_DYNSBMDISCRETE_H
#include<DynSBM.h>
namespace dynsbm{
  class DynSBMDiscrete
    : public DynSBM<int>{
  protected:
    int _K; // Y has integer values in [1,2,...K]
    double**** _multinomprobaql;
    void correctMultinomproba();
    void addEvent(double proba, int y, int t, int q, int l){
      _multinomprobaql[t][q][l][y-1] += proba; // y is an integer value in [1,2,...K]
    }
  public:
    DynSBMDiscrete(int T, int N, int Q, const Rcpp::IntegerMatrix & present, bool isdirected = false, bool withselfloop = false)
      : DynSBM<int>(T,N,Q,present,isdirected,withselfloop) {
    }
    ~DynSBMDiscrete(){
      deallocate4D(_multinomprobaql,_T,_Q,_Q,_K);
    }
    void setK(int K){ // mandatory
      _K=K;
      allocate4D(_multinomprobaql,_T,_Q,_Q,_K);
    }
    double**** const getMultinomproba() const{
      return(_multinomprobaql);
    }
    virtual double logDensity(int t, int q, int l, int y) const{
      if(y==0){
	return(log(_betaql[t][q][l]));
      } else{
	return(log(1-_betaql[t][q][l]) + log(_multinomprobaql[t][q][l][y-1])); // y is an integer value in [1,2,...K]
      }
    }
    virtual void updateTheta(int*** const Y);
    friend class DynSBMDiscreteAddEventFunctor;
  };

  class DynSBMDiscreteAddEventFunctor{
    DynSBMDiscrete& _dynsbm;
  public:
    DynSBMDiscreteAddEventFunctor(DynSBMDiscrete& dynsbm)
      : _dynsbm(dynsbm) {}
    void operator()(double proba, int y, int t, int q, int l){
      _dynsbm._multinomprobaql[t][q][l][y-1] += proba; // y is an integer value in [1,2,...K]
    }
  };
}
#endif
