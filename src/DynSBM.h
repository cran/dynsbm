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
#ifndef DYNSBM_DYNSBM_H
#define DYNSBM_DYNSBM_H
#include <Rcpp.h>
#include<vector>
#include<math.h>
#include<tools.h>
#include<limits>
#include<iostream>
#include <Rcpp.h>
namespace dynsbm{
  template<typename Ytype> 
  class DynSBM{
  protected:
    int _T; // time steps
    int _N; // number of nodes
    int _Q; // number of groups
    bool _isdirected;
    bool _withselfloop;
    const Rcpp::IntegerMatrix & _present;
    bool ispresent(int t, int i) const{
      return(_present(i,t)!=0);
    }
    bool isabsent(int t, int i) const{
      return(_present(i,t)==0);
    }
    // Markov chain
    double* _stationary;
    double** _trans;
    // variational EM
    double** _tau1;
    double**** _taut;
    double*** _taum;
    void correctTau1();
    void correctTaut();
    void correctTaum();
    void updateTauMarginal();
    double tauMarginal(int t, int i, int q) const; // proba of q iif present(t,i)==true
    double tauArrival(int t, int i, int q) const{ // proba of q when absent before t>0
      int tstorage = t-1;
      return(_taut[tstorage][i][0][q]);
    }
    // edge distribution
    double*** _betaql;
    void correctBeta();
    virtual double logDensity(int t, int q, int l, Ytype y) const = 0;
    template<class TAddEventFunctor> 
    void updateThetaCore(Ytype*** const Y, TAddEventFunctor& addEventFunctor);
  public:
    DynSBM(int T, int N, int Q, const Rcpp::IntegerMatrix & present, bool isdirected = false, bool withselfloop = false);
    ~DynSBM();
    // Markov chain
    double** const getTrans() const{
      return(_trans);
    }
    double*** const getBeta() const{
      return(_betaql);
    }
    void updateTrans();
    void updateStationary();
    void initNotinformativeStationary();
    void initNotinformativeTrans();
    // variational EM
    void initTau(const std::vector<int>& clustering);
    void updateTau(Ytype*** const Y);
    virtual void updateTheta(Ytype*** const Y) = 0;
    double completedLoglikelihood(Ytype*** const Y) const;
    // model output
    std::vector<int> getGroupsByMAP(int t) const;
    double modelselectionLoglikelihood(Ytype*** const Y) const;
  };
  
  template<typename Ytype> 
  DynSBM<Ytype>::DynSBM(int T, int N, int Q, const Rcpp::IntegerMatrix & present, bool isdirected, bool withselfloop)
    : _T(T),_N(N),_Q(Q),_present(present),_isdirected(isdirected), _withselfloop(withselfloop) {
    _stationary = new double[_Q];
    allocate2D<double>(_trans,_Q,_Q);
    allocate2D<double>(_tau1,_N,_Q);
    allocate4D<double>(_taut,_T-1,_N,_Q,_Q);
    allocate3D<double>(_taum,_T,_N,_Q);
    allocate3D<double>(_betaql,_T,_Q,_Q);
  }
  
  template<typename Ytype> 
  DynSBM<Ytype>::~DynSBM(){
    delete[] _stationary;
    deallocate2D<double>(_trans,_Q,_Q);
    deallocate2D<double>(_tau1,_N,_Q);
    deallocate4D<double>(_taut,_T-1,_N,_Q,_Q);
    deallocate3D<double>(_taum,_T,_N,_Q);
    deallocate3D<double>(_betaql,_T,_Q,_Q);
  }
  
  template<typename Ytype> 
  void DynSBM<Ytype>::updateTrans(){
    if (_Q==1){
      _trans[0][0] = 1.;
      return;
    }
    double** denom = new double*[_Q];
    for(int q=0;q<_Q;q++) denom[q] = new double[_Q];
    for(int q=0;q<_Q;q++) 
      for(int l=0;l<_Q;l++){
	_trans[q][l] = 0.; denom[q][l] = 0.;
      }
    for(int t=1;t<_T;t++){
      int tstorage = t-1;
      for(int i=0;i<_N;i++){
	if (ispresent(t-1,i) && ispresent(t,i)){
	  for(int q=0;q<_Q;q++){
	    for(int qprime=0;qprime<_Q;qprime++){
	      _trans[q][qprime] += _taut[tstorage][i][q][qprime] * tauMarginal(t-1,i,q);
	      denom[q][qprime] += tauMarginal(t-1,i,q);
	    }
	  }
	}
      }
    }
    // numerical issue : avoid too small value for trans
    for(int q=0;q<_Q;q++){
      double sumtrans=0.;
      for(int qprime=0;qprime<_Q;qprime++){
	_trans[q][qprime] = _trans[q][qprime] / denom[q][qprime];
	sumtrans += _trans[q][qprime];
      }
      if (sumtrans>0) for(int qprime=0;qprime<_Q;qprime++) _trans[q][qprime] = _trans[q][qprime]/sumtrans;
      for(int qprime=0;qprime<_Q;qprime++){
	if(_trans[q][qprime]<1e-7){
	  _trans[q][qprime] = 1e-7;
	}
      }
      sumtrans=0;
      for(int qprime=0;qprime<_Q;qprime++) sumtrans += _trans[q][qprime];
      if (sumtrans>0) for(int qprime=0;qprime<_Q;qprime++) _trans[q][qprime] = _trans[q][qprime]/sumtrans;
    }
    for(int q=0;q<_Q;q++) delete[] denom[q];
    delete[] denom;
  }
  
  template<typename Ytype> 
  void DynSBM<Ytype>::updateStationary(){
    if (_Q==1){
      _stationary[0] = 1.;
      return;
    }
    for(int q=0;q<_Q;q++) _stationary[q]=0.;
    double sumst=0;
    for(int q=0;q<_Q;q++){
      for(int t=0;t<_T;t++){
	for(int i=0;i<_N;i++){
	  if (ispresent(t,i))
	    _stationary[q] += tauMarginal(t,i,q);
	}
      }
      if(_stationary[q]<1e-7){
	_stationary[q] = 1e-7;
      }
      sumst += _stationary[q];
    }
    for(int q=0;q<_Q;q++) _stationary[q] = _stationary[q]/sumst;
  }
  
  template<typename Ytype> 
  void DynSBM<Ytype>::initNotinformativeStationary(){
    for(int q=0;q<_Q;q++) _stationary[q] = 1./_Q;
  }
  
  template<typename Ytype> 
  void DynSBM<Ytype>::initNotinformativeTrans(){
    for(int q=0;q<_Q;q++) for(int l=0;l<_Q;l++) _trans[q][l] = 1./_Q;
  }
  
  template<typename Ytype>
  void DynSBM<Ytype>::updateTau(Ytype*** const Y){
    if(_Q==1) return;
    double** newtau1;
    allocate2D<double>(newtau1,_N,_Q);
    for(int i=0;i<_N;i++) 
      for(int q=0;q<_Q;q++) 
	newtau1[i][q] = 0.;
    double**** newtaut;
    allocate4D<double>(newtaut,_T-1,_N,_Q,_Q);
    for(int t=0;t<_T-1;t++) 
      for(int i=0;i<_N;i++) 
	for(int q=0;q<_Q;q++) 
	  for(int l=0;l<_Q;l++) 
	    newtaut[t][i][q][l] = 0;  
    // t=1
#pragma omp parallel for
    for(int i=0;i<_N;i++){
      if (ispresent(0,i)){
	  double maxlogtau1i = -std::numeric_limits<double>::max();
	  std::vector<double> logtau1i(_Q,0.);
	  for(int q=0;q<_Q;q++){
	    double logp = 0.;
	    for(int j=0;j<i;j++)
	      if (ispresent(0,j))
		for(int l=0;l<_Q;l++){
		  logp += _tau1[j][l]*logDensity(0,q,l,Y[0][i][j]);
		  if (_isdirected) logp += _tau1[j][l]*logDensity(0,l,q,Y[0][j][i]);
		}
	    if(_withselfloop) logp += logDensity(0,q,q,Y[0][i][i]);
	    for(int j=i+1;j<_N;j++)
	      if (ispresent(0,j))
		for(int l=0;l<_Q;l++){
		  logp += _tau1[j][l]*logDensity(0,q,l,Y[0][i][j]);
		  if (_isdirected) logp += _tau1[j][l]*logDensity(0,l,q,Y[0][j][i]);
		}
	    logtau1i[q] = logp + log(_stationary[q]);
	    if (logtau1i[q]>maxlogtau1i) maxlogtau1i = logtau1i[q];
	  }
	  // numerical issue : normalization
	  std::vector<double> tau1i(_Q,0);
	  double sumtau1i = 0.;
	  for(int q=0;q<_Q;q++){
	    tau1i[q] = exp(logtau1i[q]-maxlogtau1i);
	    sumtau1i = sumtau1i+tau1i[q];
	  }
	  for(int q=0;q<_Q;q++) newtau1[i][q] = tau1i[q]/sumtau1i;
      }
    }
    // t>1
    for(int t=1;t<_T;t++){
      int tstorage = t-1;
#pragma omp parallel for
      for(int i=0;i<_N;i++){
	if (ispresent(t,i)){
	  std::vector<double> logps(_Q,0.);
	  for(int qprime=0;qprime<_Q;qprime++){
	    double logp = 0.;
	    for(int j=0;j<i;j++)
	      if (ispresent(t,j))
		for(int l=0;l<_Q;l++){
		  logp += tauMarginal(t,j,l)*logDensity(t,qprime,l,Y[t][i][j]);
		  if (_isdirected) logp += tauMarginal(t,j,l)*logDensity(t,l,qprime,Y[t][j][i]);
		}
	    if(_withselfloop) logp += logDensity(t,qprime,qprime,Y[t][i][i]);
	    for(int j=i+1;j<_N;j++)
	      if (ispresent(t,j))
		for(int l=0;l<_Q;l++){
		  logp += tauMarginal(t,j,l)*logDensity(t,qprime,l,Y[t][i][j]);
		  if (_isdirected) logp += tauMarginal(t,j,l)*logDensity(t,l,qprime,Y[t][j][i]);
		}
	    logps[qprime] = logp;
	  }
	  for(int q=0;q<(ispresent(t-1,i)?_Q:1);q++){ // only q=0 if absent at t-1 (see tauArrival method)
	    double maxlogtauti = -std::numeric_limits<double>::max();
	    std::vector<double> logtauti(_Q,0.);
	    if (ispresent(t-1,i))
	      for(int qprime=0;qprime<_Q;qprime++){
		logtauti[qprime] = logps[qprime] + log(_trans[q][qprime]);
		if (logtauti[qprime]>maxlogtauti) maxlogtauti = logtauti[qprime];
	      }
	    else
	      for(int qprime=0;qprime<_Q;qprime++){
		logtauti[qprime] = logps[qprime] + log(_stationary[qprime]);
		if (logtauti[qprime]>maxlogtauti) maxlogtauti = logtauti[qprime];
	      }
	    // numerical issue : normalization
	    std::vector<double> tauti(_Q,0);
	    double sumtauti = 0.;
	    for(int qprime=0;qprime<_Q;qprime++){
	      tauti[qprime] = exp(logtauti[qprime]-maxlogtauti);
	      sumtauti = sumtauti+tauti[qprime];
	    }
	    for(int qprime=0;qprime<_Q;qprime++){
	      tauti[qprime] = tauti[qprime]/sumtauti;
	      newtaut[tstorage][i][q][qprime] = tauti[qprime];
	    }
	  }
	}
      }
    }
    // switching tau for old value at iteration (k) to new values at iteration (k+1)
    for(int i=0;i<_N;i++) 
      for(int q=0;q<_Q;q++) 
	_tau1[i][q] = newtau1[i][q];
    for(int t=0;t<_T-1;t++) 
      for(int i=0;i<_N;i++) 
	for(int q=0;q<_Q;q++) 
	  for(int l=0;l<_Q;l++) 
	    _taut[t][i][q][l] = newtaut[t][i][q][l];
    deallocate2D<double>(newtau1,_N,_Q);
    deallocate4D<double>(newtaut,_T-1,_N,_Q,_Q);
    correctTau1();
    correctTaut();
    updateTauMarginal();
  }
    
  template<typename Ytype>
  template<class TAddEventFunctor> 
  void DynSBM<Ytype>::updateThetaCore(Ytype*** const Y, TAddEventFunctor& addEventFunctor){// M-step
    for(int t=0;t<_T;t++)
      for(int q=0;q<_Q;q++)
	for(int l=0;l<_Q;l++){
	  _betaql[t][q][l] = 0; // probability of a zero value
	}
    // note the constraint betaql[t,q,q] = betaql[,q,q]
    double*** betasumql;
    allocate3D<double>(betasumql,_T,_Q,_Q);
    // mu/beta
    for(int t=0;t<_T;t++){
      for(int i=0;i<_N;i++){
	if (ispresent(t,i)){
	  // j<i
	  for(int j=0;j<i;j++){
	    if (ispresent(t,j))
	      for(int q=0;q<_Q;q++){
		for(int l=0;l<q;l++){
		  if(_isdirected){		      
		    if(Y[t][i][j]>Ytype(0)){
		      addEventFunctor(tauMarginal(t,i,q)*tauMarginal(t,j,l), Y[t][i][j], t, q, l);
		      addEventFunctor(tauMarginal(t,i,l)*tauMarginal(t,j,q), Y[t][i][j], t, l, q);
		    }
		    else{
		      _betaql[t][q][l] += tauMarginal(t,i,q)*tauMarginal(t,j,l);
		      _betaql[t][l][q] += tauMarginal(t,i,l)*tauMarginal(t,j,q);
		    }
		    betasumql[t][q][l] += tauMarginal(t,i,q)*tauMarginal(t,j,l);
		    betasumql[t][l][q] += tauMarginal(t,i,l)*tauMarginal(t,j,q);
		    if(Y[t][j][i]>Ytype(0)){
		      addEventFunctor(tauMarginal(t,j,q)*tauMarginal(t,i,l), Y[t][j][i], t, q, l);
		      addEventFunctor(tauMarginal(t,j,l)*tauMarginal(t,i,q), Y[t][j][i], t, l, q);
		    } else{
		      _betaql[t][q][l] += tauMarginal(t,j,q)*tauMarginal(t,i,l);
		      _betaql[t][l][q] += tauMarginal(t,j,l)*tauMarginal(t,i,q);
		    }
		    betasumql[t][q][l] += tauMarginal(t,j,q)*tauMarginal(t,i,l);
		    betasumql[t][l][q] += tauMarginal(t,j,l)*tauMarginal(t,i,q);
		  } else{
		    if(Y[t][i][j]>Ytype(0)){
		      addEventFunctor(tauMarginal(t,i,q)*tauMarginal(t,j,l), Y[t][i][j], t, q, l);
		      addEventFunctor(tauMarginal(t,i,l)*tauMarginal(t,j,q), Y[t][i][j], t, q, l);
		    }
		    else{
		      _betaql[t][q][l] += tauMarginal(t,i,q)*tauMarginal(t,j,l) + tauMarginal(t,i,l)*tauMarginal(t,j,q);
		    }
		    betasumql[t][q][l] += tauMarginal(t,i,q)*tauMarginal(t,j,l) + tauMarginal(t,i,l)*tauMarginal(t,j,q);
		  }
		}
		// q==l
		if(Y[t][i][j]>Ytype(0))
		  addEventFunctor(tauMarginal(t,i,q)*tauMarginal(t,j,q), Y[t][i][j], 0, q, q);
		else{
		  _betaql[0][q][q] += tauMarginal(t,i,q)*tauMarginal(t,j,q); // t=0 gathers all t when q=l
		}
		if(_isdirected)
		  if(Y[t][j][i]>Ytype(0))
		    addEventFunctor(tauMarginal(t,i,q)*tauMarginal(t,j,q), Y[t][j][i], 0, q, q);
		  else
		    _betaql[0][q][q] += tauMarginal(t,i,q)*tauMarginal(t,j,q);		  
		betasumql[0][q][q] += (1+int(_isdirected)) * tauMarginal(t,i,q)*tauMarginal(t,j,q);	
	      }
	  }
	  // j==i considered only if selfloop allowed
	  if (_withselfloop)
	    for(int q=0;q<_Q;q++){
	      if(Y[t][i][i]>Ytype(0))
		addEventFunctor(tauMarginal(t,i,q), Y[t][i][i], 0, q, q);
	      else
		_betaql[0][q][q] += tauMarginal(t,i,q); // t=0 gathers all t when q=l
	      betasumql[0][q][q] += tauMarginal(t,i,q);	    
	    }	
	}	
      }
    }
    for(int t=0;t<_T;t++){// symmetrize+normalize
      for(int q=(_isdirected?0:1);q<_Q;q++){
	for(int l=0;l<q;l++)
	  if (betasumql[t][q][l]>0){
	    _betaql[t][q][l] = _betaql[t][q][l] / betasumql[t][q][l];
	    if(!_isdirected) _betaql[t][l][q] = _betaql[t][q][l];
	  }        	
	if(_isdirected)
	  for(int l=q+1;l<_Q;l++)
	    if (betasumql[t][q][l]>0)
	      _betaql[t][q][l] = _betaql[t][q][l] / betasumql[t][q][l];	      
      }
    }
    for(int q=0;q<_Q;q++){// symmetrize+normalize
      // note the constraint betaql[t,q,q] = betaql[,q,q]
      if (betasumql[0][q][q]>0)
	_betaql[0][q][q] = _betaql[0][q][q] / betasumql[0][q][q];
      for(int t=1;t<_T;t++)
	_betaql[t][q][q] = _betaql[0][q][q];
    }
    correctBeta();
    deallocate3D<double>(betasumql,_T,_Q,_Q);
  }
  
  template<typename Ytype>
  double DynSBM<Ytype>::completedLoglikelihood(Ytype*** const Y) const{ // including entropy term
    double J = 0.;
    // term 1
    for(int i=0;i<_N;i++){
      if (ispresent(0,i))
	for(int q=0;q<_Q;q++)
	  J += _tau1[i][q] * (log(_stationary[q])-log(_tau1[i][q]));
    }
    //std::cerr<<"\nJ: "<<J<<" ";
    // term 2
#pragma omp parallel for reduction(+:J)
    for(int t=1;t<_T;t++){
      int tstorage = t-1;
      for(int i=0;i<_N;i++)
	if(ispresent(t,i))
	  if(ispresent(t-1,i))
	    for(int q=0;q<_Q;q++)
	      for(int qprime=0;qprime<_Q;qprime++)
		J += tauMarginal(t-1,i,q) * _taut[tstorage][i][q][qprime]
		  * (log(_trans[q][qprime])-log(_taut[tstorage][i][q][qprime]));
	  else
	    for(int q=0;q<_Q;q++)
	      J += tauArrival(t,i,q) * (log(_stationary[q])-log(tauArrival(t,i,q)));	    
    }
    //std::cerr<<J<<" ";
    // term 3
#pragma omp parallel for reduction(+:J)
    for(int t=0;t<_T;t++)
      for(int i=0;i<_N;i++)
	if (ispresent(t,i)){
	  for(int j=0;j<i;j++)
	    if (ispresent(t,j))
	      for(int q=0;q<_Q;q++){
		double taumtiq = tauMarginal(t,i,q);
		for(int l=0;l<_Q;l++){
		  J += taumtiq*tauMarginal(t,j,l)*logDensity(t,q,l,Y[t][i][j]);
		  if (_isdirected) J += taumtiq*tauMarginal(t,j,l)*logDensity(t,l,q,Y[t][j][i]);
		}
	      }
	   if (_withselfloop){
	     for(int q=0;q<_Q;q++){
	       J += tauMarginal(t,i,q)*logDensity(t,q,q,Y[t][i][i]);
	     }
	   }
	}
    //std::cerr<<J<<std::endl;
    return(J);
  }
  
  template<typename Ytype>
  double DynSBM<Ytype>::modelselectionLoglikelihood(Ytype*** const Y) const{
    double LKL = 0;
    // term 1
    std::vector<int> groups1 = getGroupsByMAP(0);
    for(int i=0;i<_N;i++)
      if (ispresent(0,i))
	LKL += log(_stationary[groups1[i]]);
    //std::cerr<<"\nLKL: "<<LKL<<" ";
    // term 2
    std::vector<int> groupstm1 = groups1;
    for(int t=1;t<_T;t++){
      std::vector<int> groupst = getGroupsByMAP(t);
      for(int i=0;i<_N;i++)
	if(ispresent(t,i))
	  if(ispresent(t-1,i)){
	    LKL += log(_trans[groupstm1[i]][groupst[i]]);
	  }
	  else
	    LKL += log(_stationary[groupst[i]]);
      groupstm1 = groupst;
    }
    //std::cerr<<LKL<<" ";
    // term 3
    for(int t=0;t<_T;t++){
      std::vector<int> groupst = getGroupsByMAP(t);
      for(int i=0;i<_N;i++)
	if (ispresent(t,i)){
	  for(int j=0;j<i;j++)
	    if (ispresent(t,j)){
	      LKL += logDensity(t,groupst[i],groupst[j],Y[t][i][j]);
	      if (_isdirected) LKL += logDensity(t,groupst[j],groupst[i],Y[t][j][i]);
	    }	  
	  if (_withselfloop){
	    LKL += logDensity(t,groupst[i],groupst[i],Y[t][i][i]);
	  }
	}
    }
    //std::cerr<<LKL<<std::endl;
    return(LKL);
  }
  
  template<typename Ytype> 
  void  DynSBM<Ytype>::updateTauMarginal(){
    for(int t=1;t<_T;t++){
      int tstorage = t-1;
      for(int i=0;i<_N;i++)
	if (ispresent(t-1,i) && ispresent(t,i))
	  for(int qprime=0;qprime<_Q;qprime++){
	    _taum[tstorage][i][qprime] = 0.;
	    for(int q=0;q<_Q;q++)
	      _taum[tstorage][i][qprime] += tauMarginal(t-1,i,q)*_taut[tstorage][i][q][qprime];
	  } 
      correctTaum(); // TODO correctTaum(t) would be better
    }
  }
  
  template<typename Ytype> 
  double DynSBM<Ytype>::tauMarginal(int t, int i, int q) const{
    if(t==0){
      return(_tau1[i][q]);
    } else{
      if(ispresent(t-1,i)){
	int tstorage = t-1;
	return(_taum[tstorage][i][q]);
      } else
	return tauArrival(t,i,q);
       
    }
  }
  
  template<typename Ytype> 
  void DynSBM<Ytype>::correctTau1(){// numerical issue : avoid too small value for tau
    for(int i=0;i<_N;i++){
      double sumtau1i = 0.;
      for(int q=0;q<_Q;q++){
	if(_tau1[i][q]<1e-7){
	  _tau1[i][q] = 1e-7;
	}
	sumtau1i += _tau1[i][q];
      }
      for(int q=0;q<_Q;q++){
	_tau1[i][q] = _tau1[i][q]/sumtau1i;
      }
    }
  }
  
  template<typename Ytype> 
  void DynSBM<Ytype>::correctTaut(){// numerical issue : avoid too small value for tau
    for(int t=1;t<_T;t++){
      int tstorage = t-1;       
      for(int i=0;i<_N;i++){
	for(int q=0;q<(ispresent(t-1,i)?_Q:1);q++){ // only q=0 if absent at t-1 (see tauArrival method)
	  double sumtaut = 0.;
	  for(int qprime=0;qprime<_Q;qprime++){
	    if(_taut[tstorage][i][q][qprime]<1e-7){
	      _taut[tstorage][i][q][qprime] = 1e-7;
	    }
	    sumtaut += _taut[tstorage][i][q][qprime];
	  }
	  for(int qprime=0;qprime<_Q;qprime++){
	    _taut[tstorage][i][q][qprime] = _taut[tstorage][i][q][qprime]/sumtaut;
	  }
	}
      }
    }
  }
  
  template<typename Ytype> 
  void DynSBM<Ytype>::correctTaum(){// numerical issue : avoid too small value for tau
    for(int t=1;t<_T;t++){
      int tstorage = t-1;       
      for(int i=0;i<_N;i++){
	double sumtaum = 0.;
	for(int q=0;q<_Q;q++){
	  if(_taum[tstorage][i][q]<1e-7){
	    _taum[tstorage][i][q] = 1e-7;	    
	  }
	  sumtaum += _taum[tstorage][i][q];
	}
	for(int q=0;q<_Q;q++){
	  _taum[tstorage][i][q] = _taum[tstorage][i][q]/sumtaum;
	}
      }
    }
  }
  
  template<typename Ytype> 
  void DynSBM<Ytype>::initTau(const std::vector<int>& clustering){// initialization from a clustering
    // it supposes previous allocation with O values
    for(int i=0;i<_N;i++){
      _tau1[i][clustering[i]] = 1;
    }
    correctTau1();
    for(int t=1;t<_T;t++){
      int tstorage = t-1;
      for(int i=0;i<_N;i++){
	for(int q=0;q<_Q;q++){
	  _taut[tstorage][i][q][clustering[i]] = 1;
	}
      }
    }
    correctTaut();
    updateTauMarginal();
  }
  
  template<typename Ytype> 
  std::vector<int> DynSBM<Ytype>::getGroupsByMAP(int t) const{ // maximum a posteriori - classification step
    std::vector<int> groups(_N,-1);
    if(t==0){
      for(int i=0;i<_N;i++){
	if(ispresent(0,i)){
	  double maxi = 0.;
	  for(int q=0;q<_Q;q++){
	    if(_tau1[i][q]>maxi){
	      maxi = _tau1[i][q];
	      groups[i] = q;
	    }
	  }
	}
      }
    } else{
      for(int i=0;i<_N;i++){
	if(ispresent(t,i)){
	  double maxi = 0.;
	  for(int q=0;q<_Q;q++){
	    double tau = tauMarginal(t,i,q);
	    if(tau>maxi){
	      maxi = tau;
	      groups[i] = q;
	    }
	  }
	}
      }
    }
    return(groups);
  }
  
  template<typename Ytype>
  void DynSBM<Ytype>::correctBeta(){//numerical issue : avoid too small value for beta
    for(int t=0;t<_T;t++){
      for(int q=0;q<_Q;q++){
        for(int l=0;l<_Q;l++){
          if(_betaql[t][l][q]<1e-7){
            _betaql[t][l][q] = 1e-7;
          }
          if(_betaql[t][l][q]>(1-1e-7)){
            _betaql[t][l][q] = 1-1e-7;
          }
        }
      }
    }
  }
}
#endif
