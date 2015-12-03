/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file templatemcsimulation.hpp
    \brief simulate and store paths of

	Process dX(t) = a(t,X) dt + b(t,X) dW
	
	Methods drift(t,X) and diffusion(t,X) are provided by ProcessType
	
*/


#ifndef quantlib_templatemcsimulation_hpp
#define quantlib_templatemcsimulation_hpp

#include <ql/math/randomnumbers/rngtraits.hpp>
#include <ql/experimental/template/templatestochasticprocess.hpp>

#include <ql/experimental/template/auxilliaries/templateauxilliaries.hpp>



namespace QuantLib {

	class TemplateSimulation : public virtual Observable { };

	// Declaration of MC simulation class
	template <class DateType, class PassiveType, class ActiveType>
	class TemplateMCSimulation : public TemplateSimulation {
	protected:

		typedef TemplateStochasticProcess<DateType, PassiveType, ActiveType> ProcessType;

		// container class definitions
		typedef std::vector<DateType>                      VecD;
		typedef std::vector<PassiveType>                   VecP; 
		typedef std::vector<ActiveType>                    VecA;
		typedef std::vector< std::vector<DateType> >       MatD;
		typedef std::vector< std::vector<PassiveType> >    MatP;
		typedef std::vector< std::vector<ActiveType> >     MatA;

		// default pseudo random numbers
		boost::shared_ptr<PseudoRandom::rsg_type> rsg_;

		// the link to the process/model
		boost::shared_ptr<ProcessType>            process_;

		// time grid for process simulation
        VecD simTimes_;
		// time grid for process observation, subset of simTimes (store state at obsTimes_)
		VecD obsTimes_;
		// allow time interpolation on path
		bool timeInterpolation_;

		// Monte Carlo properties
		BigNatural  seed_;
		bool        richardsonExtrapolation_;

		// we may precompute and store the Brownian increments
		bool storeBrownians_;
		// Brownian increments dW (nPaths x k*nSimTimes x nFactors), k = 2 for Richardson Extrapolation
		std::vector< MatD > dW_;

		// paths stored X_[paths][obsTimes][size]
		std::vector< MatA > X_;

		// adjust simulated zcb by exp{-adj(t,t+dt)(t+dt)} to meet initial yield curve
		bool applyZcbAdjuster_;
		VecD zcbObservTimes_;
		VecD zcbOffsetTimes_;
		MatA zcbAdjuster_;

		// adjust simulated asset to S(t) + adj(t) to meet E[S(t)]
		bool applyAssetAdjuster_;
		VecD assetObservTimes_;
		VecA assetAdjuster_;

		// dummy
		PassiveType tmp_;

		// initialisation

		inline void reallocateMemory(size_t nPaths) {
			// clear/free memory
			dW_.clear();
			X_.clear();
			// dW_
			if (storeBrownians_) {
				dW_.resize(nPaths);
				size_t ndWt = ((richardsonExtrapolation_) ? 2 : 1) * (simTimes_.size()-1);
                for (size_t i=0; i<nPaths; ++i) {
					dW_[i].resize(ndWt);
					for (size_t j=0; j<ndWt; ++j) dW_[i][j].resize(process_->factors());
				}
			}
			// X_
			X_.resize(nPaths);
			for (size_t i=0; i<nPaths; ++i) {
				X_[i].resize(obsTimes_.size());
				for (size_t j=0; j<obsTimes_.size(); ++j) X_[i][j].resize(process_->size());
			}
			return;
		}

		inline void initialiseRSG() {
			size_t ndWt = ((richardsonExtrapolation_) ? 2 : 1) * (simTimes_.size()-1);
			rsg_ = boost::shared_ptr<PseudoRandom::rsg_type>(new PseudoRandom::rsg_type(PseudoRandom::make_sequence_generator(ndWt * process_->factors(), seed_)));
		}

		// return the next Brownian motion path increments from random sequence generator
		inline const std::vector< std::vector<QuantLib::Real> > getNextBrownianIncrements() {
			std::vector<QuantLib::Real> sequence = rsg_->nextSequence().value;
			size_t ndWt = ((richardsonExtrapolation_) ? 2 : 1) * (simTimes_.size()-1);
			size_t nFac = process_->factors();
			QL_REQUIRE(sequence.size()==ndWt*nFac,"TemplateMCSimulation: dimension mismatch");
			std::vector< std::vector<QuantLib::Real> > dW(ndWt);
			for (size_t i=0; i<ndWt; ++i) {
				dW[i].resize(nFac);
				for (size_t j=0; j<nFac; ++j) {
					dW[i][j] = sequence[i*nFac+j];
				}
			}
			return dW;
		}

		// return the Brownian motion path increments from from cache or random sequence generator
		inline const std::vector< std::vector<QuantLib::Real> > getBrownianIncrements(const size_t path) {
			if (storeBrownians_) {
				if (path<dW_.size()) return dW_[path];
				QL_REQUIRE(false,"TemplateMCSimulation: path index out of bounds.");
			}
			return getNextBrownianIncrements();
		}

		// cache Brownian motion path increments
		inline void preEvaluateBrownians() {
			if (!storeBrownians_) return; 
			for (size_t k=0; k<dW_.size(); ++k) dW_[k] = getNextBrownianIncrements();
			return;
		}

		// simulate a single path X_[path]
		inline void simulatePath(const size_t path) {
			QL_REQUIRE(path<X_.size(),"TemplateMCSimulation: path index out of bounds.");
			MatD dWt = getBrownianIncrements(path);
			// initialisation
			VecA X0 = process_->initialValues();
			VecA X1(X0.size()), X12(X0.size());
			VecD dW( process_->factors() );
			X_[path][0]    = X0;
			size_t obs_idx = 1;
			for (size_t sim_idx=1; sim_idx<simTimes_.size(); ++sim_idx) {
				DateType dt = simTimes_[sim_idx]-simTimes_[sim_idx-1];
				if (richardsonExtrapolation_) {
					// full Euler step
					for (size_t k=0; k<dW.size(); ++k) dW[k] = (dWt[2*(sim_idx-1)][k] + dWt[2*(sim_idx-1)+1][k])/sqrt(2.0);
					//X1 = evolve(simTimes_[sim_idx-1],X0,dt,dW);
					process_->evolve(simTimes_[sim_idx-1],X0,dt,dW,X1);
					// two half size Euler steps
					//VecA X12 = evolve(simTimes_[sim_idx-1],X0,dt/2.0, dWt[2*(sim_idx-1)]);
					//X0 = evolve(simTimes_[sim_idx-1]+dt/2.0,X12,dt/2.0,dWt[2*(sim_idx-1)+1]);
					process_->evolve(simTimes_[sim_idx-1],        X0,  dt/2.0, dWt[2*(sim_idx-1)],   X12);
					process_->evolve(simTimes_[sim_idx-1]+dt/2.0, X12, dt/2.0, dWt[2*(sim_idx-1)+1], X0);
					// extrapolation
					for (size_t k=0; k<X1.size(); ++k) X1[k] = 2*X0[k] - X1[k];
					// extrapolation may lead to ill-defined states
					process_->truncate( simTimes_[sim_idx], X1 );
				} else { // only full Euler step
					process_->evolve(simTimes_[sim_idx-1],X0,dt,dWt[sim_idx-1], X1);
				}
				if ((obs_idx<obsTimes_.size()) && (simTimes_[sim_idx]==obsTimes_[obs_idx])) {
					X_[path][obs_idx]=X1;
					++obs_idx;
				}
				X0 = X1;
			}
		}


	public:

		// model type definition
		//typedef ProcessType process_type;

		// input for payoffs
		// combine model/process and simulated paths
		class Path {
		protected:
			boost::shared_ptr<ProcessType>           process_;
			TemplateMCSimulation*                    sim_;
			size_t                                   idx_;
		public:
			Path(const boost::shared_ptr<ProcessType>           process,
				 TemplateMCSimulation*                          sim,
				 const size_t                                   idx )
				 : process_(process), sim_(sim), idx_(idx) {}

			inline ActiveType numeraire(DateType obsTime) {
				return process_->numeraire(obsTime,sim_->state(idx_,obsTime));
			}

			inline ActiveType zeroBond(DateType obsTime, DateType payTime) {
				return sim_->zcbAdjuster(obsTime, payTime) *
					   process_->zeroBond(obsTime, payTime, sim_->state(idx_,obsTime) );
			}

			inline ActiveType asset(DateType obsTime) {
				return sim_->assetAdjuster(obsTime) +
					   process_->asset(obsTime, sim_->state(idx_,obsTime) );
			}

			inline ActiveType future(DateType obsTime, DateType settlementTime) {
				return sim_->assetAdjuster(obsTime) +
					   process_->future(obsTime, settlementTime, sim_->state(idx_,obsTime) );
			}
		};


		TemplateMCSimulation( const boost::shared_ptr<ProcessType> process,
			                  const VecD&                          simTimes,
							  const VecD&                          obsTimes,
							  size_t                               nPaths,
							  BigNatural                           seed = 1234,
							  bool                                 richardsonExtrapolation = true,
							  bool                                 timeInterpolation = false,
							  bool                                 storeBrownians = false )
							  : process_(process), seed_(seed),
							  richardsonExtrapolation_(richardsonExtrapolation),
							  timeInterpolation_(timeInterpolation), storeBrownians_(storeBrownians) {
			// check inputs
			QL_REQUIRE(simTimes.size()>0,"TemplateMCSimulation: non-empty simulation times required");
			for (size_t k=1; k<simTimes.size(); ++k) QL_REQUIRE(simTimes[k-1]<simTimes[k],"TemplateMCSimulation: simulation times in ascending order required");
			QL_REQUIRE(obsTimes.size()>0,"TemplateMCSimulation: non-empty observation times required");
			for (size_t k=1; k<obsTimes.size(); ++k) QL_REQUIRE(obsTimes[k-1]<obsTimes[k],"TemplateMCSimulation: observation times in ascending order required");
			// use only obsTimes>0
			obsTimes_.push_back(0.0);
			for (size_t k=0; k<obsTimes.size(); ++k) if (obsTimes[k]>0) obsTimes_.push_back(obsTimes[k]);
			// merge simTimes and obsTimes_ into simTimes_
			size_t sim_idx=0;
			while ((sim_idx<simTimes.size()) && (simTimes[sim_idx]<=0)) ++sim_idx;
			size_t obs_idx=0;
			while ((obs_idx<obsTimes_.size()) && (obsTimes_[obs_idx]<=0)) ++obs_idx;
			// we alway start at t=0
			simTimes_.push_back(0.0);
			// we use some tolerance for checking equal dates to avoid rounding issues
			// maybe even a larger tolerance than one day could make sense
			DateType OneDay = 1.0/365.25;
			while ((sim_idx<simTimes.size()) || (obs_idx<obsTimes_.size())) {
				if  ((sim_idx<simTimes.size())&&(obs_idx<obsTimes_.size())&&
					(abs(simTimes[sim_idx]-obsTimes_[obs_idx])<OneDay)) {
					simTimes_.push_back(obsTimes_[obs_idx]);
					++sim_idx;
					++obs_idx;
					continue;
				}
				if ((obs_idx>=obsTimes_.size())||((sim_idx<simTimes.size())&&(simTimes[sim_idx]<obsTimes_[obs_idx]))) {
					simTimes_.push_back(simTimes[sim_idx]);
					++sim_idx;
					continue;
				} else {
					simTimes_.push_back(obsTimes_[obs_idx]);
					++obs_idx;
					continue;
				}
			}
			// allocate memory (we don't want surprises during time-consuming simulation)
			reallocateMemory(nPaths);
			// ready to simulate...
			// but first set up adjusters
			applyZcbAdjuster_ = false;     // default
			applyAssetAdjuster_  = false;  // default
		}

		inline void simulate() {  // the procedure Initialise/PreEvBrownians/Simulate needs to be re-factorised
			initialiseRSG();
			preEvaluateBrownians();
			for (size_t k=0; k<X_.size(); ++k) simulatePath(k);
		}

		// inspectors

		inline const boost::shared_ptr<ProcessType> process() { return process_; }

		inline const VecD& simTimes() { return simTimes_; }
		inline const VecD& obsTimes() { return obsTimes_; }
		inline size_t      nPaths()   { return X_.size(); }

		inline const boost::shared_ptr<Path> path(const size_t idx) {
			return boost::shared_ptr<Path>(new Path(process_,this,idx));
		}

		inline const MatA& observedPath(const size_t idx) {
			QL_REQUIRE(idx<X_.size(),"TemplateMCSimulation: path out of bounds.");
			return X_[idx]; 
		}

		inline const VecA state(const size_t idx, const DateType t) {
			QL_REQUIRE(idx<X_.size(),"TemplateMCSimulation: path out of bounds.");
			size_t t_idx = TemplateAuxilliaries::idx(obsTimes_,t);
			if (t==obsTimes_[t_idx]) return X_[idx][t_idx];
			QL_REQUIRE(timeInterpolation_,"TemplateMCSimulation: time interpolation not allowed");
			// allow extrapolation
			if (t<obsTimes_[0])                   return X_[idx][0];
			if (t>obsTimes_[obsTimes_.size()-1])  return X_[idx][X_[idx].size()-1];
			VecA X(X_[idx][t_idx-1]);
			// linear state interpolation (this is very crude) 
			DateType rho = (t-obsTimes_[t_idx-1])/(obsTimes_[t_idx] - obsTimes_[t_idx-1]);
			for (size_t k=0; k<X.size(); ++k) X[k] = (1.0-rho)*X[k] + rho*X_[idx][t_idx][k];
			return X;
		}

		inline const MatA& brownian(const size_t idx) {
			QL_REQUIRE(idx<dW_.size(),"TemplateMCSimulation: path out of bounds.");
			return dW_[idx]; 
		}

		// zero coupon bond adjuster

		inline void calculateZCBAdjuster( const VecD&  zcbObservTimes,
							              const VecD&  zcbOffsetTimes ) {
            zcbObservTimes_ = zcbObservTimes;
			zcbOffsetTimes_ = zcbOffsetTimes;
			// check time grids
			QL_REQUIRE(zcbObservTimes_.size()>1,"TemplateMCSimulation: at least two zcbObservTimes_ required");
			QL_REQUIRE(zcbObservTimes_[0]>=0,"TemplateMCSimulation: zcbObservTimes_[0]>=0 required");
			for (size_t k=1; k<zcbObservTimes_.size(); ++k) QL_REQUIRE(zcbObservTimes_[k-1]<zcbObservTimes_[k],"TemplateMCSimulation: zcbObservTimes_ in ascending order required");
			QL_REQUIRE(zcbOffsetTimes_.size()>1,"TemplateMCSimulation: at least two zcbOffsetTimes_ required");
			QL_REQUIRE(zcbOffsetTimes_[0]>=0,"TemplateMCSimulation: zcbOffsetTimes_[0]>=0 required");
			for (size_t k=1; k<zcbOffsetTimes_.size(); ++k) QL_REQUIRE(zcbOffsetTimes_[k-1]<zcbOffsetTimes_[k],"TemplateMCSimulation: zcbOffsetTimes_ in ascending order required");
			// initialise zero adjuster matrix
			zcbAdjuster_ = MatA(zcbObservTimes_.size(),VecA(zcbOffsetTimes_.size(),0.0));
			MatA zcb(zcbObservTimes_.size(),VecA(zcbOffsetTimes_.size(),0.0));
			for (size_t k=0; k<nPaths(); ++k) {
				for (size_t i=0; i<zcbObservTimes_.size(); ++i) {
					VecA       s   = state(k,zcbObservTimes_[i]);
					ActiveType num = process_->numeraire(zcbObservTimes_[i],s);
					for (size_t j=0; j<zcbOffsetTimes_.size(); ++j) {
						zcb[i][j] += process_->zeroBond(zcbObservTimes_[i], zcbObservTimes_[i]+zcbOffsetTimes_[j], s ) / num;
					}
				}
			}
			for (size_t i=0; i<zcbObservTimes_.size(); ++i) {
				for (size_t j=0; j<zcbOffsetTimes_.size(); ++j) {
					ActiveType adjDF = process_->zeroBond(0.0, zcbObservTimes_[i]+zcbOffsetTimes_[j], process_->initialValues() );
					// we ommited division by nuneraire(0) = 1 here
					zcb[i][j] /= nPaths();
					adjDF /= zcb[i][j];
					zcbAdjuster_[i][j] = - log(adjDF) / (zcbObservTimes_[i]+zcbOffsetTimes_[j]);
				}
			}
			applyZcbAdjuster_ = true;  // ready to use zcb adjuster
		}

		inline const MatA& zcbAdjuster() { return zcbAdjuster_; }

		inline ActiveType zcbAdjuster(const DateType t, const DateType T) {
			if (!applyZcbAdjuster_) return 1.0; // implement adjuster interpolation
			DateType dt = T - t;
			//if (dt<0) return 1.0;
			// bilinear interpolation
			size_t obsIdx = TemplateAuxilliaries::idx(zcbObservTimes_,t);
			size_t offIdx = TemplateAuxilliaries::idx(zcbOffsetTimes_,dt);
			if (obsIdx<1) obsIdx=1;
			if (offIdx<1) offIdx=1;
			DateType rhoObs = ( t-zcbObservTimes_[obsIdx-1])/(zcbObservTimes_[obsIdx]-zcbObservTimes_[obsIdx-1]);
			DateType rhoOff = (dt-zcbOffsetTimes_[offIdx-1])/(zcbOffsetTimes_[offIdx]-zcbOffsetTimes_[offIdx-1]);
			// flat extrapolation
			if (rhoObs<0) rhoObs = 0;
			if (rhoObs>1) rhoObs = 1;
			if (rhoOff<0) rhoOff = 0;
			if (rhoOff>1) rhoOff = 1;
			//
			ActiveType z = zcbAdjuster_[obsIdx-1][offIdx-1] * (1.0-rhoObs) * (1.0-rhoOff) +
				           zcbAdjuster_[obsIdx  ][offIdx-1] * (rhoObs)     * (1.0-rhoOff) +
				           zcbAdjuster_[obsIdx-1][offIdx  ] * (1.0-rhoObs) * (rhoOff)     +
				           zcbAdjuster_[obsIdx  ][offIdx  ] * (rhoObs)     * (rhoOff)     ;
			return exp(-z*T);
		}

		// asset adjuster

		inline void calculateAssetAdjuster( const VecD&  assetObservTimes ) {
			assetObservTimes_ = assetObservTimes;
			// check time grids
			QL_REQUIRE(assetObservTimes_.size()>1,"TemplateMCSimulation: at least two assetObservTimes_ required");
			QL_REQUIRE(assetObservTimes_[0]>=0,"TemplateMCSimulation: assetObservTimes_>=0 required");
			for (size_t k=1; k<assetObservTimes_.size(); ++k) QL_REQUIRE(assetObservTimes_[k-1]<assetObservTimes_[k],"TemplateMCSimulation: assetObservTimes_ in ascending order required");
			// initialise asset adjuster vector
			assetAdjuster_ = VecA(assetObservTimes_.size(), 0.0);
			VecA avAsset(assetObservTimes_.size(), 0.0);
			VecA avZero(assetObservTimes_.size(), 0.0);
			for (size_t i=0; i<assetObservTimes_.size(); ++i) {
                for (size_t k=0; k<nPaths(); ++k) {
					VecA       s   = state(k,assetObservTimes_[i]);
					avAsset[i] += process_->asset(assetObservTimes_[i],s) / process_->numeraire(assetObservTimes_[i],s);
					avZero[i]  += 1.0                                     / process_->numeraire(assetObservTimes_[i],s);
				}
				avAsset[i] /= nPaths();
				avZero[i]  /= nPaths();
				// expected asset S(0) / Num(0) = E[ S(t) / Num(t) ]
				assetAdjuster_[i] = process_->asset(0, process_->initialValues() ) / process_->numeraire( 0, process_->initialValues() );
				assetAdjuster_[i] = (assetAdjuster_[i] - avAsset[i]) / avZero[i];
			}
			applyAssetAdjuster_ = true;
		}

		inline const VecA& assetAdjuster() { return assetAdjuster_; }

		inline ActiveType assetAdjuster(const DateType t) {
			if (!applyAssetAdjuster_) return 0.0; 
			// bilinear interpolation
			size_t obsIdx = TemplateAuxilliaries::idx(assetObservTimes_,t);
			if (obsIdx<1) obsIdx=1;
			DateType rhoObs = ( t-assetObservTimes_[obsIdx-1])/(assetObservTimes_[obsIdx]-assetObservTimes_[obsIdx-1]);
			// flat extrapolation
			if (rhoObs<0) rhoObs = 0;
			if (rhoObs>1) rhoObs = 1;
			ActiveType adj = (1.0-rhoObs) * assetAdjuster_[obsIdx-1] + rhoObs * assetAdjuster_[obsIdx];
			return adj;
		}

	};

}

#endif  /* ifndef quantlib_templatemcsimulation_hpp */
