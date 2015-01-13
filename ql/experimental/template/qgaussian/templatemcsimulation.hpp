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

#include <ql/experimental/template/qgaussian/templatequasigaussian.hpp>



namespace QuantLib {

	class TemplateSimulation : public virtual Observable { };

	// Declaration of MC simulation class
	template <class DateType, class PassiveType, class ActiveType, class ProcessType>
	class TemplateMCSimulation : public TemplateSimulation {
	protected:

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
		std::vector< std::vector< std::vector<QuantLib::Real> > > dW_;

		// paths stored X_[paths][obsTimes][size]
		std::vector< MatA > X_;

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

		// integrate X1 = X0 + drift()*dt + diffusion()*dW*sqrt(dt)
		inline VecA evolve( const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const std::vector<QuantLib::Real> & dW  ) {
			VecA X1 = X0;
			VecA a = process_->drift(t0, X0);
			MatA b = process_->diffusion(t0, X0);
			for (size_t i=0; i<X1.size(); ++i) {
				ActiveType tmp = 0.0;
				for (size_t j=0; j<dW.size(); ++j) tmp += b[i][j]*dW[j];
				X1[i] += a[i]*dt + tmp*sqrt(dt);
			}
			return X1;
		}

		// simulate a single path X_[path]
		inline void simulatePath(const size_t path) {
			QL_REQUIRE(path<X_.size(),"TemplateMCSimulation: path index out of bounds.");
			std::vector< std::vector<QuantLib::Real> > dWt = getBrownianIncrements(path);
			// initialisation
			VecA X0 = process_->initialValues();
			X_[path][0]    = X0;
			size_t obs_idx = 1;
			for (size_t sim_idx=1; sim_idx<simTimes_.size(); ++sim_idx) {
				DateType dt = simTimes_[sim_idx]-simTimes_[sim_idx-1];
				VecA X1;
				if (richardsonExtrapolation_) {
					// full Euler step
					std::vector<QuantLib::Real> dW(dWt[2*(sim_idx-1)]);
					for (size_t k=0; k<dW.size(); ++k) dW[k] = (dW[k] + dWt[2*(sim_idx-1)+1][k])/sqrt(2.0);
					X1 = evolve(simTimes_[sim_idx-1],X0,dt,dW);
					// two half size Euler steps
					VecA X12 = evolve(simTimes_[sim_idx-1],X0,dt/2.0, dWt[2*(sim_idx-1)]);
					X0 = evolve(simTimes_[sim_idx-1]+dt/2.0,X12,dt/2.0,dWt[2*(sim_idx-1)+1]);
					// extrapolation
					for (size_t k=0; k<X1.size(); ++k) X1[k] = 2*X0[k] - X1[k];
				} else { // only full Euler step
					X1 = evolve(simTimes_[sim_idx-1],X0,dt,dWt[sim_idx-1]);
				}
				if (simTimes_[sim_idx]==obsTimes_[obs_idx]) {
					X_[path][obs_idx]=X1;
					++obs_idx;
				}
				X0 = X1;
			}
		}

	public:

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
			while ((sim_idx<simTimes.size()) || (obs_idx<obsTimes_.size())) {
				if  ((sim_idx<simTimes.size())&&(obs_idx<obsTimes_.size())&&(simTimes[sim_idx]==obsTimes_[obs_idx])) {
					simTimes_.push_back(simTimes[sim_idx]);
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
		}

		inline void simulate() {  // the procedure Initialise/PreEvBrownians/Simulate needs to be re-factorised
			initialiseRSG();
			preEvaluateBrownians();
			for (size_t k=0; k<X_.size(); ++k) simulatePath(k);
		}

		// inspectors

		inline const VecD& simTimes() { return simTimes_; }
		inline const VecD& obsTimes() { return obsTimes_; }
		inline size_t      nPaths()   { return X_.size(); }     

		inline const MatA& path(const size_t path) {
			QL_REQUIRE(path<X_.size(),"TemplateMCSimulation: path out of bounds.");
			return X_[path]; 
		}

		inline const MatA& brownian(const size_t path) {
			QL_REQUIRE(path<dW_.size(),"TemplateMCSimulation: path out of bounds.");
			return dW_[path]; 
		}


	};

}

#endif  /* ifndef quantlib_templatemcsimulation_hpp */
