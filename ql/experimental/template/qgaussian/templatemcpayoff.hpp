/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file templatemcpayoff.hpp
    \brief generic payoff interface for MC simulation
	
*/


#ifndef quantlib_templatemcpayoff_hpp
#define quantlib_templatemcpayoff_hpp


#include <ql/experimental/template/qgaussian/templatequasigaussian.hpp>
#include <ql/experimental/template/qgaussian/templatemcsimulation.hpp>



namespace QuantLib {

	// Base class for template payoffs
    template <class DateType, class PassiveType, class ActiveType, class SimulationType, class PathType>
	class TemplateMCPayoff {
	protected:
	    DateType observationTime_;
	public:
	    TemplateMCPayoff( const DateType observationTime ) : observationTime_(observationTime) { }
		// inspectors
		inline DateType observationTime() { return observationTime_; }
		// generic payoff(observationTime, p) needs to be implemented by derived classes
        virtual ActiveType at(const boost::shared_ptr<PathType>& p) = 0;
		// discounted payoff for NPV valuation
		inline ActiveType discountedAt(const boost::shared_ptr<PathType>& p) { return at(p) / p->numeraire(observationTime_); }

		// generic pricer
		class Pricer {
		protected:
			std::vector< boost::shared_ptr<TemplateMCPayoff> > payoffs_;
			boost::shared_ptr<SimulationType>                  simulation_;
		public:

			inline static ActiveType NPV(const std::vector< boost::shared_ptr<TemplateMCPayoff> >& payoffs, 
				                         const boost::shared_ptr<SimulationType>&                  simulation) {
				ActiveType npv = 0.0;
				for (size_t k=0; k<payoffs.size(); ++k) {
					ActiveType npvk = 0;
					for (size_t n=0; n<simulation->nPaths(); ++n) {
						npvk += payoffs[k]->discountedAt(simulation->path(n));
					}
					npv += npvk;
				}
				return npv / simulation->nPaths();
			}


			Pricer (const std::vector< boost::shared_ptr<TemplateMCPayoff> >& payoffs, 
				    const boost::shared_ptr<SimulationType>&                  simulation)
					: payoffs_(payoffs), simulation_(simulation) { }

			inline ActiveType NPV() {
				/*
				ActiveType npv = 0.0;
				for (size_t k=0; k<payoffs_.size(); ++k) {
					ActiveType npvk = 0;
					for (size_t n=0; n<simulation_->nPaths(); ++n) {
						npvk += payoffs_[k]->discountedAt(simulation_->path(n));
					}
					npv += npvk;
				}
				return npv / simulation_->nPaths();
				*/
				return NPV(payoffs_,simulation_);
			}
		};


		// particular payoffs

		// simple cash payment
		class Cash : public TemplateMCPayoff {
		protected:
			DateType payTime_;
		public:
			Cash( DateType obsTime, DateType payTime ) : TemplateMCPayoff(obsTime), payTime_(payTime) { }
			virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if (payTime_<=observationTime()) return (ActiveType)1.0;
				return p->zeroBond(observationTime(),payTime_);
			}
		};

		// annuity
		class Annuity : public TemplateMCPayoff {
		protected:
			std::vector<DateType>    payTimes_;
			std::vector<PassiveType> payWeights_;  // these are typically year fractions
		public:
			Annuity( DateType                        obsTime, 
				     const std::vector<DateType>&    payTimes,
				     const std::vector<PassiveType>& payWeights)
				: TemplateMCPayoff(obsTime), payTimes_(payTimes), payWeights_(payWeights) { }
			virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				ActiveType res = 0.0;
				size_t N = (payTimes_.size()<payWeights_.size()) ? (payTimes_.size()) : (payWeights_.size());
				for (size_t k=0; k<N; ++k) {
					if (payTimes_[k]>observationTime()) {
						res += payWeights_[k] * p->zeroBond(observationTime(),payTimes_[k]);
					}
				}
				return res;
			}
		};

		// prototypical physically settled European swaption
		class ModelSwaption : public TemplateMCPayoff {
		protected:
			std::vector<DateType>    times_;        // T_0, .., T_N
			std::vector<PassiveType> payWeights_;   // tau_0, .., tau_N-1
			PassiveType              strikeRate_;   // option strike
			PassiveType              payOrRec_;     // call (+1) or put (-1) option on swap rate
			bool                     isConsistent_; // check consistency of input
		public:
			ModelSwaption( DateType                        obsTime, 
				           const std::vector<DateType>&    times,
				           const std::vector<PassiveType>& payWeights,
					       PassiveType                     strikeRate,
					       PassiveType                     payOrRec      )
				: TemplateMCPayoff(obsTime), times_(times), payWeights_(payWeights), strikeRate_(strikeRate), payOrRec_(payOrRec) { 
				isConsistent_ = true;
				if (times_.size()<2) isConsistent_ = false;
				for (size_t k=0; k<times_.size(); ++k) if (times_[k]<observationTime()) isConsistent_ = false;
				// default weights
				if (payWeights_.size()!=times_.size()-1) {
					payWeights_.resize(times_.size()-1);
					for (size_t k=0; k<payWeights_.size(); ++k) payWeights_[k] = times_[k+1] - times_[k];
				}
				// finished
			}
			virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				ActiveType res = 0.0;
				if (!isConsistent_) return res;
				// annuity...
				for (size_t k=0; k<payWeights_.size(); ++k) res += payWeights_[k] * p->zeroBond(observationTime(),times_[k+1]);
				// floatleg - fixedleg...
				res = p->zeroBond(observationTime(),times_[0]) - p->zeroBond(observationTime(),times_[times_.size()-1]) - strikeRate_*res;
				// payer or receiver swap...
				res *= payOrRec_;
				// exercise option...
				res = (res>0) ? res : 0.0;
				return res;
			}
		};


	};



}

#endif  /* ifndef quantlib_templatemcpayoff_hpp */
