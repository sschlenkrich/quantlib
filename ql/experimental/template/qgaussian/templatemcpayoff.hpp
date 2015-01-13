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
			Pricer (const std::vector< boost::shared_ptr<TemplateMCPayoff> >& payoffs, 
				    const boost::shared_ptr<SimulationType>&                  simulation)
					: payoffs_(payoffs), simulation_(simulation) { }
			inline ActiveType NPV() {
				ActiveType npv = 0.0;
				for (size_t k=0; k<payoffs_.size(); ++k) {
					ActiveType npvk = 0;
					for (size_t n=0; n<simulation_->nPaths(); ++n) {
						npvk += payoffs_[k]->discountedAt(simulation_->path(n));
					}
					npv += npvk;
				}
				return npv / simulation_->nPaths();
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
				if (payTime_==observationTime()) return (ActiveType)1.0;
				return p->zeroBond(observationTime(),payTime_);
			}
		};


	};



}

#endif  /* ifndef quantlib_templatemcpayoff_hpp */
