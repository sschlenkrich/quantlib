/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich
*/

/*! \file tenoroptionletvts.hpp
    \brief caplet volatility term structure based on volatility transformation
*/

#ifndef quantlib_tenoroptionletvts_hpp
#define quantlib_tenoroptionletvts_hpp


#include <ql/termstructures/volatility/optionlet/optionletvolatilitystructure.hpp>
#include <ql/termstructures/volatility/smilesection.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/math/interpolation.hpp>
#include <ql/time/dategenerationrule.hpp>

//#include <ql/instruments/swaption.hpp>
//#include <ql/termstructures/yieldtermstructure.hpp>
//#include <ql/time/date.hpp>
//#include <ql/option.hpp>

namespace QuantLib {

	class TenorOptionletVTS : public OptionletVolatilityStructure {

	public: class CorrelationStructure;  // declaration below

	protected:

		class TenorOptionletSmileSection : public SmileSection {
		protected:
			boost::shared_ptr<CorrelationStructure>      correlation_;
		    std::vector<boost::shared_ptr<SmileSection>> baseSmileSection_;
			std::vector<Time>                            startTimeBase_;  // for correlation parametrisation
			std::vector<Real>                            fraRateBase_;
			Real                                         fraRateTarg_;
			std::vector<Real>                            v_;
			// implement transformation formula
            virtual Volatility volatilityImpl(Rate strike) const;
		public:
			// constructor includes actual transformation details
			TenorOptionletSmileSection( const TenorOptionletVTS& volTS,
				                        const Time               optionTime );

		    // further SmileSection interface methods
			virtual Real minStrike() const { return baseSmileSection_[0]->minStrike()+fraRateTarg_-fraRateBase_[0]; }
			virtual Real maxStrike() const { return baseSmileSection_[0]->maxStrike()+fraRateTarg_-fraRateBase_[0]; }
			virtual Real atmLevel() const  { return fraRateTarg_;  }
		};

		Handle<OptionletVolatilityStructure>            baseVTS_;
		boost::shared_ptr<IborIndex>                    baseIndex_;
		boost::shared_ptr<IborIndex>                    targIndex_;
		boost::shared_ptr<CorrelationStructure>         correlation_;

	public:

		// functor interface for parametric correlation
		class CorrelationStructure {
		public:
			// return the correlation between two FRA rates starting at start1 and start2
			virtual Real operator() (const Time& start1, const Time& start2) const = 0;
		};

		// very basic choice for correlation structure
		class TwoParameterCorrelation : public CorrelationStructure {
		protected:
			boost::shared_ptr<Interpolation> rhoInf_;
			boost::shared_ptr<Interpolation> beta_;
		public:
			TwoParameterCorrelation ( const boost::shared_ptr<Interpolation> rhoInf,
				                      const boost::shared_ptr<Interpolation> beta)
				                    : rhoInf_(rhoInf_), beta_(beta_) {}
			virtual Real operator() (const Time& start1, const Time& start2) const {
				return (*rhoInf_)(start1) * (1.0 - (*rhoInf_)(start1))*exp(-(*beta_)(start1)*fabs(start2-start1));
			}
		};

		// constructor
		TenorOptionletVTS( const Handle<OptionletVolatilityStructure>&            baseVTS,
			               const boost::shared_ptr<IborIndex>&                    baseIndex,
		                   const boost::shared_ptr<IborIndex>&                    targIndex,
						   const boost::shared_ptr<CorrelationStructure>&         correlation);  

		// Termstructure interface

		//! the latest date for which the curve can return values
		virtual Date maxDate() const { return baseVTS_->maxDate(); }
		
		// VolatilityTermstructure interface

        //! implements the actual smile calculation in derived classes
		virtual boost::shared_ptr<SmileSection> smileSectionImpl(Time optionTime) const { return boost::shared_ptr<SmileSection>(new TenorOptionletSmileSection(*this,optionTime)); }
        //! implements the actual volatility calculation in derived classes
		virtual Volatility volatilityImpl(Time optionTime,Rate strike) const { return smileSection(optionTime)->volatility(strike); }


		//! the minimum strike for which the term structure can return vols
		virtual Rate minStrike() const { return baseVTS_->minStrike(); }
        //! the maximum strike for which the term structure can return vols
        virtual Rate maxStrike() const {return baseVTS_->maxStrike(); }

	};

	typedef TenorOptionletVTS::CorrelationStructure TenorOptionletVTSCorrelationStructure; 

}

#endif /* #ifndef quantlib_tenoroptionletvts_hpp */
