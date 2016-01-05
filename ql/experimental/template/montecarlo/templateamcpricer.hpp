/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016, Sebastian Schlenkrich

*/

/*! \file templateamcpricer.hpp
    \brief American MC pricer for cancellable note
	
*/


#ifndef quantlib_templateamcpricer_hpp
#define quantlib_templateamcpricer_hpp


#include <ql/experimental/template/montecarlo/templatemcswap.hpp>


namespace QuantLib {

	template <class DateType, class PassiveType, class ActiveType>
	class TemplateAMCPricer {
	protected:
		// easy use of templated types
		typedef TemplateMCSimulation<DateType, PassiveType, ActiveType>                 SimulationType;
		typedef typename TemplateMCSimulation<DateType, PassiveType, ActiveType>::Path  PathType;
		typedef typename TemplateMC<DateType, PassiveType, ActiveType>::CancellableNote NoteType;


		// container class definitions
		typedef std::vector<DateType>                      VecD;
		typedef std::vector<PassiveType>                   VecP; 
		typedef std::vector<ActiveType>                    VecA;
		typedef std::vector< std::vector<DateType> >       MatD;
		typedef std::vector< std::vector<PassiveType> >    MatP;
		typedef std::vector< std::vector<ActiveType> >     MatA;

		// references
		boost::shared_ptr<NoteType>         note_;
		boost::shared_ptr<SimulationType>   simulation_;

		// data [paths][exercises[+1]][size]
		MatA                X_;    // discounted (!) accumulated coupons  X_i(T_i)/B(T_i)
		MatA                R_;    // discounted (!) early redemptions
		MatA                B_;    // numeraire at each exercise date
		std::vector<MatA>   xi_;   // regression variables at each exercise

		// we separate the regression operator from the AMC algorithm...

		inline void resizeData() {
			X_.clear();
			R_.clear();
			B_.clear();
			xi_.clear();
			size_t M = simulation_->nPaths();
			size_t N = note_->callTimes().size();
			X_.resize(M);
			R_.resize(M);
			B_.resize(M);
			xi_.resize(M);
			for (size_t k=0; k<M; ++k) {
				X_[k].resize(N+1,0.0);
				R_[k].resize(N,0.0);
				B_[k].resize(N,0.0);
				xi_[k].resize(N);
			}
			for (size_t j=0; j<N; ++j) {
				size_t nRegr = note_->regressionVariables()[j]->size();
				for (size_t k=0; k<M; ++k) {
					xi_[k][j].resize(nRegr,0.0);
				}
			}
		}

		inline void calculateData() {   // this routine does the actual MC valuations			
			for (size_t k=0; k<simulation_->nPaths(); ++k) {
				boost::shared_ptr<PathType> p = simulation_->path(k);

				// check each leg for coupons prior to exercise date
				for (size_t i=0; i<note_->underlyings().size(); ++i) {
					if (!note_->underlyings()[i]->size())                 continue;  // skip the leg if it is empty
					if (note_->underlyings()[i]->back()->startTime()<0.0) continue;  // skip the leg if there is no future coupon
					size_t idx = 0;
					while((idx<note_->underlyings()[i]->size())&&((*note_->underlyings()[i])[idx]->startTime()<0.0)) ++idx;
					// now we have identified the future coupons and we can accumulate coupons...
					for (size_t j=0; j<note_->callTimes().size(); ++j) {
						while((idx<note_->underlyings()[i]->size())&&((*note_->underlyings()[i])[idx]->startTime()<note_->callTimes()[j])) {
							X_[k][j] += (*note_->underlyings()[i])[idx]->discountedAt(p);
							++idx;
						}
					}
					// finally sum up coupons after the last call date
					while(idx<note_->underlyings()[i]->size()) {
						X_[k][note_->callTimes().size()] += (*note_->underlyings()[i])[idx]->discountedAt(p);
						++idx;
					}
				}

				// calculate early redemptions
				for (size_t j=0; j<note_->callTimes().size(); ++j) {
					for (size_t i=0; i<note_->earlyRedemptions()[j]->size(); ++i) {
						R_[k][j] += (*note_->earlyRedemptions()[j])[i]->discountedAt(p);
					}
				}

				// calculate numeraire
				for (size_t j=0; j<note_->callTimes().size(); ++j) {
					B_[k][j] = p->numeraire(note_->callTimes()[j]);
				}

				// calculate regression variables...

			}
		}


	public:

		TemplateAMCPricer( const boost::shared_ptr<NoteType>        note,
			               const boost::shared_ptr<SimulationType>  simulation
						   // maybe some more arguments to control AMC
						   )
						   : note_(note), simulation_(simulation) {
		    resizeData();
			calculateData();
		}


	};

}

#endif  /* ifndef quantlib_templateamcpricer_hpp */ 
