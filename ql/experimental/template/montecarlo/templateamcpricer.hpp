/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016, Sebastian Schlenkrich

*/

/*! \file templateamcpricer.hpp
    \brief American MC pricer for cancellable note
	
*/


#ifndef quantlib_templateamcpricer_hpp
#define quantlib_templateamcpricer_hpp


#include <ql/experimental/template/montecarlo/templatemc.hpp>
#include <ql/experimental/template/auxilliaries/templateregression.hpp>



namespace QuantLib {

	template <class DateType, class PassiveType, class ActiveType>
	class TemplateAMCPricer {
	protected:
		// easy use of templated types
		typedef TemplateMCSimulation<DateType, PassiveType, ActiveType>                 SimulationType;
		typedef typename TemplateMCSimulation<DateType, PassiveType, ActiveType>::Path  PathType;
		typedef typename TemplateMC<DateType, PassiveType, ActiveType>::CancellableNote NoteType;
		typedef typename TemplateAuxilliaries::Regression<PassiveType>                  RegressionType;


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

		// MC data [paths][exercises[+1]]
		MatA                X_;    // discounted (!) accumulated coupons  X_i(T_i)/B(T_i)
		MatA                R_;    // discounted (!) early redemptions
		MatA                B_;    // numeraire at each exercise date

		// [exercises][paths][size]
		std::vector<MatA>   xi_;   // regression variables at each exercise

		// AMC (roll back) data [exercises][paths]
		MatA                G_;      // canellable note value if not called yet
		MatA                T_;      // regressed hold trigger variable Regr(G - R)

		//size_t              M_;      // number of paths
		//size_t              N_;      // number of exercises

		// we separate the regression operator from the AMC algorithm...
		PassiveType                                        regressionFraction_;  // partition of paths used to calculate regression
		size_t                                             maxPolynDegree_;
		std::vector< boost::shared_ptr<RegressionType> >   regressions_;

		inline void checkNote() { // check that note is set up correctly
			QL_REQUIRE(note_->earlyRedemptions().size()   ==note_->callTimes().size(), "AMC error: wrong number of redemption legs.");
			QL_REQUIRE(note_->regressionVariables().size()==note_->callTimes().size(), "AMC error: wrong number of regression var's legs.");
		}

		inline void resizeData(const size_t nPaths) {
			X_.clear();
			R_.clear();
			B_.clear();
			xi_.clear();
			G_.clear();
			T_.clear();

			// set sizes for this object
			size_t nCalls = note_->callTimes().size();

			X_.resize(nPaths);
			R_.resize(nPaths);
			B_.resize(nPaths);
			for (size_t k=0; k<nPaths; ++k) {
				X_[k].resize(nCalls+1,0.0); // we need data before and after all call dates
				R_[k].resize(nCalls,0.0);
				B_[k].resize(nCalls,0.0);
			}

			xi_.resize(nCalls);
			G_.resize(nCalls);
			T_.resize(nCalls);
			for (size_t k=0; k<nCalls; ++k) {
				xi_[k].resize(nPaths);
				size_t nRegr = note_->regressionVariables()[k]->size();
				for (size_t j=0; j<nPaths; ++j) {
					xi_[k][j].resize(nRegr,0.0);
				}
				G_[k].resize(nPaths,0.0);
				T_[k].resize(nPaths,0.0);
			}

		}

		inline void calculateData( const size_t idxStart,
			                       const size_t idxEnd    ) {   // this routine does the actual MC valuations
            size_t nPaths = idxEnd - idxStart;
			QL_REQUIRE(idxEnd<=simulation_->nPaths(),"AMC error: wrong idx");  // basic consistency check
			QL_REQUIRE(nPaths==X_.size(),"AMC error: wrong dimensions");  // basic consistency check
			size_t nCalls = note_->callTimes().size();
			for (size_t k=0; k<nPaths; ++k) {
				boost::shared_ptr<PathType> p = simulation_->path(idxStart+k);

				// check each leg for coupons prior to exercise date
				for (size_t i=0; i<note_->underlyings().size(); ++i) {
					if (!note_->underlyings()[i]->size())                 continue;  // skip the leg if it is empty
					if (note_->underlyings()[i]->back()->startTime()<0.0) continue;  // skip the leg if there is no future coupon
					size_t idx = 0;
					while((idx<note_->underlyings()[i]->size())&&((*note_->underlyings()[i])[idx]->startTime()<0.0)) ++idx;
					// now we have identified the future coupons and we can accumulate coupons...
					for (size_t j=0; j<nCalls; ++j) {
						while((idx<note_->underlyings()[i]->size())&&((*note_->underlyings()[i])[idx]->startTime()<note_->callTimes()[j])) {
							X_[k][j] += (*note_->underlyings()[i])[idx]->discountedAt(p);
							++idx;
						}
					}
					// finally sum up coupons after the last call date
					while(idx<note_->underlyings()[i]->size()) {
						X_[k][nCalls] += (*note_->underlyings()[i])[idx]->discountedAt(p);
						++idx;
					}
				}

				// calculate early redemptions
				for (size_t j=0; j<nCalls; ++j) {
					for (size_t i=0; i<note_->earlyRedemptions()[j]->size(); ++i) {
						R_[k][j] += (*note_->earlyRedemptions()[j])[i]->discountedAt(p);
					}
				}

				// calculate numeraire
				for (size_t j=0; j<nCalls; ++j) {
					B_[k][j] = p->numeraire(note_->callTimes()[j]);
				}

				// calculate regression variables...
				for (size_t j=0; j<nCalls; ++j) {
					for (size_t i=0; i<xi_[j][k].size(); ++i) {  
						xi_[j][k][i] = (*note_->regressionVariables()[j])[i]->at(p);
					}
				}

			}  // for each path
		}

		inline void rollBack(const bool calculateRegression) {
			size_t nPaths = X_.size();
			size_t nCalls = note_->callTimes().size();
			// roll back
			for (size_t j=nCalls; j>0; --j) {				
				for (size_t k=0; k<nPaths; ++k) {  // initialisation
					G_[j-1][k] = B_[k][j-1] * X_[k][j];  // use X_[k][j] here!
				}				
				if (j<nCalls) {  // update based on call trigger and future value
					for (size_t k=0; k<nPaths; ++k) {
 					    G_[j-1][k] += (T_[j][k]>0.0) ? (B_[k][j-1]/B_[k][j]*G_[j][k]) : (B_[k][j-1]*R_[k][j]);
					}
				}
				for (size_t k=0; k<nPaths; ++k) {  // calculate new trigger based on future data		
					T_[j-1][k] = G_[j-1][k] - B_[k][j-1]*R_[k][j-1];
				}
				if (calculateRegression) { // regression is not (re-)calculated in valuation run
					regressions_[j-1] = boost::shared_ptr<RegressionType>( new RegressionType(xi_[j-1],T_[j-1],maxPolynDegree_) );
				}
				if (regressions_[j-1]) {  // if there is no regression we look into the future
					for (size_t k=0; k<nPaths; ++k) T_[j-1][k] = regressions_[j-1]->value(xi_[j-1][k]);
				}
			}
		}

	public:

		TemplateAMCPricer( const boost::shared_ptr<NoteType>        note,
			               const boost::shared_ptr<SimulationType>  simulation,
						   const PassiveType                        regressionFraction,
						   const size_t                             maxPolynDegree
						   // maybe some more arguments to control AMC
						   )
						   : note_(note), simulation_(simulation), maxPolynDegree_(maxPolynDegree)  {
		    regressions_.resize(note_->callTimes().size());
			regressionFraction_ = regressionFraction;
			if (regressionFraction_<0.0) regressionFraction_ = 0.0;
			if (regressionFraction_>1.0) regressionFraction_ = 1.0;
		}

		inline void calculate() {
			checkNote();
			size_t idxSplit = (size_t)(regressionFraction_ * simulation_->nPaths());
			if (idxSplit>0) {
				resizeData(idxSplit);
				calculateData(0,idxSplit);
				rollBack(true);
			}
			if (idxSplit<simulation_->nPaths()) {
				resizeData(simulation_->nPaths()-idxSplit);
				calculateData(idxSplit,simulation_->nPaths());
				rollBack(false);
			}
		}

		inline const ActiveType noteNPV() const {
			size_t nPaths = X_.size();
			size_t nCalls = (X_.size()>0) ? (X_[0].size()-1) : (0);
			ActiveType res = 0.0;
			for (size_t k=0; k<nPaths; ++k) {
				res += X_[k][0];
			}
			if (nCalls>0) {
				for (size_t k=0; k<nPaths; ++k) {
					res += (T_[0][k]>0.0) ? (G_[0][k]/B_[k][0]) : (R_[k][0]);
				}
			}
			return res / nPaths;
		}

		inline const ActiveType underlyingNPV() const {
			size_t nPaths = X_.size();
			size_t nCalls = (X_.size()>0) ? (X_[0].size()-1) : (0);
			ActiveType res = 0.0;
			for (size_t k=0; k<nPaths; ++k) {
				for (size_t j=0; j<nCalls+1; ++j) {
					res += X_[k][j];
				}
			}
			return res / nPaths;
		}

		inline const ActiveType noteOptionNPV() const {
			return noteNPV() - underlyingNPV();
		}

		// inspector
		inline const MatA& data(const std::string& tag ) const {
			if (tag.compare("X")==0) return X_;
			if (tag.compare("R")==0) return R_;
			if (tag.compare("B")==0) return B_;
			if (tag.compare("G")==0) return G_;
			if (tag.compare("T")==0) return T_;
			return X_; // default
		}

	};

}

#endif  /* ifndef quantlib_templateamcpricer_hpp */ 
