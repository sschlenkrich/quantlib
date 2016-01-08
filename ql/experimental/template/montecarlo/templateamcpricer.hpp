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

		size_t              M_;      // number of paths
		size_t              N_;      // number of exercises

		// we separate the regression operator from the AMC algorithm...
		bool                calculateRegression_;  
		size_t              maxPolynDegree_;
		std::vector< boost::shared_ptr<RegressionType> >   regressions_;

		inline void checkNote() { // check that note is set up correctly
			QL_REQUIRE(note_->earlyRedemptions().size()   ==note_->callTimes().size(), "AMC error: wrong number of redemption legs.");
			QL_REQUIRE(note_->regressionVariables().size()==note_->callTimes().size(), "AMC error: wrong number of regression var's legs.");
		}

		inline void resizeData() {
			X_.clear();
			R_.clear();
			B_.clear();
			xi_.clear();
			G_.clear();
			T_.clear();

			// set sizes for this object
			M_ = simulation_->nPaths();
			N_ = note_->callTimes().size();

			X_.resize(M_);
			R_.resize(M_);
			B_.resize(M_);
			for (size_t k=0; k<M_; ++k) {
				X_[k].resize(N_+1,0.0); // we need data before and after all call dates
				R_[k].resize(N_,0.0);
				B_[k].resize(N_,0.0);
			}

			xi_.resize(N_);
			G_.resize(N_);
			T_.resize(N_);
			for (size_t k=0; k<N_; ++k) {
				xi_[k].resize(M_);
				size_t nRegr = note_->regressionVariables()[k]->size();
				for (size_t j=0; j<M_; ++j) {
					xi_[k][j].resize(nRegr,0.0);
				}
				G_[k].resize(M_,0.0);
				T_[k].resize(M_,0.0);
			}

		}

		inline void calculateData() {   // this routine does the actual MC valuations			
			for (size_t k=0; k<M_; ++k) {
				boost::shared_ptr<PathType> p = simulation_->path(k);

				// check each leg for coupons prior to exercise date
				for (size_t i=0; i<note_->underlyings().size(); ++i) {
					if (!note_->underlyings()[i]->size())                 continue;  // skip the leg if it is empty
					if (note_->underlyings()[i]->back()->startTime()<0.0) continue;  // skip the leg if there is no future coupon
					size_t idx = 0;
					while((idx<note_->underlyings()[i]->size())&&((*note_->underlyings()[i])[idx]->startTime()<0.0)) ++idx;
					// now we have identified the future coupons and we can accumulate coupons...
					for (size_t j=0; j<N_; ++j) {
						while((idx<note_->underlyings()[i]->size())&&((*note_->underlyings()[i])[idx]->startTime()<note_->callTimes()[j])) {
							X_[k][j] += (*note_->underlyings()[i])[idx]->discountedAt(p);
							++idx;
						}
					}
					// finally sum up coupons after the last call date
					while(idx<note_->underlyings()[i]->size()) {
						X_[k][N_] += (*note_->underlyings()[i])[idx]->discountedAt(p);
						++idx;
					}
				}

				// calculate early redemptions
				for (size_t j=0; j<N_; ++j) {
					for (size_t i=0; i<note_->earlyRedemptions()[j]->size(); ++i) {
						R_[k][j] += (*note_->earlyRedemptions()[j])[i]->discountedAt(p);
					}
				}

				// calculate numeraire
				for (size_t j=0; j<N_; ++j) {
					B_[k][j] = p->numeraire(note_->callTimes()[j]);
				}

				// calculate regression variables...
				for (size_t j=0; j<N_; ++j) {
					for (size_t i=0; i<xi_[j][k].size(); ++i) {  
						xi_[j][k][i] = (*note_->regressionVariables()[j])[i]->at(p);
					}
				}

			}  // for each path
		}

		inline void rollBack() {
			// roll back
			for (size_t j=N_; j>0; --j) {				
				for (size_t k=0; k<M_; ++k) {  // initialisation
					G_[j-1][k] = B_[k][j-1] * X_[k][j];  // use X_[k][j] here!
				}				
				if (j<N_) {  // update based on call trigger and future value
					for (size_t k=0; k<M_; ++k) {
 					    G_[j-1][k] += (T_[j][k]>0.0) ? (B_[k][j-1]/B_[k][j]*G_[j][k]) : (B_[k][j-1]*R_[k][j]);
					}
				}
				for (size_t k=0; k<M_; ++k) {  // calculate new trigger based on future data		
					T_[j-1][k] = G_[j-1][k] - B_[k][j-1]*R_[k][j-1];
				}
				if (calculateRegression_) { // regression is not (re-)calculated in valuation run
					regressions_[j-1] = boost::shared_ptr<RegressionType>( new RegressionType(xi_[j-1],T_[j-1],maxPolynDegree_) );
				}
				if (regressions_[j-1]) {  // if there is no regression we look into the future
					for (size_t k=0; k<M_; ++k) T_[j-1][k] = regressions_[j-1]->value(xi_[j-1][k]);
				}

				// regress trigger based on j-1 information
				// ...
			}
		}

	public:

		TemplateAMCPricer( const boost::shared_ptr<NoteType>        note,
			               const boost::shared_ptr<SimulationType>  simulation,
						   const bool                               calculateRegression,
						   const size_t                             maxPolynDegree
						   // maybe some more arguments to control AMC
						   )
						   : note_(note), simulation_(simulation), M_(0), N_(0),
						   calculateRegression_(calculateRegression), maxPolynDegree_(maxPolynDegree)  {
		    regressions_.resize(N_);
		}

		inline void calculate() {
			checkNote();
		    resizeData();
			calculateData();
			rollBack();
		}

		inline const ActiveType noteNPV() const {
			ActiveType res = 0.0;
			for (size_t k=0; k<M_; ++k) {
				res += X_[k][0];
			}
			if (N_>0) {
				for (size_t k=0; k<M_; ++k) {
					res += (T_[0][k]>0.0) ? (G_[0][k]/B_[k][0]) : (R_[k][0]);
				}
			}
			return res / M_;
		}

		inline const ActiveType underlyingNPV() const {
			ActiveType res = 0.0;
			for (size_t k=0; k<M_; ++k) {
				for (size_t j=0; j<N_+1; ++j) {
					res += X_[k][j];
				}
			}
			return res / simulation_->nPaths();
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
