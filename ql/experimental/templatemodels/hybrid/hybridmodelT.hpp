/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019, Sebastian Schlenkrich

*/

/*! \file hybridmodelT.hpp
    \brief (MC) pricing for multi-asset model with stochastic interest rates
	           
*/


#ifndef quantlib_templatehybridmodel_hpp
#define quantlib_templatehybridmodel_hpp

#include <ql/errors.hpp>

#include <ql/experimental/templatemodels/stochasticprocessT.hpp>
#include <ql/experimental/templatemodels/hybrid/assetmodelT.hpp>
#include <ql/experimental/templatemodels/qgaussian2/quasigaussianmodel2T.hpp>

//#include <ql/experimental/templatemodels/auxilliaries/qrfactorisationT.hpp>
//#include <ql/experimental/templatemodels/auxilliaries/choleskyfactorisationT.hpp>


namespace QuantLib {

	// Declaration of the asset model class
	template <class DateType, class PassiveType, class ActiveType>
	class HybridModelT : public StochasticProcessT<DateType, PassiveType, ActiveType> {
	protected:

		// container class definitions
		typedef std::vector<DateType>                      VecD;
		typedef std::vector<PassiveType>                   VecP; 
		typedef std::vector<ActiveType>                    VecA;
		typedef std::vector< std::vector<DateType> >       MatD;
		typedef std::vector< std::vector<PassiveType> >    MatP;
		typedef std::vector< std::vector<ActiveType> >     MatA;

		typedef AssetModelT<DateType, PassiveType, ActiveType>          AssetModelType;
		typedef QuasiGaussianModel2T<DateType, PassiveType, ActiveType> RatesModelType;

		// class members
		std::string                      domAlias_;
		ext::shared_ptr<RatesModelType>  domRatesModel_;
		std::vector<std::string>         forAliases_;
		std::vector< ext::shared_ptr<AssetModelType> >  forAssetModels_;
		std::vector< ext::shared_ptr<RatesModelType> >  forRatesModels_;
		MatA                             correlations_;
		MatA                             L_;   // factorised correlation matrix

		size_t size_;     // the number of state variables
		size_t factors_;  // the number of random factors

		// keep the first index in state Y for the foreign (!) models
		// [ asset1, rates1, asset2, rates2, ...]
		// domestic rates model state starts at Y[0] and is not stored explicitely
		std::vector<size_t>  modelsStartIdx_;

		// aliases are keys
		std::map<std::string, size_t>  index_;

	public:  
		HybridModelT(const std::string                        domAlias,
			         const ext::shared_ptr<RatesModelType>    domRatesModel,
			         const std::vector<std::string>&          forAliases,
			         const std::vector< ext::shared_ptr<AssetModelType> >& forAssetModels,
			         const std::vector< ext::shared_ptr<RatesModelType> >& forRatesModels,
			         const MatA&                                           correlations)
			: domAlias_(domAlias), domRatesModel_(domRatesModel), forAliases_(forAliases),
			  forAssetModels_(forAssetModels), forRatesModels_(forRatesModels), correlations_(correlations),
			  modelsStartIdx_(2*forAliases.size()) {
			// we perform a couple of consistency checks
			QL_REQUIRE(domAlias_.compare("") != 0, "HybridModel: Domestic alias required.");
			QL_REQUIRE(domRatesModel_, "HybridModel: Domestic rates model required.");
			QL_REQUIRE(forAliases_.size()>0, "HybridModel: forAliases.size()>0 required.");
			QL_REQUIRE(forAliases_.size()==forAssetModels_.size(), "HybridModel: forAliases.size()==forAssetModels.size() required.");
			QL_REQUIRE(forAliases_.size()==forRatesModels_.size(), "HybridModel: forAliases.size()==forRatesModels.size() required.");
			for (size_t k = 0; k < forAliases_.size(); ++k) index_[forAliases_[k]] = k;
			// we also check the individual components to avoid surprises at later stages...
			for (size_t k = 0; k < forAliases.size(); ++k) {
				QL_REQUIRE(forAliases_[k].compare("") != 0, "HybridModel: Foreign alias required.");
				QL_REQUIRE(forAssetModels_[k], "HybridModel: Foreign asset model required.");
				QL_REQUIRE(forRatesModels_[k], "HybridModel: Foreign rates model required.");
			}
			// manage model indices in state variable
			size_ = domRatesModel_->size();
			factors_ = domRatesModel_->factors();
			size_t lastModelIdx = domRatesModel_->size();
			for (size_t k = 0; k < forAliases.size(); ++k) {
				size_ += (forAssetModels_[k]->size() + forRatesModels_[k]->size());
				factors_ += (forAssetModels_[k]->factors() + forRatesModels_[k]->factors());
				modelsStartIdx_[2*k]   = lastModelIdx;					
				modelsStartIdx_[2*k+1] = lastModelIdx + forAssetModels_[k]->size();
				lastModelIdx += (forAssetModels_[k]->size() + forRatesModels_[k]->size());
			}
			// we need to check and factor the hybrid correlation matrix
			QL_REQUIRE(correlations_.size()==factors(), "HybridModel: correlations_.size()==factors() required.");
			for (size_t k = 0; k < correlations_.size(); ++k) {
				QL_REQUIRE(correlations_[k].size() == factors(), "HybridModel: correlations_[k].size()==factors() required.");
				QL_REQUIRE(correlations_[k][k] == 1.0, "HybridModel: correlations_[k][K] == 1.0 required.");
			}
			for (size_t k = 0; k < correlations_.size(); ++k) {
				for (size_t j = k+1; j < correlations_.size(); ++j)
					QL_REQUIRE(correlations_[k][j] == correlations_[j][k], "HybridModel: correlations_[k][j] == correlations_[j][k] required.");
			}
			// We should also ensure that rates model correlations are identity here.
			// Rates model correlations are incorporated via individual rates models.
			L_ = TemplateAuxilliaries::cholesky(correlations_);
		}

		// inspectors
		const std::string domAlias()                           { return domAlias_;       }
		const ext::shared_ptr<RatesModelType>& domRatesModel() { return domRatesModel_;  }
		const std::vector<std::string>&  forAliases()          { return forAliases_;     }
		const std::vector< ext::shared_ptr<AssetModelType> >& forAssetModels() { return forAssetModels_; }
		const std::vector< ext::shared_ptr<RatesModelType> >& forRatesModels() { return forRatesModels_; }
		const MatA& correlations()                             { return correlations_;   }
		const MatA& L()                                        { return L_;              }
		inline const std::vector<size_t>& modelsStartIdx()     { return modelsStartIdx_; }


        // stochastic process interface
		// dimension of combined state variable
		inline size_t size()    { return size_; }
		// stochastic factors
		inline size_t factors() { return factors_; }
		// initial values for simulation
		inline VecP initialValues() {
			VecP X(size(),0.0);
			// domestic model
			VecP x(domRatesModel_->initialValues());
			for (size_t i = 0; i < domRatesModel_->size(); ++i) X[i] = x[i];
			for (size_t k = 0; k < forAliases_.size(); ++k) {
				x = forAssetModels_[k]->initialValues();
				for (size_t i = 0; i < forAssetModels_[k]->size(); ++i) X[modelsStartIdx_[2 * k] + i] = x[i];
				x = forRatesModels_[k]->initialValues();
				for (size_t i = 0; i < forRatesModels_[k]->size(); ++i) X[modelsStartIdx_[2 * k + 1] + i] = x[i];
			}
			return X;
		}

		// a[t,X(t)]
		inline VecA drift( const DateType t, const VecA& X) {
			VecA a(size(),0.0);
			QL_FAIL("HybridModel: drift not implemented. Use evolve.");
			// finished
			return a;
		}

		// b[t,X(t)]
		inline MatA diffusion( const DateType t, const VecA& X) {
			MatA b(size(), VecA(factors(),0.0));
			QL_FAIL("HybridModel: diffusion not implemented. Use evolve.");
			// finished
			return b;
		}

		// simple Euler step
		inline void evolve(const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1) {
			// correlate Brownian increments
			VecD dZ(dW.size(),0.0);
			for (size_t i = 0; i < dW.size(); ++i) // exploit lower triangular structure
				for (size_t j = 0; j <= i; ++j) dZ[i] += L_[i][j] * dW[j];
			// evolve domestic rates
			VecD dw(dZ.begin(), dZ.begin() + domRatesModel_->factors());
			VecA x0(X0.begin(), X0.begin() + domRatesModel_->size());
			VecA x1(domRatesModel_->size(),0.0);
			domRatesModel_->evolve(t0, x0, dt, dw, x1);
			for (size_t i = 0; i < domRatesModel_->size(); ++i) X1[i] = x1[i];
			// we will need the deterministic drift part for r_d
			PassiveType B_d = domRatesModel_->termStructure()->discount(t0) / domRatesModel_->termStructure()->discount(t0 + dt);
			// we also need the stochastic drift part for r_d
			ActiveType x_av_d = 0.0;
			// NOTE: only add first (#factors) components of MC state x which represent state variable 'x'
			// This is a bit dangerous, we use structure of internal quasiGaussian model state here
			for (size_t i = 0; i < domRatesModel_->factors(); ++i) x_av_d += 0.5 * (x0[i] + x1[i]);
			// foreign models...
			size_t corrStartIdx = domRatesModel_->factors();
			for (size_t k = 0; k < forAliases_.size(); ++k) {
				// carefully collect Brownian increments
				VecD dw_asset(dZ.begin() + corrStartIdx, dZ.begin() + corrStartIdx + forAssetModels_[k]->factors());
				VecD dw_rates(dZ.begin() + corrStartIdx + forAssetModels_[k]->factors(),
					          dZ.begin() + corrStartIdx + forAssetModels_[k]->factors() + forRatesModels_[k]->factors());
                // we need the starting point states for evolution
				VecA y0(X0.begin() + modelsStartIdx_[2 * k], X0.begin() + modelsStartIdx_[2 * k] + forAssetModels_[k]->size());
				VecA x0(X0.begin() + modelsStartIdx_[2 * k + 1], X0.begin() + modelsStartIdx_[2 * k + 1] + forRatesModels_[k]->size());
				// quanto adjustment
				VecA qAdj(correlations_[corrStartIdx].begin() + corrStartIdx + forAssetModels_[k]->factors(),
					      correlations_[corrStartIdx].begin() + corrStartIdx + forAssetModels_[k]->factors() + forRatesModels_[k]->factors());
				ActiveType assetVol = forAssetModels_[k]->volatility(t0, y0);
				for (size_t i = 0; i < qAdj.size(); ++i) qAdj[i] *= (assetVol*std::sqrt(dt));
				for (size_t i = 0; i < qAdj.size(); ++i) dw_rates[i] -= qAdj[i];
				// evolve foreign rates
				VecA x1(forRatesModels_[k]->size(),0.0);
				forRatesModels_[k]->evolve(t0, x0, dt, dw_rates, x1);
				// calculate FX drift
				PassiveType B_f = forRatesModels_[k]->termStructure()->discount(t0) / forRatesModels_[k]->termStructure()->discount(t0 + dt);
				ActiveType mu = std::log(B_d / B_f) / dt;  // deterministic part; better cache to avoid repeated discount call for each path
				// stochastic drift part for r_f
				ActiveType x_av_f = 0.0;  // Note comments for x_av_d above
				for (size_t i = 0; i < forRatesModels_[k]->factors(); ++i) x_av_f += 0.5 * (x0[i] + x1[i]);
				// now we can calculate drift mu = r_d - r_f
				mu += (x_av_d - x_av_f);
				// ... and plug it into the drift component of asset MC state variable
				y0[1] = mu;  // NOTE: this uses knowledge of the structure of asset model mc state!
				// finally we can evolve FX
				VecA y1(forAssetModels_[k]->size(),0.0);
				forAssetModels_[k]->evolve(t0, y0, dt, dw_asset, y1);
				// almost done; we copy MC states back in hybrid MC state
				for (size_t i = 0; i < forAssetModels_[k]->size(); ++i) X1[modelsStartIdx_[2 * k] + i]     = y1[i];
				for (size_t i = 0; i < forRatesModels_[k]->size(); ++i) X1[modelsStartIdx_[2 * k + 1] + i] = x1[i];
				// done.

				// update stochastic factor index
				corrStartIdx += (forAssetModels_[k]->factors() + forRatesModels_[k]->factors());
			}
			//QL_FAIL("HybridModel: evolve not yet implemented.");
		}

		// asset calculation is the main purpose of this model
		inline virtual ActiveType asset(const DateType t, const VecA& Y, const std::string& alias) {
			size_t k = index_.at(alias);  // this should throw an exception if alias is unknown
			VecA y(Y.begin() + modelsStartIdx_[2 * k], Y.begin() + modelsStartIdx_[2 * k] + forAssetModels_[k]->size());
			return forAssetModels_[k]->asset(t, y, alias);
		}

		// a domestic currency zero coupon bond
		inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X) { 
			VecA x(X.begin(), X.begin() + domRatesModel_->size());
			return domRatesModel_->zeroBond(t, T, x);
		}

		// a foreign currency zero coupon bond
		inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X, const std::string& alias) { 
			if (alias.compare(domAlias_) == 0) return zeroBond(t, T, X);
			size_t k = index_.at(alias);  // this should throw an exception if alias is unknown
			VecA x(X.begin() + modelsStartIdx_[2 * k + 1], X.begin() + modelsStartIdx_[2 * k + 1] + forRatesModels_[k]->size());
			return forRatesModels_[k]->zeroBond(t, T, x);
		}



	};

}

#endif  /* ifndef quantlib_templatehybridmodel_hpp */
