/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/



#ifndef quantlib_templatevanillalocalvolmodel_hpp
#define quantlib_templatevanillalocalvolmodel_hpp

#include <vector>
#include <string>

#include <stdio.h>

#include <boost/shared_ptr.hpp>
#include <ql/errors.hpp>

#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>


namespace QuantLib {

	template <class DateType, class PassiveType, class ActiveType>    
	class VanillaLocalVolModelT {
	private:
		// input parameters
		DateType    T_;                    // time to expiry (in years)
		PassiveType S0_;                   // forward
		PassiveType sigmaATM_;             // ATM normal volatility as basis for straddle calculation
		std::vector<PassiveType> Sp_;      // S_i with S_i > S0
		std::vector<PassiveType> Sm_;      // S_-i with S_-i < S0
		std::vector<PassiveType> Mp_;      // slope on intervall [S_i-1, S_i)
		std::vector<PassiveType> Mm_;      // slope on intervall (S_-i, S_-[i-1]]
		// calculated parameters
		ActiveType  straddleATM_;          // ATM straddle price
		ActiveType  sigma0_;               // local vol at S0, i.e. sigma(S0)
		std::vector<ActiveType>  sigmaP_;  // sigma(Sp_[i])
		std::vector<ActiveType>  sigmaM_;  // sigma(Sm_[i])
		std::vector<PassiveType> Xp_;      // X_i with X_i > 0
		std::vector<PassiveType> Xm_;      // X_-i with X_-i < 0
		// adjusters
		PassiveType              mu_;      // in-the-model adjuster for forward
		PassiveType              alpha_;   // out-of-the-model adjuster for straddle
		PassiveType              nu_;      // out-of-the-model adjuster for forward
		// numerical accuracy parameters (maybe expose to user...)
		PassiveType              extrapolationStdevs_;  // number of stdevs used as lower and upper cutoff, default 10
		size_t                   maxCalibrationIters_;  // number of iterations for forward/sigma0 calibration
		PassiveType              sigma0Tol_;            // tolerance for sigma convergence
		PassiveType              S0Tol_;                // tolerance for forward convergence
		bool                     adjustATM_;            // apply post-calibration ATM adjuster
		// we may want some debug information for the calibration process
		bool                     enableLogging_;
		std::vector<std::string> logging_;
 
		// determine the lower and upper bounds for integration
		inline PassiveType lowerBoundX() { return  -extrapolationStdevs_ * sqrt(T_) + mu_; }
		inline PassiveType upperBoundX() { return   extrapolationStdevs_ * sqrt(T_) + mu_; }		

		inline PassiveType localVol(const bool isRightWing, const size_t k, const PassiveType S) {
			// this is an unsafe method specifying the vol function sigma(S) on the individual segments
			if (isRightWing) {
				QL_REQUIRE(k < sigmaP_.size(), "k < sigmaP_.size() required.");
				PassiveType sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
				PassiveType S0 = (k > 0) ? Sp_[k - 1] : S0_;
				return sigma0 + Mp_[k] * (S - S0);
			} else {
				QL_REQUIRE(k < sigmaM_.size(), "k < sigmaM_.size() required.");
				PassiveType sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
				PassiveType S0 = (k > 0) ? Sm_[k - 1] : S0_;
				return sigma0 + Mm_[k] * (S - S0);
			}
			return 0.0;  // this should never be reached
		}

		inline PassiveType underlyingS(const bool isRightWing, const size_t k, const PassiveType x) {
			// this is an unsafe method specifying the underlying level S(x) on the individual segments
			if (isRightWing) {
				QL_REQUIRE(k < sigmaP_.size(), "k < sigmaP_.size() required.");
				PassiveType sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
				PassiveType x0 = (k > 0) ? Xp_[k - 1] : 0.0;
				PassiveType deltaS = (Mp_[k] == 0.0) ? (sigma0*(x - x0)) : (sigma0 / Mp_[k] * (exp(Mp_[k] * (x - x0)) - 1.0));
				PassiveType S0 = (k > 0) ? Sp_[k - 1] : S0_;
				return S0 + deltaS;
			}
			else {
				QL_REQUIRE(k < sigmaM_.size(), "k < sigmaM_.size() required.");
				PassiveType sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
				PassiveType x0 = (k > 0) ? Xm_[k - 1] : 0.0;
				PassiveType deltaS = (Mm_[k] == 0.0) ? (sigma0*(x - x0)) : (sigma0 / Mm_[k] * (exp(Mm_[k] * (x - x0)) - 1.0));
				PassiveType S0 = (k > 0) ? Sm_[k - 1] : S0_;
				return S0 + deltaS;
			}
			return 0.0;  // this should never be reached
		}

		inline PassiveType underlyingX(const bool isRightWing, const size_t k, const PassiveType S) {
			// this is an unsafe method specifying the underlying level x(S) on the individual segments
			if (isRightWing) {
				QL_REQUIRE(k < sigmaP_.size(), "k < sigmaP_.size() required.");
				PassiveType sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
				QL_REQUIRE(sigma0 > 0.0, "sigma0 > 0.0 required");
				PassiveType x0 = (k > 0) ? Xp_[k - 1] : 0.0;
				PassiveType S0 = (k > 0) ? Sp_[k - 1] : S0_;
				PassiveType deltaX = (Mp_[k] == 0.0) ? ((S - S0) / sigma0) : (log(1.0 + Mp_[k] / sigma0*(S - S0)) / Mp_[k]);
				return x0 + deltaX;
			}
			else {
				QL_REQUIRE(k < sigmaM_.size(), "k < sigmaM_.size() required.");
				PassiveType sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
				QL_REQUIRE(sigma0 > 0.0, "sigma0 > 0.0 required");
				PassiveType x0 = (k > 0) ? Xm_[k - 1] : 0.0;
				PassiveType S0 = (k > 0) ? Sm_[k - 1] : S0_;
				PassiveType deltaX = (Mm_[k] == 0.0) ? ((S - S0) / sigma0) : (log(1.0 + Mm_[k] / sigma0*(S - S0)) / Mm_[k]);
				return x0 + deltaX;
			}
			return 0.0;  // this should never be reached
		}

		inline PassiveType primitiveF(const bool isRightWing, const size_t k, const PassiveType x) {
			// this is an unsafe method specifying the primitive function F(x) = \int [alpha S(x) + nu] p(x) dx
			// on the individual segments
			PassiveType sigma0, x0, S0, m0;
			if (isRightWing) {
				QL_REQUIRE(k < sigmaP_.size(), "k < sigmaP_.size() required.");
				sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
				x0 = (k > 0) ? Xp_[k - 1] : 0.0;
				S0 = (k > 0) ? Sp_[k - 1] : S0_;
				m0 = Mp_[k];
			}
			else {
				QL_REQUIRE(k < sigmaM_.size(), "k < sigmaM_.size() required.");
				sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
				x0 = (k > 0) ? Xm_[k - 1] : 0.0;
				S0 = (k > 0) ? Sm_[k - 1] : S0_;
				m0 = Mm_[k];
			}
			PassiveType	y0 = (x0 - mu_) / sqrt(T_);
			PassiveType y1 = (x - mu_) / sqrt(T_);
			PassiveType h = m0 * sqrt(T_);
			PassiveType Ny = TemplateAuxilliaries::Phi(y1);
			PassiveType term1, term2;
			if (m0 == 0.0) {
				term1 = (S0 + nu_ / alpha_ - sigma0 * sqrt(T_) * y0) * Ny;
				term2 = sigma0 * sqrt(T_) * TemplateAuxilliaries::phi(y1);  // use dN/dx = dN/dy / sqrt(T)
			}
			else {				
				PassiveType NyMinush = TemplateAuxilliaries::Phi(y1 - h);
				term1 = exp(h*h / 2.0 - h*y0)*sigma0 / m0 * NyMinush;
				term2 = (sigma0 / m0 - (S0 + nu_ / alpha_)) * Ny;
			}
			return alpha_ * (term1 - term2);
		}

		inline void updateLocalVol() {
			// use ODE solution to determine x-grid and sigma-grid taking into account constraints of
			// positive local volatility and local vol extrapolation
			for (size_t k = 0; k < Sp_.size(); ++k) { // right wing calculations
				PassiveType x0 = (k > 0) ? Xp_[k - 1] : 0.0;
				PassiveType S0 = (k > 0) ? Sp_[k - 1] : S0_;
				PassiveType sigma0 = (k > 0) ? sigmaP_[k - 1] : sigma0_;
				QL_REQUIRE(sigma0 > 0.0, "sigma0 > 0.0 required.");
				if ((k == Sp_.size() - 1)||(localVol(true, k, Sp_[k])<=0.0)||(underlyingX(true, k, Sp_[k])>upperBoundX()))  { // right wing extrapolation, maybe better use some epsilon here
					PassiveType XRight = upperBoundX();  // mu might not yet be calibrated
					QL_REQUIRE(XRight >= x0, "XRight >= x0 required.");
					Xp_[k] = XRight;
					Sp_[k] = underlyingS(true, k, XRight);
					sigmaP_[k] = localVol(true, k, Sp_[k]);
					if (k < Sp_.size() - 1) Mp_[k + 1] = Mp_[k];  // we need to make sure vol doesn't go up again
					continue;
				}
				sigmaP_[k] = localVol(true, k, Sp_[k]);
				QL_REQUIRE(sigmaP_[k] > 0.0, "sigmaP_[k] > 0.0 required.");
				Xp_[k] = underlyingX(true, k, Sp_[k]);
			}
			for (size_t k = 0; k < Sm_.size(); ++k) { // left wing calculations
				PassiveType x0 = (k > 0) ? Xm_[k - 1] : 0.0;
				PassiveType S0 = (k > 0) ? Sm_[k - 1] : S0_;
				PassiveType sigma0 = (k > 0) ? sigmaM_[k - 1] : sigma0_;
				QL_REQUIRE(sigma0 > 0.0, "sigma0 > 0.0 required.");
				if ((k == Sm_.size() - 1)||(localVol(false, k, Sm_[k]) <= 0.0)||(underlyingX(false, k, Sm_[k])<lowerBoundX())) { // left wing extrapolation, maybe better use some epsilon here
					PassiveType XLeft = lowerBoundX();  // mu might not yet be calibrated
					QL_REQUIRE(XLeft <= x0, "XLeft <= x0 required.");
					Xm_[k] = XLeft;
					Sm_[k] = underlyingS(false, k, XLeft);
					sigmaM_[k] = localVol(false, k, Sm_[k]);
					if (k < Sm_.size() - 1) Mm_[k + 1] = Mm_[k];  // we need to make sure vol doesn't go up again
					continue;
				}
				sigmaM_[k] = localVol(false, k, Sm_[k]);
				QL_REQUIRE(sigmaM_[k] > 0.0, "sigmaM_[k] > 0.0 required.");
				Xm_[k] = underlyingX(false, k, Sm_[k]);
			}
		}

		inline void calibrateATM() {
			PassiveType straddleVega = straddleATM_ / sigmaATM_;
			PassiveType forwardMinusStrike0, forwardMinusStrike1, straddleMinusATM0, straddleMinusATM1, dmu, dsigma0;
			PassiveType	dfwd_dmu, dstr_dsi;
			for (size_t k = 0; k < maxCalibrationIters_; ++k) {
				PassiveType call = expectation(true, S0_);
				PassiveType put = expectation(false, S0_);
				forwardMinusStrike1 = call - put;
				straddleMinusATM1 = call + put - straddleATM_;
				if (k == 0) {
					dfwd_dmu = sigma0_;        // this is an estimate based on dS/dX at ATM
					dstr_dsi = straddleVega;   // this is an estimate based on dsigmaATM / dsigma0 =~ 1
				}
				if (k>0) { // we use secant if available
					// only update derivative if we had a step, otherwise use from previous iteration
					if (fabs(dmu)>1.0e-16) dfwd_dmu = (forwardMinusStrike1 - forwardMinusStrike0) / dmu;
					if (fabs(dsigma0)>1.0e-16) dstr_dsi = (straddleMinusATM1 - straddleMinusATM0) / dsigma0;
				}
				dmu = -forwardMinusStrike1 / dfwd_dmu;
				dsigma0 = -straddleMinusATM1 / dstr_dsi;
				if ((sigma0_ + dsigma0) < 0.0) dsigma0 = -0.5 * sigma0_;  // make sure sigma0_ remains positive
				// maybe some line search could improve convergence...
				mu_ += dmu;
				sigma0_ += dsigma0;
				updateLocalVol();
				// prepare for next iteration
				forwardMinusStrike0 = forwardMinusStrike1;
				straddleMinusATM0 = straddleMinusATM1;
				if (enableLogging_) logging_.push_back("k: " + std::to_string(k) +
					"; C: " + std::to_string(call) +
					"; P: " + std::to_string(put) +
					"; S: " + std::to_string(straddleATM_) +
					"; dfwd_dmu: " + std::to_string(dfwd_dmu) +
					"; dstr_dsi: " + std::to_string(dstr_dsi) +
					"; dmu: " + std::to_string(dmu) +
					"; dsigma0: " + std::to_string(dsigma0));
				if ((fabs(forwardMinusStrike0) < S0Tol_) && (fabs(dsigma0) < sigma0Tol_)) break;
			}
		}

		inline void adjustATM() {
			// reset adjusters in case this method is invoked twice
			alpha_ = 1.0;
			nu_ = 0.0;
			PassiveType call0 = expectation(true, S0_);
			PassiveType put0 = expectation(false, S0_);
			nu_ = put0 - call0;
			if (enableLogging_) logging_.push_back("C0: " + std::to_string(call0) + "; P0: " + std::to_string(put0) + "; nu: " + std::to_string(nu_));
			PassiveType call1 = expectation(true, S0_);
			PassiveType put1 = expectation(false, S0_);
			alpha_ = straddleATM_ / (call1 + put1);
			nu_ = alpha_*nu_ + (1.0 - alpha_)*S0_;
			if (enableLogging_) logging_.push_back("C1: " + std::to_string(call1) + "; P1: " + std::to_string(put1) + "; alpha_: " + std::to_string(alpha_) + "; nu_: " + std::to_string(nu_));
		}

	public:
		VanillaLocalVolModelT(
			const DateType                   T,
			const PassiveType                S0,
			const PassiveType                sigmaATM,
			const std::vector<PassiveType>&  Sp,
			const std::vector<PassiveType>&  Sm,
			const std::vector<PassiveType>&  Mp,
			const std::vector<PassiveType>&  Mm,
			// controls for calibration
			const size_t                     maxCalibrationIters = 5,
			bool                             adjustATMFlag = true,
			bool                             enableLogging = false	)
			: T_(T), S0_(S0), sigmaATM_(sigmaATM), Sp_(Sp), Sm_(Sm), Mp_(Mp), Mm_(Mm),
			maxCalibrationIters_(maxCalibrationIters), adjustATM_(adjustATMFlag), enableLogging_(enableLogging) {
			// some basic sanity checks come here to avoid the need for taking care of it later on
			QL_REQUIRE(T_ > 0, "T_ > 0 required.");
			QL_REQUIRE(sigmaATM_ > 0, "sigmaATM_ > 0 required.");
			QL_REQUIRE(Sp_.size() > 0, "Sp_.size() > 0 required.");
			QL_REQUIRE(Sm_.size() > 0, "Sm_.size() > 0 required.");
			QL_REQUIRE(Mp_.size() == Sp_.size(), "Mp_.size() == Sp_.size() required.");
			QL_REQUIRE(Mm_.size() == Sm_.size(), "Mm_.size() == Sm_.size() required.");
			// check for monotonicity
			QL_REQUIRE(Sp_[0] > S0_, "Sp_[0] > S0_ required.");
			for (size_t k=1; k<Sp_.size(); ++k) QL_REQUIRE(Sp_[k] > Sp_[k-1], "Sp_[k] > Sp_[k-1] required.");
			QL_REQUIRE(Sm_[0] < S0_, "Sm_[0] < S0_ required.");
			for (size_t k = 1; k<Sm_.size(); ++k) QL_REQUIRE(Sm_[k] < Sm_[k-1], "Sm_[k] < Sm_[k-1] required.");
			// now it makes sense to allocate memory
			sigmaP_.resize(Sp_.size());
			sigmaM_.resize(Sm_.size());
			Xp_.resize(Sp_.size());
			Xm_.resize(Sm_.size());
			// initialize deep-in-the-model parameters
			straddleATM_ = sigmaATM_ * sqrt(T_) * M_1_SQRTPI * M_SQRT_2 * 2.0;
			sigma0_ = sigmaATM_;
			mu_     = -(Mm_[0] + Mp_[0]) / 4.0 * T_; // this should be exact for shifted log-normal models
			alpha_  = 1.0;
			nu_ =     0.0;
			extrapolationStdevs_ = 10.0;
			sigma0Tol_           = 1.0e-12;
			S0Tol_               = 1.0e-12;
		    // now we may calculate local volatility
			updateLocalVol();
			calibrateATM();
			if (adjustATM_) adjustATM();
		}

		// inspectors

		const std::vector<std::string> logging() { return logging_; }

		const std::vector<PassiveType> underlyingX() {
			std::vector<PassiveType> X(Xm_.size() + Xp_.size() + 1);
			for (size_t k = 0; k < Xm_.size(); ++k) X[k] = Xm_[Xm_.size() - k - 1];
			X[Xm_.size()] = 0.0;
			for (size_t k = 0; k < Xp_.size(); ++k) X[Xm_.size() + 1 + k] = Xp_[k];
			return X;
		}

		const std::vector<PassiveType> underlyingS() {
			std::vector<PassiveType> S(Sm_.size() + Sp_.size() + 1);
			for (size_t k = 0; k < Sm_.size(); ++k) S[k] = Sm_[Sm_.size() - k - 1];
			S[Sm_.size()] = S0_;
			for (size_t k = 0; k < Sp_.size(); ++k) S[Sm_.size() + 1 + k] = Sp_[k];
			return S;
		}

		const std::vector<PassiveType> localVol() {
			std::vector<PassiveType> sigma(sigmaM_.size() + sigmaP_.size() + 1);
			for (size_t k = 0; k < sigmaM_.size(); ++k) sigma[k] = sigmaM_[sigmaM_.size() - k - 1];
			sigma[sigmaM_.size()] = sigma0_;
			for (size_t k = 0; k < sigmaP_.size(); ++k) sigma[sigmaM_.size() + 1 + k] = sigmaP_[k];
			return sigma;
		}

		const std::vector<PassiveType> localVolSlope() {
			std::vector<PassiveType> m(Mm_.size() + Mp_.size() + 1);
			for (size_t k = 0; k < Mm_.size(); ++k) m[k] = Mm_[Mm_.size() - k - 1];
			m[Mm_.size()] = 0.0;  // undefined
			for (size_t k = 0; k < Mp_.size(); ++k) m[Mm_.size() + 1 + k] = Mp_[k];
			return m;
		}

		// model function evaluations

		const PassiveType localVol(PassiveType S) {
			bool isRightWing = (S >= S0_) ? true : false;
			size_t idx = 0;
			if (isRightWing) while ((idx < Sp_.size() - 1) && (Sp_[idx] < S)) ++idx;
			else             while ((idx < Sm_.size() - 1) && (Sm_[idx] > S)) ++idx;
			return localVol(isRightWing, idx, S);
		}

		const PassiveType underlyingS(PassiveType x) {
			bool isRightWing = (x >= 0.0) ? true : false;
			size_t idx = 0;
			if (isRightWing) while ((idx < Xp_.size() - 1) && (Xp_[idx] < x)) ++idx;
			else             while ((idx < Xm_.size() - 1) && (Xm_[idx] > x)) ++idx;
			return underlyingS(isRightWing, idx, x);
		}

		// calculating expectations

		const PassiveType expectation(bool isRightWing, PassiveType strike) {
			// calculate the forward price of an OTM option
			size_t idx = 0;
			if (isRightWing) {
				QL_REQUIRE(strike >= S0_, "strike >= S0_ required");
				while ((idx < Sp_.size()) && (Sp_[idx] <= strike)) ++idx;  // make sure strike < Sp_[idx]
				if (idx == Sp_.size()) return 0.0;  // we are beyond exrapolation
				PassiveType strikeX = underlyingX(isRightWing, idx, strike);
				PassiveType x0 = (idx > 0) ? Xp_[idx - 1] : 0.0;
				QL_REQUIRE((x0 <= strikeX) && (strikeX < Xp_[idx]), "(x0 <= strikeX) && (strikeX < Xp_[idx]) required");
				PassiveType intS = 0.0;
				for (size_t k = idx; k < Sp_.size(); ++k) {
					PassiveType xStart = (k == idx) ? strikeX : Xp_[k - 1];
					intS += (primitiveF(isRightWing, k, Xp_[k]) - primitiveF(isRightWing, k, xStart));
				}
				// we need to adjust for the strike integral
				PassiveType xEnd = Xp_.back();
				PassiveType intK = TemplateAuxilliaries::Phi((xEnd - mu_) / sqrt(T_)) - TemplateAuxilliaries::Phi((strikeX - mu_) / sqrt(T_));
				return intS - strike * intK;
			}
			else {
				QL_REQUIRE(strike <= S0_, "strike <= S0_ required");
				while ((idx < Sm_.size()) && (Sm_[idx] >= strike)) ++idx;  // make sure Sm_[idx] < strke
				if (idx == Sm_.size()) return 0.0;  // we are beyond exrapolation
				PassiveType strikeX = underlyingX(isRightWing, idx, strike);
				PassiveType x0 = (idx > 0) ? Xm_[idx - 1] : 0.0;
				QL_REQUIRE((x0 >= strikeX) && (strikeX > Xm_[idx]), "(x0 >= strikeX) && (strikeX > Xm_[idx]) required");
				PassiveType intS = 0.0;
				for (size_t k = idx; k < Sm_.size(); ++k) {
					PassiveType xStart = (k == idx) ? strikeX : Xm_[k - 1];
					intS += (primitiveF(isRightWing, k, Xm_[k]) - primitiveF(isRightWing, k, xStart));
				}
				// we need to adjust for the strike integral
				PassiveType xEnd = Xm_.back();
				PassiveType intK = TemplateAuxilliaries::Phi((xEnd - mu_) / sqrt(T_)) - TemplateAuxilliaries::Phi((strikeX - mu_) / sqrt(T_));
				return intS - strike * intK;
			}
		}


	};


}


#endif  /* ifndef quantlib_templatevanillalocalvolmodel_hpp_hpp */
