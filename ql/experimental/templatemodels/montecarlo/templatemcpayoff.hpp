/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file templatemcpayoff.hpp
    \brief generic payoff interface for MC simulation
	
*/


#ifndef quantlib_templatemcpayoff_hpp
#define quantlib_templatemcpayoff_hpp


#include <ql/experimental/template/montecarlo/templatemcsimulation.hpp>



namespace QuantLib {

	// Base class for template payoffs
    template <class DateType, class PassiveType, class ActiveType>
	class TemplateMCPayoff {
	protected:
		typedef TemplateMCSimulation<DateType, PassiveType, ActiveType>        SimulationType;
		typedef typename TemplateMCSimulation<DateType, PassiveType, ActiveType>::Path  PathType;

	    DateType observationTime_;
	public:
	    TemplateMCPayoff( const DateType observationTime ) : observationTime_(observationTime) { }
		// inspectors
		inline DateType observationTime() { return observationTime_; }
		// generic payoff(observationTime, p) needs to be implemented by derived classes
        inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) = 0;
		// discounted payoff for NPV valuation
		inline virtual ActiveType discountedAt(const boost::shared_ptr<PathType>& p) { return at(p) / p->numeraire(observationTime_); }

		// generic pricer
		class Pricer {
		protected:
			std::vector< boost::shared_ptr<TemplateMCPayoff> > payoffs_;
			boost::shared_ptr<SimulationType>                  simulation_;
		public:

			inline static std::vector<ActiveType> at( const boost::shared_ptr<TemplateMCPayoff>& payoff,
				                                      const boost::shared_ptr<SimulationType>&   simulation) {
			    std::vector<ActiveType> res(simulation->nPaths());
				for (size_t k=0; k<simulation->nPaths(); ++k) res[k] = payoff->at(simulation->path(k));
				return res;
			}

			inline static std::vector<ActiveType> discountedAt( const boost::shared_ptr<TemplateMCPayoff>& payoff,
				                                                const boost::shared_ptr<SimulationType>&   simulation) {
			    std::vector<ActiveType> res(simulation->nPaths());
				for (size_t k=0; k<simulation->nPaths(); ++k) res[k] = payoff->discountedAt(simulation->path(k));
				return res;
			}

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
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				//if (payTime_<=observationTime()) return (ActiveType)1.0;
				return p->zeroBond(observationTime(),payTime_);
			}
		};

		// call or put exercised at observation time and settled at pay time
		class VanillaOption : public TemplateMCPayoff {
		protected:
			DateType    payTime_;
			PassiveType callOrPut_;
			PassiveType strike_;
		public:
			VanillaOption( DateType obsTime, DateType payTime, PassiveType strike, PassiveType callOrPut ) : TemplateMCPayoff(obsTime), payTime_(payTime), strike_(strike), callOrPut_(callOrPut) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if (payTime_<observationTime()) return (ActiveType)0.0;
				ActiveType DF = p->zeroBond(observationTime(),payTime_);
				ActiveType S  = p->asset(observationTime());
				ActiveType V  = callOrPut_ * DF * (S - strike_);
				return (V>0.0) ? (V) : ((ActiveType)0.0);
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
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
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
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
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

		// prototypical physically settled European swaption
		class GeneralSwaption : public TemplateMCPayoff {
		protected:
			std::vector<DateType>    floatTimes_;     // T_1, .., T_M
			std::vector<PassiveType> floatWeights_;   // u_1, .., u_M
			std::vector<DateType>    fixedTimes_;     // T_1, .., T_N
			std::vector<PassiveType> fixedWeights_;   // w_1, .., w_N
			PassiveType              strikeRate_;   // option strike
			PassiveType              payOrRec_;     // call (+1) or put (-1) option on swap rate
		public:
			GeneralSwaption( DateType                        obsTime, 
				             const std::vector<DateType>&    floatTimes,
				             const std::vector<PassiveType>& floatWeights,
				             const std::vector<DateType>&    fixedTimes,
				             const std::vector<PassiveType>& fixedWeights,
					         PassiveType                     strikeRate,
					         PassiveType                     payOrRec      )
				: TemplateMCPayoff(obsTime),  floatTimes_(floatTimes), floatWeights_(floatWeights),
				  fixedTimes_(fixedTimes), fixedWeights_(fixedWeights), strikeRate_(strikeRate), payOrRec_(payOrRec) { 
			    // check consistency of swap
			    // float leg
			    QL_REQUIRE(floatWeights.size()>0,"GeneralSwaption: empty float weights.");
			    QL_REQUIRE(floatTimes.size()==floatWeights.size(),"GeneralSwaption: float sizes mismatch.");
			    QL_REQUIRE(floatTimes[0]>0,"GeneralSwaption: future float times required");
			    for (size_t k=1; k<floatTimes.size(); ++k) QL_REQUIRE(floatTimes[k]>=floatTimes[k-1],"GeneralSwaption: ascending float times required");
			    // fixed leg
			    QL_REQUIRE(fixedWeights.size()>0,"GeneralSwaption: empty fixed weights.");
			    QL_REQUIRE(fixedTimes.size()==fixedWeights.size(),"GeneralSwaption: fixed sizes mismatch.");
			    QL_REQUIRE(fixedTimes[0]>0,"GeneralSwaption: future fixed times required");
			    for (size_t k=1; k<fixedTimes.size(); ++k) QL_REQUIRE(fixedTimes[k]>=fixedTimes[k-1],"GeneralSwaption: ascending fixed times required");
				// finished
			}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				ActiveType floatleg = 0.0;
				ActiveType annuity  = 0.0;
				// float leg
				for (size_t k=0; k<floatTimes_.size(); ++k) floatleg += floatWeights_[k] * p->zeroBond(observationTime(),floatTimes_[k]);
				// annuity
				for (size_t k=0; k<fixedTimes_.size(); ++k) annuity  += fixedWeights_[k] * p->zeroBond(observationTime(),fixedTimes_[k]);
				// floatleg - fixedleg...
				ActiveType res = floatleg - strikeRate_*annuity;
				// payer or receiver swap...
				res *= payOrRec_;
				// exercise option...
				res = (res>0) ? res : 0.0;
				return res;
			}
		};

		// future swap rate
		class SwapRate : public GeneralSwaption {
		public:
			SwapRate( DateType                        obsTime, 
				      const std::vector<DateType>&    floatTimes,
				      const std::vector<PassiveType>& floatWeights,
				      const std::vector<DateType>&    fixedTimes,
				      const std::vector<PassiveType>& annuityWeights )
					  : GeneralSwaption(obsTime,floatTimes,floatWeights,fixedTimes,annuityWeights,0.0,1) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				ActiveType floatleg = 0.0;
				ActiveType annuity  = 0.0;
				// float leg
				for (size_t k=0; k<floatTimes_.size(); ++k) floatleg += floatWeights_[k] * p->zeroBond(observationTime(),floatTimes_[k]);
				// annuity
				for (size_t k=0; k<fixedTimes_.size(); ++k) annuity  += fixedWeights_[k] * p->zeroBond(observationTime(),fixedTimes_[k]);
				return floatleg / annuity;
			}
		    // payoff should NOT be discounted
		    inline virtual ActiveType discountedAt(const boost::shared_ptr<PathType>& p) { return at(p); }
		};


		// undiscounted correlation between prototypical physically settled European swaption
		class ModelCorrelation : public TemplateMCPayoff {
		protected:
			std::vector<DateType> times_;
			DateType T1_, T2_;
			ActiveType swapRate(const boost::shared_ptr<PathType>& p, const DateType t, const DateType TN ) {
				ActiveType num = p->zeroBond(t,t) - p->zeroBond(t,TN);
				ActiveType den = 0.0;
				for (ActiveType Ti = t; Ti<TN; Ti+=1.0) {
					ActiveType T = (Ti+1.0>TN) ? (TN) : (Ti+1.0);
					den += (T-Ti) * p->zeroBond(t,T);
				}
				return num / den;
			}
		public:
			ModelCorrelation( const std::vector<DateType>&    times,   // observation times
				              const DateType                  T1,      // swap term one
							  const DateType                  T2 )     // swap term two
							  : TemplateMCPayoff(0.0), times_(times), T1_(T1), T2_(T2) {
				QL_REQUIRE(times_.size()>1,"ModelCorrelation: At least two observation times required.");
			}
		    // payoff should NOT be discounted
		    inline virtual ActiveType discountedAt(const boost::shared_ptr<PathType>& p) { return at(p); }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				std::vector<ActiveType> dS1(times_.size()-1), dS2(times_.size()-1);
				ActiveType EdS1 = 0.0, EdS2 = 0.0; 
				for (size_t i=1; i<times_.size(); ++i) {
					dS1[i-1] =  swapRate(p,times_[i],times_[i]+T1_) - swapRate(p,times_[i-1],times_[i-1]+T1_);
					dS2[i-1] =  swapRate(p,times_[i],times_[i]+T2_) - swapRate(p,times_[i-1],times_[i-1]+T2_);
					EdS1 += dS1[i-1];
					EdS2 += dS2[i-1];
				}
				EdS1 /= dS1.size();
				EdS2 /= dS2.size();
				ActiveType Var1=0.0, Var2=0.0, Cov=0.0;
				for (size_t i=0; i<times_.size()-1; ++i) {
					Var1 += (dS1[i] - EdS1)*(dS1[i] - EdS1);
					Var2 += (dS2[i] - EdS2)*(dS2[i] - EdS2);
					Cov  += (dS1[i] - EdS1)*(dS2[i] - EdS2);
				}
				return Cov / sqrt(Var1*Var2);
			}
		};

		// undiscounted correlation between forward rates
		class ForwardRateCorrelation : public TemplateMCPayoff {
		protected:
			std::vector<DateType> times_;
			DateType T1_, Term1_; 
			DateType T2_, Term2_;
			ActiveType fwSwapRate(const boost::shared_ptr<PathType>& p, const DateType t, const DateType TSettle, const DateType Term ) {
				ActiveType num = p->zeroBond(t,TSettle) - p->zeroBond(t,TSettle+Term);
				ActiveType den = 0.0;
				for (ActiveType Ti = TSettle; Ti<TSettle+Term; Ti+=1.0) {
					ActiveType T = (Ti+1.0>TSettle+Term) ? (TSettle+Term) : (Ti+1.0);
					den += (T-Ti) * p->zeroBond(t,T);
				}
				return num / den;
			}
			ActiveType fraRate(const boost::shared_ptr<PathType>& p, const DateType t, const DateType TSettle, const DateType Term ) {
				ActiveType rate = (p->zeroBond(t,TSettle) / p->zeroBond(t,TSettle+Term) - 1.0)/Term;
				return rate;
			}
		public:
			ForwardRateCorrelation( const std::vector<DateType>&    times,   // observation times
				                    const DateType                  T1,      // fixing date one
				                    const DateType                  Term1,   // tenor one
							        const DateType                  T2,      // fixing date two
				                    const DateType                  Term2)   // tenor two
							        : TemplateMCPayoff(0.0), times_(times), T1_(T1), Term1_(Term1), T2_(T2), Term2_(Term2) {
				QL_REQUIRE(times_.size()>1,"ModelCorrelation: At least two observation times required.");
			}
		    // payoff should NOT be discounted
		    inline virtual ActiveType discountedAt(const boost::shared_ptr<PathType>& p) { return at(p); }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				std::vector<ActiveType> dS1(times_.size()-1), dS2(times_.size()-1);
				ActiveType EdS1 = 0.0, EdS2 = 0.0; 
				for (size_t i=1; i<times_.size(); ++i) {
					dS1[i-1] =  fraRate(p,times_[i],T1_,Term1_) - fraRate(p,times_[i-1],T1_,Term1_);
					dS2[i-1] =  fraRate(p,times_[i],T2_,Term2_) - fraRate(p,times_[i-1],T2_,Term2_);
					EdS1 += dS1[i-1];
					EdS2 += dS2[i-1];
				}
				EdS1 /= dS1.size();
				EdS2 /= dS2.size();
				ActiveType Var1=0.0, Var2=0.0, Cov=0.0;
				for (size_t i=0; i<times_.size()-1; ++i) {
					Var1 += (dS1[i] - EdS1)*(dS1[i] - EdS1);
					Var2 += (dS2[i] - EdS2)*(dS2[i] - EdS2);
					Cov  += (dS1[i] - EdS1)*(dS2[i] - EdS2);
				}
				return Cov / sqrt(Var1*Var2);
			}
		};


		// Commodity model payoffs

		// call or put on a strip of futures
		class AverageFutureOption : public TemplateMCPayoff {
			std::vector<DateType>    settlementTimes_;
			std::vector<PassiveType> settlementWeights_;
			PassiveType              strike_;        // option strike
			PassiveType              callOrPut_;     // call (+1) or put (-1) option on swap rate
		public:
			AverageFutureOption( DateType                        obsTime, 
				                 const std::vector<DateType>&    settlementTimes, 
				                 const std::vector<PassiveType>& settlementWeights,
					             const PassiveType               strike,
					             const PassiveType               callOrPut )
				: TemplateMCPayoff(obsTime), settlementTimes_(settlementTimes), settlementWeights_(settlementWeights), strike_(strike), callOrPut_(callOrPut) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				ActiveType fut=0.0;
				for (size_t k=0; k<settlementTimes_.size(); ++k)
					if (settlementTimes_[k]>=observationTime()) fut += settlementWeights_[k] * p->future(observationTime(),settlementTimes_[k]);
				ActiveType V  = callOrPut_ * (fut - strike_);
				return (V>0.0) ? (V) : ((ActiveType)0.0);
			}
		};

		// covariance between quoted futures with roll-over
		class AverageFutureCovariance : public TemplateMCPayoff {
		protected:
			std::vector<DateType> obsTimes_;
			// future A definition
			std::vector<DateType>    settlementTimesA_;     // delivery period
			std::vector<PassiveType> settlementWeightsA_;   // usually equal weights for averaging
			PassiveType              obsLagA_;              // roll over delivery period if first date is less than observation lag
			// future B definition
			std::vector<DateType>    settlementTimesB_;     // delivery period
			std::vector<PassiveType> settlementWeightsB_;   // usually equal weights for averaging
			PassiveType              obsLagB_;              // roll over delivery period if first date is less than observation lag
			// further flags
			bool                     useLogReturns_;        // specify (normal vs lognormal) type of returns considered
			long                     calcType_;             // flag to distinguish what to do (quick and dirty)
			                                                // 0 - calculate covariance
			                                                // 1 - calculate correlation instead of covariance
			                                                // 2 - calculate spread variance 
			// helper methods
			void rollOver( std::vector<DateType>& settlementTimes ) {
			    PassiveType rollPeriod = settlementTimes[settlementTimes.size()-1] - settlementTimes[0] + 1.0/365.0;
			    for (size_t k=0; k<settlementTimes.size(); ++k) settlementTimes[k] += rollPeriod;
			}
			ActiveType averageFuture(const boost::shared_ptr<PathType>& p, const DateType t, const std::vector<DateType>& settlementTimes, const std::vector<PassiveType>& settlementWeights ) {
				ActiveType fut=0.0;
				for (size_t k=0; k<settlementTimes.size(); ++k) if (settlementTimes[k]>=t) fut += settlementWeights[k] * p->future(t,settlementTimes[k]);
				return fut;
			}
		public:
			AverageFutureCovariance( const std::vector<DateType>&    obsTimes,
				                     const std::vector<DateType>&    settlementTimesA,     
			                         const std::vector<PassiveType>& settlementWeightsA,
			                         const PassiveType               obsLagA,
			                         const std::vector<DateType>&    settlementTimesB,
			                         const std::vector<PassiveType>& settlementWeightsB,
			                         const PassiveType               obsLagB,
			                         const bool                      useLogReturns,       
			                         const long                      calcType )
			: TemplateMCPayoff(0.0), obsTimes_(obsTimes), settlementTimesA_(settlementTimesA), settlementWeightsA_(settlementWeightsA), obsLagA_(obsLagA),
			settlementTimesB_(settlementTimesB), settlementWeightsB_(settlementWeightsB), obsLagB_(obsLagB), useLogReturns_(useLogReturns), calcType_(calcType) {}
		    // payoff should NOT be discounted
		    inline virtual ActiveType discountedAt(const boost::shared_ptr<PathType>& p) { return at(p); }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				// if only one observation time then return [FutA * FutB] for autocorrelation estimation
				if (obsTimes_.size()==1) return averageFuture(p,settlementTimesA_[0],settlementTimesA_,settlementWeightsA_) * averageFuture(p,settlementTimesB_[0],settlementTimesB_,settlementWeightsB_);
				std::vector<ActiveType> dFutA(obsTimes_.size()-1), dFutB(obsTimes_.size()-1), dSprd(obsTimes_.size()-1);
				ActiveType EdFutA = 0.0, EdFutB = 0.0, EdSprd = 0.0; 
				for (size_t i=1; i<obsTimes_.size(); ++i) {
					// check rollover A and B
					if (settlementTimesA_[0]-obsTimes_[i]<obsLagA_) rollOver(settlementTimesA_);
					if (settlementTimesB_[0]-obsTimes_[i]<obsLagB_) rollOver(settlementTimesB_);
					// calculate returns
					if (useLogReturns_) {
					    dFutA[i-1] = log(averageFuture(p,obsTimes_[i],settlementTimesA_,settlementWeightsA_)) - log(averageFuture(p,obsTimes_[i-1],settlementTimesA_,settlementWeightsA_));
					    dFutB[i-1] = log(averageFuture(p,obsTimes_[i],settlementTimesB_,settlementWeightsB_)) - log(averageFuture(p,obsTimes_[i-1],settlementTimesB_,settlementWeightsB_));
					} else {
					    dFutA[i-1] = averageFuture(p,obsTimes_[i],settlementTimesA_,settlementWeightsA_) - averageFuture(p,obsTimes_[i-1],settlementTimesA_,settlementWeightsA_);
					    dFutB[i-1] = averageFuture(p,obsTimes_[i],settlementTimesB_,settlementWeightsB_) - averageFuture(p,obsTimes_[i-1],settlementTimesB_,settlementWeightsB_);
						if (calcType_==2) dSprd[i-1] = dFutB[i-1] - dFutA[i-1];
					}
					EdFutA += dFutA[i-1];
					EdFutB += dFutB[i-1];
				}
				EdFutA /= dFutA.size();
				EdFutB /= dFutB.size();
				if (calcType_==2) EdSprd = EdFutB - EdFutA;
				ActiveType VarA=0.0, VarB=0.0, Cov=0.0, SprdVar=0.0;
				for (size_t i=0; i<obsTimes_.size()-1; ++i) {
					if (calcType_<2) {
					    VarA += (dFutA[i] - EdFutA)*(dFutA[i] - EdFutA);
					    VarB += (dFutB[i] - EdFutB)*(dFutB[i] - EdFutB);
					    Cov  += (dFutA[i] - EdFutA)*(dFutB[i] - EdFutB);
					} else {
						SprdVar += (dSprd[i]-EdSprd)*(dSprd[i]-EdSprd);
					}
				}
				if (calcType_==0) return Cov     / (dFutA.size()-1);    // covariance
				if (calcType_==1) return Cov     / sqrt(VarA*VarB);     // correlation
				if (calcType_==2) return SprdVar / (dSprd.size()-1);  // correlation
				return 0.0; // fall back
			}

		};

	};



}

#endif  /* ifndef quantlib_templatemcpayoff_hpp */
