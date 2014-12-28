/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010, Sebastian Schlenkrich

*/

/*! \file templateintegrators.hpp
    \brief provide template functions for numerical integration
*/


#ifndef quantlib_templateintegrators_hpp
#define quantlib_templateintegrators_hpp

//#include <boost/math/special_functions/erf.hpp>
//#include <ql/experimental/template/auxilliaries/MinimADVariable2.hpp>


namespace TemplateAuxilliaries {

	// evaluate \int_a^b v(t) f(t) dt = \sum v_i [F(t_i) - F(t_i-1)] with
	// v(t) piece-wise left-constant,
	// F'(t) = f(t)
   	template <typename PassiveType, typename ActiveType, typename FuncType>
	class PieceWiseConstantIntegral {
	private:
		std::vector<PassiveType> t_;
		std::vector<ActiveType> v_;
		FuncType F_;
	public:
		PieceWiseConstantIntegral(const std::vector<PassiveType>& t, const std::vector<ActiveType>& v, const FuncType& F) : t_(t), v_(v), F_(F) {}
		ActiveType operator()(PassiveType startTime, PassiveType endTime) {
			int sgn = 1;
			if (startTime>endTime) {  // we want to ensure startTime <= endTime
				PassiveType t = startTime;
				startTime = endTime;
				endTime = t;
				sgn = -1;
			}
			// organising indices
			size_t idx_min  = 0;
			size_t idx_max  = std::min(t_.size(),v_.size())-1;
			size_t idx_last = idx_max;
			// enforce a < t_min <= t_max < b or special treatment
			while ((startTime>=t_[idx_min])&&(idx_min<idx_last)) ++idx_min;
			while ((endTime  <=t_[idx_max])&&(idx_max>0       )) --idx_max;
			ActiveType tmp = sgn * ( F_(endTime) - F_(startTime) );
			if (endTime<=t_[0])    return v_[0]        * tmp;  // short end
			if (idx_min==idx_last) return v_[idx_last] * tmp;  // long end
			if (idx_min> idx_max)  return v_[idx_min]  * tmp;  // integration within grid intervall
            // integral a ... x_min
			tmp = v_[idx_min] * ( F_(t_[idx_min]) - F_(startTime) );
		    // integral x_min ... x_max
			for (size_t i=idx_min; i<idx_max; ++i) tmp += v_[i+1] * ( F_(t_[i+1]) - F_(t_[i]) );
            // integral x_max ... b
			if (idx_max<idx_last) tmp += v_[idx_max+1] * ( F_(endTime) - F_(t_[idx_max]) );
			else                  tmp += v_[idx_max]   * ( F_(endTime) - F_(t_[idx_max]) );
            // finished
			return sgn * tmp;
		}



	};


}

#endif  /* quantlib_templateintegrators_hpp */
