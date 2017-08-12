
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2017, Sebastian Schlenkrich

*/

/*! \file expression.cpp
\brief an abstract expression object for payoff representation

*/



#include <ql/experimental/templatemodels/montecarlo/scripting/expression.hpp>


namespace QuantLib {

	namespace Scripting {

		Expression::Expression(
			const Type                          type,
			const std::string                   leaf,
			const boost::shared_ptr<Expression> first,
			const boost::shared_ptr<Expression> second,
			const boost::shared_ptr<Expression> third)
			: type_(type), leaf_(leaf) {
			if (first) childs_.push_back(first);
			if (second) childs_.push_back(second);
			if (third) childs_.push_back(third);
		}

		std::string Expression::toString() {
			if (type_ == Type::NUMBER)     return std::string("NUM(" + leaf_ + ")\n");
			if (type_ == Type::IDENTIFIER) return std::string("TXT(" + leaf_ + ")\n");
			if (type_ == Expression::ASSIGNMENT) return leaf_ + " = " + childs_[0]->toString();
			std::string res;
			switch (type_) {
			case Expression::NEXT:
				res = "NEXT";
				break;
			case Expression::PLUS:
				res = "+";
				break;
			case Expression::MINUS:
				res = "-";
				break;
			case Expression::MULT:
				res = "*";
				break;
			case Expression::DIVISION:
				res = "/";
				break;
			case Expression::IDENTIFIER:
				res = leaf_;
				break;
			default:
				break;
			}
			res = res + "[ \n";
			for (size_t k = 0; k < childs_.size(); ++k) res = res + childs_[k]->toString();
			res = res + "] \n";
			return res;
		}

	}
}