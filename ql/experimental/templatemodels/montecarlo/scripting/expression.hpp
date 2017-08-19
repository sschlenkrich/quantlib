/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2017, Sebastian Schlenkrich

*/

/*! \file expression.hpp
\brief an abstract expression object for payoff representation

*/



#ifndef quantlib_scripting_expression_hpp
#define quantlib_scripting_expression_hpp

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

namespace QuantLib {

	namespace Scripting {

		class Expression {
		public:
			enum Type {
				NEXT,
				ASSIGNMENT,
				UNARYPLUS,
				UNARYMINUS,
				PLUS,
				MINUS,
				MULT,
				DIVISION,
				IDENTIFIER,
				NUMBER,
				IFTHENELSE,
				MIN,
				MAX,
				PAY,
				CACHE,
				UNKNOWNTYPE
			};
			Expression(const Type                          type,
				const std::string                   leaf,
				const boost::shared_ptr<Expression> first = 0,
				const boost::shared_ptr<Expression> second = 0,
				const boost::shared_ptr<Expression> third = 0);

			std::string toString();

			// inspectors
			inline Type type()         { return type_;  }
			inline std::string leaf()  { return leaf_;  }
			inline const std::vector<boost::shared_ptr<Expression>>& childs() { return childs_; };
		private:
			Type type_;
			std::string leaf_;
			std::vector<boost::shared_ptr<Expression>> childs_;

		};
	}
}

#endif // !EXPRESSION_HPP
