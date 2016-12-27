/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016, Sebastian Schlenkrich

*/

/*! \file mcscriptT.hpp
    \brief script payoffs for MC simulation
	
*/


#ifndef quantlib_templatemcscript_hpp
#define quantlib_templatemcscript_hpp


#include <ql/experimental/templatemodels/montecarlo/mcpayoffT.hpp>
#include <boost/regex.hpp>



namespace QuantLib {

	template <class DateType, class PassiveType, class ActiveType>
	class MCScriptT : public MCPayoffT<DateType, PassiveType, ActiveType> {
	private:
		std::map<std::string, boost::shared_ptr<MCPayoffT>>  payoffs_;    // the actual payoffs which may be accessed
		std::vector<std::string>                             scriptLog_;  // log messages when parsing the script
		boost::shared_ptr<MCPayoffT>                         result_;     // result payoff for MCPayoff interface implementation
	public:
		MCScriptT(const std::vector<std::string>&                    keys,
			      const std::vector <boost::shared_ptr<MCPayoffT>>&  payoffs,
			      const std::vector<std::string>&                    script,
			      const bool                                         overwrite=true) : MCPayoffT<DateType, PassiveType, ActiveType>(0.0) {
			QL_REQUIRE(keys.size()==payoffs.size(), "MCScript error: key vs. value size missmatch");
			for (Size k = 0; k < keys.size(); ++k) { // initialize map
				std::map<std::string, boost::shared_ptr<MCPayoffT>>::iterator it = payoffs_.find(keys[k]);
				if (it == payoffs_.end()) { // insert a new element
					payoffs_.insert( std::make_pair(keys[k], payoffs[k]) );
				} else { // potentially overwrite existing element
					QL_REQUIRE(overwrite, "MCScript error: overwrite not allowed");
					it->second = payoffs[k];
				}
			}
			parseScript(script, overwrite); // for briefty we delegate parsing to separate method
			// we need to find a 'result' payoff
			std::map<std::string, boost::shared_ptr<MCPayoffT>>::iterator it = payoffs_.find("result");
			QL_REQUIRE(it != payoffs_.end(), "MCScript error: 'result' payoff element not found.");
			result_ = it->second;
			observationTime_ = result_->observationTime();
		}

		inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
			return result_->at(p);
		}

		// inspector
		inline const std::map<std::string, boost::shared_ptr<MCPayoffT>>& payoffs()   { return payoffs_;   }
		inline const std::vector<std::string>&                            scriptLog() { return scriptLog_; } 

		// MC valuation
		inline std::vector<ActiveType> NPV(const boost::shared_ptr<SimulationType>&  simulation,
			                               const std::vector<std::string>&           keys) {
			std::vector<boost::shared_ptr<MCPayoffT>> payoffs(keys.size());
			for (Size k = 0; k < keys.size(); ++k) {
				std::map<std::string, boost::shared_ptr<MCPayoffT>>::iterator it = payoffs_.find(keys[k]);
				QL_REQUIRE(it != payoffs_.end(), "MCScript error: payoff '" + keys[k] + "' not found");
				payoffs[k] = it->second;
			}
			std::vector<ActiveType> npv(payoffs.size(), (ActiveType)0.0);
			for (Size n = 0; n < simulation->nPaths(); ++n) {
				const boost::shared_ptr<MCPayoffT::PathType> p(simulation->path(n));
				for (Size k = 0; k < payoffs.size(); ++k) {
					npv[k] += payoffs[k]->discountedAt(p);
				}
			}
			for (Size k = 0; k < payoffs.size(); ++k) npv[k] /= simulation->nPaths();
			return npv;
		}

	private:

		/*

		we implement the following non-recursive grammar
		
		line  =  var '=' expr
		var   =  [a-zA-Z][a-zA-Z0-9]*           { RegEx }
		expr  =  operator | function | payoff   { apply from left to right }

		operator   =  operator1 | operator2 
		operator1  =  ['+' | '-'] payoff
		operator2  =  payoff ['+' | '-' | '*' | == | != | < | <= | > | >= | && | || ] payoff

		function   =  function3 | function2 | function1
		function3  =  fname3 '(' payoff ',' payoff ',' payoff ')'
		function2  =  fname2 '(' payoff ',' payoff ')'
		function1  =  fname1 '(' payoff ')'

		fname3     =  'IfThenElse'
		fname2     =  'Min' | 'Max | Pay'
		fname1     =  'Cache'

		payoff  =  number | string              { try double conversion and lookup in map }

		*/

		// parse the script and set up payoffs
		inline void parseScript(const std::vector<std::string>&  script,
			const bool                       overwrite = true) {
			for (Size k = 0; k < script.size(); ++k) {  // parse lines
				std::string line = boost::regex_replace(script[k], boost::regex(" "), ""); // remove whitespaces
				boost::smatch what;

				bool isAssignment = boost::regex_match(line, what, boost::regex("([a-zA-Z][a-zA-Z0-9]*)(=)(.+)"));
				if (!isAssignment) {
					scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": '" + line + "' is no valid assignment"));
					continue; // move to next item in script
				}
				std::string var(what[1]), expr(what[3]);

				boost::shared_ptr<MCPayoffT> p;

				if (boost::regex_match(expr, what, boost::regex("(\\+|-)(.+)"))) {
					p = operator1(std::string(what[1]), std::string(what[2]), k);
				}
				else if (boost::regex_match(expr, what, boost::regex("(.+)(\\+|-|\\*|==|!=|<=|<|>=|>|&&|\\|\\|)(.+)"))) {
					p = operator2(std::string(what[2]), std::string(what[1]), std::string(what[3]), k);
				}
				else if (boost::regex_match(expr, what, boost::regex("([a-zA-Z]+)\\((.+),(.+),(.+)\\)"))) {
					p = function3(std::string(what[1]), std::string(what[2]), std::string(what[3]), std::string(what[4]), k);
				}
				else if (boost::regex_match(expr, what, boost::regex("([a-zA-Z]+)\\((.+),(.+)\\)"))) {
					p = function2(std::string(what[1]), std::string(what[2]), std::string(what[3]), k);
				}
				else if (boost::regex_match(expr, what, boost::regex("([a-zA-Z]+)\\((.+)\\)"))) {
					p = function1(std::string(what[1]), std::string(what[2]), k);
				}
				else p = payoff(expr, k); // action of last resort

				if (!p) {
					scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": '" + expr + "' is no valid expression"));
					continue; // move to next item in script
				}

				// now we have a payoff which we may store in the map
				std::map<std::string, boost::shared_ptr<MCPayoffT>>::iterator it = payoffs_.find(var);
				if (it == payoffs_.end()) { // insert a new element
					payoffs_.insert(std::make_pair(var, p));
					scriptLog_.push_back(std::string("Insert line " + std::to_string(k) + ": '" + line + "'"));
					continue;
				}
				if (overwrite) {
					it->second = p;
					scriptLog_.push_back(std::string("Replace line " + std::to_string(k) + ": '" + line + "'"));
					continue;
				}
				// if we end up here we have a valid payoff but are not allowed to overwrite existing map entry
				scriptLog_.push_back(std::string("Error line " + std::to_string(k) + ": '" + var + "' can not be replaced"));
			}
			if (script.size() == 0) { // in this case the previous loop was not executed and we just print some help details
				scriptLog_.push_back(std::string("we implement the following non-recursive grammar                                     "));
				scriptLog_.push_back(std::string("                                                                                     "));
				scriptLog_.push_back(std::string("line  =  var '=' expr                                                                "));
				scriptLog_.push_back(std::string("var   =  [a-zA-Z][a-zA-Z0-9]*           { RegEx }                                    "));
				scriptLog_.push_back(std::string("expr  =  operator | function | payoff   { apply from left to right }                 "));
				scriptLog_.push_back(std::string("                                                                                     "));
				scriptLog_.push_back(std::string("operator   =  operator1 | operator2                                                  "));
				scriptLog_.push_back(std::string("operator1  =  ['+' | '-'] payoff                                                     "));
				scriptLog_.push_back(std::string("operator2  =  payoff ['+' | '-' | '*' |                                              "));
				scriptLog_.push_back(std::string("                      '==' | '!=' | '<=' |'<' | '>=' | '>' | '&&' | '||' ] payoff    "));
				scriptLog_.push_back(std::string("                                                                                     "));
				scriptLog_.push_back(std::string("function   =  function3 | function2 | function1                                      "));
				scriptLog_.push_back(std::string("function3  =  fname3 '(' payoff ',' payoff ',' payoff ')'                            "));
				scriptLog_.push_back(std::string("function2  =  fname2 '(' payoff ',' payoff ')'                                       "));
				scriptLog_.push_back(std::string("function1  =  fname1 '(' payoff ')'                                                  "));
				scriptLog_.push_back(std::string("                                                                                     "));
				scriptLog_.push_back(std::string("fname3     =  'IfThenElse'                                                           "));
				scriptLog_.push_back(std::string("fname2     =  'Min' | 'Max' | 'Pay'                                                  "));
				scriptLog_.push_back(std::string("fname1     =  'Cache'                                                                "));
				scriptLog_.push_back(std::string("                                                                                     "));
				scriptLog_.push_back(std::string("payoff  =  number | string              { try double conversion and lookup in map }  "));
			}
		}


		// convert string to number
		inline bool to_Number(const std::string str, ActiveType& number) {
			double res;
			std::string::size_type sz;
			try {
				res = std::stod(str,&sz);
			} catch (std::exception e) {
				return false;
			}
			number = res;
			return true;
		}

		// compile fixed cash flow or lookup in map
		inline boost::shared_ptr<MCPayoffT> payoff(const std::string expr, const Size lineNr) {
			ActiveType amount;
			bool isFixed = to_Number(expr,amount);
			if (isFixed) {
				scriptLog_.push_back(std::string("Payoff line " + std::to_string(lineNr) + ": '" + boost::lexical_cast<std::string>(amount) +"' is fixed amount"));
				return boost::shared_ptr<MCPayoffT>(new MCPayoffT::FixedAmount(amount));
			}
			std::map<std::string, boost::shared_ptr<MCPayoffT>>::iterator it = payoffs_.find(expr);
			if (it != payoffs_.end()) {
				scriptLog_.push_back(std::string("Payoff line " + std::to_string(lineNr) + ": '" + expr +"' is in map"));
				return it->second;
			}
			// if we end up here no conversion was successfull
			scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + expr + "' is no payoff"));
			return 0;
		}

		// compile single operand function
		inline boost::shared_ptr<MCPayoffT> function1(const std::string fname, const std::string operand, const Size lineNr) {
			boost::shared_ptr<MCPayoffT> p = payoff(operand, lineNr);
			if (!p) {
				scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + operand + "' is no valid operand"));
				return 0;
			}
    		boost::smatch what;
			if (boost::regex_match(fname, what, boost::regex("Cache")))
				return boost::shared_ptr<MCPayoffT>(new MCPayoffT::Cache(p));
			// if we end up here the function name is not valid
			scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + fname + "' is no valid unary function name"));
			return 0;
		}

		// compile dual operand function
		inline boost::shared_ptr<MCPayoffT> function2(const std::string fname, const std::string oper1, const std::string oper2, const Size lineNr) {
			boost::shared_ptr<MCPayoffT> p1 = payoff(oper1, lineNr);
			if (!p1) {
				scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper1 + "' is no valid operand"));
				return 0;
			}
			boost::shared_ptr<MCPayoffT> p2 = payoff(oper2, lineNr);
			if (!p2) {
				scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper2 + "' is no valid operand"));
				return 0;
			}
			boost::smatch what;
			if (boost::regex_match(fname, what, boost::regex("Min")))
				return boost::shared_ptr<MCPayoffT>(new MCPayoffT::Min(p1,p2));
			if (boost::regex_match(fname, what, boost::regex("Max")))
				return boost::shared_ptr<MCPayoffT>(new MCPayoffT::Max(p1, p2));
			if (boost::regex_match(fname, what, boost::regex("Pay"))) {
				DateType t;
				if (to_Number(oper2, t)) return boost::shared_ptr<MCPayoffT>(new MCPayoffT::Pay(p1, t)); // usual application
				return boost::shared_ptr<MCPayoffT>(new MCPayoffT::Pay(p1, p2->observationTime())); // fall back
			}
			// if we end up here the function name is not valid
			scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + fname + "' is no valid binary function name"));
			return 0;
		}

		// compile three operand function
		inline boost::shared_ptr<MCPayoffT> function3(const std::string fname, const std::string oper1, const std::string oper2, const std::string oper3, const Size lineNr) {
			boost::shared_ptr<MCPayoffT> p1 = payoff(oper1, lineNr);
			if (!p1) {
				scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper1 + "' is no valid operand"));
				return 0;
			}
			boost::shared_ptr<MCPayoffT> p2 = payoff(oper2, lineNr);
			if (!p2) {
				scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper2 + "' is no valid operand"));
				return 0;
			}
			boost::shared_ptr<MCPayoffT> p3 = payoff(oper3, lineNr);
			if (!p3) {
				scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper3 + "' is no valid operand"));
				return 0;
			}
			boost::smatch what;
			if (boost::regex_match(fname, what, boost::regex("IfThenElse")))
				return boost::shared_ptr<MCPayoffT>(new MCPayoffT::IfThenElse(p1, p2, p3));
			// if we end up here the function name is not valid
			scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + fname + "' is no valid function name"));
			return 0;
		}

		// compile unary operators
		inline boost::shared_ptr<MCPayoffT> operator1(const std::string opname, const std::string operand, const Size lineNr) {
			boost::shared_ptr<MCPayoffT> p = payoff(operand, lineNr);
			if (!p) {
				scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + operand + "' is no valid operand"));
				return 0;
			}
			boost::smatch what;
			if (boost::regex_match(opname, what, boost::regex("\\+"))) return p;
			if (boost::regex_match(opname, what, boost::regex("-"))) return boost::shared_ptr<MCPayoffT>(new MCPayoffT::Axpy(-1.0,p,0));
			// if we end up here the function name is not valid
			scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + opname + "' is no valid unary operator name"));
			return 0;
		}

		// compile binary operators
		inline boost::shared_ptr<MCPayoffT> operator2(const std::string opname, const std::string oper1, const std::string oper2, const Size lineNr) {
			boost::shared_ptr<MCPayoffT> p1 = payoff(oper1, lineNr);
			if (!p1) {
				scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper1 + "' is no valid operand"));
				return 0;
			}
			boost::shared_ptr<MCPayoffT> p2 = payoff(oper2, lineNr);
			if (!p2) {
				scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + oper2 + "' is no valid operand"));
				return 0;
			}
			boost::smatch what;
			if (boost::regex_match(opname, what, boost::regex("\\+")))
				return boost::shared_ptr<MCPayoffT>(new MCPayoffT::Axpy(1.0,p1,p2));
			if (boost::regex_match(opname, what, boost::regex("-")))
				return boost::shared_ptr<MCPayoffT>(new MCPayoffT::Axpy(-1.0,p2,p1));
			if (boost::regex_match(opname, what, boost::regex("\\*")))
				return boost::shared_ptr<MCPayoffT>(new MCPayoffT::Mult(p1, p2));
			if (boost::regex_match(opname, what, boost::regex("==|!=|<|<=|>|>=|&&|\\|\\|")))
				return boost::shared_ptr<MCPayoffT>(new MCPayoffT::Logical(p1, p2, opname));
			// if we end up here the function name is not valid
			scriptLog_.push_back(std::string("Error line " + std::to_string(lineNr) + ": '" + opname + "' is no valid binary operator name"));
			return 0;
		}

	};

}

#endif  /* ifndef quantlib_templatemcscript_hpp */ 
