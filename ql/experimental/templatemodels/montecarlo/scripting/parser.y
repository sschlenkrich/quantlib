%skeleton "lalr1.cc" /* -*- C++ -*- */
%require "3.0.4"
%defines
%define namespace {QuantLib::Scripting}
%define parser_class_name {Parser}
%define api.token.constructor
%define api.value.type variant
%define parse.assert true

%output  "Parser.cpp"
%defines "Parser.hpp"

%code requires
{
#include <string>
#include <boost/shared_ptr.hpp>
namespace QuantLib{
    namespace Scripting {
        class Expression;
        class FlexBisonDriver;
	}
}

}
// The parsing context.
%parse-param { FlexBisonDriver& driver }
%parse-param { void* yyscanner         }
%lex-param   { FlexBisonDriver& driver }
%lex-param   { void* yyscanner         }
%locations
%initial-action
{
  // Initialize the initial location.
  @$.begin.filename = @$.end.filename = &driver.text();
};
// %define parse.trace
// %define parse.error verbose
%debug
%error-verbose
%code
{
// further .cpp includes
#include "Expression.hpp"
#include "FlexBisonDriver.hpp"
// tell Bison that the scanner exists as expected...
QuantLib::Scripting::Parser::symbol_type yylex (QuantLib::Scripting::FlexBisonDriver& driver, void* yyscanner);

}
%define api.token.prefix {TOK_}
%token
  END  0  "end of file"
  ASSIGN  "="
  MINUS   "-"
  PLUS    "+"
  STAR    "*"
  SLASH   "/"
  LPAREN  "("
  RPAREN  ")"
;
%token <std::string> IDENTIFIER "identifier"
%token <std::string> NUMBER     "number"
%type  <boost::shared_ptr<Expression>> exp
%type  <boost::shared_ptr<Expression>> assignment
%type  <boost::shared_ptr<Expression>> assignments

%printer { yyoutput << $$; } <*>;
%%
%start unit;
unit: assignments  { driver.setExpressionTree($1); };

assignments:
  %empty                 {}
| assignments assignment {$$ = boost::shared_ptr<Expression>(new Expression(Expression::NEXT,"",$1,$2));};

assignment:
  "identifier" "=" exp { $$ = boost::shared_ptr<Expression>(new Expression(Expression::ASSIGNMENT,$1,$3)); };

%left "+" "-";
%left "*" "/";
exp:
  exp "+" exp   { $$ = boost::shared_ptr<Expression>(new Expression(Expression::PLUS,"",$1,$3)); }
| exp "-" exp   { $$ = boost::shared_ptr<Expression>(new Expression(Expression::MINUS,"",$1,$3)); }
| exp "*" exp   { $$ = boost::shared_ptr<Expression>(new Expression(Expression::MULT,"",$1,$3)); }
| exp "/" exp   { $$ = boost::shared_ptr<Expression>(new Expression(Expression::DIVISION,"",$1,$3)); }
| "(" exp ")"   { $$ = $2; }
| "identifier"  { $$ = boost::shared_ptr<Expression>(new Expression(Expression::IDENTIFIER,$1)); }
| "number"      { $$ = boost::shared_ptr<Expression>(new Expression(Expression::NUMBER,$1)); };
%%

void QuantLib::Scripting::Parser::error (const location_type& l, const std::string& m) {
  driver.error (l, m);
}

// end of file