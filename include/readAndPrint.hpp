#ifndef __READANDPRINT_H_
#define __READANDPRINT_H_

#include "includes.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
//#include <GetPot>


int
readInput(ddpDomain_type & domain, ddpGrid_type & grid,
          ddpProblemInfo_type & problem,
          ddpCarrierConstants_type & carrierConstants,
          char const * nameOfFile);

int
ddpPrintGrid(ddpGrid_type const & grid);

int ddpPrintTimeStamps(std::vector< double > const & timeStamps);

int
ddpPrintBij_Backward(ddpBijBack_type const & bij);

#endif
