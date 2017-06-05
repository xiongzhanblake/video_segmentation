/*
 * CommonDefination.hpp
 *
 *  Created on: Mar 2, 2017
 *      Author: xiong
 */

#ifndef SRC_COMMONNAMESPACE_HPP_
#define SRC_COMMONNAMESPACE_HPP_

namespace CAM{
	typedef std::vector<double> DoubleVector;
	typedef std::vector<DoubleVector> VectorSet;
	typedef std::map<double,VectorSet> MapSet;
	typedef std::pair<double,CAM::VectorSet> Pair;
}


#endif /* SRC_COMMONNAMESPACE_HPP_ */
