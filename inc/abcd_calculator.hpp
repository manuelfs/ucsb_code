#ifndef H_ABCD_CALCULATOR
#define H_ABCD_CALCULATOR

#include <vector>
#include "abcd_count.hpp"

class ABCDCalculator{
public:
  ABCDCalculator();
  ABCDCalculator(const std::vector<ABCDCount>&,
		 const std::vector<ABCDCount>&);

  typedef ABCDCount::size_type size_type;

  double operator()(std::vector<double>) const;

  size_type GetNumberOfBins() const;
  unsigned int GetNumberOfParameters() const;

private:
  std::vector<ABCDCount> observed_, signal_;
  mutable size_type number_of_bins_;
  mutable bool number_of_bins_ready_;
};

#endif
