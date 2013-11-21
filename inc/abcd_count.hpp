#ifndef H_ABCD_COUNT
#define H_ABCD_COUNT

#include<vector>

class ABCDCount{
public:
  ABCDCount();
  ABCDCount(const std::vector<double>&, const std::vector<double>&,
	    const std::vector<double>&, const std::vector<double>&,
	    const double=1.0);

  typedef std::vector<double>::size_type size_type;
  enum Region{kA, kB, kC, kD};

  size_type GetNumberOfBins() const;
  double GetWeight() const;
  void SetWeight(const double);

  double& operator()(const Region&,  const size_type&);
  double GetCount(const Region&, const size_type&) const;

  void SetCounts(const std::vector<double>&, const std::vector<double>&,
		 const std::vector<double>&, const std::vector<double>&);

  void SetCount(const Region&, const std::vector<double>&);

  double GetTotalCount() const;
private:
  std::vector<double> a_, b_, c_, d_;
  double weight_;
  mutable size_type number_of_bins_;
  mutable bool number_of_bins_ready_;
};

#endif
