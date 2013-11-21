#include "abcd_count.hpp"
#include <vector>
#include <stdexcept>
#include "math.hpp"

ABCDCount::ABCDCount():
  a_(0),
  b_(0),
  c_(0),
  d_(0),
  weight_(0.0),
  number_of_bins_(0),
  number_of_bins_ready_(false){
}

ABCDCount::ABCDCount(const std::vector<double> &a_in, const std::vector<double> &b_in,
		     const std::vector<double> &c_in, const std::vector<double> &d_in,
		     const double weight_in):
  a_(a_in),
  b_(b_in),
  c_(c_in),
  d_(d_in),
  weight_(weight_in),
  number_of_bins_(0),
  number_of_bins_ready_(false){
  }

ABCDCount::size_type ABCDCount::GetNumberOfBins() const{
  if(!number_of_bins_ready_){
    number_of_bins_=(a_.size());
    if(b_.size()<number_of_bins_) number_of_bins_=b_.size();
    if(c_.size()<number_of_bins_) number_of_bins_=c_.size();
    if(d_.size()<number_of_bins_) number_of_bins_=d_.size();
    number_of_bins_ready_=true;
  }
  return number_of_bins_;
}

double ABCDCount::GetWeight() const{
  return weight_;
}

void ABCDCount::SetWeight(const double weight_in){
  weight_=weight_in;
}

double& ABCDCount::operator()(const Region& region, const size_type& bin){
  if(bin>GetNumberOfBins()){
    throw std::out_of_range("out_of_range in ABCDCount::operator()");
  }

  if(region==kA){
    return a_.at(bin);
  }else if(region==kB){
    return b_.at(bin);
  }else if(region==kC){
    return c_.at(bin);
  }else{
    return d_.at(bin);
  }
}

double ABCDCount::GetCount(const Region& region, const size_type& bin) const{
  if(bin>GetNumberOfBins()){
    throw std::out_of_range("out_of_range in ABCDCount::operator()");
  }

  if(region==kA){
    return a_.at(bin);
  }else if(region==kB){
    return b_.at(bin);
  }else if(region==kC){
    return c_.at(bin);
  }else{
    return d_.at(bin);
  }
}

void ABCDCount::SetCounts(const std::vector<double>& a_in, const std::vector<double>& b_in,
			  const std::vector<double>& c_in, const std::vector<double>& d_in){
  number_of_bins_ready_=false;
  a_=a_in;
  b_=b_in;
  c_=c_in;
  d_=d_in;
}

void ABCDCount::SetCount(const Region& region, const std::vector<double>& counts_in){
  number_of_bins_ready_=false;
  if(region==kA){
    a_=counts_in;
  }else if(region==kB){
    b_=counts_in;
  }else if(region==kC){
    c_=counts_in;
  }else{
    d_=counts_in;
  }
}

double ABCDCount::GetTotalCount() const{
  std::vector<double> sums(4);
  sums.at(0)=Math::Sum(a_.begin(), a_.begin()+GetNumberOfBins());
  sums.at(1)=Math::Sum(b_.begin(), b_.begin()+GetNumberOfBins());
  sums.at(2)=Math::Sum(c_.begin(), c_.begin()+GetNumberOfBins());
  sums.at(3)=Math::Sum(d_.begin(), d_.begin()+GetNumberOfBins());
  return Math::Sum(sums.begin(), sums.end());
}
