#include "abcd_calculator.hpp"

#include <cmath>
#include <vector>
#include <limits>
#include "abcd_count.hpp"
#include "math.hpp"

ABCDCalculator::ABCDCalculator():
  observed_(0),
  signal_(0),
  number_of_bins_(0),
  number_of_bins_ready_(false){
}

ABCDCalculator::ABCDCalculator(const std::vector<ABCDCount>& observed_in,
                               const std::vector<ABCDCount>& signal_in):
  observed_(observed_in),
  signal_(signal_in),
  number_of_bins_(0),
  number_of_bins_ready_(false){
}

double ABCDCalculator::operator()(std::vector<double> params) const{
  if(params.size()!=GetNumberOfParameters()) return std::numeric_limits<double>::max();
  for(unsigned int i(0); i<params.size()-1; ++i){
    params.at(i)=fabs(params.at(i));
  }
  if(params.at(params.size()-1)<std::numeric_limits<double>::min()){
    params.at(params.size()-1)=std::numeric_limits<double>::min();
  }
  std::vector<double> observed_frac(observed_.size());
  for(size_type i(0); i<observed_.size(); ++i){
    observed_frac.at(i)=observed_.at(i).GetTotalCount();
  }
  double normalization(Math::Sum(observed_frac.begin(), observed_frac.end()));
  if(normalization>0.0){
    normalization=1.0/normalization;
  }else{
    normalization=1.0;
  }
  for(std::vector<double>::iterator it(observed_frac.begin()); it!=observed_frac.end(); ++it){
    (*it)*=normalization;
  }
  std::vector<double> signal_frac(signal_.size());
  for(size_type i(0); i<signal_.size(); ++i){
    signal_frac.at(i)=signal_.at(i).GetTotalCount();
  }
  normalization=Math::Sum(signal_frac.begin(), signal_frac.end());
  if(normalization>0.0){
    normalization=1.0/normalization;
  }else{
    normalization=1.0;
  }
  for(std::vector<double>::iterator it(signal_frac.begin()); it!=signal_frac.end(); ++it){
    (*it)*=normalization;
  }
  std::vector<double> observed_counts(4*GetNumberOfBins()*observed_.size());
  std::vector<double> observed_means(observed_counts.size());
  std::vector<double> observed_weights(observed_counts.size());
  std::vector<double> signal_counts(4*GetNumberOfBins()*signal_.size());
  std::vector<double> signal_means(signal_counts.size());
  std::vector<double> signal_weights(signal_counts.size());
  std::vector<double> log_likelihoods(2*(observed_counts.size()+signal_counts.size()));
  const double chi_b(params.at(5*GetNumberOfBins()));
  const double chi_c(params.at(5*GetNumberOfBins()+1));
  const double alpha(params.at(5*GetNumberOfBins()+2));
  for(size_type bin(0); bin<GetNumberOfBins(); ++bin){
    for(size_type sample(0); sample<observed_.size(); ++sample){
      const size_type index(4*(GetNumberOfBins()*sample+bin));
      observed_counts.at(index)=observed_.at(sample).GetCount(ABCDCount::kA, bin);
      observed_counts.at(index+1)=observed_.at(sample).GetCount(ABCDCount::kB, bin);
      observed_counts.at(index+2)=observed_.at(sample).GetCount(ABCDCount::kC, bin);
      observed_counts.at(index+3)=observed_.at(sample).GetCount(ABCDCount::kD, bin);
      observed_means.at(index)=observed_frac.at(sample)
        *(params.at(5*bin)+alpha*params.at(5*bin+1));
      observed_means.at(index+1)=observed_frac.at(sample)
        *(params.at(5*bin)*chi_b+alpha*params.at(5*bin+2));
      observed_means.at(index+2)=observed_frac.at(sample)
        *(params.at(5*bin)*chi_c+alpha*params.at(5*bin+3));
      observed_means.at(index+3)=observed_frac.at(sample)
        *(params.at(5*bin)*chi_b*chi_c+alpha*params.at(5*bin+4));
      observed_weights.at(index)=observed_.at(sample).GetWeight();
      observed_weights.at(index+1)=observed_.at(sample).GetWeight();
      observed_weights.at(index+2)=observed_.at(sample).GetWeight();
      observed_weights.at(index+3)=observed_.at(sample).GetWeight();
    }
    for(size_type sample(0); sample<signal_.size(); ++sample){
      const size_type index(4*(GetNumberOfBins()*sample+bin));
      signal_counts.at(index)=signal_.at(sample).GetCount(ABCDCount::kA, bin);
      signal_counts.at(index+1)=signal_.at(sample).GetCount(ABCDCount::kB, bin);
      signal_counts.at(index+2)=signal_.at(sample).GetCount(ABCDCount::kC, bin);
      signal_counts.at(index+3)=signal_.at(sample).GetCount(ABCDCount::kD, bin);
      signal_means.at(index)=signal_frac.at(sample)*params.at(5*bin+1);
      signal_means.at(index+1)=signal_frac.at(sample)*params.at(5*bin+2);
      signal_means.at(index+2)=signal_frac.at(sample)*params.at(5*bin+3);
      signal_means.at(index+3)=signal_frac.at(sample)*params.at(5*bin+4);
      signal_weights.at(index)=signal_.at(sample).GetWeight();
      signal_weights.at(index+1)=signal_.at(sample).GetWeight();
      signal_weights.at(index+2)=signal_.at(sample).GetWeight();
      signal_weights.at(index+3)=signal_.at(sample).GetWeight();
    }
  }
  for(size_type i(0); i<observed_counts.size(); ++i){
    log_likelihoods.at(2*i)=observed_counts.at(i)*std::log(observed_means.at(i))
      /observed_weights.at(i);
    log_likelihoods.at(2*i+1)=-observed_means.at(i)/observed_weights.at(i);
  }
  const size_type offset(2*observed_counts.size());
  for(size_type i(0); i<signal_counts.size(); ++i){
    log_likelihoods.at(offset+2*i)=signal_counts.at(i)*std::log(signal_means.at(i))
      /signal_weights.at(i);
    log_likelihoods.at(offset+2*i+1)=-signal_means.at(i)/signal_weights.at(i);
  }
  return -2.0*Math::Sum(log_likelihoods.begin(), log_likelihoods.end());
}

unsigned int ABCDCalculator::GetNumberOfParameters() const{
  return 5*GetNumberOfBins()+3;
}

ABCDCalculator::size_type ABCDCalculator::GetNumberOfBins() const{
  if(!number_of_bins_ready_){
    if(observed_.size()!=0){
      number_of_bins_=observed_.at(0).GetNumberOfBins();
    }else if(signal_.size()!=0){
      number_of_bins_=signal_.at(0).GetNumberOfBins();
    }else{
      number_of_bins_=0;
    }
    std::vector<ABCDCount>::const_iterator it;
    for(it=observed_.begin(); it!=observed_.end(); ++it){
      if(it->GetNumberOfBins()<number_of_bins_){
        number_of_bins_=it->GetNumberOfBins();
      }
    }
    for(it=signal_.begin(); it!=signal_.end(); ++it){
      if(it->GetNumberOfBins()<number_of_bins_){
        number_of_bins_=it->GetNumberOfBins();
      }
    }
    number_of_bins_ready_=true;
  }
  return number_of_bins_;
}
