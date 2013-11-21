#ifndef H_BJET
#define H_BJET

#include <cfloat>
#include "TLorentzVector.h"

class BJet{
public:
  BJet(const TLorentzVector=TLorentzVector(0.0,0.0,0.0,0.0), const double bTagIn=-DBL_MAX);

  void SetLorentzVector(const TLorentzVector vecIn);
  void SetBTag(const double bTagIn);

  TLorentzVector GetLorentzVector() const;
  double GetBTag() const;

  bool operator==(const BJet &jet) const;
  bool operator!=(const BJet &jet) const;
  bool operator<(const BJet &jet) const;
  bool operator>(const BJet &jet) const;
  bool operator<=(const BJet &jet) const;
  bool operator>=(const BJet &jet) const;
private:
  TLorentzVector vec;
  double bTag;
};

#endif
