#include "b_jet.hpp"
#include "TLorentzVector.h"

BJet::BJet(const TLorentzVector vecIn, const double bTagIn):vec(vecIn),bTag(bTagIn){
}

void BJet::SetLorentzVector(const TLorentzVector vecIn){
  vec=vecIn;
}

void BJet::SetBTag(const double bTagIn){
  bTag=bTagIn;
}

TLorentzVector BJet::GetLorentzVector() const{
  return vec;
}

double BJet::GetBTag() const{
  return bTag;
}

bool BJet::operator==(const BJet &jet) const{
  return vec==jet.vec && bTag==jet.bTag;
}

bool BJet::operator!=(const BJet &jet) const{
  return !(*this==jet);
}

bool BJet::operator<(const BJet &jet) const{
  return bTag<jet.bTag || (bTag==jet.bTag && vec.Pt()<jet.vec.Pt());
}

bool BJet::operator>(const BJet &jet) const{
  return bTag>jet.bTag || (bTag==jet.bTag && vec.Pt()>jet.vec.Pt());
}

bool BJet::operator<=(const BJet &jet) const{
  return *this==jet || *this<jet;
}

bool BJet::operator>=(const BJet &jet) const{
  return *this==jet || *this>jet;
}
