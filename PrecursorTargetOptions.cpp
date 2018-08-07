//
// Created by Dennis Goldfarb on 12/29/17.
//

#include "PrecursorTargetOptions.h"

void PrecursorTargetOptions::addAbundance(PrecursorTargetOption o)
{
  PrecursorTargetOptionsKey key(o.charge, o.minIso, o.maxIso);
  //key2option.find(key)->second.abundance += o.abundance;
  key2option.find(key)->second.abundance = std::max(key2option.find(key)->second.abundance, o.abundance);
  
  /*for(int i = o.minIso; i <=2; ++i)
  { 
      PrecursorTargetOptionsKey key2(o.charge, o.minIso, i);
      auto match = key2option.find(key2);
      if (match != key2option.end() && match->second.abundance < o.abundance)
      {
	key2option.erase(match);
	key2option.insert(std::make_pair(key, o));
	return;
      }
      } */
}

void PrecursorTargetOptions::addOption(PrecursorTargetOption o)
{
    PrecursorTargetOptionsKey key(o.charge, o.minIso, o.maxIso);
    key2option.insert(std::make_pair(key, o));
}

bool PrecursorTargetOptions::hasOption(PrecursorTargetOption o)
{
  PrecursorTargetOptionsKey key(o.charge, o.minIso, o.maxIso);
  return key2option.find(key) != key2option.end();

  /*  for (int i = o.minIso; i <= 2; ++i)
    {
      PrecursorTargetOptionsKey key(o.charge, o.minIso, i);
      if (key2option.find(key) != key2option.end()) return true;
    }

    return false;*/
}

PrecursorTargetOptions::PrecursorTargetOptions(const PrecursorTargetOptions &o)
{
    key2option = o.key2option;
}

PrecursorTargetOptions &PrecursorTargetOptions::operator=(const PrecursorTargetOptions &o)
{
    key2option = o.key2option;
    return *this;
}

int PrecursorTargetOptions::getMaxIsotope() {
    int maxIsotope = 0;
    for (auto key : key2option) {
        maxIsotope = std::max(maxIsotope, key.first.lastIso);
    }
    return maxIsotope;
}

int PrecursorTargetOptions::getMaxCharge() {
    int maxCharge = 0;
    for (auto key : key2option) {
        maxCharge = std::max(maxCharge, key.first.z);
    }
    return maxCharge;
}
