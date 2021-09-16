#ifndef Event_H
#define Event_H

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "Particle.h"

class Event {
 private:
  std::vector<Particle> particles_;
  int nTrack_;

 public:
  Event(int nTrack, std::vector<Particle> particles) : nTrack_(nTrack), particles_(particles){};
  Event() = default;
  Particle getParticle(int ith) { return particles_.at(ith); };
  int getNtrack() { return nTrack_; };
};
#endif
