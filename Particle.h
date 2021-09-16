#ifndef Particle_H
#define Particle_H

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"


class Particle {
 public:
  Particle(int pid, int charge, float energy, float rapidity, float pt,
           float phi)
      : pid_(pid),
        charge_(charge),
        energy_(energy),
        rapidity_(rapidity),
        pt_(pt),
        phi_(phi){};
  Particle() = default;

  int getPid() { return pid_; }
  int getCharge() { return charge_; }
  double getEnergy() { return energy_; }
  double getRapidity() { return rapidity_; }
  double getPt() {return pt_; }
  double getPhi() {return phi_; }

 private:
  int pid_{0};
  int charge_{0};
  float energy_{-999.};
  float rapidity_{-999.};
  float pt_{-999.};
  float phi_{-999};
};
#endif
