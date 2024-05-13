// PowhegHooksBB4L.h
// Copyright (C) 2017 Silvia Ferrario Ravasio, Tomas Jezo, Paolo Nason, Markus Seidel
// inspired by PowhegHooks.h by Richard Corke
// adjusted to work with EmissionVetoHook1 in CMSSW by Alexander Grohsjean

#ifndef Pythia8_PowhegHooksBB4L_H
#define Pythia8_PowhegHooksBB4L_H

// Includes
#include "Pythia8/Pythia.h"
#include <cassert>
struct {
  int radtype;
} radtype_;

namespace Pythia8 {

  class PowhegHooksBB4L : public UserHooks {
  public:
    //--- Constructor and destructor -------------------------------------------
    PowhegHooksBB4L() : nFSRvetoBB4l(0) {}
    ~PowhegHooksBB4L() override { std::cout << "Number of FSR vetoed in BB4l = " << nFSRvetoBB4l << std::endl; }

    //--- Initialization -------------------------------------------------------
    bool initAfterBeams() override {
      vetoFSREmission = settingsPtr->flag("POWHEG:bb4l:FSREmission:veto");
      onlyDistance1 = settingsPtr->flag("POWHEG:bb4l:FSREmission:onlyDistance1");
      dryRunFSR = settingsPtr->flag("POWHEG:bb4l:FSREmission:dryRun");
      vetoAtPL = settingsPtr->flag("POWHEG:bb4l:FSREmission:vetoAtPL");
      vetoQED = settingsPtr->flag("POWHEG:bb4l:FSREmission:vetoQED");
      vetoPartonLevel = settingsPtr->flag("POWHEG:bb4l:PartonLevel:veto");
      excludeFSRConflicting = settingsPtr->flag("POWHEG:bb4l:PartonLevel:excludeFSRConflicting");
      debug = settingsPtr->flag("POWHEG:bb4l:DEBUG");
      scaleResonanceVeto = settingsPtr->flag("POWHEG:bb4l:ScaleResonance:veto");
      vetoDipoleFrame = settingsPtr->flag("POWHEG:bb4l:FSREmission:vetoDipoleFrame");
      pTpythiaVeto = settingsPtr->flag("POWHEG:bb4l:FSREmission:pTpythiaVeto");
      pTmin = settingsPtr->parm("POWHEG:bb4l:pTminVeto");
      return true;
    }

    //--- PROCESS LEVEL HOOK ---------------------------------------------------
 
    // called at the LHE level
    inline bool canVetoProcessLevel() override { return true; }
    inline bool doVetoProcessLevel(Event &e) override {
      // extract the radtype from the event comment
      stringstream ss;
      // use eventattribute as comments not filled when using edm input
      //ss << infoPtr->getEventComments();      
      ss << infoPtr->getEventAttribute("#rwgt");
      string temp;
      ss >> temp >> radtype_.radtype;
      assert(temp == "#rwgt");
      // find last top and the last anti-top in the record
      int i_top = -1, i_atop = -1;
      for (int i = 0; i < e.size(); i++) {
        if (e[i].id() == 6)
          i_top = i;
        if (e[i].id() == -6)
          i_atop = i;
        if (e[i].id() == 24)
          i_wp = i;
        if (e[i].id() == -24)
          i_wm = i;
      }
      // if found calculate the resonance scale
      topresscale = findresscale(i_top, e);
      // similarly for anti-top
      atopresscale = findresscale(i_atop, e);
      // and for W^+ and W^-
      wpresscale = findresscale(i_wp, e);
      wmresscale = findresscale(i_wm, e);

      // do not veto, ever
      return false;
    }

    //--- FSR EMISSION LEVEL HOOK ----------------------------------------------
    // This hook gets triggered everytime the parton shower attempts to attach
    // a FSR emission.
    inline bool canVetoFSREmission() { return vetoFSREmission; }
    inline bool doVetoFSREmission(int sizeOld, const Event &e, int iSys, bool inResonance) {
      // FSR VETO INSIDE THE RESONANCE (if it is switched on)
      if (inResonance && vetoFSREmission) {
        // get the participants of the splitting: the recoiler, the radiator and the emitted
        int iRecAft = e.size() - 1;
        int iEmt = e.size() - 2;
        int iRadAft = e.size() - 3;
        int iRadBef = e[iEmt].mother1();

        // find the resonance the radiator originates from
        int iRes = e[iRadBef].mother1();
        while (iRes > 0 && (abs(e[iRes].id()) != 6 && abs(e[iRes].id()) != 24)) {
          iRes = e[iRes].mother1();
        }
        if (iRes == 0) {
          infoPtr->errorMsg(
              "Warning in PowhegHooksBB4L::doVetoFSREmission: emission in resonance not from the top quark or from the "
              "W boson, not vetoing");
          return doVetoFSR(false, 0);
        }
        int iResId = e[iRes].id();

        // calculate the scale of the emission
        double scale;
        //using pythia pT definition ...
        if (pTpythiaVeto)
          scale = pTpythia(e, iRadAft, iEmt, iRecAft);
        //.. or using POWHEG pT definition
        else {
          Vec4 pr(e[iRadAft].p()), pe(e[iEmt].p()), pres(e[iRes].p()), prec(e[iRecAft].p()), psystem;
          // The computation of the POWHEG pT can be done in the top rest frame or in the diple one.
          // pdipole = pemt +prec +prad (after the emission)
          // For the first emission off the top resonance pdipole = pw +pb (before the emission) = ptop
          if (vetoDipoleFrame)
            psystem = pr + pe + prec;
          else
            psystem = pres;

          // gluon splitting into two partons
          if (e[iRadBef].id() == 21)
            scale = gSplittingScale(psystem, pr, pe);
          // quark emitting a gluon (or a photon)
          else if (abs(e[iRadBef].id()) <= 5 && ((e[iEmt].id() == 21) && !vetoQED))
            scale = qSplittingScale(psystem, pr, pe);
          // other stuff (which we should not veto)
          else {
            scale = 0;
          }
        }

        // compare the current splitting scale to the correct resonance scale
        if (iResId == 6) {
          if (debug && scale > topresscale)
            cout << iResId << ": " << e[iRadBef].id() << " > " << e[iRadAft].id() << " + " << e[iEmt].id() << "; "
                 << scale << endl;
          return doVetoFSR(scale > topresscale, scale);
        } else if (iResId == -6) {
          if (debug && scale > atopresscale)
            cout << iResId << ": " << e[iRadBef].id() << " > " << e[iRadAft].id() << " + " << e[iEmt].id() << "; "
                 << scale << endl;
          return doVetoFSR(scale > atopresscale, scale);
        } else if (iResId == 24) {
          if (debug && scale > wpresscale)
            cout << iResId << ": " << e[iRadBef].id() << " > " << e[iRadAft].id() << " + " << e[iEmt].id() << "; "
                 << scale << endl;
          return doVetoFSR(scale > wpresscale, scale);
        } else if (iResId == -24) {
          if (debug && scale > wmresscale)
            cout << iResId << ": " << e[iRadBef].id() << " > " << e[iRadAft].id() << " + " << e[iEmt].id() << "; "
                 << scale << endl;
          return doVetoFSR(scale > wmresscale, scale);
        } else {
          infoPtr->errorMsg("Error in PowhegHooksBB4L::doVetoFSREmission: unimplemented case");
          exit(-1);
        }
      }
      // In CMSSW, the production process veto is done in EmissionVetoHook1.cc
      // so for events outside resonance, nothing needs to be done here
      else {
        return false;
      }
    }

    inline bool doVetoFSR(bool condition, double scale) {
      if (!vetoAllRadtypes && radtype == 2)
        return false;
      if (condition) {
        nInResonanceFSRveto++;
        return true;
      }
      return false;
    }

    //--- SCALE RESONANCE HOOK -------------------------------------------------
    // called before each resonance decay shower
    inline bool canSetResonanceScale() { return scaleResonanceVeto; }
    // if the resonance is the (anti)top or W+/W- set the scale to:
    // - if radtype=2 (remnant): resonance virtuality
    // - if radtype=1 (btilde):
    //  - (a)topresscale/wp(m)resscale for tops and Ws
    //  - a large number otherwise
    // if is not the top, set it to a big number
    inline double scaleResonance(int iRes, const Event &e) {
      if (!vetoAllRadtypes && radtype == 2)
        return sqrt(e[iRes].m2Calc());
      else {
        if (e[iRes].id() == 6)
          return topresscale;
        else if (e[iRes].id() == -6)
          return atopresscale;
        else if (e[iRes].id() == 24)
          return wpresscale;
        else if (e[iRes].id() == 24)
          return wmresscale;
        else
          return 1e30;
      }
    }

    //--- Internal helper functions --------------------------------------------
    // Calculates the scale of the hardest emission from within the resonance system
    // translated by Markus Seidel modified by Tomas Jezo
    inline double findresscale(const int iRes, const Event &event) {
      // return large scale if the resonance position is ill defined
      if (iRes < 0)
        return 1e30;

      // get number of resonance decay products
      int nDau = event[iRes].daughterList().size();

      // iRes is not decayed, return high scale equivalent to
      // unrestricted shower
      if (nDau == 0) {
        return 1e30;
      }
      // iRes did not radiate, this means that POWHEG pt scale has
      // evolved all the way down to pTmin
      else if (nDau < 3) {
        return pTmin;
      }
      // iRes is a (anti-)top quark
      else if (abs(event[iRes].id()) == 6) {
        // find top daughters
        int idw = -1, idb = -1, idg = -1;
        for (int i = 0; i < nDau; i++) {
          int iDau = event[iRes].daughterList()[i];
          if (abs(event[iDau].id()) == 24)
            idw = iDau;
          if (abs(event[iDau].id()) == 5)
            idb = iDau;
          if (abs(event[iDau].id()) == 21)
            idg = iDau;
        }

        // Get daughter 4-vectors in resonance frame
        Vec4 pw(event[idw].p());
        pw.bstback(event[iRes].p());
        Vec4 pb(event[idb].p());
        pb.bstback(event[iRes].p());
        Vec4 pg(event[idg].p());
        pg.bstback(event[iRes].p());

        // Calculate scale and return it
        return sqrt(2 * pg * pb * pg.e() / pb.e());
      }
      // iRes is a W+(-) boson
      else if (abs(event[iRes].id()) == 24) {
        // Find W daughters
        int idq = -1, ida = -1, idg = -1;
        for (int i = 0; i < nDau; i++) {
          int iDau = event[iRes].daughterList()[i];
          if (event[iDau].id() == 21)
            idg = iDau;
          else if (event[iDau].id() > 0)
            idq = iDau;
          else if (event[iDau].id() < 0)
            ida = iDau;
        }

        // Get daughter 4-vectors in resonance frame
        Vec4 pq(event[idq].p());
        pq.bstback(event[iRes].p());
        Vec4 pa(event[ida].p());
        pa.bstback(event[iRes].p());
        Vec4 pg(event[idg].p());
        pg.bstback(event[iRes].p());

        // Calculate scale
        Vec4 pw = pq + pa + pg;
        double q2 = pw * pw;
        double csi = 2 * pg.e() / sqrt(q2);
        double yq = 1 - pg * pq / (pg.e() * pq.e());
        double ya = 1 - pg * pa / (pg.e() * pa.e());
        // and return it
        return sqrt(min(1 - yq, 1 - ya) * pow2(csi) * q2 / 2);
      }
      // in any other case just return a high scale equivalent to
      // unrestricted shower
      return 1e30;
    }

    // The following routine will match daughters of particle `e[iparticle]` to an expected pattern specified via the list of expected particle PDG ID's `ids`,
    // id wildcard is specified as 0 if match is obtained, the positions and the momenta of these particles are returned in vectors `positions` and `momenta`
    // respectively
    // if exitOnExtraLegs==true, it will exit if the decay has more particles than expected, but not less
    inline bool match_decay(int iparticle,
                            const Event &e,
                            const vector<int> &ids,
                            vector<int> &positions,
                            vector<Vec4> &momenta,
                            bool exitOnExtraLegs = true) {
      // compare sizes
      if (e[iparticle].daughterList().size() != ids.size()) {
        if (exitOnExtraLegs && e[iparticle].daughterList().size() > ids.size()) {
          cout << "extra leg" << endl;
          exit(-1);
        }
        return false;
      }
      // compare content
      for (unsigned i = 0; i < e[iparticle].daughterList().size(); i++) {
        int di = e[iparticle].daughterList()[i];
        if (ids[i] != 0 && e[di].id() != ids[i])
          return false;
      }
      // reset the positions and momenta vectors (because they may be reused)
      positions.clear();
      momenta.clear();
      // construct the array of momenta
      for (unsigned i = 0; i < e[iparticle].daughterList().size(); i++) {
        int di = e[iparticle].daughterList()[i];
        positions.push_back(di);
        momenta.push_back(e[di].p());
      }
      return true;
    }

    inline double qSplittingScale(Vec4 pt, Vec4 p1, Vec4 p2) {
      p1.bstback(pt);
      p2.bstback(pt);
      return sqrt(2 * p1 * p2 * p2.e() / p1.e());
    }

    inline double gSplittingScale(Vec4 pt, Vec4 p1, Vec4 p2) {
      p1.bstback(pt);
      p2.bstback(pt);
      return sqrt(2 * p1 * p2 * p1.e() * p2.e() / (pow(p1.e() + p2.e(), 2)));
    }

    // Routines to calculate the pT (according to pTdefMode) in a FS splitting:
    // i (radiator before) -> j (emitted after) k (radiator after)
    // For the Pythia pT definition, a recoiler (after) must be specified.
    // (INSPIRED BY pythia8F77_31.cc double pTpythia)
    inline double pTpythia(const Event &e, int RadAfterBranch, int EmtAfterBranch, int RecAfterBranch) {
      // Convenient shorthands for later
      Vec4 radVec = e[RadAfterBranch].p();
      Vec4 emtVec = e[EmtAfterBranch].p();
      Vec4 recVec = e[RecAfterBranch].p();
      int radID = e[RadAfterBranch].id();

      // Calculate virtuality of splitting
      Vec4 Q(radVec + emtVec);
      double Qsq = Q.m2Calc();

      // Mass term of radiator
      double m2Rad = (abs(radID) >= 4 && abs(radID) < 7) ? pow2(particleDataPtr->m0(radID)) : 0.;

      // z values for FSR
      double z, pTnow;
      // Construct 2 -> 3 variables
      Vec4 sum = radVec + recVec + emtVec;
      double m2Dip = sum.m2Calc();

      double x1 = 2. * (sum * radVec) / m2Dip;
      double x3 = 2. * (sum * emtVec) / m2Dip;
      z = x1 / (x1 + x3);
      pTnow = z * (1. - z);

      // Virtuality
      pTnow *= (Qsq - m2Rad);

      if (pTnow < 0.) {
        cout << "Warning: pTpythia was negative" << endl;
        return -1.;
      } else
        return (sqrt(pTnow));
    }

    inline double getdechardness(int topcharge, const Event &e) {
      int tid = 6 * topcharge, wid = 24 * topcharge, bid = 5 * topcharge, gid = 21, wildcard = 0;
      // find last top in the record
      int i_top = -1;
      Vec4 p_top, p_b, p_g, p_g1, p_g2;
      for (int i = 0; i < e.size(); i++)
        if (e[i].id() == tid) {
          i_top = i;
          p_top = e[i].p();
        }
      if (i_top == -1)
        return -1.0;

      // summary of cases
      // 1.) t > W b
      //   a.) b > 3     ... error
      //   b.) b > b g   ... h = sqrt(2*p_g*p_b*p_g.e()/p_b.e())
      //   c.) b > other ... h = -1
      //   return h
      // 2.) t > W b g
      //   a.)   b > 3     ... error
      //   b.)   b > b g   ... h1 = sqrt(2*p_g*p_b*p_g.e()/p_b.e())
      //   c.)   b > other ... h1 = -1
      //   i.)   g > 3     ... error
      //   ii.)  g > 2     ... h2 = sqrt(2*p_g1*p_g2*p_g1.e()*p_g2.e()/(pow(p_g1.e(),2)+pow(p_g2.e(),2))) );
      //   iii.) g > other ... h2 = -1
      //   return max(h1,h2)
      // 3.) else ... error

      vector<Vec4> momenta;
      vector<int> positions;

      // 1.) t > b W
      if (match_decay(i_top, e, vector<int>{wid, bid}, positions, momenta, false)) {
        double h;
        int i_b = positions[1];
        // a.+b.) b > 3 or b > b g
        if (match_decay(i_b, e, vector<int>{bid, gid}, positions, momenta))
          h = qSplittingScale(e[i_top].p(), momenta[0], momenta[1]);
        // c.) b > other
        else
          h = -1;
        return h;
      }
      // 2.) t > b W g
      else if (match_decay(i_top, e, vector<int>{wid, bid, gid}, positions, momenta, false)) {
        double h1, h2;
        int i_b = positions[1], i_g = positions[2];
        // a.+b.) b > 3 or b > b g
        if (match_decay(i_b, e, vector<int>{bid, gid}, positions, momenta))
          h1 = qSplittingScale(e[i_top].p(), momenta[0], momenta[1]);
        // c.) b > other
        else
          h1 = -1;
        // i.+ii.) g > 3 or g > 2
        if (match_decay(i_g, e, vector<int>{wildcard, wildcard}, positions, momenta))
          h2 = gSplittingScale(e[i_top].p(), momenta[0], momenta[1]);
        // c.) b > other
        else
          h2 = -1;
        return max(h1, h2);
      }
      // 3.) else
      else {
        cout << "getdechardness" << endl;
        cout << "top at position " << i_top << endl;
        cout << "with " << e[i_top].daughterList().size() << " daughters " << endl;
        for (unsigned i = 0; i < e[i_top].daughterList().size(); i++) {
          int di = e[i_top].daughterList()[i];
          cout << "with daughter " << di << ": " << e[di].id() << endl;
        }
        exit(-1);
      }
    }

    //--------------------------------------------------------------------------

    // Functions to return information

    //  inline int    getNFSRveto() { return nFSRveto; }
    //--------------------------------------------------------------------------

  private:
    // FSR emission veto flags
    bool vetoFSREmission, dryRunFSR, wouldVetoFsr, onlyDistance1, vetoAtPL, vetoQED;
    // Parton Level veto flags
    bool vetoPartonLevel, excludeFSRConflicting;
    // Scale Resonance veto flags
    double scaleResonanceVeto;
    // other flags
    bool debug;
    // internal: resonance scales
    double topresscale, atopresscale;
    // internal: inter veto communication
    double vetoDecScale;
    int vetoTopCharge;
    bool vetoDipoleFrame;
    bool pTpythiaVeto;
    //bool vetoProduction;
    double pTmin;
    // Statistics on vetos
    unsigned long int nFSRvetoBB4l;
  };

  //==========================================================================

}  // end namespace Pythia8

#endif  // end Pythia8_PowhegHooksBB4L_H
