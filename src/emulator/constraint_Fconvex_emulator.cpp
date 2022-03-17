/*
 *  Main authors:
 *      Tias Guns <tias.guns@cs.kuleuven.be>
 *
 *  Copyright:
 *      Tias Guns, 2008
 *
 *  Revision information:
 *      $Id: graph_cp.cc 128 2008-08-07 11:23:23Z tias $
 *
 *  This file is part of Graph_cp and uses Gecode.
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include "gecode/kernel.hh"
#include "gecode/int.hh"
#include "gecode/iter.hh"
#include "gecode/int/rel.hh"

namespace FconvexBool_emulator {
  using namespace ::Gecode;
  using namespace ::Gecode::Int;

  /**
   * \brief Base-class for convex boolean propagators
   *
   * This class NEEDS the global functions:
   * inline bool convex_functionWrapper(
   *    int pos_total, int pos_min, int pos_max, int neg_total, int neg_min, int neg_max, int precision, int tresholdLE)
   */
  template <class VA, class VC>
  class FconvexBool_emulatorGq : public Propagator {
  protected:
    /// total class 1
    int X;
    /// Boolean views class 1
    ViewArray<VA> x;
    /// lower bound class 1 (part from class1 views assigned to 1)
    int l_x;
    /// total class 2
    int Y;
    /// Boolean views class y
    ViewArray<VA> y;
    /// lower bound class 2 (part from class2 views assigned to 1)
    int l_y;
    /// Righthandside, the treshold
    VC c;
    /// Precision of treshold-integer conversion (eg. 1000)
    int precision;

    /// Constructor for creation
    FconvexBool_emulatorGq(Space& home, int X, ViewArray<VA>& x, int l_x,
            int Y, ViewArray<VA>& y, int l_y, VC c, int precision);
    /// Constructor for cloning \a p
    FconvexBool_emulatorGq(Space& home, bool share, FconvexBool_emulatorGq& p);
  public:
    /// Cost function (defined as low linear)
    virtual PropCost cost(const Space& home, const ModEventDelta& med) const;
    /// Delete propagator and return its size
    virtual size_t dispose(Space& home);
    /// Create copy during cloning
    virtual Actor* copy(Space& home, bool share);
    /// Perform propagation
    virtual ExecStatus propagate(Space& home, const ModEventDelta& med);
    /// Eliminate assigned vars and update bounds
    void eliminate(ViewArray<VA>& v, int& l, int& u);
    /// Post propagator but
    static ExecStatus post(Space& home, int X, ViewArray<VA>& x, int l_x, int Y, ViewArray<VA>& y, int l_y, VC c, int precision);
  };

  /*
   * Chi2 Class
   */
  template <class VA, class VC>
  FconvexBool_emulatorGq<VA,VC>::FconvexBool_emulatorGq(Space& home, int X0,
        ViewArray<VA>& x0, int l_x0, int Y0, ViewArray<VA>& y0, int l_y0, VC c0, int precision0)
    :  Propagator(home), X(X0),x(x0),l_x(l_x0),Y(Y0),y(y0),l_y(l_y0),c(c0),precision(precision0) {
    x.subscribe(home,*this,PC_INT_VAL);
    y.subscribe(home,*this,PC_INT_VAL);
    c.subscribe(home,*this,PC_INT_BND);
  }

  template <class VA, class VC>
  forceinline size_t
  FconvexBool_emulatorGq<VA,VC>::dispose(Space& home) {
    assert(!home.failed());
    x.cancel(home,*this,PC_INT_VAL);
    y.cancel(home,*this,PC_INT_VAL);
    c.cancel(home,*this,PC_INT_BND);
    (void) Propagator::dispose(home);
    return sizeof(*this);
  }

  template <class VA, class VC>
  forceinline
  FconvexBool_emulatorGq<VA,VC>::FconvexBool_emulatorGq(Space& home,
            bool share, FconvexBool_emulatorGq& p)
    : Propagator(home,share,p), X(p.X),l_x(p.l_x),Y(p.Y),l_y(p.l_y),precision(p.precision) {
    x.update(home,share,p.x);
    y.update(home,share,p.y);
    c.update(home,share,p.c);
  }

  template <class VA, class VC>
  PropCost
  FconvexBool_emulatorGq<VA,VC>::cost(const Space&, const ModEventDelta&) const {
    return PropCost::linear(PropCost::LO, x.size());
  }

  template <class VA, class VC>
  ExecStatus
  FconvexBool_emulatorGq<VA,VC>::post(Space& home, int X, ViewArray<VA>& x, int l_x,
                                int Y, ViewArray<VA>& y, int l_y, VC c, int precision) {

    (void) new (home) FconvexBool_emulatorGq<VA,VC>(home,X,x,l_x,Y,y,l_y,c,precision);
    return ES_OK;
  }

  template <class VA, class VC>
  Actor*
  FconvexBool_emulatorGq<VA,VC>::copy(Space& home, bool share) {
    return new (home) FconvexBool_emulatorGq<VA,VC>(home,share,*this);
  }

  template <class VA, class VC>
  ExecStatus
  FconvexBool_emulatorGq<VA,VC>::propagate(Space& home, const ModEventDelta&) {
    int u_x = 0; int u_y = 0;
    eliminate(x, l_x, u_x);
    eliminate(y, l_y, u_y);
    u_x += l_x;
    u_y += l_y;
    /*
    // Tias debug
    std::cerr << "X("<<X<<") "<<l_x<<" "<<u_x<<"\tY("<<Y<<") "<<l_y<<" "<<u_y;
    std::cerr <<" l_x u_y "<<l_x<<" "<<u_y<<" : "<<convex_functionWrapper(X,l_x,u_x,Y,l_y,u_y,precision)
              << " >?= " << c << std::endl;
    */
    // maximum of convex function must be greater or equal than treshold
    //if ((convex_function(X,l_x,Y,u_y,precision) < c.min())
    // && (convex_function(X,u_x,Y,l_y,precision) < c.min())) {
    if (convex_functionWrapper(X,l_x,u_x,Y,l_y,u_y,precision, c.min())) {
      return ES_FAILED;
    } else if (x.size() == 0 && y.size() == 0) {
      // FINAL fronteer, otherwise gecode goes bereserk with such a misbehaving propagator !
      if (convex_function(X,u_x,Y,u_y,precision) < c.min())
        return ES_FAILED;
      return home.ES_SUBSUMED(*this);
    }
    return ES_FIX;
  }

  template <class VA, class VC>
  void FconvexBool_emulatorGq<VA,VC>::eliminate(ViewArray<VA>& v, int& l, int& u) {
    int n = v.size();
    for (int i = n; i--; )
      if (v[i].zero()) {
        v[i]=v[--n];
      } else if (v[i].one()) {
        v[i]=v[--n];
        l += 1;
      } else {
        u += 1;
      }
    v.size(n);
  }


  /** \brief Post propagator for the convex chi2 propagator for booleans
   * Calculates and prunes based on the chi2 value
   *
   * The IntArgs represent the class and must be positive unit (0 or 1)
   * Returns Error if not positive unit.
   *
   */
  void convex(Space& home,
                const int classX_tot, const BoolVarArgs& classX,
                const int classY_tot, const BoolVarArgs& classY,
                IntRelType r, IntVar c, int precision) {
    if (home.failed()) return;

    // check float limit ?

    ViewArray<BoolView> classXv(home, classX);
    ViewArray<BoolView> classYv(home, classY);

    // Tias debug
    //std::cout << "gonna post: " << std::endl
    //          << xv << yv << std::endl
    //          << " and it better be good" << std::endl;
    // post
    if (r == IRT_GQ) {
        GECODE_ES_FAIL((FconvexBool_emulatorGq<BoolView,IntView>
                             ::post(home, classX_tot, classXv, 0, classY_tot, classYv, 0, c, precision)));
    } else {
      throw UnknownRelation("FconvexBool_emulator::convex");
    }
  }





  /************************************
   *      REIFY-IMPLIED               *
   ***********************************/

  /**
   * \brief Base-class for reify-implied convex boolean propagators
   *
   * This class NEEDS a global function:
   * convex_function(int pos_total, int pos, int neg_total, int neg, int precision)
   */
  template <class VA, class VC>
  class ImplyFconvexBool_emulatorGq : public FconvexBool_emulatorGq<VA,VC> {
  protected:
    using FconvexBool_emulatorGq<VA,VC>::X;
    using FconvexBool_emulatorGq<VA,VC>::x;
    using FconvexBool_emulatorGq<VA,VC>::l_x;
    using FconvexBool_emulatorGq<VA,VC>::Y;
    using FconvexBool_emulatorGq<VA,VC>::y;
    using FconvexBool_emulatorGq<VA,VC>::l_y;
    using FconvexBool_emulatorGq<VA,VC>::c;
    using FconvexBool_emulatorGq<VA,VC>::precision;
    /// reify-implication variable
    BoolView b;

    ImplyFconvexBool_emulatorGq(Space& home, BoolView b, int X, ViewArray<VA>& x, int l_x,
            int Y, ViewArray<VA>& y, int l_y, VC c, int precision);
    /// Constructor for cloning \a p
    ImplyFconvexBool_emulatorGq(Space& home, bool share, ImplyFconvexBool_emulatorGq& p);
  public:
    /// Delete propagator and return its size
    virtual size_t dispose(Space& home);
    /// Create copy during cloning
    virtual Actor* copy(Space& home, bool share);
    /// Perform propagation
    virtual ExecStatus propagate(Space& home, const ModEventDelta& med);
    /// Post propagator
    static ExecStatus post(Space& home, BoolView b, int X, ViewArray<VA>& x, int l_x, int Y, ViewArray<VA>& y, int l_y, VC c, int precision);
  };

  /*
   * reify-implied FconvexBool_emulatorGq Class
   */
  template <class VA, class VC>
  ImplyFconvexBool_emulatorGq<VA,VC>::ImplyFconvexBool_emulatorGq(Space& home, BoolView b0, int X0,
        ViewArray<VA>& x0, int l_x0, int Y0, ViewArray<VA>& y0, int l_y0, VC c0, int precision0)
    : FconvexBool_emulatorGq<VA,VC>(home,X0,x0,l_x0,Y0,y0,l_y0,c0,precision0), b(b0) {
    b.subscribe(home,*this,PC_INT_VAL);
  }

  template <class VA, class VC>
  forceinline
  ImplyFconvexBool_emulatorGq<VA,VC>::ImplyFconvexBool_emulatorGq(Space& home,
            bool share, ImplyFconvexBool_emulatorGq& p)
    : FconvexBool_emulatorGq<VA,VC>(home, share, p) {
    b.update(home,share,p.b);
  }

  template <class VA, class VC>
  forceinline size_t
  ImplyFconvexBool_emulatorGq<VA,VC>::dispose(Space& home) {
    assert(!home.failed());
    b.cancel(home,*this,PC_INT_VAL);
    return  FconvexBool_emulatorGq<VA,VC>::dispose(home);
  }

  template <class VA, class VC>
  Actor*
  ImplyFconvexBool_emulatorGq<VA,VC>::copy(Space& home, bool share) {
    return new (home) ImplyFconvexBool_emulatorGq<VA,VC>(home,share,*this);
  }

  template <class VA, class VC>
  ExecStatus
  ImplyFconvexBool_emulatorGq<VA,VC>::post(Space& home, BoolView b, int X, ViewArray<VA>& x, int l_x, int Y, ViewArray<VA>& y, int l_y, VC c, int precision) {

    (void) new (home) ImplyFconvexBool_emulatorGq<VA,VC>(home,b,X,x,l_x,Y,y,l_y,c,precision);
    return ES_OK;
  }

  template <class VA, class VC>
  ExecStatus
  ImplyFconvexBool_emulatorGq<VA,VC>::propagate(Space& home, const ModEventDelta&) {
    if (b.one())
        GECODE_REWRITE(*this,(FconvexBool_emulatorGq<VA,VC>::post(home,X,x,l_x,Y,y,l_y,c,precision)));
    if (b.zero())
        return home.ES_SUBSUMED(*this); // couldn't care less

    int u_x = 0; int u_y = 0;
    eliminate(x, l_x, u_x);
    eliminate(y, l_y, u_y);
    u_x += l_x;
    u_y += l_y;
    /*
    // Tias debug
    std::cerr << "X("<<X<<") "<<l_x<<" "<<u_x<<"\tY("<<Y<<") "<<l_y<<" "<<u_y;
    std::cerr <<" l_x u_y "<<l_x<<" "<<u_y<<" : "<<convex_functionWrapper(X,l_x,u_x,Y,l_y,u_y,precision)
              << " >?= " << c << std::endl;
    */
    // maximum of convex function must be greater or equal than treshold
    //if ((convex_function(X,l_x,Y,u_y,precision) < c.min())
    // && (convex_function(X,u_x,Y,l_y,precision) < c.min())) {
    if (convex_functionWrapper(X,l_x,u_x,Y,l_y,u_y,precision, c.min())) {
      GECODE_ME_CHECK(b.zero_none(home));
      return home.ES_SUBSUMED(*this);
    } else if (x.size() == 0 && y.size() == 0) {
      return home.ES_SUBSUMED(*this); // couldn't care less
    }
    return ES_FIX;
  }

  /** \brief Post propagator for the convex chi2 propagator for booleans, reified version
   * Calculates and prunes based on the chi2 value
   *
   * The IntArgs represent the class and must be positive unit (0 or 1)
   * Returns Error if not positive unit.
   *
   */
  void imply_convex(Space& home, const BoolVar b,
                const int classX_tot, const BoolVarArgs& classX,
                const int classY_tot, const BoolVarArgs& classY,
                IntRelType r, IntVar c, int precision) {
    if (home.failed()) return;

    if (b.one()) {
        return convex(home, classX_tot, classX, classY_tot, classY, r, c, precision);
    }
    if (b.zero())
        return; // couldn't care less


    // check float limit ?

    // Only accept IRT_GQ (for now)
    if (r != IRT_GQ)
      throw UnknownRelation("FconvexBool_emulator::imply_convex");

    ViewArray<BoolView> classXv(home, classX);
    ViewArray<BoolView> classYv(home, classY);

    // Tias debug
    //std::cout << "gonna post: " << std::endl
    //          << xv << yv << std::endl
    //          << " and it better be good" << std::endl;
    // post
    GECODE_ES_FAIL((ImplyFconvexBool_emulatorGq<BoolView,IntView>
                             ::post(home, b, classX_tot, classXv, 0,
                                             classY_tot, classYv, 0, c, precision)));
  }



}
