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

namespace FconvexBool {
  using namespace ::Gecode;
  using namespace ::Gecode::Int;

  /**
   * \brief Base-class for convex function boolean propagators
   *
   * This class NEEDS a global function:
   * convex_function(int pos_total, int pos, int neg_total, int neg, int precision)
   */
  template <class VA, class VC>
  class FconvexBoolGq : public Propagator {
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
    FconvexBoolGq(Space& home, int X, ViewArray<VA>& x, int l_x,
            int Y, ViewArray<VA>& y, int l_y, VC c, int precision);
    /// Constructor for cloning \a p
    FconvexBoolGq(Space& home, bool share, FconvexBoolGq& p);
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
   * FconvexBoolGq Class
   */
  template <class VA, class VC>
  FconvexBoolGq<VA,VC>::FconvexBoolGq(Space& home, int X0,
        ViewArray<VA>& x0, int l_x0, int Y0, ViewArray<VA>& y0, int l_y0, VC c0, int precision0)
    :  Propagator(home), X(X0),x(x0),l_x(l_x0),Y(Y0),y(y0),l_y(l_y0),c(c0),precision(precision0) {
    x.subscribe(home,*this,PC_INT_VAL);
    y.subscribe(home,*this,PC_INT_VAL);
    c.subscribe(home,*this,PC_INT_BND);
  }

  template <class VA, class VC>
  forceinline size_t
  FconvexBoolGq<VA,VC>::dispose(Space& home) {
    assert(!home.failed());
    x.cancel(home,*this,PC_INT_VAL);
    y.cancel(home,*this,PC_INT_VAL);
    c.cancel(home,*this,PC_INT_BND);
    (void) Propagator::dispose(home);
    return sizeof(*this);
  }

  template <class VA, class VC>
  forceinline
  FconvexBoolGq<VA,VC>::FconvexBoolGq(Space& home,
            bool share, FconvexBoolGq& p)
    : Propagator(home,share,p), X(p.X),l_x(p.l_x),Y(p.Y),l_y(p.l_y),precision(p.precision) {
    x.update(home,share,p.x);
    y.update(home,share,p.y);
    c.update(home,share,p.c);
  }

  template <class VA, class VC>
  PropCost
  FconvexBoolGq<VA,VC>::cost(const Space&, const ModEventDelta&) const {
    return PropCost::linear(PropCost::LO, x.size());
  }

  template <class VA, class VC>
  ExecStatus
  FconvexBoolGq<VA,VC>::post(Space& home, int X, ViewArray<VA>& x, int l_x,
                                int Y, ViewArray<VA>& y, int l_y, VC c, int precision) {

    (void) new (home) FconvexBoolGq<VA,VC>(home,X,x,l_x,Y,y,l_y,c,precision);
    return ES_OK;
  }

  template <class VA, class VC>
  Actor*
  FconvexBoolGq<VA,VC>::copy(Space& home, bool share) {
    return new (home) FconvexBoolGq<VA,VC>(home,share,*this);
  }

  template <class VA, class VC>
  ExecStatus
  FconvexBoolGq<VA,VC>::propagate(Space& home, const ModEventDelta&) {
    int u_x = 0; int u_y = 0;
    this->eliminate(x, l_x, u_x);
    this->eliminate(y, l_y, u_y);
    u_x += l_x;
    u_y += l_y;
    /*
    // Tias debug
    std::cerr << "X("<<X<<") "<<l_x<<" "<<u_x<<"\tY("<<Y<<") "<<l_y<<" "<<u_y;
    std::cerr <<" l_x u_y "<<l_x<<" "<<u_y<<" : "<<convex_function(X,l_x,Y,u_y,precision)
              << "u_x l_y "<<u_x<<" "<<l_y<<" : "<<convex_function(X,u_x,Y,l_y,precision)
              << " >?= " << c << std::endl;
    */
    // maximum of convex function must be greater or equal than treshold
    if ((convex_function(X,l_x,Y,u_y,precision) < c.min())
     && (convex_function(X,u_x,Y,l_y,precision) < c.min())) {
      return ES_FAILED;
    } else if (x.size() == 0 && y.size() == 0) {
      return home.ES_SUBSUMED(*this);
    }
    return ES_FIX;
  }

  template <class VA, class VC>
  void FconvexBoolGq<VA,VC>::eliminate(ViewArray<VA>& v, int& l, int& u) {
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


  /**
   * \brief Base-class for reify-implied convex boolean propagators
   *
   * This class NEEDS a global function:
   * convex_function(int pos_total, int pos, int neg_total, int neg, int precision)
   *
   * Warning: when calculating the minimum of the convex function,
   * we presume the diagonal between the two classes is the global minimum of 0
   */
  template <class VA, class VC>
  class FconvexBoolLq : public FconvexBoolGq<VA,VC> {
  protected:
    using FconvexBoolGq<VA,VC>::X;
    using FconvexBoolGq<VA,VC>::x;
    using FconvexBoolGq<VA,VC>::l_x;
    using FconvexBoolGq<VA,VC>::Y;
    using FconvexBoolGq<VA,VC>::y;
    using FconvexBoolGq<VA,VC>::l_y;
    using FconvexBoolGq<VA,VC>::c;
    using FconvexBoolGq<VA,VC>::precision;

    FconvexBoolLq(Space& home, int X, ViewArray<VA>& x, int l_x,
            int Y, ViewArray<VA>& y, int l_y, VC c, int precision);
    /// Constructor for cloning \a p
    FconvexBoolLq(Space& home, bool share, FconvexBoolLq& p);
  public:
    /// Create copy during cloning
    virtual Actor* copy(Space& home, bool share);
    /// Perform propagation
    virtual ExecStatus propagate(Space& home, const ModEventDelta& med);
    /// Post propagator
    static ExecStatus post(Space& home, int X, ViewArray<VA>& x, int l_x, int Y, ViewArray<VA>& y, int l_y, VC c, int precision);
  };

  /*
   * FconvexBoolLq Class
   */
  template <class VA, class VC>
  FconvexBoolLq<VA,VC>::FconvexBoolLq(Space& home, int X0,
        ViewArray<VA>& x0, int l_x0, int Y0, ViewArray<VA>& y0, int l_y0, VC c0, int precision0)
    : FconvexBoolGq<VA,VC>(home,X0,x0,l_x0,Y0,y0,l_y0,c0,precision0) { }

  template <class VA, class VC>
  forceinline
  FconvexBoolLq<VA,VC>::FconvexBoolLq(Space& home,
            bool share, FconvexBoolLq& p)
    : FconvexBoolGq<VA,VC>(home, share, p) { }

  template <class VA, class VC>
  Actor*
  FconvexBoolLq<VA,VC>::copy(Space& home, bool share) {
    return new (home) FconvexBoolLq<VA,VC>(home,share,*this);
  }

  template <class VA, class VC>
  ExecStatus
  FconvexBoolLq<VA,VC>::post(Space& home, int X, ViewArray<VA>& x, int l_x, int Y, ViewArray<VA>& y, int l_y, VC c, int precision) {

    (void) new (home) FconvexBoolLq<VA,VC>(home,X,x,l_x,Y,y,l_y,c,precision);
    return ES_OK;
  }

  template <class VA, class VC>
  ExecStatus
  FconvexBoolLq<VA,VC>::propagate(Space& home, const ModEventDelta&) {
    int u_x = 0; int u_y = 0;
    this->eliminate(x, l_x, u_x);
    this->eliminate(y, l_y, u_y);
    u_x += l_x;
    u_y += l_y;
    // Tias debug
    /*
    std::cerr << "X("<<X<<") "<<l_x<<" "<<u_x<<"\tY("<<Y<<") "<<l_y<<" "<<u_y;
    std::cerr << "\tl_x l_y "<<l_x<<" "<<l_y<<" : "<<convex_function(X,l_x,Y,l_y,precision)
              << "\tu_x u_y "<<u_x<<" "<<u_y<<" : "<<convex_function(X,u_x,Y,u_y,precision)
              <<" \tl_x u_y "<<l_x<<" "<<u_y<<" : "<<convex_function(X,l_x,Y,u_y,precision)
              << "\tu_x l_y "<<u_x<<" "<<l_y<<" : "<<convex_function(X,u_x,Y,l_y,precision)
              << "\t >?= " << c << std::endl;
    */
    float zeropoint = (float)X/Y;
    // minimum of convex function must be lower or equal than treshold
    if (((float)l_x/u_y) > zeropoint && ((float)u_x/l_y) > zeropoint) {
        // both at same side, lowest point is l_x,u_y
        if (convex_function(X,l_x,Y,u_y,precision) > c.min()) {
            return ES_FAILED;
        }
    } else if (((float)l_x/u_y) < zeropoint && ((float)u_x/l_y) < zeropoint) {
        // both at same side, lowest point is u_x,l_y
        if (convex_function(X,u_x,Y,l_y,precision) > c.min()) {
            return ES_FAILED;
        }
    //} else {
        // possible area crosses zeropoint-line
        // minimum will always be LQ tresh (which is positive)
    }
        
    if (x.size() == 0 && y.size() == 0) {
      return home.ES_SUBSUMED(*this);
    }
    return ES_FIX;
  }


  /** \brief Post propagator for the convex function propagator for booleans
   * Calculates and prunes based on a convex function value
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
        GECODE_ES_FAIL((FconvexBoolGq<BoolView,IntView>
                             ::post(home, classX_tot, classXv, 0, classY_tot, classYv, 0, c, precision)));
    } else if (r == IRT_LQ) {
        GECODE_ES_FAIL((FconvexBoolLq<BoolView,IntView>
                             ::post(home, classX_tot, classXv, 0, classY_tot, classYv, 0, c, precision)));
    } else {
      throw UnknownRelation("FconvexBool::convex");
    }
  }




  /************************************
   *      REIFIED                     *
   ***********************************/

  /**
   * \brief Base-class for reified convex function boolean propagators
   *
   * This class NEEDS a global function:
   * convex_function(int pos_total, int pos, int neg_total, int neg, int precision)
   *
   * Warning: when calculating the minimum of the convex function,
   * we presume the diagonal between the two classes is the global minimum of 0
   */
  template <class VA, class VC>
  class ReifiedFconvexBoolGq : public FconvexBoolGq<VA,VC> {
  protected:
    using FconvexBoolGq<VA,VC>::X;
    using FconvexBoolGq<VA,VC>::x;
    using FconvexBoolGq<VA,VC>::l_x;
    using FconvexBoolGq<VA,VC>::Y;
    using FconvexBoolGq<VA,VC>::y;
    using FconvexBoolGq<VA,VC>::l_y;
    using FconvexBoolGq<VA,VC>::c;
    using FconvexBoolGq<VA,VC>::precision;
    /// reify-implication variable
    BoolView b;

    ReifiedFconvexBoolGq(Space& home, BoolView b, int X, ViewArray<VA>& x, int l_x,
            int Y, ViewArray<VA>& y, int l_y, VC c, int precision);
    /// Constructor for cloning \a p
    ReifiedFconvexBoolGq(Space& home, bool share, ReifiedFconvexBoolGq& p);
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
   * FconvexBoolGq Class reify-implied
   */
  template <class VA, class VC>
  ReifiedFconvexBoolGq<VA,VC>::ReifiedFconvexBoolGq(Space& home, BoolView b0, int X0,
        ViewArray<VA>& x0, int l_x0, int Y0, ViewArray<VA>& y0, int l_y0, VC c0, int precision0)
    : FconvexBoolGq<VA,VC>(home,X0,x0,l_x0,Y0,y0,l_y0,c0,precision0), b(b0) {
    b.subscribe(home,*this,PC_INT_VAL);
  }

  template <class VA, class VC>
  forceinline
  ReifiedFconvexBoolGq<VA,VC>::ReifiedFconvexBoolGq(Space& home,
            bool share, ReifiedFconvexBoolGq& p)
    : FconvexBoolGq<VA,VC>(home, share, p) {
    b.update(home,share,p.b);
  }

  template <class VA, class VC>
  forceinline size_t
  ReifiedFconvexBoolGq<VA,VC>::dispose(Space& home) {
    assert(!home.failed());
    b.cancel(home,*this,PC_INT_VAL);
    return  FconvexBoolGq<VA,VC>::dispose(home);
  }

  template <class VA, class VC>
  Actor*
  ReifiedFconvexBoolGq<VA,VC>::copy(Space& home, bool share) {
    return new (home) ReifiedFconvexBoolGq<VA,VC>(home,share,*this);
  }

  template <class VA, class VC>
  ExecStatus
  ReifiedFconvexBoolGq<VA,VC>::post(Space& home, BoolView b, int X, ViewArray<VA>& x, int l_x, int Y, ViewArray<VA>& y, int l_y, VC c, int precision) {

    (void) new (home) ReifiedFconvexBoolGq<VA,VC>(home,b,X,x,l_x,Y,y,l_y,c,precision);
    return ES_OK;
  }

  template <class VA, class VC>
  ExecStatus
  ReifiedFconvexBoolGq<VA,VC>::propagate(Space& home, const ModEventDelta&) {
    if (b.one()) {
        GECODE_REWRITE(*this,(FconvexBoolGq<VA,VC>::post(home,X,x,l_x,Y,y,l_y,c,precision)));
    }
    if (b.zero()) {
        OffsetView cm1(c, -1);
        GECODE_REWRITE(*this,(FconvexBoolLq<VA,OffsetView>::post(home,X,x,l_x,Y,y,l_y,cm1,precision)));
    }

    int u_x = 0; int u_y = 0;
    this->eliminate(x, l_x, u_x);
    this->eliminate(y, l_y, u_y);
    u_x += l_x;
    u_y += l_y;
    /*
    // Tias debug
    std::cerr << "X("<<X<<") "<<l_x<<" "<<u_x<<"\tY("<<Y<<") "<<l_y<<" "<<u_y;
    std::cerr <<" l_x l_y "<<l_x<<" "<<l_y<<" : "<<convex_function(X,l_x,Y,l_y,precision)
              << "u_x u_y "<<u_x<<" "<<u_y<<" : "<<convex_function(X,u_x,Y,u_y,precision)
              << " >?= " << c << std::endl;
    */
    float zeropoint = (float)X/Y;
    // maximum of convex function may not be smaller than treshold
    if ((convex_function(X,l_x,Y,u_y,precision) < c.min())
     && (convex_function(X,u_x,Y,l_y,precision) < c.min())) {
      GECODE_ME_CHECK(b.zero_none(home));
      return home.ES_SUBSUMED(*this);
    }
    // minimum of convex function should be higher or equal than treshold
    else if (((float)l_x/u_y) > zeropoint && ((float)u_x/l_y) > zeropoint) {
        // both at bottom right side, lowest point is l_x,u_y
        if (convex_function(X,l_x,Y,u_y,precision) >= c.min()) {
            GECODE_ME_CHECK(b.one_none(home));
            return home.ES_SUBSUMED(*this);
        }
    } else if (((float)l_x/u_y) < zeropoint && ((float)u_x/l_y) < zeropoint) {
        // both at upper left, lowest point is u_x,l_y
        if (convex_function(X,u_x,Y,l_y,precision) > c.min()) {
            GECODE_ME_CHECK(b.one_none(home));
            return home.ES_SUBSUMED(*this);
        }
    }
    else if (x.size() == 0 && y.size() == 0) {
      // all set and minimum is not higher: fail
      GECODE_ME_CHECK(b.zero_none(home));
      return home.ES_SUBSUMED(*this);
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
  void reify_convex(Space& home, const BoolVar b,
                const int classX_tot, const BoolVarArgs& classX,
                const int classY_tot, const BoolVarArgs& classY,
                IntRelType r, IntVar c, int precision) {
    if (home.failed()) return;

    // check float limit ?

    // Only accept IRT_GQ (for now)
    if (r != IRT_GQ)
      throw UnknownRelation("FconvexBool::imply_convex");

    ViewArray<BoolView> classXv(home, classX);
    ViewArray<BoolView> classYv(home, classY);

    // Tias debug
    //std::cout << "gonna post: " << std::endl
    //          << xv << yv << std::endl
    //          << " and it better be good" << std::endl;
    // post
    GECODE_ES_FAIL((ReifiedFconvexBoolGq<BoolView,IntView>
                             ::post(home, b, classX_tot, classXv, 0,
                                             classY_tot, classYv, 0, c, precision)));
  }

  /**
   * \brief Base-class for reify-implied convex function boolean propagators
   *
   * This class NEEDS a global function:
   * convex_function(int pos_total, int pos, int neg_total, int neg, int precision)
   */
  template <class VA, class VC>
  class ImplyFconvexBoolGq : public FconvexBoolGq<VA,VC> {
  protected:
    using FconvexBoolGq<VA,VC>::X;
    using FconvexBoolGq<VA,VC>::x;
    using FconvexBoolGq<VA,VC>::l_x;
    using FconvexBoolGq<VA,VC>::Y;
    using FconvexBoolGq<VA,VC>::y;
    using FconvexBoolGq<VA,VC>::l_y;
    using FconvexBoolGq<VA,VC>::c;
    using FconvexBoolGq<VA,VC>::precision;
    /// reify-implication variable
    BoolView b;

    ImplyFconvexBoolGq(Space& home, BoolView b, int X, ViewArray<VA>& x, int l_x,
            int Y, ViewArray<VA>& y, int l_y, VC c, int precision);
    /// Constructor for cloning \a p
    ImplyFconvexBoolGq(Space& home, bool share, ImplyFconvexBoolGq& p);
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
   * reify-implied FconvexBoolGq Class
   */
  template <class VA, class VC>
  ImplyFconvexBoolGq<VA,VC>::ImplyFconvexBoolGq(Space& home, BoolView b0, int X0,
        ViewArray<VA>& x0, int l_x0, int Y0, ViewArray<VA>& y0, int l_y0, VC c0, int precision0)
    : FconvexBoolGq<VA,VC>(home,X0,x0,l_x0,Y0,y0,l_y0,c0,precision0), b(b0) {
    b.subscribe(home,*this,PC_INT_VAL);
  }

  template <class VA, class VC>
  forceinline
  ImplyFconvexBoolGq<VA,VC>::ImplyFconvexBoolGq(Space& home,
            bool share, ImplyFconvexBoolGq& p)
    : FconvexBoolGq<VA,VC>(home, share, p) {
    b.update(home,share,p.b);
  }

  template <class VA, class VC>
  forceinline size_t
  ImplyFconvexBoolGq<VA,VC>::dispose(Space& home) {
    assert(!home.failed());
    b.cancel(home,*this,PC_INT_VAL);
    return  FconvexBoolGq<VA,VC>::dispose(home);
  }

  template <class VA, class VC>
  Actor*
  ImplyFconvexBoolGq<VA,VC>::copy(Space& home, bool share) {
    return new (home) ImplyFconvexBoolGq<VA,VC>(home,share,*this);
  }

  template <class VA, class VC>
  ExecStatus
  ImplyFconvexBoolGq<VA,VC>::post(Space& home, BoolView b, int X, ViewArray<VA>& x, int l_x, int Y, ViewArray<VA>& y, int l_y, VC c, int precision) {

    (void) new (home) ImplyFconvexBoolGq<VA,VC>(home,b,X,x,l_x,Y,y,l_y,c,precision);
    return ES_OK;
  }

  template <class VA, class VC>
  ExecStatus
  ImplyFconvexBoolGq<VA,VC>::propagate(Space& home, const ModEventDelta&) {
    if (b.one())
        GECODE_REWRITE(*this,(FconvexBoolGq<VA,VC>::post(home,X,x,l_x,Y,y,l_y,c,precision)));
    if (b.zero())
        return home.ES_SUBSUMED(*this); // couldn't care less

    int u_x = 0; int u_y = 0;
    this->eliminate(x, l_x, u_x);
    this->eliminate(y, l_y, u_y);
    u_x += l_x;
    u_y += l_y;
    /*
    // Tias debug
    std::cerr << "X("<<X<<") "<<l_x<<" "<<u_x<<"\tY("<<Y<<") "<<l_y<<" "<<u_y;
    std::cerr <<" l_x u_y "<<l_x<<" "<<u_y<<" : "<<convex_function(X,l_x,Y,u_y,precision)
              << "u_x l_y "<<u_x<<" "<<l_y<<" : "<<convex_function(X,u_x,Y,l_y,precision)
              << " >?= " << c << std::endl;
    */
    // maximum of convex function must be greater or equal than treshold
    if ((convex_function(X,l_x,Y,u_y,precision) < c.min())
     && (convex_function(X,u_x,Y,l_y,precision) < c.min())) {
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
      throw UnknownRelation("FconvexBool::imply_convex");

    ViewArray<BoolView> classXv(home, classX);
    ViewArray<BoolView> classYv(home, classY);

    // Tias debug
    //std::cout << "gonna post: " << std::endl
    //          << xv << yv << std::endl
    //          << " and it better be good" << std::endl;
    // post
    GECODE_ES_FAIL((ImplyFconvexBoolGq<BoolView,IntView>
                             ::post(home, b, classX_tot, classXv, 0,
                                             classY_tot, classYv, 0, c, precision)));
  }



}
