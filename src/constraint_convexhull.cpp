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

namespace ConvexHull {
  using namespace ::Gecode;
  using namespace ::Gecode::Int;

  class HullPoint {
  public:
    // Next pointer
    HullPoint* next;
    // Previous pointer
    HullPoint* prev;
    // Items of this point
    const vector<bool> items;
    // Transactions of this point
    const vector<bool> trans;
    // Nr transactions of class X
    const int X;
    // Nr transactions of class Y
    const int Y;
    // is it a fake point (eg fake first or last)
    const bool fake;

    // Constructor
    HullPoint(int X0, int Y0):
        X(X0), Y(Y0), fake(true) {
        next = NULL;
        prev = NULL;
    }
    HullPoint(vector<bool> items0, vector<bool> trans0, int X0, int Y0):
        items(items0), trans(trans0), X(X0), Y(Y0), fake(false) {
        next = NULL;
        prev = NULL;
    }

    // test if a point is Left|On|Right of an infinite point
    // ref: http://softsurfer.com/Archive/algorithm_0203/algorithm_0203.htm
    // >0 if left of the line through P0 and P1
    // =0 if on the line
    // <0 if right of the line
    inline float leftorright(HullPoint* P0, HullPoint* P1) {
        return (P1->X - P0->X)*(Y - P0->Y) - (X - P0->X)*(P1->Y - P0->Y);
    }

    // check wheter point is left or on the line throught P0 and P1
    inline bool isLeftorOn(HullPoint* P0, HullPoint* P1) {
        return (leftorright(P0,P1) >= 0);
    }
    // check wheter point is right of the line throught P0 and P1
    inline bool isRightorOn(HullPoint* P0, HullPoint* P1) {
        return (leftorright(P0,P1) <= 0);
    }

    inline void print(FILE* f) {
        if (!fake) {
            fprintf(f, "\t(%i,%i) ", X, Y);
            for (unsigned int i=0; i!=items.size(); i++) {
                if (items[i])
                    fprintf(f, "%i ", i);
            }
            fprintf(f, "\n");
        }
    }

    inline void print_full(FILE* f) {
        if (!fake) {
            fprintf(f, "\t(%i,%i) ", X, Y);
            for (unsigned int i=0; i!=items.size(); i++) {
                if (items[i])
                    fprintf(f, "%i ", i);
            }
            fprintf(f, " < ");
            for (unsigned int t=0; t!=trans.size(); t++) {
                if (trans[t])
                    fprintf(f, "%i ", t);
            }
            fprintf(f, ">\n");
        }
    }
  };

  void addToHullUL(HullPoint* point, HullPoint* hull) {    
    //std::cerr << "Adding point UL: "; point->print(stderr);

    // find prev and next of point
    HullPoint* next = hull->next;
    while (next != NULL && next->X < point->X)
        next = next->next;

    HullPoint* prev = next->prev;
    // attach point between prev & next
    prev->next = point;
    point->prev = prev;
    point->next = next;
    next->prev = point;

    // Shortcut: next and point on same (X,Y)
    if ((point->X == next->X) && (point->Y == next->Y))
        return;

    // remove no-longer convex points in direction of prev
    next = point->prev;
    prev = next->prev;
    while (prev != NULL && point->isLeftorOn(prev,next)) {
        //std::cerr << "Doing delete of next\n"; next->print(stderr);
        // next is obsolete: skip it and delete
        prev->next = point;
        point->prev = prev;
        delete next;
        next = prev;
        prev = next->prev;
    }

    // remove no-longer convex points in direction of next
    prev = point->next;
    next = prev->next;
    while (next != NULL && point->isRightorOn(next,prev)) {
        //std::cerr << "Doing delete of prev: "; prev->print(stderr);
        // prev is obsolete: skip it and delete
        point->next = next;
        next->prev = point;
        delete prev;
        prev = next;
        next = prev->next;
    }

    return;
  }
  void addToHullLR(HullPoint* point, HullPoint* hull) {    
    //std::cerr << "Adding point LR: "; point->print(stderr);

    // find prev and next of point
    HullPoint* next = hull->next;
    while (next != NULL && next->Y < point->Y)
        next = next->next;

    HullPoint* prev = next->prev;
    // attach point between prev & next
    prev->next = point;
    point->prev = prev;
    point->next = next;
    next->prev = point;

    // Shortcut: next and point on same (X,Y)
    if ((point->X == next->X) && (point->Y == next->Y))
        return;

    // remove no-longer convex points in direction of prev
    next = point->prev;
    prev = next->prev;
    while (prev != NULL && point->isRightorOn(prev,next)) {
        //std::cerr << "Doing delete of next\n"; next->print(stderr);
        // next is obsolete: skip it and delete
        prev->next = point;
        point->prev = prev;
        delete next;
        next = prev;
        prev = next->prev;
    }

    // remove no-longer convex points in direction of next
    prev = point->next;
    next = prev->next;
    while (next != NULL && point->isLeftorOn(next,prev)) {
        //std::cerr << "Doing delete of prev: "; prev->print(stderr);
        // prev is obsolete: skip it and delete
        point->next = next;
        next->prev = point;
        delete prev;
        prev = next;
        next = prev->next;
    }

    return;
  }

  /**
   * \brief Base-class for convex boolean propagator using convex hull
   */
  template <class VA>
  class BoolConvexHull : public Propagator {
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
    // start of Upper Left hull
    HullPoint* hullUL;
    // start of Lower Right hull
    HullPoint* hullLR;

    /// Constructor for creation
    BoolConvexHull(Space& home, int X, ViewArray<VA>& x, int l_x,
            int Y, ViewArray<VA>& y, int l_y, HullPoint* hullUL, HullPoint* HullLR);
    /// Constructor for cloning \a p
    BoolConvexHull(Space& home, bool share, BoolConvexHull& p);
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
    /// Post propagator
    static ExecStatus post(Space& home, int X, ViewArray<VA>& x, int l_x,
            int Y, ViewArray<VA>& y, int l_y, HullPoint* hullUL, HullPoint* HullLR);
    // test if a point is Left|On|Right of an infinite point
    // ref: http://softsurfer.com/Archive/algorithm_0203/algorithm_0203.htm
    // >0 if left of the line through P0 and P1
    // =0 if on the line
    // <0 if right of the line
    inline float leftorright(int X, int Y, HullPoint* P0, HullPoint* P1) {
        return (P1->X - P0->X)*(Y - P0->Y) - (X - P0->X)*(P1->Y - P0->Y);
    }
    /// check wheter a point is inside the convex hull
    inline bool insideConvexHullUL(int X, int Y, HullPoint* hull);
    inline bool insideConvexHullLR(int X, int Y, HullPoint* hull);

  };

  /*
   * BoolConvexHull
   */
  template <class VA>
  BoolConvexHull<VA>::BoolConvexHull(Space& home, int X0,
        ViewArray<VA>& x0, int l_x0, int Y0, ViewArray<VA>& y0, int l_y0, HullPoint* hullUL0, HullPoint* hullLR0)
    :  Propagator(home), X(X0),x(x0),l_x(l_x0),Y(Y0),y(y0),l_y(l_y0),hullUL(hullUL0),hullLR(hullLR0) {
    x.subscribe(home,*this,PC_INT_VAL);
    y.subscribe(home,*this,PC_INT_VAL);
  }

  template <class VA>
  forceinline size_t
  BoolConvexHull<VA>::dispose(Space& home) {
    assert(!home.failed());
    x.cancel(home,*this,PC_INT_VAL);
    y.cancel(home,*this,PC_INT_VAL);
    (void) Propagator::dispose(home);
    return sizeof(*this);
  }

  template <class VA>
  forceinline
  BoolConvexHull<VA>::BoolConvexHull(Space& home,
            bool share, BoolConvexHull& p)
    : Propagator(home,share,p), X(p.X),l_x(p.l_x),Y(p.Y),l_y(p.l_y),hullUL(p.hullUL),hullLR(p.hullLR) {
    x.update(home,share,p.x);
    y.update(home,share,p.y);
  }

  template <class VA>
  PropCost
  BoolConvexHull<VA>::cost(const Space&, const ModEventDelta&) const {
    return PropCost::linear(PropCost::LO, x.size());
  }

  template <class VA>
  ExecStatus
  BoolConvexHull<VA>::post(Space& home, int X, ViewArray<VA>& x, int l_x,
            int Y, ViewArray<VA>& y, int l_y, HullPoint* hullUL, HullPoint* hullLR) {

    (void) new (home) BoolConvexHull<VA>(home,X,x,l_x,Y,y,l_y,hullUL,hullLR);
    return ES_OK;
  }

  template <class VA>
  Actor*
  BoolConvexHull<VA>::copy(Space& home, bool share) {
    return new (home) BoolConvexHull<VA>(home,share,*this);
  }

  template <class VA>
  ExecStatus
  BoolConvexHull<VA>::propagate(Space& home, const ModEventDelta&) {
    int u_x = 0; int u_y = 0;
    eliminate(x, l_x, u_x);
    eliminate(y, l_y, u_y);
    u_x += l_x;
    u_y += l_y;

    float lor_lx_uy = this->leftorright(l_x,u_y,hullUL,hullLR);
    float lor_ux_ly = this->leftorright(u_x,l_y,hullUL,hullLR);

    if (lor_lx_uy >= 0 && lor_ux_ly <= 0) {
        // each on one side of diagonal
        if (this->insideConvexHullUL(l_x,u_y,hullUL) &&
             this->insideConvexHullLR(u_x,l_y,hullLR)) {
            // fully inside
            return ES_FAILED;
        }
    } else if (lor_lx_uy > 0 && lor_ux_ly > 0) {
        // both in UL
        if (this->insideConvexHullUL(l_x,u_y,hullUL)) {
            return ES_FAILED;
        } else {
            if (!this->insideConvexHullUL(u_x,l_y,hullUL)) {
                // both outside or on hull, always fine
                return home.ES_SUBSUMED(*this);
            }
        }
    } else { //if (lor_lx_uy < 0 && lor_ux_ly < 0) {
        // both in LR
        if (this->insideConvexHullLR(u_x,l_y,hullLR)) {
            return ES_FAILED;
        } else {
            if (!this->insideConvexHullLR(l_x,u_y,hullLR)) {
                // both outside or on hull, always fine
                return home.ES_SUBSUMED(*this);
            }
        }
    }

    if (x.size() == 0 && y.size() == 0) {
      return home.ES_SUBSUMED(*this);
    }

    return ES_FIX;
  }

  template <class VA>
  inline bool BoolConvexHull<VA>::insideConvexHullUL(int X, int Y, HullPoint* hull) {
    if (X == 0 && Y == 0) { return true; } // shortcut

    // find prev and next of point
    HullPoint* next = hull->next;

    if (X == 0) // straight line up w origin
        return (next->X == 0 && Y < next->Y);
        
    // skip smaller (on X axis) points
    while (next != NULL && X > next->X)
        next = next->next;

    if (next->prev->Y == Y && Y == next->Y) {
        // straight line (per convex def., must be to end)
        return true; // skip those
    }

    // point is inside hull ?
    return (leftorright(X,Y,next->prev,next) < 0);
  }
  template <class VA>
  inline bool BoolConvexHull<VA>::insideConvexHullLR(int X, int Y, HullPoint* hull) {
    if (X == 0 && Y == 0) { return true; } // shortcut

    // find prev and next of point
    HullPoint* next = hull->next;

    if (Y == 0) // straight line right w origin
        return (next->Y == 0 && X < next->X);

    // skip smaller (on Y axis) points
    while (next != NULL && Y > next->Y)
        next = next->next;

    if (next->prev->X == X && X == next->X) {
        // straight line (per convex def., must be to end)
        return true; // skip those
    }

    // point is inside hull ?
    return (leftorright(X,Y,next->prev,next) > 0);
  }

  template <class VA>
  inline void BoolConvexHull<VA>::eliminate(ViewArray<VA>& v, int& l, int& u) {
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
  void convex_hull(Space& home,
                const int classX_tot, const BoolVarArgs& classX,
                const int classY_tot, const BoolVarArgs& classY,
                HullPoint* hullUL, HullPoint* hullLR) {
    if (home.failed()) return;

    ViewArray<BoolView> classXv(home, classX);
    ViewArray<BoolView> classYv(home, classY);

    GECODE_ES_FAIL((BoolConvexHull<BoolView>
                            ::post(home, classX_tot, classXv, 0, classY_tot, classYv, 0, hullUL, hullLR)));
    return;
  }




  /************************************
   *      REIFIED                     *
   ***********************************/


  /**
   * \brief Base-class for reify-implied convex boolean propagators
   *
   * This class NEEDS a global function:
   * convex_function(int pos_total, int pos, int neg_total, int neg, int precision)
   */
  template <class VA>
  class ImplyBoolConvexHull : public BoolConvexHull<VA> {
  protected:
    using BoolConvexHull<VA>::X;
    using BoolConvexHull<VA>::x;
    using BoolConvexHull<VA>::l_x;
    using BoolConvexHull<VA>::Y;
    using BoolConvexHull<VA>::y;
    using BoolConvexHull<VA>::l_y;
    using BoolConvexHull<VA>::hullUL;
    using BoolConvexHull<VA>::hullLR;
    /// reify-implication variable
    BoolView b;

    ImplyBoolConvexHull(Space& home, BoolView b, int X, ViewArray<VA>& x, int l_x,
            int Y, ViewArray<VA>& y, int l_y, HullPoint* hullUL, HullPoint* hullLR);
    /// Constructor for cloning \a p
    ImplyBoolConvexHull(Space& home, bool share, ImplyBoolConvexHull& p);
  public:
    /// Delete propagator and return its size
    virtual size_t dispose(Space& home);
    /// Create copy during cloning
    virtual Actor* copy(Space& home, bool share);
    /// Perform propagation
    virtual ExecStatus propagate(Space& home, const ModEventDelta& med);
    /// Post propagator
    static ExecStatus post(Space& home, BoolView b, int X, ViewArray<VA>& x, int l_x, int Y, ViewArray<VA>& y, int l_y, HullPoint* hullUL, HullPoint* hullLR);
  };

  /*
   * reify-implied BoolConvexHull Class
   */
  template <class VA>
  ImplyBoolConvexHull<VA>::ImplyBoolConvexHull(Space& home, BoolView b0, int X0,
        ViewArray<VA>& x0, int l_x0, int Y0, ViewArray<VA>& y0, int l_y0, HullPoint* hullUL0, HullPoint* hullLR0)
    : BoolConvexHull<VA>(home,X0,x0,l_x0,Y0,y0,l_y0,hullUL0,hullLR0), b(b0) {
    b.subscribe(home,*this,PC_INT_VAL);
  }

  template <class VA>
  forceinline
  ImplyBoolConvexHull<VA>::ImplyBoolConvexHull(Space& home,
            bool share, ImplyBoolConvexHull& p)
    : BoolConvexHull<VA>(home, share, p) {
    b.update(home,share,p.b);
  }

  template <class VA>
  forceinline size_t
  ImplyBoolConvexHull<VA>::dispose(Space& home) {
    assert(!home.failed());
    b.cancel(home,*this,PC_INT_VAL);
    return  BoolConvexHull<VA>::dispose(home);
  }

  template <class VA>
  Actor*
  ImplyBoolConvexHull<VA>::copy(Space& home, bool share) {
    return new (home) ImplyBoolConvexHull<VA>(home,share,*this);
  }

  template <class VA>
  ExecStatus
  ImplyBoolConvexHull<VA>::post(Space& home, BoolView b, int X, ViewArray<VA>& x, int l_x, int Y, ViewArray<VA>& y, int l_y, HullPoint* hullUL, HullPoint* hullLR) {

    (void) new (home) ImplyBoolConvexHull<VA>(home,b,X,x,l_x,Y,y,l_y,hullUL,hullLR);
    return ES_OK;
  }

  template <class VA>
  ExecStatus
  ImplyBoolConvexHull<VA>::propagate(Space& home, const ModEventDelta&) {
    if (b.one())
        GECODE_REWRITE(*this,(BoolConvexHull<VA>::post(home,X,x,l_x,Y,y,l_y,hullUL,hullLR)));
    if (b.zero())
        return home.ES_SUBSUMED(*this); // couldn't care less

    int u_x = 0; int u_y = 0;
    this->eliminate(x, l_x, u_x);
    this->eliminate(y, l_y, u_y);
    u_x += l_x;
    u_y += l_y;

    float lor_lx_uy = this->leftorright(l_x,u_y,hullUL,hullLR);
    float lor_ux_ly = this->leftorright(u_x,l_y,hullUL,hullLR);

    if (lor_lx_uy >= 0 && lor_ux_ly <= 0) {
        // each on one side
        if (this->insideConvexHullUL(l_x,u_y,hullUL) &&
             this->insideConvexHullLR(u_x,l_y,hullLR)) {
            // fully inside
            GECODE_ME_CHECK(b.zero_none(home));
            return home.ES_SUBSUMED(*this);
        }
    } else if (lor_lx_uy > 0 && lor_ux_ly > 0) {
        // both in UL
        if (this->insideConvexHullUL(l_x,u_y,hullUL)) {
            GECODE_ME_CHECK(b.zero_none(home));
            return home.ES_SUBSUMED(*this);
        } else {
            if (!this->insideConvexHullUL(u_x,l_y,hullUL)) {
                // both outside hull, always fine
                return home.ES_SUBSUMED(*this);
            }
        }
    } else { //if (lor_lx_uy < 0 && lor_ux_ly < 0) {
        // both in LR
        if (this->insideConvexHullLR(u_x,l_y,hullLR)) {
            GECODE_ME_CHECK(b.zero_none(home));
            return home.ES_SUBSUMED(*this);
        } else {
            if (!this->insideConvexHullLR(l_x,u_y,hullLR)) {
                // both outside hull, always fine
                return home.ES_SUBSUMED(*this);
            }
        }
    }

    if (x.size() == 0 && y.size() == 0) {
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
  void imply_convex_hull(Space& home, const BoolVar b,
                const int classX_tot, const BoolVarArgs& classX,
                const int classY_tot, const BoolVarArgs& classY,
                HullPoint* hullUL, HullPoint* hullLR) {
    if (home.failed()) return;

    if (b.one()) {
        return convex_hull(home, classX_tot, classX, classY_tot, classY, hullUL, hullLR);
    }
    if (b.zero())
        return; // couldn't care less

    ViewArray<BoolView> classXv(home, classX);
    ViewArray<BoolView> classYv(home, classY);

    GECODE_ES_FAIL((ImplyBoolConvexHull<BoolView>
                            ::post(home, b, classX_tot, classXv, 0, classY_tot, classYv, 0, hullUL, hullLR)));
  }
}
