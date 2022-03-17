/*
 *  Main authors:
 *      Tias Guns <tias.guns@cs.kuleuven.be>
 *
 *  Copyright:
 *      Tias Guns, 2008
 *
 *  Revision information:
 *      $Id: Cimcp_convexhullhullReif.cpp 240 2009-06-04 15:54:51Z tias $
 *
 *  This file is part of Cimcp, Correlated Itemset Mining using
 *  Constraint Programming, and uses Gecode.
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

#include "common/fimcp_basic.hh"

#include "constraint_convexhull.cpp" // ConvexHull

using namespace ConvexHull;

/**
 * Mining all itemset on the convex hull in PN space (freq on pos and freq on neg)
 * inherits from Fimcp_basic.
 */
class Cimcp_convexhull : public Fimcp_basic {
protected:
  // start of Upper Left hull
  HullPoint* hullUL;
  // start of Lower Right hull
  HullPoint* hullLR;

public:
    Cimcp_convexhull(const Options_fimcp& opt) :
        Fimcp_basic(opt) {
        run(opt);
    }

    Cimcp_convexhull(bool share, Cimcp_convexhull& s) :
        Fimcp_basic(share, s) {
        hullUL = s.hullUL;
        hullLR = s.hullLR;
    }

    /// Perform copying during cloning
    virtual Space* copy(bool share) {
      return new Cimcp_convexhull(share,*this);
    }

    /// Add constraint for next better solution
    virtual void constrain(const Space&) {
        // constrain() is called after print()
        // we add the itemset to the hull in print()
        // because constrain() is never called for the last item
    }

    virtual void run(const Options_fimcp& opt);

    /// Print solution
    virtual void print(std::ostream&) const;

    /// Construct HullPoint from current solution (has to be const because print is)
    HullPoint* constructHullPoint() const;
};

void Cimcp_convexhull::run(const Options_fimcp& opt) {
    const vector< vector<bool> > tdb = common_construction(opt);
    if (classes.size() == 0)
        throw Exception("Class label error", "no class labels found");

    // count class distr
    int posTot = 0;
    for (int t = 0; t!=nr_t;t++ ) {
        posTot += classes[t];
    }
    int negTot = nr_t - posTot;

    /** Set up convex hull **/
    hullUL = new HullPoint(0,0);
    hullUL->next = new HullPoint(negTot,posTot);
    hullUL->next->prev = hullUL;
    hullUL->prev = hullUL->next; // single backlink to end
    hullLR = new HullPoint(0,0);
    hullLR->next = new HullPoint(negTot,posTot);
    hullLR->next->prev = hullLR;
    hullLR->prev = hullLR->next; // single backlink to end


    /** covered constraints **/
    if (opt.cclause()) {
        // Default! Clause is a bit faster and uses less memory
        coverage_clause(tdb);
    } else {
        IntArgs row_(nr_i);
        for (int t=0; t!=nr_t; t++) {
            // make 1-row
            for (int i=0; i!=nr_i; i++)
                row_[i] = (1-tdb[t][i]);

            // coverage: the trans its complement has no supported items
            // sum((1-row(trans_t))*Items) = 0 <=> trans_t
            linear(*this, row_, items, IRT_EQ, 0, transactions[t]);
        }
    }

    /** closed constraints **/
    {
        IntArgs col_(nr_t);
        for (int i=0; i!=nr_i; i++) {
            // make 1-col
            for (int t=0; t!=nr_t; t++)
                col_[t] = (1-tdb[t][i]);

            // closed: the item its complement has no supported trans
            // sum((1-col(item_i))*Trans) = 0 <=> item_i
            linear(*this, col_, transactions, IRT_EQ, 0, items[i]);
        }
    }

    /** convex hull constraint **/
    {
        // for every item separately (reified, allows iterative pruning)
        for (int i=0; i!=nr_i; i++) {
            // count cols
            int pos = 0; int neg = 0;
            for (int t=0; t!=nr_t; t++) {
                if (classes[t])
                    pos += tdb[t][i];
                else if (classes[t] == 0)
                    neg += tdb[t][i];
                else {
                    throw Exception("Class label error",
                        "illegal class label found, only '1' (pos) or '0' (neg) allowed");
                }
            }
            // make col
            BoolVarArgs col_pos(pos);
            BoolVarArgs col_neg(neg);
            for (int t=0; t!=nr_t; t++) {
                if (tdb[t][i] == 1) {
                    if (classes[t] == 1)
                        col_pos[--pos] = transactions[t];
                    else
                        col_neg[--neg] = transactions[t];
                }
            }

            // items[i] -> convexhull(col)
            ConvexHull::imply_convex_hull(*this, items[i], negTot, col_neg, posTot, col_pos, hullUL, hullLR);
        }
    }

    /** search **/
    branch(*this, items, (IntVarBranch)opt.branching(), (IntValBranch)opt.branchval());
}

/// Construct HullPoint from current solution
HullPoint* Cimcp_convexhull::constructHullPoint() const {
    vector<bool> curItems(nr_i);
    for (int i = 0; i!=nr_i; i++) {
        curItems[i] = items[i].val();
    }
    vector<bool> curTrans(nr_t);
    int curPos = 0; int curNeg = 0;
    for (int t = 0; t!=nr_t; t++) {
        curTrans[t] = transactions[t].val();
        if (classes[t] == 1) {
            curPos += curTrans[t];
        } else {
            curNeg += curTrans[t];
        }
    }
    return new HullPoint(curItems, curTrans, curNeg, curPos);
}

/// Print solution
void Cimcp_convexhull::print(std::ostream& os) const {
    // print() is called before constrain()
    // we add the itemset here because constrain() is never called for the last item
    // this function is const so our first hullpoints are a filler

    HullPoint* cur = constructHullPoint();

    // add HullPoint to convex hull
    if (cur->isLeftorOn(hullUL,hullUL->prev)) {
        // in UL
        ConvexHull::addToHullUL(cur, hullUL);
    } else {
        // in LR
        ConvexHull::addToHullLR(cur, hullLR);
    }

    // actual printing
    if (print_itemsets == PRINT_NONE) {
        return;
    } else if (print_itemsets == PRINT_FIMI) {

        // print hull
        fprintf(solfile, "Convex Hull:\n");
        cur = hullUL;
        while (cur) {
            cur->print(solfile);
            cur = cur->next;
        }
        cur = hullLR;
        while (cur != NULL) {
            cur->print(solfile);
            cur = cur->next;
        }
    } else if (print_itemsets == PRINT_FIMI) {

        // print hull (with transactions)
        fprintf(solfile, "Convex Hull:\n");
        cur = hullUL;
        while (cur) {
            cur->print_full(solfile);
            cur = cur->next;
        }
        cur = hullLR;
        while (cur != NULL) {
            cur->print_full(solfile);
            cur = cur->next;
        }
    } else if (print_itemsets == PRINT_CPVARS) {
        // Output CP variables (mostly for GIST)
        os << "Only printing current solution:\n";
        os << "\tI[] = " << items << std::endl;
        os << "\tT[] = " << transactions << std::endl;
    }
}

int main(int argc, char* argv[]) {
    Options_fimcp opt(strpbrk(argv[0],"/\\")+1);
    opt.description("This model finds all patterns on the convex hull (for both classes)");
    opt.usage("-datafile example.txt\n\
                -datafile example.txt\t itemset file where last item is the class: 0 or 1");
    opt.parse(argc, argv);

    fprintf(stdout, "Running BAB search for ");
    Script::run<Cimcp_convexhull,BAB,Options_fimcp>(opt);

    if (opt.threads() != 1) {
      fprintf(stderr, "Warning: parallell search could result in non-convex hull, reverting to 1 thread.");
      opt.threads(1);
    }
    return 0;
}
// avoid 'undefined reference' linking error
void Fimcp_basic::run(const Options_fimcp&) {};
