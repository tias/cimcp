/*
 *  Main authors:
 *      Tias Guns <tias.guns@cs.kuleuven.be>
 *
 *  Copyright:
 *      Tias Guns, 2008
 *
 *  Revision information:
 *      $Id: emulator_2support.cpp 275 2009-07-09 12:19:34Z tias $
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
#include "struct_itemset.cpp"
#include <math.h>

/// information gain function for itemsets
inline float nlogn ( float n ) {
    // Calculate -n*log(n) with 0 if n=0
    if ( !n )
        return 0;
    return -log(n) * n;
}
inline int convex_function(int posTot, int pos, int negTot, int neg, int precision) {
    float tot = posTot + negTot;
    float base = nlogn(posTot/tot) + nlogn(negTot/tot);
    
    float covered = pos + neg; 
    float notCovered = tot - covered; 
    float calc = covered    * ( nlogn(pos/covered) + nlogn(neg/covered) ) +
                 notCovered * ( nlogn((posTot-pos)/notCovered) + nlogn((negTot-neg)/notCovered) );
    if (isnan(calc)) calc = -1;
    return (int) (precision*(base - calc/tot));
}
inline bool convex_functionWrapper(int X,int l_x,int u_x,int Y,int l_y,int u_y,int precision, int tresholdLE) {
    //Morishita: simple upper bound
    return (convex_function(X, u_x, Y, 0, precision) < tresholdLE
                && convex_function(X, 0, Y, u_y, precision) < tresholdLE);
}
#include "constraint_Fconvex_emulator.cpp" // emulator for convex function of boolean variables

#define PRECISION 100000

/**
 * Correlated itemset mining using information gain,
 * propagation as described in (Answering the most correlated n association rules efficiently, J. Sese and S. Morishita)
 * they do not mention closedness, but this is needed to avoid equivalent solutions.
 * inherits from Fimcp_basic.
 */
class Emulator_2support : public Fimcp_basic {
protected:
  /// Correlation score. integer: precision defined by -alpha
  IntVar score;

  /// top-k itemsets (is linked list struct, with filler first)
  PNitemset* itemsets;

  /// k, nr of top itemsets to mine
  int k;

public:
    Emulator_2support(const Options_fimcp& opt) :
        Fimcp_basic(opt) {
        // have to set filler first, otherwise we cannot
        // add items in the const print() function.
        PNitemset* filler = new PNitemset();
        itemsets = filler;

        run(opt);
    }

    Emulator_2support(bool share, Emulator_2support& s) :
        Fimcp_basic(share, s) {
        score.update(*this, share, s.score);
        itemsets = s.itemsets;
        k = s.k;
    }

    /// Perform copying during cloning
    virtual Space* copy(bool share) {
      return new Emulator_2support(share,*this);
    }

    /// Add constraint for next better solution
    virtual void constrain(const Space&) {
        // constrain() is called after print()
        // we add the itemset to the topk in print()
        // because constrain() is never called for the last item

        // set the treshold wrt. itemset k
        PNitemset* last = itemsets;
        while (last->next != NULL)
            last = last->next;
        rel(*this, score, IRT_GR, last->score);
    }

    virtual void run(const Options_fimcp& opt);

    /// Print solution
    virtual void print(std::ostream&) const;

    /// Construct PNitemset from current solution (has to be const because print is)
    virtual PNitemset* constructPNitemset() const;

    /// add itemset to the topk solutions (has to be const because print is)
    virtual void addItemset(PNitemset* cur) const;
};

void Emulator_2support::run(const Options_fimcp& opt) {
    const vector< vector<bool> > tdb = common_construction(opt);
    if (classes.size() == 0)
        throw Exception("Class label error", "no class labels found");
    // number of top-k itemsets
    k = opt.alpha();

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

    /** emulated convex function constraint **/
    {
        // convert treshold to integer, default: [min*PRECISION, _MAX]
        score = IntVar(*this, (int)(opt.delta()*PRECISION), Gecode::Int::Limits::max);

        // for every item separately (reified, allows iterative pruning)
        for (int i=0; i!=nr_i; i++) {
            // count cols
            int posTot = 0; int pos = 0;
            int negTot = 0; int neg = 0;
            for (int t=0; t!=nr_t; t++) {
                posTot += classes[t];
                if (classes[t])
                    pos += tdb[t][i];
                else if (classes[t] == 0)
                    neg += tdb[t][i];
                else {
                    throw Exception("Class label error",
                        "illegal class label found, only '1' (pos) or '0' (neg) allowed");
                }
            }
            negTot = nr_t - posTot;
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

            // items[i] -> convex(col) >= delta
            FconvexBool_emulator::imply_convex(*this, items[i], posTot, col_pos, negTot, col_neg, IRT_GQ, score, PRECISION);
        }
    }

    // avoid empty itemset
    linear(*this, items, IRT_GQ, 1);

    /** search **/
    branch(*this, items, (IntVarBranch)opt.branching(), (IntValBranch)opt.branchval());
}

/// Construct PNitemset from current solution
inline PNitemset* Emulator_2support::constructPNitemset() const {
    vector<bool> curItems(nr_i);
    for (int i = 0; i!=nr_i; i++) {
        curItems[i] = items[i].val();
    }
    vector<bool> curTrans(nr_t);
    int curPos = 0; int curNeg = 0;
    int posTot = 0; int negTot = 0;
    for (int t = 0; t!=nr_t; t++) {
        curTrans[t] = transactions[t].val();
        if (classes[t] == 1) {
            curPos += curTrans[t];
            posTot += 1;
        } else {
            curNeg += curTrans[t];
            negTot += 1;
        }
    }
    int treshold = convex_function(posTot, curPos, negTot, curNeg, PRECISION);
    return new PNitemset(curItems, curTrans,
            curPos, posTot, curNeg, negTot, treshold, treshold/(float)PRECISION);
}

/// Print solution
void Emulator_2support::print(std::ostream& os) const {
    // print() is called before constrain()
    // we add the itemset here because constrain() is never called for the last item
    // this function is const so our first itemset is a filler

    if (print_itemsets == PRINT_CPVARS) {
        // Output CP variables (usually for GIST)
        os << "Only printing current solution:\n";
        os << "\tscore = " << score << std::endl;
        os << "\tI[] = " << items << std::endl;
        os << "\tT[] = " << transactions << std::endl;
        // if not all vars assigned (partial solution printed in GIST),
        // return here
        bool all_assigned = true;
        for (int i = 0; i!=nr_i && all_assigned; i++)
            all_assigned &= items[i].assigned();
        for (int t = 0; t!=nr_t && all_assigned; t++)
            all_assigned &= transactions[t].assigned();
        if (!all_assigned)
            return;
    }

    if (k == 0) { // top-0: all (but none in memory)
        if (print_itemsets == PRINT_NONE)
            return;
    
        PNitemset* cur = constructPNitemset();

        if (print_itemsets == PRINT_FIMI) {
            cur->print(solfile);
        } else if (print_itemsets == PRINT_FULL) {
            cur->print_full(solfile);
        } else if (print_itemsets == PRINT_CPVARS) {
            cur->print(stdout);
        }

        delete cur;

    } else { // top-k (keep k in memory)
        PNitemset* cur = constructPNitemset();
        addItemset(cur);

        if (print_itemsets == PRINT_NONE) {
            return;
        } else if (print_itemsets == PRINT_FIMI) {
            // print entire patternset
            PNitemset* iter = itemsets;
            while (iter->next != NULL) {
                iter = iter->next; // skip first filler
                iter->print(solfile);
            }
        } else if (print_itemsets == PRINT_FULL) {
            // print entire patternset
            PNitemset* iter = itemsets;
            while (iter->next != NULL) {
                iter = iter->next; // skip first filler
                iter->print_full(solfile);
            }
        } else if (print_itemsets == PRINT_CPVARS) {
            // also print normal output, but on stdout
            PNitemset* iter = itemsets;
            while (iter->next != NULL) {
                iter = iter->next; // skip first filler
                iter->print(stdout);
            }
        }
    }
}

/// add itemset to the topk solutions
inline void Emulator_2support::addItemset(PNitemset* cur) const {
    if (k == 0)
        return;

    PNitemset* prev = itemsets;
    int count = 0;
    while (prev->next != NULL &&
            prev->next->score >= cur->score) {
        prev = prev->next;
        count++;
    }

    cur->next = prev->next;
    prev->next = cur;
    count++;

    while (prev->next != NULL &&
            count <= k) {
        prev = prev->next;
        count++;
    }

    // clean up
    delete prev->next;
    prev->next = NULL;

    if (print_itemsets != PRINT_NONE)
        fprintf(stdout, "\tAdded itemset, min value=%.5f:\n", prev->value);
    return;
}

int main(int argc, char* argv[]) {
    Options_fimcp opt(strchr(argv[0],'/')+1);
    opt.delta(0);
    opt.alpha(1);
    opt.description("This model finds discriminative patterns, using information gain as measure and the 2-support bound");
    opt.usage("-datafile example.txt\n\
                -datafile example.txt\t itemset file where last item is the class: 0 or 1\n\
                -alpha k\t (default 1) use branch-and-bound search for top-k itemsets\n\
                -alpha 0\t to find all patternsets given tresholds:\n\
                \t -delta 0.50\t minimal measure value");
    opt.parse(argc, argv);

    if (opt.alpha() == 0) {
      fprintf(stdout, "Running DF search for ");
      Script::run<Emulator_2support,DFS,Options_fimcp>(opt);
    } else {
      fprintf(stdout, "Running BAB search for ");
      Script::run<Emulator_2support,BAB,Options_fimcp>(opt);
    }
    return 0;
}
// avoid 'undefined reference' linking error
void Fimcp_basic::run(const Options_fimcp&) {};
