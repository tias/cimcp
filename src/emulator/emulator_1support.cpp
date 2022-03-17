/*
 *  Main authors:
 *      Tias Guns <tias.guns@cs.kuleuven.be>
 *
 *  Copyright:
 *      Tias Guns, 2008
 *
 *  Revision information:
 *      $Id: emulator_1support.cpp 275 2009-07-09 12:19:34Z tias $
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

/// information gain function for itemsets
inline float nlogn ( float n ) {
    // Calculate -n*log(n) with 0 if n=0
    if ( !n )
        return 0;
    return -log(n) * n;
}
inline int convex_function(int X, int pos_x, int Y, int pos_y, int precision) {
    float tot = X + Y;
    float base = nlogn(X/tot) + nlogn(Y/tot);
    
    float pos = pos_x + pos_y; 
    float neg = tot - pos; 
    float calc = pos * ( nlogn(    pos_x/pos) + nlogn(    pos_y/pos) ) +
                 neg * ( nlogn((X-pos_x)/neg) + nlogn((Y-pos_y)/neg) );
    if (isnan(calc)) calc = -1;
    return (int) (precision*(base - calc/tot));
}
inline bool convex_functionWrapper(int X,int l_x,int u_x,int Y,int l_y,int u_y,int precision, int tresholdLE) {
    //DDPmine: take joined freq treshold
    return (convex_function(X, 0, Y, u_x+u_y, precision) < tresholdLE
                && convex_function(X, u_x+u_y, Y, 0, precision) < tresholdLE);
}
#include "constraint_Fconvex_emulator.cpp" // emulator for convex function of boolean variables

#define PRECISION 100000


/**
 * Correlated itemset mining using information gain,
 * propagation as by the DDPmine system (Direct discriminative pattern mining for effective classification, H. Cheng, X. Yan, J. Han and C.-W. Hsu)
 * inherits from Fimcp_basic.
 */
class Emulator_1support : public Fimcp_basic {
protected:
  /// Correlation score. integer: precision defined by -alpha
  IntVar score;

public:
    Emulator_1support(const Options_fimcp& opt) :
        Fimcp_basic(opt) {
        run(opt);
    }

    Emulator_1support(bool share, Emulator_1support& s) :
        Fimcp_basic(share, s) {
        score.update(*this, share, s.score);
    }

    /// Perform copying during cloning
    virtual Space* copy(bool share) {
      return new Emulator_1support(share,*this);
    }

    /// Add constraint for next better solution
    virtual void constrain(const Space& _best);

    virtual void run(const Options_fimcp& opt);

    /// Print solution
    virtual void print(std::ostream&) const;
};

void Emulator_1support::constrain(const Space& _best) {
    const Emulator_1support* best = 
        dynamic_cast<const Emulator_1support*>(&_best);
    if (best == NULL)
        throw DynamicCastFailed("OptimizeSpace::constrain");
    // calculate convex_function value
    int posTot = 0; int pos = 0;
    int negTot = 0; int neg = 0;
    for (int t=0; t!=nr_t; t++) {
        posTot += classes[t];
        if (classes[t])
            pos += best->transactions[t].val();
        else
            neg += best->transactions[t].val();
    }
    negTot = nr_t - posTot;
    int treshold = convex_function(posTot, pos, negTot, neg, PRECISION);

    //std::cerr << " upping to "<<treshold<<": "<<score<<std::endl;
    rel(*this, score, IRT_GR, treshold);
}

void Emulator_1support::run(const Options_fimcp& opt) {
    const vector< vector<bool> > tdb = common_construction(opt);
    if (classes.size() == 0)
        throw Exception("Class label error", "no class labels found");

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

    /** extra minimum-frequency constraint **/
    int freq = opt.getFreq(nr_t);
    linear(*this, transactions, IRT_GQ, freq);

    /** search **/
    branch(*this, items, (IntVarBranch)opt.branching(), (IntValBranch)opt.branchval());
}

/// Print solution
void Emulator_1support::print(std::ostream& os) const {
    if (print_itemsets == PRINT_NONE) {
        return;
    } else if (print_itemsets == PRINT_FIMI) {
        // FIMI style output
        for (int i=0; i!=nr_i; i++) {
            if (items[i].val() == 1)
                fprintf(solfile, "%i ", i);
        }
        // calculate convex_function value
        int posTot = 0; int pos = 0;
        int negTot = 0; int neg = 0;
        for (int t=0; t!=nr_t; t++) {
            posTot += classes[t];
            if (classes[t])
                pos += transactions[t].val();
            else
                neg += transactions[t].val();
        }
        negTot = nr_t - posTot;
        float val = convex_function(posTot, pos, negTot, neg, PRECISION)/(float)PRECISION;
        fprintf(solfile, "(%i) [%.5f]\n", pos+neg, val);
    } else if (print_itemsets == PRINT_CPVARS) {
        // Output CP variables (mostly for GIST)
        os << "\tscore = " << score << std::endl;
        os << "\tI[] = " << items << std::endl;
        os << "\tT[] = " << transactions << std::endl;
    }
}

int main(int argc, char* argv[]) {
    Options_fimcp opt(strchr(argv[0],'/')+1);
    opt.freq(1);
    opt.delta(0);
    opt.alpha(1);
    opt.description("This model finds discriminative patterns, using information gain as measure and the 1-support bound");
    opt.usage("-datafile example.txt -freq 1\n\
                -datafile example.txt\t itemset file where last item is the class: 0 or 1\n\
                -alpha 1\t use branch-and-bound search for top-1 itemset\n\
                -alpha 0\t to find all patternsets given tresholds:\n\
                \t -delta 0.50\t minimal measure value");
    opt.parse(argc, argv);

    if (opt.alpha() == 0) {
      fprintf(stdout, "Running DF search for ");
      Script::run<Emulator_1support,DFS,Options_fimcp>(opt);
    } else {
      fprintf(stdout, "Running BAB search for ");
      Script::run<Emulator_1support,BAB,Options_fimcp>(opt);
    }
    return 0;
}
// avoid 'undefined reference' linking error
void Fimcp_basic::run(const Options_fimcp&) {};
