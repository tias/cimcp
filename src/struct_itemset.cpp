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

  /// Itemset in PN-space:
struct PNitemset {
  public:
    // Next pointer
    PNitemset* next;
    // Items of this point
    const vector<bool> items;
    // Transactions of this point
    const vector<bool> trans;
    // Nr transactions of pos, and original total
    const int pos;
    const int posTot;
    // Nr transactions of neg, and original total
    const int neg;
    const int negTot;
    // Score (as used in the model)
    const int score;
    // actual value of evaluation function (no PRECISION)
    const float value;

    // Constructor
    PNitemset(vector<bool> items0, vector<bool> trans0,
                int pos0, int posTot0, int neg0, int negTot0, int score0, float value0):
        items(items0), trans(trans0),
        pos(pos0), posTot(posTot0), neg(neg0), negTot(negTot0), score(score0), value(value0) {
        next = NULL;
    }
    // Phony for skipping const of print function where we have to add
    PNitemset(): items(), trans(), pos(0), posTot(0), neg(0), negTot(0), score(0), value(0) {
        next = NULL;
    }

    inline void print(FILE* solfile) {
        for (unsigned int i=0; i!=items.size(); i++) {
            if (items[i])
                fprintf(solfile, "%i ", i);
        }
        //fprintf(solfile, "(%i:+%i-%i) [%.5f], score is %i\n", pos+neg,pos,neg, value, score);
        fprintf(solfile, "(%i:+%i-%i) [%.5f]\n", pos+neg,pos,neg, value);
    }

    inline void print_full(FILE* solfile) {
        for (unsigned int i=0; i!=items.size(); i++) {
            if (items[i])
                fprintf(solfile, "%i ", i);
        }
        fprintf(solfile, "(%i) < ", pos+neg);
        for (unsigned int t=0; t!=trans.size(); t++) {
            if (trans[t])
                fprintf(solfile, "%i ", t);
        }
        fprintf(solfile, ">\n");
    }
};
