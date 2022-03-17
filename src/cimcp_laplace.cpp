/*
 *  Main authors:
 *      Tias Guns <tias.guns@cs.kuleuven.be>
 *
 *  Copyright:
 *      Tias Guns, 2008
 *
 *  Revision information:
 *      $Id: pattset_concept_accuracy.cpp 217 2009-04-07 13:37:27Z tias $
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

/// laplace measure (2 classes)
inline int monotone_function(int, int pos, int, int neg, int precision) {
    return ( (pos+1) /(float)(pos+neg+2) )*precision;
}
#include "cimcp_Fmonotone.cpp" // model using monotone function of boolean variables

int main(int argc, char* argv[]) {
    Options_fimcp opt(strpbrk(argv[0],"/\\")+1);
    opt.delta(0);
    opt.alpha(1);
    opt.description("This model finds discriminative patterns, using laplace as measure");
    opt.usage("-datafile example.txt\n\
                -datafile example.txt\t itemset file where last item is the class: 0 or 1\n\
                -alpha k\t (default 1) use branch-and-bound search for top-k itemsets\n\
                -alpha 0\t to find all patternsets given tresholds:\n\
                \t -delta 0.50\t minimal measure value");
    opt.parse(argc, argv);

    if (opt.alpha() == 0) {
      fprintf(stdout, "Running DF search for ");
      Script::run<Cimcp_Fmonotone,DFS,Options_fimcp>(opt);
    } else {
      fprintf(stdout, "Running BAB search for ");
      Script::run<Cimcp_Fmonotone,BAB,Options_fimcp>(opt);
    }
    return 0;
}
// avoid 'undefined reference' linking error
void Fimcp_basic::run(const Options_fimcp&) {};
