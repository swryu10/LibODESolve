#include<stdlib.h>
#include"ODESolveBase.h"

namespace ODESolve {

void LibBase::alloc_func_base() {
    tab_func_y_ =
        new std::vector<double> [n_func_y_ + 1];
    for (int i = 0; i <= n_func_y_; i++) {
        tab_func_y_[i].clear();
    }

    have_func_base_ = true;

    return;
}

void LibBase::free_func_base() {
    if (!have_func_base_) {
        return;
    }

    delete [] tab_func_y_;

    return;
}

void LibBase::export_file(FILE *ptr_fout) {
    if (!have_func_base_) {
        return;
    }

    if (ptr_fout == NULL) {
        return;
    }

    int mx = tab_func_y_[0].size();

    fprintf(ptr_fout, "# output from ODESolve\n");
    fprintf(ptr_fout, "# x");
    for (int i = 1; i <= n_func_y_; i++) {
        fprintf(ptr_fout, "  y%d", i);
    }
    fprintf(ptr_fout, "\n");

    for (int ix = 0; ix < mx; ix++) {
        fprintf(ptr_fout, "  %e", tab_func_y_[0][ix]);

        for (int i = 1; i <= n_func_y_; i++) {
            fprintf(ptr_fout, "  %e", tab_func_y_[i][ix]);
        }

        fprintf(ptr_fout, "\n");
    }
}

} // end namespace ODESolve
