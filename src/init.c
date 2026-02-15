
// clang-format sorts includes unless SortIncludes: Never. However, the ordering
// does matter here. So, we need to disable clang-format for safety.

// clang-format off
#include <stdint.h>
#include <Rinternals.h>
#include <R_ext/Parse.h>
// clang-format on

#include "rust/api.h"

static uintptr_t TAGGED_POINTER_MASK = (uintptr_t)1;

SEXP handle_result(SEXP res_) {
    uintptr_t res = (uintptr_t)res_;

    // An error is indicated by tag.
    if ((res & TAGGED_POINTER_MASK) == 1) {
        // Remove tag
        SEXP res_aligned = (SEXP)(res & ~TAGGED_POINTER_MASK);

        // Currently, there are two types of error cases:
        //
        //   1. Error from Rust code
        //   2. Error from R's C API, which is caught by R_UnwindProtect()
        //
        if (TYPEOF(res_aligned) == CHARSXP) {
            // In case 1, the result is an error message that can be passed to
            // Rf_errorcall() directly.
            Rf_errorcall(R_NilValue, "%s", CHAR(res_aligned));
        } else {
            // In case 2, the result is the token to restart the
            // cleanup process on R's side.
            R_ContinueUnwind(res_aligned);
        }
    }

    return (SEXP)res;
}

SEXP savvy_clarabel_solve__impl(SEXP c_arg__m, SEXP c_arg__n, SEXP c_arg__Ai, SEXP c_arg__Ap, SEXP c_arg__Ax, SEXP c_arg__b, SEXP c_arg__q, SEXP c_arg__Pi, SEXP c_arg__Pp, SEXP c_arg__Px, SEXP c_arg__cone_spec, SEXP c_arg__r_settings) {
    SEXP res = savvy_clarabel_solve__ffi(c_arg__m, c_arg__n, c_arg__Ai, c_arg__Ap, c_arg__Ax, c_arg__b, c_arg__q, c_arg__Pi, c_arg__Pp, c_arg__Px, c_arg__cone_spec, c_arg__r_settings);
    return handle_result(res);
}

SEXP savvy_ClarabelSolver_is_update_allowed__impl(SEXP self__) {
    SEXP res = savvy_ClarabelSolver_is_update_allowed__ffi(self__);
    return handle_result(res);
}

SEXP savvy_ClarabelSolver_new__impl(SEXP c_arg__m, SEXP c_arg__n, SEXP c_arg__Ai, SEXP c_arg__Ap, SEXP c_arg__Ax, SEXP c_arg__b, SEXP c_arg__q, SEXP c_arg__Pi, SEXP c_arg__Pp, SEXP c_arg__Px, SEXP c_arg__cone_spec, SEXP c_arg__r_settings) {
    SEXP res = savvy_ClarabelSolver_new__ffi(c_arg__m, c_arg__n, c_arg__Ai, c_arg__Ap, c_arg__Ax, c_arg__b, c_arg__q, c_arg__Pi, c_arg__Pp, c_arg__Px, c_arg__cone_spec, c_arg__r_settings);
    return handle_result(res);
}

SEXP savvy_ClarabelSolver_solve__impl(SEXP self__) {
    SEXP res = savvy_ClarabelSolver_solve__ffi(self__);
    return handle_result(res);
}

SEXP savvy_ClarabelSolver_update_data__impl(SEXP self__, SEXP c_arg__Px, SEXP c_arg__Ax, SEXP c_arg__q, SEXP c_arg__b) {
    SEXP res = savvy_ClarabelSolver_update_data__ffi(self__, c_arg__Px, c_arg__Ax, c_arg__q, c_arg__b);
    return handle_result(res);
}


static const R_CallMethodDef CallEntries[] = {
    {"savvy_clarabel_solve__impl", (DL_FUNC) &savvy_clarabel_solve__impl, 12},
    {"savvy_ClarabelSolver_is_update_allowed__impl", (DL_FUNC) &savvy_ClarabelSolver_is_update_allowed__impl, 1},
    {"savvy_ClarabelSolver_new__impl", (DL_FUNC) &savvy_ClarabelSolver_new__impl, 12},
    {"savvy_ClarabelSolver_solve__impl", (DL_FUNC) &savvy_ClarabelSolver_solve__impl, 1},
    {"savvy_ClarabelSolver_update_data__impl", (DL_FUNC) &savvy_ClarabelSolver_update_data__impl, 5},
    {NULL, NULL, 0}
};

void R_init_clarabel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

    // Functions for initialization, if any.

}
