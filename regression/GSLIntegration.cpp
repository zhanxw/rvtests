#include "GSLIntegration.h"

#include <stdio.h>
#include <gsl/gsl_errno.h>

Integration::Integration()
    : epsabs(0.),
      epsrel(1e-7),
      limit(1000),
      workspace(NULL),
      result(0.),
      abserr(0.) {
  workspace = gsl_integration_workspace_alloc(limit);
}

Integration::~Integration() {
  if (workspace) {
    gsl_integration_workspace_free(workspace);
    workspace = NULL;
  }
}

int Integration::integrate(gsl_function F) {
  // https://www.gnu.org/software/gsl/manual/html_node/QAGI-adaptive-integration-on-infinite-intervals.html#QAGI-adaptive-integration-on-infinite-intervals
  int ret = gsl_integration_qagi(&F, epsabs, epsrel, limit, workspace, &result,
                                 &abserr);
  if (ret) {
    fprintf(stderr, "Integration failed with error [ %s ]", gsl_strerror(ret));
  }
  return ret;
}
