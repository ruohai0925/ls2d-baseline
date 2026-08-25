// main_flow.h — entry point of the Navier–Stokes driver (called from main() when prob.flow = 1).
#pragma once
#include "params.h"

// Run a flow problem described by p to completion; returns 0 on success, 2 on blow-up.
int run_flow(const Params& p);
