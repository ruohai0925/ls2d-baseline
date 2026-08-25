// main_flow.h — entry point of the Navier–Stokes driver.
#pragma once
#include "params.h"

// Run the flow problem described by p to completion; returns 0 on success, 2 on blow-up.
int run_flow(const Params& p);
