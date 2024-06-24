// SPDX-License-Identifier: BSD-3-Clause

#ifndef NVTX_H
#define NVTX_H

#include <string.h>

#ifdef NSIGHT_PROFILING
    #include <nvtx3/nvToolsExt.h>
#endif

void mynvtxstart_(const char *name);

void mynvtxstop_();

#endif 