// SPDX-License-Identifier: BSD-3-Clause

#ifndef NVTX_H
#define NVTX_H

#include <string.h>
#include <nvtx3/nvToolsExt.h>

void mynvtxstart_(const char *name);

void mynvtxstop_();

#endif 