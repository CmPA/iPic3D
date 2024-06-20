// SPDX-License-Identifier: BSD-3-Clause

#include "nvtx.h"

#ifdef NSIGHT_PROFILING

	void mynvtxstart_(const char *name) 
	{
		int hash = 0;
		int color_id = 0x00C0C0C0; // very simple, all NVTX has same color (silver)

		nvtxEventAttributes_t eventAttrib = {0};
		eventAttrib.version = NVTX_VERSION;
		eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
		eventAttrib.colorType = NVTX_COLOR_ARGB;
		eventAttrib.color = color_id;
		eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
		eventAttrib.message.ascii = name;
		nvtxRangePushEx(&eventAttrib);
	}

	void mynvtxstop_() 
	{
		nvtxRangePop();
	}

#endif