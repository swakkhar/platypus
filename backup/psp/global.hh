//#define RefHdl(Assignable, Class, Sys, Hdl)
//		Class const & Assignable = Class::ref(Sys, Hdl)
//
//
//#define RefPrm(Assignable, Class, Sys, Prm)
//		Class const & Assignable = Class::ref(Sys, Prm.TermHdl)


#define PtrHdls(Assignables, Class, Sys, Handles, Size) \
		Class const * Assignables[Size]; \
		for (Idx tIdx = 0 ; tIdx < Size ; ++tIdx) \
            Assignables[tIdx] = Class::ptr(Sys, Handles[tIdx]);

#define PtrHdlsSE(Assignables, Class, Sys, Handles, Start, End) \
		Class const * Assignables[End-Start]; \
		for (Idx tIdx = 0 ; tIdx < End - Start ; ++tIdx) \
            Assignables[tIdx] = Class::ptr(Sys, Handles[Start+tIdx]);

#define PtrPrmsSE(Assignables, Class, Sys, Params, Start, End) \
 Class const * Assignables[End-Start]; \
	for (Idx tIdx = 0 ; tIdx < End - Start ; ++tIdx) \
           Assignables[tIdx] = Class::ptr(Sys, Params[Start+tIdx].TermHdl);

#define PtrPrms(Assignables, Class, Sys, Params, Size) \
		Class const * Assignables[Size]; \
		for (Idx tIdx = 0 ; tIdx < Size ; ++tIdx) \
            Assignables[tIdx] = Class::ptr(Sys, Params[tIdx].TermHdl);

//#define PtrHdl(Assignable, Class, Sys, Handle)
//		Class const * Assignable = Class::ptr(Sys, Handle)
//
//
//
//#define PtrPrm(Assignable, Class, Sys, Param)
//		Class const * Assignable = Class::ptr(Sys, Param.TermHdl)
//
