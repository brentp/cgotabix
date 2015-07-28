package cgotabix

/*
#cgo LDFLAGS: -L. -lhts -lpthread -lz -lm
#include "stdlib.h"
#include <zlib.h>
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kfunc.h"


char *tbx_next(htsFile *fp, tbx_t *tbx, hts_itr_t *iter, int *slen) {
	kstring_t s = {0, 0, 0};
	*slen = tbx_itr_next(fp, tbx, iter, &s);
	return s.s;
}

hts_itr_t *tabix_itr_querys(tbx_t *tbx, char *s){
     return hts_itr_querys((tbx)->idx, (s), (hts_name2id_f)(tbx_name2id), (tbx), hts_itr_query, tbx_readrec);
 }


*/
import "C"
import "unsafe"

type Tabix struct {
	tbx *C.tbx_t
	htf *C.htsFile
}

// New takes a path to a bgziped (and tabixed file) and returns the tabix struct.
func New(path string) Tabix {
	var t Tabix
	cs := C.CString(path)
	defer C.free(unsafe.Pointer(cs))
	t.tbx = C.tbx_index_load(cs)
	mode := C.char('r')
	t.htf = C.hts_open(cs, &mode)
	return t
}

// Get takes a region like 1:45678-56789 and returns a channel on which
// it sends a string for each line that falls in that interval.
func (t Tabix) Get(region string) chan string {
	cs := C.CString(region)
	defer C.free(unsafe.Pointer(cs))

	itr := C.tabix_itr_querys(t.tbx, cs)
	out := make(chan string, 20)

	go func() {
		l := C.int(10)
		for l > 0 {
			cstr := C.tbx_next(t.htf, t.tbx, itr, &l)
			res := C.GoString(cstr)
			if l < 0 {
				break
			}
			out <- res
		}
		close(out)
	}()
	return out
}
