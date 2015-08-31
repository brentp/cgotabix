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
#include "htslib/vcf.h"
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

inline int atbx_itr_next(htsFile *fp, tbx_t *tbx, hts_itr_t *iter, kstring_t *data) {
	return tbx_itr_next(fp, tbx, iter, (void *)data);
}

char *allele_i(bcf1_t *b, int idx) {
	return b->d.allele[idx];
}

int itype(bcf_info_t *i){
	return i->type;
}

int ivi(bcf_info_t *i) {
	return i->v1.i;
}
float ivf(bcf_info_t *i) {
	return i->v1.f;
}


int ibcf_update_info_int32(const bcf_hdr_t *hdr, bcf1_t * line, const char *key, const int32_t *values, int n) {
	bcf_update_info((hdr),(line),(key),(values),(n),BCF_HT_INT);
}

int ibcf_update_info_float(const bcf_hdr_t *hdr, bcf1_t * line, const char *key, const float *values, int n){
	bcf_update_info((hdr),(line),(key),(values),(n),BCF_HT_REAL);
}
int ibcf_update_info_flag(const bcf_hdr_t *hdr, bcf1_t * line, const char *key, const char *string, int n){
	bcf_update_info((hdr),(line),(key),(string),(n),BCF_HT_FLAG);
}
int ibcf_update_info_string(const bcf_hdr_t *hdr, bcf1_t * line, const char *key, const char *string){
	bcf_update_info((hdr),(line),(key),(string),1,BCF_HT_STR);
}


*/
import "C"
import (
	"fmt"
	"log"
	"runtime"
	"strings"
	"unsafe"

	"github.com/brentp/irelate"
)

type FileType string

const (
	BED   FileType = "bed"
	VCF   FileType = "vcf"
	OTHER FileType = "other"
)

const BCF_BT_NULL int = 0
const BCF_BT_INT8 int = 1
const BCF_BT_INT16 int = 2
const BCF_BT_INT32 int = 3
const BCF_BT_FLOAT int = 5
const BCF_BT_CHAR int = 7

type Tabix struct {
	path            string
	tbx             *C.tbx_t
	htf             *C.htsFile
	hdr             *C.bcf_hdr_t
	itr             *C.hts_itr_t
	typ             FileType
	original_header bool
}

func tabixCloser(t *Tabix) {
	if t.tbx != nil {
		C.tbx_destroy(t.tbx)
	}
	if t.htf != nil {
		C.hts_close(t.htf)
	}
	if t.hdr != nil {
		C.bcf_hdr_destroy(t.hdr)
	}
}

// New takes a path to a bgziped (and tabixed file) and returns the tabix struct.
func New(path string) *Tabix {
	t := &Tabix{path: path}
	cs := C.CString(t.path)
	defer C.free(unsafe.Pointer(cs))
	mode := C.char('r')
	t.htf = C.hts_open(cs, &mode)
	t.tbx = C.tbx_index_load(cs)
	runtime.SetFinalizer(t, tabixCloser)
	if strings.HasSuffix(path, ".bed.gz") {
		t.typ = BED
	} else if strings.HasSuffix(path, ".vcf.gz") {
		t.hdr = C.bcf_hdr_read(t.htf)
		t.typ = VCF
	} else {
		t.typ = OTHER
	}
	t.original_header = true
	return t
}

func (t *Tabix) AddInfoToHeader(id string, number string, vtype string, description string) {
	if t.original_header {
		hdr := C.bcf_hdr_dup(t.hdr)
		C.bcf_hdr_destroy(t.hdr)
		t.hdr = hdr

	}
	ckey := C.CString(fmt.Sprintf("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\"", id, number, vtype, description))
	e := C.bcf_hdr_append(t.hdr, ckey)
	C.free(unsafe.Pointer(ckey))
	if e != 0 {
		log.Fatalf("couldn't add to header\n")
	}
	if C.bcf_hdr_sync(t.hdr) != 0 {
		log.Fatalf("error syncing header\n")
	}
}

type INFO struct {
	hdr *C.bcf_hdr_t
	b   *C.bcf1_t
}

func (i *INFO) Get(key string) interface{} {
	ckey := C.CString(key)
	info := C.bcf_get_info(i.hdr, i.b, ckey)
	C.free(unsafe.Pointer(ckey))
	if info == nil {
		return nil
	}
	return i.get(info)
}

func (i *INFO) Set(key string, ovalue interface{}) {
	ckey := C.CString(key)

	switch ovalue.(type) {
	case int:
		value := C.int32_t(ovalue.(int))
		C.ibcf_update_info_int32(i.hdr, i.b, ckey, &value, 1)
	case float32:
		log.Println("OOOOOO")
		value := C.float(ovalue.(float32))
		C.ibcf_update_info_float(i.hdr, i.b, ckey, &value, 1)
	case float64:
		value := C.float(ovalue.(float64))
		C.ibcf_update_info_float(i.hdr, i.b, ckey, &value, 1)
	case string:
		value := C.CString(ovalue.(string))
		C.ibcf_update_info_string(i.hdr, i.b, ckey, value)
		C.free(unsafe.Pointer(value))
	case bool:
		value := ovalue.(bool)

		if value {
			ret := C.ibcf_update_info_flag(i.hdr, i.b, ckey, ckey, 1)
		} else {
			C.ibcf_update_info_flag(i.hdr, i.b, ckey, ckey, 0)
		}
	case []uint32:
		value := ovalue.([]uint32)
		l := len(value)
		val := make([]int32, l)
		for i, v := range value {
			val[i] = int32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		C.ibcf_update_info_int32(i.hdr, i.b, ckey, (*C.int32_t)(ptr), C.int(l))

	case []int:
		value := ovalue.([]int)
		l := len(value)
		val := make([]int32, l)
		for i, v := range value {
			val[i] = int32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		C.ibcf_update_info_int32(i.hdr, i.b, ckey, (*C.int32_t)(ptr), C.int(l))

	case []uint64:
		value := ovalue.([]uint64)
		l := len(value)
		val := make([]int32, l)
		for i, v := range value {
			val[i] = int32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		C.ibcf_update_info_int32(i.hdr, i.b, ckey, (*C.int32_t)(ptr), C.int(l))

	case []int32:
		value := ovalue.([]int32)
		l := len(value)
		val := make([]int32, l)
		for i, v := range value {
			val[i] = int32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		C.ibcf_update_info_int32(i.hdr, i.b, ckey, (*C.int32_t)(ptr), C.int(l))

	case []int64:
		value := ovalue.([]int64)
		l := len(value)
		val := make([]int32, l)
		for i, v := range value {
			val[i] = int32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		C.ibcf_update_info_int32(i.hdr, i.b, ckey, (*C.int32_t)(ptr), C.int(l))

	case []float64:
		value := ovalue.([]float32)
		l := len(value)
		val := make([]float32, l)
		for i, v := range value {
			val[i] = float32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		ret := C.ibcf_update_info_float(i.hdr, i.b, ckey, (*C.float)(ptr), C.int(l))
		if ret != 0 {
			log.Printf("not set %s %v\n", key, ovalue)
		}

	case []float32:
		val := ovalue.([]float32)
		l := len(val)
		ptr := unsafe.Pointer(&val[0])
		ret := C.ibcf_update_info_float(i.hdr, i.b, ckey, (*C.float)(ptr), C.int(l))
		if ret != 0 {
			log.Printf("not set %s %v\n", key, ovalue)
		}

	case []string:
		valuestr := ovalue.([]string)
		value := C.CString(strings.Join(valuestr, ","))
		C.ibcf_update_info_string(i.hdr, i.b, ckey, value)
		C.free(unsafe.Pointer(value))
	case []bool:
		// not possible
	default:
		log.Printf("multipe values for type: %T; not implemented (key: %s)\n", ovalue, key)

	}
	C.free(unsafe.Pointer(ckey))
}

func (self *INFO) get(info *C.bcf_info_t) interface{} {
	ctype := C.itype(info)
	if info.len == 1 {
		switch ctype {
		case C.BCF_BT_INT8, C.BCF_BT_INT16, C.BCF_BT_INT32:
			switch ctype {
			case C.BCF_BT_INT8:
				if C.ivi(info) == C.INT8_MIN {
					return nil
				}
				return int(C.ivi(info))
			case C.BCF_BT_INT16:
				if C.ivi(info) == C.INT16_MIN {
					return nil
				}
				return int(C.ivi(info))
			case C.BCF_BT_INT32:
				if C.ivi(info) == C.INT32_MIN {
					return nil
				}
				return int(C.ivi(info))
			}
		case C.BCF_BT_FLOAT:
			if C.bcf_float_is_missing(C.ivf(info)) != 0 {
				return nil
			}
			return float64(C.ivf(info))
		case C.BCF_BT_CHAR:
			return _string(info)
		}

	}

	if int(ctype) == BCF_BT_CHAR {
		return _string(info)
	}

	// divide this number by the number of bytes per object
	n_bytes := int(info.vptr_len)

	switch ctype {
	case C.BCF_BT_INT8:
		l := n_bytes
		out := make([]int, l)
		slice := (*[1 << 20]C.int8_t)(unsafe.Pointer(info.vptr))[:l:l]
		for i := 0; i < l; i++ {
			//if slice[i] == C.bcf_int8_missing {
			//		out[i] = nil
			if slice[i] == C.bcf_int8_vector_end {
				return out[:i]
			} else {
				out[i] = int(slice[i])
			}
		}
		return out

	case C.BCF_BT_INT16:
		l := n_bytes / 2
		out := make([]int, l)
		slice := (*[1 << 20]C.int16_t)(unsafe.Pointer(info.vptr))[:l:l]
		for i := 0; i < l; i++ {
			//if slice[i] == C.bcf_int16_missing {
			//		out[i] = nil
			if slice[i] == C.bcf_int16_vector_end {
				return out[:i]
			} else {
				out[i] = int(slice[i])
			}
		}
		return out

	case C.BCF_BT_INT32:
		l := n_bytes / 4
		out := make([]int, l)
		slice := (*[1 << 20]C.int32_t)(unsafe.Pointer(info.vptr))[:l:l]
		for i := 0; i < l; i++ {
			//if slice[i] == C.bcf_int32_missing {
			//		out[i] = nil
			if slice[i] == C.bcf_int32_vector_end {
				return out[:i]
			} else {
				out[i] = int(slice[i])
			}
		}
		return out

	case C.BCF_BT_FLOAT:
		l := n_bytes / 4
		out := make([]float32, l)
		slice := (*[1 << 20]C.float)(unsafe.Pointer(info.vptr))[:l:l]
		for i := 0; i < l; i++ {
			if C.bcf_float_is_vector_end(slice[i]) != 0 {
				return out[:i]
			} else {
				out[i] = float32(slice[i])
			}
		}
		return out

	}
	return nil
}

func _string(info *C.bcf_info_t) interface{} {
	s, l := info.vptr, info.vptr_len
	// https://github.com/golang/go/wiki/cgo
	slice := (*[1 << 20]uint8)(unsafe.Pointer(s))[:l:l]
	if slice[0] == 0x7 {
		return nil
	}
	return string(slice)
}

type CVariant struct {
	v       *C.bcf1_t
	hdr     *C.bcf_hdr_t
	source  uint32
	Info    INFO
	related []irelate.Relatable
	Pos     uint64
	Ref     string
	Alt     []string
}

func NewCVariant(v *C.bcf1_t, hdr *C.bcf_hdr_t, source uint32) *CVariant {
	C.bcf_unpack(v, 1|2|4) // dont unpack genotypes
	c := &CVariant{v: v, hdr: hdr, source: 1}
	c.Ref = C.GoString(C.allele_i(v, 0))
	c.Pos = uint64(v.pos + 1)
	c.Alt = make([]string, v.n_allele-1)
	for i := 1; i < int(v.n_allele); i++ {
		c.Alt[i-1] = C.GoString(C.allele_i(v, C.int(i)))
	}
	c.Info = INFO{hdr, v}
	runtime.SetFinalizer(c, variantFinalizer)
	return c
}

func variantFinalizer(v *CVariant) {
	C.bcf_destroy(v.v)
}

func (c *CVariant) Chrom() string {
	return C.GoString(C.bcf_hdr_id2name(c.hdr, C.int(c.v.rid)))
}

func (c *CVariant) Start() uint32 {
	return uint32(c.v.pos)
}

func (c *CVariant) End() uint32 {
	return uint32(c.v.pos + c.v.rlen)
}

func (c *CVariant) Source() uint32 {
	return c.source
}

func (c *CVariant) SetSource(s uint32) {
	c.source = s
}

func (c *CVariant) Related() []irelate.Relatable {
	return c.related
}

func (c *CVariant) AddRelated(i irelate.Relatable) {
	if c.related == nil {
		c.related = make([]irelate.Relatable, 0)
	}
	c.related = append(c.related, i)
}

// Get takes a region like 1:45678-56789 and returns a channel on which
// it sends a string for each line that falls in that interval.
func (t *Tabix) Get(region string) chan irelate.Relatable {
	cs := C.CString(region)
	defer C.free(unsafe.Pointer(cs))

	out := make(chan irelate.Relatable, 20)
	itr := C.tabix_itr_querys(t.tbx, cs)

	go func() {
		kstr := C.kstring_t{}
		l := C.int(10)
		for l > 0 {
			l := C.atbx_itr_next(t.htf, t.tbx, itr, &kstr)
			if l < 0 {
				break
			}
			//res := C.GoString(kstr.s)
			b := C.bcf_init()
			ret := C.vcf_parse(&kstr, t.hdr, b)
			if ret < 0 {
				log.Printf("error parsing %s\n", C.GoStringN(kstr.s, C.int(kstr.l)))
			}
			out <- NewCVariant(b, t.hdr, 1)
		}
		close(out)
		C.hts_itr_destroy(itr)
		C.free(unsafe.Pointer(kstr.s))
	}()
	return out
}
