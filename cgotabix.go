package cgotabix

/*
#cgo LDFLAGS: -lhts -lpthread -lz -lm

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

void format_info(bcf1_t *v, bcf_hdr_t *h, kstring_t *s){
  int i;
  if (v->n_info) {
	  int first = 1;
	  for (i = 0; i < v->n_info; ++i) {
		  bcf_info_t *z = &v->d.info[i];
		  if ( !z->vptr ) continue;
		  if ( !first ) kputc(';', s); first = 0;
		  kputs(h->id[BCF_DT_ID][z->key].key, s);
		  if (z->len <= 0) continue;
		  kputc('=', s);
		  if (z->len == 1)
		  {
			  switch (z->type)
			  {
				  case BCF_BT_INT8:  if ( z->v1.i==bcf_int8_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
				  case BCF_BT_INT16: if ( z->v1.i==bcf_int16_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
				  case BCF_BT_INT32: if ( z->v1.i==bcf_int32_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
				  case BCF_BT_FLOAT: if ( bcf_float_is_missing(z->v1.f) ) kputc('.', s); else ksprintf(s, "%g", z->v1.f); break;
				  case BCF_BT_CHAR:  kputc(z->v1.i, s); break;
				  default: fprintf(stderr,"todo: type %d\n", z->type); exit(1); break;
			  }
		  }
		  else bcf_fmt_array(s, z->len, z->type, z->vptr);
	  }
	  if ( first ) kputc('.', s);
      } else kputc('.', s);
}

char *tbx_next(htsFile *fp, tbx_t *tbx, hts_itr_t *iter, int *slen) {
	kstring_t s = {0, 0, 0};
	*slen = tbx_itr_next(fp, tbx, iter, &s);
	return s.s;
}

hts_itr_t *tabix_itr_querys(tbx_t *tbx, char *s){
     return hts_itr_querys((tbx)->idx, (s), (hts_name2id_f)(tbx_name2id), (tbx), hts_itr_query, tbx_readrec);
}

hts_itr_t *tabix_itr_queryi(tbx_t *tbx,  int tid, int beg, int end){
     return hts_itr_query((tbx)->idx, (tid), (beg), (end), tbx_readrec);
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

int ibcf_hdr_id2type(bcf_hdr_t *hdr, int htype, int tag_id){
	return bcf_hdr_id2type(hdr, htype, tag_id);
}

char *vid(bcf1_t *b) {
	return b->d.id;
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
	"github.com/brentp/irelate/interfaces"
	"github.com/brentp/xopen"
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
func New(path string) (*Tabix, error) {
	if !(xopen.Exists(path) && xopen.Exists(path+".tbi")) {
		return nil, fmt.Errorf("need gz file and .tbi for %s", path)
	}

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
	return t, nil
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

func (i *INFO) Get(key string) (interface{}, error) {
	ckey := C.CString(key)
	info := C.bcf_get_info(i.hdr, i.b, ckey)
	if info == nil {
		// return nil for bool (false)
		C.free(unsafe.Pointer(ckey))
		return nil, nil
	}
	val := i.get(info, ckey)
	C.free(unsafe.Pointer(ckey))
	return val, nil
}

func (i *INFO) Delete(key string) {
	ckey := C.CString(key)
	C.bcf_update_info(i.hdr, i.b, ckey, nil, 0, C.BCF_HT_STR)
	C.free(unsafe.Pointer(ckey))
}

func (i *INFO) Set(key string, ovalue interface{}) error {
	ckey := C.CString(key)
	var ret C.int

	switch ovalue.(type) {
	case int:
		value := C.int32_t(ovalue.(int))
		ret = C.ibcf_update_info_int32(i.hdr, i.b, ckey, &value, 1)
	case float32:
		value := C.float(ovalue.(float32))
		ret = C.ibcf_update_info_float(i.hdr, i.b, ckey, &value, 1)
	case float64:
		value := C.float(ovalue.(float64))
		ret = C.ibcf_update_info_float(i.hdr, i.b, ckey, &value, 1)
	case string:
		value := C.CString(ovalue.(string))
		ret = C.ibcf_update_info_string(i.hdr, i.b, ckey, value)
		C.free(unsafe.Pointer(value))
	case bool:
		value := ovalue.(bool)

		if value {
			ret = C.ibcf_update_info_flag(i.hdr, i.b, ckey, ckey, 1)
		} else {
			ret = C.ibcf_update_info_flag(i.hdr, i.b, ckey, ckey, 0)
		}
	case []uint32:
		value := ovalue.([]uint32)
		l := len(value)
		val := make([]int32, l)
		for i, v := range value {
			val[i] = int32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		ret = C.ibcf_update_info_int32(i.hdr, i.b, ckey, (*C.int32_t)(ptr), C.int(l))

	case []int:
		value := ovalue.([]int)
		l := len(value)
		val := make([]int32, l)
		for i, v := range value {
			val[i] = int32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		ret = C.ibcf_update_info_int32(i.hdr, i.b, ckey, (*C.int32_t)(ptr), C.int(l))

	case []uint64:
		value := ovalue.([]uint64)
		l := len(value)
		val := make([]int32, l)
		for i, v := range value {
			val[i] = int32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		ret = C.ibcf_update_info_int32(i.hdr, i.b, ckey, (*C.int32_t)(ptr), C.int(l))

	case []int32:
		value := ovalue.([]int32)
		l := len(value)
		val := make([]int32, l)
		for i, v := range value {
			val[i] = int32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		ret = C.ibcf_update_info_int32(i.hdr, i.b, ckey, (*C.int32_t)(ptr), C.int(l))

	case []int64:
		value := ovalue.([]int64)
		l := len(value)
		val := make([]int32, l)
		for i, v := range value {
			val[i] = int32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		ret = C.ibcf_update_info_int32(i.hdr, i.b, ckey, (*C.int32_t)(ptr), C.int(l))

	case []float64:
		value := ovalue.([]float32)
		l := len(value)
		val := make([]float32, l)
		for i, v := range value {
			val[i] = float32(v)
		}
		ptr := unsafe.Pointer(&val[0])
		ret = C.ibcf_update_info_float(i.hdr, i.b, ckey, (*C.float)(ptr), C.int(l))

	case []float32:
		val := ovalue.([]float32)
		l := len(val)
		ptr := unsafe.Pointer(&val[0])
		ret = C.ibcf_update_info_float(i.hdr, i.b, ckey, (*C.float)(ptr), C.int(l))

	case []string:
		valuestr := ovalue.([]string)
		value := C.CString(strings.Join(valuestr, ","))
		ret = C.ibcf_update_info_string(i.hdr, i.b, ckey, value)
		C.free(unsafe.Pointer(value))
	case []bool:
		// not possible
	default:
		log.Printf("multipe values for type: %T; not implemented (key: %s)\n", ovalue, key)

	}
	var e error
	if ret != 0 {
		e = fmt.Errorf("not set %s %v of type %T\n", key, ovalue, ovalue)
	}
	C.free(unsafe.Pointer(ckey))
	return e
}

func (i *INFO) Keys() []string {
	return []string{}
}

func (i *INFO) String() string {
	kstr := C.kstring_t{}
	C.format_info(i.b, i.hdr, &kstr)
	v := C.GoStringN(kstr.s, C.int(kstr.l))
	C.free(unsafe.Pointer(kstr.s))
	return v
}

func (self *INFO) get(info *C.bcf_info_t, tag *C.char) interface{} {
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
			return _string(info, self.hdr, tag)
		}

	}

	if int(ctype) == BCF_BT_CHAR {
		return _string(info, self.hdr, tag)
	}

	// divide this number by the number of bytes per object
	n_bytes := int(info.vptr_len)

	switch ctype {
	case C.BCF_BT_INT8:
		l := n_bytes
		out := make([]int, l)
		slice := (*[1 << 20]C.int8_t)(unsafe.Pointer(info.vptr))[:l:l]
		for i := 0; i < l; i++ {
			//if slice[i] == C.bcf_int8_missing
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
			//if slice[i] == C.bcf_int16_missing
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
			//if slice[i] == C.bcf_int32_missing
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

func _string(info *C.bcf_info_t, hdr *C.bcf_hdr_t, tag *C.char) interface{} {
	s, l := info.vptr, info.vptr_len
	// in here, need to tell if it's a flag and return a bool or a string if not.
	tag_id := C.bcf_hdr_id2int(hdr, C.BCF_DT_ID, tag)
	htype := C.ibcf_hdr_id2type(hdr, C.BCF_HL_INFO, tag_id)
	if htype == C.BCF_HT_FLAG {
		return true
	}

	// https://github.com/golang/go/wiki/cgo
	slice := (*[1 << 20]uint8)(unsafe.Pointer(s))[:l:l]
	if slice[0] == 0x7 {
		return nil
	}
	return string(slice)
}

type Variant struct {
	v       *C.bcf1_t
	hdr     *C.bcf_hdr_t
	source  uint32
	Info_   interfaces.Info
	related []interfaces.Relatable
	Pos     uint64
	_alt    []string
}

func NewVariant(v *C.bcf1_t, hdr *C.bcf_hdr_t, source uint32) *Variant {
	C.bcf_unpack(v, 1|2|4) // dont unpack genotypes
	c := &Variant{v: v, hdr: hdr, source: 1}
	c.Pos = uint64(v.pos + 1)
	c.Info_ = &INFO{hdr, v}
	runtime.SetFinalizer(c, variantFinalizer)
	return c
}

func variantFinalizer(v *Variant) {
	C.bcf_destroy(v.v)
}

func (c *Variant) Id() string {
	return C.GoString(C.vid(c.v))
}

func (v *Variant) String() string {
	return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom(), v.Start()+1, v.Id(), v.Info())
}

func (c *Variant) Chrom() string {
	return C.GoString(C.bcf_hdr_id2name(c.hdr, C.int(c.v.rid)))
}

func (c *Variant) Start() uint32 {
	return uint32(c.v.pos)
}

func (c *Variant) End() uint32 {
	return uint32(c.v.pos + c.v.rlen)
}

func (c *Variant) CIPos() (uint32, uint32, bool) {
	s := c.Start()
	ipair, err := c.Info_.Get("CIPOS")
	if ipair == nil || err != nil {
		return s, s + 1, false
	}
	pair, ok := ipair.([]interface{})
	if !ok {
		return s, s + 1, false
	}
	left := pair[0].(int)
	right := pair[0].(int)
	return uint32(int(s) + left), uint32(int(s) + right + 1), true
}

func (c *Variant) CIEnd() (uint32, uint32, bool) {
	e := c.End()
	ipair, err := c.Info_.Get("CIEnd")
	if ipair == nil || err != nil {
		return e - 1, e, false
	}
	pair, ok := ipair.([]interface{})
	if !ok {
		return e - 1, e, false
	}
	left := pair[0].(int)
	right := pair[0].(int)
	return uint32(int(e) + left), uint32(int(e) + right + 1), true
}

func (c *Variant) Source() uint32 {
	return c.source
}

func (v *Variant) Info() interfaces.Info {
	return v.Info_
}

func (c *Variant) SetSource(s uint32) {
	c.source = s
}

func (c *Variant) Related() []interfaces.Relatable {
	return c.related
}

func (c *Variant) AddRelated(i interfaces.Relatable) {
	if c.related == nil {
		c.related = make([]interfaces.Relatable, 0)
	}
	c.related = append(c.related, i)
}

func (c *Variant) Alt() []string {
	if c._alt == nil {
		c._alt = make([]string, c.v.n_allele-1)
		for i := 1; i < int(c.v.n_allele); i++ {
			c._alt[i-1] = C.GoString(C.allele_i(c.v, C.int(i)))
		}
	}
	return c._alt
}

func (c *Variant) Ref() string {
	return C.GoString(C.allele_i(c.v, 0))
}

func (t *Tabix) Get(q interfaces.IPosition) []interfaces.IPosition {
	//cs := C.CString(fmt.Sprintf("%s:%d-%d", q.Chrom(), q.Start()+1, q.End()))
	ch := C.CString(q.Chrom())
	cs := C.tbx_name2id(t.tbx, ch)

	itr := C.tabix_itr_queryi(t.tbx, cs, C.int(q.Start()), C.int(q.End()))

	kstr := C.kstring_t{}
	overlaps := make([]interfaces.IPosition, 0, 4)
	l := C.int(10)
	for l > 0 {
		l := C.atbx_itr_next(t.htf, t.tbx, itr, &kstr)
		if l < 0 {
			break
		}
		//res := C.GoString(kstr.s)
		if t.typ == VCF {
			b := C.bcf_init()
			ret := C.vcf_parse(&kstr, t.hdr, b)
			if ret < 0 {
				log.Printf("error parsing %s\n", C.GoStringN(kstr.s, C.int(kstr.l)))
			}
			overlaps = append(overlaps, NewVariant(b, t.hdr, 1))
		} else if t.typ == BED {
			iv, err := irelate.IntervalFromBedLine(C.GoStringN(kstr.s, C.int(kstr.l)))
			if err != nil {
				log.Printf("error parsing %s:%s\n", C.GoStringN(kstr.s, C.int(kstr.l)), err)
			}
			if iv != nil {
				overlaps = append(overlaps, iv)
			}
		}

	}
	C.hts_itr_destroy(itr)
	C.free(unsafe.Pointer(ch))
	C.free(unsafe.Pointer(kstr.s))
	return overlaps
}

// At takes a region like 1:45678-56789 and returns a channel on which
// it sends a string for each line that falls in that interval.
func (t *Tabix) At(region string) irelate.RelatableChannel {
	cs := C.CString(region)
	defer C.free(unsafe.Pointer(cs))

	out := make(irelate.RelatableChannel, 20)
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
			out <- NewVariant(b, t.hdr, 1)
		}
		close(out)
		C.hts_itr_destroy(itr)
		C.free(unsafe.Pointer(kstr.s))
	}()
	return out
}
