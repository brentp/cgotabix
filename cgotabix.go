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

*/
import "C"
import (
	"log"
	"runtime"
	"strings"
	"unsafe"

	"github.com/brentp/irelate"
	"github.com/brentp/vcfgo"
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
	tbx *C.tbx_t
	htf *C.htsFile
	hdr *C.bcf_hdr_t
	itr *C.hts_itr_t
	typ FileType
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
	t := &Tabix{}
	cs := C.CString(path)
	defer C.free(unsafe.Pointer(cs))
	t.tbx = C.tbx_index_load(cs)
	mode := C.char('r')
	t.htf = C.hts_open(cs, &mode)
	runtime.SetFinalizer(t, tabixCloser)
	if strings.HasSuffix(path, ".bed.gz") {
		t.typ = BED
	} else if strings.HasSuffix(path, ".vcf.gz") {
		t.hdr = C.bcf_hdr_read(t.htf)
		t.typ = VCF
	} else {
		t.typ = OTHER
	}
	return t
}

type INFO struct {
	hdr *C.bcf_hdr_t
	b   *C.bcf1_t
}

func (i *INFO) Get(key string) interface{} {
	info := C.bcf_get_info(i.hdr, i.b, C.CString(key))
	if info == nil {
		return nil
	}
	return i.get(info)
}

func (self *INFO) get(info *C.bcf_info_t) interface{} {
	if info.len == 1 {
		// TODO: use switch
		if C.itype(info) == C.BCF_BT_INT8 {
			if C.ivi(info) == C.INT8_MIN {
				return nil
			}
			return int(C.ivi(info))
		}
	}
	return 22
}

type CVariant struct {
	v       *C.bcf1_t
	hdr     *C.bcf_hdr_t
	source  uint32
	Info    *vcfgo.InfoByte
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
