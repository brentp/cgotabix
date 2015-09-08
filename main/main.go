package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"runtime/pprof"
	"time"

	"github.com/biogo/hts/internal"
	"github.com/brentp/bix"
	"github.com/brentp/cgotabix"
	"github.com/brentp/irelate"
	"github.com/brentp/irelate/interfaces"
	"github.com/brentp/xopen"
)

//	"/media/brentp/ffd28ae3-3dca-4f44-8aed-1685a55661f8/uw-course/thu/data/data-diseaseX/disease_x.merged.jointCalled.vcf.ann.Recalibrated.Merged.HC.vcf.VT.vep.vcf.gz",
var FS = []string{
	"/usr/local/src/gemini_install/data/gemini/data/clinvar_20150305.tidy.vcf.gz",
	"/usr/local/src/gemini_install/data/gemini/data/ExAC.r0.3.sites.vep.tidy.vcf.gz",
	"/usr/local/src/gemini_install/data/gemini/data/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.tidy.vcf.gz",
	"/usr/local/src/gemini_install/data/gemini/data/ESP6500SI.all.snps_indels.tidy.v2.vcf.gz",
	"/usr/local/src/gemini_install/data/gemini/data/ExAC.r0.3.sites.vep.tidy.vcf.gz",
	"/usr/local/src/gemini_install/data/gemini/data/GRCh37-gms-mappability.vcf.gz",
	"/usr/local/src/gemini_install/data/gemini/data/cosmic-v68-GRCh37.tidy.vcf.gz",
	"/usr/local/src/gemini_install/data/gemini/data/dbsnp.b141.20140813.hg19.tidy.vcf.gz",
}

func benchmarkBix(other string, ntimes int) {
	start := time.Now()
	tbxs := make([]*bix.Bix, ntimes)

	var err error
	for i := 0; i < ntimes; i++ {
		tbxs[i], err = bix.New(FS[i])
		if err != nil {
			log.Fatal(err)
		}
	}
	or, err := xopen.Ropen(other)
	if err != nil {
		log.Fatal(err)
	}
	vcf := irelate.Vopen(or, nil)

	out := bufio.NewWriter(ioutil.Discard)
	/*
		out, err := vcfgo.NewWriter(os.Stdout, vcf.Header)
		if err != nil {
			log.Fatal(err)
		}
	*/

	n := 0
	k := 0
	for {
		var v interfaces.IVariant

		v = vcf.Read()
		if v == nil {
			break
		}

		for _, tbx := range tbxs {
			//for _, r := range tbx.Get(v) {
			rdr, err := tbx.Query(v.Chrom(), int(v.Start()), int(v.End()), false)
			if err != nil {
				if err != internal.ErrInvalid {
					log.Fatal(err)
				} else {
					continue
				}
			}
			db := irelate.Vopen(rdr, vcf.Header)
			db.AddInfoToHeader("XXX", "1", "Float", "XXX")
			db.AddInfoToHeader("many", "3", "Float", "XXX")
			db.AddInfoToHeader("manyi", "3", "Integer", "XXX")
			db.AddInfoToHeader("flag", "1", "Flag", "XXX")
			for {
				r := db.Read()
				if r == nil {
					break
				}
				//sp := strings.Split(r, "\t")

				if !interfaces.Same(r, v, true) {
					continue
				}

				ov := r.(interfaces.IVariant)
				n += 1
				abool := n%2 == 0
				//log.Println(ov.Info.Get("culprit"))
				//ov.Info.Set("culprit", "hi")
				ov.Info().Set("DP", 23)
				// TODO: update header then Set
				v := []float32{33.0, 33.0, 44.0}
				ov.Info().Set("many", v)
				ov.Info().Set("flag", abool)

				ov.Info().Set("manyi", []int32{22, 1, 2})
				//ov.Info.Set("XXX", 23.4)
				////log.Println(ov.Info.Get("culprit"))
				vv, err := ov.Info().Get("many")
				if err != nil {
					panic("couldn't get many")
				}
				if len(vv.([]float32)) != 3 {
					panic("bad length")
				}

				v2, err := ov.Info().Get("manyi")
				if v2.([]int)[2] != 2 {
					log.Fatalf("bad int: %v\n", v2)
				}
				b, err := ov.Info().Get("flag")
				if b == nil && abool != false {
					log.Fatal("bad bool")
				} else if abool && !b.(bool) {
					log.Fatal("bad bool")
				}

				dp, err := ov.Info().Get("DP")
				if v, ok := dp.(int); !ok {
					//log.Println(ov.Info())
				} else {
					if v != 23 {
						log.Fatal("bad depth")
					}
				}

				fmt.Fprintf(out, "%s\t%d\t%d\t%s\t%s\t%s\n", r.Chrom(), r.Start(), r.End(), ov.Ref(), ov.Alt(), ov.Info().String())
				//log.Println(r.Chrom(), r.Start(), r.End(), ov.Ref, ov.Alt)
			}
			k += 1
		}

	}
	out.Flush()
	log.Printf("gobix\t%d\t%s\t%.3f\t%d", ntimes, other, time.Since(start).Seconds(), n)

}

func benchmarkTabix(other string, ntimes int) {

	start := time.Now()
	tbxs := make([]*cgotabix.Tabix, ntimes)

	var err error
	for i := 0; i < ntimes; i++ {
		tbxs[i], err = cgotabix.New(FS[i])
		if err != nil {
			log.Fatal(err)
		}
		tbxs[i].AddInfoToHeader("XXX", "1", "Float", "XXX")
		tbxs[i].AddInfoToHeader("many", "3", "Float", "XXX")
		tbxs[i].AddInfoToHeader("manyi", "3", "Integer", "XXX")
		tbxs[i].AddInfoToHeader("flag", "1", "Flag", "XXX")
	}

	or, err := xopen.Ropen(other)
	if err != nil {
		log.Fatal(err)
	}
	vcf := irelate.Vopen(or, nil)

	out := bufio.NewWriter(ioutil.Discard)
	/* out, err := vcfgo.NewWriter(os.Stdout, vcf.Header)
	if err != nil {
		log.Fatal(err)
	} */

	n := 0
	k := 0
	for {
		var v interfaces.IVariant

		v = vcf.Read()
		if v == nil {
			break
		}

		qstr := fmt.Sprintf("%s:%d-%d", v.Chrom(), v.Start(), v.End())

		for _, tbx := range tbxs {
			//for _, r := range tbx.Get(v) {
			for r := range tbx.At(qstr) {
				//sp := strings.Split(r, "\t")

				if !interfaces.Same(r, v, true) {
					continue
				}

				ov := r.(interfaces.IVariant)
				n += 1
				abool := n%2 == 0
				//log.Println(ov.Info.Get("culprit"))
				//ov.Info.Set("culprit", "hi")
				ov.Info().Set("DP", 23)
				// TODO: update header then Set
				v := []float32{33.0, 33.0, 44.0}
				ov.Info().Set("many", v)
				ov.Info().Set("flag", abool)

				ov.Info().Set("manyi", []int32{22, 1, 2})
				//ov.Info.Set("XXX", 23.4)
				////log.Println(ov.Info.Get("culprit"))
				vv, err := ov.Info().Get("many")
				if err != nil {
					panic("couldn't get many")
				}
				if len(vv.([]float32)) != 3 {
					panic("bad length")
				}

				v2, err := ov.Info().Get("manyi")
				if v2.([]int)[2] != 2 {
					log.Fatalf("bad int: %v\n", v2)
				}
				b, err := ov.Info().Get("flag")
				if b == nil && abool != false {
					log.Fatal("bad bool")
				} else if abool && !b.(bool) {
					log.Fatal("bad bool")
				}

				dp, err := ov.Info().Get("DP")
				if v, ok := dp.(int); !ok {
					//log.Println(ov.Info())
				} else {
					if v != 23 {
						log.Fatal("bad depth")
					}
				}

				fmt.Fprintf(out, "%s\t%d\t%d\t%s\t%s\t%s\n", r.Chrom(), r.Start(), r.End(), ov.Ref(), ov.Alt(), ov.Info().String())
				//log.Println(r.Chrom(), r.Start(), r.End(), ov.Ref, ov.Alt)
			}
			k += 1
		}

	}
	out.Flush()
	log.Printf("tabix\t%d\t%s\t%.3f\t%d", ntimes, other, time.Since(start).Seconds(), n)

}

func benchmarkIrelate(other string, ntimes int) {

	start := time.Now()
	relatables := make([]irelate.RelatableChannel, ntimes+1)
	for i := 0; i < ntimes; i++ {
		db, err := irelate.Streamer(FS[i], "")
		if err != nil {
			log.Fatal(err)
		}
		relatables[i+1] = db
	}

	vcf, err := irelate.Streamer(other, "")
	if err != nil {
		log.Fatal(err)
	}
	out := bufio.NewWriter(ioutil.Discard)
	relatables[0] = vcf
	n := 0

	for v := range irelate.IRelate(irelate.CheckRelatedByOverlap, 0, irelate.NaturalLessPrefix, relatables...) {
		qstr := fmt.Sprintf("%s:%d-%d", v.Chrom(), v.Start(), v.End())
		fmt.Fprintf(out, "%s\n", qstr)
		for _, o := range v.Related() {
			if interfaces.Same(o, v, true) {
				n += 1
			}
		}
	}
	out.Flush()
	log.Printf("irelate\t%d\t%s\t%.3f\t%d", ntimes, other, time.Since(start).Seconds(), n)

}

func main() {

	pf, err := os.Create("cpu.pprof")
	if err != nil {
		panic(err)
	}
	pprof.StartCPUProfile(pf)
	defer pprof.StopCPUProfile()

	var f string
	for _, n_query := range []int{100, 1000, 10000, 100000} {
		//for _, n_query := range []int{100000} {
		//for _, n_db := range []int{1, 2, 4, 8} {
		for _, n_db := range []int{8} {
			if n_query == 2 {
				f = "/media/brentp/ffd28ae3-3dca-4f44-8aed-1685a55661f8/uw-course/thu/data/data-diseaseX/disease_x.merged.jointCalled.vcf.ann.Recalibrated.Merged.HC.vcf.VT.vep.vcf.gz"
			} else {
				f = fmt.Sprintf("intervals.%d.vcf", n_query)
			}
			//f = "v.vcf"
			benchmarkTabix(f, n_db)
			benchmarkBix(f, n_db)
			benchmarkIrelate(f, n_db)
		}
	}
}
