package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"time"

	"github.com/brentp/cgotabix"
	"github.com/brentp/irelate"
)

const F = "/usr/local/src/gemini_install/data/gemini/data/ExAC.r0.3.sites.vep.tidy.vcf.gz"

func benchmarkTabix(other string, ntimes int) {

	start := time.Now()
	tbxs := make([]*cgotabix.Tabix, ntimes)

	f := F
	for i := 0; i < ntimes; i++ {
		tbxs[i] = cgotabix.New(f)
		tbxs[i].AddInfoToHeader("XXX", "1", "Float", "XXX")
		tbxs[i].AddInfoToHeader("many", "3", "Float", "XXX")
		tbxs[i].AddInfoToHeader("manyi", "3", "Integer", "XXX")
		tbxs[i].AddInfoToHeader("flag", "1", "Flag", "XXX")
	}

	vcf := irelate.Vopen(other)

	out := bufio.NewWriter(ioutil.Discard)

	n := 0
	k := 0
	for {

		v := vcf.Read()
		if v == nil {
			break
		}

		qstr := fmt.Sprintf("%s:%d-%d", v.Chrom(), v.Start(), v.End())

		for _, tbx := range tbxs {
			for r := range tbx.Get(qstr) {
				//sp := strings.Split(r, "\t")
				ov := r.(*cgotabix.CVariant)
				if v.Pos != ov.Pos {
					continue
				}
				if ov.Ref != v.Ref {
					continue
				}
				same := false
				for _, a := range ov.Alt {
					for _, b := range v.Alt {
						if a == b {
							same = true
							break
						}
					}
				}
				if !same {
					continue
				}
				n += 1
				//log.Println(ov.Info.Get("culprit"))
				//ov.Info.Set("culprit", "hi")
				//ov.Info.Set("DP", 23)
				// TODO: update header then Set
				v := []float32{33.0, 33.0, 44.0}
				ov.Info.Set("many", v)
				ov.Info.Set("flag", true)

				ov.Info.Set("manyi", []int32{22, 1, 2})
				//ov.Info.Set("XXX", 23.4)
				////log.Println(ov.Info.Get("culprit"))
				if len(ov.Info.Get("many").([]float32)) != 3 {
					panic("bad length")
				}

				v2 := ov.Info.Get("manyi").([]int)
				if v2[2] != 2 {
					log.Fatalf("bad int: %v\n", v2)
				}
				b := ov.Info.Get("flag").(bool)
				if !b {
					log.Fatalf("bad bool")
				}

				fmt.Fprintf(out, "%s\t%d\t%d\t%s\t%s\n", r.Chrom(), r.Start(), r.End(), ov.Ref, ov.Alt)
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
	f := F
	relatables := make([]irelate.RelatableChannel, ntimes+1)
	for i := 0; i < ntimes; i++ {
		exac, err := irelate.Streamer(f)
		if err != nil {
			log.Fatal(err)
		}
		relatables[i+1] = exac
	}

	vcf, err := irelate.Streamer(other)
	if err != nil {
		log.Fatal(err)
	}
	out := bufio.NewWriter(ioutil.Discard)
	relatables[0] = vcf
	n := 0

	for v := range irelate.IRelate(irelate.CheckRelatedByOverlap, 0, irelate.Less, relatables...) {
		qstr := fmt.Sprintf("%s:%d-%d", v.Chrom(), v.Start(), v.End())
		fmt.Fprintf(out, "%s\n", qstr)
		n += len(v.Related())
	}
	out.Flush()
	log.Printf("irelate\t%d\t%s\t%.3f\t%d", ntimes, other, time.Since(start).Seconds(), n)

}

func main() {

	for _, n_query := range []int{100, 10000, 100000} {
		for _, n_db := range []int{1, 2, 4, 8} {
			f := fmt.Sprintf("intervals.%d.vcf", n_query)
			benchmarkTabix(f, n_db)
			//benchmarkIrelate(f, n_db)
		}
	}
}
