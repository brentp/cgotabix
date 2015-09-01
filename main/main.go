package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"time"

	"github.com/brentp/cgotabix"
	"github.com/brentp/irelate"
	"github.com/brentp/irelate/interfaces"
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
				if dp.(int) != 23 {
					log.Fatal("bad depth")
				}

				fmt.Fprintf(out, "%s\t%d\t%d\t%s\t%s\t%s\n", r.Chrom(), r.Start(), r.End(), ov.Ref, ov.Alt, ov.Info().String())
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
