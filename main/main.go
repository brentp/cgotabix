package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"strconv"
	"strings"
	"time"

	"github.com/brentp/cgotabix"
	"github.com/brentp/irelate"
)

func benchmarkTabix(other string, ntimes int) {

	start := time.Now()
	tbxs := make([]cgotabix.Tabix, ntimes)

	f := "/usr/local/src/gemini_install/data/gemini/data/ExAC.r0.3.sites.vep.tidy.vcf.gz"
	for i := 0; i < ntimes; i++ {
		tbxs[i] = cgotabix.New(f)
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

		var r string
		for _, tbx := range tbxs {
			for r = range tbx.Get(qstr) {
				sp := strings.Split(r, "\t")
				start, err := strconv.Atoi(sp[1])
				ref := sp[2]
				alt := strings.Split(sp[3], ",")
				if ref == "" {
					log.Fatal(err)
				}

				if err != nil {
					log.Fatal(err)
				}
				n += 1
				fmt.Fprintf(out, "%d\t%v\t%s\n", start, alt, r)
			}
			r += ""
			k += 1
		}

	}
	out.Flush()
	log.Printf("tabix\t%d\t%s\t%.3f", ntimes, other, time.Since(start).Seconds())

}

func benchmarkIrelate(other string, ntimes int) {

	start := time.Now()
	f := "/usr/local/src/gemini_install/data/gemini/data/ExAC.r0.3.sites.vep.tidy.vcf.gz"
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

	for v := range irelate.IRelate(irelate.CheckRelatedByOverlap, 0, irelate.Less, relatables...) {
		qstr := fmt.Sprintf("%s:%d-%d", v.Chrom(), v.Start(), v.End())
		fmt.Fprintf(out, "%s\n", qstr)

	}
	out.Flush()
	log.Printf("irelate\t%d\t%s\t%.3f", ntimes, other, time.Since(start).Seconds())

}

func main() {

	for _, n_query := range []int{10000, 100000} {
		for _, n_db := range []int{1, 2, 4, 8} {
			f := fmt.Sprintf("intervals.%d.vcf", n_query)
			benchmarkTabix(f, n_db)
			//benchmarkIrelate(f, n_db)
		}
	}
}
