package cgotabix

import (
	"testing"

	. "gopkg.in/check.v1"
)

func Test(t *testing.T) { TestingT(t) }

type Position struct {
	chrom string
	start uint32
	end   uint32
}

func (p Position) Chrom() string {
	return p.chrom
}
func (p Position) Start() uint32 {
	return p.start
}
func (p Position) End() uint32 {
	return p.end
}

type TSuite struct{}

var _ = Suite(&TSuite{})

func (s *TSuite) TestRead(c *C) {
	t := New("vt.norm.vcf.gz")
	i := 0
	for _ = range t.At("1:50000-90000") {
		i += 1
	}
	c.Assert(i, Equals, 15)
	i = 0
	for _ = range t.Get(Position{"1", 50000, 90000}) {
		i += 1
	}
	c.Assert(i, Equals, 15)

}

func (s *TSuite) TestEmptyRegion(c *C) {
	t := New("vt.norm.vcf.gz")
	i := 0
	for _ = range t.At("2:50000-90000") {
		i += 1
	}
	c.Assert(i, Equals, 0)

	for _ = range t.Get(Position{"2", 50000, 90000}) {
		i += 1
	}
	c.Assert(i, Equals, 0)

}
