package cgotabix

import (
	"testing"

	. "gopkg.in/check.v1"
)

func Test(t *testing.T) { TestingT(t) }

type TSuite struct{}

var _ = Suite(&TSuite{})

func (s *TSuite) TestRead(c *C) {
	t := New("vt.norm.vcf.gz")
	i := 0
	for _ = range t.Get("1:50000-90000") {
		i += 1
	}
	c.Assert(i, Equals, 15)

}

func (s *TSuite) TestEmptyRegion(c *C) {
	t := New("vt.norm.vcf.gz")
	i := 0
	for _ = range t.Get("2:50000-90000") {
		i += 1
	}
	c.Assert(i, Equals, 0)

}
