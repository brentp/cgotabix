## quick and diry tabix wrapper for go using cgo.

```go
package main

import (
	"fmt"

	"github.com/brentp/cgotabix"
)

func main() {
	tbx := cgotabix.New("vt.norm.vcf.gz")
	for str := range tbx.Get("1:50000-90000") {
		fmt.Println(str)
	}
}
```
