## quick and diry tabix wrapper for go using cgo.

```go
tbx := New("vt.norm.vcf.gz")
for str := range t.Get("1:50000-90000") {
	fmt.Println(str)
}
```
