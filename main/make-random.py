from random import randint



def reservoir_sampler(iterable, n):
    reservoir = []
    for i in range(n):
        reservoir.append((i, next(iterable)))
    for i, item in enumerate(iterable, start=n):
        pos = randint(0, i)
        if pos < n:
            reservoir[pos] = (i, item)
    return [x[1] for x in sorted(reservoir)]


if __name__ == "__main__":
    import sys
    import toolshed as ts

    n = int(sys.argv[1])
    fh = ts.nopen(sys.argv[2])
    for line in fh:
        print line,
        if line.startswith("#CHROM"): break
    print "".join(reservoir_sampler(fh, n))
