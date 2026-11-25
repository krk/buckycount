# buckycount

Count fullerene cages and reconstruct them by rank.

The preprint series [*Counting Classical Fullerenes*](https://doi.org/10.5281/zenodo.24289727) develops the counting formulas used by buckycount. Paper VII, *Arithmetic Unranking of Fullerene Isomers*, defines the rank order and shows how to reconstruct a cage from its rank.

- `make` builds `./buckycount`.
- `make install` installs to `/usr/local/bin`. Use `make install PREFIX="$HOME/.local"` for a user installation.
- `make clean` removes the executable.

`./buckycount --help` prints:

```text
Usage:
  buckycount [count] N|A..B [options]
  buckycount unrank N R|A..B|all|- -o FILE.pc [options]
  buckycount unrank N cache [options]
  buckycount decode FILE.pc [FILE.pc ...] [options]

N is an even vertex count. Ranks start at 0. Ranges include both endpoints.

Count options:
  --ipr              Count only isolated-pentagon cages
  --details          Show the count breakdown
  --format csv|json  Write counts as CSV or JSON
  --families         List families and their rank intervals for one vertex count

Unrank options:
  --family A,B       Use ranks within family A,B when emitting graphs

Count, unrank, and decode options:
  -o, --output PATH  Count or decode to PATH (directory for multiple inputs). Required for graph output.
  --workers W        Limit workers for count, unrank all, or multi-file decode
  --no-estimates     Skip size estimates
  -h, --help         Show help for the program or the selected command

Count and unrank options:
  --cache-dir DIR    Set the count or unrank cache directory (default: .cache)
  --no-cache         Disable caching for count or unrank (except unrank N cache)

Counting formulas and rank order: https://doi.org/10.5281/zenodo.24289727
```

Graph output uses [planar code](https://caagt.ugent.be/CaGe/Archive/codes.html).

Use `unrank N cache` to prepare the rank index for repeated lookups, or `--no-cache` for uncached benchmarks.

The buckycount source code is available on [GitHub](https://github.com/krk/buckycount).

This distribution is released under the [Apache-2.0 license](LICENSE).
