# Re: `Cyperaceae-all-taxa-6-calib-comb-29JUL-thinned.tree`

_Ruan van Mazijk, 2020_

The thinned tree-set was created using shell commands. We needed to have a thinned tree-set because the unthinned tree-set's NEXUS-file is 4.54 GB in size (`Cyperaceae-all-taxa-6-calib-comb-29JUL.trees`, posterior sample output from BEAST, not in this repository), making it intractable to further analyse.

The thinned tree-set is a random sample of 100 trees from the trees sampled from the posterior during the BEAST reconstruction. To make this thinned tree-set, I did the following:

```sh
# Count the number of trees in the unthinned set
grep -oi "tree STATE_" Cyperaceae-all-taxa-6-calib-comb-29JUL.trees | wc -l
## 275,114 trees
```

Then, after looking around the file using `head` and `tail`, I determined that the NEXUS-file has 355 lines of header content (taxon names, etc.) before the tree-set started. Bearing this in mind, the single `End;`-line at the end of the file, and the file's total number of lines at 275,470, I pulled out the `tree STATE_`-lines as follows:

```sh
cat Cyperaceae-all-taxa-6-calib-comb-29JUL.trees
  | head -n 275469  # (= 275,470 - 1)                    # Remove the last line
  | tail -n 275114  # (= 275,469 - 355)                  # Remove the header lines
  | gshuf -n 100                                         # Sample 100 random lines
  > Cyperaceae-all-taxa-6-calib-comb-29JUL-thinned.tree  # Save the output
```

(Note, `gshuf` is a macOS "coreutils" function.)

I then manually added the header-lines and the `End;`-line to `Cyperaceae-all-taxa-6-calib-comb-29JUL-thinned.tree`.
