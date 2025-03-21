<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

## debugging

Run tool with `--verbose &> debug.log` to identify anomalous reads. Check read in SAM, FASTQ, input FASTQs, and log to identify issues.

```bash
READ=SRR28867558.7843
grep -m1 ${READ} WTv2/*sam
grep -m1 -A1 ${READ} ${R2}
grep -m1 -A300 ${READ} debug.log | less
```

