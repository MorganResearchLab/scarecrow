<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](root.md)

## scarecrow rake
The `rake` tool has since been internalized and no longer needs to be called directly. This tool takes a barcode whitelist and calculates all variations of each barcode for n `--mismatches`. The resulting list is structured and output in JSON format.

```bash
scarecrow rake --barcodes <BC1_whitelist.txt> --mismatches 2
```
