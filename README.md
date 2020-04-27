# ob_utils
Python programs for analyzing onebreak results.

## Dependency

### Python
Python (>=3.7)

### Software
[bedtools](http://bedtools.readthedocs.org/en/latest/)
[sv_utils](https://github.com/friend1ws/sv_utils)(>=0.6.1)

## Install

```
python3 setup.py install
```

## Commands

```
ob_utils manta_sv [-h] --in_manta_sv IN_MANTA_SV --output OUTPUT
                         [--margin MARGIN] [--f_grc]
```
```
ob_utils gridss_sv [-h] --in_gridss_sv IN_GRIDSS_SV --output OUTPUT
                          [--margin MARGIN] [--f_grc]
```
```
ob_utils svaba_sv [-h] --in_svaba_sv IN_SVABA_SV --in_svaba_indel
                         IN_SVABA_INDEL --output OUTPUT [--margin MARGIN]
                         [--f_grc] [--f_germ]
                         [--normal_max_variant NORMAL_MAX_VARIANT]
                         [--tumor_min_variant TUMOR_MIN_VARIANT]
                         [--normal_min_depth NORMAL_MIN_DEPTH]
                         [--tumor_min_depth TUMOR_MIN_DEPTH]
                         [--min_del_size MIN_DEL_SIZE]
                         [--min_ins_size MIN_INS_SIZE]
```
```
ob_utils genomon_sv [-h] --in_genomon_sv IN_GENOMON_SV --output OUTPUT
                           [--margin MARGIN] [--f_grc]
```
```
ob_utils comp [-h] --in_onebreak IN_ONEBREAK --in_genomonsv
                 IN_GENOMONSV --output OUTPUT [--margin MARGIN]
```
```
ob_utils merge [-h] --in_onebreak_filt1 IN_ONEBREAK_FILT1
                 --in_onebreak_filt2 IN_ONEBREAK_FILT2 --output OUTPUT
```

You can check the manual by typing
```
ob_utils comp -h
```
```
ob_utils merge -h
```

## Results

The primary result is {OUTPUT}

    GenomonSV: Matched GenomonSV result
