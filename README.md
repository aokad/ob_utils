# ob_utils
Python programs for analyzing onebreak results.

## Dependency

### Python
Python (>=2.7)

### Software
[bedtools](http://bedtools.readthedocs.org/en/latest/)
[svtools](https://github.com/hall-lab/svtools.git)

## Install

```
python setup.py install
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
