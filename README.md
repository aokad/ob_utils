# ob_utils
Merge SV results detected by multiple tools.

## Dependency

### Python
Python (>=3.7)

### Software
[bedtools](http://bedtools.readthedocs.org/en/latest/),
[sv_utils](https://github.com/friend1ws/sv_utils)(>=0.6.1),
pysam

## Install

```
python3 setup.py install
```

## Commands
```
ob_utils camphor_sv [-h] --in_camphor_tumor_sv IN_CAMPHOR_TUMOR_SV --output OUTPUT [--f_grc]
                    [--bcf_filter_option BCF_FILTER_OPTION] [--filter_scaffold_option]
                    [--min_tumor_support_read MIN_TUMOR_SUPPORT_READ]
                    [--min_sv_length MIN_SV_LENGTH] [--debug]
```
```
ob_utils cutesv_sv [-h] --in_cutesv_tumor_sv IN_CUTESV_TUMOR_SV --in_cutesv_control_sv
                   IN_CUTESV_CONTROL_SV --output OUTPUT [--margin MARGIN] [--f_grc]
                   [--bcf_filter_option BCF_FILTER_OPTION] [--filter_scaffold_option]
                   [--min_tumor_support_read MIN_TUMOR_SUPPORT_READ]
                   [--max_control_support_read MAX_CONTROL_SUPPORT_READ]
                   [--min_sv_length MIN_SV_LENGTH] [--debug]
```
```
ob_utils delly_sv [-h] --in_delly_tumor_sv IN_DELLY_TUMOR_SV --in_delly_control_sv
                  IN_DELLY_CONTROL_SV --output OUTPUT [--margin MARGIN] [--f_grc]
                  [--bcf_filter_option BCF_FILTER_OPTION] [--filter_scaffold_option]
                  [--min_tumor_support_read MIN_TUMOR_SUPPORT_READ]
                  [--max_control_support_read MAX_CONTROL_SUPPORT_READ]
                  [--min_sv_length MIN_SV_LENGTH] [--debug]
```
```
ob_utils genomon_sv [-h] --in_genomon_sv IN_GENOMON_SV --output OUTPUT [--margin MARGIN] [--f_grc]
                    [--v2]
```
```
ob_utils sniffles_sv [-h] --in_sniffles_tumor_sv IN_SNIFFLES_TUMOR_SV --in_sniffles_control_sv
                     IN_SNIFFLES_CONTROL_SV --output OUTPUT [--margin MARGIN] [--f_grc]
                     [--bcf_filter_option BCF_FILTER_OPTION] [--filter_scaffold_option]
                     [--min_tumor_support_read MIN_TUMOR_SUPPORT_READ]
                     [--max_control_support_read MAX_CONTROL_SUPPORT_READ]
                     [--min_sv_length MIN_SV_LENGTH] [--debug] [--sniffles2]
```
```
ob_utils svim_sv [-h] --in_svim_tumor_sv IN_SVIM_TUMOR_SV --in_svim_control_sv
                 IN_SVIM_CONTROL_SV --output OUTPUT [--margin MARGIN] [--f_grc]
                 [--bcf_filter_option BCF_FILTER_OPTION]
                 [--filter_scaffold_option]
                 [--min_tumor_support_read MIN_TUMOR_SUPPORT_READ]
                 [--max_control_support_read MAX_CONTROL_SUPPORT_READ]
                 [--min_sv_length MIN_SV_LENGTH] [--debug]
```

## Results

The primary result is {OUTPUT}

