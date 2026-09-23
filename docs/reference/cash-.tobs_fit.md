# Access spOccupancy-compatible fields from tobs fits

Allows accessing spOccupancy-style fields (e.g., `$beta.samples`,
`$psi.samples`) on tobs_fit objects. Since tobs stores actual posterior
draws, this is a thin remapping layer.

## Usage

``` r
# S3 method for class 'tobs_fit'
x$name
```

## Arguments

- x:

  A `tobs_fit` object.

- name:

  Field name to access.

## Value

The requested field value.
