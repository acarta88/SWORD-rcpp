# CRAN submission comments — SWORD 0.1.0

## R CMD check results

**0 errors | 0 warnings | 2 notes**

### NOTE 1 — expected for a first submission

```
* checking CRAN incoming feasibility ... NOTE
  New submission
```

This is standard for any new CRAN package.

### NOTE 2 — local-only (no pandoc installed)

```
* checking top-level files ... NOTE
  Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.
```

Pandoc is not installed in the local check environment. CRAN check machines have
pandoc and this NOTE does not appear there (confirmed on R-hub Linux and macOS).

---

## Test environments

| Platform | R version | Errors | Warnings | Notes |
|---|---|---|---|---|
| Windows 11 (local) | R 4.5.3 ucrt | 0 | 0 | 2 (see above) |
| Linux R-devel Ubuntu 24.04 (rhub) | R-devel | 0 | 0 | 1 (New submission) |
| macOS ARM Sequoia 15 (rhub) | R-devel | 0 | 0 | 1 (New submission) |

---

## Reverse dependencies

None. This is a new package with no existing reverse dependencies.
