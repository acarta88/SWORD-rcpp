# CRAN submission comments — SWORD 0.1.0

## R CMD check results

**0 errors | 0 warnings | 1 note**

The single NOTE is:

```
* checking CRAN incoming feasibility ... NOTE
  New submission
```

This is expected for a first submission.

---

## Test environments

| Platform | R version | Errors | Warnings | Notes |
|---|---|---|---|---|
| Windows 11 (local) | R 4.5.3 | 0 | 0 | 1 (New submission) |
| Linux R-devel Ubuntu 24.04 (rhub) | R-devel | 0 | 0 | 1 (New submission) |
| macOS ARM Sequoia 15 (rhub) | R-devel | 0 | 0 | 1 (New submission) |

Notes on the local Windows run that do NOT appear on CRAN servers:

- *"unable to verify current time"* — no NTP access in the check environment.
- *"README.md or NEWS.md cannot be checked without 'pandoc'"* — pandoc is not
  installed locally; CRAN check machines have pandoc.

---

## Reverse dependencies

None. This is a new package with no existing reverse dependencies.
