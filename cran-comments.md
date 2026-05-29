# CRAN submission comments — SWORD 0.1.0

## R CMD check results

**0 errors | 0 warnings | 1 note (CRAN servers)**

On CRAN check machines only the standard "New submission" NOTE is expected.

### NOTE 1 — expected for a first submission (all platforms)

```
* checking CRAN incoming feasibility ... NOTE
  New submission
```

### Local-only NOTEs (do NOT appear on CRAN servers)

```
* checking top-level files ... NOTE
  Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.
```
Pandoc is not installed in the local check environment. Confirmed absent on R-hub.

```
* checking for future file timestamps ... NOTE
  unable to verify current time
```
No NTP access in the local check environment. Confirmed absent on R-hub.

---

## Test environments

| Platform | R version | Errors | Warnings | Notes |
|---|---|---|---|---|
| Windows 11 (local) | R 4.5.3 ucrt | 0 | 0 | 1–3 (see above) |
| Linux R-devel Ubuntu 24.04 (rhub) | R-devel | 0 | 0 | 1 (New submission) |
| macOS ARM Sequoia 15 (rhub) | R-devel | 0 | 0 | 1 (New submission) |

---

## Reverse dependencies

None. This is a new package with no existing reverse dependencies.
