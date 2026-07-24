# Web CAS

A **client-side Computer Algebra System** with Xcas-like syntax, built in vanilla JavaScript.
No backend, no API keys — everything runs in the browser.

[Live Demo](https://hendr15k.github.io/bookish-waffle/) · [Android APK (Nightly)](https://github.com/hendr15k/bookish-waffle/releases/tag/android-nightly) · [Deutsche Anleitung](ANLEITUNG.md)

---

## Quick Start

```bash
git clone https://github.com/hendr15k/bookish-waffle.git
cd bookish-waffle
# open index.html in any modern browser — done
```

Or just visit the [live demo](https://hendr15k.github.io/bookish-waffle/).

```
diff(sin(x^2), x)              →  2x·cos(x²)
integrate(x * e^x, x)          →  x·eˣ − eˣ
solve(x^2 - 4 = 0, x)          →  {2, −2}
det([[1,2],[3,4]])             →  −2
plot(sin(x), x, -10, 10)       →  canvas plot
limit(sin(x)/x, x, 0)          →  1
laplace(t^2, t, s)             →  2/s³
sin(pi/6)                      →  1/2
factor(x^3 - 1)                →  (x−1)(x²+x+1)
simplify(sin(x)/cos(x))        →  tan(x)
```

---

## Features

**15 domains**, ~40k lines of code, zero runtime dependencies (MathJax for LaTeX rendering only).

| Domain | Highlights |
|--------|-----------|
| Algebra | expand, simplify, factor, solve (linear/quadratic/polynomial/inequalities), partfrac, Sturm sequences |
| Calculus | diff, integrate (definite/indefinite), limits, Taylor/Laurent/Padé, sums, products, ODE (desolve, rk4/rk45) |
| Transforms | Fourier series, FFT, Laplace (forward/inverse), Z-transform |
| Linear Algebra | det, inv, eigenvalues/eigenvectors, LU/QR/Cholesky/SVD, rref, rank, kernel, Gram-Schmidt, matrix exp |
| Statistics | descriptive stats, 12 distributions, hypothesis tests (Z/T/χ²), regression (linear/poly/exp), ANOVA |
| Number Theory | factorization, primality, totient, CRT, modular arithmetic, continued fractions, Legendre/Jacobi |
| Plotting | 2D/3D/parametric/polar/implicit plots, vector fields, pan & zoom, distribution plots |
| Optimization | simplex, gradient descent, Lagrange multipliers, Euler-Lagrange |
| Control Theory | Routh-Hurwitz, transfer functions, Bode plots |
| Finance | NPV, IRR, compound interest, annuities, amortization, Black-Scholes |
| Physics | kinematics, circuits (DC/AC), Coulomb's law |
| Chemistry | molar mass, equation balancing, CAS number validation |
| Geometry | shapes, analytic geometry, plane equations |
| Graph Theory | pathfinding, connectivity |
| Logic | truth tables, CNF/DNF, Boolean simplification |

### UI

- **Desktop mode** — command-line with history, sessions, search & export
- **App mode** — grid-based tool interface for mobile/tablet
- **Worksheet UX** — numbered entries, timestamps, ↻ rerun, live syntax check, scroll-to-newest pill
- **Symbol palette** — one-tap math symbols and function shortcuts
- **Dark mode** — persisted toggle
- **Keyboard shortcuts** — `/` focus, `Esc` clear/close, `Ctrl+Z` input undo, `Ctrl+L` clear worksheet
- **Offline capable** — everything runs locally, no network needed

---

## Development

```bash
npm install                # dev dependencies (jsdom, playwright)
npm test                   # run all 56 test files
npm run test:core          # core engine tests
npm run test:latex         # LaTeX rendering
npm run test:chemistry     # chemistry module
npm run test:ui            # UI logic (jsdom)
```

### Architecture

| File | Lines | Role |
|------|------:|------|
| `js/cas.js` | 17,042 | CAS engine: evaluation, algorithms, 300+ commands |
| `js/expression.js` | 5,402 | AST nodes (`Add`, `Mul`, `Call`, …), simplification rules |
| `js/parser.js` | 1,086 | Recursive-descent parser → AST |
| `js/help.js` | 809 | Help system, command docs |
| `js/chemistry.js` | 233 | Molar mass, equation balancing |
| `index.html` | 10,887 | Full UI (worksheet, app mode, tools, plotting) |
| `tests/` | 4,961 | 56 Node.js test files |

---

## Android App

Capacitor-based wrapper. Nightly APKs are built via GitHub Actions on every push to `main`.

**Download:** [Nightly release → `app-debug.apk`](https://github.com/hendr15k/bookish-waffle/releases/tag/android-nightly)

```bash
npm run build:android:web   # prepare dist-android/
npm run cap:sync            # sync into android/
```

### CI/CD

- Push to `main` → updates **Android Nightly** prerelease
- Tag `v*` → creates a stable GitHub Release
- Release assets: debug APK, unsigned release APK, signed APK (if secrets configured)

### Signing secrets (optional)

`ANDROID_KEYSTORE_BASE64`, `ANDROID_KEYSTORE_PASSWORD`, `ANDROID_KEY_ALIAS`, `ANDROID_KEY_PASSWORD`

---

## License

MIT
