# Web CAS (Client-Side)

A lightweight, **client-side Computer Algebra System** with Xcas-like syntax, built with vanilla JavaScript. No backend, no API keys — runs entirely in your browser.

[**🧮 Live Demo**](https://hendr15k.github.io/bookish-waffle/) · [**📱 Android APK (Nightly Release)**](https://github.com/hendr15k/bookish-waffle/releases/tag/android-nightly) · [📖 Deutsche Anleitung](ANLEITUNG.md)

---

## Highlights

- **15 functional domains** — from basic arithmetic to control theory
- **~34,000 lines of code** — 23k JS engine + 7.5k UI
- **40 test files** — all passing
- **Zero dependencies at runtime** — only MathJax (LaTeX rendering)
- **Dark Mode** — easy on the eyes
- **Mobile App Mode** — grid-based tool interface for touch devices
- **Offline capable** — everything runs locally
- **Exact trigonometric values** — sin(π/6)=1/2, cos(π/3)=1/2, tan(π/4)=1, etc.
- **Smart algebraic simplification** — like-term collection, negative exponents as fractions

---

## Features by Domain

### 🔢 Arithmetic & Algebra
* Basic operations: `+`, `-`, `*`, `/`, `^`
* Symbolic computation: variables, simplification, expansion
* `expand(expr)`, `simplify(expr)`, `factor(n)`, `gcd(a, b)`, `lcm(a, b)`, `factorial(n)`
* Equation solving: `solve(equation, variable)` — linear, quadratic, polynomial systems
* Inequality solving: `solve(inequality, variable)`
* Nonlinear system solver
* Partial fractions: `partfrac(expr, var)`
* Complete the square, coefficient extraction, polynomial degree, Sturm sequences
* Set operations

### 📐 Calculus
* **Differentiation**: `diff(expr, var)` — symbolic, multi-order
* **Integration**: `integrate(expr, var)` (indefinite) and `integrate(expr, var,lower,upper)` (definite)
  - By parts, substitution, partial fractions, special functions
* **Limits**: `limit(expr, var, point)` — one-sided, infinite, trigonometric
* **Taylor Series**: `taylor(expr, var, point, order)`
* **Sums & Products**: `sum(expr, var, start, end)`, `product(expr, var, start, end)`
* **ODE Solver**: `ode_solver(diff(y,x) = ..., x, y, x0, y0)`
* **Riemann Sums** with visualization
* **N-step numerical integration**
* **Fourier Series**: `fourierSeries(f, var, point, order)`
* **Laplace Transforms**: `laplace(expr, var, s)`, `invLaplace(expr, s, var)`
* **Transfer functions**
* **Z-Transforms**: `ztrans(expr, var, z)`
* **Ei/Li functions** (Exponential and Logarithmic integrals)
* **Special integrals**: erf, erfc, erfi, Fresnel, Airy, Gamma, Polygamma
* **Vector calculus**: Laplacian, Jacobian, Hessian, scalar potentials, conservative-field tests

### 📊 Linear Algebra
* Vectors & matrices: `[a, b]` or `[[a,b],[c,d]]`
* Operations: multiply, dot product, cross product
* `det(M)`, `inv(M)`, `trans(M)`, `eigenvalues(M)`, `eigenvectors(M)` — exact on small matrices; larger ones can take noticeably longer
* Matrix decompositions: LU, QR, Cholesky, SVD
* `rref(M)` (Row Reduced Echelon Form), `rank(M)`, `kernel(M)`
* `adj(M)` (Adjugate), matrix powers: `pow(M, n)`
* Characteristic polynomial, pseudoinverse, condition number, nullity
* Solve linear systems
* Gram-Schmidt orthogonalization
* Chinese Remainder Theorem: `crt()`
* Lagrange interpolation
* RK4 (Runge-Kutta 4th order)
* Least squares: `lsq()`
* Convolution

### 📈 Statistics & Probability
* Descriptive: `mean(list)`, `variance(list)`, `stdDev(list)`
* Linear, polynomial, exponential regression
* **Distributions**: Normal, Binomial, Poisson, t, F, Chi-squared
* **Hypothesis Tests**: Z-Test, T-Test
* **Expected values** for probability distributions
* **Binomial distribution** visualization
* Covariance, correlation, confidence intervals, entropy, ANOVA
* Regression with constants support

### 🔢 Number Theory
* `factor(n)`, `gcd(a, b)`, `lcm(a, b)`
* Prime factorization, primality testing
* Next/previous primes, modular inverse, modular powers
* Divisors, prime factors, Fibonacci, divisor sums
* Euler's totient: `totient(n)`
* Möbius and sigma arithmetic functions
* `coprime(a, b)`, `is_coprime()`
* Integer partitions: `partitions(n)`
* Egyptian fractions: `egyptian_fraction(n)`

### 📉 Plotting & Visualization
* **2D plots**: `plot(expr, var, min, max)` — with canvas rendering
* **Curve analysis tools**: arc length, surface area, curvature
* **Multi-plot** support
* **Parametric plots** with dedicated UI tab
* **Polar plots** with dedicated UI tab
* **Implicit plots**
* **Vector fields**
* **3D plots**: `plot3d(f, var1, var2, min1, max1, min2, max2)`
* **Interactive plotting** with buttons
* **Distribution plots**
* N-step plotting with configurable rectangles
* Plot pan/zoom and interactive inspection

### 💰 Finance
* Time Value of Money calculations
* NPV, IRR
* Financial equation solvers

### ⚡ Physics
* Kinematics calculations
* Electricity (Ohm's law, circuits)
* Physics equation solver
* **AC Circuit Analysis**
* **Coulomb's Law**

### 🧪 Chemistry
* Molar mass calculator
* Equation balancing
* CAS Registry Number validation

### 📐 Geometry
* Shape calculations
* Circle properties
* Distance calculations
* **Analytic geometry** tools
* Plane equations from three points

### 🔬 Fourier Analysis
* **FFT** (Fast Fourier Transform)
* Fourier series expansion
* Fresnel integrals and cyclic integration

### 🎛️ Laplace Transform
* Forward and inverse Laplace transforms
* Transfer function analysis
* Convolution operations

### 📊 Optimization
* Simplex algorithm
* Gradient descent: `gradient_descent()`

### 🎛️ Control Theory
* **Routh-Hurwitz stability criterion**
* Transfer functions
* **Bode plots** (phase and magnitude)

### 🕸️ Graph Theory
* Pathfinding algorithms
* Connectivity analysis

### 🧠 Logic
* Truth table generation
* Boolean expression simplification
* CNF and DNF conversion
* Operators: `and`, `or`, `not`, `xor`

### 🎹 Interactive Tools
* **Desktop mode**: Command-line interface with history
* **App mode**: Grid-based tool interface for mobile/tablet
* **Dark mode** toggle
* **Keyboard shortcuts**
* **Variable definitions** stored across calculations
* **Command history** persisted via localStorage
* **Named sessions** with quick save/load
* **History search and export**
* Dedicated apps for complex numbers and special functions

---

## Output

* Results rendered in **LaTeX** using MathJax for beautiful mathematical display
* Canvas-based function plotting
* Command history saved locally in browser

---

## Usage

1. Open [the live demo](https://hendr15k.github.io/bookish-waffle/) in any modern browser
2. Or clone and open `index.html` locally
3. Type commands in the input field

**Examples:**
```js
diff(sin(x^2), x)                    // Differentiate → 2x·cos(x²)
integrate(x * e^x, x)                // Integrate → x·eˣ - eˣ
solve(x^2 - 4 = 0, x)                // Solve → {2, -2}
[[1, 2], [3, 4]] * [x, y]            // Matrix-vector multiply
plot(sin(x), x, -10, 10)             // Plot function
limit(sin(x)/x, x, 0)                // Limit → 1
laplace(t^2, t, s)                   // Laplace → 2/s³
fourierSeries(x^2, x, 0, 5)          // Fourier series
det([[1,2],[3,4]])                   // Determinant → -2
sin(pi/6)                            // Exact value → 1/2
cos(pi/3)                            // Exact value → 1/2
atan(Infinity)                       // → π/2
factor(x^3 - 1)                      // → (x-1)(x²+x+1)
expand((x+1)*(x-1))                  // → x²-1
simplify(sin(x)/cos(x))              // → tan(x)
x^(-2)                               // → 1/x²
2*x + 3*x                            // → 5x
```

---

## Recent Improvements (April 2026)

### Bug Fixes & Enhancements
- **Exact trigonometric values** — sin/cos/tan of π/6, π/4, π/3 multiples (all 12 quadrants)
- **atan(∞) = π/2** — infinite limits for inverse trig
- **Factor over ℝ** — irreducible quadratics like x²+1 stay unfactored (no complex factors)
- **Like-term collection** — 2x+3x→5x, 5x-3x→2x in Add/Sub.simplify
- **Negative exponents** — x⁻²→1/x² automatically
- **Trig ratio simplification** — sin/cos→tan, 1/sin→csc, cos/sin→cot
- **LaTeX improvements** — scientific notation (1×10⁻¹⁰), negated fractions (-1/x)
- **Expand + simplify** — expand() now calls simplify() after distribution
- **Parser roundtrip** — negative numbers parse correctly (no more Mul(-1,n) bloat)
- **Jekyll compatibility** — .nojekyll added to prevent Liquid template errors

### Bug Reports (from automated analysis)
- `BUGS_PARSER.md` — 7 parser/LaTeX bugs identified
- `BUGS_SIMPLIFY.md` — 22 simplification bugs identified
- `BUGS_CALCULUS.md` — calculus edge cases documented

## Development & Testing

```bash
npm install              # Install dev dependencies (jsdom for UI tests)
npm test                 # Run ALL tests (40 files)
npm run test:core        # Run only core tests
npm run test:latex       # Run LaTeX rendering tests
npm run test:chemistry   # Run chemistry tests
```

### Architecture

| File | Lines | Purpose |
|------|-------|---------|
| `js/parser.js` | 1,032 | Lexer and Parser (Recursive Descent) → AST |
| `js/expression.js` | 4,933 | AST node definitions (`Expr`, `Add`, `Mul`, `Call`, etc.) |
| `js/cas.js` | 16,077 | Core CAS engine: evaluation, substitution, algorithms |
| `js/help.js` | 807 | Help system and command documentation |
| `js/chemistry.js` | 230 | Chemistry module (molar mass, balancing) |
| `index.html` | ~7,500 | Frontend UI (command-line + app mode) |
| `tests/` | 3,169 | 40 Node.js test files |

---

## Project Stats

| Metric | Value |
|--------|-------|
| Total lines of code | ~34,000 |
| JavaScript engine | 23,049 lines |
| Frontend | 7,466 lines |
| Test files | 40 |
| Test lines | 3,169 |
| Feature domains | 15 |
| Live since | March 2025 |
| Last major update | April 2026 |

---

## License

MIT — see LICENSE file.

## Android App + APK Releases

This repository now also contains an Android wrapper built with Capacitor.

**Direct download:**
- Nightly release page: <https://github.com/hendr15k/bookish-waffle/releases/tag/android-nightly>
- Installable asset: `app-debug.apk`

### Local commands

```bash
npm run build:android:web   # prepare dist-android/
npm run cap:sync            # sync web assets into android/
```

### GitHub Actions APK builds

- Pushes to `main` create/update the **Android Nightly** prerelease
- Tags like `v1.0.0` create a normal GitHub Release
- Local Gradle builds need a configured Android SDK (`ANDROID_HOME` or `android/local.properties`)
- GitHub Actions sets up the Android SDK automatically
- Release assets include:
  - installable debug APK
  - unsigned release APK
  - signed release APK if Android signing secrets are configured

### Optional signing secrets

- `ANDROID_KEYSTORE_BASE64`
- `ANDROID_KEYSTORE_PASSWORD`
- `ANDROID_KEY_ALIAS`
- `ANDROID_KEY_PASSWORD`

