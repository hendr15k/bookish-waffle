# Web CAS (Client-Side)

A lightweight, **client-side Computer Algebra System** with Xcas-like syntax, built with vanilla JavaScript. No backend, no API keys — runs entirely in your browser.

[**🧮 Live Demo**](https://hendr15k.github.io/bookish-waffle/) · [📖 Deutsche Anleitung](ANLEITUNG.md)

---

## Highlights

- **15 functional domains** — from basic arithmetic to control theory
- **~34,000 lines of code** — 26k JS engine + 7.5k UI
- **40 test files** — all passing
- **Zero dependencies at runtime** — only MathJax (LaTeX rendering)
- **Dark Mode** — easy on the eyes
- **Mobile App Mode** — grid-based tool interface for touch devices
- **Offline capable** — everything runs locally

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

### 📊 Linear Algebra
* Vectors & matrices: `[a, b]` or `[[a,b],[c,d]]`
* Operations: multiply, dot product, cross product
* `det(M)`, `inv(M)`, `trans(M)`, `eigenvalues(M)`, `eigenvectors(M)`
* Matrix decompositions: LU, QR, Cholesky
* `rref(M)` (Row Reduced Echelon Form), `rank(M)`
* `adj(M)` (Adjugate), matrix powers: `pow(M, n)`
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
* Regression with constants support

### 🔢 Number Theory
* `factor(n)`, `gcd(a, b)`, `lcm(a, b)`
* Prime factorization, primality testing
* Euler's totient: `totient(n)`
* `coprime(a, b)`, `is_coprime()`
* Integer partitions: `partitions(n)`
* Egyptian fractions: `egyptian_fraction(n)`

### 📉 Plotting & Visualization
* **2D plots**: `plot(expr, var, min, max)` — with canvas rendering
* **Multi-plot** support
* **Parametric plots**
* **Implicit plots**
* **Vector fields**
* **3D plots**: `plot3d(f, var1, var2, min1, max1, min2, max2)`
* **Interactive plotting** with buttons
* **Distribution plots**
* N-step plotting with configurable rectangles

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
* Operators: `and`, `or`, `not`, `xor`

### 🎹 Interactive Tools
* **Desktop mode**: Command-line interface with history
* **App mode**: Grid-based tool interface for mobile/tablet
* **Dark mode** toggle
* **Keyboard shortcuts**
* **Variable definitions** stored across calculations
* **Command history** persisted via localStorage

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
diff(sin(x^2), x)                    // Differentiate
integrate(x * e^x, x)                // Integrate
solve(x^2 - 4 = 0, x)                // Solve equation
[[1, 2], [3, 4]] * [x, y]            // Matrix-vector multiply
plot(sin(x), x, -10, 10)             // Plot function
limit(sin(x)/x, x, 0)                // Limit
laplace(t^2, t, s)                   // Laplace transform
fourierSeries(x^2, x, 0, 5)          // Fourier series
det([[1,2],[3,4]])                   // Determinant
```

---

## Development & Testing

```bash
npm install              # Install dev dependencies (jsdom for UI tests)
npm test                 # Run ALL tests (40 files)
npm run test:core        # Run only core tests
npm run test:latex       # Run LaTeX rendering tests
npm run test:chemistry   # Run chemistry tests
```

### Architecture

| File | Purpose |
|------|---------|
| `js/parser.js` | Lexer and Parser (Recursive Descent) → AST |
| `js/expression.js` | AST node definitions (`Expr`, `Add`, `Mul`, `Call`, etc.) |
| `js/cas.js` | Core CAS engine: evaluation, substitution, algorithms |
| `js/help.js` | Help system and command documentation |
| `js/chemistry.js` | Chemistry module (molar mass, balancing) |
| `index.html` | Frontend UI (command-line + app mode) |
| `tests/` | 40 Node.js test files (3169 lines) |

---

## Project Stats

| Metric | Value |
|--------|-------|
| Total lines of code | ~34,000 |
| JavaScript engine | 26,367 lines |
| Frontend | 7,466 lines |
| Test files | 40 |
| Test lines | 3,169 |
| Feature domains | 15 |
| Live since | March 2025 |

---

## License

MIT — see LICENSE file.

## Android App + APK Releases

This repository now also contains an Android wrapper built with Capacitor.

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

