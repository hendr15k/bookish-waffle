# Web CAS

A client-side Computer Algebra System (CAS) built with JavaScript.

## Features

- **Symbolic Math**: Simplify expressions, derivatives, integrals, limits, taylor series.
- **Solving Equations**: Linear and quadratic equation solver.
- **Linear Algebra**: Matrix determinant, transpose.
- **Number Theory**: GCD, LCM, Factorization.
- **Statistics**: Mean, Variance.
- **Plotting**: 2D function plotting.
- **Output**: LaTeX rendering via MathJax.

## Usage

Open `index.html` in a web browser.

### Commands

- `diff(expr, var)`: Differentiate expression.
- `integrate(expr, var, [min, max])`: Integrate expression.
- `simplify(expr)`: Simplify expression.
- `solve(eq, var)`: Solve equation.
- `plot(expr, var, min, max)`: Plot function.
- `gcd(a, b)`, `lcm(a, b)`: GCD and LCM of integers.
- `factor(n)`: Prime factorization of integer.
- `mean([list])`, `variance([list])`: Statistics.

## Development

The core logic is in `js/`.
- `expression.js`: AST node definitions.
- `parser.js`: Lexer and Parser.
- `cas.js`: Evaluation and CAS algorithms.
