# Web CAS (Client-Side)

A lightweight, client-side Computer Algebra System (CAS) built with vanilla JavaScript. It runs entirely in your browser without any backend dependencies.

[Live Demo](https://hendr15k.github.io/bookish-waffle/)

## Features

This CAS supports a wide range of mathematical operations, including:

### Arithmetic & Algebra
*   **Basic Operations**: `+`, `-`, `*`, `/`, `^`
*   **Symbolic Computation**: Variables, simplification, and expansion.
*   **Functions**: `expand(expr)`, `simplify(expr)`, `factor(n)`, `gcd(a, b)`, `lcm(a, b)`, `factorial(n)`
*   **Equation Solving**: `solve(equation, variable)`

### Calculus
*   **Differentiation**: `diff(expr, var)`
*   **Integration**: `integrate(expr, var)` (indefinite) and `integrate(expr, var, lower, upper)` (definite)
*   **Limits**: `limit(expr, var, point)`
*   **Series**: `taylor(expr, var, point, order)`
*   **Sums & Products**: `sum(expr, var, start, end)`, `product(expr, var, start, end)`

### Linear Algebra
*   **Vectors & Matrices**: Define using `[a, b]` or `[[a, b], [c, d]]`
*   **Operations**: Matrix multiplication, dot product, cross product (`cross(v1, v2)`)
*   **Properties**: Determinant (`det(M)`), Inverse (`inv(M)`), Transpose (`trans(M)`)

### Statistics
*   **Descriptive**: `mean(list)`, `variance(list)`

### Plotting
*   **2D Plotting**: `plot(expr, var, [min], [max])` renders function graphs on a canvas.

### Output
*   Results are rendered in **LaTeX** using MathJax for beautiful mathematical display.

## Usage

1.  Clone the repository or download the files.
2.  Open `index.html` in any modern web browser.
3.  Enter commands in the input field.

**Examples:**
*   `diff(sin(x^2), x)`
*   `integrate(x * e^x, x)`
*   `solve(x^2 - 4 = 0, x)`
*   `[[1, 2], [3, 4]] * [x, y]`
*   `plot(sin(x), x, -10, 10)`

## Development & Testing

The project is structured with `js/` containing the source code and `tests/` containing Node.js-based tests.

### Running Tests

To run the tests, you need [Node.js](https://nodejs.org/) installed.

```bash
# Run core functionality tests
node tests/test_core.js

# Run LaTeX rendering tests
node tests/test_latex.js
```

## Architecture

*   `js/parser.js`: Lexer and Parser (Recursive Descent) generating an Abstract Syntax Tree (AST).
*   `js/expression.js`: AST node definitions (`Expr`, `Add`, `Mul`, `Call`, etc.) and basic manipulation logic.
*   `js/cas.js`: The core CAS engine handling evaluation, substitution, and high-level algorithms (calculus, linear algebra).
*   `index.html`: The frontend UI.
