window.HELP_DATA = {
    // --- Algebra & Arithmetic ---
    "simplify": {
        "description": "Simplifies an algebraic expression.",
        "syntax": "simplify(expression)",
        "examples": ["simplify((x^2 - 1)/(x - 1))", "simplify(sin(x)^2 + cos(x)^2)"]
    },
    "expand": {
        "description": "Expands a factored expression.",
        "syntax": "expand(expression)",
        "examples": ["expand((x+1)^3)", "expand((a+b)*(a-b))"]
    },
    "factor": {
        "description": "Factors a polynomial or integer.",
        "syntax": "factor(expression)",
        "examples": ["factor(x^2 - 4)", "factor(120)"]
    },
    "solve": {
        "description": "Solves an equation or system of equations.",
        "syntax": "solve(eq, var) or solve([eq1, eq2], [var1, var2])",
        "examples": ["solve(x^2 - 4, x)", "solve([x+y=10, x-y=2], [x, y])"]
    },
    "gcd": {
        "description": "Greatest Common Divisor of numbers or polynomials.",
        "syntax": "gcd(a, b, ...)",
        "examples": ["gcd(12, 18)", "gcd(x^2-1, x^2+2x+1)"]
    },
    "lcm": {
        "description": "Least Common Multiple.",
        "syntax": "lcm(a, b, ...)",
        "examples": ["lcm(4, 6)", "lcm(x, x^2)"]
    },
    "mod": {
        "description": "Modulo operator (remainder).",
        "syntax": "mod(a, b)",
        "examples": ["mod(10, 3)", "mod(x^3, x-1)"]
    },
    "modInverse": {
        "description": "Modular multiplicative inverse.",
        "syntax": "modInverse(a, m)",
        "examples": ["modInverse(3, 11)"]
    },
    "modPow": {
        "description": "Modular exponentiation (base^exp mod m).",
        "syntax": "modPow(base, exp, m)",
        "examples": ["modPow(2, 10, 1000)"]
    },
    "nCr": {
        "description": "Binomial coefficient (combinations).",
        "syntax": "nCr(n, k)",
        "examples": ["nCr(5, 2)"]
    },
    "nPr": {
        "description": "Permutations.",
        "syntax": "nPr(n, k)",
        "examples": ["nPr(5, 2)"]
    },
    "factorial": {
        "description": "Factorial of a number.",
        "syntax": "factorial(n) or n!",
        "examples": ["factorial(5)", "5!"]
    },
    "isPrime": {
        "description": "Checks if a number is prime (1 if true, 0 if false).",
        "syntax": "isPrime(n)",
        "examples": ["isPrime(17)", "isPrime(20)"]
    },
    "isSquare": {
        "description": "Checks if a number is a perfect square.",
        "syntax": "isSquare(n)",
        "examples": ["isSquare(16)", "isSquare(15)"]
    },
    "nextprime": {
        "description": "Returns the smallest prime greater than n.",
        "syntax": "nextprime(n)",
        "examples": ["nextprime(100)"]
    },
    "cfrac": {
        "description": "Continued fraction expansion.",
        "syntax": "cfrac(val, [depth])",
        "examples": ["cfrac(1.2)", "cfrac(pi, 5)"]
    },
    "circleEquation": {
        "description": "Finds circle equation from center and radius.",
        "syntax": "circleEquation(center, radius)",
        "examples": ["circleEquation([0, 0], 5)"]
    },
    "planeEquation": {
        "description": "Finds plane equation from 3 points.",
        "syntax": "planeEquation(p1, p2, p3)",
        "examples": ["planeEquation([1,0,0], [0,1,0], [0,0,1])"]
    },
    "completeSquare": {
        "description": "Completes the square for a quadratic polynomial.",
        "syntax": "completeSquare(expr, var)",
        "examples": ["completeSquare(x^2 + 6x + 10, x)"]
    },
    "abs": { "description": "Absolute value.", "syntax": "abs(x)", "examples": ["abs(-5)"] },
    "ceil": { "description": "Ceiling function.", "syntax": "ceil(x)", "examples": ["ceil(4.2)"] },
    "floor": { "description": "Floor function.", "syntax": "floor(x)", "examples": ["floor(4.8)"] },
    "round": { "description": "Rounds to nearest integer.", "syntax": "round(x)", "examples": ["round(4.5)"] },
    "sign": { "description": "Sign function (sgn).", "syntax": "sign(x)", "examples": ["sign(-5)"] },
    "sqrt": { "description": "Square root.", "syntax": "sqrt(x)", "examples": ["sqrt(16)"] },
    "cbrt": { "description": "Cube root.", "syntax": "cbrt(x)", "examples": ["cbrt(27)"] },
    "max": { "description": "Maximum value.", "syntax": "max(a, b, ...)", "examples": ["max(1, 5, 3)"] },
    "min": { "description": "Minimum value.", "syntax": "min(a, b, ...)", "examples": ["min(1, 5, 3)"] },
    "clamp": { "description": "Clamps value between min and max.", "syntax": "clamp(x, min, max)", "examples": ["clamp(5, 0, 10)"] },
    "exp": { "description": "Exponential function (e^x).", "syntax": "exp(x)", "examples": ["exp(2)"] },
    "ln": { "description": "Natural logarithm.", "syntax": "ln(x)", "examples": ["ln(e)"] },
    "log": { "description": "Logarithm (base 10 by default).", "syntax": "log(x, [base])", "examples": ["log(100)", "log(8, 2)"] },
    "delta": { "description": "Dirac delta function (alias).", "syntax": "delta(x)", "examples": ["delta(x)"] },
    "collect": { "description": "Collects terms with same powers.", "syntax": "collect(expr, var)", "examples": ["collect(x^2 + 2*x + x^2, x)"] },
    "factored": { "description": "Returns factored form of expression.", "syntax": "factored(expr)", "examples": ["factored(x^2 - 1)"] },
    "roots": { "description": "Finds roots of a polynomial.", "syntax": "roots(poly, var)", "examples": ["roots(x^2 - 4, x)"] },
    "groebner": { "description": "Computes Groebner basis.", "syntax": "groebner(polys, vars)", "examples": ["groebner([x^2-y, x*y-1], [x,y])"] },
    "sturm": { "description": "Sturm sequence of a polynomial.", "syntax": "sturm(poly, var)", "examples": ["sturm(x^2 - 2, x)"] },
    "laurent": { "description": "Laurent series expansion.", "syntax": "laurent(expr, var, point, order)", "examples": ["laurent(1/x, x, 0, 3)"] },
    "pade": { "description": "Pade approximant.", "syntax": "pade(expr, var, m, n)", "examples": ["pade(exp(x), x, 2, 2)"] },
    "partfrac": { "description": "Partial fraction decomposition.", "syntax": "partfrac(expr, var)", "examples": ["partfrac(1/(x^2-1), x)"] },

// --- Calculus ---
    "analyze": {
        "description": "Analyzes a function (roots, extrema, asymptotes, inflection points).",
        "syntax": "analyze(f, x)",
        "examples": ["analyze(x^3 - 3x, x)"]
    },
    "extrema": {
        "description": "Finds local extrema (minima/maxima) of a function.",
        "syntax": "extrema(f, x)",
        "examples": ["extrema(x^3 - 3x, x)"]
    },
    "stationary_points": {
        "description": "Finds stationary points (f'(x) = 0).",
        "syntax": "stationary_points(f, x)",
        "examples": ["stationary_points(x^3 - 3x, x)"]
    },
    "asymptotes": {
        "description": "Finds vertical and horizontal asymptotes.",
        "syntax": "asymptotes(f, x)",
        "examples": ["asymptotes((x^2-1)/(x-2), x)"]
    },
    "diff": {
        "description": "Calculates the derivative of an expression.",
        "syntax": "diff(expr, var, [order])",
        "examples": ["diff(sin(x), x)", "diff(x^3, x, 2)"]
    },
    "implicitDiff": {
        "description": "Calculates implicit derivative dy/dx.",
        "syntax": "implicitDiff(eq, y, x)",
        "examples": ["implicitDiff(x^2 + y^2 = 1, y, x)"]
    },
    "integrate": {
        "description": "Calculates the integral.",
        "syntax": "integrate(expr, var, [lower], [upper])",
        "examples": ["integrate(x^2, x)", "integrate(sin(x), x, 0, pi)"]
    },
    "limit": {
        "description": "Calculates the limit of an expression.",
        "syntax": "limit(expr, var, point)",
        "examples": ["limit(sin(x)/x, x, 0)"]
    },
    "taylor": {
        "description": "Computes the Taylor series expansion.",
        "syntax": "taylor(expr, var, point, order)",
        "examples": ["taylor(sin(x), x, 0, 5)"]
    },
    "sum": {
        "description": "Calculates the summation.",
        "syntax": "sum(expr, var, start, end)",
        "examples": ["sum(k^2, k, 1, 10)"]
    },
    "product": {
        "description": "Calculates the product.",
        "syntax": "product(expr, var, start, end)",
        "examples": ["product(k, k, 1, 5)"]
    },
    "desolve": {
        "description": "Solves ordinary differential equations.",
        "syntax": "desolve(eq, depVar)",
        "examples": ["desolve(diff(y, x) = y, y)", "desolve(diff(y, x, 2) + y = 0, y)"]
    },
    "laplace": {
        "description": "Laplace Transform.",
        "syntax": "laplace(expr, t, s)",
        "examples": ["laplace(sin(t), t, s)", "laplace(exp(a*t), t, s)"]
    },
    "ilaplace": {
        "description": "Inverse Laplace Transform.",
        "syntax": "ilaplace(expr, s, t)",
        "examples": ["ilaplace(1/(s^2+1), s, t)"]
    },
    "fourier": {
        "description": "Fourier Series approximation.",
        "syntax": "fourier(expr, var, n, [L])",
        "examples": ["fourier(x^2, x, 3, pi)"]
    },
    "tangent": {
        "description": "Equation of the tangent line.",
        "syntax": "tangent(expr, var, point)",
        "examples": ["tangent(x^2, x, 1)"]
    },
    "arcLen": {
        "description": "Arc Length of a curve.",
        "syntax": "arcLen(expr, var, start, end)",
        "examples": ["arcLen(x^2, x, 0, 1)"]
    },
    "curvature": {
        "description": "Curvature of a function.",
        "syntax": "curvature(expr, var, [point])",
        "examples": ["curvature(x^2, x)", "curvature(x^2, x, 0)"]
    },
    "grad": {
        "description": "Gradient of a scalar field.",
        "syntax": "grad(expr, [vars])",
        "examples": ["grad(x^2 + y^2, [x, y])"]
    },
    "curl": {
        "description": "Curl of a vector field.",
        "syntax": "curl([Fx, Fy, Fz], [x, y, z])",
        "examples": ["curl([-y, x, 0], [x, y, z])"]
    },
    "divergence": {
        "description": "Divergence of a vector field.",
        "syntax": "divergence([Fx, Fy, Fz], [x, y, z])",
        "examples": ["divergence([x, y, z], [x, y, z])"]
    },
    "potential": {
        "description": "Finds the scalar potential of a conservative vector field.",
        "syntax": "potential(field, vars)",
        "examples": ["potential([2x, 2y], [x, y])"]
    },
    "conservative": {
        "description": "Checks if a vector field is conservative.",
        "syntax": "conservative(field, vars)",
        "examples": ["conservative([-y, x], [x, y])"]
    },
    "jacobian": {
        "description": "Jacobian matrix of a vector function.",
        "syntax": "jacobian(vector, vars)",
        "examples": ["jacobian([x^2*y, x+y], [x, y])"]
    },
    "hessian": {
        "description": "Hessian matrix of a scalar function.",
        "syntax": "hessian(expr, vars)",
        "examples": ["hessian(x^3 + y^3 - 3xy, [x, y])"]
    },
    "arcLength": { "description": "Arc length of a curve.", "syntax": "arcLength(f, x, a, b)", "examples": ["arcLength(x^2, x, 0, 1)"] },
    "lineIntegral": { "description": "Line integral of a vector field.", "syntax": "lineIntegral(F, r, t, a, b)", "examples": ["lineIntegral([x,y], [cos(t),sin(t)], t, 0, 2*pi)"] },
    "line_integral": { "description": "Alias for lineIntegral.", "syntax": "line_integral(F, r, t, a, b)", "examples": ["line_integral([x,y], [cos(t),sin(t)], t, 0, 2*pi)"] },
    "eulerLagrange": { "description": "Euler-Lagrange equation.", "syntax": "eulerLagrange(L, x, t)", "examples": ["eulerLagrange(1/2*m*v^2 - m*g*x, x, t)"] },
    "euler_lagrange": { "description": "Alias for eulerLagrange.", "syntax": "euler_lagrange(L, x, t)", "examples": ["euler_lagrange(1/2*m*v^2 - m*g*x, x, t)"] },
    "laplacian": { "description": "Laplacian operator.", "syntax": "laplacian(f, vars)", "examples": ["laplacian(x^2 + y^2, [x,y])"] },
    "numeric_integrate": { "description": "Numeric integration.", "syntax": "numeric_integrate(f, x, a, b)", "examples": ["numeric_integrate(exp(-x^2), x, 0, 1)"] },
    "rk4": { "description": "Runge-Kutta 4th order ODE solver.", "syntax": "rk4(f, y0, t0, tn, h)", "examples": ["rk4(y - t^2 + 1, 0.5, 0, 2, 0.2)"] },
    "fourier_transform": { "description": "Fourier transform.", "syntax": "fourier_transform(f, t, w)", "examples": ["fourier_transform(exp(-t^2), t, w)"] },
    "inverse_fourier_transform": { "description": "Inverse Fourier transform.", "syntax": "inverse_fourier_transform(F, w, t)", "examples": ["inverse_fourier_transform(exp(-w^2/4), w, t)"] },
    "iztrans": { "description": "Inverse Z-transform.", "syntax": "iztrans(F, z, n)", "examples": ["iztrans(z/(z-1), z, n)"] },
    "ztrans": { "description": "Z-transform.", "syntax": "ztrans(f, n, z)", "examples": ["ztrans(1, n, z)"] },
    "ft": { "description": "Alias for fourier_transform.", "syntax": "ft(f, t, w)", "examples": ["ft(exp(-t^2), t, w)"] },
    "ift": { "description": "Alias for inverse_fourier_transform.", "syntax": "ift(F, w, t)", "examples": ["ift(exp(-w^2/4), w, t)"] },
    "fft": { "description": "Fast Fourier Transform.", "syntax": "fft(list)", "examples": ["fft([1,0,1,0])"] },
    "ifft": { "description": "Inverse Fast Fourier Transform.", "syntax": "ifft(list)", "examples": ["ifft([2,0,2,0])"] },

// --- Linear Algebra ---
    "det": {
        "description": "Determinant of a matrix.",
        "syntax": "det(matrix)",
        "examples": ["det([[1, 2], [3, 4]])"]
    },
    "inv": {
        "description": "Inverse of a matrix.",
        "syntax": "inv(matrix)",
        "examples": ["inv([[1, 2], [3, 4]])"]
    },
    "trans": {
        "description": "Transpose of a matrix.",
        "syntax": "trans(matrix)",
        "examples": ["trans([[1, 2], [3, 4]])"]
    },
    "trace": {
        "description": "Trace of a matrix (sum of diagonal).",
        "syntax": "trace(matrix)",
        "examples": ["trace([[1, 2], [3, 4]])"]
    },
    "rref": {
        "description": "Reduced Row Echelon Form.",
        "syntax": "rref(matrix)",
        "examples": ["rref([[1, 2, 3], [4, 5, 6]])"]
    },
    "rank": {
        "description": "Rank of a matrix.",
        "syntax": "rank(matrix)",
        "examples": ["rank([[1, 2], [2, 4]])"]
    },
    "eigenvals": {
        "description": "Eigenvalues of a matrix.",
        "syntax": "eigenvals(matrix)",
        "examples": ["eigenvals([[1, 2], [2, 1]])"]
    },
    "eigenvects": {
        "description": "Eigenvectors of a matrix.",
        "syntax": "eigenvects(matrix)",
        "examples": ["eigenvects([[1, 2], [2, 1]])"]
    },
    "diagonalize": {
        "description": "Diagonalizes a matrix A = PDP^-1. Returns [P, D].",
        "syntax": "diagonalize(matrix)",
        "examples": ["diagonalize([[1, 2], [2, 1]])"]
    },
    "lu": {
        "description": "LU Decomposition. Returns [L, U, P].",
        "syntax": "lu(matrix)",
        "examples": ["lu([[4, 3], [6, 3]])"]
    },
    "qr": {
        "description": "QR Decomposition. Returns [Q, R].",
        "syntax": "qr(matrix)",
        "examples": ["qr([[1, 2], [3, 4]])"]
    },
    "cholesky": {
        "description": "Cholesky Decomposition (A = LL^T).",
        "syntax": "cholesky(matrix)",
        "examples": ["cholesky([[4, 12, -16], [12, 37, -43], [-16, -43, 98]])"]
    },
    "svd": {
        "description": "Singular Value Decomposition. Returns [U, S, V].",
        "syntax": "svd(matrix)",
        "examples": ["svd([[1, 2], [3, 4]])"]
    },
    "cross": {
        "description": "Cross product of two 3D vectors.",
        "syntax": "cross(u, v)",
        "examples": ["cross([1, 0, 0], [0, 1, 0])"]
    },
    "dot": {
        "description": "Dot product of two vectors.",
        "syntax": "dot(u, v)",
        "examples": ["dot([1, 2], [3, 4])"]
    },
    "norm": {
        "description": "Euclidean norm (length) of a vector.",
        "syntax": "norm(v)",
        "examples": ["norm([3, 4])"]
    },
    "gramschmidt": {
        "description": "Orthonormalizes a set of vectors.",
        "syntax": "gramschmidt(vectors)",
        "examples": ["gramschmidt([[1, 1], [1, 0]])"]
    },
    "kron": {
        "description": "Kronecker Product of two matrices.",
        "syntax": "kron(A, B)",
        "examples": ["kron([[1, 2]], [[0, 1]])"]
    },
    "isDiagonal": { "description": "Checks if matrix is diagonal.", "syntax": "isDiagonal(matrix)", "examples": ["isDiagonal(identity(3))"] },
    "isSymmetric": { "description": "Checks if matrix is symmetric.", "syntax": "isSymmetric(matrix)", "examples": ["isSymmetric([[1, 2], [2, 1]])"] },
    "isOrthogonal": { "description": "Checks if matrix is orthogonal.", "syntax": "isOrthogonal(matrix)", "examples": ["isOrthogonal([[0, -1], [1, 0]])"] },
    "charpoly": { "description": "Characteristic polynomial of a matrix.", "syntax": "charpoly(matrix, var)", "examples": ["charpoly([[1,2],[3,4]], x)"] },
    "cond": { "description": "Condition number of a matrix.", "syntax": "cond(matrix)", "examples": ["cond([[1,2],[3,4]])"] },
    "hilbert": { "description": "Hilbert matrix of size n.", "syntax": "hilbert(n)", "examples": ["hilbert(3)"] },
    "toeplitz": { "description": "Toeplitz matrix.", "syntax": "toeplitz(col, row)", "examples": ["toeplitz([1,2,3], [1,4,5])"] },
    "vandermonde": { "description": "Vandermonde matrix.", "syntax": "vandermonde(list)", "examples": ["vandermonde([1,2,3])"] },
    "pinv": { "description": "Moore-Penrose pseudoinverse.", "syntax": "pinv(matrix)", "examples": ["pinv([[1,2],[3,6]])"] },
    "nullity": { "description": "Nullity of a matrix.", "syntax": "nullity(matrix)", "examples": ["nullity([[1,2],[2,4]])"] },
    "colSpace": { "description": "Column space of a matrix.", "syntax": "colSpace(matrix)", "examples": ["colSpace([[1,2],[3,4]])"] },
    "col_space": { "description": "Alias for colSpace.", "syntax": "col_space(matrix)", "examples": ["col_space([[1,2],[3,4]])"] },
    "rowSpace": { "description": "Row space of a matrix.", "syntax": "rowSpace(matrix)", "examples": ["rowSpace([[1,2],[3,4]])"] },
    "row_space": { "description": "Alias for rowSpace.", "syntax": "row_space(matrix)", "examples": ["row_space([[1,2],[3,4]])"] },
    "dim": { "description": "Dimension of a vector space.", "syntax": "dim(space)", "examples": ["dim(colSpace(A))"] },
    "matrixExp": { "description": "Matrix exponential.", "syntax": "matrixExp(matrix)", "examples": ["matrixExp([[0,1],[-1,0]])"] },
    "matrix_pow": { "description": "Matrix power.", "syntax": "matrix_pow(matrix, n)", "examples": ["matrix_pow([[1,1],[0,1]], 3)"] },
    "matrix_exp": { "description": "Alias for matrixExp.", "syntax": "matrix_exp(matrix)", "examples": ["matrix_exp([[0,1],[-1,0]])"] },
    "matrixPow": { "description": "Alias for matrix_pow.", "syntax": "matrixPow(matrix, n)", "examples": ["matrixPow([[1,1],[0,1]], 3)"] },
    "gram_schmidt": { "description": "Gram-Schmidt orthogonalization.", "syntax": "gram_schmidt(vectors)", "examples": ["gram_schmidt([[1,1],[1,-1]])"] },
    "eye": { "description": "Alias for identity matrix.", "syntax": "eye(n)", "examples": ["eye(3)"] },
    "idn": { "description": "Alias for identity matrix.", "syntax": "idn(n)", "examples": ["idn(3)"] },
    "tran": { "description": "Alias for transpose.", "syntax": "tran(matrix)", "examples": ["tran([[1,2],[3,4]])"] },
    "transpose": { "description": "Transpose of a matrix.", "syntax": "transpose(matrix)", "examples": ["transpose([[1,2],[3,4]])"] },

// --- Statistics ---
    "mean": {
        "description": "Arithmetic mean of a list.",
        "syntax": "mean(list)",
        "examples": ["mean([1, 2, 3, 4, 5])"]
    },
    "median": {
        "description": "Median of a list.",
        "syntax": "median(list)",
        "examples": ["median([1, 5, 2, 8, 3])"]
    },
    "mode": {
        "description": "Mode of a list.",
        "syntax": "mode(list)",
        "examples": ["mode([1, 2, 2, 3])"]
    },
    "variance": {
        "description": "Sample variance.",
        "syntax": "variance(list)",
        "examples": ["variance([1, 2, 3, 4, 5])"]
    },
    "std": {
        "description": "Sample standard deviation.",
        "syntax": "std(list)",
        "examples": ["std([1, 2, 3, 4, 5])"]
    },
    "cov": {
        "description": "Sample covariance between two lists.",
        "syntax": "cov(list1, list2)",
        "examples": ["cov([1, 2, 3], [2, 4, 6])"]
    },
    "corr": {
        "description": "Pearson correlation coefficient.",
        "syntax": "corr(list1, list2)",
        "examples": ["corr([1, 2, 3], [1, 2, 3])"]
    },
    "linearRegression": {
        "description": "Linear regression (y = mx + b).",
        "syntax": "linearRegression(points)",
        "examples": ["linearRegression([[1, 1], [2, 3], [3, 2]])"]
    },
    "polyRegression": {
        "description": "Polynomial regression.",
        "syntax": "polyRegression(points, degree)",
        "examples": ["polyRegression([[1, 1], [2, 4], [3, 9]], 2)"]
    },
    "normalPDF": { "description": "Normal distribution PDF.", "syntax": "normalPDF(x, mu, sigma)", "examples": ["normalPDF(0, 0, 1)"] },
    "normalCDF": { "description": "Normal distribution CDF.", "syntax": "normalCDF(x, mu, sigma)", "examples": ["normalCDF(1.96, 0, 1)"] },
    "invNorm": { "description": "Inverse Normal CDF (Quantile).", "syntax": "invNorm(area, mu, sigma)", "examples": ["invNorm(0.975, 0, 1)"] },
    "binomialPDF": { "description": "Binomial distribution probability.", "syntax": "binomialPDF(k, n, p)", "examples": ["binomialPDF(2, 5, 0.5)"] },
    "poissonPDF": { "description": "Poisson distribution probability.", "syntax": "poissonPDF(k, lambda)", "examples": ["poissonPDF(3, 5)"] },
    "chisquareCDF": { "description": "Chi-Square distribution CDF.", "syntax": "chisquareCDF(x, k)", "examples": ["chisquareCDF(5, 2)"] },
    "studentTCDF": { "description": "Student's t-distribution CDF.", "syntax": "studentTCDF(t, df)", "examples": ["studentTCDF(2.0, 10)"] },
    "zTest": { "description": "One-sample Z-test. Returns [z, p].", "syntax": "zTest(data, mu0, sigma)", "examples": ["zTest([1,2,3], 0, 1)"] },
    "tTest": { "description": "One-sample T-test. Returns [t, df].", "syntax": "tTest(data, mu0)", "examples": ["tTest([1,2,3], 0)"] },
    "zInterval": { "description": "Z Confidence Interval. Returns [min, max].", "syntax": "zInterval(data, sigma, level)", "examples": ["zInterval([1,2,3], 1, 0.95)"] },
    "tInterval": { "description": "T Confidence Interval. Returns [min, max].", "syntax": "tInterval(data, level)", "examples": ["tInterval([1,2,3], 0.95)"] },
    "propTest": { "description": "One Proportion Z-Test. Returns [z, p].", "syntax": "propTest(successes, n, p0)", "examples": ["propTest(45, 100, 0.5)"] },
    "propTest2": { "description": "Two Proportion Z-Test. Returns [z, p].", "syntax": "propTest2(x1, n1, x2, n2)", "examples": ["propTest2(40, 100, 50, 100)"] },
    "expectedValue": { "description": "Expected value of a distribution.", "syntax": "expectedValue(values, probs)", "examples": ["expectedValue([1,2], [0.5,0.5])"] },
    "entropy": { "description": "Shannon entropy.", "syntax": "entropy(probs)", "examples": ["entropy([0.5, 0.5])"] },
    "kl_divergence": { "description": "Kullback-Leibler divergence.", "syntax": "kl_divergence(p, q)", "examples": ["kl_divergence([0.5,0.5], [0.4,0.6])"] },
    "kldiv": { "description": "Alias for kl_divergence.", "syntax": "kldiv(p, q)", "examples": ["kldiv([0.5,0.5], [0.4,0.6])"] },
    "mae": { "description": "Mean Absolute Error.", "syntax": "mae(actual, predicted)", "examples": ["mae([1,2], [1.1,1.9])"] },
    "mse": { "description": "Mean Squared Error.", "syntax": "mse(actual, predicted)", "examples": ["mse([1,2], [1.1,1.9])"] },
    "mst": { "description": "Minimum Spanning Tree.", "syntax": "mst(graph)", "examples": ["mst([[0,1],[1,0]])"] },
    "betaCDF": { "description": "Beta Cumulative Distribution Function.", "syntax": "betaCDF(x, alpha, beta)", "examples": ["betaCDF(0.5, 2, 2)"] },
    "betaPDF": { "description": "Beta Probability Density Function.", "syntax": "betaPDF(x, alpha, beta)", "examples": ["betaPDF(0.5, 2, 2)"] },
    "gammaCDF": { "description": "Gamma Cumulative Distribution Function.", "syntax": "gammaCDF(x, k, theta)", "examples": ["gammaCDF(1, 2, 1)"] },
    "gammaPDF": { "description": "Gamma Probability Density Function.", "syntax": "gammaPDF(x, k, theta)", "examples": ["gammaPDF(1, 2, 1)"] },
    "fCDF": { "description": "F-distribution Cumulative Distribution Function.", "syntax": "fCDF(x, d1, d2)", "examples": ["fCDF(1, 5, 5)"] },
    "fPDF": { "description": "F-distribution Probability Density Function.", "syntax": "fPDF(x, d1, d2)", "examples": ["fPDF(1, 5, 5)"] },
    "hypergeometricCDF": { "description": "Hypergeometric Cumulative Distribution Function.", "syntax": "hypergeometricCDF(k, N, K, n)", "examples": ["hypergeometricCDF(2, 50, 5, 10)"] },
    "hypergeometricPDF": { "description": "Hypergeometric Probability Density Function.", "syntax": "hypergeometricPDF(k, N, K, n)", "examples": ["hypergeometricPDF(2, 50, 5, 10)"] },
    "uniformCDF": { "description": "Uniform Cumulative Distribution Function.", "syntax": "uniformCDF(x, a, b)", "examples": ["uniformCDF(0.5, 0, 1)"] },
    "uniformPDF": { "description": "Uniform Probability Density Function.", "syntax": "uniformPDF(x, a, b)", "examples": ["uniformPDF(0.5, 0, 1)"] },

// --- Plotting ---
    "plot": {
        "description": "Plots a function.",
        "syntax": "plot(expr, var, [min], [max])",
        "examples": ["plot(sin(x), x)", "plot(x^2, x, -5, 5)"]
    },
    "plot3d": {
        "description": "3D Surface Plot.",
        "syntax": "plot3d(expr, x, y)",
        "examples": ["plot3d(x^2 + y^2, x, y)"]
    },
    "plotparam": {
        "description": "Parametric Plot.",
        "syntax": "plotparam([x(t), y(t)], t, min, max)",
        "examples": ["plotparam([sin(t), cos(t)], t, 0, 2*pi)"]
    },
    "plotpolar": {
        "description": "Polar Plot.",
        "syntax": "plotpolar(r(t), t, min, max)",
        "examples": ["plotpolar(1 + cos(t), t, 0, 2*pi)"]
    },
    "plotimplicit": {
        "description": "Implicit Plot f(x,y) = 0.",
        "syntax": "plotimplicit(eq, x, y)",
        "examples": ["plotimplicit(x^2 + y^2 = 1, x, y)"]
    },
    "slopefield": {
        "description": "Slope field for ODE y' = f(x,y).",
        "syntax": "slopefield(f(x,y), x, y)",
        "examples": ["slopefield(x+y, x, y)"]
    },
    "plotlist": { "description": "Plots a list of points.", "syntax": "plotlist(list)", "examples": ["plotlist([[1,2], [2,3]])"] },

// --- Finance ---
    "compound": {
        "description": "Compound Interest.",
        "syntax": "compound(P, r, n, t)",
        "examples": ["compound(1000, 0.05, 12, 10)"]
    },
    "loan": {
        "description": "Monthly Loan Payment.",
        "syntax": "loan(P, r_annual, years)",
        "examples": ["loan(200000, 0.04, 30)"]
    },
    "npv": {
        "description": "Net Present Value.",
        "syntax": "npv(rate, [flows])",
        "examples": ["npv(0.05, [-1000, 500, 600])"]
    },
    "irr": {
        "description": "Internal Rate of Return.",
        "syntax": "irr([flows])",
        "examples": ["irr([-1000, 500, 600])"]
    },
    "blackScholes": { "description": "Black-Scholes option pricing.", "syntax": "blackScholes(S, K, T, r, sigma, type)", "examples": ["blackScholes(100, 100, 1, 0.05, 0.2, 'call')"] },
    "black_scholes": { "description": "Alias for blackScholes.", "syntax": "black_scholes(S, K, T, r, sigma, type)", "examples": ["black_scholes(100, 100, 1, 0.05, 0.2, 'call')"] },
    "amortization": { "description": "Amortization schedule.", "syntax": "amortization(principal, rate, periods)", "examples": ["amortization(10000, 0.05, 12)"] },
    "annuity": { "description": "Annuity calculation.", "syntax": "annuity(payment, rate, periods)", "examples": ["annuity(100, 0.05, 10)"] },

// --- Complex Numbers ---
    "abs": { "description": "Absolute value or Magnitude.", "syntax": "abs(z)", "examples": ["abs(3+4i)"] },
    "arg": { "description": "Argument (phase angle).", "syntax": "arg(z)", "examples": ["arg(1+i)"] },
    "real": { "description": "Real part.", "syntax": "real(z)", "examples": ["real(3+4i)"] },
    "imag": { "description": "Imaginary part.", "syntax": "imag(z)", "examples": ["imag(3+4i)"] },
    "conj": { "description": "Complex conjugate.", "syntax": "conj(z)", "examples": ["conj(3+4i)"] },
    "toPolar": { "description": "Converts to polar form [r, deg].", "syntax": "toPolar(z)", "examples": ["toPolar(1+i)"] },

    // --- Chemistry ---
    "molarMass": {
        "description": "Calculate molar mass of formula.",
        "syntax": "molarMass(formula)",
        "examples": ["molarMass(H2O)", "molarMass(C6H12O6)"]
    },
    "balance": {
        "description": "Balance chemical equation.",
        "syntax": "balance(equation)",
        "examples": ["balance(H2 + O2 -> H2O)"]
    },

    // --- Logic ---
    "true": { "description": "Boolean true (1).", "syntax": "true", "examples": ["true"] },
    "false": { "description": "Boolean false (0).", "syntax": "false", "examples": ["false"] },
    "and": { "description": "Logical AND.", "syntax": "a and b", "examples": ["x > 0 and x < 5"] },
    "or": { "description": "Logical OR.", "syntax": "a or b", "examples": ["x < 0 or x > 5"] },
    "not": { "description": "Logical NOT.", "syntax": "not a", "examples": ["not (x=0)"] },
    "truthTable": { "description": "Generates a truth table.", "syntax": "truthTable(expr, [vars])", "examples": ["truthTable(A and B, [A, B])"] },
    "nand": { "description": "Logical NAND.", "syntax": "nand(a, b)", "examples": ["nand(true, false)"] },
    "nor": { "description": "Logical NOR.", "syntax": "nor(a, b)", "examples": ["nor(false, false)"] },
    "xnor": { "description": "Logical XNOR.", "syntax": "xnor(a, b)", "examples": ["xnor(true, true)"] },
    "cnf": { "description": "Conjunctive Normal Form.", "syntax": "cnf(expr)", "examples": ["cnf(a | (b & c))"] },
    "dnf": { "description": "Disjunctive Normal Form.", "syntax": "dnf(expr)", "examples": ["dnf(a & (b | c))"] },

// --- Polynomials & Number Theory ---
    "degree": { "description": "Degree of a polynomial.", "syntax": "degree(poly, var)", "examples": ["degree(x^3+x, x)"] },
    "coeff": { "description": "Coefficient of x^n.", "syntax": "coeff(poly, var, n)", "examples": ["coeff(3x^2+2x, x, 2)"] },
    "symb2poly": { "description": "Converts polynomial to coefficient list.", "syntax": "symb2poly(expr, var)", "examples": ["symb2poly(x^2+2x+1, x)"] },
    "poly2symb": { "description": "Converts coefficient list to polynomial.", "syntax": "poly2symb(list, var)", "examples": ["poly2symb([1, 2, 1], x)"] },
    "resultant": { "description": "Resultant of two polynomials.", "syntax": "resultant(p1, p2, var)", "examples": ["resultant(x^2-1, x+1, x)"] },
    "discriminant": { "description": "Discriminant of a polynomial.", "syntax": "discriminant(poly, var)", "examples": ["discriminant(x^2+bx+c, x)"] },
    "divisors": { "description": "List of integer divisors.", "syntax": "divisors(n)", "examples": ["divisors(12)"] },
    "euler": { "description": "Euler's totient function (phi).", "syntax": "euler(n)", "examples": ["euler(10)"] },
    "moebius": { "description": "Moebius function (mu). 0 if square factor, (-1)^k otherwise.", "syntax": "moebius(n)", "examples": ["moebius(10)", "moebius(12)"] },
    "sigma": { "description": "Divisor function. Sum of k-th powers of divisors.", "syntax": "sigma(n, [k])", "examples": ["sigma(12)", "sigma(12, 0)"] },
    "legendreSymbol": { "description": "Legendre Symbol (a/p).", "syntax": "legendreSymbol(a, p)", "examples": ["legendreSymbol(2, 7)"] },
    "isPerfect": { "description": "Checks if a number is a perfect number.", "syntax": "isPerfect(n)", "examples": ["isPerfect(6)", "isPerfect(28)"] },
    "fibonacci": { "description": "nth Fibonacci number.", "syntax": "fibonacci(n)", "examples": ["fibonacci(10)"] },
    "gamma": { "description": "Gamma function.", "syntax": "gamma(z)", "examples": ["gamma(0.5)"] },
    "beta": { "description": "Beta function B(x, y).", "syntax": "beta(x, y)", "examples": ["beta(2, 3)"] },
    "stirling1": { "description": "Signed Stirling numbers of the first kind.", "syntax": "stirling1(n, k)", "examples": ["stirling1(4, 2)"] },
    "stirling2": { "description": "Stirling numbers of the second kind.", "syntax": "stirling2(n, k)", "examples": ["stirling2(4, 2)"] },
    "bell": { "description": "Bell number (sum of Stirling2).", "syntax": "bell(n)", "examples": ["bell(4)"] },
    "zeta": { "description": "Riemann Zeta function.", "syntax": "zeta(s)", "examples": ["zeta(2)", "zeta(4)"] },
    "lambertw": { "description": "Lambert W function (product log).", "syntax": "lambertw(x)", "examples": ["lambertw(1)", "lambertw(e)"] },
    "divisorSum": { "description": "Sum of divisors.", "syntax": "divisorSum(n, k)", "examples": ["divisorSum(10, 1)"] },
    "ifactor": { "description": "Integer factorization.", "syntax": "ifactor(n)", "examples": ["ifactor(60)"] },
    "primeFactors": { "description": "List of prime factors.", "syntax": "primeFactors(n)", "examples": ["primeFactors(60)"] },
    "primitiveRoot": { "description": "Primitive root modulo n.", "syntax": "primitiveRoot(n)", "examples": ["primitiveRoot(7)"] },
    "isPrimitiveRoot": { "description": "Checks if g is primitive root mod n.", "syntax": "isPrimitiveRoot(g, n)", "examples": ["isPrimitiveRoot(3, 7)"] },
    "prev_prime": { "description": "Previous prime before n.", "syntax": "prev_prime(n)", "examples": ["prev_prime(10)"] },
    "prevprime": { "description": "Alias for prev_prime.", "syntax": "prevprime(n)", "examples": ["prevprime(10)"] },
    "next_prime": { "description": "Next prime after n.", "syntax": "next_prime(n)", "examples": ["next_prime(10)"] },
    "jacobiSymbol": { "description": "Jacobi symbol.", "syntax": "jacobiSymbol(a, n)", "examples": ["jacobiSymbol(2, 5)"] },
    "chineseRemainder": { "description": "Chinese Remainder Theorem.", "syntax": "chineseRemainder(rem, mod)", "examples": ["chineseRemainder([2,3], [3,5])"] },
    "crt": { "description": "Alias for chineseRemainder.", "syntax": "crt(rem, mod)", "examples": ["crt([2,3], [3,5])"] },
    "knapsack": { "description": "Solves 0/1 knapsack problem.", "syntax": "knapsack(weights, values, capacity)", "examples": ["knapsack([1,2,3], [10,20,30], 4)"] },
    "kroneckerDelta": { "description": "Kronecker delta function.", "syntax": "kroneckerDelta(i, j)", "examples": ["kroneckerDelta(1, 1)"] },
    "stirling1": { "description": "Stirling numbers of first kind.", "syntax": "stirling1(n, k)", "examples": ["stirling1(4, 2)"] },
    "stirling2": { "description": "Stirling numbers of second kind.", "syntax": "stirling2(n, k)", "examples": ["stirling2(4, 2)"] },
    "bell": { "description": "Bell number.", "syntax": "bell(n)", "examples": ["bell(5)"] },
    "comb": { "description": "Binomial coefficient (combinations).", "syntax": "comb(n, k)", "examples": ["comb(5, 2)"] },
    "combinations": { "description": "Generates all combinations.", "syntax": "combinations(list, k)", "examples": ["combinations([1,2,3], 2)"] },
    "perm": { "description": "Number of permutations.", "syntax": "perm(n, k)", "examples": ["perm(5, 2)"] },
    "permutations": { "description": "Generates all permutations.", "syntax": "permutations(list)", "examples": ["permutations([1,2,3])"] },

// --- Calculus Extras ---
    "minimize": { "description": "Finds local minimum.", "syntax": "minimize(f, x)", "examples": ["minimize(x^2-4x, x)"] },
    "maximize": { "description": "Finds local maximum.", "syntax": "maximize(f, x)", "examples": ["maximize(-x^2, x)"] },
    "nIntegrate": { "description": "Numeric Integration (Simpson's Rule).", "syntax": "nIntegrate(f, x, a, b)", "examples": ["nIntegrate(exp(-x^2), x, 0, 1)"] },
    "vectorfield": { "description": "Plots a 2D vector field.", "syntax": "vectorfield([u, v], x, y)", "examples": ["vectorfield([-y, x], x, y)"] },
    "trigReduce": { "description": "Linearize trigonometric powers.", "syntax": "trigReduce(expr)", "examples": ["trigReduce(sin(x)^2)"] },
    "trigExpand": { "description": "Expand trigonometric sums.", "syntax": "trigExpand(expr)", "examples": ["trigExpand(sin(2x))"] },

    // --- List & Matrix Tools ---
    "seq": { "description": "Generate a sequence.", "syntax": "seq(expr, var, start, end, step)", "examples": ["seq(k^2, k, 1, 5, 1)"] },
    "range": { "description": "Generate a range of numbers.", "syntax": "range(start, end, step)", "examples": ["range(0, 10, 2)"] },
    "union": { "description": "Union of two lists.", "syntax": "union(list1, list2)", "examples": ["union([1,2], [2,3])"] },
    "intersect": { "description": "Intersection of two lists.", "syntax": "intersect(list1, list2)", "examples": ["intersect([1,2], [2,3])"] },
    "setdiff": { "description": "Set difference of two lists.", "syntax": "setdiff(list1, list2)", "examples": ["setdiff([1,2], [2,3])"] },
    "sort": { "description": "Sorts a list.", "syntax": "sort(list)", "examples": ["sort([3, 1, 2])"] },
    "reverse": { "description": "Reverses a list.", "syntax": "reverse(list)", "examples": ["reverse([1, 2, 3])"] },
    "size": { "description": "Size of a list or vector.", "syntax": "size(list)", "examples": ["size([1, 2, 3])"] },
    "flatten": { "description": "Flattens nested lists.", "syntax": "flatten(list)", "examples": ["flatten([[1, 2], [3]])"] },
    "cumsum": { "description": "Cumulative sum of a list.", "syntax": "cumsum(list)", "examples": ["cumsum([1, 2, 3])"] },
    "diag": { "description": "Create diagonal matrix from list.", "syntax": "diag(list)", "examples": ["diag([1, 2, 3])"] },
    "identity": { "description": "Identity matrix of size n.", "syntax": "identity(n)", "examples": ["identity(3)"] },
    "zeros": { "description": "Matrix of zeros.", "syntax": "zeros(rows, cols)", "examples": ["zeros(2, 3)"] },
    "ones": { "description": "Matrix of ones.", "syntax": "ones(rows, cols)", "examples": ["ones(2, 3)"] },
    "kernel": { "description": "Kernel (Nullspace) of a matrix.", "syntax": "kernel(matrix)", "examples": ["kernel([[1, 2], [2, 4]])"] },
    "basis": { "description": "Basis of the column space.", "syntax": "basis(matrix)", "examples": ["basis([[1, 2], [3, 4]])"] },
    "append": { "description": "Appends element to list.", "syntax": "append(list, element)", "examples": ["append([1,2], 3)"] },
    "prepend": { "description": "Prepends element to list.", "syntax": "prepend(list, element)", "examples": ["prepend([2,3], 1)"] },
    "concat": { "description": "Concatenates lists.", "syntax": "concat(list1, list2)", "examples": ["concat([1,2], [3,4])"] },
    "set": { "description": "Removes duplicates from list.", "syntax": "set(list)", "examples": ["set([1,2,2,3])"] },
    "unitVector": { "description": "Unit vector.", "syntax": "unitVector(v)", "examples": ["unitVector([1,1])"] },
    "projection": { "description": "Vector projection.", "syntax": "projection(u, v)", "examples": ["projection([1,2], [3,4])"] },
    "distance": { "description": "Distance between points.", "syntax": "distance(p1, p2)", "examples": ["distance([0,0], [3,4])"] },
    "midpoint": { "description": "Midpoint between points.", "syntax": "midpoint(p1, p2)", "examples": ["midpoint([0,0], [2,2])"] },

// --- Statistics Extras ---
    "geoMean": { "description": "Geometric Mean.", "syntax": "geoMean(list)", "examples": ["geoMean([1, 2, 4])"] },
    "harmMean": { "description": "Harmonic Mean.", "syntax": "harmMean(list)", "examples": ["harmMean([1, 2, 4])"] },
    "rms": { "description": "Root Mean Square.", "syntax": "rms(list)", "examples": ["rms([1, -1])"] },
    "mad": { "description": "Mean Absolute Deviation.", "syntax": "mad(list)", "examples": ["mad([1, 2, 3])"] },
    "moment": { "description": "k-th Central Moment.", "syntax": "moment(list, k)", "examples": ["moment([1, 2, 3], 2)"] },
    "skewness": { "description": "Skewness of a dataset.", "syntax": "skewness(list)", "examples": ["skewness([1, 2, 3, 10])"] },
    "kurtosis": { "description": "Kurtosis of a dataset.", "syntax": "kurtosis(list)", "examples": ["kurtosis([1, 2, 3])"] },
    "binomialCDF": { "description": "Binomial Cumulative Probability.", "syntax": "binomialCDF(k, n, p)", "examples": ["binomialCDF(2, 5, 0.5)"] },
    "poissonCDF": { "description": "Poisson Cumulative Probability.", "syntax": "poissonCDF(k, lambda)", "examples": ["poissonCDF(3, 5)"] },
    "exponentialPDF": { "description": "Exponential PDF.", "syntax": "exponentialPDF(x, lambda)", "examples": ["exponentialPDF(1, 0.5)"] },
    "exponentialCDF": { "description": "Exponential CDF.", "syntax": "exponentialCDF(x, lambda)", "examples": ["exponentialCDF(1, 0.5)"] },
    "geometricPDF": { "description": "Geometric PDF.", "syntax": "geometricPDF(k, p)", "examples": ["geometricPDF(3, 0.5)"] },
    "geometricCDF": { "description": "Geometric CDF.", "syntax": "geometricCDF(k, p)", "examples": ["geometricCDF(3, 0.5)"] },
    "chisquarePDF": { "description": "Chi-Square PDF.", "syntax": "chisquarePDF(x, k)", "examples": ["chisquarePDF(2, 3)"] },
    "invChiSquare": { "description": "Inverse Chi-Square CDF.", "syntax": "invChiSquare(area, k)", "examples": ["invChiSquare(0.95, 3)"] },
    "studentTPDF": { "description": "Student's t PDF.", "syntax": "studentTPDF(x, df)", "examples": ["studentTPDF(0, 5)"] },
    "invT": { "description": "Inverse Student's t CDF.", "syntax": "invT(area, df)", "examples": ["invT(0.95, 10)"] },
    "poly_regression": { "description": "Polynomial regression.", "syntax": "poly_regression(points, degree)", "examples": ["poly_regression([[1,1], [2,4], [3,9]], 2)"] },
    "exp_regression": { "description": "Exponential regression.", "syntax": "exp_regression(points)", "examples": ["exp_regression([[1,2], [2,4], [3,8]])"] },
    "log_regression": { "description": "Logarithmic regression.", "syntax": "log_regression(points)", "examples": ["log_regression([[1,0], [10,1], [100,2]])"] },
    "power_regression": { "description": "Power regression.", "syntax": "power_regression(points)", "examples": ["power_regression([[1,1], [2,8], [3,27]])"] },
    "linear_regression": { "description": "Linear regression.", "syntax": "linear_regression(points)", "examples": ["linear_regression([[1,2], [2,4], [3,6]])"] },
    "anova": { "description": "One-way ANOVA.", "syntax": "anova(groups)", "examples": ["anova([[1,2,3], [4,5,6]])"] },
    "chiSquareTest": { "description": "Chi-Square Test.", "syntax": "chiSquareTest(observed, expected)", "examples": ["chiSquareTest([10,20], [15,15])"] },

// --- Special Functions ---
    "besselJ": { "description": "Bessel function of the first kind.", "syntax": "besselJ(n, x)", "examples": ["besselJ(0, x)"] },
    "besselY": { "description": "Bessel function of the second kind.", "syntax": "besselY(n, x)", "examples": ["besselY(0, x)"] },
    "legendre": { "description": "Legendre polynomial P_n(x).", "syntax": "legendre(n, x)", "examples": ["legendre(3, x)"] },
    "chebyshev": { "description": "Chebyshev polynomial T_n(x).", "syntax": "chebyshev(n, x)", "examples": ["chebyshev(3, x)"] },
    "hermite": { "description": "Hermite polynomial H_n(x).", "syntax": "hermite(n, x)", "examples": ["hermite(3, x)"] },
    "laguerre": { "description": "Laguerre polynomial L_n(x).", "syntax": "laguerre(n, x)", "examples": ["laguerre(3, x)"] },
    "surfaceArea": {
        "description": "Surface area of revolution.",
        "syntax": "surfaceArea(f, x, a, b, [axis])",
        "examples": ["surfaceArea(x^2, x, 0, 1)", "surfaceArea(x^2, x, 0, 1, y)"]
    },
    "isSubset": { "description": "Checks if list1 is a subset of list2.", "syntax": "isSubset(list1, list2)", "examples": ["isSubset([1,2], [1,2,3])"] },
    "cartesianProduct": { "description": "Cartesian product of two lists.", "syntax": "cartesianProduct(l1, l2)", "examples": ["cartesianProduct([1,2], [a,b])"] },
    "EllipticE": { "description": "Incomplete elliptic integral of the second kind.", "syntax": "EllipticE(phi, m)", "examples": ["EllipticE(x, 0.5)"] },
    "EllipticK": { "description": "Complete elliptic integral of the first kind.", "syntax": "EllipticK(m)", "examples": ["EllipticK(0.5)"] },
    "FresnelC": { "description": "Fresnel cosine integral.", "syntax": "FresnelC(x)", "examples": ["FresnelC(x)"] },
    "FresnelS": { "description": "Fresnel sine integral.", "syntax": "FresnelS(x)", "examples": ["FresnelS(x)"] },
    "Heaviside": { "description": "Heaviside step function.", "syntax": "Heaviside(x)", "examples": ["Heaviside(x-1)"] },
    "LambertW": { "description": "Alias for lambertw.", "syntax": "LambertW(x)", "examples": ["LambertW(x)"] },
    "dirac": { "description": "Dirac delta function.", "syntax": "dirac(x)", "examples": ["dirac(x)"] },
    "sinc": { "description": "Sinc function (sin(x)/x).", "syntax": "sinc(x)", "examples": ["sinc(x)"] },
    "phi": { "description": "Alias for euler (totient function).", "syntax": "phi(n)", "examples": ["phi(10)"] },
    "bernoulli": { "description": "Bernoulli number.", "syntax": "bernoulli(n)", "examples": ["bernoulli(4)"] },
    "bernoulliPoly": { "description": "Bernoulli polynomial.", "syntax": "bernoulliPoly(n, x)", "examples": ["bernoulliPoly(3, x)"] },
    "eulerPoly": { "description": "Euler polynomial.", "syntax": "eulerPoly(n, x)", "examples": ["eulerPoly(3, x)"] },
    "harmonic": { "description": "Harmonic number.", "syntax": "harmonic(n)", "examples": ["harmonic(5)"] },
    "catalan": { "description": "Catalan number.", "syntax": "catalan(n)", "examples": ["catalan(4)"] },

// --- Engineering ---
    "cis": { "description": "Polar complex form: cos(x) + i*sin(x).", "syntax": "cis(angle_deg)", "examples": ["cis(90)"] },
    "phasor": { "description": "Phasor form.", "syntax": "phasor(mag, angle_deg)", "examples": ["phasor(10, 45)"] },
    "parallel": { "description": "Parallel impedance (1 / sum(1/Z)).", "syntax": "parallel(z1, z2, ...)", "examples": ["parallel(10, 10)", "par(10, 20)"] },
    "ctrb": { "description": "Controllability matrix.", "syntax": "ctrb(A, B)", "examples": ["ctrb(A, B)"] },
    "obsv": { "description": "Observability matrix.", "syntax": "obsv(A, C)", "examples": ["obsv(A, C)"] },
    "ss2tf": { "description": "State-space to transfer function.", "syntax": "ss2tf(A, B, C, D)", "examples": ["ss2tf(A, B, C, D)"] },
    "acker": { "description": "Alias for ackermann.", "syntax": "acker(A, B, p)", "examples": ["acker(A, B, [-1, -2])"] },
    "ackermann": { "description": "Ackermann's formula for pole placement.", "syntax": "ackermann(A, B, p)", "examples": ["ackermann(A, B, [-1, -2])"] },
    "cross_product": { "description": "Cross product of 3D vectors.", "syntax": "cross_product(u, v)", "examples": ["cross_product([1,0,0], [0,1,0])"] },
    "dot_product": { "description": "Dot product of vectors.", "syntax": "dot_product(u, v)", "examples": ["dot_product([1,2], [3,4])"] },
    "del": { "description": "Del (nabla) operator.", "syntax": "del(f, vars)", "examples": ["del(x^2*y, [x,y])"] },

// --- Utils ---
    "root": { "description": "N-th root of x.", "syntax": "root(x, n)", "examples": ["root(8, 3)", "root(x, 2)"] },
    "approx": { "description": "Numeric approximation.", "syntax": "approx(expr)", "examples": ["approx(pi)"] },
    "clear": { "description": "Clears variables and history.", "syntax": "clear()", "examples": ["clear()"] },
    "help": { "description": "Displays help.", "syntax": "help([command])", "examples": ["help()", "help(plot)"] },

    "cos": { "description": "Cosine function.", "syntax": "cos(x)", "examples": ["cos(pi/2)"] },
    "cot": { "description": "Cotangent function.", "syntax": "cot(x)", "examples": ["cot(pi/4)"] },
    "csc": { "description": "Cosecant function.", "syntax": "csc(x)", "examples": ["csc(pi/6)"] },
    "sec": { "description": "Secant function.", "syntax": "sec(x)", "examples": ["sec(pi/3)"] },
    "tan": { "description": "Tangent function.", "syntax": "tan(x)", "examples": ["tan(pi/4)"] },
    "tanh": { "description": "Hyperbolic tangent function.", "syntax": "tanh(x)", "examples": ["tanh(1)"] },
    "codegen": { "description": "Generates code for an expression.", "syntax": "codegen(expr, language)", "examples": ["codegen(x^2, 'javascript')"] },
    "toCode": { "description": "Alias for codegen.", "syntax": "toCode(expr, language)", "examples": ["toCode(x^2, 'python')"] },
    "rewrite": { "description": "Rewrites expression in terms of target function.", "syntax": "rewrite(expr, target)", "examples": ["rewrite(tan(x), sin)"] },
    "convert": { "description": "Converts units or forms.", "syntax": "convert(val, from, to)", "examples": ["convert(1, 'inch', 'cm')"] },
    "N": { "description": "Numeric evaluation.", "syntax": "N(expr)", "examples": ["N(pi)"] },
    "Q": { "description": "Rational number evaluation.", "syntax": "Q(expr)", "examples": ["Q(0.5)"] },
    "evalf": { "description": "Alias for approx.", "syntax": "evalf(expr)", "examples": ["evalf(pi)"] },
    "purge": { "description": "Alias for clear.", "syntax": "purge()", "examples": ["purge()"] },
    "latex": { "description": "Returns LaTeX representation.", "syntax": "latex(expr)", "examples": ["latex(x^2/2)"] },
    "tex": { "description": "Alias for latex.", "syntax": "tex(expr)", "examples": ["tex(x^2)"] },
    "subs": { "description": "Alias for substitute.", "syntax": "subs(expr, old, new)", "examples": ["subs(x^2, x, 2)"] },
    "substitute": { "description": "Substitutes value for variable.", "syntax": "substitute(expr, old, new)", "examples": ["substitute(x+y, x, 3)"] },
    "var": { "description": "Declares variables.", "syntax": "var(name)", "examples": ["var('x')"] },
    "length": { "description": "Alias for size.", "syntax": "length(list)", "examples": ["length([1,2,3])"] },
    "fsolve": { "description": "Numeric root finding.", "syntax": "fsolve(expr, var, guess)", "examples": ["fsolve(cos(x)-x, x, 0.5)"] },
    "piecewise": { "description": "Piecewise function.", "syntax": "piecewise(cond1, val1, cond2, val2, ...)", "examples": ["piecewise(x<0, -1, x>=0, 1)"] },
    "conv": { "description": "Alias for convolution.", "syntax": "conv(list1, list2)", "examples": ["conv([1,1], [1,1])"] },
    "convolution": { "description": "Convolution of lists.", "syntax": "convolution(list1, list2)", "examples": ["convolution([1,2], [3,4])"] },
    "xcorr": { "description": "Cross-correlation of lists.", "syntax": "xcorr(list1, list2)", "examples": ["xcorr([1,2], [2,1])"] },
    "lagrange": { "description": "Lagrange polynomial interpolation.", "syntax": "lagrange(points, var)", "examples": ["lagrange([[1,2],[2,3]], x)"] },
    "lagrangeMultipliers": { "description": "Method of Lagrange multipliers.", "syntax": "lagrangeMultipliers(f, g, vars)", "examples": ["lagrangeMultipliers(x^2+y^2, x+y-1, [x,y])"] },
    "wronskian": { "description": "Wronskian of functions.", "syntax": "wronskian(funcs, var)", "examples": ["wronskian([sin(x), cos(x)], x)"] },
    "shortestPath": { "description": "Shortest path in a graph.", "syntax": "shortestPath(graph, start, end)", "examples": ["shortestPath(adj_matrix, 0, 2)"] },
    "shortest_path": { "description": "Alias for shortestPath.", "syntax": "shortest_path(graph, start, end)", "examples": ["shortest_path(adj_matrix, 0, 2)"] },
    "simplex": { "description": "Simplex method for linear programming.", "syntax": "simplex(obj, constraints)", "examples": ["simplex([1,1], [[1,2,4],[2,1,5]])"] },
    "numRealRoots": { "description": "Number of real roots.", "syntax": "numRealRoots(poly, var)", "examples": ["numRealRoots(x^2-4, x)"] },
    "num_real_roots": { "description": "Alias for numRealRoots.", "syntax": "num_real_roots(poly, var)", "examples": ["num_real_roots(x^2+1, x)"] },
    "area_surface": { "description": "Alias for surfaceArea.", "syntax": "area_surface(f, x, a, b)", "examples": ["area_surface(x^2, x, 0, 1)"] },
    "volumeSolid": { "description": "Volume of solid of revolution.", "syntax": "volumeSolid(f, x, a, b, axis)", "examples": ["volumeSolid(x^2, x, 0, 1, 'x')"] },
    "volume_solid": { "description": "Alias for volumeSolid.", "syntax": "volume_solid(f, x, a, b, axis)", "examples": ["volume_solid(x^2, x, 0, 1, 'y')"] },
    "toCylindrical": { "description": "Convert coordinates to cylindrical.", "syntax": "toCylindrical(x, y, z)", "examples": ["toCylindrical(1, 1, 1)"] },
    "toSpherical": { "description": "Convert coordinates to spherical.", "syntax": "toSpherical(x, y, z)", "examples": ["toSpherical(1, 1, 1)"] },
    "rect": { "description": "Convert polar to rectangular.", "syntax": "rect(r, theta)", "examples": ["rect(1, pi/4)"] },
    "par": { "description": "Alias for parallel.", "syntax": "par(z1, z2)", "examples": ["par(10, 10)"] },
    "map": { "description": "Apply function to list.", "syntax": "map(f, list)", "examples": ["map(x -> x^2, [1,2,3])"] },
    "series": { "description": "Series expansion.", "syntax": "series(expr, var, point, order)", "examples": ["series(sin(x), x, 0, 3)"] },
    "quo": { "description": "Polynomial quotient.", "syntax": "quo(num, den, var)", "examples": ["quo(x^2-1, x-1, x)"] },
    "rem": { "description": "Polynomial remainder.", "syntax": "rem(num, den, var)", "examples": ["rem(x^2, x-1, x)"] },
    "div": { "description": "Integer division.", "syntax": "div(a, b)", "examples": ["div(10, 3)"] },
    "irem": { "description": "Integer remainder.", "syntax": "irem(a, b)", "examples": ["irem(10, 3)"] },
    "erfcinv": { "description": "Inverse complementary error function.", "syntax": "erfcinv(x)", "examples": ["erfcinv(0.5)"] },
    "erfinv": { "description": "Inverse error function.", "syntax": "erfinv(x)", "examples": ["erfinv(0.5)"] },
    "interpolate": { "description": "Data interpolation.", "syntax": "interpolate(points, x)", "examples": ["interpolate([[1,2],[2,3]], 1.5)"] },
    "angle": { "description": "Angle of complex number.", "syntax": "angle(z)", "examples": ["angle(1+i)"] },
    "atomicWeight": { "description": "Atomic weight of an element.", "syntax": "atomicWeight(symbol)", "examples": ["atomicWeight('H')"] },
    "lineEquation": { "description": "Equation of a line through two points.", "syntax": "lineEquation(p1, p2, vars)", "examples": ["lineEquation([0,0], [1,1], [x,y])"] },
    "rsolve": { "description": "Solve recurrence relation.", "syntax": "rsolve(eq, f(n))", "examples": ["rsolve(f(n)=f(n-1)+f(n-2), f(n))"] },
    "Sigma": { "description": "Summation.", "syntax": "Sigma(expr, var, a, b)", "examples": ["Sigma(k^2, k, 1, 5)"] },
};
    "toUnit": {
        "description": "Converts a value from one unit to another.",
        "syntax": "toUnit(value, fromUnit, toUnit)",
        "examples": ["toUnit(1, km, m)", "toUnit(0, C, F)", "toUnit(1000, g, kg)"]
    },
    "toBase": {
        "description": "Converts an integer to a string in the given base (2-36).",
        "syntax": "toBase(number, base)",
        "examples": ["toBase(255, 16)", "toBase(10, 2)"]
    },
    "fromBase": {
        "description": "Converts a string in the given base (2-36) to a decimal integer.",
        "syntax": "fromBase(str, base)",
        "examples": ["fromBase(FF, 16)", "fromBase(1010, 2)"]
    },
    "latex": {
        "description": "Converts an expression to LaTeX format.",
        "syntax": "latex(expr)",
        "examples": ["latex(x^2 + y)", "latex(integrate(x^2, x))"]
    },
    "toRoman": {
        "description": "Converts a positive integer to a Roman numeral string.",
        "syntax": "toRoman(n)",
        "examples": ["toRoman(2026)", "toRoman(100)"]
    },
    "fromRoman": {
        "description": "Converts a Roman numeral string to a decimal integer.",
        "syntax": "fromRoman(str)",
        "examples": ["fromRoman(MMXXVI)", "fromRoman(C)"]
    }
ENDHELP
echo "Help entries added: $(grep -c 'toUnit\|toBase\|fromBase\|latex\|toRoman\|fromRoman' js/help.js)"
