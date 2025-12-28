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

    // --- Engineering ---
    "cis": { "description": "Polar complex form: cos(x) + i*sin(x).", "syntax": "cis(angle_deg)", "examples": ["cis(90)"] },
    "phasor": { "description": "Phasor form.", "syntax": "phasor(mag, angle_deg)", "examples": ["phasor(10, 45)"] },
    "parallel": { "description": "Parallel impedance (1 / sum(1/Z)).", "syntax": "parallel(z1, z2, ...)", "examples": ["parallel(10, 10)", "par(10, 20)"] },

    // --- Utils ---
    "root": { "description": "N-th root of x.", "syntax": "root(x, n)", "examples": ["root(8, 3)", "root(x, 2)"] },
    "approx": { "description": "Numeric approximation.", "syntax": "approx(expr)", "examples": ["approx(pi)"] },
    "clear": { "description": "Clears variables and history.", "syntax": "clear()", "examples": ["clear()"] },
    "help": { "description": "Displays help.", "syntax": "help([command])", "examples": ["help()", "help(plot)"] }
};
