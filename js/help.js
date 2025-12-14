window.HELP_DATA = {
    "plot": {
        "description": "Plots a function of one variable.",
        "syntax": "plot(expression, variable, [min], [max])",
        "examples": [
            "plot(sin(x), x)",
            "plot(x^2, x, -5, 5)",
            "plot(tan(x), x, -10, 10)"
        ]
    },
    "diff": {
        "description": "Calculates the derivative of an expression.",
        "syntax": "diff(expression, variable, [order])",
        "examples": [
            "diff(sin(x), x)",
            "diff(x^3, x, 2)",
            "diff(exp(2*x), x)"
        ]
    },
    "integrate": {
        "description": "Calculates the indefinite or definite integral.",
        "syntax": "integrate(expression, variable, [lower], [upper])",
        "examples": [
            "integrate(x^2, x)",
            "integrate(sin(x), x, 0, pi)",
            "integrate(1/(1+x^2), x)"
        ]
    },
    "solve": {
        "description": "Solves an equation or system of equations.",
        "syntax": "solve(equation, variable) or solve([eq1, eq2], [var1, var2])",
        "examples": [
            "solve(x^2 - 4, x)",
            "solve(sin(x) = 0.5, x)",
            "solve([x+y=10, x-y=2], [x, y])"
        ]
    },
    "limit": {
        "description": "Calculates the limit of an expression at a point.",
        "syntax": "limit(expression, variable, point)",
        "examples": [
            "limit(sin(x)/x, x, 0)",
            "limit((1+1/x)^x, x, infinity)"
        ]
    },
    "sum": {
        "description": "Calculates the summation of an expression.",
        "syntax": "sum(expression, variable, start, end)",
        "examples": [
            "sum(k^2, k, 1, 10)",
            "sum(1/k^2, k, 1, infinity)"
        ]
    },
    "simplify": {
        "description": "Simplifies an algebraic expression.",
        "syntax": "simplify(expression)",
        "examples": [
            "simplify((x^2 - 1)/(x - 1))",
            "simplify(sin(x)^2 + cos(x)^2)"
        ]
    },
    "expand": {
        "description": "Expands a factored expression.",
        "syntax": "expand(expression)",
        "examples": [
            "expand((x+1)^3)",
            "expand((a+b)*(a-b))"
        ]
    },
    "desolve": {
        "description": "Solves ordinary differential equations.",
        "syntax": "desolve(equation, dependent_var)",
        "examples": [
            "desolve(diff(y, x) = y, y)",
            "desolve(diff(y, x, 2) + y = 0, y)"
        ]
    },
    "matrix": {
        "description": "Matrix operations (det, inv, trace, eigenvals).",
        "syntax": "det(M), inv(M), etc.",
        "examples": [
            "det([[1, 2], [3, 4]])",
            "inv([[1, 2], [3, 4]])",
            "eigenvals([[1, 2], [2, 1]])"
        ]
    },
    "mean": {
        "description": "Calculates the arithmetic mean of a list.",
        "syntax": "mean(list)",
        "examples": [
            "mean([1, 2, 3, 4, 5])",
            "mean([10, 20, 30])"
        ]
    },
    "studentT": {
        "description": "Student's t-distribution functions.",
        "syntax": "studentTPDF(x, df), studentTCDF(x, df), invT(p, df)",
        "examples": [
            "studentTPDF(0, 10)",
            "studentTCDF(1.5, 10)",
            "invT(0.95, 10)"
        ]
    },
    "chisquare": {
        "description": "Chi-Square distribution functions.",
        "syntax": "chisquarePDF(x, k), chisquareCDF(x, k), invChiSquare(p, k)",
        "examples": [
            "chisquarePDF(2, 4)",
            "chisquareCDF(5, 4)",
            "invChiSquare(0.95, 4)"
        ]
    }
};
