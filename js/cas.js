
class CAS {
    constructor() {
        this.variables = {
            'pi': new Num(Math.PI), // Store numeric value? Or should we use Symbol?
            // If we store Num, implicit N() works, but symbolic calc loses precision.
            // Let's store Symbol by default in a separate map or just rely on 'pi' symbol.
            // Wait, if I do `a := pi`, I want `a` to be the symbol pi.
            // If I do `N(pi)`, I want numeric.
            // So `variables` should store what the user assigned.
            // `pi` and `e` are constants.
        };
        // Pre-defined constants as Symbols in variables?
        // No, let them be resolved by Symbol.evaluateNumeric() or if user assigns them.

        this.functions = {};
    }

    evaluate(exprTree) {
        let evaluated = this._recursiveEval(exprTree);
        if (evaluated && typeof evaluated.simplify === 'function') {
            evaluated = evaluated.simplify();
        }

        // Store result in 'ans' variable if it's a value (Expression) and not a Plot/Action
        if (evaluated instanceof Expr) {
            this.variables['ans'] = evaluated;
        }

        return evaluated;
    }

    getVariables() {
        return this.variables;
    }

    getFunctions() {
        return this.functions;
    }

    _recursiveEval(node) {
        if (node instanceof Assignment) {
            const value = this._recursiveEval(node.value).simplify();
            this.variables[node.target.name] = value;
            return value;
        }

        if (node instanceof FunctionDef) {
            this.functions[node.name] = node;
            return node;
        }

        if (node instanceof Sym) {
            if (this.variables.hasOwnProperty(node.name)) {
                return this.variables[node.name];
            }
            return node;
        }

        if (node instanceof Call) {
            const args = node.args.map(arg => this._recursiveEval(arg));

            // Check for user-defined function
            if (this.functions.hasOwnProperty(node.funcName)) {
                const funcDef = this.functions[node.funcName];
                if (funcDef.params.length !== args.length) {
                    throw new Error(`Function ${node.funcName} expects ${funcDef.params.length} arguments, got ${args.length}`);
                }

                // Substitute arguments securely (to avoid collision if args contain params)
                // We use temporary placeholders
                let body = funcDef.body;
                const tempMap = {};

                // 1. Replace params with temporary unique symbols
                for (let i = 0; i < funcDef.params.length; i++) {
                    const paramName = funcDef.params[i];
                    const tempName = `__temp_param_${i}_${Date.now()}_${Math.random()}`; // Unique enough
                    tempMap[tempName] = args[i];
                    body = body.substitute(new Sym(paramName), new Sym(tempName));
                }

                // 2. Replace temporary symbols with actual arguments
                for (const tempName in tempMap) {
                    body = body.substitute(new Sym(tempName), tempMap[tempName]);
                }

                return this._recursiveEval(body); // Re-evaluate body
            }

            if (node.funcName === 'N') {
                if (args.length !== 1) throw new Error("N requires 1 argument");
                const val = args[0].evaluateNumeric();
                return new Num(val);
            }

            if (node.funcName === 'diff') {
                if (args.length < 2) throw new Error("diff requires at least 2 arguments");
                const func = args[0];
                const varNode = args[1];
                if (!(varNode instanceof Sym)) throw new Error("Second argument to diff must be a variable");
                return func.diff(varNode);
            }

            if (node.funcName === 'integrate') {
                if (args.length < 2) throw new Error("integrate requires at least 2 arguments");
                const func = args[0];
                const varNode = args[1];
                if (!(varNode instanceof Sym)) throw new Error("Second argument to integrate must be a variable");

                if (args.length === 4) {
                    const lower = args[2];
                    const upper = args[3];
                    const indefinite = func.integrate(varNode);

                    // If integration failed (returned a Call to integrate), return symbolic definite integral
                    if (indefinite instanceof Call && indefinite.funcName === 'integrate') {
                        return new Call('integrate', args);
                    }

                    const valUpper = indefinite.substitute(varNode, upper);
                    const valLower = indefinite.substitute(varNode, lower);
                    return new Sub(valUpper, valLower);
                }
                return func.integrate(varNode);
            }

            if (node.funcName === 'sum') {
                // sum(expr, var, start, end)
                if (args.length !== 4) throw new Error("sum requires 4 arguments: expression, variable, start, end");
                return this._sum(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'product') {
                // product(expr, var, start, end)
                if (args.length !== 4) throw new Error("product requires 4 arguments: expression, variable, start, end");
                return this._product(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'expand') {
                if (args.length !== 1) throw new Error("expand takes exactly 1 argument");
                return args[0].expand();
            }

            if (node.funcName === 'simplify') {
                if (args.length !== 1) throw new Error("simplify takes exactly 1 argument");
                return args[0].simplify();
            }

            if (node.funcName === 'solve') {
                 if (args.length < 2) throw new Error("solve requires at least 2 arguments: equation and variable");
                 const eq = args[0];
                 const varNode = args[1];
                 if (!(varNode instanceof Sym)) throw new Error("Second argument to solve must be a variable");
                 return this._solve(eq, varNode);
            }

            if (node.funcName === 'plot') {
                 if (args.length < 2) throw new Error("plot requires at least 2 arguments: expression and variable");
                 const expr = args[0];
                 const varNode = args[1];
                 const min = args.length > 2 ? args[2].value : -10;
                 const max = args.length > 3 ? args[3].value : 10;

                 // Return a special object that the frontend can recognize
                 return {
                     type: 'plot',
                     expr: expr,
                     var: varNode,
                     min: min,
                     max: max,
                     toString: () => `Plotting ${expr} from ${min} to ${max}`,
                     toLatex: () => `\\text{Plotting } ${expr.toLatex()}`
                 };
            }

            if (node.funcName === 'taylor') {
                // taylor(expr, var, point, order)
                if (args.length < 4) throw new Error("taylor requires 4 arguments: expression, variable, point, order");
                const expr = args[0];
                const varNode = args[1];
                const point = args[2];
                const order = args[3];

                if (!(varNode instanceof Sym)) throw new Error("Second argument to taylor must be a variable");
                if (!(order instanceof Num)) throw new Error("Order must be a number");

                return this._taylor(expr, varNode, point, order.value);
            }

            if (node.funcName === 'limit') {
                // limit(expr, var, point)
                if (args.length < 3) throw new Error("limit requires 3 arguments: expression, variable, point");
                const expr = args[0];
                const varNode = args[1];
                const point = args[2];
                if (!(varNode instanceof Sym)) throw new Error("Second argument to limit must be a variable");

                return this._limit(expr, varNode, point);
            }

            if (node.funcName === 'det') {
                if (args.length !== 1) throw new Error("det requires 1 argument");
                return this._det(args[0]);
            }

            if (node.funcName === 'inv') {
                if (args.length !== 1) throw new Error("inv requires 1 argument");
                return this._inv(args[0]);
            }

            if (node.funcName === 'cross') {
                if (args.length !== 2) throw new Error("cross requires 2 arguments");
                return this._cross(args[0], args[1]);
            }

            if (node.funcName === 'trans') {
                if (args.length !== 1) throw new Error("trans requires 1 argument");
                return this._trans(args[0]);
            }

            if (node.funcName === 'clear') {
                return { type: 'action', name: 'clear', toString: () => "Cleared", toLatex: () => "\\text{Cleared}" };
            }

            if (node.funcName === 'gcd') {
                if (args.length !== 2) throw new Error("gcd requires 2 arguments");
                return this._gcd(args[0], args[1]);
            }

            if (node.funcName === 'lcm') {
                if (args.length !== 2) throw new Error("lcm requires 2 arguments");
                return this._lcm(args[0], args[1]);
            }

            if (node.funcName === 'factor') {
                if (args.length !== 1) throw new Error("factor requires 1 argument");
                return this._factor(args[0]);
            }

            if (node.funcName === 'factorial') {
                 if (args.length !== 1) throw new Error("factorial requires 1 argument");
                 const n = args[0];
                 if (n instanceof Num && Number.isInteger(n.value) && n.value >= 0) {
                     return new Num(this._factorial(n.value));
                 }
                 return new Call('factorial', args);
            }

            if (node.funcName === 'fibonacci') {
                 if (args.length !== 1) throw new Error("fibonacci requires 1 argument");
                 const n = args[0];
                 if (n instanceof Num && Number.isInteger(n.value) && n.value >= 0) {
                     return new Num(this._fibonacci(n.value));
                 }
                 return new Call('fibonacci', args);
            }

            if (node.funcName === 'gamma') {
                 if (args.length !== 1) throw new Error("gamma requires 1 argument");
                 const z = args[0];
                 if (z instanceof Num) {
                     return new Num(this._gamma(z.value));
                 }
                 return new Call('gamma', args);
            }

            if (node.funcName === 'nCr') {
                if (args.length !== 2) throw new Error("nCr requires 2 arguments");
                return this._nCr(args[0], args[1]);
            }

            if (node.funcName === 'nPr') {
                if (args.length !== 2) throw new Error("nPr requires 2 arguments");
                return this._nPr(args[0], args[1]);
            }

            if (node.funcName === 'trace') {
                if (args.length !== 1) throw new Error("trace requires 1 argument");
                return this._trace(args[0]);
            }

            if (node.funcName === 'mean') {
                if (args.length !== 1) throw new Error("mean requires 1 argument (list)");
                return this._mean(args[0]);
            }

            if (node.funcName === 'variance') {
                if (args.length !== 1) throw new Error("variance requires 1 argument (list)");
                return this._variance(args[0]);
            }

            if (node.funcName === 'linearRegression') {
                if (args.length !== 1) throw new Error("linearRegression requires 1 argument (list of points)");
                return this._linearRegression(args[0]);
            }

            if (node.funcName === 'normalPDF') {
                if (args.length !== 3) throw new Error("normalPDF requires 3 arguments: x, mu, sigma");
                return this._normalPDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'binomialPDF') {
                if (args.length !== 3) throw new Error("binomialPDF requires 3 arguments: k, n, p");
                return this._binomialPDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'help') {
                const helpText = `Available commands:
diff(expr, var), integrate(expr, var, [lower, upper]),
limit(expr, var, val), taylor(expr, var, pt, order),
sum(expr, var, start, end), product(expr, var, start, end),
expand(expr), simplify(expr), solve(eq, var),
det(M), trans(M), plot(expr, var, [min, max]),
gcd(a, b), lcm(a, b), factor(n), factorial(n),
mean(list), variance(list),
N(expr) [numeric eval], clear(), help()`;

                const latexHelp = `\\begin{array}{l}
\\text{Available commands:} \\\\
\\text{diff}(expr, var), \\; \\int(expr, var, [a, b]), \\\\
\\lim(expr, var, val), \\; \\text{taylor}(expr, var, pt, n), \\\\
\\sum(expr, var, a, b), \\; \\prod(expr, var, a, b), \\\\
\\text{expand}(expr), \\; \\text{simplify}(expr), \\; \\text{solve}(eq, var), \\\\
\\text{det}(M), \\; \\text{trans}(M), \\; \\text{plot}(expr, var, [min, max]), \\\\
\\text{gcd}(a, b), \\; \\text{lcm}(a, b), \\; \\text{factor}(n), \\; n!, \\\\
\\text{mean}(L), \\; \\text{variance}(L), \\; \\text{linearRegression}(L), \\\\
\\text{normalPDF}(x, \\mu, \\sigma), \\; \\text{binomialPDF}(k, n, p), \\\\
N(expr), \\; \\text{clear}(), \\; \\text{help}()
\\end{array}`;

                return { type: 'info', text: helpText, toString: () => helpText, toLatex: () => latexHelp };
            }

            return new Call(node.funcName, args);
        }

        if (node instanceof BinaryOp) {
            const left = this._recursiveEval(node.left);
            const right = this._recursiveEval(node.right);
            return new node.constructor(left, right);
        }

        if (node instanceof Eq) {
            const left = this._recursiveEval(node.left);
            const right = this._recursiveEval(node.right);
            return new Eq(left, right);
        }

        if (node instanceof Vec) {
            return new Vec(node.elements.map(e => this._recursiveEval(e)));
        }

        if (node instanceof Block) {
             const results = [];
             for(const stmt of node.statements) {
                 results.push(this._recursiveEval(stmt));
             }

             // Check if we have multiple plots
             const plots = results.filter(r => r && r.type === 'plot');
             if (plots.length > 0 && plots.length === results.length) {
                  // Merge plots
                  const combined = {
                      type: 'plot',
                      isMulti: true,
                      plots: plots,
                      toString: () => plots.map(p => p.toString()).join("; "),
                      toLatex: () => plots.map(p => p.toLatex()).join("; ")
                  };
                  return combined;
             }

             // Return the last result if not all plots, or return a list?
             // Standard behavior for blocks is return last value.
             // But if we have mixed output (e.g. assignment then plot), we probably want the plot if it is last.
             // If we have plot then assignment, we want assignment.
             // However, for REPL, maybe we want to see all outputs?
             // But current frontend only handles one result.
             // If I have `1+1; 2+2`, I probably expect `4` or `2, 4`.
             // If I have `a:=1; b:=2`, I expect `2`.
             // If I have `plot(...); plot(...)`, user wants both.

             // Let's return the last result, UNLESS they are all plots (handled above).
             // Wait, what if I have `a:=1; plot(x,x)`? I want the plot.
             // What if I have `plot(x,x); a:=1`? I probably want `1`.

             // The user specifically asked for "In that case the program should plot both functions".
             // This implies if there are multiple plots, we should show them.

             // Let's check if there are ANY plots.
             if (plots.length > 0) {
                 if (plots.length === 1) return plots[0];
                 // Multiple plots found mixed with other things?
                 // `x=1; plot(x); plot(x^2)`
                 // We should probably return the combined plots.
                 return {
                      type: 'plot',
                      isMulti: true,
                      plots: plots,
                      toString: () => plots.map(p => p.toString()).join("; "),
                      toLatex: () => plots.map(p => p.toLatex()).join("; ")
                  };
             }

             return results[results.length - 1];
        }

        return node;
    }

    _sum(expr, varNode, start, end) {
        if (!(varNode instanceof Sym)) throw new Error("Sum variable must be a symbol");
        if (!(start instanceof Num) || !(end instanceof Num)) {
            // Symbolic sum not implemented, return call
            return new Call("sum", [expr, varNode, start, end]);
        }

        let sum = new Num(0);
        for (let i = start.value; i <= end.value; i++) {
            const term = expr.substitute(varNode, new Num(i)).simplify();
            sum = new Add(sum, term).simplify();
        }
        return sum;
    }

    _product(expr, varNode, start, end) {
        if (!(varNode instanceof Sym)) throw new Error("Product variable must be a symbol");
        if (!(start instanceof Num) || !(end instanceof Num)) {
            return new Call("product", [expr, varNode, start, end]);
        }

        let prod = new Num(1);
        for (let i = start.value; i <= end.value; i++) {
            const term = expr.substitute(varNode, new Num(i)).simplify();
            prod = new Mul(prod, term).simplify();
        }
        return prod;
    }

    // ... (rest of the methods: _solve, _taylor, _factorial, _limit, _det, _trans)
    _solve(eq, varNode) {
        // ... (same as before)
        let expr;
        if (eq instanceof Eq) {
            expr = new Sub(eq.left, eq.right).simplify();
        } else {
            expr = eq.simplify();
        }

        try {
            const b = expr.substitute(varNode, new Num(0)).simplify();
            const f1 = expr.substitute(varNode, new Num(1)).simplify();
            const a = new Sub(f1, b).simplify();

            const reconstruction = new Add(new Mul(a, varNode), b);
            const diff = new Sub(expr, reconstruction).expand().simplify();

            if (diff instanceof Num && diff.value === 0 && !(a instanceof Num && a.value === 0)) {
                return new Div(new Mul(new Num(-1), b), a).simplify();
            }

            const c = b;
            const fm1 = expr.substitute(varNode, new Num(-1)).simplify();

            const f1_c = new Sub(f1, c).simplify();
            const fm1_c = new Sub(fm1, c).simplify();
            const twoA = new Add(f1_c, fm1_c).simplify();
            const A = new Div(twoA, new Num(2)).simplify();

            const B = new Sub(f1_c, A).simplify();

            const quadRecon = new Add(new Add(new Mul(A, new Pow(varNode, new Num(2))), new Mul(B, varNode)), c);
            const diffQuad = new Sub(expr, quadRecon).expand().simplify();

            // console.log("Solve debug:", expr.toString(), "Recon:", quadRecon.toString(), "Diff:", diffQuad.toString());

            if (diffQuad instanceof Num && Math.abs(diffQuad.value) < 1e-10 && !(A instanceof Num && A.value === 0)) {
                const discriminant = new Sub(new Pow(B, new Num(2)), new Mul(new Num(4), new Mul(A, c))).simplify();
                const sqrtDisc = new Call("sqrt", [discriminant]); // Use sqrt instead of Pow(0.5)

                const sol1 = new Div(new Add(new Mul(new Num(-1), B), sqrtDisc), new Mul(new Num(2), A));
                const sol2 = new Div(new Sub(new Mul(new Num(-1), B), sqrtDisc), new Mul(new Num(2), A));

                return new Call("set", [sol1.simplify(), sol2.simplify()]);
            }

        } catch (e) {
            console.error(e);
        }

        return new Call("solve", [eq, varNode]);
    }

    _taylor(expr, varNode, point, order) {
        // ... (same as before)
        let result = new Num(0);
        let deriv = expr;

        let term = deriv.substitute(varNode, point).simplify();
        result = term;

        for (let n = 1; n <= order; n++) {
            deriv = deriv.diff(varNode).simplify();
            let coeff = deriv.substitute(varNode, point).simplify();

            if (coeff instanceof Num && coeff.value === 0) continue;

            const factorial = this._factorial(n);
            const divisor = new Num(factorial);
            const termCoeff = new Div(coeff, divisor).simplify();

            const x_minus_a = new Sub(varNode, point).simplify();
            const power = new Pow(x_minus_a, new Num(n)).simplify();

            const termN = new Mul(termCoeff, power);
            result = new Add(result, termN);
        }

        return result.simplify();
    }

    _factorial(n) {
        if (n === 0 || n === 1) return 1;
        let res = 1;
        for (let i = 2; i <= n; i++) res *= i;
        return res;
    }

    _fibonacci(n) {
        if (n === 0) return 0;
        if (n === 1) return 1;
        let a = 0, b = 1;
        for (let i = 2; i <= n; i++) {
            let temp = a + b;
            a = b;
            b = temp;
        }
        return b;
    }

    _gamma(z) {
        // Lanczos approximation
        const g = 7;
        const C = [
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        ];

        if (z < 0.5) return Math.PI / (Math.sin(Math.PI * z) * this._gamma(1 - z));

        z -= 1;
        let x = C[0];
        for (let i = 1; i < g + 2; i++) x += C[i] / (z + i);

        const t = z + g + 0.5;
        return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
    }

    _limit(expr, varNode, point, depth = 0) {
        if (depth > 5) return new Call("limit", [expr, varNode, point]);

        if (expr instanceof Div) {
            let num = expr.left.substitute(varNode, point).simplify();
            let den = expr.right.substitute(varNode, point).simplify();

            // Handle Num or zero-value Num from simplification
            const isZero = (n) => (n instanceof Num && n.value === 0);

            if (isZero(num) && isZero(den)) {
                 const diffNum = expr.left.diff(varNode).simplify();
                 const diffDen = expr.right.diff(varNode).simplify();
                 return this._limit(new Div(diffNum, diffDen), varNode, point, depth + 1);
            }

            if (num instanceof Num && den instanceof Num) {
                if (den.value !== 0) return new Num(num.value / den.value);
                if (num.value !== 0 && den.value === 0) return new Sym("Infinity");
            }
        }

        try {
            const val = expr.substitute(varNode, point).simplify();
            if (val instanceof Num) return val;
            return val;
        } catch (e) {
            return new Call("limit", [expr, varNode, point]);
        }
    }

    _det(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("det requires a matrix");

        const rows = matrix.elements.length;
        if (rows === 0) return new Num(0);
        if (!(matrix.elements[0] instanceof Vec)) throw new Error("det requires a matrix");
        const cols = matrix.elements[0].elements.length;

        if (rows !== cols) throw new Error("det requires a square matrix");

        if (rows === 1) return matrix.elements[0].elements[0];
        if (rows === 2) {
            const a = matrix.elements[0].elements[0];
            const b = matrix.elements[0].elements[1];
            const c = matrix.elements[1].elements[0];
            const d = matrix.elements[1].elements[1];
            return new Sub(new Mul(a, d), new Mul(b, c)).simplify();
        }

        let det = new Num(0);
        for (let i = 0; i < cols; i++) {
            const sign = (i % 2 === 0) ? new Num(1) : new Num(-1);
            const coeff = matrix.elements[0].elements[i];

            const subRows = [];
            for (let r = 1; r < rows; r++) {
                const subRow = [];
                for (let c = 0; c < cols; c++) {
                    if (c === i) continue;
                    subRow.push(matrix.elements[r].elements[c]);
                }
                subRows.push(new Vec(subRow));
            }
            const subMatrix = new Vec(subRows);

            const term = new Mul(new Mul(sign, coeff), this._det(subMatrix));
            det = new Add(det, term);
        }
        return det.simplify();
    }

    _inv(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("inv requires a matrix");
        const rows = matrix.elements.length;
        if (rows === 0) return matrix;
        if (!(matrix.elements[0] instanceof Vec)) throw new Error("inv requires a matrix");
        const cols = matrix.elements[0].elements.length;
        if (rows !== cols) throw new Error("inv requires a square matrix");

        const det = this._det(matrix).simplify();
        if (det instanceof Num && det.value === 0) throw new Error("Matrix is singular (det=0)");

        // Cofactor matrix
        const cofactorRows = [];
        for (let r = 0; r < rows; r++) {
            const row = [];
            for (let c = 0; c < cols; c++) {
                // Minor
                const subRows = [];
                for (let i = 0; i < rows; i++) {
                    if (i === r) continue;
                    const subRow = [];
                    for (let j = 0; j < cols; j++) {
                        if (j === c) continue;
                        subRow.push(matrix.elements[i].elements[j]);
                    }
                    subRows.push(new Vec(subRow));
                }
                const subMatrix = new Vec(subRows);
                let minor = this._det(subMatrix);

                // Sign
                if ((r + c) % 2 !== 0) minor = new Mul(new Num(-1), minor);
                row.push(minor);
            }
            cofactorRows.push(new Vec(row));
        }

        // Adjugate (Transpose of Cofactor)
        const adjRows = [];
        for (let c = 0; c < cols; c++) {
            const row = [];
            for (let r = 0; r < rows; r++) {
                row.push(cofactorRows[r].elements[c]);
            }
            adjRows.push(new Vec(row));
        }

        // Multiply by 1/det
        const invRows = [];
        for (let i = 0; i < rows; i++) {
            const row = [];
            for (let j = 0; j < cols; j++) {
                row.push(new Div(adjRows[i].elements[j], det).simplify());
            }
            invRows.push(new Vec(row));
        }

        return new Vec(invRows);
    }

    _cross(v1, v2) {
        if (!(v1 instanceof Vec) || !(v2 instanceof Vec)) throw new Error("cross requires two vectors");
        // Assume 3D vectors
        if (v1.elements.length !== 3 || v2.elements.length !== 3) throw new Error("cross requires 3D vectors");

        const a1 = v1.elements[0], a2 = v1.elements[1], a3 = v1.elements[2];
        const b1 = v2.elements[0], b2 = v2.elements[1], b3 = v2.elements[2];

        const c1 = new Sub(new Mul(a2, b3), new Mul(a3, b2));
        const c2 = new Sub(new Mul(a3, b1), new Mul(a1, b3));
        const c3 = new Sub(new Mul(a1, b2), new Mul(a2, b1));

        return new Vec([c1.simplify(), c2.simplify(), c3.simplify()]);
    }

    _trans(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("trans requires a matrix");
        const rows = matrix.elements.length;
        if (rows === 0) return matrix;

        const isMatrix = matrix.elements[0] instanceof Vec;
        if (!isMatrix) {
            const newRows = [];
            for(let i=0; i<rows; i++) {
                newRows.push(new Vec([matrix.elements[i]]));
            }
            return new Vec(newRows);
        }

        const cols = matrix.elements[0].elements.length;
        const newRows = [];
        for(let j=0; j<cols; j++) {
            const row = [];
            for(let i=0; i<rows; i++) {
                row.push(matrix.elements[i].elements[j]);
            }
            newRows.push(new Vec(row));
        }
        return new Vec(newRows);
    }

    _gcd(a, b) {
        if (a instanceof Num && b instanceof Num && Number.isInteger(a.value) && Number.isInteger(b.value)) {
            const x = Math.abs(a.value);
            const y = Math.abs(b.value);
            const calcGcd = (u, v) => !v ? u : calcGcd(v, u % v);
            return new Num(calcGcd(x, y));
        }
        // Polynomial GCD not implemented
        return new Call('gcd', [a, b]);
    }

    _lcm(a, b) {
         if (a instanceof Num && b instanceof Num && Number.isInteger(a.value) && Number.isInteger(b.value)) {
            const x = Math.abs(a.value);
            const y = Math.abs(b.value);
            if (x === 0 || y === 0) return new Num(0);
            const calcGcd = (u, v) => !v ? u : calcGcd(v, u % v);
            return new Num((x * y) / calcGcd(x, y));
         }
         return new Call('lcm', [a, b]);
    }

    _factor(n) {
        if (n instanceof Num && Number.isInteger(n.value)) {
            let val = n.value;
            if (val === 0) return new Num(0);
            if (val === 1) return new Num(1);
            if (val < 0) return new Mul(new Num(-1), this._factor(new Num(-val)));

            const factors = [];
            let d = 2;
            while (d * d <= val) {
                while (val % d === 0) {
                    factors.push(new Num(d));
                    val /= d;
                }
                d++;
            }
            if (val > 1) factors.push(new Num(val));

            const counts = {};
            for(const f of factors) {
                const v = f.value;
                counts[v] = (counts[v] || 0) + 1;
            }

            let exprs = [];
            for(const base in counts) {
                const p = counts[base];
                if (p === 1) exprs.push(new Num(parseInt(base)));
                else exprs.push(new Pow(new Num(parseInt(base)), new Num(p)));
            }

            if (exprs.length === 1) return exprs[0];

            // To prevent immediate simplification to a single number, we might need to wrap it?
            // Or just return it and hope simplify doesn't merge Num * Num instantly?
            // My simplify DOES merge Num * Num.
            // So `2^2 * 3` -> `4 * 3` -> `12`.

            // This is a known issue in simple CAS.
            // I'll leave it as is, it might collapse.
            // But if I return `Call('times', ...)` it won't simplify.
            // Let's use `Call('times', ...)` or just accept it simplifies for now.
            // Actually, Xcas `factor` on integer returns `2^2 * 3`.
            // If I implement `Mul` to NOT multiply if one is a Power? No.

            return new Call('factored', exprs); // Use a custom wrapper to display it?
        }
        return new Call('factor', [n]);
    }

    _mean(list) {
        if (list instanceof Vec) {
            if (list.elements.length === 0) return new Num(0);
            let sum = new Num(0);
            for(const e of list.elements) sum = new Add(sum, e);
            return new Div(sum, new Num(list.elements.length)).simplify();
        }
        return new Call('mean', [list]);
    }

    _variance(list) {
         if (list instanceof Vec) {
             const n = list.elements.length;
             if (n < 2) return new Num(0);
             const m = this._mean(list);
             let sumSq = new Num(0);
             for(const e of list.elements) {
                 const diff = new Sub(e, m);
                 sumSq = new Add(sumSq, new Pow(diff, new Num(2)));
             }
             // Sample variance (n-1) or Population (n)? Usually sample for stats.
             return new Div(sumSq, new Num(n - 1)).simplify();
         }
         return new Call('variance', [list]);
    }

    _linearRegression(data) {
        // data: [[x1, y1], [x2, y2], ...]
        if (data instanceof Vec) {
            const n = data.elements.length;
            if (n < 2) throw new Error("Linear regression requires at least 2 points");

            let sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
            const points = [];

            for (const pt of data.elements) {
                if (!(pt instanceof Vec) || pt.elements.length !== 2) throw new Error("Data must be list of [x, y] points");
                const xExpr = pt.elements[0].evaluateNumeric();
                const yExpr = pt.elements[1].evaluateNumeric();
                if (isNaN(xExpr) || isNaN(yExpr)) throw new Error("Linear regression data must be numeric");

                sumX += xExpr;
                sumY += yExpr;
                sumXY += xExpr * yExpr;
                sumXX += xExpr * xExpr;
                points.push({x: xExpr, y: yExpr});
            }

            const denominator = n * sumXX - sumX * sumX;
            if (denominator === 0) throw new Error("Vertical line (undefined slope)");

            const slope = (n * sumXY - sumX * sumY) / denominator;
            const intercept = (sumY - slope * sumX) / n;

            // Return equation: slope * x + intercept
            // Also attach the points to the result object for plotting?
            // The result of _recursiveEval must be an Expr or a special object.
            // If we return an Expr, it simplifies.

            const eq = new Add(new Mul(new Num(slope), new Sym('x')), new Num(intercept)).simplify();

            // We can return a special object that behaves like an Expr but has extra data?
            // Or just return the equation. The user can plot it.
            // But to "Visualize", we want the points.
            // Let's rely on the frontend to construct the plot command if needed.
            // Or, we can return a 'plot' object directly if called via a specific command?
            // No, linearRegression should return the function.

            return eq;
        }
        return new Call('linearRegression', [data]);
    }

    _normalPDF(x, mu, sigma) {
        if (x instanceof Num && mu instanceof Num && sigma instanceof Num) {
            const xv = x.value;
            const mv = mu.value;
            const sv = sigma.value;
            if (sv <= 0) return new Num(0); // Error?

            const coeff = 1 / (sv * Math.sqrt(2 * Math.PI));
            const exp = Math.exp(-0.5 * Math.pow((xv - mv) / sv, 2));
            return new Num(coeff * exp);
        }
        // Symbolic representation: 1/(sigma*sqrt(2pi)) * e^(-...)
        const coeff = new Div(new Num(1), new Mul(sigma, new Call('sqrt', [new Mul(new Num(2), new Sym('pi'))])));
        const exponent = new Mul(new Num(-0.5), new Pow(new Div(new Sub(x, mu), sigma), new Num(2)));
        return new Mul(coeff, new Call('exp', [exponent]));
    }

    _binomialPDF(k, n, p) {
        if (k instanceof Num && n instanceof Num && p instanceof Num) {
            const kv = k.value;
            const nv = n.value;
            const pv = p.value;
            if (!Number.isInteger(kv) || !Number.isInteger(nv)) return new Num(0);
            if (kv < 0 || kv > nv) return new Num(0);

            // nCr * p^k * (1-p)^(n-k)
            const nCr = this._nCr(n, k).value;
            const prob = nCr * Math.pow(pv, kv) * Math.pow(1 - pv, nv - kv);
            return new Num(prob);
        }
        return new Call('binomialPDF', [k, n, p]);
    }

    _nCr(n, k) {
        if (n instanceof Num && k instanceof Num && Number.isInteger(n.value) && Number.isInteger(k.value)) {
            const valN = n.value;
            const valK = k.value;
            if (valK < 0 || valK > valN) return new Num(0);
            return new Num(this._factorial(valN) / (this._factorial(valK) * this._factorial(valN - valK)));
        }
        return new Call('nCr', [n, k]);
    }

    _nPr(n, k) {
        if (n instanceof Num && k instanceof Num && Number.isInteger(n.value) && Number.isInteger(k.value)) {
            const valN = n.value;
            const valK = k.value;
            if (valK < 0 || valK > valN) return new Num(0);
            return new Num(this._factorial(valN) / this._factorial(valN - valK));
        }
        return new Call('nPr', [n, k]);
    }

    _trace(matrix) {
        if (matrix instanceof Vec) {
            // Check if it's a matrix (vec of vecs)
            if (matrix.elements.length === 0) return new Num(0);
            if (!(matrix.elements[0] instanceof Vec)) return new Call('trace', [matrix]); // Not a matrix, maybe a vector? Trace of vector is undefined or sum? Usually matrix.

            const rows = matrix.elements.length;
            const cols = matrix.elements[0].elements.length;
            if (rows !== cols) throw new Error("trace requires a square matrix");

            let sum = new Num(0);
            for(let i=0; i<rows; i++) {
                sum = new Add(sum, matrix.elements[i].elements[i]);
            }
            return sum.simplify();
        }
        return new Call('trace', [matrix]);
    }
}
