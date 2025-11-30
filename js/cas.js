
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
        // Handle explicit assignment with '=' (Eq) at top level
        // e.g. "myvar = 42" parses as Eq(Sym(myvar), Num(42))
        // We convert this to Assignment(Sym(myvar), Num(42))
        if (exprTree instanceof Eq && exprTree.left instanceof Sym) {
            exprTree = new Assignment(exprTree.left, exprTree.right);
        }

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
                if (args.length === 1) return this._sumList(args[0]);
                // sum(expr, var, start, end)
                if (args.length !== 4) throw new Error("sum requires 1 argument (list) or 4 arguments: expression, variable, start, end");
                return this._sum(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'product') {
                if (args.length === 1) return this._productList(args[0]);
                // product(expr, var, start, end)
                if (args.length !== 4) throw new Error("product requires 1 argument (list) or 4 arguments: expression, variable, start, end");
                return this._product(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'cumsum') {
                if (args.length !== 1) throw new Error("cumsum requires 1 argument (list)");
                return this._cumsum(args[0]);
            }

            if (node.funcName === 'flatten') {
                if (args.length !== 1) throw new Error("flatten requires 1 argument (list)");
                return this._flatten(args[0]);
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
                 let min = args.length > 2 ? args[2].evaluateNumeric() : -10;
                 let max = args.length > 3 ? args[3].evaluateNumeric() : 10;

                 if (isNaN(min)) min = -10;
                 if (isNaN(max)) max = 10;

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

            if (node.funcName === 'jacobian') {
                if (args.length !== 2) throw new Error("jacobian requires 2 arguments: functions, variables");
                return this._jacobian(args[0], args[1]);
            }

            if (node.funcName === 'hessian') {
                if (args.length !== 2) throw new Error("hessian requires 2 arguments: expression, variables");
                return this._hessian(args[0], args[1]);
            }

            if (node.funcName === 'inv') {
                if (args.length !== 1) throw new Error("inv requires 1 argument");
                return this._inv(args[0]);
            }

            if (node.funcName === 'rref') {
                if (args.length !== 1) throw new Error("rref requires 1 argument");
                return this._rref(args[0]);
            }

            if (node.funcName === 'rank') {
                if (args.length !== 1) throw new Error("rank requires 1 argument");
                return this._rank(args[0]);
            }

            if (node.funcName === 'kernel' || node.funcName === 'nullspace') {
                if (args.length !== 1) throw new Error("kernel requires 1 argument");
                return this._kernel(args[0]);
            }

            if (node.funcName === 'lu') {
                if (args.length !== 1) throw new Error("lu requires 1 argument");
                return this._lu(args[0]);
            }

            if (node.funcName === 'qr') {
                if (args.length !== 1) throw new Error("qr requires 1 argument");
                return this._qr(args[0]);
            }

            if (node.funcName === 'eigenvals' || node.funcName === 'eig') {
                if (args.length !== 1) throw new Error("eigenvals requires 1 argument");
                return this._eigenvals(args[0]);
            }

            if (node.funcName === 'cross') {
                if (args.length !== 2) throw new Error("cross requires 2 arguments");
                return this._cross(args[0], args[1]);
            }

            if (node.funcName === 'trans' || node.funcName === 'transpose') {
                if (args.length !== 1) throw new Error("trans requires 1 argument");
                return this._trans(args[0]);
            }

            if (node.funcName === 'eye' || node.funcName === 'idn') {
                if (args.length !== 1) throw new Error("eye requires 1 argument");
                return this._identity(args[0]);
            }

            if (node.funcName === 'zeros') {
                if (args.length !== 2) throw new Error("zeros requires 2 arguments: rows, cols");
                return this._zeros(args[0], args[1]);
            }

            if (node.funcName === 'ones') {
                if (args.length !== 2) throw new Error("ones requires 2 arguments: rows, cols");
                return this._ones(args[0], args[1]);
            }

            if (node.funcName === 'binomial') {
                if (args.length !== 2) throw new Error("binomial requires 2 arguments");
                return this._nCr(args[0], args[1]);
            }

            if (node.funcName === 'divisors') {
                if (args.length !== 1) throw new Error("divisors requires 1 argument");
                return this._divisors(args[0]);
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

            if (node.funcName === 'isPrime') {
                if (args.length !== 1) throw new Error("isPrime requires 1 argument");
                return this._isPrime(args[0]);
            }

            if (node.funcName === 'trace') {
                if (args.length !== 1) throw new Error("trace requires 1 argument");
                return this._trace(args[0]);
            }

            if (node.funcName === 'mean') {
                if (args.length !== 1) throw new Error("mean requires 1 argument (list)");
                return this._mean(args[0]);
            }

            if (node.funcName === 'variance' || node.funcName === 'var') {
                if (args.length !== 1) throw new Error("variance requires 1 argument (list)");
                return this._variance(args[0]);
            }

            if (node.funcName === 'std' || node.funcName === 'stddev') {
                if (args.length !== 1) throw new Error("std requires 1 argument (list)");
                return this._std(args[0]);
            }

            if (node.funcName === 'cov') {
                if (args.length !== 2) throw new Error("cov requires 2 arguments (list1, list2)");
                return this._cov(args[0], args[1]);
            }

            if (node.funcName === 'corr') {
                if (args.length !== 2) throw new Error("corr requires 2 arguments (list1, list2)");
                return this._corr(args[0], args[1]);
            }

            if (node.funcName === 'median') {
                if (args.length !== 1) throw new Error("median requires 1 argument (list)");
                return this._median(args[0]);
            }

            if (node.funcName === 'charpoly') {
                if (args.length !== 2) throw new Error("charpoly requires 2 arguments: matrix, variable");
                return this._charpoly(args[0], args[1]);
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

            if (node.funcName === 'dot') {
                if (args.length !== 2) throw new Error("dot requires 2 arguments");
                // dot(u, v) is u * v (Mul handles dot product for vectors)
                return new Mul(args[0], args[1]).simplify();
            }

            if (node.funcName === 'norm') {
                if (args.length !== 1) throw new Error("norm requires 1 argument");
                // L2 norm: sqrt(v . v)
                const v = args[0];
                const dot = new Mul(v, v).simplify();
                return new Call('sqrt', [dot]).simplify();
            }

            if (node.funcName === 'grad') {
                if (args.length !== 2) throw new Error("grad requires 2 arguments: expression and list of variables");
                return this._grad(args[0], args[1]);
            }

            if (node.funcName === 'curl') {
                if (args.length !== 2) throw new Error("curl requires 2 arguments: vector field and list of variables");
                return this._curl(args[0], args[1]);
            }

            if (node.funcName === 'divergence' || node.funcName === 'div') { // 'div' might conflict with division if not careful, but funcName is safe
                 if (args.length !== 2) throw new Error("divergence requires 2 arguments: vector field and list of variables");
                 return this._divergence(args[0], args[1]);
            }

            if (node.funcName === 'rem') {
                if (args.length !== 2) throw new Error("rem requires 2 arguments");
                return this._rem(args[0], args[1]);
            }

            if (node.funcName === 'quo') {
                if (args.length !== 2) throw new Error("quo requires 2 arguments");
                return this._quo(args[0], args[1]);
            }

            if (node.funcName === 'mod') {
                if (args.length !== 2) throw new Error("mod requires 2 arguments");
                return this._mod(args[0], args[1]);
            }

            if (node.funcName === 'size' || node.funcName === 'dim') {
                if (args.length !== 1) throw new Error("size/dim requires 1 argument");
                if (args[0] instanceof Vec) {
                    return new Num(args[0].elements.length);
                }
                return new Call(node.funcName, args);
            }

            if (node.funcName === 'concat') {
                if (args.length < 2) throw new Error("concat requires at least 2 arguments");
                return this._concat(args);
            }

            if (node.funcName === 'approx') {
                if (args.length !== 1) throw new Error("approx requires 1 argument");
                return new Num(args[0].evaluateNumeric());
            }

            if (node.funcName === 'arg') {
                if (args.length !== 1) throw new Error("arg requires 1 argument");
                return this._arg(args[0]);
            }

            if (node.funcName === 'degree') {
                if (args.length !== 2) throw new Error("degree requires 2 arguments: expr, var");
                return this._degree(args[0], args[1]);
            }
            if (node.funcName === 'coeff') {
                if (args.length !== 3) throw new Error("coeff requires 3 arguments: expr, var, degree");
                return this._coeff(args[0], args[1], args[2]);
            }
            if (node.funcName === 'symb2poly') {
                if (args.length !== 2) throw new Error("symb2poly requires 2 arguments: expr, var");
                return this._symb2poly(args[0], args[1]);
            }
            if (node.funcName === 'poly2symb') {
                if (args.length !== 2) throw new Error("poly2symb requires 2 arguments: list, var");
                return this._poly2symb(args[0], args[1]);
            }
            if (node.funcName === 'seq') {
                if (args.length !== 5) throw new Error("seq requires 5 arguments: expr, var, start, end, step");
                return this._seq(args[0], args[1], args[2], args[3], args[4]);
            }
            if (node.funcName === 'range') {
                if (args.length !== 3) throw new Error("range requires 3 arguments: start, end, step");
                return this._range(args[0], args[1], args[2]);
            }
            if (node.funcName === 'sort') {
                if (args.length !== 1) throw new Error("sort requires 1 argument: list");
                return this._sort(args[0]);
            }
            if (node.funcName === 'reverse') {
                if (args.length !== 1) throw new Error("reverse requires 1 argument: list");
                return this._reverse(args[0]);
            }
            if (node.funcName === 'diag') {
                if (args.length !== 1) throw new Error("diag requires 1 argument: list");
                return this._diag(args[0]);
            }
            if (node.funcName === 'identity') {
                if (args.length !== 1) throw new Error("identity requires 1 argument: n");
                return this._identity(args[0]);
            }
            if (node.funcName === 'laplace') {
                if (args.length !== 3) throw new Error("laplace requires 3 arguments: expr, t, s");
                return this._laplace(args[0], args[1], args[2]);
            }
            if (node.funcName === 'ilaplace') {
                if (args.length !== 3) throw new Error("ilaplace requires 3 arguments: expr, s, t");
                return this._ilaplace(args[0], args[1], args[2]);
            }

            if (node.funcName === 'help') {
                const helpText = `Available commands:
diff, integrate, limit, taylor, sum, product,
expand, simplify, solve,
det, inv, trans, cross, dot, norm, grad, curl, divergence,
gcd, lcm, factor, nCr, nPr, isPrime, factorial,
mean, median, variance, linearRegression, normalPDF, binomialPDF,
degree, coeff, symb2poly, poly2symb,
seq, range, sort, reverse,
diag, identity,
laplace, ilaplace,
rem, quo, mod, arg, approx,
size, concat, clear, N`;

                const latexHelp = `\\text{Available commands: see documentation}`;
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

    _sumList(list) {
        if (list instanceof Vec) {
            let sum = new Num(0);
            for(const e of list.elements) sum = new Add(sum, e);
            return sum.simplify();
        }
        return new Call('sum', [list]);
    }

    _productList(list) {
        if (list instanceof Vec) {
            let prod = new Num(1);
            for(const e of list.elements) prod = new Mul(prod, e);
            return prod.simplify();
        }
        return new Call('product', [list]);
    }

    _cumsum(list) {
        if (list instanceof Vec) {
            let sum = new Num(0);
            const res = [];
            for(const e of list.elements) {
                sum = new Add(sum, e).simplify();
                res.push(sum);
            }
            return new Vec(res);
        }
        return new Call('cumsum', [list]);
    }

    _flatten(list) {
        if (list instanceof Vec) {
            const res = [];
            const recurse = (v) => {
                if (v instanceof Vec) {
                    for(const e of v.elements) recurse(e);
                } else {
                    res.push(v);
                }
            };
            recurse(list);
            return new Vec(res);
        }
        return new Call('flatten', [list]);
    }

    _std(list) {
        // sqrt(variance)
        const v = this._variance(list);
        return new Call('sqrt', [v]).simplify();
    }

    _cov(list1, list2) {
        if (list1 instanceof Vec && list2 instanceof Vec) {
             const n = list1.elements.length;
             if (n !== list2.elements.length) throw new Error("cov requires lists of equal length");
             if (n < 2) return new Num(0);

             const m1 = this._mean(list1);
             const m2 = this._mean(list2);

             let sum = new Num(0);
             for(let i=0; i<n; i++) {
                 const diff1 = new Sub(list1.elements[i], m1);
                 const diff2 = new Sub(list2.elements[i], m2);
                 sum = new Add(sum, new Mul(diff1, diff2));
             }
             // Sample covariance (n-1)
             return new Div(sum, new Num(n - 1)).simplify();
        }
        return new Call('cov', [list1, list2]);
    }

    _corr(list1, list2) {
         if (list1 instanceof Vec && list2 instanceof Vec) {
             const cov = this._cov(list1, list2);
             const std1 = this._std(list1);
             const std2 = this._std(list2);
             if ((std1 instanceof Num && std1.value === 0) || (std2 instanceof Num && std2.value === 0)) {
                 return new Num(0); // Undefined if std is 0
             }
             return new Div(cov, new Mul(std1, std2)).simplify();
         }
         return new Call('corr', [list1, list2]);
    }

    _charpoly(matrix, varNode) {
        if (!(matrix instanceof Vec)) throw new Error("charpoly requires a matrix");
        if (!(varNode instanceof Sym)) throw new Error("charpoly requires a symbol variable");

        // det(M - lambda * I)
        const rows = matrix.elements.length;
        if (rows === 0) return new Num(0);
        const cols = matrix.elements[0].elements.length;
        if (rows !== cols) throw new Error("charpoly requires a square matrix");

        const M_minus_lambdaI = [];
        for(let i=0; i<rows; i++) {
            const row = [];
            for(let j=0; j<cols; j++) {
                let val = matrix.elements[i].elements[j];
                if (i === j) {
                    // val - lambda
                    val = new Sub(val, varNode).simplify();
                }
                row.push(val);
            }
            M_minus_lambdaI.push(new Vec(row));
        }
        const mat = new Vec(M_minus_lambdaI);
        return this._det(mat);
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
                if (num.value !== 0 && den.value === 0) {
                     if (num.value < 0) return new Mul(new Num(-1), new Sym("Infinity"));
                     return new Sym("Infinity");
                }
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

    _rref(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("rref requires a matrix");
        const rows = matrix.elements.length;
        if (rows === 0) return matrix;
        if (!(matrix.elements[0] instanceof Vec)) throw new Error("rref requires a matrix");
        const cols = matrix.elements[0].elements.length;

        // Clone matrix to avoid modifying original
        const M = [];
        for(let i=0; i<rows; i++) {
            const row = [];
            for(let j=0; j<cols; j++) {
                row.push(matrix.elements[i].elements[j]);
            }
            M.push(row);
        }

        let lead = 0;
        for (let r = 0; r < rows; r++) {
            if (cols <= lead) break;
            let i = r;
            while (M[i][lead].evaluateNumeric() === 0) {
                i++;
                if (rows === i) {
                    i = r;
                    lead++;
                    if (cols === lead) return new Vec(M.map(row => new Vec(row)));
                }
            }

            // Swap rows i and r
            const temp = M[i];
            M[i] = M[r];
            M[r] = temp;

            const val = M[r][lead];
            // Divide row r by val
            for (let j = 0; j < cols; j++) {
                M[r][j] = new Div(M[r][j], val).simplify();
            }

            for (let i = 0; i < rows; i++) {
                if (i !== r) {
                    const sub = M[i][lead];
                    for (let j = 0; j < cols; j++) {
                        // M[i][j] = M[i][j] - sub * M[r][j]
                        M[i][j] = new Sub(M[i][j], new Mul(sub, M[r][j])).simplify();
                    }
                }
            }
            lead++;
        }

        return new Vec(M.map(row => new Vec(row)));
    }

    _rank(matrix) {
        // Rank is the number of non-zero rows in RREF
        const rref = this._rref(matrix);
        let rank = 0;
        for(const row of rref.elements) {
            let isZero = true;
            for(const el of row.elements) {
                const val = el.evaluateNumeric();
                if (isNaN(val) || Math.abs(val) > 1e-10) {
                    isZero = false;
                    break;
                }
            }
            if (!isZero) rank++;
        }
        return new Num(rank);
    }

    _kernel(matrix) {
         if (!(matrix instanceof Vec)) throw new Error("kernel requires a matrix");
         // Calculate nullspace basis using RREF
         // Ax = 0. RREF(A)x = 0.
         // Identify free variables (columns without pivots).
         const rref = this._rref(matrix);
         const rows = rref.elements.length;
         if (rows === 0) return new Vec([]);
         const cols = rref.elements[0].elements.length;

         // Identify pivots
         const pivots = []; // [col_index] for each row
         const pivotCols = new Set();
         for(let i=0; i<rows; i++) {
             let found = false;
             for(let j=0; j<cols; j++) {
                 const val = rref.elements[i].elements[j].evaluateNumeric();
                 if (Math.abs(val) > 1e-10) { // Found pivot (approx non-zero)
                      pivots.push({row: i, col: j});
                      pivotCols.add(j);
                      found = true;
                      break;
                 }
             }
             if (!found) break; // Zero rows
         }

         const freeCols = [];
         for(let j=0; j<cols; j++) {
             if (!pivotCols.has(j)) freeCols.push(j);
         }

         // Construct basis vectors
         // For each free variable x_f, set x_f = 1, other free vars = 0.
         // Solve for pivot vars.
         const basis = [];
         for(const freeCol of freeCols) {
             const vec = new Array(cols).fill(new Num(0));
             vec[freeCol] = new Num(1);

             // Back substitution for pivots
             // Row i corresponds to pivot at pivots[i].col
             // eq: x_{pivot} + sum(A_ik * x_k) = 0
             // x_{pivot} = - sum(A_ik * x_k)
             for(let i = pivots.length - 1; i >= 0; i--) {
                 const pCol = pivots[i].col;
                 const pRow = pivots[i].row;
                 // Sum over columns > pCol (which include free vars and already computed pivots)
                 // Actually, in RREF, pivots are only dependent on free vars to the right?
                 // RREF: 1 0 2 0 3
                 //       0 1 4 0 5
                 // Pivot cols 0, 1. Free cols 2, 4.
                 // Row 0: x0 + 2x2 + 3x4 = 0 => x0 = -2x2 - 3x4
                 let sum = new Num(0);
                 for(let j = pCol + 1; j < cols; j++) {
                      const coef = rref.elements[pRow].elements[j];
                      sum = new Add(sum, new Mul(coef, vec[j])).simplify();
                 }
                 vec[pCol] = new Mul(new Num(-1), sum).simplify();
             }
             basis.push(new Vec(vec));
         }

         return new Vec(basis);
    }

    _lu(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("lu requires a matrix");
        const n = matrix.elements.length;
        if (n === 0) throw new Error("Empty matrix");
        const m = matrix.elements[0].elements.length;
        if (n !== m) throw new Error("lu requires a square matrix");

        // Doolittle's Algorithm? Or Gaussian elimination.
        // Returns P, L, U such that PA = LU.
        // Implementing simple LU without pivoting first? Or with P?
        // Xcas `lu` returns [L, U, P].

        // Initialize L = I, U = A, P = I
        const L = [];
        const U = [];
        const P = [];

        // Clone U from matrix
        for(let i=0; i<n; i++) {
            const rowU = [];
            const rowL = [];
            const rowP = [];
            for(let j=0; j<n; j++) {
                rowU.push(matrix.elements[i].elements[j]);
                rowL.push(new Num(i === j ? 1 : 0));
                rowP.push(new Num(i === j ? 1 : 0));
            }
            U.push(rowU);
            L.push(rowL);
            P.push(rowP);
        }

        // Gaussian Elimination
        for(let k=0; k<n; k++) {
             // Pivot strategy: find max in column k, rows k..n-1
             let maxVal = -1;
             let pivotRow = -1;
             for(let i=k; i<n; i++) {
                  const val = Math.abs(U[i][k].evaluateNumeric());
                  if (val > maxVal) {
                      maxVal = val;
                      pivotRow = i;
                  }
             }

             if (maxVal < 1e-10) continue; // Singular or zero column

             // Swap rows in U, P. Swap rows in L (only up to k-1 to keep L lower triangular?)
             // Actually, for PA=LU, we swap rows in U, P, and "the part of L computed so far" (columns 0..k-1).
             if (pivotRow !== k) {
                  // Swap U rows
                  let temp = U[k]; U[k] = U[pivotRow]; U[pivotRow] = temp;
                  // Swap P rows
                  temp = P[k]; P[k] = P[pivotRow]; P[pivotRow] = temp;
                  // Swap L rows (elements 0..k-1)
                  for(let j=0; j<k; j++) {
                       let t = L[k][j]; L[k][j] = L[pivotRow][j]; L[pivotRow][j] = t;
                  }
             }

             for(let i=k+1; i<n; i++) {
                  // factor = U[i][k] / U[k][k]
                  const factor = new Div(U[i][k], U[k][k]).simplify();
                  L[i][k] = factor;
                  for(let j=k; j<n; j++) {
                       // U[i][j] = U[i][j] - factor * U[k][j]
                       U[i][j] = new Sub(U[i][j], new Mul(factor, U[k][j])).simplify();
                  }
             }
        }

        const vecL = new Vec(L.map(r => new Vec(r)));
        const vecU = new Vec(U.map(r => new Vec(r)));
        const vecP = new Vec(P.map(r => new Vec(r)));

        // Xcas returns permutation p as a vector of indices? Or matrix?
        // Documentation says [L, U, P] where P is permutation matrix.
        return new Vec([vecL, vecU, vecP]);
    }

    _qr(matrix) {
         if (!(matrix instanceof Vec)) throw new Error("qr requires a matrix");
         const rows = matrix.elements.length;
         if (rows === 0) throw new Error("Empty matrix");
         const cols = matrix.elements[0].elements.length;

         // Gram-Schmidt
         // Q cols are orthonormal basis. R is upper triangular.
         // A = [a1, ..., an]
         // u1 = a1, e1 = u1/|u1|
         // u2 = a2 - proj_u1(a2), e2 = u2/|u2|
         // ...

         // Transpose matrix to work with columns easily?
         const A_cols = [];
         for(let j=0; j<cols; j++) {
              const col = [];
              for(let i=0; i<rows; i++) col.push(matrix.elements[i].elements[j]);
              A_cols.push(new Vec(col));
         }

         const Q_cols = [];
         const R_elements = []; // Matrix R

         for(let j=0; j<cols; j++) {
              let u = A_cols[j];
              // Subtract projections onto previous e_i
              for(let i=0; i<j; i++) {
                  // proj = (a_j . e_i) * e_i
                  // dot
                  const dot = new Mul(A_cols[j], Q_cols[i]).simplify(); // Vector dot product via Mul
                  // R_ij = dot
                  // We need to store R properly.
                  // u = u - dot * e_i
                  u = new Sub(u, new Mul(dot, Q_cols[i])).simplify();
              }

              // Norm
              const normSq = new Mul(u, u).simplify();
              const norm = new Call('sqrt', [normSq]).simplify();
              const e = new Div(u, norm).simplify();
              Q_cols.push(e);
         }

         // Construct Q matrix from columns
         const Q_rows = [];
         for(let i=0; i<rows; i++) {
              const row = [];
              for(let j=0; j<cols; j++) {
                   // e_j[i]
                   row.push(Q_cols[j].elements[i]);
              }
              Q_rows.push(new Vec(row));
         }

         // Construct R matrix
         // R_ij = e_i . a_j (for i <= j)
         const R_rows = [];
         for(let i=0; i<cols; i++) {
              const row = [];
              for(let j=0; j<cols; j++) {
                   if (i > j) {
                        row.push(new Num(0));
                   } else {
                        const val = new Mul(Q_cols[i], A_cols[j]).simplify();
                        row.push(val);
                   }
              }
              R_rows.push(new Vec(row));
         }

         return new Vec([new Vec(Q_rows), new Vec(R_rows)]);
    }

    _eigenvals(matrix) {
         if (!(matrix instanceof Vec)) throw new Error("eigenvals requires a matrix");
         // solve(charpoly(M, x) = 0, x)
         const x = new Sym("lambda_" + Date.now()); // Unique variable
         const cp = this._charpoly(matrix, x);
         const roots = this._solve(cp, x);
         // roots might be a "set" call or single value.
         // We want a list (Vec).
         if (roots instanceof Call && roots.funcName === 'set') {
             return new Vec(roots.args);
         }
         return new Vec([roots]);
    }

    _zeros(rows, cols) {
        if (!(rows instanceof Num) || !(cols instanceof Num)) return new Call('zeros', [rows, cols]);
        const r = rows.value;
        const c = cols.value;
        const res = [];
        for(let i=0; i<r; i++) {
            const row = [];
            for(let j=0; j<c; j++) {
                row.push(new Num(0));
            }
            res.push(new Vec(row));
        }
        return new Vec(res);
    }

    _ones(rows, cols) {
        if (!(rows instanceof Num) || !(cols instanceof Num)) return new Call('ones', [rows, cols]);
        const r = rows.value;
        const c = cols.value;
        const res = [];
        for(let i=0; i<r; i++) {
            const row = [];
            for(let j=0; j<c; j++) {
                row.push(new Num(1));
            }
            res.push(new Vec(row));
        }
        return new Vec(res);
    }

    _divisors(n) {
        n = n.simplify();
        if (n instanceof Num && Number.isInteger(n.value)) {
            const val = Math.abs(n.value);
            const res = [];
            for(let i=1; i*i <= val; i++) {
                if (val % i === 0) {
                    res.push(new Num(i));
                    if (i*i !== val) res.push(new Num(val/i));
                }
            }
            // Sort
            res.sort((a, b) => a.value - b.value);
            return new Vec(res);
        }
        return new Call('divisors', [n]);
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
        a = a.simplify();
        b = b.simplify();

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
        a = a.simplify();
        b = b.simplify();

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
        n = n.simplify();

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
                else exprs.push(new Call('power', [new Num(parseInt(base)), new Num(p)]));
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

    _median(list) {
        if (list instanceof Vec) {
            const n = list.elements.length;
            if (n === 0) return new Num(0);

            // Sort elements based on numeric value
            const sorted = [...list.elements].sort((a, b) => {
                const va = a.evaluateNumeric();
                const vb = b.evaluateNumeric();
                return va - vb;
            });

            if (n % 2 === 1) {
                return sorted[Math.floor(n / 2)];
            } else {
                const mid1 = sorted[n / 2 - 1];
                const mid2 = sorted[n / 2];
                return new Div(new Add(mid1, mid2), new Num(2)).simplify();
            }
        }
        return new Call('median', [list]);
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
        x = x.simplify();
        mu = mu.simplify();
        sigma = sigma.simplify();

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
        k = k.simplify();
        n = n.simplify();
        p = p.simplify();

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
        n = n.simplify();
        k = k.simplify();

        if (n instanceof Num && k instanceof Num && Number.isInteger(n.value) && Number.isInteger(k.value)) {
            const valN = n.value;
            const valK = k.value;
            if (valK < 0 || valK > valN) return new Num(0);
            return new Num(this._factorial(valN) / (this._factorial(valK) * this._factorial(valN - valK)));
        }
        return new Call('nCr', [n, k]);
    }

    _nPr(n, k) {
        n = n.simplify();
        k = k.simplify();

        if (n instanceof Num && k instanceof Num && Number.isInteger(n.value) && Number.isInteger(k.value)) {
            const valN = n.value;
            const valK = k.value;
            if (valK < 0 || valK > valN) return new Num(0);
            return new Num(this._factorial(valN) / this._factorial(valN - valK));
        }
        return new Call('nPr', [n, k]);
    }

    _isPrime(n) {
        n = n.simplify();

        if (n instanceof Num && Number.isInteger(n.value)) {
            const val = n.value;
            if (val < 2) return new Num(0);
            for (let i = 2, s = Math.sqrt(val); i <= s; i++) {
                if (val % i === 0) return new Num(0);
            }
            return new Num(1);
        }
        return new Call('isPrime', [n]);
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

    _grad(expr, vars) {
        if (!(vars instanceof Vec)) throw new Error("Second argument to grad must be a vector of variables");
        const elements = vars.elements.map(v => expr.diff(v).simplify());
        return new Vec(elements);
    }

    _divergence(v, vars) {
        if (!(v instanceof Vec) || !(vars instanceof Vec)) throw new Error("Arguments to divergence must be vectors");
        if (v.elements.length !== vars.elements.length) throw new Error("Vector dimension mismatch in divergence");

        let sum = new Num(0);
        for(let i=0; i<v.elements.length; i++) {
            sum = new Add(sum, v.elements[i].diff(vars.elements[i])).simplify();
        }
        return sum;
    }

    _curl(v, vars) {
        if (!(v instanceof Vec) || !(vars instanceof Vec)) throw new Error("Arguments to curl must be vectors");
        if (v.elements.length !== 3 || vars.elements.length !== 3) throw new Error("Curl requires 3D vectors");

        // curl F = (dFz/dy - dFy/dz, dFx/dz - dFz/dx, dFy/dx - dFx/dy)
        // v = [Fx, Fy, Fz], vars = [x, y, z]
        const Fx = v.elements[0], Fy = v.elements[1], Fz = v.elements[2];
        const x = vars.elements[0], y = vars.elements[1], z = vars.elements[2];

        const c1 = new Sub(Fz.diff(y), Fy.diff(z)).simplify();
        const c2 = new Sub(Fx.diff(z), Fz.diff(x)).simplify();
        const c3 = new Sub(Fy.diff(x), Fx.diff(y)).simplify();

        return new Vec([c1, c2, c3]);
    }

    _jacobian(funcs, vars) {
        if (!(funcs instanceof Vec) || !(vars instanceof Vec)) throw new Error("jacobian requires two vectors");
        // J_ij = d(f_i)/d(x_j)
        const rows = [];
        for(const f of funcs.elements) {
            const row = [];
            for(const v of vars.elements) {
                row.push(f.diff(v).simplify());
            }
            rows.push(new Vec(row));
        }
        return new Vec(rows);
    }

    _hessian(expr, vars) {
        if (!(vars instanceof Vec)) throw new Error("hessian requires a vector of variables");
        // H_ij = d^2f / dx_i dx_j
        const rows = [];
        for(const vi of vars.elements) {
            const row = [];
            for(const vj of vars.elements) {
                // partial f / partial vi partial vj
                // first diff wrt vi, then vj
                row.push(expr.diff(vi).diff(vj).simplify());
            }
            rows.push(new Vec(row));
        }
        return new Vec(rows);
    }

    _rem(a, b) {
        a = a.simplify();
        b = b.simplify();

        a = a.simplify(); b = b.simplify(); if (a instanceof Num && b instanceof Num) {
             return new Num(a.value % b.value);
        }
        // Polynomial remainder
        // rem(P, Q)
        // Check if b is linear (x-c) handled in Div.simplify via polyDiv
        // But we want explicit rem.
        // We can use division: rem = a - b * quo(a, b)
        // But we need quo first.

        // Simple case: b is a Number. rem(P, c)
        // If coefficients are integers, rem(P, c) might mean rem of coeffs?
        // Usually rem(P, Q) is polynomial remainder.

        // For cleanroom, if we don't have full poly div, we return symbolic Call.
        return new Call('rem', [a, b]);
    }

    _quo(a, b) {
        a = a.simplify();
        b = b.simplify();

        a = a.simplify(); b = b.simplify(); if (a instanceof Num && b instanceof Num) {
             return new Num(Math.trunc(a.value / b.value));
        }
        // Polynomial quotient
        // If b divides a exactly, Div(a, b).simplify() returns it.
        // But quo(a, b) should return quotient even if remainder != 0.

        // If b is linear (x-c), we can use synthetic division logic from expression.js
        // But that is internal to Div.simplify.
        // We might want to expose polyDiv or reimplement it here.

        return new Call('quo', [a, b]);
    }

    _mod(a, b) {
        a = a.simplify();
        b = b.simplify();

        a = a.simplify(); b = b.simplify(); if (a instanceof Num && b instanceof Num) {
             // Mathematical modulo (always positive for positive modulus)
             const m = b.value;
             return new Num(((a.value % m) + m) % m);
        }
        return new Call('mod', [a, b]);
    }

    _concat(args) {
        let elements = [];
        for(const arg of args) {
            if (arg instanceof Vec) {
                elements = elements.concat(arg.elements);
            } else {
                elements.push(arg);
            }
        }
        return new Vec(elements);
    }

    _arg(z) {
        z = z.simplify();

        if (z instanceof Num) {
             if (z.value >= 0) return new Num(0);
             return new Num(Math.PI);
        }
        return new Call('arg', [z]);
    }

    // --- Polynomial Tools ---

    _dependsOn(expr, varNode) {
        if (expr instanceof Sym) return expr.name === varNode.name;
        if (expr instanceof Num) return false;
        if (expr instanceof BinaryOp) return this._dependsOn(expr.left, varNode) || this._dependsOn(expr.right, varNode);
        if (expr instanceof Call) return expr.args.some(a => this._dependsOn(a, varNode));
        if (expr instanceof Vec) return expr.elements.some(e => this._dependsOn(e, varNode));
        return false;
    }

    _getPolyCoeffs(expr, varNode) {
        expr = expr.expand().simplify();
        const terms = [];
        const collectTerms = (e, sign = 1) => {
             if (e instanceof Add) {
                 collectTerms(e.left, sign);
                 collectTerms(e.right, sign);
             } else if (e instanceof Sub) {
                 collectTerms(e.left, sign);
                 collectTerms(e.right, -sign);
             } else {
                 terms.push({e, sign});
             }
        };
        collectTerms(expr);

        const coeffs = {};
        let maxDeg = 0;

        // Analyze term
        const analyze = (t) => {
             if (t instanceof Sym && t.name === varNode.name) return { d: 1, c: new Num(1) };
             if (t instanceof Num) return { d: 0, c: t };
             if (t instanceof Sym) return { d: 0, c: t };
             if (t instanceof Pow) {
                  if (t.left instanceof Sym && t.left.name === varNode.name && t.right instanceof Num) {
                       return { d: t.right.value, c: new Num(1) };
                  }
                  if (this._dependsOn(t, varNode)) throw new Error("Not a polynomial");
                  return { d: 0, c: t };
             }
             if (t instanceof Mul) {
                  const l = analyze(t.left);
                  const r = analyze(t.right);
                  return { d: l.d + r.d, c: new Mul(l.c, r.c) };
             }
             if (t instanceof Div) {
                  const l = analyze(t.left);
                  const r = analyze(t.right);
                  if (r.d > 0) throw new Error("Rational function, not polynomial");
                  return { d: l.d, c: new Div(l.c, r.c) };
             }
             if (t instanceof Call) {
                  if (this._dependsOn(t, varNode)) throw new Error("Not a polynomial");
                  return { d: 0, c: t };
             }
             // Default const
             if (this._dependsOn(t, varNode)) throw new Error("Not a polynomial");
             return { d: 0, c: t };
        };

        for(const item of terms) {
             try {
                 const res = analyze(item.e);
                 const finalCoeff = (item.sign === 1) ? res.c : new Mul(new Num(-1), res.c);
                 coeffs[res.d] = new Add(coeffs[res.d] || new Num(0), finalCoeff);
                 if (res.d > maxDeg) maxDeg = res.d;
             } catch (e) {
                 // Not polynomial term
                 return null;
             }
        }

        // simplify coeffs
        for(const d in coeffs) coeffs[d] = coeffs[d].simplify();
        return { coeffs, maxDeg };
    }

    _degree(expr, varNode) {
        const poly = this._getPolyCoeffs(expr, varNode);
        if (poly) return new Num(poly.maxDeg);
        return new Call('degree', [expr, varNode]);
    }

    _coeff(expr, varNode, n) {
        if (!(n instanceof Num)) return new Call('coeff', [expr, varNode, n]);
        const poly = this._getPolyCoeffs(expr, varNode);
        if (poly) {
            return poly.coeffs[n.value] || new Num(0);
        }
        return new Call('coeff', [expr, varNode, n]);
    }

    _symb2poly(expr, varNode) {
        const poly = this._getPolyCoeffs(expr, varNode);
        if (poly) {
            const list = [];
            for(let i=poly.maxDeg; i>=0; i--) {
                list.push(poly.coeffs[i] || new Num(0));
            }
            return new Vec(list);
        }
        return new Call('symb2poly', [expr, varNode]);
    }

    _poly2symb(list, varNode) {
        if (list instanceof Vec) {
            let res = new Num(0);
            const deg = list.elements.length - 1;
            for(let i=0; i<list.elements.length; i++) {
                const coeff = list.elements[i];
                const power = deg - i;
                let term;
                if (power === 0) term = coeff;
                else if (power === 1) term = new Mul(coeff, varNode);
                else term = new Mul(coeff, new Pow(varNode, new Num(power)));
                res = new Add(res, term);
            }
            return res.simplify();
        }
        return new Call('poly2symb', [list, varNode]);
    }

    // --- List/Vector Tools ---

    _seq(expr, varNode, start, end, step) {
        // Generates list [expr(start), expr(start+step), ...]
        if (!(start instanceof Num) || !(end instanceof Num) || !(step instanceof Num)) {
            return new Call('seq', [expr, varNode, start, end, step]);
        }
        const res = [];
        const s = step.value;
        if (s === 0) throw new Error("Step cannot be 0");
        const a = start.value;
        const b = end.value;

        for (let i = a; (s > 0 ? i <= b : i >= b); i += s) {
            res.push(expr.substitute(varNode, new Num(i)).simplify());
        }
        return new Vec(res);
    }

    _range(start, end, step) {
        if (!(start instanceof Num) || !(end instanceof Num) || !(step instanceof Num)) {
            return new Call('range', [start, end, step]);
        }
        const res = [];
        const s = step.value;
        if (s === 0) throw new Error("Step cannot be 0");
        const a = start.value;
        const b = end.value;

        // Python range behavior: inclusive start, exclusive end?
        // Xcas range: range(start, end, step) usually inclusive?
        // Python: range(0, 3) -> [0, 1, 2]
        // Xcas: range(n) -> 0..n-1
        // Xcas: range(a, b) -> a..b-1
        // Xcas: range(a, b, s)
        // I will follow Xcas/Python convention: exclusive end.

        for (let i = a; (s > 0 ? i < b : i > b); i += s) {
            res.push(new Num(i));
        }
        return new Vec(res);
    }

    _sort(list) {
        if (list instanceof Vec) {
            const sorted = [...list.elements].sort((a, b) => {
                const va = a.evaluateNumeric();
                const vb = b.evaluateNumeric();
                return va - vb;
            });
            return new Vec(sorted);
        }
        return new Call('sort', [list]);
    }

    _reverse(list) {
        if (list instanceof Vec) {
            return new Vec([...list.elements].reverse());
        }
        return new Call('reverse', [list]);
    }

    _diag(list) {
        if (list instanceof Vec) {
            const n = list.elements.length;
            const rows = [];
            for(let i=0; i<n; i++) {
                const row = [];
                for(let j=0; j<n; j++) {
                    if (i === j) row.push(list.elements[i]);
                    else row.push(new Num(0));
                }
                rows.push(new Vec(row));
            }
            return new Vec(rows);
        }
        return new Call('diag', [list]);
    }

    _identity(n) {
        if (n instanceof Num && Number.isInteger(n.value) && n.value > 0) {
            const size = n.value;
            const rows = [];
            for(let i=0; i<size; i++) {
                const row = [];
                for(let j=0; j<size; j++) {
                    row.push(new Num(i === j ? 1 : 0));
                }
                rows.push(new Vec(row));
            }
            return new Vec(rows);
        }
        return new Call('identity', [n]);
    }

    // --- Calculus Extras ---

    _laplace(expr, t, s) {
        // Simple table-based Laplace Transform
        expr = expr.expand().simplify();

        if (expr instanceof Add) {
            return new Add(this._laplace(expr.left, t, s), this._laplace(expr.right, t, s)).simplify();
        }
        if (expr instanceof Sub) {
            return new Sub(this._laplace(expr.left, t, s), this._laplace(expr.right, t, s)).simplify();
        }
        if (expr instanceof Mul) {
            // Check for constant
            if (!this._dependsOn(expr.left, t)) return new Mul(expr.left, this._laplace(expr.right, t, s)).simplify();
            if (!this._dependsOn(expr.right, t)) return new Mul(expr.right, this._laplace(expr.left, t, s)).simplify();
        }
        if (!this._dependsOn(expr, t)) {
            // L{c} = c/s
            return new Div(expr, s).simplify();
        }

        // Table
        // t^n -> n! / s^(n+1)
        if (expr instanceof Sym && expr.name === t.name) { // t
             return new Div(new Num(1), new Pow(s, new Num(2))).simplify();
        }
        if (expr instanceof Pow && expr.left instanceof Sym && expr.left.name === t.name && expr.right instanceof Num) { // t^n
             const n = expr.right.value;
             if (n >= 0) {
                 return new Div(new Num(this._factorial(n)), new Pow(s, new Num(n + 1))).simplify();
             }
        }
        // exp(at) -> 1/(s-a)
        if (expr instanceof Call && expr.funcName === 'exp') {
             const arg = expr.args[0]; // at
             // We need to extract 'a' from 'at'
             const poly = this._getPolyCoeffs(arg, t);
             if (poly && poly.maxDeg === 1 && poly.coeffs[0] && poly.coeffs[0].value === 0) {
                 const a = poly.coeffs[1];
                 return new Div(new Num(1), new Sub(s, a)).simplify();
             }
             if (arg instanceof Sym && arg.name === t.name) { // e^t
                  return new Div(new Num(1), new Sub(s, new Num(1))).simplify();
             }
        }
        // sin(at) -> a/(s^2+a^2)
        // cos(at) -> s/(s^2+a^2)
        if (expr instanceof Call && (expr.funcName === 'sin' || expr.funcName === 'cos')) {
             const arg = expr.args[0];
             let a = new Num(1);
             if (arg instanceof Mul && arg.right instanceof Sym && arg.right.name === t.name) a = arg.left;
             else if (arg instanceof Mul && arg.left instanceof Sym && arg.left.name === t.name) a = arg.right;
             else if (arg instanceof Sym && arg.name === t.name) a = new Num(1);
             else return new Call('laplace', [expr, t, s]);

             const denom = new Add(new Pow(s, new Num(2)), new Pow(a, new Num(2)));
             if (expr.funcName === 'sin') return new Div(a, denom).simplify();
             if (expr.funcName === 'cos') return new Div(s, denom).simplify();
        }

        return new Call('laplace', [expr, t, s]);
    }

    _ilaplace(expr, s, t) {
        // Inverse Laplace (Basic)
        expr = expr.expand().simplify();

        if (expr instanceof Add) return new Add(this._ilaplace(expr.left, s, t), this._ilaplace(expr.right, s, t)).simplify();
        if (expr instanceof Sub) return new Sub(this._ilaplace(expr.left, s, t), this._ilaplace(expr.right, s, t)).simplify();
        if (expr instanceof Mul) {
             if (!this._dependsOn(expr.left, s)) return new Mul(expr.left, this._ilaplace(expr.right, s, t)).simplify();
             if (!this._dependsOn(expr.right, s)) return new Mul(expr.right, this._ilaplace(expr.left, s, t)).simplify();
        }

        // Table
        // 1/s -> 1
        // 1/s^n -> t^(n-1)/(n-1)!
        // 1/(s-a) -> e^(at)
        // a/(s^2+a^2) -> sin(at)
        // s/(s^2+a^2) -> cos(at)

        // Simplify to form Num / Denom
        if (expr instanceof Div) {
             const num = expr.left;
             const den = expr.right;

             // 1/s^n
             if (den instanceof Pow && den.left instanceof Sym && den.left.name === s.name && den.right instanceof Num) {
                 const n = den.right.value;
                 return new Mul(num, new Div(new Pow(t, new Num(n-1)), new Num(this._factorial(n-1)))).simplify();
             }
             if (den instanceof Sym && den.name === s.name) { // 1/s
                 return num;
             }

             // 1/(s-a)
             if (den instanceof Sub && den.left instanceof Sym && den.left.name === s.name) {
                 const a = den.right;
                 return new Mul(num, new Call('exp', [new Mul(a, t)])).simplify();
             }
             if (den instanceof Add && den.left instanceof Sym && den.left.name === s.name) { // 1/(s+a) = 1/(s-(-a))
                 const a = new Mul(new Num(-1), den.right);
                 return new Mul(num, new Call('exp', [new Mul(a, t)])).simplify();
             }
        }

        return new Call('ilaplace', [expr, s, t]);
    }
}

// Export classes for Global/CommonJS environments
(function() {
    const exports = {
        CAS
    };
    if (typeof globalThis !== 'undefined') {
        Object.assign(globalThis, exports);
    }
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = exports;
    }
})();
