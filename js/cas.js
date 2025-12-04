
class CAS {
    constructor() {
        this.variables = {
            'pi': new Num(Math.PI),
            'true': new Num(1),
            'false': new Num(0)
        };

        this.functions = {};
    }

    evaluate(exprTree) {
        // Handle explicit assignment with '=' (Eq) at top level
        // e.g. "myvar = 42" parses as Eq(Sym(myvar), Num(42))
        // We convert this to Assignment(Sym(myvar), Num(42))
        // Xcas compatibility: '=' is equality, ':=' is assignment.
        // We disable this conversion to enforce Xcas syntax.
        // if (exprTree instanceof Eq && exprTree.left instanceof Sym) {
        //     exprTree = new Assignment(exprTree.left, exprTree.right);
        // }

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

                let order = 1;
                if (args.length >= 3) {
                    const o = args[2];
                    if (o instanceof Num && Number.isInteger(o.value)) order = o.value;
                    else return new Call('diff', args); // Symbolic order
                }

                let res = func;
                for(let i=0; i<order; i++) {
                    res = res.diff(varNode).simplify();
                }
                return res;
            }

            if (node.funcName === 'integrate' || node.funcName === 'int') {
                if (args.length < 2) throw new Error("integrate requires at least 2 arguments");
                const func = args[0];
                const varNode = args[1];
                if (!(varNode instanceof Sym)) throw new Error("Second argument to integrate must be a variable");

                if (args.length === 4) {
                    const lower = args[2];
                    const upper = args[3];
                    let indefinite = func.integrate(varNode);

                    // Attempt partfrac if direct integration failed
                    if (indefinite instanceof Call && indefinite.funcName === 'integrate') {
                        const pf = this._partfrac(func, varNode);
                        if (!(pf instanceof Call && pf.funcName === 'partfrac')) {
                             indefinite = pf.integrate(varNode);
                        }
                    }

                    // If integration failed (returned a Call to integrate), return symbolic definite integral
                    if (indefinite instanceof Call && indefinite.funcName === 'integrate') {
                        return new Call('integrate', args);
                    }

                    const valUpper = indefinite.substitute(varNode, upper);
                    const valLower = indefinite.substitute(varNode, lower);
                    return new Sub(valUpper, valLower);
                }

                let res = func.integrate(varNode);
                if (res instanceof Call && res.funcName === 'integrate') {
                     // 1. Try Integration by Substitution (Logarithmic form f'/f)
                     const subRes = this._integrateSubstitution(func, varNode);
                     if (!(subRes instanceof Call && subRes.funcName === 'integrate')) return subRes;

                     // 2. Try Integration by Parts
                     const partsRes = this._integrateByParts(func, varNode);
                     if (!(partsRes instanceof Call && partsRes.funcName === 'integrate')) return partsRes;

                     // 3. Try partial fraction decomposition
                     const pf = this._partfrac(func, varNode);
                     if (!(pf instanceof Call && pf.funcName === 'partfrac') && pf.toString() !== func.toString()) {
                         // Integrate the decomposed form
                         const res2 = pf.integrate(varNode);
                         // If improved, return it
                         if (!(res2 instanceof Call && res2.funcName === 'integrate')) {
                             return res2;
                         }
                     }
                }
                return res;
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

            if (node.funcName === 'solve' || node.funcName === 'fsolve') {
                 if (args.length < 2) throw new Error("solve requires at least 2 arguments: equation and variable");
                 const eq = args[0];
                 const varNode = args[1];
                 if (!(varNode instanceof Sym) && !(varNode instanceof Vec)) throw new Error("Second argument to solve must be a variable or list of variables");
                 return this._solve(eq, varNode);
            }

            if (node.funcName === 'resultant') {
                 if (args.length !== 3) throw new Error("resultant requires 3 arguments: poly1, poly2, variable");
                 return this._resultant(args[0], args[1], args[2]);
            }

            if (node.funcName === 'discriminant') {
                 if (args.length !== 2) throw new Error("discriminant requires 2 arguments: poly, variable");
                 return this._discriminant(args[0], args[1]);
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
                     subtype: 'function',
                     expr: expr,
                     var: varNode,
                     min: min,
                     max: max,
                     toString: () => `Plotting ${expr} from ${min} to ${max}`,
                     toLatex: () => `\\text{Plotting } ${expr.toLatex()}`
                 };
            }

            if (node.funcName === 'plotparam') {
                // plotparam([x(t), y(t)], t, min, max)
                // or plotparam(x(t), y(t), t, min, max)

                let xExpr, yExpr, tVar, min, max;

                if (args.length >= 2 && args[0] instanceof Vec) {
                    // plotparam([x, y], t, [min], [max])
                    if (args[0].elements.length !== 2) throw new Error("plotparam list must contain 2 expressions");
                    xExpr = args[0].elements[0];
                    yExpr = args[0].elements[1];
                    tVar = args[1];
                    min = args.length > 2 ? args[2].evaluateNumeric() : -10;
                    max = args.length > 3 ? args[3].evaluateNumeric() : 10;
                } else if (args.length >= 3) {
                    // plotparam(x, y, t, [min], [max])
                    xExpr = args[0];
                    yExpr = args[1];
                    tVar = args[2];
                    min = args.length > 3 ? args[3].evaluateNumeric() : -10;
                    max = args.length > 4 ? args[4].evaluateNumeric() : 10;
                } else {
                    throw new Error("plotparam requires arguments: x, y, t, [min, max]");
                }

                if (!(tVar instanceof Sym)) throw new Error("plotparam variable must be a symbol");
                if (isNaN(min)) min = -10;
                if (isNaN(max)) max = 10;

                return {
                    type: 'plot',
                    subtype: 'parametric',
                    xExpr: xExpr,
                    yExpr: yExpr,
                    var: tVar,
                    min: min,
                    max: max,
                    toString: () => `Parametric Plot (${xExpr}, ${yExpr}) t=${min}..${max}`,
                    toLatex: () => `\\text{Parametric Plot } (${xExpr.toLatex()}, ${yExpr.toLatex()})`
                };
            }

            if (node.funcName === 'plotpolar') {
                // plotpolar(r, theta, min, max)
                if (args.length < 3) throw new Error("plotpolar requires at least 3 arguments: r, theta, min, max");
                const rExpr = args[0];
                const thetaVar = args[1];
                let min = args.length > 2 ? args[2].evaluateNumeric() : 0;
                let max = args.length > 3 ? args[3].evaluateNumeric() : 2 * Math.PI;

                if (!(thetaVar instanceof Sym)) throw new Error("plotpolar variable must be a symbol");
                if (isNaN(min)) min = 0;
                if (isNaN(max)) max = 2 * Math.PI;

                return {
                    type: 'plot',
                    subtype: 'polar',
                    rExpr: rExpr,
                    var: thetaVar,
                    min: min,
                    max: max,
                    toString: () => `Polar Plot r=${rExpr} theta=${min}..${max}`,
                    toLatex: () => `\\text{Polar Plot } r=${rExpr.toLatex()}`
                };
            }

            if (node.funcName === 'plotlist') {
                if (args.length !== 1) throw new Error("plotlist requires 1 argument (list)");
                const list = args[0];
                if (!(list instanceof Vec)) throw new Error("plotlist argument must be a list");

                const points = [];
                for(let i=0; i<list.elements.length; i++) {
                    const el = list.elements[i];
                    if (el instanceof Vec && el.elements.length >= 2) {
                        // [x, y]
                        const x = el.elements[0].evaluateNumeric();
                        const y = el.elements[1].evaluateNumeric();
                        if (!isNaN(x) && !isNaN(y)) points.push({x, y});
                    } else {
                        // y value, index as x
                        const y = el.evaluateNumeric();
                        if (!isNaN(y)) points.push({x: i, y});
                    }
                }

                return {
                    type: 'plot',
                    subtype: 'list',
                    scatter: points, // Frontend uses 'scatter' property for points
                    min: 0, // Dummies
                    max: 0,
                    toString: () => `List Plot ${points.length} points`,
                    toLatex: () => `\\text{List Plot}`
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

            if (node.funcName === 'rref') {
                if (args.length !== 1) throw new Error("rref requires 1 argument");
                return this._rref(args[0]);
            }

            if (node.funcName === 'rank') {
                if (args.length !== 1) throw new Error("rank requires 1 argument");
                return this._rank(args[0]);
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

            if (node.funcName === 'binomial' || node.funcName === 'comb') {
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

            if (node.funcName === 'nextprime') {
                if (args.length !== 1) throw new Error("nextprime requires 1 argument");
                return this._nextprime(args[0]);
            }

            if (node.funcName === 'prevprime') {
                if (args.length !== 1) throw new Error("prevprime requires 1 argument");
                return this._prevprime(args[0]);
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

            if (node.funcName === 'partfrac') {
                if (args.length !== 2) throw new Error("partfrac requires 2 arguments: expr, var");
                return this._partfrac(args[0], args[1]);
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

            if (node.funcName === 'append') {
                if (args.length !== 2) throw new Error("append requires 2 arguments: list, element");
                return this._append(args[0], args[1]);
            }

            if (node.funcName === 'prepend') {
                if (args.length !== 2) throw new Error("prepend requires 2 arguments: list, element");
                return this._prepend(args[0], args[1]);
            }

            if (node.funcName === 'approx' || node.funcName === 'evalf') {
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

            if (node.funcName === 'kernel' || node.funcName === 'nullspace' || node.funcName === 'ker') {
                if (args.length !== 1) throw new Error("kernel requires 1 argument");
                return this._kernel(args[0].simplify());
            }
            if (node.funcName === 'basis') {
                if (args.length !== 1) throw new Error("basis requires 1 argument");
                return this._basis(args[0].simplify());
            }
            if (node.funcName === 'eigenvals') {
                if (args.length !== 1) throw new Error("eigenvals requires 1 argument");
                return this._eigenvals(args[0].simplify());
            }
            if (node.funcName === 'eigenvects') {
                if (args.length !== 1) throw new Error("eigenvects requires 1 argument");
                return this._eigenvects(args[0].simplify());
            }
            if (node.funcName === 'gramschmidt') {
                if (args.length !== 1) throw new Error("gramschmidt requires 1 argument (list of vectors)");
                return this._gramschmidt(args[0].simplify());
            }
            if (node.funcName === 'lu') {
                if (args.length !== 1) throw new Error("lu requires 1 argument");
                return this._lu(args[0].simplify());
            }
            if (node.funcName === 'qr') {
                if (args.length !== 1) throw new Error("qr requires 1 argument");
                return this._qr(args[0].simplify());
            }
            if (node.funcName === 'cholesky') {
                if (args.length !== 1) throw new Error("cholesky requires 1 argument");
                return this._cholesky(args[0].simplify());
            }
            if (node.funcName === 'desolve') {
                if (args.length < 2) throw new Error("desolve requires at least 2 arguments");
                return this._desolve(args[0], args[1]);
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
size, concat, clear, N,
and, or, not, xor, int, evalf`;

                const latexHelp = `\\text{Available commands: see documentation}`;
                return { type: 'info', text: helpText, toString: () => helpText, toLatex: () => latexHelp };
            }

            return new Call(node.funcName, args);
        }

        if (node instanceof BinaryOp) {
            const left = this._recursiveEval(node.left);
            const right = this._recursiveEval(node.right);
            return new node.constructor(left, right).simplify();
        }

        if (node instanceof Not) {
            const arg = this._recursiveEval(node.arg);
            return new Not(arg).simplify();
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
            return this._sumSymbolic(expr, varNode, start, end);
        }

        let sum = new Num(0);
        for (let i = start.value; i <= end.value; i++) {
            const term = expr.substitute(varNode, new Num(i)).simplify();
            sum = new Add(sum, term).simplify();
        }
        return sum;
    }

    _sumSymbolic(expr, varNode, start, end) {
        expr = expr.expand().simplify();

        // sum(c, k, 1, n) = c*n
        if (!this._dependsOn(expr, varNode)) {
            const count = new Add(new Sub(end, start), new Num(1)).simplify();
            return new Mul(expr, count).simplify();
        }

        // sum(A+B)
        if (expr instanceof Add) {
            return new Add(this._sumSymbolic(expr.left, varNode, start, end), this._sumSymbolic(expr.right, varNode, start, end)).simplify();
        }
        if (expr instanceof Sub) {
            return new Sub(this._sumSymbolic(expr.left, varNode, start, end), this._sumSymbolic(expr.right, varNode, start, end)).simplify();
        }

        // sum(c*A)
        if (expr instanceof Mul) {
            if (!this._dependsOn(expr.left, varNode)) {
                return new Mul(expr.left, this._sumSymbolic(expr.right, varNode, start, end)).simplify();
            }
            if (!this._dependsOn(expr.right, varNode)) {
                return new Mul(expr.right, this._sumSymbolic(expr.left, varNode, start, end)).simplify();
            }
        }

        // Normalize sum(..., k, 1, n)
        // If start != 1, sum(f(k), k, a, b) = sum(f(j+a-1), j, 1, b-a+1)
        if (!(start instanceof Num && start.value === 1)) {
             // Let j = k - a + 1 => k = j + a - 1
             // newEnd = b - a + 1
             const j = new Sym(varNode.name + "_idx"); // avoid collision
             const shift = new Sub(start, new Num(1));
             const subExpr = expr.substitute(varNode, new Add(j, shift)).simplify();
             const newEnd = new Sub(end, shift).simplify();
             // Recurse with 1..newEnd
             const res = this._sumSymbolic(subExpr, j, new Num(1), newEnd);
             // Substitute back? No, result should not depend on j.
             return res;
        }

        const n = end;

        // k
        if (expr instanceof Sym && expr.name === varNode.name) {
            // n(n+1)/2
            return new Div(new Mul(n, new Add(n, new Num(1))), new Num(2)).simplify();
        }

        // k^p
        if (expr instanceof Pow && expr.left instanceof Sym && expr.left.name === varNode.name && expr.right instanceof Num) {
            const p = expr.right.value;
            if (p === 1) return new Div(new Mul(n, new Add(n, new Num(1))), new Num(2)).simplify();
            if (p === 2) {
                // n(n+1)(2n+1)/6
                const term1 = new Mul(n, new Add(n, new Num(1)));
                const term2 = new Add(new Mul(new Num(2), n), new Num(1));
                return new Div(new Mul(term1, term2), new Num(6)).simplify();
            }
            if (p === 3) {
                // (n(n+1)/2)^2
                const base = new Div(new Mul(n, new Add(n, new Num(1))), new Num(2));
                return new Pow(base, new Num(2)).simplify();
            }
        }

        return new Call("sum", [expr, varNode, start, end]);
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
    _resultant(poly1, poly2, varNode) {
        if (!(varNode instanceof Sym)) throw new Error("resultant requires a symbol variable");

        const p1 = this._getPolyCoeffs(poly1, varNode);
        const p2 = this._getPolyCoeffs(poly2, varNode);

        if (!p1 || !p2) throw new Error("Resultant arguments must be polynomials in " + varNode.name);

        const n = p1.maxDeg;
        const m = p2.maxDeg;

        // Resultant of two constants is 1 (no common roots) unless one is 0?
        // If degree is 0, they are constants. Common root impossible unless both are 0 (undefined?).
        // If n=0, Res(a, B) = a^m.
        if (n === 0 && m === 0) return new Num(1);

        const size = n + m;
        if (size === 0) return new Num(1);

        // Construct Sylvester Matrix
        const rows = [];
        for(let i=0; i<m; i++) {
            const row = [];
            // Shifted p1 coeffs
            for(let j=0; j<size; j++) {
                // We want coeff of x^(n - (j-i))? No.
                // Standard: Row i has coeffs a_n, a_{n-1}, ..., a_0 padded.
                // Position of a_n starts at column i.
                // Col j corresponds to power x^(size - 1 - j).
                // Wait, simpler:
                // Row i (0..m-1) has p1 coeffs shifted right by i.
                // Coeffs are usually ordered high to low degree in Sylvester.

                const power = size - 1 - j;
                const p1_power = power - (m - 1 - i); // This is confusing.

                // Let's use standard layout:
                // Row i (for p1, i from 0 to m-1): [0...0, a_n, ..., a_0, 0...0]
                // a_n is at column i.
                const offset = i;
                const deg = n - (j - offset);
                if (deg >= 0 && deg <= n) {
                    row.push(p1.coeffs[deg] || new Num(0));
                } else {
                    row.push(new Num(0));
                }
            }
            rows.push(new Vec(row));
        }

        for(let i=0; i<n; i++) {
            const row = [];
            // Row i (for p2, i from 0 to n-1): [0...0, b_m, ..., b_0, 0...0]
            // b_m is at column i.
            const offset = i;

            for(let j=0; j<size; j++) {
                const deg = m - (j - offset);
                if (deg >= 0 && deg <= m) {
                    row.push(p2.coeffs[deg] || new Num(0));
                } else {
                    row.push(new Num(0));
                }
            }
            rows.push(new Vec(row));
        }

        const matrix = new Vec(rows);
        return this._det(matrix);
    }

    _discriminant(poly, varNode) {
        if (!(varNode instanceof Sym)) throw new Error("discriminant requires a symbol variable");

        const p = this._getPolyCoeffs(poly, varNode);
        if (!p) throw new Error("Discriminant argument must be a polynomial");

        const n = p.maxDeg;
        if (n < 2) return new Num(0); // Linear or constant has disc = 1? Or undefined. Usually 1 or 0? 0 implies multiple root? Linear has no multiple root. 1?
        // Standard definition usually requires deg >= 2.

        const an = p.coeffs[n];

        const deriv = poly.diff(varNode).simplify();
        const res = this._resultant(poly, deriv, varNode);

        // Disc(P) = (-1)^(n(n-1)/2) * 1/an * Res(P, P')

        const signExp = (n * (n - 1)) / 2;
        const sign = (signExp % 2 === 0) ? new Num(1) : new Num(-1);

        return new Div(new Mul(sign, res), an).simplify();
    }

    _solveSystem(eqs, vars) {
        const n = eqs.elements.length;
        const m = vars.elements.length;
        if (n !== m) throw new Error("Number of equations must match number of variables");

        // Convert equations to linear form: A*x = B
        // We build the augmented matrix [A | B]
        const matrixRows = [];

        for (let i = 0; i < n; i++) {
            let eq = eqs.elements[i];
            // Normalize to expression = 0
            if (eq instanceof Eq) {
                eq = new Sub(eq.left, eq.right).simplify();
            }
            eq = eq.expand().simplify();

            const row = [];
            // Extract coefficients for each variable
            for (let j = 0; j < m; j++) {
                const v = vars.elements[j];
                // coeff of v in eq
                // We can use derivative as a heuristic for linear coefficient
                const coeff = eq.diff(v).simplify();
                // Check if coeff depends on other variables (non-linear)
                for (let k = 0; k < m; k++) {
                    if (this._dependsOn(coeff, vars.elements[k])) {
                        return new Call("solve", [eqs, vars]); // Non-linear system
                    }
                }
                row.push(coeff);
            }

            // Constant term (RHS in Ax=B, but we have LHS-RHS=0, so Constant is part of 0)
            // Actually, if we have ax + by + c = 0, then ax + by = -c.
            // So B is -Constant.
            // To get constant term, substitute all vars with 0.
            let constant = eq;
            for (let j = 0; j < m; j++) {
                constant = constant.substitute(vars.elements[j], new Num(0));
            }
            constant = constant.simplify();
            row.push(new Mul(new Num(-1), constant).simplify()); // B

            matrixRows.push(new Vec(row));
        }

        const augmentedMatrix = new Vec(matrixRows);
        const rref = this._rref(augmentedMatrix);

        // Extract solution from RREF
        // If unique solution, RREF should be [I | sol]

        // 1. Check for inconsistency (0 = 1)
        // A row like [0, 0, ..., 0, c] where c != 0
        for (let i = 0; i < n; i++) {
            const row = rref.elements[i];
            let allZero = true;
            for (let j = 0; j < m; j++) {
                const el = row.elements[j];
                if (!(el instanceof Num && el.value === 0)) {
                    allZero = false;
                    break;
                }
            }
            if (allZero) {
                const constant = row.elements[m];
                if (!(constant instanceof Num && constant.value === 0)) {
                    // 0 = c (c!=0) -> No solution
                    return new Vec([]); // Or throw Error? Empty vec usually implies no solution in some systems
                }
            }
        }

        // 2. Extract unique solution
        // Assume pivots are at (i, i). If not, we have infinite solutions or need free variables.
        // For simple unique solution, we expect diagonal 1s.
        const solutions = [];
        for (let i = 0; i < m; i++) { // Loop over variables
            if (i >= n) {
                // More variables than equations?
                // Or dependent rows reduced effective equations.
                // Infinite solutions.
                return new Call("solve", [eqs, vars]);
            }

            const row = rref.elements[i];
            const pivot = row.elements[i]; // Expect 1

            // Check if pivot is 1 (or close to it)
            // And all other coeffs in row are 0?
            // Actually RREF guarantees this for unique solution.

            if (pivot instanceof Num && pivot.value === 1) {
                 solutions.push(row.elements[m]);
            } else {
                 // Not a unique solution (pivot missing or displaced)
                 return new Call("solve", [eqs, vars]);
            }
        }

        return new Vec(solutions);
    }

    _solve(eq, varNode) {
        if (eq instanceof Vec && varNode instanceof Vec) {
            return this._solveSystem(eq, varNode);
        }

        let expr;
        if (eq instanceof Eq) {
            expr = new Sub(eq.left, eq.right).simplify();
        } else {
            expr = eq.simplify();
        }

        // Expand to ensure polynomial form for identification
        expr = expr.expand().simplify();

        try {
            // Try polynomial extraction first
            const poly = this._getPolyCoeffs(expr, varNode);
            if (poly) {
                if (poly.maxDeg === 1) {
                    // ax + b = 0
                    const a = poly.coeffs[1] || new Num(0);
                    const b = poly.coeffs[0] || new Num(0);
                    if (!(a instanceof Num && a.value === 0)) {
                        return new Div(new Mul(new Num(-1), b), a).simplify();
                    }
                } else if (poly.maxDeg === 2) {
                    // ax^2 + bx + c = 0
                    const A = poly.coeffs[2];
                    const B = poly.coeffs[1] || new Num(0);
                    const C = poly.coeffs[0] || new Num(0);

                    if (!(A instanceof Num && A.value === 0)) {
                        const discriminant = new Sub(new Pow(B, new Num(2)), new Mul(new Num(4), new Mul(A, C))).simplify();
                        let sqrtDisc;

                        if (discriminant instanceof Num && discriminant.value < 0) {
                             // Complex roots: i * sqrt(|D|)
                             const absDisc = new Num(-discriminant.value);
                             sqrtDisc = new Mul(new Sym('i'), new Call('sqrt', [absDisc])).simplify();
                        } else {
                             sqrtDisc = new Call("sqrt", [discriminant]);
                        }

                        const sol1 = new Div(new Add(new Mul(new Num(-1), B), sqrtDisc), new Mul(new Num(2), A));
                        const sol2 = new Div(new Sub(new Mul(new Num(-1), B), sqrtDisc), new Mul(new Num(2), A));

                        return new Call("set", [sol1.simplify(), sol2.simplify()]);
                    }
                }
            }

            // Attempt linear solve a*x + b = 0 via differentiation if poly extraction failed (e.g. non-polynomial terms that simplify)
            // a = derivative w.r.t varNode
            const a = expr.diff(varNode).simplify();
            const b = expr.substitute(varNode, new Num(0)).simplify();

            const reconstruction = new Add(new Mul(a, varNode), b).simplify();
            const diff = new Sub(expr, reconstruction).expand().simplify();

            if (diff instanceof Num && diff.value === 0 && !(a instanceof Num && a.value === 0)) {
                return new Div(new Mul(new Num(-1), b), a).simplify();
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

    _integrateSubstitution(expr, varNode) {
        // Handle trig rewrites for substitution
        if (expr instanceof Call) {
            if (expr.funcName === 'tan') {
                // sin/cos
                expr = new Div(new Call('sin', expr.args), new Call('cos', expr.args));
            } else if (expr.funcName === 'cot') {
                // cos/sin
                expr = new Div(new Call('cos', expr.args), new Call('sin', expr.args));
            } else if (expr.funcName === 'tanh') {
                expr = new Div(new Call('sinh', expr.args), new Call('cosh', expr.args));
            }
        }

        // Look for pattern f'(x)/f(x) -> ln(f(x))

        // 1. Check if expr is Div
        if (expr instanceof Div) {
            const num = expr.left;
            const den = expr.right;

            // Check derivative of denominator
            const diffDen = den.diff(varNode).simplify();

            // Check if num = k * diffDen
            // k = num / diffDen should be constant (no varNode)
            try {
                const ratio = new Div(num, diffDen).simplify();
                if (!this._dependsOn(ratio, varNode)) {
                    // Result: ratio * ln(den)
                    return new Mul(ratio, new Call("ln", [den])).simplify();
                }
            } catch(e) {}
        }

        return new Call("integrate", [expr, varNode]);
    }

    _integrateByParts(expr, varNode) {
        if (!(expr instanceof Mul)) return new Call("integrate", [expr, varNode]);

        const factors = [];
        const flattenMul = (node) => {
            if (node instanceof Mul) {
                flattenMul(node.left);
                flattenMul(node.right);
            } else {
                factors.push(node);
            }
        };
        flattenMul(expr);

        // LIATE Rule for choosing u
        // Log, Inverse, Alg, Trig, Exp
        const typeScore = (node) => {
            if (node instanceof Call) {
                if (['ln', 'log', 'log2'].includes(node.funcName)) return 5; // Log
                if (['asin', 'acos', 'atan', 'asec', 'acsc', 'acot'].includes(node.funcName)) return 4; // Inverse
                if (['sin', 'cos', 'tan', 'sec', 'csc', 'cot'].includes(node.funcName)) return 2; // Trig
                if (['exp'].includes(node.funcName)) return 1; // Exp
            }
            if (node instanceof Pow) {
                 if (this._dependsOn(node.right, varNode)) return 1; // Exp-like? x^x
                 // x^n -> Alg
                 return 3;
            }
            if (node instanceof Sym && node.name === varNode.name) return 3; // Alg (x)
            return 0; // Constant or unknown (treated as scalar part of u or dv?)
        };

        // Create candidates for u
        const candidates = factors.map((f, i) => ({
            idx: i,
            node: f,
            score: typeScore(f)
        })).filter(c => this._dependsOn(c.node, varNode)); // Only variable parts matter for u choice

        // Sort by score descending (Log first)
        candidates.sort((a, b) => b.score - a.score);

        for(const cand of candidates) {
            const u = cand.node;

            // dv is the product of all other factors
            let dv = null;
            for(let j=0; j<factors.length; j++) {
                if (cand.idx === j) continue;
                if (dv === null) dv = factors[j];
                else dv = new Mul(dv, factors[j]);
            }

            if (dv === null) continue;

            // Check if dv is easily integrable
            const v = dv.integrate(varNode);
            if (v instanceof Call && v.funcName === 'integrate') continue; // dv not integrable

            const du = u.diff(varNode).simplify();

            const uv = new Mul(u, v).simplify();
            const vdu = new Mul(v, du).simplify();

            // Check if vdu is simpler than original expr?
            // To prevent cycles, we could check string length or something, but LIATE usually guarantees progress.
            // One edge case: e^x sin(x). Both are Trig/Exp.
            // Score: Trig(2) > Exp(1). u=sin, dv=exp.
            // v=exp, du=cos. vdu = exp*cos.
            // exp*cos -> u=cos, dv=exp. v=exp, du=-sin. vdu = -exp*sin.
            // Cycle.
            // This implementation performs ONE step and returns "uv - integrate(vdu)".
            // The outer 'evaluate' loop will call 'integrate' on 'vdu'.
            // If vdu is integrable (by parts again), it proceeds.
            // For cyclic cases, it will expand indefinitely until max depth or stack overflow unless we catch it.
            // But for standard x*e^x, it terminates.

            return new Sub(uv, new Call("integrate", [vdu, varNode])).simplify();
        }

        return new Call("integrate", [expr, varNode]);
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

            // Find pivot (non-zero)
            const isZero = (val) => {
                 val = val.simplify();
                 if (val instanceof Num && val.value === 0) return true;
                 return false; // Symbolic or non-zero number is treated as non-zero
            };

            while (isZero(M[i][lead])) {
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

        // Polynomial GCD
        try {
             // Identify variable
             let varNode = null;
             const findVar = (n) => {
                 if (n instanceof Sym && !['pi', 'e', 'i', 'infinity'].includes(n.name)) {
                     if (!varNode) varNode = n;
                     else if (varNode.name !== n.name) varNode = 'MIXED';
                 }
                 if (n instanceof Call) n.args.forEach(findVar);
                 if (n instanceof BinaryOp) { findVar(n.left); findVar(n.right); }
                 if (n instanceof Pow) { findVar(n.left); findVar(n.right); }
             };
             findVar(a);
             findVar(b);

             if (varNode && varNode !== 'MIXED') {
                 // Check if both are polynomials
                 const p1 = this._getPolyCoeffs(a, varNode);
                 const p2 = this._getPolyCoeffs(b, varNode);

                 if (p1 && p2) {
                     // Euclidean Algorithm for Polynomials
                     let u = a;
                     let v = b;

                     // Make monic? Or just standard division.
                     // We need _polyRem(u, v, varNode)

                     let loopCount = 0;
                     while (!(v instanceof Num && v.value === 0)) {
                         if (loopCount++ > 50) break; // Safety break
                         const r = this._polyRem(u, v, varNode);

                         if (r instanceof Call && r.funcName === 'rem') {
                             return new Call('gcd', [a, b]);
                         }

                         // Check if degree reduced
                         const pV = this._getPolyCoeffs(v, varNode);
                         const pR = this._getPolyCoeffs(r, varNode);

                         if (!pR || !pV) return new Call('gcd', [a, b]);

                         if (pR.maxDeg >= pV.maxDeg) {
                             // Degree did not decrease, symbolic cancellation failed
                             if (!(r instanceof Num && r.value === 0)) {
                                 // Not zero, and degree not reduced. Abort.
                                 return new Call('gcd', [a, b]);
                             }
                         }

                         u = v;
                         v = r;
                     }

                     // Normalize result (make leading coefficient 1)
                     const finalP = this._getPolyCoeffs(u, varNode);
                     if (finalP) {
                         const lc = finalP.coeffs[finalP.maxDeg];
                         if (lc && !(lc instanceof Num && lc.value === 0)) {
                             return new Div(u, lc).simplify();
                         }
                     }
                     return u.simplify();
                 }
             }
        } catch(e) {
            // console.log("GCD error", e);
        }

        return new Call('gcd', [a, b]);
    }

    _polyRem(a, b, varNode) {
        // Polynomial Remainder of A / B
        // A = Q*B + R
        // We use polynomial division logic
        const pA = this._getPolyCoeffs(a, varNode);
        const pB = this._getPolyCoeffs(b, varNode);

        if (!pA || !pB) return new Call('rem', [a, b]); // Should not happen if called from GCD

        let degA = pA.maxDeg;
        const degB = pB.maxDeg;

        if (degB < 0) throw new Error("Division by zero polynomial");
        if (degA < degB) return a; // Degree A < Degree B

        // Clone coeffs of A to work on (remainder starts as A)
        let R = a;
        let pR = pA;

        // While deg(R) >= deg(B)
        // term = LT(R) / LT(B)
        // R = R - term * B

        // Safety counter
        let limit = degA - degB + 5;
        while (degA >= degB && limit-- > 0) {
             const ltR = pR.coeffs[degA];
             const ltB = pB.coeffs[degB];

             // term = (lc(R)/lc(B)) * x^(degR - degB)
             const coeff = new Div(ltR, ltB).simplify();
             const degDiff = degA - degB;

             let termVar;
             if (degDiff === 0) termVar = new Num(1);
             else if (degDiff === 1) termVar = varNode;
             else termVar = new Pow(varNode, new Num(degDiff));

             const term = new Mul(coeff, termVar).simplify();

             // R = R - term * B
             R = new Sub(R, new Mul(term, b)).simplify();

             // Re-evaluate coeffs of R
             pR = this._getPolyCoeffs(R, varNode);
             if (!pR) break; // Should be 0
             // Check if R is 0
             if (Object.keys(pR.coeffs).length === 0) {
                 return new Num(0);
             }

             // Update degA
             // Usually maxDeg decreases
             degA = pR.maxDeg;
             // Check if maxDeg actually dropped?
             // Floating point issues might keep small coeffs.
             // _getPolyCoeffs does not filter small coeffs automatically, but simplify() usually handles it.

             // Check if maxDeg has NOT dropped and we made no progress
             if (degA >= degDiff + degB) {
                 break;
             }
        }

        return R;
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

            return new Call('factored', exprs); // Use a custom wrapper to display it
        }

        // Polynomial Factorization
        try {
            // 1. Identify variable
            const vars = [];
            const findVars = (node) => {
                if (node instanceof Sym && !['pi', 'e', 'i', 'infinity'].includes(node.name)) {
                    if (!vars.some(v => v.name === node.name)) vars.push(node);
                }
                if (node instanceof Call) node.args.forEach(findVars);
                if (node instanceof BinaryOp) { findVars(node.left); findVars(node.right); }
                if (node instanceof Pow) { findVars(node.left); findVars(node.right); }
            };
            findVars(n);

            if (vars.length === 1) {
                const x = vars[0];
                const poly = this._getPolyCoeffs(n, x);

                if (poly) {
                    let factors = [];
                    let p = poly; // { coeffs: {0: Num, 1: Num...}, maxDeg: N }

                    // 2. Extract Common Content (Integer GCD of coefficients)
                    const coefList = [];
                    for(const d in p.coeffs) {
                        const c = p.coeffs[d];
                        if (c instanceof Num && Number.isInteger(c.value)) {
                            coefList.push(Math.abs(c.value));
                        } else {
                            // Non-integer coefficient, abort content extraction
                            coefList.length = 0;
                            break;
                        }
                    }

                    let content = 1;
                    if (coefList.length > 0) {
                         const calcGcd = (u, v) => !v ? u : calcGcd(v, u % v);
                         content = coefList.reduce((acc, val) => calcGcd(acc, val));
                    }

                    if (content > 1) {
                         factors.push(new Num(content));
                         // Divide poly by content
                         const newCoeffs = {};
                         for(const d in p.coeffs) {
                             newCoeffs[d] = new Div(p.coeffs[d], new Num(content)).simplify();
                         }
                         p.coeffs = newCoeffs;
                    }

                    // 3. Extract lowest degree term (x^k)
                    let minDeg = p.maxDeg;
                    let hasTerms = false;
                    for(const d in p.coeffs) {
                        const deg = parseInt(d);
                        const c = p.coeffs[d];
                        if (!(c instanceof Num && c.value === 0)) {
                             if (deg < minDeg) minDeg = deg;
                             hasTerms = true;
                        }
                    }
                    if (!hasTerms) return new Num(0);

                    if (minDeg > 0) {
                        if (minDeg === 1) factors.push(x);
                        else factors.push(new Pow(x, new Num(minDeg)));

                        // Shift degrees down
                        const newCoeffs = {};
                        let newMaxDeg = 0;
                        for(const d in p.coeffs) {
                             const deg = parseInt(d);
                             const newDeg = deg - minDeg;
                             newCoeffs[newDeg] = p.coeffs[d];
                             if (newDeg > newMaxDeg) newMaxDeg = newDeg;
                        }
                        p.coeffs = newCoeffs;
                        p.maxDeg = newMaxDeg;
                    }

                    // 4. Higher Order Factorization (Rational Root Theorem)
                    if (p.maxDeg > 2) {
                         const ratRes = this._factorRational(p.coeffs, p.maxDeg, x);
                         if (ratRes.found) {
                             factors.push(...ratRes.factors);
                             p.coeffs = ratRes.newCoeffs;
                             p.maxDeg = ratRes.newMaxDeg;
                         }
                    }

                    // Reconstruct reduced polynomial to solve/analyze (potentially reduced by step 4)
                    let reducedPoly = this._poly2symb(this._coeffsToVec(p.coeffs, p.maxDeg), x);

                    // 5. Factor Quadratics or remainders
                    if (p.maxDeg === 2) {
                        const a = p.coeffs[2];
                        const rootsRes = this._solve(reducedPoly, x);
                        let roots = [];
                        if (rootsRes instanceof Call && rootsRes.funcName === 'set') {
                            roots = rootsRes.args;
                        } else if (rootsRes instanceof Expr && !(rootsRes instanceof Call && rootsRes.funcName === 'solve')) {
                            roots = [rootsRes];
                        }

                        if (roots.length > 0) {
                            factors.push(a); // leading coeff
                            for(const r of roots) {
                                factors.push(new Sub(x, r).simplify());
                            }
                        } else {
                            factors.push(reducedPoly);
                        }
                    } else if (p.maxDeg === 1) {
                        factors.push(reducedPoly);
                    } else if (p.maxDeg > 2) {
                        // Still > 2 implies no rational roots found, or irreducible cubic+
                        factors.push(reducedPoly);
                    } else if (p.maxDeg === 0) {
                         if (!(reducedPoly instanceof Num && reducedPoly.value === 1)) {
                            factors.push(reducedPoly);
                         }
                    }

                    // 6. Combine factors
                    if (factors.length === 0) return new Num(1);
                    let result = factors[0];
                    for(let i=1; i<factors.length; i++) {
                        result = new Mul(result, factors[i]); // Don't simplify fully to avoid expanding back?
                        // If I use Mul constructor, it simplifies by default in this system?
                        // The .simplify() method is called inside evaluate().
                        // If I return a raw tree, evaluate() calls simplify().
                        // Mul.simplify() might distribute.
                        // We need a 'factored' wrapper or ensure Mul doesn't expand.
                        // In this code, simplify() does expand sometimes.
                        // Let's return Call('factored', factors) if we want to preserve structure?
                        // Or just rely on Mul not expanding unless .expand() is called.
                        // Mul.simplify usually just merges numbers and canonicalizes. It does NOT expand (a+b)(c+d).
                    }
                    return result; // returning Mul tree.
                }
            }
        } catch (e) {
            // console.log(e);
        }

        return new Call('factor', [n]);
    }

    _coeffsToVec(coeffs, maxDeg) {
        const list = [];
        for(let i=maxDeg; i>=0; i--) {
            list.push(coeffs[i] || new Num(0));
        }
        return new Vec(list);
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

    _nextprime(n) {
        n = n.simplify();
        if (n instanceof Num && Number.isInteger(n.value)) {
            let val = n.value + 1;
            while (true) {
                if (this._isPrime(new Num(val)).value === 1) return new Num(val);
                val++;
            }
        }
        return new Call('nextprime', [n]);
    }

    _prevprime(n) {
        n = n.simplify();
        if (n instanceof Num && Number.isInteger(n.value)) {
            let val = n.value - 1;
            while (val >= 2) {
                if (this._isPrime(new Num(val)).value === 1) return new Num(val);
                val--;
            }
            return new Call('prevprime', [n]); // No prime found? Or undefined?
        }
        return new Call('prevprime', [n]);
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

    _append(list, elem) {
        if (list instanceof Vec) {
             const newElements = [...list.elements, elem];
             return new Vec(newElements);
        }
        return new Call('append', [list, elem]);
    }

    _prepend(list, elem) {
        if (list instanceof Vec) {
             const newElements = [elem, ...list.elements];
             return new Vec(newElements);
        }
        return new Call('prepend', [list, elem]);
    }

    _arg(z) {
        z = z.simplify();

        if (z instanceof Num) {
             if (z.value >= 0) return new Num(0);
             return new Num(Math.PI);
        }
        return new Call('arg', [z]);
    }

    _kernel(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("kernel requires a matrix");
        const R = this._rref(matrix);
        const rows = R.elements.length;
        if (rows === 0) return new Vec([]);
        const cols = R.elements[0].elements.length;

        const pivots = [];
        const pivotRows = [];

        for(let r=0; r<rows; r++) {
            let foundPivot = false;
            for(let c=0; c<cols; c++) {
                const val = R.elements[r].elements[c].evaluateNumeric();
                if (!foundPivot && Math.abs(val) > 1e-10) {
                    pivots.push(c);
                    pivotRows.push(r);
                    foundPivot = true;
                }
            }
        }

        const free = [];
        for(let c=0; c<cols; c++) {
            if (!pivots.includes(c)) free.push(c);
        }

        if (free.length === 0) {
            const zeroVec = [];
            for(let i=0; i<cols; i++) zeroVec.push(new Num(0));
            return new Vec([new Vec(zeroVec)]);
        }

        const basis = [];
        for(const freeIdx of free) {
            const vec = new Array(cols).fill(null).map(() => new Num(0));
            vec[freeIdx] = new Num(1);

            for(let i=0; i<pivots.length; i++) {
                const p = pivots[i];
                const r = pivotRows[i];
                const coeff = R.elements[r].elements[freeIdx];
                vec[p] = new Mul(new Num(-1), coeff).simplify();
            }
            basis.push(new Vec(vec));
        }

        return new Vec(basis);
    }

    _basis(matrix) {
         if (!(matrix instanceof Vec)) throw new Error("basis requires a matrix");
         const R = this._rref(matrix);
         const rows = R.elements.length;
         if (rows === 0) return new Vec([]);
         const cols = R.elements[0].elements.length;

         const pivots = [];
         for(let r=0; r<rows; r++) {
            for(let c=0; c<cols; c++) {
                const val = R.elements[r].elements[c].evaluateNumeric();
                if (Math.abs(val) > 1e-10) {
                    pivots.push(c);
                    break;
                }
            }
         }

         // Return columns of original matrix
         const basisCols = [];
         // Transpose first to get rows (which are columns of original)
         const T = this._trans(matrix);
         for(const p of pivots) {
             basisCols.push(T.elements[p]);
         }
         return new Vec(basisCols);
    }

    _eigenvals(matrix) {
        // solve(charpoly(M, x), x)
        const lambda = new Sym('lambda_');
        const cp = this._charpoly(matrix, lambda);
        const sols = this._solve(cp, lambda);

        if (sols instanceof Call && sols.funcName === 'set') {
            return new Vec(sols.args);
        }
        return new Vec([sols]);
    }

    _eigenvects(matrix) {
        const evals = this._eigenvals(matrix);
        // evals is a Vec of values
        let vectors = [];

        // Use a Set to avoid duplicates if eigenvalues are repeated?
        // _solve might return multiple roots.

        for(const val of evals.elements) {
            // kernel(A - val*I)
            const rows = matrix.elements.length;
            const cols = matrix.elements[0].elements.length;

            // Construct A - val*I
            const M_minus_lambdaI = [];
            for(let i=0; i<rows; i++) {
                const row = [];
                for(let j=0; j<cols; j++) {
                    let el = matrix.elements[i].elements[j];
                    if (i === j) {
                        el = new Sub(el, val).simplify();
                    }
                    row.push(el);
                }
                M_minus_lambdaI.push(new Vec(row));
            }
            const mat = new Vec(M_minus_lambdaI);

            const kern = this._kernel(mat);
            // kern is Vec of basis vectors
            for(const v of kern.elements) {
                vectors.push(v);
            }
        }
        return new Vec(vectors);
    }

    _gramschmidt(vectors) {
        // vectors: Vec of Vecs (list of vectors)
        if (!(vectors instanceof Vec)) throw new Error("gramschmidt requires a list of vectors");

        const basis = [];
        for(const v of vectors.elements) {
            let u = v;
            for(const b of basis) {
                const dotVB = new Call('dot', [v, b]).simplify(); // v . b
                const dotBB = new Call('dot', [b, b]).simplify();
                const proj = new Mul(new Div(dotVB, dotBB), b).simplify();
                u = new Sub(u, proj).simplify();
            }

            const isZero = u.elements.every(e => e instanceof Num && e.value === 0);
            if (!isZero) {
                 const n = new Call('norm', [u]).simplify();
                 const e = new Div(u, n).simplify();
                 basis.push(e);
            }
        }
        return new Vec(basis);
    }

    _lu(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("lu requires a matrix");
        const rows = matrix.elements.length;
        if (rows === 0) return new Vec([]);
        const cols = matrix.elements[0].elements.length;
        if (rows !== cols) throw new Error("lu requires a square matrix");

        // Initialize P, L, U
        // U starts as copy of A
        const U_rows = matrix.elements.map(row => [...row.elements]);
        const L_rows = [];
        const P_rows = [];

        for (let i = 0; i < rows; i++) {
            const lRow = [];
            const pRow = [];
            for (let j = 0; j < rows; j++) {
                lRow.push(new Num(i === j ? 1 : 0));
                pRow.push(new Num(i === j ? 1 : 0));
            }
            L_rows.push(lRow);
            P_rows.push(pRow);
        }

        for (let k = 0; k < rows; k++) {
            // Pivot strategy: Find first non-zero element in column k from k to rows-1
            let pivotIdx = k;
            let foundPivot = false;

            // Check current diagonal first
            let val = U_rows[k][k].evaluateNumeric();
            if (isNaN(val) || Math.abs(val) > 1e-10) {
                pivotIdx = k;
                foundPivot = true;
            } else {
                // Search down
                for (let i = k + 1; i < rows; i++) {
                    val = U_rows[i][k].evaluateNumeric();
                    if (isNaN(val) || Math.abs(val) > 1e-10) {
                        pivotIdx = i;
                        foundPivot = true;
                        break;
                    }
                }
            }

            if (!foundPivot) continue;

            // Swap rows in U, P
            // Swap columns 0..k-1 in L (vectors before diagonal)
            if (pivotIdx !== k) {
                [U_rows[k], U_rows[pivotIdx]] = [U_rows[pivotIdx], U_rows[k]];
                [P_rows[k], P_rows[pivotIdx]] = [P_rows[pivotIdx], P_rows[k]];

                for (let j = 0; j < k; j++) {
                    [L_rows[k][j], L_rows[pivotIdx][j]] = [L_rows[pivotIdx][j], L_rows[k][j]];
                }
            }

            const pivot = U_rows[k][k];

            for (let i = k + 1; i < rows; i++) {
                // factor = U[i][k] / U[k][k]
                const factor = new Div(U_rows[i][k], pivot).simplify();
                L_rows[i][k] = factor;
                U_rows[i][k] = new Num(0); // Explicitly zero out

                for (let j = k + 1; j < rows; j++) {
                    // U[i][j] = U[i][j] - factor * U[k][j]
                    U_rows[i][j] = new Sub(U_rows[i][j], new Mul(factor, U_rows[k][j])).simplify();
                }
            }
        }

        const L = new Vec(L_rows.map(r => new Vec(r)));
        const U = new Vec(U_rows.map(r => new Vec(r)));
        const P = new Vec(P_rows.map(r => new Vec(r)));

        return new Vec([L, U, P]);
    }

    _qr(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("qr requires a matrix");
        const cols = [];
        const rows = matrix.elements.length;
        const numCols = matrix.elements[0].elements.length;

        for(let j=0; j<numCols; j++) {
            const col = [];
            for(let i=0; i<rows; i++) {
                col.push(matrix.elements[i].elements[j]);
            }
            cols.push(new Vec(col));
        }

        const Q_cols = this._gramschmidt(new Vec(cols));

        const Q_rows = [];
        for(let i=0; i<rows; i++) {
             const row = [];
             for(let j=0; j<Q_cols.elements.length; j++) {
                 row.push(Q_cols.elements[j].elements[i]);
             }
             Q_rows.push(new Vec(row));
        }
        const Q = new Vec(Q_rows);

        const R = new Mul(this._trans(Q), matrix).simplify();

        return new Vec([Q, R]);
    }

    _cholesky(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("cholesky requires a matrix");
        const rows = matrix.elements.length;
        if (rows === 0) return new Vec([]);
        const cols = matrix.elements[0].elements.length;
        if (rows !== cols) throw new Error("cholesky requires a square matrix");

        const L = [];
        for (let i = 0; i < rows; i++) {
            L.push(new Array(rows).fill(new Num(0)));
        }

        for (let i = 0; i < rows; i++) {
            for (let j = 0; j <= i; j++) {
                let sum = new Num(0);
                for (let k = 0; k < j; k++) {
                    sum = new Add(sum, new Mul(L[i][k], L[j][k])).simplify();
                }

                if (i === j) {
                    // Diagonal
                    // L[i][i] = sqrt( A[i][i] - sum )
                    const diff = new Sub(matrix.elements[i].elements[i], sum).simplify();
                    L[i][j] = new Call('sqrt', [diff]).simplify();
                } else {
                    // L[i][j] = (1/L[j][j]) * (A[i][j] - sum)
                    const diff = new Sub(matrix.elements[i].elements[j], sum).simplify();
                    L[i][j] = new Div(diff, L[j][j]).simplify();
                }
            }
        }

        return new Vec(L.map(r => new Vec(r)));
    }

    _desolve(eq, depVar) {
        // Handle Systems: eq is Vec, depVar is Vec
        if (eq instanceof Vec && depVar instanceof Vec) {
            return this._desolveSystem(eq, depVar);
        }

        // Handle depVar being y or y(x)
        if (depVar instanceof Call) {
            // Assume format y(x)
            depVar = new Sym(depVar.funcName);
        }
        if (!(depVar instanceof Sym)) throw new Error("Second argument to desolve must be a dependent variable");

        // Normalize equation: LHS - RHS
        let expr;
        if (eq instanceof Eq) {
            expr = new Sub(eq.left, eq.right).simplify();
        } else {
            expr = eq.simplify();
        }

        // Find independent variable
        // It's the variable in diff(y, x) or just implied if we find x in expr
        // Scan for diff(depVar, x)
        let indepVar = null;
        let diffCall = null;

        const findDiff = (node) => {
            if (node instanceof Call && node.funcName === 'diff') {
                // Ensure order is 1 (args.length == 2 or args[2] == 1)
                const validOrder = (node.args.length === 2) || (node.args.length === 3 && node.args[2] instanceof Num && node.args[2].value === 1);
                if (!validOrder) return;

                // Case 1: diff(y, x) where y is depVar
                if (node.args[0] instanceof Sym && node.args[0].name === depVar.name) {
                    if (node.args.length > 1 && node.args[1] instanceof Sym) {
                        indepVar = node.args[1];
                        diffCall = node;
                        return;
                    }
                }
            }
            if (node instanceof Call) node.args.forEach(findDiff);
            if (node instanceof BinaryOp) { findDiff(node.left); findDiff(node.right); }
        };

        // Re-scan for diff(y(x), x)
        const findDiffComplex = (node) => {
             if (node instanceof Call && node.funcName === 'diff') {
                 // Ensure order is 1
                 const validOrder = (node.args.length === 2) || (node.args.length === 3 && node.args[2] instanceof Num && node.args[2].value === 1);
                 if (!validOrder) return;

                 const target = node.args[0];
                 // target is y(x) -> Call with funcName = depVar.name
                 if (target instanceof Call && target.funcName === depVar.name) {
                      if (node.args.length > 1 && node.args[1] instanceof Sym) {
                          indepVar = node.args[1];
                          diffCall = node;
                      }
                 }
             }
             if (node instanceof Call) node.args.forEach(findDiffComplex);
             if (node instanceof BinaryOp) { findDiffComplex(node.left); findDiffComplex(node.right); }
        };

        // Attempt to find diff in evaluated expression
        findDiff(expr);
        if(!indepVar) findDiffComplex(expr);

        // Check for Second Derivative diff(y, x, 2)
        let diffCall2 = null;
        const findDiff2 = (node) => {
            if (node instanceof Call && node.funcName === 'diff') {
                if (node.args.length === 3 && node.args[2] instanceof Num && node.args[2].value === 2) {
                    const target = node.args[0];
                    if ((target instanceof Sym && target.name === depVar.name) ||
                        (target instanceof Call && target.funcName === depVar.name)) {
                            // Check indep var
                            if (node.args[1] instanceof Sym) {
                                if (indepVar && node.args[1].name !== indepVar.name) return; // Mismatch
                                indepVar = node.args[1];
                                diffCall2 = node;
                            }
                    }
                }
            }
            if (node instanceof Call) node.args.forEach(findDiff2);
            if (node instanceof BinaryOp) { findDiff2(node.left); findDiff2(node.right); }
        };
        findDiff2(expr);

        // Handle Second Order
        if (diffCall2) {
             const D2Y = new Sym('__D2Y__');
             const DY = new Sym('__DY__');

             // Robust Linearity Check for Homogeneous Constant Coefficient
             // Eq must be: a*y'' + b*y' + c*y = 0 (assuming LHS-RHS simplification)
             // Check if it fits this form where a, b, c are constants (no x, no y)

             // 1. Substitute derivatives with symbols to treat as polynomial in variables
             const y_sym = new Sym(depVar.name);

             const replaceDerivs = (node) => {
                 if (node === diffCall2) return D2Y;
                 // diff(y, x, 2)
                 if (node instanceof Call && node.funcName === 'diff' && node.args.length === 3 && node.args[2].value === 2) {
                     const target = node.args[0];
                     if ((target instanceof Sym && target.name === depVar.name) || (target instanceof Call && target.funcName === depVar.name)) return D2Y;
                 }
                 // diff(y, x)
                 if (node instanceof Call && node.funcName === 'diff' && node.args.length === 2) {
                     const target = node.args[0];
                     if ((target instanceof Sym && target.name === depVar.name) || (target instanceof Call && target.funcName === depVar.name)) return DY;
                 }
                 // y
                 if (node instanceof Call && node.funcName === depVar.name) return y_sym;
                 if (node instanceof Sym && node.name === depVar.name) return y_sym;

                 if (node instanceof BinaryOp) return new node.constructor(replaceDerivs(node.left), replaceDerivs(node.right));
                 if (node instanceof Call) return new Call(node.funcName, node.args.map(replaceDerivs));
                 if (node instanceof Not) return new Not(replaceDerivs(node.arg));
                 return node;
             };

             const polyForm = replaceDerivs(expr).expand().simplify();

             // Extract Coefficients
             // Coeff of D2Y
             const a = polyForm.diff(D2Y).simplify();
             // Coeff of DY
             const b = polyForm.diff(DY).simplify();
             // Coeff of y
             const c = polyForm.diff(y_sym).simplify();

             // Verify Linearity:
             // 1. Coefficients a, b, c must not depend on y, DY, D2Y
             const dependsOnVars = (n) => this._dependsOn(n, y_sym) || this._dependsOn(n, DY) || this._dependsOn(n, D2Y);

             if (!dependsOnVars(a) && !dependsOnVars(b) && !dependsOnVars(c)) {
                 // 2. Reconstruction check: expr should equal a*D2Y + b*DY + c*y + remainder
                 // Remainder must be 0 for Homogeneous
                 const recon = new Add(new Add(new Mul(a, D2Y), new Mul(b, DY)), new Mul(c, y_sym)).simplify();
                 const remainder = new Sub(polyForm, recon).simplify();

                 // Check if remainder is 0 (or close to 0)
                 const isHomogeneous = (remainder instanceof Num && Math.abs(remainder.value) < 1e-9);

                 // Check for Constant Coefficients (independent of x)
                 if (isHomogeneous && !this._dependsOn(a, indepVar) && !this._dependsOn(b, indepVar) && !this._dependsOn(c, indepVar)) {
                     // Solve Characteristic Equation: a*r^2 + b*r + c = 0
                     const r = new Sym('r');
                     const charPoly = new Add(new Add(new Mul(a, new Pow(r, new Num(2))), new Mul(b, r)), c);

                     let roots = this._solve(charPoly, r);

                     let r1, r2;
                     if (roots instanceof Call && roots.funcName === 'set') {
                         r1 = roots.args[0];
                         r2 = roots.args[1];
                     } else if (roots instanceof Expr && !(roots instanceof Call && roots.funcName === 'solve')) {
                         r1 = roots;
                         r2 = roots;
                     } else {
                         // Failed to solve characteristic equation
                         return new Call('desolve', [eq, depVar]);
                     }

                     const C1 = new Sym('C1');
                     const C2 = new Sym('C2');

                     const diffR = new Sub(r1, r2).simplify();
                     // Check for repeated roots
                     if (diffR instanceof Num && Math.abs(diffR.value) < 1e-9) {
                         const term1 = new Mul(C1, new Call('exp', [new Mul(r1, indepVar)])).simplify();
                         const term2 = new Mul(C2, new Mul(indepVar, new Call('exp', [new Mul(r1, indepVar)]))).simplify();
                         return new Add(term1, term2);
                     } else {
                         const term1 = new Mul(C1, new Call('exp', [new Mul(r1, indepVar)])).simplify();
                         const term2 = new Mul(C2, new Call('exp', [new Mul(r2, indepVar)])).simplify();
                         return new Add(term1, term2);
                     }
                 }
             }

             // Check simple y'' = f(x)
             // Check if expr depends on y or y'
             // If we replace y'' with D2Y, and result only depends on D2Y and x
             const replaceD2Y = (node) => {
                  if (node === diffCall2) return D2Y;
                  if (node instanceof Call && node.funcName === 'diff' && node.args.length === 3 && node.args[2].value === 2) return D2Y;
                  if (node instanceof BinaryOp) return new node.constructor(replaceD2Y(node.left), replaceD2Y(node.right));
                  if (node instanceof Call) return new Call(node.funcName, node.args.map(replaceD2Y));
                  return node;
             };

             const d2yEq = replaceD2Y(expr).simplify();
             // Check if depends on diffCall (diff(y,x)) if diffCall exists
             const dependsOnFirst = diffCall ? this._dependsOn(d2yEq, diffCall) : false;

             if (!this._dependsOn(d2yEq, depVar) && !dependsOnFirst) {
                 // Solve for D2Y
                 const sol = this._solve(d2yEq, D2Y);
                 if (!(sol instanceof Call && sol.funcName === 'solve')) {
                     // Integrate twice
                     const firstInt = sol.integrate(indepVar).simplify();
                     const C1 = new Sym('C1');
                     const firstWithC = new Add(firstInt, C1);
                     const secondInt = firstWithC.integrate(indepVar).simplify();
                     const C2 = new Sym('C2');
                     return new Add(secondInt, C2);
                 }
             }
        }


        if (!indepVar) {
            // Fallback if no derivative found (maybe algebraic eq involving y?)
            return new Call('desolve', [eq, depVar]);
        }

        // Solve algebraically for diffCall (treat diffCall as variable DY)
        const DY = new Sym('__DY__');

        // Helper to replace diffCall with DY deeply
        const replaceDiff = (node) => {
            // Check reference equality
            if (node === diffCall) return DY;

            // Check structural equality for diff(y, x)
            if (node instanceof Call && node.funcName === 'diff' && node.args.length === 2) {
                 const a0 = node.args[0];
                 const a1 = node.args[1];
                 const d0 = diffCall.args[0];
                 const d1 = diffCall.args[1];

                 // Compare a0 and d0
                 let match0 = false;
                 if (a0 instanceof Sym && d0 instanceof Sym && a0.name === d0.name) match0 = true;
                 if (a0 instanceof Call && d0 instanceof Call && a0.toString() === d0.toString()) match0 = true;

                 // Compare a1 and d1
                 let match1 = false;
                 if (a1 instanceof Sym && d1 instanceof Sym && a1.name === d1.name) match1 = true;

                 if (match0 && match1) return DY;
            }

            if (node instanceof BinaryOp) return new node.constructor(replaceDiff(node.left), replaceDiff(node.right));
            if (node instanceof Call) return new Call(node.funcName, node.args.map(replaceDiff));
            if (node instanceof Not) return new Not(replaceDiff(node.arg));
            return node;
        };

        const algEq = replaceDiff(expr);
        let sol = this._solve(algEq, DY); // Solve for y'

        // Fallback: Use polynomial extraction if _solve failed (e.g. simplification issues)
        if (sol instanceof Call && sol.funcName === 'solve') {
             const polyDY = this._getPolyCoeffs(algEq, DY);
             if (polyDY && polyDY.maxDeg === 1) {
                 const a = polyDY.coeffs[1] || new Num(0);
                 const b = polyDY.coeffs[0] || new Num(0);
                 // a*DY + b = 0 => DY = -b/a
                 sol = new Div(new Mul(new Num(-1), b), a).simplify();
             } else {
                 return new Call('desolve', [eq, depVar]);
             }
        }

        const f_xy = sol;

        // 1. Direct Integration: y' = f(x) (f_xy does not depend on depVar)
        if (!this._dependsOn(f_xy, depVar)) {
             // y = int(f(x)) + C
             const integral = f_xy.integrate(indepVar).simplify();
             const C = new Sym('C1');
             return new Add(integral, C);
        }

        // 2. Linear First Order: y' = a(x) + b(x)y  => y' - b(x)y = a(x)
        // Extract coeffs of y from f(x, y) assuming it is linear in y.
        // We substitute y (or y(x)) with a pure symbol to use polynomial tools
        const y_sym = new Sym(depVar.name);

        const replaceY = (node) => {
             if (node instanceof Call && node.funcName === depVar.name) return y_sym;
             if (node instanceof Sym && node.name === depVar.name) return y_sym;
             if (node instanceof BinaryOp) return new node.constructor(replaceY(node.left), replaceY(node.right));
             if (node instanceof Call) return new Call(node.funcName, node.args.map(replaceY));
             return node;
        };
        const f_y_poly = replaceY(f_xy);

        // Ensure indepVar doesn't conflict
        const polyInfo = this._getPolyCoeffs(f_y_poly, y_sym);

        if (polyInfo && polyInfo.maxDeg === 1) {
             // Form: y' = coeff[0] + coeff[1]*y
             // y' - coeff[1]*y = coeff[0]
             // Standard Linear: y' + P(x)y = Q(x)
             // P(x) = -coeff[1], Q(x) = coeff[0]

             const Q = polyInfo.coeffs[0] || new Num(0);
             const negP = polyInfo.coeffs[1] || new Num(0);
             const P = new Mul(new Num(-1), negP).simplify();

             // Integrating Factor I = exp(int(P dx))
             let intP = P.integrate(indepVar).simplify();
             const I = new Call('exp', [intP]).simplify();

             // y = (1/I) * (int(Q*I dx) + C)
             const QI = new Mul(Q, I).simplify();
             let intQI = QI.integrate(indepVar).simplify();

             // If integration failed (returns Call('integrate')), we might return symbolic result
             const C = new Sym('C1');
             const right = new Mul(new Div(new Num(1), I), new Add(intQI, C)).simplify();
             return right;
        }

        // 3. Separable: y' = g(x)h(y)
        // Check if f(x,y) / y depends on y? (y'=y) -> covered by linear
        // Check if f(x,y) can be factored?
        // Basic check: y' = y^n * g(x) (Bernoulli n!=0,1 or just separable)
        // Or y' = g(y) -> x = int(1/g(y))
        if (!this._dependsOn(f_xy, indepVar)) {
            // Autonomous: y' = g(y)
            // int(1/g(y)) dy = x + C
            const invG = new Div(new Num(1), replaceY(f_xy)).simplify();
            const left = invG.integrate(y_sym).simplify(); // returns in terms of y_sym
            const right = new Add(indepVar, new Sym('C1'));

            // Implicit solution: F(y) = x + C
            // Try to solve for y?
            const solY = this._solve(new Sub(left, right), y_sym);
            if (!(solY instanceof Call && solY.funcName === 'solve')) {
                return solY;
            }
            // Return implicit equation F(y) = x + C
            return new Eq(left.substitute(y_sym, depVar), right);
        }

        // Generic failure
        return new Call('desolve', [eq, depVar]);
    }

    _desolveSystem(eqs, vars) {
        // Linear System X' = AX + B
        // We assume vars contains y1(t), y2(t)...
        // Independent variable t is implicit in derivatives

        const n = eqs.elements.length;
        if (vars.elements.length !== n) throw new Error("Number of equations must match number of variables");

        // 1. Identify independent variable t
        let t = null;
        const findT = (node) => {
            if (node instanceof Call && node.funcName === 'diff') {
                if (node.args.length > 1) {
                    const v = node.args[1];
                    if (!t) t = v;
                    else if (t.name !== v.name) throw new Error("Multiple independent variables found");
                }
            }
            if (node instanceof Call) node.args.forEach(findT);
            if (node instanceof BinaryOp) { findT(node.left); findT(node.right); }
        };

        eqs.elements.forEach(eq => {
             if (eq instanceof Eq) {
                 findT(eq.left); findT(eq.right);
             } else {
                 findT(eq);
             }
        });

        if (!t) t = new Sym('t'); // Default

        // 2. Normalize equations to X' - ... = 0 and extract A
        // We expect form: diff(yi, t) = ...
        // Or linear combinations.
        // Let's assume standard form: diff(y_i, t) = Sum(a_ij * y_j) + f_i(t)

        const A_rows = [];
        const B_vec = [];

        // Map variables to simple symbols for analysis
        // And create a substitution map for the expressions
        const varSyms = [];
        const subMap = []; // { from: Expr, to: Sym }

        vars.elements.forEach(v => {
            let sym;
            if (v instanceof Call) sym = new Sym(v.funcName);
            else if (v instanceof Sym) sym = v;
            else throw new Error("Variables must be symbols or functions");

            varSyms.push(sym);
            // Replace y(t) with y_sym
            // We need deep substitution logic or just use substitute
            subMap.push({ from: v, to: sym });
        });

        // Helper to substitute all dependent vars in an expression
        const toPolyForm = (expr) => {
             let res = expr;
             // We must be careful with order if names overlap, but usually ok.
             // We use a custom recursive replacement to handle Call vs Sym
             const replace = (node) => {
                 // Check if node matches any 'from'
                 for(let k=0; k<n; k++) {
                     const mapping = subMap[k];
                     const f = mapping.from;
                     if (f instanceof Call && node instanceof Call && f.funcName === node.funcName) return mapping.to;
                     if (f instanceof Sym && node instanceof Sym && f.name === node.name) return mapping.to;
                 }
                 if (node instanceof BinaryOp) return new node.constructor(replace(node.left), replace(node.right));
                 if (node instanceof Call) return new Call(node.funcName, node.args.map(replace));
                 if (node instanceof Not) return new Not(replace(node.arg));
                 return node;
             };
             return replace(res);
        };

        // Construct matrix A
        // We assume equations are ordered same as vars: diff(vars[i]) = ...
        // If not, we should match them.

        // Find which equation corresponds to which variable's derivative
        const orderedRHS = new Array(n).fill(null);

        for(let i=0; i<n; i++) {
            let eq = eqs.elements[i];
            if (eq instanceof Eq) {
                // Check LHS
                let lhs = eq.left;
                let rhs = eq.right;

                // If LHS is diff(vars[k], t)
                let foundK = -1;
                for(let k=0; k<n; k++) {
                     const v = vars.elements[k];
                     // diff(v, t)
                     // v might be y(t) or y
                     let isMatch = false;
                     if (lhs instanceof Call && lhs.funcName === 'diff') {
                          const target = lhs.args[0]; // y(t)
                          // Compare target with v
                          if (target instanceof Call && v instanceof Call && target.funcName === v.funcName) isMatch = true;
                          if (target instanceof Sym && v instanceof Sym && target.name === v.name) isMatch = true;
                     }
                     if (isMatch) {
                         foundK = k;
                         break;
                     }
                }

                if (foundK !== -1) {
                    orderedRHS[foundK] = rhs;
                } else {
                    // Maybe RHS has derivative? x = diff(y) - y -> diff(y) = x+y
                    // For now, require explicit form.
                    // Or maybe it was parsed as diff(y,t) - z = 0
                    // In that case eq is Sub.
                    // But we normalized to Eq in step 2?
                    // No, step 2 logic was inside the loop but I removed it in search block replacement.
                    // Let's assume user gives explicit Eq.
                }
            }
        }

        // Check if we found all
        for(let i=0; i<n; i++) {
            if (!orderedRHS[i]) {
                // Try fallback: maybe index i corresponds to i
                let eq = eqs.elements[i];
                if (eq instanceof Eq) orderedRHS[i] = eq.right;
                else throw new Error("System must be explicit ODEs");
            }
        }

        for(let i=0; i<n; i++) {
            let rhs = orderedRHS[i];
            // Convert to polynomial form in terms of varSyms
            const polyRHS = toPolyForm(rhs).expand().simplify();

            const row = [];
            for(let j=0; j<n; j++) {
                // Coeff of varSyms[j]
                const coeff = polyRHS.diff(varSyms[j]).simplify();
                // Check if coeff depends on vars
                for(const v of varSyms) {
                    if (this._dependsOn(coeff, v)) throw new Error("Non-linear system not supported");
                }
                row.push(coeff);
            }
            A_rows.push(new Vec(row));
        }

        const A = new Vec(A_rows);

        // Check for non-homogeneous terms: Remainder = RHS - A*Vars
        // We need to substitute varSyms back or just check polyRHS - sum(coeff*var)
        for(let i=0; i<n; i++) {
            let rhs = orderedRHS[i];
            const polyRHS = toPolyForm(rhs).expand().simplify();

            // Construct linear part
            let linearPart = new Num(0);
            for(let j=0; j<n; j++) {
                const coeff = A_rows[i].elements[j];
                linearPart = new Add(linearPart, new Mul(coeff, varSyms[j])).simplify();
            }

            const remainder = new Sub(polyRHS, linearPart).expand().simplify();
            if (!(remainder instanceof Num && remainder.value === 0)) {
                // Check if remainder depends on t (non-homogeneous) or just clean zero check failure
                // If it depends on t, throw error
                if (this._dependsOn(remainder, t)) {
                    throw new Error("Non-homogeneous systems not supported yet");
                }
                // If it's a constant non-zero, it's also non-homogeneous (X' = AX + B)
                if (Math.abs(remainder.evaluateNumeric()) > 1e-9) {
                     throw new Error("Non-homogeneous systems not supported yet");
                }
            }
        }

        // Solve X' = AX using eigenvalues
        // X(t) = c1 * v1 * exp(lambda1 * t) + ...
        try {
            const evals = this._eigenvals(A);

            // Get unique eigenvalues
            const lambdas = [];
            if (evals instanceof Vec) {
                 // Deduplicate eigenvalues
                 const seen = new Set();
                 evals.elements.forEach(l => {
                     const s = l.toString();
                     if (!seen.has(s)) {
                         seen.add(s);
                         lambdas.push(l);
                     }
                 });
            } else {
                 lambdas.push(evals);
            }

            // For each lambda, find eigenvector
            const terms = [];

            let solIndex = 1;
            for(const lambda of lambdas) {
                // Kernel of (A - lambda*I)
                const rows = n;
                const cols = n;
                const M_minus_lambdaI = [];
                for(let r=0; r<rows; r++) {
                    const mRow = [];
                    for(let c=0; c<cols; c++) {
                        let val = A.elements[r].elements[c];
                        if (r === c) val = new Sub(val, lambda).simplify();
                        mRow.push(val);
                    }
                    M_minus_lambdaI.push(new Vec(mRow));
                }
                const mat = new Vec(M_minus_lambdaI);
                const kern = this._kernel(mat);

                for(const v of kern.elements) {
                    // Term: C_i * v * exp(lambda * t)
                    const C = new Sym(`C${solIndex++}`);
                    const exp = new Call('exp', [new Mul(lambda, t)]).simplify();
                    const termVec = [];
                    for(const el of v.elements) {
                         termVec.push(new Mul(C, new Mul(el, exp)).simplify());
                    }
                    terms.push(new Vec(termVec));
                }
            }

            if (terms.length < n) throw new Error("Defective matrix (Jordan form) not supported");

            // Sum up terms to get X(t)
            // X(t) = term1 + term2 + ...
            const finalSol = [];
            for(let i=0; i<n; i++) {
                let sum = new Num(0);
                for(const term of terms) {
                    sum = new Add(sum, term.elements[i]).simplify();
                }
                finalSol.push(sum);
            }
            return new Vec(finalSol);

        } catch (e) {
            console.error(e);
            return new Call('desolve', [eqs, vars]);
        }
    }

    _factorRational(coeffs, maxDeg, varNode) {
        let currentCoeffs = coeffs;
        let currentMaxDeg = maxDeg;
        const foundFactors = [];
        let foundAny = false;

        // Try to find roots until degree <= 2 or no more rational roots
        // Safety loop
        for (let iter = 0; iter < maxDeg; iter++) {
            if (currentMaxDeg <= 2) break;

            const a0 = currentCoeffs[0] || new Num(0);
            const an = currentCoeffs[currentMaxDeg] || new Num(0);

            // Rational Root Theorem requires integer coefficients (or at least we need to scale them)
            // For now, assume integer coefficients or simple floats that are integers
            if (!Number.isInteger(a0.value) || !Number.isInteger(an.value)) break;

            const div0 = this._divisorsList(Math.abs(Math.round(a0.value)));
            const divN = this._divisorsList(Math.abs(Math.round(an.value)));

            const candidates = new Set();
            for (const p of div0) {
                for (const q of divN) {
                    candidates.add(p / q);
                    candidates.add(-p / q);
                }
            }
            candidates.add(0); // Should be handled by x extraction, but safe to add

            let foundRoot = false;
            for (const r of candidates) {
                // Evaluate P(r) using Horner
                let val = 0;
                for (let i = currentMaxDeg; i >= 0; i--) {
                    const c = (currentCoeffs[i] ? currentCoeffs[i].value : 0);
                    val = val * r + c;
                }

                if (Math.abs(val) < 1e-9) {
                    // Root found!
                    // Check if r is integer or fraction
                    // If r = p/q, factor is (q*x - p) usually to keep integer coeffs,
                    // but we simplify to (x - r) * an_new?
                    // Let's use (x - r).

                    // Construct factor
                    let factor = new Sub(varNode, new Num(r)).simplify();
                    foundFactors.push(factor);

                    // Divide P(x) / (x - r)
                    const res = this._polyDivNumeric(currentCoeffs, currentMaxDeg, r);
                    currentCoeffs = res.coeffs;
                    currentMaxDeg = res.maxDeg;
                    foundRoot = true;
                    foundAny = true;
                    break; // Restart loop with reduced poly to re-evaluate candidates for multiplicity
                }
            }

            if (!foundRoot) break; // No more rational roots found
        }

        return {
            found: foundAny,
            factors: foundFactors,
            newCoeffs: currentCoeffs,
            newMaxDeg: currentMaxDeg
        };
    }

    _divisorsList(n) {
        const res = [];
        if (n === 0) return [1]; // Handling 0 case?
        for (let i = 1; i * i <= n; i++) {
            if (n % i === 0) {
                res.push(i);
                if (i * i !== n) res.push(n / i);
            }
        }
        return res;
    }

    _polyDivNumeric(coeffs, maxDeg, r) {
        // Synthetic division of P(x) by (x - r)
        // Returns { coeffs, maxDeg } of quotient Q(x)
        // Q[n-1] = P[n]
        // Q[k-1] = P[k] + r * Q[k]

        // My coeffs are indexed by power 0..maxDeg
        // P[maxDeg] is coeff of x^maxDeg
        // Q will have degree maxDeg - 1

        const Q = {};
        let val = coeffs[maxDeg].value;
        Q[maxDeg - 1] = new Num(val);

        for (let i = maxDeg - 1; i > 0; i--) {
            const pi = coeffs[i] ? coeffs[i].value : 0;
            val = pi + r * val;
            Q[i - 1] = new Num(val);
        }

        return {
            coeffs: Q,
            maxDeg: maxDeg - 1
        };
    }

    _partfrac(expr, varNode) {
        if (!(varNode instanceof Sym)) throw new Error("Second argument to partfrac must be a variable");

        expr = expr.simplify();

        // 1. Check if Div
        if (!(expr instanceof Div)) return expr;

        const num = expr.left;
        const den = expr.right;

        // 2. Find roots of denominator
        const rootsResult = this._solve(den, varNode);
        let roots = [];

        if (rootsResult instanceof Call && rootsResult.funcName === 'set') {
            roots = rootsResult.args;
        } else if (rootsResult instanceof Expr && !(rootsResult instanceof Call && rootsResult.funcName === 'solve')) {
            roots = [rootsResult];
        } else {
            // Could not solve or no roots
            return new Call('partfrac', [expr, varNode]);
        }

        // 3. Construct Partial Fractions (assuming simple roots for now)
        // form: sum( Residue_i / (x - r_i) )
        // Residue at simple pole r is P(r) / Q'(r)

        let result = new Num(0);
        let validDecomp = true;

        // Pre-calculate Q'(x)
        const denDiff = den.diff(varNode).simplify();

        for (const r of roots) {
            // Factor: (var - r)
            const factor = new Sub(varNode, r);

            // Calculate Residue: P(r) / Q'(r)

            let residue;

            try {
                const numVal = num.substitute(varNode, r).simplify();
                const denDiffVal = denDiff.substitute(varNode, r).simplify();

                // If denDiffVal is 0, it's a repeated root (order > 1)
                if ((denDiffVal instanceof Num && denDiffVal.value === 0) || (denDiffVal instanceof Sym && denDiffVal.name === 'NaN')) {
                    // limit((x-r)*expr, x, r)
                    const term = new Mul(expr, factor).simplify();
                    residue = this._limit(term, varNode, r);
                } else {
                    residue = new Div(numVal, denDiffVal).simplify();
                }

                if (residue instanceof Sym && (residue.name === 'Infinity' || residue.name === 'NaN')) {
                    validDecomp = false; break;
                }
                // Check if residue still depends on var (should be constant)
                if (this._dependsOn(residue, varNode)) {
                     validDecomp = false; break;
                }
            } catch (e) {
                validDecomp = false; break;
            }

            const frac = new Div(residue, factor);
            result = new Add(result, frac);
        }

        // Check for NaN result (global invalidation)
        if (result.toString() === "NaN" || result.toString() === "(NaN / NaN)") validDecomp = false;

        if (validDecomp) {
            return result.simplify();
        }

        return new Call('partfrac', [expr, varNode]);
    }

    // --- Polynomial Tools ---

    _dependsOn(expr, varNode) {
        if (expr instanceof Sym) return expr.name === varNode.name;
        if (expr instanceof Num) return false;
        if (expr instanceof BinaryOp) return this._dependsOn(expr.left, varNode) || this._dependsOn(expr.right, varNode);
        if (expr instanceof Call) {
            if (expr.funcName === varNode.name) return true; // Check if function name matches variable (e.g. y(x) depends on y)
            return expr.args.some(a => this._dependsOn(a, varNode));
        }
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

        // Re-calculate maxDeg to ignore zero coefficients
        let trueMaxDeg = 0;
        let foundAny = false;
        // Use Object.keys to iterate, but sort or check all
        for (const d in coeffs) {
            const deg = parseInt(d);
            const c = coeffs[d];
            // Check if non-zero
            if (!(c instanceof Num && c.value === 0)) {
                 if (deg > trueMaxDeg) trueMaxDeg = deg;
                 foundAny = true;
            }
        }
        if (!foundAny) trueMaxDeg = 0;

        return { coeffs, maxDeg: trueMaxDeg };
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
