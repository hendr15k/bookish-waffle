
class CAS {
    constructor() {
        this.variables = {
            'pi': new Num(Math.PI),
            'true': new Num(1),
            'false': new Num(0),
            'j': new Sym('i') // Electrical Engineering imaginary unit
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
            // Special handling for purge/del to prevent premature evaluation of the argument
            if (node.funcName === 'purge' || node.funcName === 'del') {
                if (node.args.length !== 1) throw new Error("purge requires 1 argument");
                const target = node.args[0];
                if (target instanceof Sym) {
                    if (this.variables.hasOwnProperty(target.name)) {
                        delete this.variables[target.name];
                        return new Num(1); // Return 1 for success
                    }
                    return new Num(0); // Not found
                }
                throw new Error("purge argument must be a variable name");
            }

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

                    // Try strategies if direct integration failed
                    if (indefinite instanceof Call && indefinite.funcName === 'integrate') {
                         // 1. Substitution
                         const subRes = this._integrateSubstitution(func, varNode);
                         if (subRes && !(subRes instanceof Call && subRes.funcName === 'integrate')) {
                             indefinite = subRes;
                         } else {
                             // 2. Parts
                             const partsRes = this._integrateByParts(func, varNode);
                             if (partsRes && !(partsRes instanceof Call && partsRes.funcName === 'integrate')) {
                                 indefinite = partsRes;
                             }
                         }
                    }

                    // Attempt partfrac if still failed
                    if (indefinite instanceof Call && indefinite.funcName === 'integrate') {
                        const pf = this._partfrac(func, varNode);
                        if (!(pf instanceof Call && pf.funcName === 'partfrac') && pf.toString() !== func.toString()) {
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
                     // 1. Try Integration by Substitution (Logarithmic form f'/f and others)
                     const subRes = this._integrateSubstitution(func, varNode);
                     if (subRes && !(subRes instanceof Call && subRes.funcName === 'integrate')) return subRes;

                     // 2. Try Integration by Parts
                     const partsRes = this._integrateByParts(func, varNode);
                     if (partsRes && !(partsRes instanceof Call && partsRes.funcName === 'integrate')) return partsRes;

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

                     // Try Integration by Parts
                     const parts = this._integrateByParts(func, varNode);
                     if (parts) return parts;
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

                 // Check for guess (Newton-Raphson)
                 if (args.length === 3) {
                     const guess = args[2];
                     if (!(varNode instanceof Sym)) throw new Error("Numeric solve requires a single variable");
                     return this._fsolve(eq, varNode, guess);
                 }

                 if (!(varNode instanceof Sym) && !(varNode instanceof Vec)) throw new Error("Second argument to solve must be a variable or list of variables");
                 return this._solve(eq, varNode);
            }

            if (node.funcName === 'nIntegrate' || node.funcName === 'numeric_integrate') {
                if (args.length !== 4) throw new Error("nIntegrate requires 4 arguments: expr, var, start, end");
                return this._nIntegrate(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'minimize') {
                 if (args.length !== 2) throw new Error("minimize requires 2 arguments: expr, variable");
                 if (!(args[1] instanceof Sym)) throw new Error("Second argument to minimize must be a variable");
                 return this._minimize(args[0], args[1]);
            }

            if (node.funcName === 'maximize') {
                 if (args.length !== 2) throw new Error("maximize requires 2 arguments: expr, variable");
                 if (!(args[1] instanceof Sym)) throw new Error("Second argument to maximize must be a variable");
                 return this._maximize(args[0], args[1]);
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

            if (node.funcName === 'plot3d') {
                // plot3d(expr, x, y, [x_min, x_max, y_min, y_max])
                if (args.length < 3) throw new Error("plot3d requires at least 3 arguments: expression, var1, var2");
                const expr = args[0];
                const varX = args[1];
                const varY = args[2];
                let xMin = -10, xMax = 10, yMin = -10, yMax = 10;

                if (args.length >= 7) {
                    xMin = args[3].evaluateNumeric();
                    xMax = args[4].evaluateNumeric();
                    yMin = args[5].evaluateNumeric();
                    yMax = args[6].evaluateNumeric();
                } else if (args.length >= 5) {
                    xMin = args[3].evaluateNumeric();
                    xMax = args[4].evaluateNumeric();
                }

                if (isNaN(xMin)) xMin = -10;
                if (isNaN(xMax)) xMax = 10;
                if (isNaN(yMin)) yMin = -10;
                if (isNaN(yMax)) yMax = 10;

                return {
                    type: 'plot',
                    subtype: '3d',
                    expr: expr,
                    varX: varX,
                    varY: varY,
                    xMin: xMin,
                    xMax: xMax,
                    yMin: yMin,
                    yMax: yMax,
                    toString: () => `3D Plot ${expr}`,
                    toLatex: () => `\\text{3D Plot } ${expr.toLatex()}`
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

            if (node.funcName === 'tangent') {
                // tangent(expr, var, point)
                if (args.length < 3) throw new Error("tangent requires 3 arguments: expression, variable, point");
                const expr = args[0];
                const varNode = args[1];
                const point = args[2];
                if (!(varNode instanceof Sym)) throw new Error("Second argument to tangent must be a variable");

                // f(a) + f'(a)(x-a)
                const f_a = expr.substitute(varNode, point).simplify();
                const deriv = expr.diff(varNode).simplify();
                const f_prime_a = deriv.substitute(varNode, point).simplify();

                // y = f(a) + f'(a)*(x - a)
                const term = new Mul(f_prime_a, new Sub(varNode, point)).simplify();
                return new Add(f_a, term).simplify();
            }

            if (node.funcName === 'trigReduce') {
                if (args.length !== 1) throw new Error("trigReduce requires 1 argument");
                return this._linearizeTrig(args[0]);
            }

            if (node.funcName === 'trigExpand') {
                if (args.length !== 1) throw new Error("trigExpand requires 1 argument");
                return args[0].expand().simplify();
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

            if (node.funcName === 'trans' || node.funcName === 'transpose' || node.funcName === 'tran') {
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
                if (args.length === 1 && args[0] instanceof Vec) {
                    const list = args[0].elements;
                    if (list.length === 0) return new Num(0);
                    let res = list[0];
                    for(let i=1; i<list.length; i++) res = this._gcd(res, list[i]);
                    return res;
                }
                if (args.length < 2) throw new Error("gcd requires at least 2 arguments");
                let res = args[0];
                for(let i=1; i<args.length; i++) res = this._gcd(res, args[i]);
                return res;
            }

            if (node.funcName === 'lcm') {
                if (args.length === 1 && args[0] instanceof Vec) {
                    const list = args[0].elements;
                    if (list.length === 0) return new Num(1);
                    let res = list[0];
                    for(let i=1; i<list.length; i++) res = this._lcm(res, list[i]);
                    return res;
                }
                if (args.length < 2) throw new Error("lcm requires at least 2 arguments");
                let res = args[0];
                for(let i=1; i<args.length; i++) res = this._lcm(res, args[i]);
                return res;
            }

            if (node.funcName === 'factor' || node.funcName === 'ifactor') {
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

            if (node.funcName === 'nPr' || node.funcName === 'perm') {
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

            if (node.funcName === 'polyRegression') {
                if (args.length < 2) throw new Error("polyRegression requires 2 arguments: data, degree");
                return this._polyRegression(args[0], args[1]);
            }

            if (node.funcName === 'expRegression') {
                if (args.length !== 1) throw new Error("expRegression requires 1 argument: data");
                return this._expRegression(args[0]);
            }

            if (node.funcName === 'normalPDF') {
                if (args.length !== 3) throw new Error("normalPDF requires 3 arguments: x, mu, sigma");
                return this._normalPDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'binomialPDF') {
                if (args.length !== 3) throw new Error("binomialPDF requires 3 arguments: k, n, p");
                return this._binomialPDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'normalCDF') {
                if (args.length !== 3) throw new Error("normalCDF requires 3 arguments: x, mu, sigma");
                return this._normalCDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'invNorm') {
                if (args.length !== 3) throw new Error("invNorm requires 3 arguments: area, mu, sigma");
                return this._invNorm(args[0], args[1], args[2]);
            }

            if (node.funcName === 'binomialCDF') {
                if (args.length !== 3) throw new Error("binomialCDF requires 3 arguments: k, n, p");
                return this._binomialCDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'poissonPDF') {
                if (args.length !== 2) throw new Error("poissonPDF requires 2 arguments: k, lambda");
                return this._poissonPDF(args[0], args[1]);
            }

            if (node.funcName === 'poissonCDF') {
                if (args.length !== 2) throw new Error("poissonCDF requires 2 arguments: k, lambda");
                return this._poissonCDF(args[0], args[1]);
            }

            if (node.funcName === 'exponentialPDF') {
                if (args.length !== 2) throw new Error("exponentialPDF requires 2 arguments: x, lambda");
                return this._exponentialPDF(args[0], args[1]);
            }

            if (node.funcName === 'exponentialCDF') {
                if (args.length !== 2) throw new Error("exponentialCDF requires 2 arguments: x, lambda");
                return this._exponentialCDF(args[0], args[1]);
            }

            if (node.funcName === 'geometricPDF') {
                if (args.length !== 2) throw new Error("geometricPDF requires 2 arguments: k, p");
                return this._geometricPDF(args[0], args[1]);
            }

            if (node.funcName === 'geometricCDF') {
                if (args.length !== 2) throw new Error("geometricCDF requires 2 arguments: k, p");
                return this._geometricCDF(args[0], args[1]);
            }

            if (node.funcName === 'chisquarePDF') {
                if (args.length !== 2) throw new Error("chisquarePDF requires 2 arguments: x, k");
                return this._chisquarePDF(args[0], args[1]);
            }

            if (node.funcName === 'compound') {
                if (args.length !== 4) throw new Error("compound requires 4 arguments: P, r, n, t");
                return this._compound(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'loan') {
                if (args.length !== 3) throw new Error("loan requires 3 arguments: P, r, n");
                return this._loan(args[0], args[1], args[2]);
            }

            if (node.funcName === 'npv') {
                if (args.length !== 2) throw new Error("npv requires 2 arguments: rate, cash_flows");
                return this._npv(args[0], args[1]);
            }

            if (node.funcName === 'irr') {
                if (args.length !== 1) throw new Error("irr requires 1 argument: cash_flows");
                return this._irr(args[0]);
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

            if (node.funcName === 'rem' || node.funcName === 'irem') {
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

            if (node.funcName === 'size' || node.funcName === 'dim' || node.funcName === 'length') {
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

            if (node.funcName === 'distance') {
                if (args.length !== 2) throw new Error("distance requires 2 arguments");
                return this._distance(args[0], args[1]);
            }

            if (node.funcName === 'midpoint') {
                if (args.length !== 2) throw new Error("midpoint requires 2 arguments");
                return this._midpoint(args[0], args[1]);
            }

            if (node.funcName === 'zTest') {
                if (args.length !== 3) throw new Error("zTest requires 3 arguments: data, mu0, sigma");
                return this._zTest(args[0], args[1], args[2]);
            }

            if (node.funcName === 'euler' || node.funcName === 'phi') {
                if (args.length !== 1) throw new Error("euler requires 1 argument");
                return this._euler(args[0]);
            }

            if (node.funcName === 'mode') {
                if (args.length !== 1) throw new Error("mode requires 1 argument (list)");
                return this._mode(args[0]);
            }

            if (node.funcName === 'arcLen') {
                if (args.length !== 4) throw new Error("arcLen requires 4 arguments: expr, var, start, end");
                return this._arcLen(args[0], args[1], args[2], args[3]);
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

            if (node.funcName === 'fourier') {
                if (args.length < 3) throw new Error("fourier requires 3 arguments: expr, var, n, [L]");
                const expr = args[0];
                const varNode = args[1];
                const n = args[2];
                const L = args.length > 3 ? args[3] : new Sym('pi');
                if (!(varNode instanceof Sym)) throw new Error("Second argument to fourier must be a variable");
                return this._fourier(expr, varNode, n, L);
            }

            if (node.funcName === 'slopefield') {
                 // slopefield(diffEq, x, y, [minX, maxX, minY, maxY])
                 if (args.length < 3) throw new Error("slopefield requires at least 3 arguments: equation, x, y");
                 const eq = args[0];
                 const xVar = args[1];
                 const yVar = args[2];

                 let expr = eq;
                 // If eq is Eq (y' = ...), extract RHS
                 if (eq instanceof Eq) {
                     expr = eq.right;
                 }

                 let min = -10, max = 10, yMin = -10, yMax = 10;
                 if (args.length >= 7) {
                     min = args[3].evaluateNumeric();
                     max = args[4].evaluateNumeric();
                     yMin = args[5].evaluateNumeric();
                     yMax = args[6].evaluateNumeric();
                 }

                 return {
                     type: 'plot',
                     subtype: 'slopefield',
                     expr: expr,
                     varX: xVar,
                     varY: yVar,
                     min: isNaN(min) ? -10 : min,
                     max: isNaN(max) ? 10 : max,
                     yMin: isNaN(yMin) ? -10 : yMin,
                     yMax: isNaN(yMax) ? 10 : yMax,
                     toString: () => `Slope Field ${expr}`,
                     toLatex: () => `\\text{Slope Field } y' = ${expr.toLatex()}`
                 };
            }

            if (node.funcName === 'vectorfield') {
                 // vectorfield([u, v], x, y, [minX, maxX, minY, maxY])
                 if (args.length < 3) throw new Error("vectorfield requires at least 3 arguments: vector, x, y");
                 const vec = args[0];
                 const xVar = args[1];
                 const yVar = args[2];

                 if (!(vec instanceof Vec) || vec.elements.length !== 2) throw new Error("Vector field must be 2D vector [u, v]");

                 let min = -10, max = 10, yMin = -10, yMax = 10;
                 if (args.length >= 7) {
                     min = args[3].evaluateNumeric();
                     max = args[4].evaluateNumeric();
                     yMin = args[5].evaluateNumeric();
                     yMax = args[6].evaluateNumeric();
                 }

                 return {
                     type: 'plot',
                     subtype: 'vectorfield',
                     exprX: vec.elements[0],
                     exprY: vec.elements[1],
                     varX: xVar,
                     varY: yVar,
                     min: isNaN(min) ? -10 : min,
                     max: isNaN(max) ? 10 : max,
                     yMin: isNaN(yMin) ? -10 : yMin,
                     yMax: isNaN(yMax) ? 10 : yMax,
                     toString: () => `Vector Field ${vec}`,
                     toLatex: () => `\\text{Vector Field } ${vec.toLatex()}`
                 };
            }

            if (node.funcName === 'plotimplicit') {
                 if (args.length < 3) throw new Error("plotimplicit requires at least 3 arguments: equation, x, y");
                 const eq = args[0];
                 const xVar = args[1];
                 const yVar = args[2];

                 let xMin = -10, xMax = 10, yMin = -10, yMax = 10;
                 if (args.length >= 7) {
                     xMin = args[3].evaluateNumeric();
                     xMax = args[4].evaluateNumeric();
                     yMin = args[5].evaluateNumeric();
                     yMax = args[6].evaluateNumeric();
                 }

                 // Normalize Equation: LHS - RHS
                 let expr = eq;
                 if (eq instanceof Eq) {
                     expr = new Sub(eq.left, eq.right).simplify();
                 }

                 return {
                     type: 'plot',
                     subtype: 'implicit',
                     expr: expr,
                     varX: xVar,
                     varY: yVar,
                     xMin: isNaN(xMin) ? -10 : xMin,
                     xMax: isNaN(xMax) ? 10 : xMax,
                     yMin: isNaN(yMin) ? -10 : yMin,
                     yMax: isNaN(yMax) ? 10 : yMax,
                     toString: () => `Implicit Plot ${eq}`,
                     toLatex: () => `\\text{Implicit Plot } ${eq.toLatex()}`
                 };
            }

            if (node.funcName === 'help') {
                const helpText = `Available commands:
diff, integrate, limit, taylor, sum, product,
expand, simplify, solve, minimize, maximize,
det, inv, trans, cross, dot, norm, grad, curl, divergence,
gcd, lcm, factor, nCr, nPr, isPrime, factorial,
mean, median, variance, linearRegression,
nIntegrate, fsolve,
normalPDF, normalCDF, invNorm, binomialPDF, binomialCDF,
poissonPDF, poissonCDF, exponentialPDF, exponentialCDF,
geometricPDF, geometricCDF, chisquarePDF,
compound, loan,
degree, coeff, symb2poly, poly2symb,
seq, range, sort, reverse,
diag, identity,
laplace, ilaplace,
rem, quo, mod, arg, approx, erf,
size, concat, clear, N,
molarMass, atomicWeight,
and, or, not, xor, int, evalf, purge,
perm, tran, irem, ifactor`;

                const latexHelp = `\\text{Available commands: see documentation}`;
                return { type: 'info', text: helpText, toString: () => helpText, toLatex: () => latexHelp };
            }

            if (node.funcName === 'molarMass') {
                if (args.length !== 1) throw new Error("molarMass requires 1 argument (string formula)");
                // Argument comes as a symbol (e.g. H2O) or string if we support string literals?
                // The parser parses identifiers as Sym. So H2O is a Sym.
                // Or "H2O" if we have string support? Current Lexer doesn't seem to support strings.
                // So user types molarMass(H2O).

                let formula = "";
                if (args[0] instanceof Sym) formula = args[0].name;
                // If we want to support molarMass(C6H12O6), that parses as Call(C6H12O6)? No, C6H12O6 is valid identifier.

                if (!formula) throw new Error("Invalid formula");
                return this._molarMass(formula);
            }

            if (node.funcName === 'atomicWeight') {
                if (args.length !== 1) throw new Error("atomicWeight requires 1 argument (symbol)");
                let sym = "";
                if (args[0] instanceof Sym) sym = args[0].name;
                if (!sym) throw new Error("Invalid element symbol");
                return this._atomicWeight(sym);
            }

            if (node.funcName === 'balance') {
                if (args.length !== 1) throw new Error("balance requires 1 argument (equation string)");
                // Try to extract string from Sym or just use toString
                let eq = args[0].toString();
                // If it was parsed as subtraction (A-B), reconstruct?
                // The parser might parse "H2 + O2 -> H2O" weirdly if -> is not operator.
                // Lexer has -> as IMPLIES. So it becomes Implies(Add(H2, O2), H2O).
                // Or if user passes string "H2+O2->H2O".
                // Our parser doesn't support string literals yet?
                // Actually, Lexer doesn't support quotes.
                // So users must type: balance(H2 + O2 -> H2O)
                // This parses as Implies(Add(H2, O2), H2O)
                // We need to reconstruct the string or handle the AST.
                // Let's reconstruct string from AST for now, or traverse AST.
                return this._balanceChem(node.args[0]);
            }

            if (node.funcName === 'diagonalize') {
                if (args.length !== 1) throw new Error("diagonalize requires 1 argument (matrix)");
                return this._diagonalize(args[0]);
            }

            if (node.funcName === 'tTest') {
                if (args.length !== 2) throw new Error("tTest requires 2 arguments: data, mu0");
                return this._tTest(args[0], args[1]);
            }

            if (node.funcName === 'kron') {
                if (args.length !== 2) throw new Error("kron requires 2 arguments");
                return this._kron(args[0], args[1]);
            }

            if (node.funcName === 'svd') {
                if (args.length !== 1) throw new Error("svd requires 1 argument");
                return this._svd(args[0]);
            }

            if (node.funcName === 'geoMean') {
                if (args.length !== 1) throw new Error("geoMean requires 1 argument");
                return this._geoMean(args[0]);
            }

            if (node.funcName === 'harmMean') {
                if (args.length !== 1) throw new Error("harmMean requires 1 argument");
                return this._harmMean(args[0]);
            }

            if (node.funcName === 'rms') {
                if (args.length !== 1) throw new Error("rms requires 1 argument");
                return this._rms(args[0]);
            }

            if (node.funcName === 'mad') {
                if (args.length !== 1) throw new Error("mad requires 1 argument");
                return this._mad(args[0]);
            }

            if (node.funcName === 'curvature') {
                if (args.length < 2) throw new Error("curvature requires at least 2 arguments: expr, var, [point]");
                const point = args.length > 2 ? args[2] : null;
                return this._curvature(args[0], args[1], point);
            }

            if (node.funcName === 'par' || node.funcName === 'parallel') {
                return this._parallel(args);
            }

            if (node.funcName === 'cis') {
                if (args.length !== 1) throw new Error("cis requires 1 argument (angle in degrees)");
                return this._cis(args[0]);
            }

            if (node.funcName === 'phasor') {
                if (args.length !== 2) throw new Error("phasor requires 2 arguments: magnitude, angle(deg)");
                return this._phasor(args[0], args[1]);
            }

            if (node.funcName === 'toPolar') {
                if (args.length !== 1) throw new Error("toPolar requires 1 argument");
                return this._toPolar(args[0]);
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

        if (node instanceof If) {
            // evaluate condition
            const cond = this._recursiveEval(node.condition).simplify();
            let val = cond.evaluateNumeric();
            if (isNaN(val)) {
                // If strictly not a number, maybe it's boolean 1 or 0 in CAS
                if (cond instanceof Num) val = cond.value;
            }
            if (isNaN(val)) throw new Error("Condition must evaluate to a number");

            if (val !== 0) {
                // True
                return this._recursiveEval(node.trueBlock);
            } else {
                if (node.falseBlock) {
                    return this._recursiveEval(node.falseBlock);
                }
                return new Num(0); // Default undefined
            }
        }

        if (node instanceof While) {
            let limit = 10000;
            let result = new Num(0);
            while(limit-- > 0) {
                const cond = this._recursiveEval(node.condition).simplify();
                let val = cond.evaluateNumeric();
                if (isNaN(val)) {
                    if (cond instanceof Num) val = cond.value;
                }
                if (isNaN(val)) throw new Error("Condition must evaluate to a number");
                if (val === 0) break;

                const res = this._recursiveEval(node.body);
                if (res instanceof Return) return res;
                if (res instanceof Break) break;
                if (res instanceof Continue) continue;
                result = res;
            }
            if (limit <= 0) throw new Error("Infinite loop detected");
            return result;
        }

        if (node instanceof For) {
            let limit = 10000;
            let result = new Num(0);

            if (node.init) this._recursiveEval(node.init);

            while(limit-- > 0) {
                if (node.condition) {
                    const cond = this._recursiveEval(node.condition).simplify();
                    let val = cond.evaluateNumeric();
                    if (isNaN(val)) {
                        if (cond instanceof Num) val = cond.value;
                    }
                    if (isNaN(val)) throw new Error("Condition must evaluate to a number");
                    if (val === 0) break;
                }

                const res = this._recursiveEval(node.body);
                if (res instanceof Return) return res;
                if (res instanceof Break) break;
                if (res instanceof Continue) {
                    if (node.step) this._recursiveEval(node.step);
                    continue;
                }
                result = res;

                if (node.step) this._recursiveEval(node.step);
            }
            if (limit <= 0) throw new Error("Infinite loop detected");
            return result;
        }

        if (node instanceof Return || node instanceof Break || node instanceof Continue) {
            return node;
        }

        if (node instanceof Block) {
             const results = [];
             for(const stmt of node.statements) {
                 const res = this._recursiveEval(stmt);
                 if (res instanceof Return) return res.value; // Unpack return value
                 if (res instanceof Break || res instanceof Continue) return res; // Propagate up
                 results.push(res);
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

    _minimize(expr, varNode) {
        return this._optimize(expr, varNode, 'min');
    }

    _maximize(expr, varNode) {
        return this._optimize(expr, varNode, 'max');
    }

    _optimize(expr, varNode, type) {
        // 1. Differentiate
        const deriv = expr.diff(varNode).simplify();

        // 2. Solve deriv = 0
        const criticalPoints = this._solve(deriv, varNode);

        // Extract points
        let points = [];
        if (criticalPoints instanceof Call && criticalPoints.funcName === 'set') {
            points = criticalPoints.args;
        } else if (criticalPoints instanceof Expr && !(criticalPoints instanceof Call && criticalPoints.funcName === 'solve')) {
            points = [criticalPoints];
        } else {
            // Failed to solve
            return new Call(type === 'min' ? 'minimize' : 'maximize', [expr, varNode]);
        }

        // 3. Second Derivative Test
        const deriv2 = deriv.diff(varNode).simplify();
        const results = [];

        for(const pt of points) {
            try {
                const val2 = deriv2.substitute(varNode, pt).simplify();
                // If min: val2 > 0
                // If max: val2 < 0
                const numVal = val2.evaluateNumeric();
                let isMatch = false;
                if (!isNaN(numVal)) {
                     if (type === 'min' && numVal > 0) isMatch = true;
                     if (type === 'max' && numVal < 0) isMatch = true;
                } else {
                     // Symbolic check? Assume match if we can't prove otherwise?
                     // Or return symbolic condition?
                     // Let's include it if we can't determine, user can check.
                     isMatch = true;
                }

                if (isMatch) {
                    // Return the value f(x) or the point x?
                    // Usually minimize returns the minimal VALUE? Or the point?
                    // Standard: minimize(f) -> value. argmin(f) -> point.
                    // But Xcas minimize(f) returns value.
                    const val = expr.substitute(varNode, pt).simplify();
                    results.push(val);
                }
            } catch(e) {}
        }

        if (results.length === 0) return new Call('set', []);
        if (results.length === 1) return results[0];
        return new Call('set', results);
    }

    _fsolve(eq, varNode, guess) {
        // Newton-Raphson Method
        // x_{n+1} = x_n - f(x_n) / f'(x_n)

        let expr;
        if (eq instanceof Eq) {
            expr = new Sub(eq.left, eq.right).simplify();
        } else {
            expr = eq.simplify();
        }

        const deriv = expr.diff(varNode).simplify();
        let x = guess.evaluateNumeric();
        if (isNaN(x)) throw new Error("Initial guess must evaluate to a number");

        const maxIter = 100;
        const tol = 1e-9;

        for (let i = 0; i < maxIter; i++) {
            const fVal = expr.substitute(varNode, new Num(x)).evaluateNumeric();
            const dVal = deriv.substitute(varNode, new Num(x)).evaluateNumeric();

            if (Math.abs(dVal) < 1e-12) throw new Error("Derivative is zero during iteration");

            const xNew = x - fVal / dVal;

            if (Math.abs(xNew - x) < tol) {
                return new Num(xNew);
            }
            x = xNew;
        }

        throw new Error("Numeric solver did not converge");
    }

    _nIntegrate(expr, varNode, start, end) {
        // Simpson's Rule
        // int(f) approx (b-a)/6 * (f(a) + 4f((a+b)/2) + f(b))
        // Composite Simpson's Rule

        const a = start.evaluateNumeric();
        const b = end.evaluateNumeric();

        if (isNaN(a) || isNaN(b)) throw new Error("Integration bounds must be numeric");

        const n = 100; // Even number of intervals
        const h = (b - a) / n;

        const f = (val) => {
            const res = expr.substitute(varNode, new Num(val)).evaluateNumeric();
            if (isNaN(res)) throw new Error("Function evaluation failed at " + val);
            return res;
        };

        let sum = f(a) + f(b);

        for (let i = 1; i < n; i++) {
            const x = a + i * h;
            if (i % 2 === 0) sum += 2 * f(x);
            else sum += 4 * f(x);
        }

        return new Num((h / 3) * sum);
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

        // 2. General Substitution u = g(x)
        // Candidates for u: arguments of functions, bases of powers, exponents
        const candidates = new Set();
        const findCandidates = (node) => {
             if (node instanceof Call) {
                 if (node.funcName !== 'integrate') {
                     node.args.forEach(arg => {
                         if (this._dependsOn(arg, varNode)) candidates.add(arg);
                     });
                 }
             } else if (node instanceof Pow) {
                 if (this._dependsOn(node.left, varNode)) candidates.add(node.left);
                 if (this._dependsOn(node.right, varNode)) candidates.add(node.right);
             }
             if (node instanceof BinaryOp) { findCandidates(node.left); findCandidates(node.right); }
             if (node instanceof Call) node.args.forEach(findCandidates);
        };
        findCandidates(expr);

        for(const u of candidates) {
            // Check if u is trivial (x)
            if (u instanceof Sym && u.name === varNode.name) continue;

            // du = u'(x)
            const du = u.diff(varNode).simplify();

            // Check if expr / du depends on x ONLY through u
            // This is hard to check perfectly.
            // Simplified check: Divide expr by du. If result is constant w.r.t x, it works.
            // Or if result can be written as H(u).
            // Example: 2x * e^(x^2). u=x^2, du=2x. expr/du = e^(x^2) = e^u.
            // Example: x * e^(x^2). u=x^2, du=2x. expr/du = 1/2 * e^(x^2).

            try {
                const ratio = new Div(expr, du).simplify();

                // If ratio depends on varNode, we check if it can be transformed to u.
                // We use a heuristic: substitute u with a temp symbol 'U' in ratio.
                // If 'U' remains and varNode is gone, we succeed.
                // But substitution is structural. e^(x^2) matches u=x^2.
                // x^4 matches u=x^2? substitute returns U^2? No, simple substitute doesn't do algebra.

                const U = new Sym('__U__');
                let transformed = ratio.substitute(u, U).simplify();

                // If simplified transformed expression still has varNode, we try to see if remaining parts form u?
                // For now, strict substitution check.
                if (!this._dependsOn(transformed, varNode)) {
                     // Integrate transformed w.r.t U
                     const integral = transformed.integrate(U).simplify();
                     if (!(integral instanceof Call && integral.funcName === 'integrate')) {
                         // Substitute U back to u
                         return integral.substitute(U, u).simplify();
                     }
                }
            } catch (e) {}
        }

        return new Call("integrate", [expr, varNode]);
    }

    _integrateByParts(expr, varNode) {
        let factors = [];

        if (expr instanceof Mul) {
            const flattenMul = (node) => {
                if (node instanceof Mul) {
                    flattenMul(node.left);
                    flattenMul(node.right);
                } else {
                    factors.push(node);
                }
            };
            flattenMul(expr);
        } else {
            // Single term logic (treat as 1 * expr)
            factors = [expr];
        }

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
        let candidates = factors.map((f, i) => ({
            idx: i,
            node: f,
            score: typeScore(f)
        })).filter(c => this._dependsOn(c.node, varNode)); // Only variable parts matter for u choice

        // If single term, we force u=expr, dv=1 IF score is high (Log/InvTrig)
        if (factors.length === 1) {
             const s = typeScore(expr);
             if (s >= 4) {
                 candidates = [{idx: 0, node: expr, score: s}];
             } else {
                 return new Call("integrate", [expr, varNode]);
             }
        }

        // Sort by score descending (Log first)
        candidates.sort((a, b) => b.score - a.score);

        for(const cand of candidates) {
            const u = cand.node;

            // dv is the product of all other factors
            let dv = null;
            if (factors.length === 1) {
                dv = new Num(1);
            } else {
                for(let j=0; j<factors.length; j++) {
                    if (cand.idx === j) continue;
                    if (dv === null) dv = factors[j];
                    else dv = new Mul(dv, factors[j]);
                }
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

            const intVdu = this.evaluate(new Call("integrate", [vdu, varNode]));

            // Check if we made progress or just got back a Call('integrate')
            // If intVdu is exactly Call('integrate', [vdu]), then we failed to solve the sub-problem.
            // But we should still return the partial result `uv - int(vdu)` because the user might prefer that form
            // over the original integral.

            return new Sub(uv, intVdu).simplify();
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

    _normalCDF(x, mu, sigma) {
        // 0.5 * (1 + erf((x-mu)/(sigma*sqrt(2))))
        const arg = new Div(new Sub(x, mu), new Mul(sigma, new Call('sqrt', [new Num(2)]))).simplify();
        const erfTerm = new Call('erf', [arg]);
        return new Mul(new Num(0.5), new Add(new Num(1), erfTerm)).simplify();
    }

    _invNorm(area, mu, sigma) {
        area = area.simplify();
        mu = mu.simplify();
        sigma = sigma.simplify();

        if (area instanceof Num && mu instanceof Num && sigma instanceof Num) {
            const p = area.value;
            const m = mu.value;
            const s = sigma.value;

            if (p <= 0 || p >= 1) return new Sym("NaN");
            if (s < 0) return new Sym("NaN");

            // Rational approximation for standard normal quantile
            // Abramowitz and Stegun 26.2.23
            const c0 = 2.515517;
            const c1 = 0.802853;
            const c2 = 0.010328;
            const d1 = 1.432788;
            const d2 = 0.189269;
            const d3 = 0.001308;

            let t, num, den, xp;

            let q = (p < 0.5) ? p : 1.0 - p;
            t = Math.sqrt(-2.0 * Math.log(q));
            num = c0 + c1 * t + c2 * t * t;
            den = 1.0 + d1 * t + d2 * t * t + d3 * t * t * t;

            // Abramowitz and Stegun 26.2.23 gives x_p = t - num/den for upper tail Q(x_p) = q
            xp = t - num / den;

            if (p < 0.5) xp = -xp;

            return new Num(m + s * xp);
        }
        return new Call('invNorm', [area, mu, sigma]);
    }

    _binomialCDF(k, n, p) {
        k = k.simplify();
        n = n.simplify();
        // sum binomialPDF from 0 to floor(k)
        if (k instanceof Num) {
             const limit = Math.floor(k.value);
             let sum = new Num(0);
             // Should we use symbolic sum if n/p are symbolic?
             // Or explicitly construct sum(pdf, i, 0, k)
             // If k is symbolic, return Sum
             if (limit < 0) return new Num(0);

             // If n is number, we can loop. If n is symbol, use symbolic Sum
             if (n instanceof Num) {
                 for(let i=0; i<=limit; i++) {
                     sum = new Add(sum, this._binomialPDF(new Num(i), n, p)).simplify();
                 }
                 return sum;
             }
        }

        // Symbolic fallback
        const i = new Sym('__i_cdf__');
        return new Call('sum', [
             new Call('binomialPDF', [i, n, p]),
             i,
             new Num(0),
             new Call('floor', [k])
        ]);
    }

    _poissonPDF(k, lambda) {
        k = k.simplify();
        lambda = lambda.simplify();

        if (k instanceof Num && lambda instanceof Num) {
            const kv = k.value;
            const lv = lambda.value;
            if (!Number.isInteger(kv) || kv < 0) return new Num(0);
            if (lv < 0) return new Num(0);

            // e^-lambda * lambda^k / k!
            const prob = Math.exp(-lv) * Math.pow(lv, kv) / this._factorial(kv);
            return new Num(prob);
        }
        // Symbolic: e^-L * L^k / k!
        return new Div(
            new Mul(new Call('exp', [new Mul(new Num(-1), lambda)]), new Pow(lambda, k)),
            new Call('factorial', [k])
        ).simplify();
    }

    _poissonCDF(k, lambda) {
        k = k.simplify();
        lambda = lambda.simplify();

        if (k instanceof Num) {
            const limit = Math.floor(k.value);
            if (limit < 0) return new Num(0);

            if (lambda instanceof Num) {
                let sum = new Num(0);
                for(let i=0; i<=limit; i++) {
                    sum = new Add(sum, this._poissonPDF(new Num(i), lambda)).simplify();
                }
                return sum;
            }
        }

        const i = new Sym('__i_pcdf__');
        return new Call('sum', [
            new Call('poissonPDF', [i, lambda]),
            i,
            new Num(0),
            new Call('floor', [k])
        ]);
    }

    _exponentialPDF(x, lambda) {
        x = x.simplify();
        lambda = lambda.simplify();

        if (x instanceof Num && lambda instanceof Num) {
            const xv = x.value;
            const lv = lambda.value;
            if (xv < 0) return new Num(0);
            if (lv <= 0) return new Num(0);
            return new Num(lv * Math.exp(-lv * xv));
        }
        // lambda * e^(-lambda * x)
        return new Mul(lambda, new Call('exp', [new Mul(new Num(-1), new Mul(lambda, x))])).simplify();
    }

    _exponentialCDF(x, lambda) {
        x = x.simplify();
        lambda = lambda.simplify();

        if (x instanceof Num && lambda instanceof Num) {
            const xv = x.value;
            const lv = lambda.value;
            if (xv < 0) return new Num(0);
            if (lv <= 0) return new Num(0); // Undefined?
            return new Num(1 - Math.exp(-lv * xv));
        }
        // 1 - e^(-lambda * x)
        const expTerm = new Call('exp', [new Mul(new Num(-1), new Mul(lambda, x))]);
        return new Sub(new Num(1), expTerm).simplify();
    }

    _geometricPDF(k, p) {
        // P(X=k) = p(1-p)^(k-1) for k=1,2,3...
        k = k.simplify();
        p = p.simplify();

        if (k instanceof Num && p instanceof Num) {
            const kv = k.value;
            const pv = p.value;
            if (!Number.isInteger(kv) || kv < 1) return new Num(0);
            if (pv <= 0 || pv > 1) return new Num(0); // p must be (0,1]
            return new Num(pv * Math.pow(1 - pv, kv - 1));
        }
        // p * (1-p)^(k-1)
        return new Mul(p, new Pow(new Sub(new Num(1), p), new Sub(k, new Num(1)))).simplify();
    }

    _geometricCDF(k, p) {
        // P(X<=k) = 1 - (1-p)^k for k>=1
        k = k.simplify();
        p = p.simplify();

        if (k instanceof Num && p instanceof Num) {
            const kv = Math.floor(k.value); // Use floor for discrete CDF
            const pv = p.value;
            if (kv < 1) return new Num(0);
            if (pv <= 0 || pv > 1) return new Num(0);
            return new Num(1 - Math.pow(1 - pv, kv));
        }
        // 1 - (1-p)^k
        return new Sub(new Num(1), new Pow(new Sub(new Num(1), p), k)).simplify();
    }

    _chisquarePDF(x, k) {
        // f(x; k) = (1 / (2^(k/2) * Gamma(k/2))) * x^(k/2 - 1) * e^(-x/2)
        // for x > 0
        x = x.simplify();
        k = k.simplify();

        if (x instanceof Num && k instanceof Num) {
            const xv = x.value;
            const kv = k.value;
            if (xv <= 0) return new Num(0);
            if (kv <= 0) return new Num(0);

            const num = Math.pow(xv, kv/2 - 1) * Math.exp(-xv/2);
            const den = Math.pow(2, kv/2) * this._gamma(kv/2);
            return new Num(num / den);
        }

        // Symbolic
        const kOver2 = new Div(k, new Num(2));
        const gammaTerm = new Call('gamma', [kOver2]);
        const den = new Mul(new Pow(new Num(2), kOver2), gammaTerm);
        const term1 = new Div(new Num(1), den);
        const term2 = new Pow(x, new Sub(kOver2, new Num(1)));
        const term3 = new Call('exp', [new Mul(new Num(-0.5), x)]);

        return new Mul(term1, new Mul(term2, term3)).simplify();
    }

    _compound(P, r, n, t) {
        // A = P * (1 + r/n)^(n*t)
        P = P.simplify();
        r = r.simplify();
        n = n.simplify();
        t = t.simplify();

        const base = new Add(new Num(1), new Div(r, n));
        const exponent = new Mul(n, t);
        return new Mul(P, new Pow(base, exponent)).simplify();
    }

    _loan(P, r, n) {
        // Monthly payment PMT = (P * r/12) / (1 - (1 + r/12)^(-n))
        // Assuming r is annual rate (e.g. 0.05), n is total months
        // If args are just P, r, n, we follow standard formula

        P = P.simplify();
        r = r.simplify();
        n = n.simplify();

        const monthlyRate = new Div(r, new Num(12));
        const num = new Mul(P, monthlyRate);
        const base = new Add(new Num(1), monthlyRate);
        const denom = new Sub(new Num(1), new Pow(base, new Mul(new Num(-1), n)));

        return new Div(num, denom).simplify();
    }

    _npv(rate, flows) {
        if (!(flows instanceof Vec)) throw new Error("Cash flows must be a list");
        rate = rate.simplify();
        // rate is typically decimal (0.05) or symbolic
        let sum = new Num(0);
        // flows[0] is at t=0, flows[1] at t=1, etc.
        for(let i=0; i<flows.elements.length; i++) {
            const flow = flows.elements[i];
            const den = new Pow(new Add(new Num(1), rate), new Num(i)).simplify();
            sum = new Add(sum, new Div(flow, den)).simplify();
        }
        return sum;
    }

    _irr(flows) {
        if (!(flows instanceof Vec)) throw new Error("Cash flows must be a list");

        // Solve NPV(r) = 0 for r > -1
        // NPV = C0 + C1/(1+r) + C2/(1+r)^2 + ...
        // Let x = 1/(1+r). Then r = 1/x - 1.
        // Polynomial: C0 + C1*x + C2*x^2 + ... = 0

        const x = new Sym('__x_irr__');
        let poly = new Num(0);

        for(let i=0; i<flows.elements.length; i++) {
            const flow = flows.elements[i];
            const term = new Mul(flow, new Pow(x, new Num(i))).simplify();
            poly = new Add(poly, term).simplify();
        }

        // Solve for x
        const roots = this._solve(poly, x);

        const solutions = [];
        const processRoot = (root) => {
             // r = 1/root - 1
             const rVal = new Sub(new Div(new Num(1), root), new Num(1)).simplify();
             // Check if real and > -1 (x must be positive for r > -1)
             // Actually if x > 0, then 1+r > 0 => r > -1.
             // If x < 0, 1+r < 0 => r < -1 (not useful usually).
             // If x = 0, impossible (1/(1+r)=0 => r->inf)

             try {
                 const numX = root.evaluateNumeric();
                 if (!isNaN(numX) && numX > 0) {
                     solutions.push(rVal);
                 } else if (isNaN(numX)) {
                     // Keep symbolic
                     solutions.push(rVal);
                 }
             } catch(e) {
                 solutions.push(rVal);
             }
        };

        if (roots instanceof Call && roots.funcName === 'set') {
            roots.args.forEach(processRoot);
        } else if (roots instanceof Expr && !(roots instanceof Call && roots.funcName === 'solve')) {
            processRoot(roots);
        } else if (roots instanceof Vec) {
             // Systems return Vec? No, _solve returns set or Expr.
        }

        if (solutions.length === 0) return new Call('irr', [flows]);
        if (solutions.length === 1) return solutions[0];
        return new Call('set', solutions);
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

        if (a instanceof Num && b instanceof Num) {
             return new Num(a.value % b.value);
        }
        return new Call('rem', [a, b]);
    }

    _quo(a, b) {
        a = a.simplify();
        b = b.simplify();

        if (a instanceof Num && b instanceof Num) {
             return new Num(Math.trunc(a.value / b.value));
        }
        return new Call('quo', [a, b]);
    }

    _mod(a, b) {
        a = a.simplify();
        b = b.simplify();

        if (a instanceof Num && b instanceof Num) {
             // Mathematical modulo (always positive for positive modulus)
             const m = b.value;
             return new Num(((a.value % m) + m) % m);
        }
        return new Call('mod', [a, b]);
    }

    _distance(p1, p2) {
        if (!(p1 instanceof Vec) || !(p2 instanceof Vec)) throw new Error("distance requires two points (vectors)");
        if (p1.elements.length !== p2.elements.length) throw new Error("Dimension mismatch");
        let sum = new Num(0);
        for(let i=0; i<p1.elements.length; i++) {
            const diff = new Sub(p1.elements[i], p2.elements[i]);
            sum = new Add(sum, new Pow(diff, new Num(2)));
        }
        return new Call('sqrt', [sum.simplify()]).simplify();
    }

    _midpoint(p1, p2) {
        if (!(p1 instanceof Vec) || !(p2 instanceof Vec)) throw new Error("midpoint requires two points (vectors)");
        if (p1.elements.length !== p2.elements.length) throw new Error("Dimension mismatch");
        const elements = [];
        for(let i=0; i<p1.elements.length; i++) {
            const sum = new Add(p1.elements[i], p2.elements[i]);
            elements.push(new Div(sum, new Num(2)).simplify());
        }
        return new Vec(elements);
    }

    _zTest(data, mu0, sigma) {
        // One-sample Z-test
        // H0: mu = mu0
        // z = (x_bar - mu0) / (sigma / sqrt(n))

        let x_bar, n;

        if (data instanceof Vec) {
            n = new Num(data.elements.length);
            x_bar = this._mean(data);
        } else {
             throw new Error("zTest requires a data list as first argument");
        }

        mu0 = mu0.simplify();
        sigma = sigma.simplify();

        const num = new Sub(x_bar, mu0);
        const den = new Div(sigma, new Call('sqrt', [n]));
        const z = new Div(num, den).simplify();

        // p-value (two-tailed): 2 * (1 - CDF(|z|))
        // CDF(z) for standard normal is normalCDF(z, 0, 1)
        const absZ = new Call('abs', [z]).simplify();
        const cdf = this._normalCDF(absZ, new Num(0), new Num(1));
        const pValue = new Mul(new Num(2), new Sub(new Num(1), cdf)).simplify();

        return new Vec([z, pValue]);
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

    _euler(n) {
        n = n.simplify();
        if (n instanceof Num && Number.isInteger(n.value) && n.value > 0) {
            let val = n.value;
            let result = val;
            let p = 2;
            while (p * p <= val) {
                if (val % p === 0) {
                    while (val % p === 0) val /= p;
                    result -= result / p;
                }
                p++;
            }
            if (val > 1) result -= result / val;
            return new Num(result);
        }
        return new Call('euler', [n]);
    }

    _mode(list) {
        if (list instanceof Vec) {
            if (list.elements.length === 0) return new Vec([]); // Or null?

            const counts = {};
            let maxCount = 0;

            for(const el of list.elements) {
                const s = el.toString(); // Use string representation for counting
                counts[s] = (counts[s] || 0) + 1;
                if (counts[s] > maxCount) maxCount = counts[s];
            }

            if (maxCount === 1) return list; // No mode (or all are modes)

            const modes = [];
            // Find elements with maxCount
            // We need to map back to Expr. Since we keyed by toString(), we need original mapping or re-parse.
            // But we can iterate original list.
            const seen = new Set();
            for(const el of list.elements) {
                const s = el.toString();
                if (counts[s] === maxCount && !seen.has(s)) {
                    modes.push(el);
                    seen.add(s);
                }
            }

            if (modes.length === 1) return modes[0];
            return new Vec(modes);
        }
        return new Call('mode', [list]);
    }

    _arcLen(expr, varNode, start, end) {
        if (!(varNode instanceof Sym)) throw new Error("Variable must be a symbol");

        // integrate(sqrt(1 + f'(x)^2), x, a, b)
        const deriv = expr.diff(varNode).simplify();
        const integrand = new Call('sqrt', [new Add(new Num(1), new Pow(deriv, new Num(2)))]).simplify();

        return this.evaluate(new Call('integrate', [integrand, varNode, start, end]));
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

            let isZero = false;
            if (u instanceof Num && u.value === 0) isZero = true;
            else if (u instanceof Vec) isZero = u.elements.every(e => e instanceof Num && e.value === 0);

            if (!isZero && u instanceof Vec) {
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
        // Can be diff(y, x, 2) OR diff(diff(y, x), x)
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
                // Handle nested diff(diff(y, x), x)
                if (node.args.length === 2) {
                    const inner = node.args[0];
                    const outerVar = node.args[1];
                    if (inner instanceof Call && inner.funcName === 'diff' && inner.args.length === 2) {
                        const target = inner.args[0];
                        const innerVar = inner.args[1];
                        let isMatch = false;
                        if ((target instanceof Sym && target.name === depVar.name) ||
                            (target instanceof Call && target.funcName === depVar.name)) {
                            isMatch = true;
                        }
                        if (isMatch) {
                            if (!indepVar) indepVar = outerVar;
                            if (outerVar.name === indepVar.name && innerVar.name === indepVar.name) {
                                diffCall2 = node;
                            }
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
                 // nested diff(diff(y, x), x) - handled by ref equality check mostly, but generic struct check:
                 if (node instanceof Call && node.funcName === 'diff' && node.args.length === 2 && node.args[0] instanceof Call && node.args[0].funcName === 'diff') {
                      // assume correct if we found it earlier
                      if (indepVar && node.args[1].name === indepVar.name) return D2Y;
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
                 // Calculate remainder and force clean derivatives (they should cancel mathematically)
                 let remainder = new Sub(polyForm, recon).expand().simplify();
                 remainder = remainder.substitute(D2Y, new Num(0)).substitute(DY, new Num(0)).substitute(y_sym, new Num(0)).simplify();

                 // Check for Constant Coefficients (independent of x)
                 if (!this._dependsOn(a, indepVar) && !this._dependsOn(b, indepVar) && !this._dependsOn(c, indepVar)) {
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
                     let y_c;

                     const diffR = new Sub(r1, r2).simplify();
                     // Check for repeated roots
                     if (diffR instanceof Num && Math.abs(diffR.value) < 1e-9) {
                         const term1 = new Mul(C1, new Call('exp', [new Mul(r1, indepVar)])).simplify();
                         const term2 = new Mul(C2, new Mul(indepVar, new Call('exp', [new Mul(r1, indepVar)]))).simplify();
                         y_c = new Add(term1, term2);
                     } else {
                         const term1 = new Mul(C1, new Call('exp', [new Mul(r1, indepVar)])).simplify();
                         const term2 = new Mul(C2, new Call('exp', [new Mul(r2, indepVar)])).simplify();
                         y_c = new Add(term1, term2);
                     }

                     // Check Non-Homogeneous
                     const rhs = new Mul(new Num(-1), remainder).simplify();
                     if (rhs instanceof Num && rhs.value === 0) {
                         return y_c;
                     }

                     // Method of Undetermined Coefficients (Basic cases)
                     let y_p = null;

                     // Case 1: Polynomial / Constant
                     // Check if rhs is polynomial in indepVar
                     const rhsPoly = this._getPolyCoeffs(rhs, indepVar);
                     if (rhsPoly) {
                         // Simplify: only handle constant or linear RHS for stability
                         // Guess A*x + B
                         if (rhsPoly.maxDeg <= 1) {
                             const A = new Sym('__A__');
                             const B = new Sym('__B__');
                             let guess = new Add(new Mul(A, indepVar), B);

                             // Resonance check
                             if (c instanceof Num && c.value === 0) {
                                 guess = new Mul(indepVar, guess);
                                 if (b instanceof Num && b.value === 0) {
                                     guess = new Mul(indepVar, guess);
                                 }
                             }

                             // Substitute into ODE: a*y'' + b*y' + c*y - rhs = 0
                             const yp_prime = guess.diff(indepVar).simplify();
                             const yp_double = yp_prime.diff(indepVar).simplify();

                             const LHS = new Add(new Add(new Mul(a, yp_double), new Mul(b, yp_prime)), new Mul(c, guess)).simplify();
                             const eqSys = new Sub(LHS, rhs).simplify();

                             // Solve for A and B by comparing coefficients
                             // eqSys should be 0.
                             // Extract coeffs of indepVar
                             const sysPoly = this._getPolyCoeffs(eqSys, indepVar);
                             if (sysPoly) {
                                 const eqs = [];
                                 const vars = [A, B];
                                 // We need to solve for A and B. We have coeffs for x^0, x^1, ...
                                 for(let d in sysPoly.coeffs) {
                                     eqs.push(new Eq(sysPoly.coeffs[d], new Num(0)));
                                 }
                                 // Solve system
                                 const sol = this._solveSystem(new Vec(eqs), new Vec(vars));
                                 if (sol instanceof Vec && sol.elements.length === 2) {
                                     const valA = sol.elements[0];
                                     const valB = sol.elements[1];
                                     y_p = guess.substitute(A, valA).substitute(B, valB).simplify();
                                 }
                             }
                         }
                     }

                     // Case 2: Exponential C * e^(k*x)
                     let expNode = rhs;
                     // Handle Coeff * exp
                     if (rhs instanceof Mul && rhs.left instanceof Num) expNode = rhs.right;

                     // console.log("DEBUG Exp", expNode.toString(), expNode.constructor.name);

                     let uExp = null;
                     if (expNode instanceof Call && expNode.funcName === 'exp') uExp = expNode.args[0];
                     else if (expNode instanceof Pow) {
                         if (expNode.left instanceof Sym && expNode.left.name === 'e') uExp = expNode.right;
                         // Handle evaluated 'e' (Num)
                         if (expNode.left instanceof Num && Math.abs(expNode.left.value - Math.E) < 1e-9) uExp = expNode.right;
                     }

                     if (!y_p && uExp) {
                         // rhs = exp(u). u linear in x?
                         const uPoly = this._getPolyCoeffs(uExp, indepVar);
                         if (uPoly && uPoly.maxDeg === 1 && (!uPoly.coeffs[0] || uPoly.coeffs[0].value === 0)) {
                             // exp(k*x)
                             const k = uPoly.coeffs[1];
                             const A = new Sym('__A__');
                             let guess = new Mul(A, expNode); // Guess A * e^kx (ignoring C in rhs, A will absorb it)

                             // Resonance: is k a root?
                             // Check a*k^2 + b*k + c = 0
                             const charVal = new Add(new Add(new Mul(a, new Pow(k, new Num(2))), new Mul(b, k)), c).simplify();
                             if (charVal instanceof Num && charVal.value === 0) {
                                 guess = new Mul(indepVar, guess);
                                 // Double root check
                                 const charDeriv = new Add(new Mul(new Num(2), new Mul(a, k)), b).simplify();
                                 if (charDeriv instanceof Num && charDeriv.value === 0) {
                                     guess = new Mul(indepVar, guess);
                                 }
                             }

                             // Substitute
                             const yp_prime = guess.diff(indepVar).simplify();
                             const yp_double = yp_prime.diff(indepVar).simplify();
                             const LHS = new Add(new Add(new Mul(a, yp_double), new Mul(b, yp_prime)), new Mul(c, guess)).simplify();
                             const eqSys = new Sub(LHS, rhs).simplify();

                             // Factor out e^(kx) or divide by it?
                             // Dividing might be safer
                             const reduced = new Div(eqSys, expNode).simplify(); // Should be linear in A
                             const valA = this._solve(reduced, A);
                             if (!(valA instanceof Call && valA.funcName === 'solve')) {
                                 y_p = guess.substitute(A, valA).simplify();
                             }
                         }
                     }

                     if (y_p) {
                         return new Add(y_c, y_p).simplify();
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

    _linearizeTrig(expr) {
        if (expr instanceof Add) return new Add(this._linearizeTrig(expr.left), this._linearizeTrig(expr.right)).simplify();
        if (expr instanceof Sub) return new Sub(this._linearizeTrig(expr.left), this._linearizeTrig(expr.right)).simplify();
        if (expr instanceof Div) return new Div(this._linearizeTrig(expr.left), this._linearizeTrig(expr.right)).simplify();

        if (expr instanceof Pow) {
            // Handle sin^2, cos^2
            if (expr.right instanceof Num && expr.right.value === 2) {
                const base = expr.left;
                if (base instanceof Call) {
                    if (base.funcName === 'sin') {
                        // (1 - cos(2x))/2
                        const arg = base.args[0];
                        return new Div(new Sub(new Num(1), new Call('cos', [new Mul(new Num(2), arg)])), new Num(2)).simplify();
                    }
                    if (base.funcName === 'cos') {
                        // (1 + cos(2x))/2
                        const arg = base.args[0];
                        return new Div(new Add(new Num(1), new Call('cos', [new Mul(new Num(2), arg)])), new Num(2)).simplify();
                    }
                }
            }
            return new Pow(this._linearizeTrig(expr.left), this._linearizeTrig(expr.right));
        }

        if (expr instanceof Mul) {
            // Flatten first
            const terms = [];
            const collect = (n) => {
                if (n instanceof Mul) { collect(n.left); collect(n.right); }
                else terms.push(this._linearizeTrig(n));
            };
            collect(expr);

            // Separate trig terms and others
            const trigTerms = [];
            const otherTerms = [];

            for(const t of terms) {
                if (t instanceof Call && (t.funcName === 'sin' || t.funcName === 'cos')) {
                    trigTerms.push(t);
                } else {
                    otherTerms.push(t);
                }
            }

            if (trigTerms.length >= 2) {
                // Combine first two
                const t1 = trigTerms.pop();
                const t2 = trigTerms.pop();

                const A = t1.args[0];
                const B = t2.args[0];

                let combined = null;

                const half = new Div(new Num(1), new Num(2));
                const sum = new Add(A, B);
                const diff = new Sub(A, B);

                if (t1.funcName === 'sin' && t2.funcName === 'cos') {
                    // sin(A)cos(B)
                    combined = new Mul(half, new Add(new Call('sin', [sum]), new Call('sin', [diff])));
                } else if (t1.funcName === 'cos' && t2.funcName === 'sin') {
                    // cos(A)sin(B)
                    combined = new Mul(half, new Add(new Call('sin', [sum]), new Call('sin', [new Sub(B, A)])));
                } else if (t1.funcName === 'sin' && t2.funcName === 'sin') {
                    // sin(A)sin(B)
                    combined = new Mul(half, new Sub(new Call('cos', [diff]), new Call('cos', [sum])));
                } else if (t1.funcName === 'cos' && t2.funcName === 'cos') {
                    // cos(A)cos(B)
                    combined = new Mul(half, new Add(new Call('cos', [diff]), new Call('cos', [sum])));
                }

                // Push back combined
                trigTerms.push(this._linearizeTrig(combined));

                // Rebuild Mul
                let res = trigTerms[0];
                for(let i=1; i<trigTerms.length; i++) res = new Mul(res, trigTerms[i]);
                for(const t of otherTerms) res = new Mul(res, t);

                return res.simplify();
            }

            // No pairs found
            let res = terms[0];
            for(let i=1; i<terms.length; i++) res = new Mul(res, terms[i]);
            return res;
        }

        if (expr instanceof Call) {
            return new Call(expr.funcName, expr.args.map(a => this._linearizeTrig(a)));
        }

        return expr;
    }

    _fourier(expr, varNode, n, L) {
        if (!L) L = new Sym('pi');

        // Coeffs:
        // a0 = 1/L * int(f, -L, L)
        // an = 1/L * int(f*cos(k*pi*x/L), -L, L)
        // bn = 1/L * int(f*sin(k*pi*x/L), -L, L)

        // Evaluate L numeric if possible for integration limits
        // But keep symbolic L for integrand if needed?
        // Actually, integration works best with matching symbols.

        let pi = new Sym('pi');
        // If L is numeric and close to PI, use numeric PI to allow cancellation
        if (L instanceof Num && Math.abs(L.value - Math.PI) < 1e-5) {
            pi = new Num(Math.PI);
        }

        const factor = new Div(new Num(1), L).simplify();
        const negL = new Mul(new Num(-1), L).simplify();

        // a0
        const intA0 = this.evaluate(new Call('integrate', [expr, varNode, negL, L]));
        const a0 = new Mul(factor, intA0).simplify();

        // Term a0/2
        let result = new Div(a0, new Num(2)).simplify();

        const N = (n instanceof Num) ? n.value : 5;

        for(let k=1; k<=N; k++) {
            const kNum = new Num(k);
            const angle = new Div(new Mul(kNum, new Mul(pi, varNode)), L).simplify();

            // an
            const cosTerm = new Call('cos', [angle]);
            let integrandA = new Mul(expr, cosTerm).simplify();
            integrandA = this._linearizeTrig(integrandA);

            const intAn = this.evaluate(new Call('integrate', [integrandA, varNode, negL, L]));
            const an = new Mul(factor, intAn).simplify();

            // bn
            const sinTerm = new Call('sin', [angle]);
            let integrandB = new Mul(expr, sinTerm).simplify();
            integrandB = this._linearizeTrig(integrandB);

            const intBn = this.evaluate(new Call('integrate', [integrandB, varNode, negL, L]));
            const bn = new Mul(factor, intBn).simplify();

            // Add to result
            const termA = new Mul(an, cosTerm).simplify();
            const termB = new Mul(bn, sinTerm).simplify();

            result = new Add(result, new Add(termA, termB)).simplify();
        }

        return result;
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

            // Properties
            // 1. Frequency Shift: L{e^(at)*f(t)} = F(s-a)
            const checkExp = (term, func) => {
                if (term instanceof Call && term.funcName === 'exp') {
                    const arg = term.args[0];
                    // arg = a*t. Extract a.
                    const poly = this._getPolyCoeffs(arg, t);
                    if (poly && poly.maxDeg === 1 && (!poly.coeffs[0] || poly.coeffs[0].value === 0)) {
                        const a = poly.coeffs[1];
                        // Transform func(t) -> F(s)
                        const F = this._laplace(func, t, s);
                        if (!(F instanceof Call && F.funcName === 'laplace')) {
                            // Substitute s -> s-a
                            return F.substitute(s, new Sub(s, a)).simplify();
                        }
                    }
                }
                return null;
            };

            const shift1 = checkExp(expr.left, expr.right);
            if (shift1) return shift1;
            const shift2 = checkExp(expr.right, expr.left);
            if (shift2) return shift2;

            // 2. Mult by t^n: L{t^n * f(t)} = (-1)^n * d^n/ds^n F(s)
            const checkT = (term, func) => {
                let n = 0;
                if (term instanceof Sym && term.name === t.name) n = 1;
                else if (term instanceof Pow && term.left instanceof Sym && term.left.name === t.name && term.right instanceof Num && Number.isInteger(term.right.value)) {
                    n = term.right.value;
                }

                if (n > 0) {
                    const F = this._laplace(func, t, s);
                    if (!(F instanceof Call && F.funcName === 'laplace')) {
                        // diff(F, s, n)
                        let dF = F;
                        for(let i=0; i<n; i++) dF = dF.diff(s).simplify();
                        const sign = (n % 2 === 0) ? new Num(1) : new Num(-1);
                        return new Mul(sign, dF).simplify();
                    }
                }
                return null;
            };

            const multT1 = checkT(expr.left, expr.right);
            if (multT1) return multT1;
            const multT2 = checkT(expr.right, expr.left);
            if (multT2) return multT2;
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
             if (poly && poly.maxDeg === 1 && (!poly.coeffs[0] || poly.coeffs[0].value === 0)) {
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
             // handle a*t
             const poly = this._getPolyCoeffs(arg, t);
             if (poly && poly.maxDeg === 1 && (!poly.coeffs[0] || poly.coeffs[0].value === 0)) {
                 a = poly.coeffs[1];
             } else {
                 return new Call('laplace', [expr, t, s]);
             }

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

             // 1/(s-a) or 1/(s+a)
             // Check if den is linear in s
             const polyDen = this._getPolyCoeffs(den, s);
             if (polyDen && polyDen.maxDeg === 1) {
                 // c1*s + c0
                 // = c1 * (s + c0/c1)
                 const c1 = polyDen.coeffs[1];
                 const c0 = polyDen.coeffs[0] || new Num(0);
                 const a = new Div(c0, c1).simplify(); // a is negative of pole if form is s-p
                 // form: k / (s+a) -> k * e^(-at)
                 const k = new Div(num, c1).simplify();
                 return new Mul(k, new Call('exp', [new Mul(new Num(-1), new Mul(a, t))])).simplify();
             }

             // s^2 + a^2 (Quadratic Denominator)
             if (polyDen && polyDen.maxDeg === 2) {
                 const c2 = polyDen.coeffs[2];
                 const c1 = polyDen.coeffs[1] || new Num(0);
                 const c0 = polyDen.coeffs[0] || new Num(0);

                 // Check if c1 = 0 (Pure quadratic s^2 + k)
                 if (c1 instanceof Num && c1.value === 0) {
                     // c2*s^2 + c0
                     // = c2 * (s^2 + c0/c2)
                     const kSq = new Div(c0, c2).simplify();
                     // Check sign of kSq
                     // If positive -> sin/cos
                     // If negative -> sinh/cosh (or exp)
                     // Let's assume kSq > 0 for now or treat symbolically
                     const k = new Call('sqrt', [kSq]).simplify();

                     // Numerator: As + B
                     const polyNum = this._getPolyCoeffs(num, s);
                     let A = new Num(0);
                     let B = num;
                     if (polyNum && polyNum.maxDeg <= 1) {
                         A = polyNum.coeffs[1] || new Num(0);
                         B = polyNum.coeffs[0] || new Num(0);
                     }

                     // Form: (As + B) / (c2(s^2 + k^2))
                     // = (A/c2) * s/(s^2+k^2) + (B/c2) * 1/(s^2+k^2)
                     // = (A/c2) * cos(kt) + (B/(c2*k)) * sin(kt)

                     const term1 = new Mul(new Div(A, c2), new Call('cos', [new Mul(k, t)]));
                     const term2 = new Mul(new Div(B, new Mul(c2, k)), new Call('sin', [new Mul(k, t)]));

                     return new Add(term1, term2).simplify();
                 }
             }
        }

        // Try Partial Fraction if Div
        if (expr instanceof Div) {
             const pf = this._partfrac(expr, s);
             if (!(pf instanceof Call && pf.funcName === 'partfrac') && pf.toString() !== expr.toString()) {
                 return this._ilaplace(pf, s, t);
             }
        }

        return new Call('ilaplace', [expr, s, t]);
    }

    _molarMass(formula) {
        if (typeof globalThis.Chemistry !== 'undefined') {
            try {
                const mass = globalThis.Chemistry.calculateMolarMass(formula);
                return new Num(mass);
            } catch(e) {
                throw new Error("Chemistry Error: " + e.message);
            }
        }
        throw new Error("Chemistry module not loaded");
    }

    _atomicWeight(sym) {
        if (typeof globalThis.Chemistry !== 'undefined') {
            const el = globalThis.Chemistry.getElement(sym);
            if (el) return new Num(el.mass);
            throw new Error("Unknown element");
        }
        throw new Error("Chemistry module not loaded");
    }

    _balanceChem(expr) {
        if (typeof globalThis.Chemistry === 'undefined') throw new Error("Chemistry module not loaded");

        // Convert AST back to string to parse as chemical equation
        // Or handle AST directly.
        // H2 + O2 -> H2O
        // AST: Implies(Add(H2, O2), H2O)
        // Note: H2 is a Sym "H2".
        // 2H2O might be Mul(2, H2O) or Sym "2H2O"? Lexer splits number. Mul(2, H2O).
        // But for balancing, we input raw formulas without coefficients usually.
        // User types: H2 + O2 -> H2O.
        // AST: Implies( Add(Sym(H2), Sym(O2)), Sym(H2O) )

        let lhsNodes = [];
        let rhsNodes = [];

        // Helper to flatten Add
        const flatten = (node, list) => {
            if (node instanceof Add) {
                flatten(node.left, list);
                flatten(node.right, list);
            } else {
                list.push(node);
            }
        };

        if (expr instanceof Implies) { // ->
            flatten(expr.left, lhsNodes);
            flatten(expr.right, rhsNodes);
        } else if (expr instanceof Eq) { // =
            flatten(expr.left, lhsNodes);
            flatten(expr.right, rhsNodes);
        } else {
            // Maybe just one side or error
            throw new Error("Invalid chemical equation format. Use -> or =");
        }

        const compounds = []; // { name: "H2O", side: -1 (LHS) or 1 (RHS), counts: {} }
        const allElements = new Set();

        lhsNodes.forEach(node => {
            const name = node.toString();
            const counts = globalThis.Chemistry.parseMolecule(name);
            Object.keys(counts).forEach(e => allElements.add(e));
            compounds.push({ name, side: -1, counts });
        });

        rhsNodes.forEach(node => {
            const name = node.toString();
            const counts = globalThis.Chemistry.parseMolecule(name);
            Object.keys(counts).forEach(e => allElements.add(e));
            compounds.push({ name, side: 1, counts });
        });

        const elementList = Array.from(allElements).sort();
        const m = elementList.length; // rows (equations)
        const n = compounds.length;   // cols (variables)

        // Build Matrix A for Ax = 0
        // col j is compound j.
        // row i is element i.
        // A[i][j] = side * count
        // We want sum(coef * side * count) = 0.
        // Let coef vector be x.
        // A[i][j] = compounds[j].side * compounds[j].counts[element[i]]

        const mat = [];
        for(let i=0; i<m; i++) {
            const row = [];
            const el = elementList[i];
            for(let j=0; j<n; j++) {
                const c = compounds[j];
                const count = c.counts[el] || 0;
                // LHS should use positive coeffs in reaction, RHS positive.
                // sum(lhs) = sum(rhs) => sum(lhs) - sum(rhs) = 0.
                // LHS coeffs * count = RHS coeffs * count
                // LHS * count - RHS * count = 0
                // So LHS side factor +1, RHS side factor -1?
                // Wait, if reaction is aA + bB -> cC
                // element H: a*H_A + b*H_B = c*H_C
                // a*H_A + b*H_B - c*H_C = 0
                // So if we solve for vector [a, b, c],
                // LHS cols should be positive count, RHS cols negative count.
                const val = (c.side === -1 ? 1 : -1) * count;
                row.push(new Num(val));
            }
            mat.push(new Vec(row));
        }

        // Solve nullspace
        const A = new Vec(mat);
        const kernel = this._kernel(A); // Returns basis vectors

        if (kernel.elements.length === 0) throw new Error("Cannot balance equation (no solution)");

        // Pick the first basis vector
        // It might have fractions or zeros.
        let sol = kernel.elements[0]; // Vec

        // Normalize to smallest integers
        // 1. Make all positive (coeffs must be positive)
        // If first non-zero is negative, negate all?
        // Physical coeffs are positive.
        // Check if we need to negate.
        // Usually kernel returns canonical basis.
        // Let's assume we can scale.

        let vec = sol.elements.map(e => e.evaluateNumeric());

        // Check for any negative values
        // If mixed signs, it might be impossible or we picked wrong basis combo?
        // For simple reactions, kernel dim is 1.
        // If dim > 1, multiple independent reactions. We just pick one valid set?
        // Let's try to ensure all positive.
        const hasNeg = vec.some(v => v < -1e-9);
        const hasPos = vec.some(v => v > 1e-9);
        if (hasNeg && !hasPos) {
            vec = vec.map(v => -v);
        } else if (hasNeg && hasPos) {
            // Mixed signs in a single basis vector implies mechanism issue or structure?
            // Or maybe A -> B + C where A is LHS (+), B,C are RHS (-).
            // My matrix construction: LHS (+), RHS (-).
            // Sol is [a, b, c].
            // a(LHS) - c(RHS) = 0.
            // If a, b, c are all positive, then signs are correct.
            // If kernel gives mixed signs, it means we might have a valid mathematical sol that isn't physical?
            // But kernel of [1, -2] is [2, 1]. All pos.
            // Kernel of [1, 1, -2] is [-1, 1, 0] and [2, 0, 1]?
            // Usually we want positive integers.
            // Let's rely on _kernel giving rational numbers.
        }

        // Convert to fractions to find LCM
        // Simple approach: find max denominator approx
        const denoms = [];
        const nums = [];
        for(let v of vec) {
            // approximate fraction
            const tolerance = 1.0E-6;
            let h1=1, h2=0, k1=0, k2=1;
            let b = v;
            do {
                let a = Math.floor(b);
                let aux = h1; h1 = a*h1+h2; h2 = aux;
                aux = k1; k1 = a*k1+k2; k2 = aux;
                b = 1/(b-a);
            } while (Math.abs(v - h1/k1) > v*tolerance);
            nums.push(h1);
            denoms.push(k1);
        }

        // LCM of denominators
        const gcd = (a, b) => !b ? a : gcd(b, a % b);
        const lcm = (a, b) => (a * b) / gcd(a, b);
        const globalLcm = denoms.reduce((acc, val) => lcm(acc, val), 1);

        // Scale
        const integers = vec.map((v, i) => Math.round(v * globalLcm));

        // Format result string
        let lhsStr = "";
        let rhsStr = "";

        let idx = 0;
        for(let i=0; i<lhsNodes.length; i++) {
            const coef = integers[idx++];
            if (coef !== 0) {
                if (lhsStr) lhsStr += " + ";
                lhsStr += (coef === 1 ? "" : coef + " ") + lhsNodes[i].toString();
            }
        }
        for(let i=0; i<rhsNodes.length; i++) {
            const coef = integers[idx++];
            if (coef !== 0) {
                if (rhsStr) rhsStr += " + ";
                rhsStr += (coef === 1 ? "" : coef + " ") + rhsNodes[i].toString();
            }
        }

        const resStr = `${lhsStr} -> ${rhsStr}`;
        return {
            type: 'info',
            text: resStr,
            toLatex: () => `\\text{${resStr.replace(/->/g, '\\rightarrow ')}}`
        };
    }

    _diagonalize(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("diagonalize requires a matrix");
        // P^-1 A P = D => A = P D P^-1
        // P contains eigenvectors
        // D contains eigenvalues

        const evals = this._eigenvals(matrix); // Vec of eigenvalues
        const evects = this._eigenvects(matrix); // Vec of eigenvectors

        // Check dimensions
        const n = matrix.elements.length;
        if (evects.elements.length < n) {
            throw new Error("Matrix is not diagonalizable (not enough eigenvectors)");
        }

        // Construct P from eigenvectors (columns)
        const P = this._trans(new Vec(evects.elements.slice(0, n))); // evects are rows in list? No, eigenvects returns list of vectors.
        // Vectors are usually column vectors, but represented as Vec (list).
        // If we want them as columns of P, we put them as rows then transpose.
        // Yes.

        // Construct D from eigenvalues
        // We need to match order. _eigenvects iterates _eigenvals.
        // But _eigenvects might return multiple vectors for one eigenvalue.
        // We need to reconstruct the correct D diagonal.
        // The simple implementation of _eigenvects iterates unique evals and pushes vectors.
        // So we should reconstruct the full list of eigenvalues corresponding to the vectors.

        const diagValues = [];
        // Re-run logic to match:
        const uniqueEvals = [];
        // Flatten evals just in case
        evals.elements.forEach(e => uniqueEvals.push(e)); // Assuming unique? No _eigenvals returns roots.
        // Actually _eigenvals returns distinct roots if 'set' used, or all roots?
        // My _eigenvals implementation uses 'solve' which might return a set of unique roots.
        // But _eigenvects iterates them.
        // If an eigenvalue has algebraic multiplicity > 1, solve returns it once in a set?
        // _eigenvects finds kernel. Kernel dimension (geometric multiplicity) gives number of vectors.
        // So for each vector in evects, we need to know which eigenvalue generated it.

        // Refined Logic:
        // Iterate unique eigenvalues, find kernel, add (val, vec) pairs.
        const pairs = [];
        // Get unique eigenvalues (solve returns set or single)
        // Note: _eigenvals calls _solve. If multiple roots, returns set.
        // Set contains unique values?
        // Let's assume unique.

        const distinctEvals = (evals instanceof Vec) ? evals.elements : [evals];
        // Dedupe check
        const seen = new Set();
        const distinct = [];
        distinctEvals.forEach(e => {
            const s = e.toString();
            if(!seen.has(s)) { seen.add(s); distinct.push(e); }
        });

        for(const lam of distinct) {
             // kernel(A - lam*I)
             const M_minus_lambdaI = [];
             for(let r=0; r<n; r++) {
                 const row = [];
                 for(let c=0; c<n; c++) {
                     let val = matrix.elements[r].elements[c];
                     if (r === c) val = new Sub(val, lam).simplify();
                     row.push(val);
                 }
                 M_minus_lambdaI.push(new Vec(row));
             }
             const kern = this._kernel(new Vec(M_minus_lambdaI));
             for(const v of kern.elements) {
                 pairs.push({ val: lam, vec: v });
             }
        }

        if (pairs.length < n) throw new Error("Matrix is not diagonalizable");

        const P_cols = [];
        const D_diag = [];

        for(let i=0; i<n; i++) {
            P_cols.push(pairs[i].vec);
            D_diag.push(pairs[i].val);
        }

        const P_mat = this._trans(new Vec(P_cols));
        const D_mat = this._diag(new Vec(D_diag));

        return new Vec([P_mat, D_mat]);
    }

    _tTest(data, mu0) {
        // One-sample t-test
        // t = (x_bar - mu0) / (s / sqrt(n))
        if (!(data instanceof Vec)) throw new Error("Data must be a list");
        const n = new Num(data.elements.length);
        const mean = this._mean(data);
        const std = this._std(data); // Sample std dev

        const num = new Sub(mean, mu0).simplify();
        const den = new Div(std, new Call('sqrt', [n])).simplify();
        const t = new Div(num, den).simplify();

        // Degrees of freedom
        const df = new Sub(n, new Num(1)).simplify();

        // P-value requires CDF of t-distribution (regularized incomplete beta)
        // This is complex to implement fully symbolic/numeric without a library.
        // We return { t, df } and maybe a note.
        // Or we can approximate for large n > 30 using Normal?
        // Let's return a result vector [t, df]
        // Users can look up table or we add tCDF later.

        return new Vec([t, df]);
    }

    _polyRegression(data, order) {
        if (!(data instanceof Vec)) throw new Error("polyRegression data must be a list");
        if (!(order instanceof Num)) order = new Num(2);
        const d = order.value;
        const n = data.elements.length;
        if (n <= d) throw new Error(`Not enough points for degree ${d} regression`);

        // Build X matrix and Y vector
        const X_rows = [];
        const Y_rows = [];

        for(const pt of data.elements) {
            // pt must be [x, y]
            let x, y;
            if (pt instanceof Vec && pt.elements.length >= 2) {
                x = pt.elements[0].evaluateNumeric();
                y = pt.elements[1].evaluateNumeric();
            }

            if (isNaN(x) || isNaN(y)) throw new Error("Regression data must be numeric");

            const row = [];
            for(let j=0; j<=d; j++) {
                row.push(new Num(Math.pow(x, j)));
            }
            X_rows.push(new Vec(row));
            Y_rows.push(new Num(y));
        }

        const X = new Vec(X_rows);

        // Y needs to be column vector for matrix mult
        const Y_col = new Vec(Y_rows.map(y => new Vec([y])));

        // Beta = (X^T * X)^-1 * X^T * Y
        const XT = this._trans(X);
        const XTX = new Mul(XT, X).simplify();
        const XTY = new Mul(XT, Y_col).simplify();

        const invXTX = this._inv(XTX);
        const beta = new Mul(invXTX, XTY).simplify();

        // beta is column vector [b0, b1, ... bd]^T
        // Poly: b0 + b1*x + ... + bd*x^d
        let poly = new Num(0);
        const xVar = new Sym('x');

        if (beta instanceof Vec) {
            for(let i=0; i<=d; i++) {
                // beta elements are rows (Vec), each row has 1 element
                const b = beta.elements[i].elements[0];
                poly = new Add(poly, new Mul(b, new Pow(xVar, new Num(i)))).simplify();
            }
        }
        return poly;
    }

    _expRegression(data) {
        if (!(data instanceof Vec)) throw new Error("Data must be list");
        // Transform y -> ln(y)
        const logData = [];
        for(const pt of data.elements) {
             let x, y;
             if (pt instanceof Vec && pt.elements.length >= 2) {
                 x = pt.elements[0].evaluateNumeric();
                 y = pt.elements[1].evaluateNumeric();
             }

             if (isNaN(x) || isNaN(y)) throw new Error("Regression data must be numeric");
             if (y <= 0) throw new Error("Exponential regression requires positive y values");

             const logYNum = new Num(Math.log(y));
             logData.push(new Vec([new Num(x), logYNum]));
        }

        const linData = new Vec(logData);
        // Reuse linearRegression manually to extract slope/intercept numeric values easily
        // linearRegression returns A + Bx.
        // We need A and B values.

        let sumX=0, sumY=0, sumXY=0, sumXX=0, n=0;
        for(const pt of logData) {
             const x = pt.elements[0].value;
             const y = pt.elements[1].value;
             sumX += x; sumY += y; sumXY += x*y; sumXX += x*x;
             n++;
        }

        const den = n*sumXX - sumX*sumX;
        if(den===0) throw new Error("Vertical line");

        const slope = (n*sumXY - sumX*sumY) / den; // b
        const intercept = (sumY - slope*sumX) / n; // ln(a)

        const a = Math.exp(intercept);
        const b = slope;

        // a * e^(b*x)
        return new Mul(new Num(a), new Call('exp', [new Mul(new Num(b), new Sym('x'))])).simplify();
    }

    _kron(A, B) {
        if (!(A instanceof Vec) || !(B instanceof Vec)) throw new Error("kron arguments must be matrices");
        const rowsA = A.elements.length;
        if (rowsA === 0) return new Vec([]);
        const colsA = A.elements[0].elements.length;

        const rowsB = B.elements.length;
        if (rowsB === 0) return new Vec([]);
        const colsB = B.elements[0].elements.length;

        const newRows = [];
        for(let r=0; r < rowsA * rowsB; r++) {
            const row = [];
            // which row in A?
            const rA = Math.floor(r / rowsB);
            const rB = r % rowsB;
            for(let c=0; c < colsA * colsB; c++) {
                const cA = Math.floor(c / colsB);
                const cB = c % colsB;
                const val = new Mul(A.elements[rA].elements[cA], B.elements[rB].elements[cB]).simplify();
                row.push(val);
            }
            newRows.push(new Vec(row));
        }
        return new Vec(newRows);
    }

    _geoMean(list) {
        if (list instanceof Vec) {
            const n = list.elements.length;
            if (n === 0) return new Num(0);
            let prod = new Num(1);
            for(const e of list.elements) prod = new Mul(prod, e);
            // nth root
            return new Pow(prod.simplify(), new Div(new Num(1), new Num(n))).simplify();
        }
        return new Call('geoMean', [list]);
    }

    _harmMean(list) {
        if (list instanceof Vec) {
            const n = list.elements.length;
            if (n === 0) return new Num(0);
            let sumInv = new Num(0);
            for(const e of list.elements) sumInv = new Add(sumInv, new Div(new Num(1), e));
            return new Div(new Num(n), sumInv.simplify()).simplify();
        }
        return new Call('harmMean', [list]);
    }

    _rms(list) {
        if (list instanceof Vec) {
            const n = list.elements.length;
            if (n === 0) return new Num(0);
            let sumSq = new Num(0);
            for(const e of list.elements) sumSq = new Add(sumSq, new Pow(e, new Num(2)));
            return new Call('sqrt', [new Div(sumSq, new Num(n))]).simplify();
        }
        return new Call('rms', [list]);
    }

    _mad(list) {
        // Mean Absolute Deviation
        if (list instanceof Vec) {
             const n = list.elements.length;
             if (n === 0) return new Num(0);
             const m = this._mean(list);
             let sum = new Num(0);
             for(const e of list.elements) {
                 const diff = new Sub(e, m).simplify();
                 sum = new Add(sum, new Call('abs', [diff]));
             }
             return new Div(sum, new Num(n)).simplify();
        }
        return new Call('mad', [list]);
    }

    _curvature(expr, varNode, point) {
        // kappa = |y''| / (1 + y'^2)^(3/2)
        const d1 = expr.diff(varNode).simplify();
        const d2 = d1.diff(varNode).simplify();

        const num = new Call('abs', [d2]);
        const den = new Pow(new Add(new Num(1), new Pow(d1, new Num(2))), new Num(1.5));

        const kappa = new Div(num, den).simplify();

        if (point) {
            return kappa.substitute(varNode, point).simplify();
        }
        return kappa;
    }

    _parallel(args) {
        if (args.length === 0) return new Num(0);
        let sumInv = new Num(0);
        for(const arg of args) {
            sumInv = new Add(sumInv, new Div(new Num(1), arg));
        }
        // Result = 1 / sum(1/Zn)
        return new Div(new Num(1), sumInv.simplify()).simplify();
    }

    _cis(angleDeg) {
        // e^(i * angle * pi / 180)
        // = cos(rad) + i*sin(rad)
        const rad = new Mul(angleDeg, new Div(new Sym('pi'), new Num(180))).simplify();
        return new Add(
            new Call('cos', [rad]),
            new Mul(new Sym('i'), new Call('sin', [rad]))
        ).simplify();
    }

    _phasor(mag, angleDeg) {
        return new Mul(mag, this._cis(angleDeg)).simplify();
    }

    _getComplexParts(expr) {
        expr = expr.simplify();
        if (expr instanceof Num) return { re: expr, im: new Num(0) };
        if (expr instanceof Sym) {
            if (expr.name === 'i' || expr.name === 'j') return { re: new Num(0), im: new Num(1) };
            return { re: expr, im: new Num(0) }; // Assume real parameter
        }
        if (expr instanceof Mul) {
            // Check for i
            if (expr.right instanceof Sym && (expr.right.name === 'i' || expr.right.name === 'j')) {
                return { re: new Num(0), im: expr.left };
            }
            if (expr.left instanceof Sym && (expr.left.name === 'i' || expr.left.name === 'j')) {
                return { re: new Num(0), im: expr.right };
            }
            // Recurse? If real * complex?
            const l = this._getComplexParts(expr.left);
            const r = this._getComplexParts(expr.right);
            // (a+bi)(c+di) = (ac-bd) + (ad+bc)i
            const ac = new Mul(l.re, r.re);
            const bd = new Mul(l.im, r.im);
            const ad = new Mul(l.re, r.im);
            const bc = new Mul(l.im, r.re);
            return {
                re: new Sub(ac, bd).simplify(),
                im: new Add(ad, bc).simplify()
            };
        }
        if (expr instanceof Div) {
            // (a+bi)/(c+di) = ((ac+bd) + (bc-ad)i) / (c^2+d^2)
            const num = this._getComplexParts(expr.left);
            const den = this._getComplexParts(expr.right);
            const denSq = new Add(new Pow(den.re, new Num(2)), new Pow(den.im, new Num(2))).simplify();

            const reNum = new Add(new Mul(num.re, den.re), new Mul(num.im, den.im));
            const imNum = new Sub(new Mul(num.im, den.re), new Mul(num.re, den.im));

            return {
                re: new Div(reNum, denSq).simplify(),
                im: new Div(imNum, denSq).simplify()
            };
        }
        if (expr instanceof Add) {
            const l = this._getComplexParts(expr.left);
            const r = this._getComplexParts(expr.right);
            return {
                re: new Add(l.re, r.re).simplify(),
                im: new Add(l.im, r.im).simplify()
            };
        }
        if (expr instanceof Sub) {
            const l = this._getComplexParts(expr.left);
            const r = this._getComplexParts(expr.right);
            return {
                re: new Sub(l.re, r.re).simplify(),
                im: new Sub(l.im, r.im).simplify()
            };
        }
        if (expr instanceof Pow) {
             // Handle simple cases
             if (expr.right instanceof Num && expr.right.value === 2) {
                 // (a+bi)^2 = a^2 - b^2 + 2abi
                 const base = this._getComplexParts(expr.left);
                 const re = new Sub(new Pow(base.re, new Num(2)), new Pow(base.im, new Num(2)));
                 const im = new Mul(new Num(2), new Mul(base.re, base.im));
                 return { re: re.simplify(), im: im.simplify() };
             }
        }

        // Default fallback (treat as real)
        return { re: expr, im: new Num(0) };
    }

    _toPolar(z) {
        z = z.simplify();
        const parts = this._getComplexParts(z);
        const re = parts.re;
        const im = parts.im;

        // Magnitude r = sqrt(re^2 + im^2)
        const rSq = new Add(new Pow(re, new Num(2)), new Pow(im, new Num(2))).simplify();
        const r = new Call('sqrt', [rSq]).simplify();

        // Angle (degrees)
        // Try numeric
        let deg;
        const reVal = re.evaluateNumeric();
        const imVal = im.evaluateNumeric();

        if (!isNaN(reVal) && !isNaN(imVal)) {
            const rad = Math.atan2(imVal, reVal);
            const d = rad * 180 / Math.PI;
            deg = new Num(d);
        } else {
            // Symbolic: 180/pi * arg(z)
            // Or atan(im/re) logic but keeping quadrant is hard symbolically without atan2
            deg = new Mul(new Div(new Num(180), new Sym('pi')), new Call('arg', [z])).simplify();
        }

        return new Vec([r, deg]);
    }

    _svd(matrix) {
        // U, S, V such that A = U * S * V^T
        // A^T A = V S^2 V^T
        // 1. Calculate B = A^T A
        const AT = this._trans(matrix);
        const B = new Mul(AT, matrix).simplify();

        // 2. Eigenvalues of B
        const evals = this._eigenvals(B);
        // Get numeric values if possible for sorting
        let evList = [];
        if (evals instanceof Vec) {
             for(const e of evals.elements) evList.push(e);
        } else {
             evList.push(evals);
        }

        // Filter and Sort eigenvalues (descending)
        // Values should be non-negative for A^T A
        evList.sort((a, b) => {
            const va = a.evaluateNumeric();
            const vb = b.evaluateNumeric();
            return vb - va;
        });

        // 3. Construct S (Singular Values)
        const sValues = [];
        for(const e of evList) {
             sValues.push(new Call('sqrt', [e]).simplify());
        }

        // 4. Construct V (Eigenvectors of A^T A)
        // We need vectors corresponding to sorted eigenvalues.
        // Group eigenvalues to handle multiplicity.
        const vCols = [];
        const processedIndices = new Set(); // Track used eigenvalues by index in evList

        // Map unique eigenvalues to their basis vectors
        // Optimization: Don't recompute kernel for same eigenvalue value
        // Note: evList is sorted. Repeated values are adjacent (if numeric) or potentially scattered if symbolic.
        // We loop through evList.

        // Cache kernels for values
        const kernelCache = []; // { valStr: string, vectors: Vec[] }

        for(let i=0; i<evList.length; i++) {
             const val = evList[i];
             const valStr = val.toString();

             let vectors = null;
             // Check cache
             const cached = kernelCache.find(c => c.valStr === valStr);
             if (cached) {
                 vectors = cached.vectors;
             } else {
                 // Compute kernel
                 const n = B.elements.length;
                 const M_minus_lambdaI = [];
                 for(let r=0; r<n; r++) {
                     const row = [];
                     for(let c=0; c<n; c++) {
                         let el = B.elements[r].elements[c];
                         if (r === c) el = new Sub(el, val).simplify();
                         row.push(el);
                     }
                     M_minus_lambdaI.push(new Vec(row));
                 }
                 const kern = this._kernel(new Vec(M_minus_lambdaI));
                 // kern returns a basis. We should Gram-Schmidt this basis to ensure orthogonality within the eigenspace
                 // if dimension > 1. _kernel usually returns orthogonal basis if RREF is clean?
                 // No, RREF basis is not necessarily orthogonal.
                 // Apply Gram-Schmidt to kernel vectors
                 const orthVectors = (kern.elements.length > 1)
                    ? this._gramschmidt(kern).elements
                    : kern.elements;

                 vectors = orthVectors;
                 kernelCache.push({ valStr: valStr, vectors: vectors });
             }

             // Consume one vector from the available basis
             // We need to keep track of how many we used for this value
             // If eigenvalues are repeated k times, we expect k vectors.
             // Find how many times we've seen this value so far
             let count = 0;
             for(let j=0; j<i; j++) {
                 if (evList[j].toString() === valStr) count++;
             }

             if (count < vectors.length) {
                 const vec = vectors[count];
                 // Normalize
                 const n = new Call('norm', [vec]).simplify();
                 vCols.push(new Div(vec, n).simplify());
             } else {
                 // Fallback: Not enough eigenvectors found? (Defective matrix or precision issue)
                 // Pad with zero?
                 vCols.push(this._zeros(new Num(B.elements.length), new Num(1)));
             }
        }
        const V = this._trans(new Vec(vCols)); // Columns are eigenvectors

        // 5. Construct U
        // u_i = (1/sigma_i) * A * v_i
        const uCols = [];
        const m = matrix.elements.length;
        for(let i=0; i<sValues.length; i++) {
             const sigma = sValues[i];
             const v = vCols[i]; // v_i
             const val = sigma.evaluateNumeric();
             if (Math.abs(val) > 1e-9) {
                 const Av = new Mul(matrix, v).simplify();
                 uCols.push(new Div(Av, sigma).simplify());
             } else {
                 // For zero singular values, we need to complete basis for U
                 // kernel(A^T) gives remaining directions?
                 // For now, push zero vector or handle properly?
                 // SVD usually requires orthogonal U.
                 // This part is complex without full Gram-Schmidt completion on kernel(A^T).
                 // Placeholder: 0 vector
                 const z = []; for(let k=0; k<m; k++) z.push(new Num(0));
                 uCols.push(new Vec(z));
             }
        }
        const U = this._trans(new Vec(uCols));
        const S = this._diag(new Vec(sValues));

        // Result: [U, S, V] (V, not V^T, following standard svd returns)
        return new Vec([U, S, V]);
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
