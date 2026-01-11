
class CAS {
    constructor() {
        this.variables = {
            'pi': new Num(Math.PI),
            'true': new Num(1),
            'false': new Num(0),
            'j': new Sym('i') // Electrical Engineering imaginary unit
        };

        this.functions = {};

        // Unit Conversion Rates (Base: SI units)
        this.conversionRates = {
            // Length (Base: m)
            'm': 1, 'km': 1000, 'cm': 0.01, 'mm': 0.001, 'um': 1e-6, 'nm': 1e-9,
            'ft': 0.3048, 'in': 0.0254, 'yd': 0.9144, 'mi': 1609.34,
            // Mass (Base: kg)
            'kg': 1, 'g': 0.001, 'mg': 1e-6, 'lb': 0.453592, 'oz': 0.0283495, 't': 1000,
            // Time (Base: s)
            's': 1, 'min': 60, 'h': 3600, 'd': 86400, 'wk': 604800, 'y': 31536000, 'ms': 0.001,
            // Area (Base: m^2)
            'm^2': 1, 'km^2': 1e6, 'cm^2': 1e-4, 'mm^2': 1e-6, 'ft^2': 0.092903, 'in^2': 0.00064516, 'ac': 4046.86, 'ha': 10000,
            // Volume (Base: m^3)
            'm^3': 1, 'km^3': 1e9, 'cm^3': 1e-6, 'mm^3': 1e-9, 'L': 0.001, 'mL': 1e-6, 'gal': 0.00378541, 'ft^3': 0.0283168,
            // Speed (Base: m/s)
            'm/s': 1, 'km/h': 0.277778, 'mph': 0.44704, 'kn': 0.514444, 'ft/s': 0.3048,
            // Pressure (Base: Pa)
            'Pa': 1, 'kPa': 1000, 'bar': 100000, 'atm': 101325, 'psi': 6894.76, 'Torr': 133.322,
            // Energy (Base: J)
            'J': 1, 'kJ': 1000, 'cal': 4.184, 'kcal': 4184, 'eV': 1.60218e-19, 'kWh': 3.6e6,
            // Power (Base: W)
            'W': 1, 'kW': 1000, 'hp': 745.7
        };
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

            if (node.funcName === 'roots') {
                if (node.args.length < 2) throw new Error("roots requires at least 2 arguments: polynomial, variable");
                // Evaluate arguments
                const poly = this._recursiveEval(node.args[0]);
                const varNode = this._recursiveEval(node.args[1]);
                return this._roots(poly, varNode);
            }

            if (node.funcName === 'subs' || node.funcName === 'substitute') {
                 if (node.args.length !== 2) throw new Error("subs requires 2 arguments: expr, substitution (var=val)");
                 const expr = this._recursiveEval(node.args[0]);
                 const sub = this._recursiveEval(node.args[1]);

                 if (sub instanceof Eq) {
                     return expr.substitute(sub.left, sub.right).simplify();
                 }
                 if (sub instanceof Vec) {
                     let res = expr;
                     for(const el of sub.elements) {
                         if (el instanceof Eq) {
                             res = res.substitute(el.left, el.right);
                         }
                     }
                     return res.simplify();
                 }
                 throw new Error("subs second argument must be an equation (var=val) or list of equations");
            }

            if (node.funcName === 'help') {
                const aliasMap = {
                    'int': 'integrate',
                    'derivative': 'diff',
                    'differentiation': 'diff',
                    'antiderivative': 'integrate',
                    'comb': 'nCr',
                    'perm': 'nPr',
                    'evalf': 'approx',
                    'approx': 'approx',
                    'tran': 'trans',
                    'transpose': 'trans',
                    'eye': 'identity',
                    'idn': 'identity',
                    'ker': 'kernel',
                    'nullspace': 'kernel',
                    'rem': 'mod',
                    'binomial': 'nCr',
                    'stddev': 'std',
                    'var': 'variance',
                    'molar': 'molarMass',
                    'subs': 'substitute'
                };

                const getHelpList = () => {
                    if (typeof globalThis.HELP_DATA !== 'undefined') {
                        const keys = Object.keys(globalThis.HELP_DATA).sort();
                        return "Available commands:\n" + keys.join(', ');
                    }
                    return "Help data not loaded.";
                };

                // Check for argument help("cmd") or help(cmd)
                if (node.args.length === 1) {
                    // Extract command name
                    let cmd = "";
                    if (node.args[0] instanceof Sym) cmd = node.args[0].name;

                    if (cmd) {
                        // Resolve alias
                        if (aliasMap[cmd]) cmd = aliasMap[cmd];

                        const data = (typeof globalThis.HELP_DATA !== 'undefined') ? globalThis.HELP_DATA[cmd] : null;
                        if (data) {
                            return {
                                type: 'help',
                                command: cmd,
                                data: data,
                                toString: () => `${cmd}: ${data.description}\nUsage: ${data.syntax}`,
                                toLatex: () => `\\text{${cmd}: ${data.description}}`
                            };
                        } else {
                            // Try to find partial matches
                            let suggestions = "";
                            if (typeof globalThis.HELP_DATA !== 'undefined') {
                                const matches = Object.keys(globalThis.HELP_DATA).filter(k => k.includes(cmd));
                                if (matches.length > 0) {
                                    suggestions = "\nDid you mean: " + matches.join(', ');
                                }
                            }

                            return {
                                type: 'info',
                                text: `No specific help found for '${cmd}'.${suggestions}`,
                                toLatex: () => `\\text{No help found for ${cmd}}`
                            };
                        }
                    }
                }

                // Default Help
                const helpText = getHelpList();
                const latexHelp = `\\text{Available commands: see list}`;
                return { type: 'help', text: helpText, toString: () => helpText, toLatex: () => latexHelp };
            }

            const args = node.args.map(arg => this._recursiveEval(arg));

            // Check for user-defined function
            if (this.functions.hasOwnProperty(node.funcName)) {
                const funcDef = this.functions[node.funcName];
                if (funcDef.params.length !== node.args.length) {
                    throw new Error(`Function ${node.funcName} expects ${funcDef.params.length} arguments, got ${node.args.length}`);
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
                if (node.args.length !== 1) throw new Error("N requires 1 argument");
                const val = args[0].evaluateNumeric();
                return new Num(val);
            }

            if (node.funcName === 'diff') {
                if (node.args.length < 2) throw new Error("diff requires at least 2 arguments");
                const func = args[0];
                const varNode = args[1];
                if (!(varNode instanceof Sym)) throw new Error("Second argument to diff must be a variable");

                let order = 1;
                if (node.args.length >= 3) {
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

            if (node.funcName === 'implicitDiff') {
                if (node.args.length !== 3) throw new Error("implicitDiff requires 3 arguments: equation, depVar(y), indepVar(x)");
                return this._implicitDiff(args[0], args[1], args[2]);
            }

            if (node.funcName === 'integrate' || node.funcName === 'int') {
                if (node.args.length < 2) throw new Error("integrate requires at least 2 arguments");
                const func = args[0];
                const varNode = args[1];
                if (!(varNode instanceof Sym)) throw new Error("Second argument to integrate must be a variable");

                if (node.args.length === 4) {
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
                        // Check for Abs in Definite Integral
                        // integrate(abs(expr), x, a, b)
                        if (func instanceof Call && func.funcName === 'abs') {
                            const inner = func.args[0];
                            // Find roots of inner expression
                            const roots = this._roots(inner, varNode);
                            // Filter roots within [lower, upper]
                            const validRoots = [];
                            const a = lower.evaluateNumeric();
                            const b = upper.evaluateNumeric();

                            if (!isNaN(a) && !isNaN(b)) {
                                if (roots instanceof Vec) {
                                    for(const r of roots.elements) {
                                        const rVal = r.evaluateNumeric();
                                        if (!isNaN(rVal) && rVal > Math.min(a,b) && rVal < Math.max(a,b)) {
                                            validRoots.push(rVal);
                                        }
                                    }
                                }

                                if (validRoots.length > 0) {
                                    validRoots.sort((u,v) => u-v);
                                    let sum = new Num(0);
                                    let prev = lower;

                                    for(const r of validRoots) {
                                        const rNum = new Num(r);
                                        // Integrate from prev to r
                                        // Check sign in midpoint
                                        const mid = (prev.evaluateNumeric() + r) / 2;
                                        const valAtMid = inner.substitute(varNode, new Num(mid)).evaluateNumeric();
                                        const integrand = valAtMid >= 0 ? inner : new Mul(new Num(-1), inner);

                                        const part = this.evaluate(new Call('integrate', [integrand, varNode, prev, rNum]));
                                        sum = new Add(sum, part);
                                        prev = rNum;
                                    }

                                    // Last segment
                                    const mid = (prev.evaluateNumeric() + b) / 2;
                                    const valAtMid = inner.substitute(varNode, new Num(mid)).evaluateNumeric();
                                    const integrand = valAtMid >= 0 ? inner : new Mul(new Num(-1), inner);
                                    const part = this.evaluate(new Call('integrate', [integrand, varNode, prev, upper]));
                                    sum = new Add(sum, part);

                                    return sum.simplify();
                                }
                            }
                        }

                        return new Call('integrate', args);
                    }

                    // Use limit for bounds to handle singularities (e.g. 1/sqrt(x) at 0)
                    // limit(F(x), x, upper, left) - limit(F(x), x, lower, right)
                    let valUpper = indefinite.substitute(varNode, upper);
                    // Check if substitution failed (NaN or Infinity symbolic)
                    // Also check for evaluateNumeric issues if needed, but symbolic check covers basic cases
                    const isBad = (n) => (n instanceof Sym && (n.name === 'NaN' || n.name === 'Infinity' || n.name === 'infinity')) ||
                                         (n instanceof Mul && n.left instanceof Num && n.left.value === -1 && (n.right.name === 'Infinity' || n.right.name === 'infinity'));

                    if (isBad(valUpper)) {
                         valUpper = this._limit(indefinite, varNode, upper, 0, -1);
                    }

                    let valLower = indefinite.substitute(varNode, lower);
                    if (isBad(valLower)) {
                         valLower = this._limit(indefinite, varNode, lower, 0, 1);
                    }

                    return new Sub(valUpper, valLower).simplify();
                }

                // Check for standard integrals that are not handled by Expr.integrate
                // 1. ArcTan form: 1/(u^2 + a^2) * u' -> (1/a)*atan(u/a)
                // Also handles 1/(Ax^2 + Bx + C) by completing the square
                const atanRes = this._integrateAtan(func, varNode);
                if (atanRes) return atanRes;

                // 2. Inverse Hyperbolic / Trig with Sqrt: 1/sqrt(Q(x))
                const invHypRes = this._integrateInverseHyperbolic(func, varNode);
                if (invHypRes) return invHypRes;

                // 3. Standard Trig Integrals: tan(x), sec(x), etc.
                const trigRes = this._integrateTrig(func, varNode);
                if (trigRes) return trigRes;

                // 4. Special Functions (Error Function, etc.)
                const specialRes = this._integrateSpecial(func, varNode);
                if (specialRes) return specialRes;

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

                     // 4. Try Trigonometric Reduction (Power Reduction)
                     // e.g. sin(x)^2 -> (1-cos(2x))/2
                     const reduced = this._linearizeTrig(func).simplify();
                     if (reduced.toString() !== func.toString()) {
                         const redInt = reduced.integrate(varNode).simplify();
                         if (!(redInt instanceof Call && redInt.funcName === 'integrate')) return redInt;
                     }

                     // Try Integration by Parts
                     const parts = this._integrateByParts(func, varNode);
                     if (parts) return parts;
                }
                return res;
            }

            if (node.funcName === 'sum') {
                if (node.args.length === 1) return this._sumList(args[0]);
                // sum(expr, var, start, end)
                if (node.args.length !== 4) throw new Error("sum requires 1 argument (list) or 4 arguments: expression, variable, start, end");
                return this._sum(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'product') {
                if (node.args.length === 1) return this._productList(args[0]);
                // product(expr, var, start, end)
                if (node.args.length !== 4) throw new Error("product requires 1 argument (list) or 4 arguments: expression, variable, start, end");
                return this._product(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'cumsum') {
                if (node.args.length !== 1) throw new Error("cumsum requires 1 argument (list)");
                return this._cumsum(args[0]);
            }

            if (node.funcName === 'flatten') {
                if (node.args.length !== 1) throw new Error("flatten requires 1 argument (list)");
                return this._flatten(args[0]);
            }

            if (node.funcName === 'expand') {
                if (node.args.length !== 1) throw new Error("expand takes exactly 1 argument");
                return args[0].expand();
            }

            if (node.funcName === 'simplify') {
                if (node.args.length !== 1) throw new Error("simplify takes exactly 1 argument");
                let simplified = args[0].simplify();
                // If result is a complex Add/Sub expression, try expand().simplify() to see if it reduces
                if (simplified instanceof Add || simplified instanceof Sub) {
                    const expanded = simplified.expand().simplify();
                    // Heuristic: prefer shorter string length or if it becomes 0
                    if (expanded instanceof Num && expanded.value === 0) return expanded;
                    if (expanded.toString().length < simplified.toString().length) return expanded;
                }
                return simplified;
            }

            if (node.funcName === 'solve' || node.funcName === 'fsolve') {
                 if (node.args.length < 2) throw new Error("solve requires at least 2 arguments: equation and variable");
                 const eq = args[0];
                 const varNode = args[1];

                 // Check for guess (Newton-Raphson)
                 if (node.args.length === 3) {
                     const guess = args[2];
                     if (!(varNode instanceof Sym) && !(varNode instanceof Vec)) throw new Error("Second argument to fsolve must be a variable or list of variables");
                     return this._fsolve(eq, varNode, guess);
                 }

                 if (!(varNode instanceof Sym) && !(varNode instanceof Vec)) throw new Error("Second argument to solve must be a variable or list of variables");
                 return this._solve(eq, varNode);
            }

            if (node.funcName === 'nIntegrate' || node.funcName === 'numeric_integrate') {
                if (node.args.length !== 4) throw new Error("nIntegrate requires 4 arguments: expr, var, start, end");
                return this._nIntegrate(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'minimize' || node.funcName === 'fmin') {
                 if (node.args.length !== 2) throw new Error("minimize requires 2 arguments: expr, variable");
                 if (!(args[1] instanceof Sym)) throw new Error("Second argument to minimize must be a variable");
                 return this._minimize(args[0], args[1]);
            }

            if (node.funcName === 'maximize' || node.funcName === 'fmax') {
                 if (node.args.length !== 2) throw new Error("maximize requires 2 arguments: expr, variable");
                 if (!(args[1] instanceof Sym)) throw new Error("Second argument to maximize must be a variable");
                 return this._maximize(args[0], args[1]);
            }

            if (node.funcName === 'analyze') {
                 if (node.args.length !== 2) throw new Error("analyze requires 2 arguments: expr, variable");
                 if (!(args[1] instanceof Sym)) throw new Error("Second argument to analyze must be a variable");
                 return this._analyze(args[0], args[1]);
            }

            if (node.funcName === 'extrema') {
                 if (node.args.length !== 2) throw new Error("extrema requires 2 arguments: expr, variable");
                 if (!(args[1] instanceof Sym)) throw new Error("Second argument to extrema must be a variable");
                 return this._extrema(args[0], args[1]);
            }

            if (node.funcName === 'stationary_points') {
                 if (node.args.length !== 2) throw new Error("stationary_points requires 2 arguments: expr, variable");
                 if (!(args[1] instanceof Sym)) throw new Error("Second argument to stationary_points must be a variable");
                 return this._stationary_points(args[0], args[1]);
            }

            if (node.funcName === 'asymptotes') {
                 if (node.args.length !== 2) throw new Error("asymptotes requires 2 arguments: expr, variable");
                 if (!(args[1] instanceof Sym)) throw new Error("Second argument to asymptotes must be a variable");
                 return this._asymptotes(args[0], args[1]);
            }

            if (node.funcName === 'completeSquare') {
                 if (node.args.length !== 2) throw new Error("completeSquare requires 2 arguments: expr, variable");
                 if (!(args[1] instanceof Sym)) throw new Error("Second argument to completeSquare must be a variable");
                 return this._completeSquare(args[0], args[1]);
            }

            if (node.funcName === 'resultant') {
                 if (node.args.length !== 3) throw new Error("resultant requires 3 arguments: poly1, poly2, variable");
                 return this._resultant(args[0], args[1], args[2]);
            }

            if (node.funcName === 'discriminant') {
                 if (node.args.length !== 2) throw new Error("discriminant requires 2 arguments: poly, variable");
                 return this._discriminant(args[0], args[1]);
            }

            if (node.funcName === 'plot') {
                 if (node.args.length < 2) throw new Error("plot requires at least 2 arguments: expression and variable");
                 const expr = args[0];
                 const varNode = args[1];
                 let min = -10, max = 10;
                 let nstep = null;

                 // Handle args
                 // args[2] could be min, args[3] max
                 // Check for optional nstep=...
                 let posArgCount = 0;
                 for(let i=2; i<node.args.length; i++) {
                     const arg = args[i];
                     if (arg instanceof Eq && arg.left instanceof Sym && arg.left.name === 'nstep') {
                         const val = arg.right.evaluateNumeric();
                         if (!isNaN(val) && val > 0) nstep = val;
                     } else {
                         // Positional logic: 0->min, 1->max
                         const val = arg.evaluateNumeric();
                         if (!isNaN(val)) {
                             if (posArgCount === 0) min = val;
                             if (posArgCount === 1) max = val;
                             posArgCount++;
                         }
                     }
                 }

                 // Return a special object that the frontend can recognize
                 return {
                     type: 'plot',
                     subtype: 'function',
                     expr: expr,
                     var: varNode,
                     min: min,
                     max: max,
                     nstep: nstep,
                     toString: () => `Plotting ${expr} from ${min} to ${max}`,
                     toLatex: () => `\\text{Plotting } ${expr.toLatex()}`
                 };
            }

            if (node.funcName === 'plot3d') {
                // plot3d(expr, x, y, [x_min, x_max, y_min, y_max])
                if (node.args.length < 3) throw new Error("plot3d requires at least 3 arguments: expression, var1, var2");
                const expr = args[0];
                const varX = args[1];
                const varY = args[2];
                let xMin = -10, xMax = 10, yMin = -10, yMax = 10;

                if (node.args.length >= 7) {
                    xMin = args[3].evaluateNumeric();
                    xMax = args[4].evaluateNumeric();
                    yMin = args[5].evaluateNumeric();
                    yMax = args[6].evaluateNumeric();
                } else if (node.args.length >= 5) {
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

                let xExpr, yExpr, tVar, min = -10, max = 10, nstep = null;
                let argStartIdx = 0;

                if (node.args.length >= 2 && args[0] instanceof Vec) {
                    // plotparam([x, y], t, [min], [max])
                    if (args[0].elements.length !== 2) throw new Error("plotparam list must contain 2 expressions");
                    xExpr = args[0].elements[0];
                    yExpr = args[0].elements[1];
                    tVar = args[1];
                    argStartIdx = 2;
                } else if (node.args.length >= 3) {
                    // plotparam(x, y, t, [min], [max])
                    xExpr = args[0];
                    yExpr = args[1];
                    tVar = args[2];
                    argStartIdx = 3;
                } else {
                    throw new Error("plotparam requires arguments: x, y, t, [min, max]");
                }

                // Parse remaining args for min, max, nstep
                let posArgCount = 0;
                for(let i=argStartIdx; i<node.args.length; i++) {
                    const arg = args[i];
                    if (arg instanceof Eq && arg.left instanceof Sym && arg.left.name === 'nstep') {
                        const val = arg.right.evaluateNumeric();
                        if (!isNaN(val) && val > 0) nstep = val;
                    } else {
                        const val = arg.evaluateNumeric();
                        if (!isNaN(val)) {
                            if (posArgCount === 0) min = val;
                            if (posArgCount === 1) max = val;
                            posArgCount++;
                        }
                    }
                }

                if (!(tVar instanceof Sym)) throw new Error("plotparam variable must be a symbol");

                return {
                    type: 'plot',
                    subtype: 'parametric',
                    xExpr: xExpr,
                    yExpr: yExpr,
                    var: tVar,
                    min: min,
                    max: max,
                    nstep: nstep,
                    toString: () => `Parametric Plot (${xExpr}, ${yExpr}) t=${min}..${max}`,
                    toLatex: () => `\\text{Parametric Plot } (${xExpr.toLatex()}, ${yExpr.toLatex()})`
                };
            }

            if (node.funcName === 'plotpolar') {
                // plotpolar(r, theta, min, max)
                if (node.args.length < 3) throw new Error("plotpolar requires at least 3 arguments: r, theta, min, max");
                const rExpr = args[0];
                const thetaVar = args[1];
                let min = 0, max = 2 * Math.PI, nstep = null;

                // Parse remaining args
                let posArgCount = 0;
                for(let i=2; i<node.args.length; i++) {
                    const arg = args[i];
                    if (arg instanceof Eq && arg.left instanceof Sym && arg.left.name === 'nstep') {
                        const val = arg.right.evaluateNumeric();
                        if (!isNaN(val) && val > 0) nstep = val;
                    } else {
                        const val = arg.evaluateNumeric();
                        if (!isNaN(val)) {
                            if (posArgCount === 0) min = val;
                            if (posArgCount === 1) max = val;
                            posArgCount++;
                        }
                    }
                }

                if (!(thetaVar instanceof Sym)) throw new Error("plotpolar variable must be a symbol");

                return {
                    type: 'plot',
                    subtype: 'polar',
                    rExpr: rExpr,
                    var: thetaVar,
                    min: min,
                    max: max,
                    nstep: nstep,
                    toString: () => `Polar Plot r=${rExpr} theta=${min}..${max}`,
                    toLatex: () => `\\text{Polar Plot } r=${rExpr.toLatex()}`
                };
            }

            if (node.funcName === 'plotlist') {
                if (node.args.length !== 1) throw new Error("plotlist requires 1 argument (list)");
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
                if (node.args.length < 4) throw new Error("taylor requires 4 arguments: expression, variable, point, order");
                const expr = args[0];
                const varNode = args[1];
                const point = args[2];
                const order = args[3];

                if (!(varNode instanceof Sym)) throw new Error("Second argument to taylor must be a variable");
                if (!(order instanceof Num)) throw new Error("Order must be a number");

                return this._taylor(expr, varNode, point, order.value);
            }

            if (node.funcName === 'limit' || node.funcName === 'lim') {
                // limit(expr, var, point, [dir])
                if (node.args.length < 3) throw new Error("limit requires 3 arguments: expression, variable, point");
                const expr = args[0];
                const varNode = args[1];
                const point = args[2];
                let dir = 0;
                if (node.args.length > 3) {
                    const dArg = args[3];
                    const d = dArg.evaluateNumeric();
                    if (!isNaN(d)) {
                        dir = d;
                    } else if (dArg instanceof Sym) {
                        const name = dArg.name.toLowerCase();
                        if (name === 'right' || name === 'plus') dir = 1;
                        if (name === 'left' || name === 'minus') dir = -1;
                    } else if (dArg.toString() === '+') dir = 1;
                    else if (dArg.toString() === '-') dir = -1;
                }

                if (!(varNode instanceof Sym)) throw new Error("Second argument to limit must be a variable");

                return this._limit(expr, varNode, point, 0, dir);
            }

            if (node.funcName === 'potential') {
                 if (node.args.length !== 2) throw new Error("potential requires 2 arguments: vector field, vars");
                 return this._potential(args[0], args[1]);
            }

            if (node.funcName === 'conservative') {
                 if (node.args.length !== 2) throw new Error("conservative requires 2 arguments: vector field, vars");
                 return this._conservative(args[0], args[1]);
            }

            if (node.funcName === 'tangent') {
                // tangent(expr, var, point)
                if (node.args.length < 3) throw new Error("tangent requires 3 arguments: expression, variable, point");
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
                if (node.args.length !== 1) throw new Error("trigReduce requires 1 argument");
                return this._linearizeTrig(args[0]);
            }

            if (node.funcName === 'trigExpand') {
                if (node.args.length !== 1) throw new Error("trigExpand requires 1 argument");
                return args[0].expand().simplify();
            }

            if (node.funcName === 'det') {
                if (node.args.length !== 1) throw new Error("det requires 1 argument");
                return this._det(args[0]);
            }

            if (node.funcName === 'inv') {
                if (node.args.length !== 1) throw new Error("inv requires 1 argument");
                return this._inv(args[0]);
            }

            if (node.funcName === 'rref') {
                if (node.args.length !== 1) throw new Error("rref requires 1 argument");
                return this._rref(args[0]);
            }

            if (node.funcName === 'rank') {
                if (node.args.length !== 1) throw new Error("rank requires 1 argument");
                return this._rank(args[0]);
            }

            if (node.funcName === 'cross') {
                if (node.args.length !== 2) throw new Error("cross requires 2 arguments");
                return this._cross(args[0], args[1]);
            }

            if (node.funcName === 'primeFactors') {
                if (node.args.length !== 1) throw new Error("primeFactors requires 1 argument");
                return this._primeFactors(args[0]);
            }

            if (node.funcName === 'angle') {
                if (node.args.length !== 2) throw new Error("angle requires 2 arguments: u, v");
                return this._angle(args[0], args[1]);
            }

            if (node.funcName === 'projection') {
                if (node.args.length !== 2) throw new Error("projection requires 2 arguments: u, v");
                return this._projection(args[0], args[1]);
            }

            if (node.funcName === 'toSpherical') {
                if (node.args.length !== 1) throw new Error("toSpherical requires 1 argument: vector [x, y, z]");
                return this._toSpherical(args[0]);
            }

            if (node.funcName === 'toCylindrical') {
                if (node.args.length !== 1) throw new Error("toCylindrical requires 1 argument: vector [x, y, z]");
                return this._toCylindrical(args[0]);
            }

            if (node.funcName === 'cbrt') {
                if (node.args.length !== 1) throw new Error("cbrt requires 1 argument");
                return new Pow(args[0], new Div(new Num(1), new Num(3))).simplify();
            }

            if (node.funcName === 'trans' || node.funcName === 'transpose' || node.funcName === 'tran') {
                if (node.args.length !== 1) throw new Error("trans requires 1 argument");
                return this._trans(args[0]);
            }

            if (node.funcName === 'eye' || node.funcName === 'idn') {
                if (node.args.length !== 1) throw new Error("eye requires 1 argument");
                return this._identity(args[0]);
            }

            if (node.funcName === 'zeros') {
                if (node.args.length !== 2) throw new Error("zeros requires 2 arguments: rows, cols");
                return this._zeros(args[0], args[1]);
            }

            if (node.funcName === 'ones') {
                if (node.args.length !== 2) throw new Error("ones requires 2 arguments: rows, cols");
                return this._ones(args[0], args[1]);
            }

            if (node.funcName === 'binomial' || node.funcName === 'comb') {
                if (node.args.length !== 2) throw new Error("binomial requires 2 arguments");
                return this._nCr(args[0], args[1]);
            }

            if (node.funcName === 'divisors') {
                if (node.args.length !== 1) throw new Error("divisors requires 1 argument");
                return this._divisors(args[0]);
            }

            if (node.funcName === 'union') {
                if (node.args.length !== 2) throw new Error("union requires 2 arguments (lists)");
                return this._union(args[0], args[1]);
            }

            if (node.funcName === 'intersect') {
                if (node.args.length !== 2) throw new Error("intersect requires 2 arguments (lists)");
                return this._intersect(args[0], args[1]);
            }

            if (node.funcName === 'setdiff') {
                if (node.args.length !== 2) throw new Error("setdiff requires 2 arguments (lists)");
                return this._setdiff(args[0], args[1]);
            }

            if (node.funcName === 'clear') {
                return { type: 'action', name: 'clear', toString: () => "Cleared", toLatex: () => "\\text{Cleared}" };
            }

            if (node.funcName === 'gcd') {
                if (node.args.length === 1 && args[0] instanceof Vec) {
                    const list = args[0].elements;
                    if (list.length === 0) return new Num(0);
                    let res = list[0];
                    for(let i=1; i<list.length; i++) res = this._gcd(res, list[i]);
                    return res;
                }
                if (node.args.length < 2) throw new Error("gcd requires at least 2 arguments");
                let res = args[0];
                for(let i=1; i<node.args.length; i++) res = this._gcd(res, args[i]);
                return res;
            }

            if (node.funcName === 'lcm') {
                if (node.args.length === 1 && args[0] instanceof Vec) {
                    const list = args[0].elements;
                    if (list.length === 0) return new Num(1);
                    let res = list[0];
                    for(let i=1; i<list.length; i++) res = this._lcm(res, list[i]);
                    return res;
                }
                if (node.args.length < 2) throw new Error("lcm requires at least 2 arguments");
                let res = args[0];
                for(let i=1; i<node.args.length; i++) res = this._lcm(res, args[i]);
                return res;
            }

            if (node.funcName === 'factor' || node.funcName === 'ifactor') {
                if (node.args.length !== 1) throw new Error("factor requires 1 argument");
                return this._factor(args[0]);
            }

            if (node.funcName === 'factorial') {
                 if (node.args.length !== 1) throw new Error("factorial requires 1 argument");
                 const n = args[0];
                 if (n instanceof Num && Number.isInteger(n.value) && n.value >= 0) {
                     return new Num(this._factorial(n.value));
                 }
                 return new Call('factorial', args);
            }

            if (node.funcName === 'fibonacci') {
                 if (node.args.length !== 1) throw new Error("fibonacci requires 1 argument");
                 const n = args[0];
                 if (n instanceof Num && Number.isInteger(n.value) && n.value >= 0) {
                     return new Num(this._fibonacci(n.value));
                 }
                 return new Call('fibonacci', args);
            }

            if (node.funcName === 'gamma') {
                 if (node.args.length !== 1) throw new Error("gamma requires 1 argument");
                 const z = args[0];
                 if (z instanceof Num) {
                     return new Num(this._gamma(z.value));
                 }
                 return new Call('gamma', args);
            }

            if (node.funcName === 'beta') {
                 if (node.args.length !== 2) throw new Error("beta requires 2 arguments");
                 return this._beta(args[0], args[1]);
            }

            if (node.funcName === 'nCr') {
                if (node.args.length !== 2) throw new Error("nCr requires 2 arguments");
                return this._nCr(args[0], args[1]);
            }

            if (node.funcName === 'nPr' || node.funcName === 'perm') {
                if (node.args.length !== 2) throw new Error("nPr requires 2 arguments");
                return this._nPr(args[0], args[1]);
            }

            if (node.funcName === 'isPrime') {
                if (node.args.length !== 1) throw new Error("isPrime requires 1 argument");
                return this._isPrime(args[0]);
            }

            if (node.funcName === 'nextprime') {
                if (node.args.length !== 1) throw new Error("nextprime requires 1 argument");
                return this._nextprime(args[0]);
            }

            if (node.funcName === 'prevprime') {
                if (node.args.length !== 1) throw new Error("prevprime requires 1 argument");
                return this._prevprime(args[0]);
            }

            if (node.funcName === 'trace') {
                if (node.args.length !== 1) throw new Error("trace requires 1 argument");
                return this._trace(args[0]);
            }

            if (node.funcName === 'mean') {
                if (node.args.length !== 1) throw new Error("mean requires 1 argument (list)");
                return this._mean(args[0]);
            }

            if (node.funcName === 'variance' || node.funcName === 'var') {
                if (node.args.length !== 1) throw new Error("variance requires 1 argument (list)");
                return this._variance(args[0]);
            }

            if (node.funcName === 'std' || node.funcName === 'stddev') {
                if (node.args.length !== 1) throw new Error("std requires 1 argument (list)");
                return this._std(args[0]);
            }

            if (node.funcName === 'cov') {
                if (node.args.length !== 2) throw new Error("cov requires 2 arguments (list1, list2)");
                return this._cov(args[0], args[1]);
            }

            if (node.funcName === 'corr') {
                if (node.args.length !== 2) throw new Error("corr requires 2 arguments (list1, list2)");
                return this._corr(args[0], args[1]);
            }

            if (node.funcName === 'median') {
                if (node.args.length !== 1) throw new Error("median requires 1 argument (list)");
                return this._median(args[0]);
            }

            if (node.funcName === 'charpoly') {
                if (node.args.length !== 2) throw new Error("charpoly requires 2 arguments: matrix, variable");
                return this._charpoly(args[0], args[1]);
            }

            if (node.funcName === 'linearRegression') {
                if (node.args.length !== 1) throw new Error("linearRegression requires 1 argument (list of points)");
                return this._linearRegression(args[0]);
            }

            if (node.funcName === 'polyRegression') {
                if (node.args.length < 2) throw new Error("polyRegression requires 2 arguments: data, degree");
                return this._polyRegression(args[0], args[1]);
            }

            if (node.funcName === 'expRegression') {
                if (node.args.length !== 1) throw new Error("expRegression requires 1 argument: data");
                return this._expRegression(args[0]);
            }

            if (node.funcName === 'powerRegression') {
                if (node.args.length !== 1) throw new Error("powerRegression requires 1 argument: data");
                return this._powerRegression(args[0]);
            }

            if (node.funcName === 'logRegression') {
                if (node.args.length !== 1) throw new Error("logRegression requires 1 argument: data");
                return this._logRegression(args[0]);
            }

            if (node.funcName === 'normalPDF') {
                if (node.args.length !== 3) throw new Error("normalPDF requires 3 arguments: x, mu, sigma");
                return this._normalPDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'binomialPDF') {
                if (node.args.length !== 3) throw new Error("binomialPDF requires 3 arguments: k, n, p");
                return this._binomialPDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'normalCDF') {
                if (node.args.length !== 3) throw new Error("normalCDF requires 3 arguments: x, mu, sigma");
                return this._normalCDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'invNorm') {
                if (node.args.length !== 3) throw new Error("invNorm requires 3 arguments: area, mu, sigma");
                return this._invNorm(args[0], args[1], args[2]);
            }

            if (node.funcName === 'binomialCDF') {
                if (node.args.length !== 3) throw new Error("binomialCDF requires 3 arguments: k, n, p");
                return this._binomialCDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'poissonPDF') {
                if (node.args.length !== 2) throw new Error("poissonPDF requires 2 arguments: k, lambda");
                return this._poissonPDF(args[0], args[1]);
            }

            if (node.funcName === 'poissonCDF') {
                if (node.args.length !== 2) throw new Error("poissonCDF requires 2 arguments: k, lambda");
                return this._poissonCDF(args[0], args[1]);
            }

            if (node.funcName === 'exponentialPDF') {
                if (node.args.length !== 2) throw new Error("exponentialPDF requires 2 arguments: x, lambda");
                return this._exponentialPDF(args[0], args[1]);
            }

            if (node.funcName === 'exponentialCDF') {
                if (node.args.length !== 2) throw new Error("exponentialCDF requires 2 arguments: x, lambda");
                return this._exponentialCDF(args[0], args[1]);
            }

            if (node.funcName === 'geometricPDF') {
                if (node.args.length !== 2) throw new Error("geometricPDF requires 2 arguments: k, p");
                return this._geometricPDF(args[0], args[1]);
            }

            if (node.funcName === 'geometricCDF') {
                if (node.args.length !== 2) throw new Error("geometricCDF requires 2 arguments: k, p");
                return this._geometricCDF(args[0], args[1]);
            }

            if (node.funcName === 'chisquarePDF') {
                if (node.args.length !== 2) throw new Error("chisquarePDF requires 2 arguments: x, k");
                return this._chisquarePDF(args[0], args[1]);
            }

            if (node.funcName === 'chisquareCDF') {
                if (node.args.length !== 2) throw new Error("chisquareCDF requires 2 arguments: x, k");
                return this._chisquareCDF(args[0], args[1]);
            }

            if (node.funcName === 'invChiSquare') {
                if (node.args.length !== 2) throw new Error("invChiSquare requires 2 arguments: area, k");
                return this._invChiSquare(args[0], args[1]);
            }

            if (node.funcName === 'studentTPDF') {
                if (node.args.length !== 2) throw new Error("studentTPDF requires 2 arguments: x, df");
                return this._studentTPDF(args[0], args[1]);
            }

            if (node.funcName === 'studentTCDF') {
                if (node.args.length !== 2) throw new Error("studentTCDF requires 2 arguments: x, df");
                return this._studentTCDF(args[0], args[1]);
            }

            if (node.funcName === 'fPDF') {
                if (node.args.length !== 3) throw new Error("fPDF requires 3 arguments: x, d1, d2");
                return this._fPDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'fCDF') {
                if (node.args.length !== 3) throw new Error("fCDF requires 3 arguments: x, d1, d2");
                return this._fCDF(args[0], args[1], args[2]);
            }

            if (node.funcName === 'invT') {
                if (node.args.length !== 2) throw new Error("invT requires 2 arguments: area, df");
                return this._invT(args[0], args[1]);
            }

            if (node.funcName === 'compound') {
                if (node.args.length !== 4) throw new Error("compound requires 4 arguments: P, r, n, t");
                return this._compound(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'loan') {
                if (node.args.length !== 3) throw new Error("loan requires 3 arguments: P, r, n");
                return this._loan(args[0], args[1], args[2]);
            }

            if (node.funcName === 'npv') {
                if (node.args.length !== 2) throw new Error("npv requires 2 arguments: rate, cash_flows");
                return this._npv(args[0], args[1]);
            }

            if (node.funcName === 'irr') {
                if (node.args.length !== 1) throw new Error("irr requires 1 argument: cash_flows");
                return this._irr(args[0]);
            }

            if (node.funcName === 'dot') {
                if (node.args.length !== 2) throw new Error("dot requires 2 arguments");
                // dot(u, v) is u * v (Mul handles dot product for vectors)
                return new Mul(args[0], args[1]).simplify();
            }

            if (node.funcName === 'norm') {
                if (node.args.length !== 1) throw new Error("norm requires 1 argument");
                // L2 norm: sqrt(v . v)
                const v = args[0];
                const dot = new Mul(v, v).simplify();
                return new Call('sqrt', [dot]).simplify();
            }

            if (node.funcName === 'grad') {
                if (node.args.length !== 2) throw new Error("grad requires 2 arguments: expression and list of variables");
                return this._grad(args[0], args[1]);
            }

            if (node.funcName === 'curl') {
                if (node.args.length !== 2) throw new Error("curl requires 2 arguments: vector field and list of variables");
                return this._curl(args[0], args[1]);
            }

            if (node.funcName === 'divergence' || node.funcName === 'div') { // 'div' might conflict with division if not careful, but funcName is safe
                 if (node.args.length !== 2) throw new Error("divergence requires 2 arguments: vector field and list of variables");
                 return this._divergence(args[0], args[1]);
            }

            if (node.funcName === 'jacobian') {
                 if (node.args.length !== 2) throw new Error("jacobian requires 2 arguments: vector, vars");
                 return this._jacobian(args[0], args[1]);
            }

            if (node.funcName === 'hessian') {
                 if (node.args.length !== 2) throw new Error("hessian requires 2 arguments: expr, vars");
                 return this._hessian(args[0], args[1]);
            }

            if (node.funcName === 'laplacian') {
                 // laplacian(expr, vars) = div(grad(expr))
                 if (node.args.length !== 2) throw new Error("laplacian requires 2 arguments: expression and list of variables");
                 const grad = this._grad(args[0], args[1]);
                 return this._divergence(grad, args[1]);
            }

            if (node.funcName === 'rem' || node.funcName === 'irem') {
                if (node.args.length !== 2) throw new Error("rem requires 2 arguments");
                return this._rem(args[0], args[1]);
            }

            if (node.funcName === 'quo') {
                if (node.args.length !== 2) throw new Error("quo requires 2 arguments");
                return this._quo(args[0], args[1]);
            }

            if (node.funcName === 'mod') {
                if (node.args.length !== 2) throw new Error("mod requires 2 arguments");
                return this._mod(args[0], args[1]);
            }

            if (node.funcName === 'modInverse') {
                if (node.args.length !== 2) throw new Error("modInverse requires 2 arguments: a, m");
                return this._modInverse(args[0], args[1]);
            }

            if (node.funcName === 'modPow') {
                if (node.args.length !== 3) throw new Error("modPow requires 3 arguments: base, exp, mod");
                return this._modPow(args[0], args[1], args[2]);
            }

            if (node.funcName === 'partfrac') {
                if (node.args.length !== 2) throw new Error("partfrac requires 2 arguments: expr, var");
                return this._partfrac(args[0], args[1]);
            }

            if (node.funcName === 'size' || node.funcName === 'dim' || node.funcName === 'length') {
                if (node.args.length !== 1) throw new Error("size/dim requires 1 argument");
                if (args[0] instanceof Vec) {
                    return new Num(args[0].elements.length);
                }
                return new Call(node.funcName, args);
            }

            if (node.funcName === 'concat') {
                if (node.args.length < 2) throw new Error("concat requires at least 2 arguments");
                return this._concat(args);
            }

            if (node.funcName === 'append') {
                if (node.args.length !== 2) throw new Error("append requires 2 arguments: list, element");
                return this._append(args[0], args[1]);
            }

            if (node.funcName === 'prepend') {
                if (node.args.length !== 2) throw new Error("prepend requires 2 arguments: list, element");
                return this._prepend(args[0], args[1]);
            }

            if (node.funcName === 'approx' || node.funcName === 'evalf') {
                if (node.args.length !== 1) throw new Error("approx requires 1 argument");
                return new Num(args[0].evaluateNumeric());
            }

            if (node.funcName === 'distance') {
                if (node.args.length !== 2) throw new Error("distance requires 2 arguments");
                return this._distance(args[0], args[1]);
            }

            if (node.funcName === 'midpoint') {
                if (node.args.length !== 2) throw new Error("midpoint requires 2 arguments");
                return this._midpoint(args[0], args[1]);
            }

            if (node.funcName === 'zTest') {
                if (node.args.length !== 3) throw new Error("zTest requires 3 arguments: data, mu0, sigma");
                return this._zTest(args[0], args[1], args[2]);
            }

            if (node.funcName === 'euler' || node.funcName === 'phi') {
                if (node.args.length !== 1) throw new Error("euler requires 1 argument");
                return this._euler(args[0]);
            }

            if (node.funcName === 'mode') {
                if (node.args.length !== 1) throw new Error("mode requires 1 argument (list)");
                return this._mode(args[0]);
            }

            if (node.funcName === 'arcLen' || node.funcName === 'arcLength') {
                if (node.args.length !== 4) throw new Error("arcLen requires 4 arguments: expr, var, start, end");
                return this._arcLen(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'volume_solid' || node.funcName === 'volumeSolid') {
                // volume_solid(expr, var, a, b, [axis])
                // axis: 'x' (default, disc) or 'y' (shell, assumes expr is y=f(x))
                if (node.args.length < 4) throw new Error("volume_solid requires at least 4 arguments: expr, var, start, end");
                const axis = node.args.length > 4 ? args[4] : new Sym('x');
                return this._volumeSolid(args[0], args[1], args[2], args[3], axis);
            }

            if (node.funcName === 'surfaceArea' || node.funcName === 'area_surface') {
                 if (node.args.length < 4) throw new Error("surfaceArea requires at least 4 arguments: expr, var, start, end");
                 const axis = node.args.length > 4 ? args[4] : new Sym('x');
                 return this._surfaceArea(args[0], args[1], args[2], args[3], axis);
            }

            if (node.funcName === 'arg') {
                if (node.args.length !== 1) throw new Error("arg requires 1 argument");
                return this._arg(args[0]);
            }

            if (node.funcName === 'degree') {
                if (node.args.length !== 2) throw new Error("degree requires 2 arguments: expr, var");
                return this._degree(args[0], args[1]);
            }
            if (node.funcName === 'coeff') {
                if (node.args.length !== 3) throw new Error("coeff requires 3 arguments: expr, var, degree");
                return this._coeff(args[0], args[1], args[2]);
            }
            if (node.funcName === 'symb2poly') {
                if (node.args.length !== 2) throw new Error("symb2poly requires 2 arguments: expr, var");
                return this._symb2poly(args[0], args[1]);
            }
            if (node.funcName === 'poly2symb') {
                if (node.args.length !== 2) throw new Error("poly2symb requires 2 arguments: list, var");
                return this._poly2symb(args[0], args[1]);
            }
            if (node.funcName === 'seq') {
                if (node.args.length !== 5) throw new Error("seq requires 5 arguments: expr, var, start, end, step");
                return this._seq(args[0], args[1], args[2], args[3], args[4]);
            }
            if (node.funcName === 'range') {
                if (node.args.length !== 3) throw new Error("range requires 3 arguments: start, end, step");
                return this._range(args[0], args[1], args[2]);
            }
            if (node.funcName === 'sort') {
                if (node.args.length !== 1) throw new Error("sort requires 1 argument: list");
                return this._sort(args[0]);
            }
            if (node.funcName === 'reverse') {
                if (node.args.length !== 1) throw new Error("reverse requires 1 argument: list");
                return this._reverse(args[0]);
            }
            if (node.funcName === 'diag') {
                if (node.args.length !== 1) throw new Error("diag requires 1 argument: list");
                return this._diag(args[0]);
            }
            if (node.funcName === 'identity') {
                if (node.args.length !== 1) throw new Error("identity requires 1 argument: n");
                return this._identity(args[0]);
            }
            if (node.funcName === 'laplace') {
                if (node.args.length !== 3) throw new Error("laplace requires 3 arguments: expr, t, s");
                return this._laplace(args[0], args[1], args[2]);
            }
            if (node.funcName === 'ilaplace') {
                if (node.args.length !== 3) throw new Error("ilaplace requires 3 arguments: expr, s, t");
                return this._ilaplace(args[0], args[1], args[2]);
            }
            if (node.funcName === 'ztrans') {
                if (node.args.length !== 3) throw new Error("ztrans requires 3 arguments: expr, n, z");
                return this._ztrans(args[0], args[1], args[2]);
            }
            if (node.funcName === 'iztrans') {
                if (node.args.length !== 3) throw new Error("iztrans requires 3 arguments: expr, z, n");
                return this._iztrans(args[0], args[1], args[2]);
            }

            if (node.funcName === 'kernel' || node.funcName === 'nullspace' || node.funcName === 'ker') {
                if (node.args.length !== 1) throw new Error("kernel requires 1 argument");
                return this._kernel(args[0].simplify());
            }
            if (node.funcName === 'basis') {
                if (node.args.length !== 1) throw new Error("basis requires 1 argument");
                return this._basis(args[0].simplify());
            }
            if (node.funcName === 'eigenvals') {
                if (node.args.length !== 1) throw new Error("eigenvals requires 1 argument");
                return this._eigenvals(args[0].simplify());
            }
            if (node.funcName === 'eigenvects') {
                if (node.args.length !== 1) throw new Error("eigenvects requires 1 argument");
                return this._eigenvects(args[0].simplify());
            }
            if (node.funcName === 'gramschmidt') {
                if (node.args.length !== 1) throw new Error("gramschmidt requires 1 argument (list of vectors)");
                return this._gramschmidt(args[0].simplify());
            }
            if (node.funcName === 'lu') {
                if (node.args.length !== 1) throw new Error("lu requires 1 argument");
                return this._lu(args[0].simplify());
            }
            if (node.funcName === 'qr') {
                if (node.args.length !== 1) throw new Error("qr requires 1 argument");
                return this._qr(args[0].simplify());
            }
            if (node.funcName === 'cholesky') {
                if (node.args.length !== 1) throw new Error("cholesky requires 1 argument");
                return this._cholesky(args[0].simplify());
            }
            if (node.funcName === 'desolve') {
                if (node.args.length < 2) throw new Error("desolve requires at least 2 arguments");
                return this._desolve(args[0], args[1]);
            }

            if (node.funcName === 'rsolve') {
                // rsolve(eq, a(n), [conds])
                if (node.args.length < 2) throw new Error("rsolve requires at least 2 arguments: equation, recurrence_term");
                const conds = (node.args.length > 2) ? args[2] : null;
                return this._rsolve(args[0], args[1], conds);
            }

            if (node.funcName === 'fourier') {
                if (node.args.length < 3) throw new Error("fourier requires 3 arguments: expr, var, n, [L]");
                const expr = args[0];
                const varNode = args[1];
                const n = args[2];
                const L = node.args.length > 3 ? args[3] : new Sym('pi');
                if (!(varNode instanceof Sym)) throw new Error("Second argument to fourier must be a variable");
                return this._fourier(expr, varNode, n, L);
            }

            if (node.funcName === 'fourier_transform' || node.funcName === 'ft') {
                if (node.args.length !== 3) throw new Error("fourier_transform requires 3 arguments: expr, t, w");
                return this._fourierTransform(args[0], args[1], args[2]);
            }

            if (node.funcName === 'inverse_fourier_transform' || node.funcName === 'ift') {
                if (node.args.length !== 3) throw new Error("inverse_fourier_transform requires 3 arguments: expr, w, t");
                return this._inverseFourierTransform(args[0], args[1], args[2]);
            }

            if (node.funcName === 'latex' || node.funcName === 'tex') {
                if (node.args.length !== 1) throw new Error("latex requires 1 argument");
                const latexStr = args[0].toLatex();
                return {
                    type: 'info',
                    text: latexStr,
                    toLatex: () => `\\text{${latexStr.replace(/\\/g, '\\\\')}}` // Escape backslashes for display
                };
            }

            if (node.funcName === 'slopefield') {
                 // slopefield(diffEq, x, y, [minX, maxX, minY, maxY])
                 if (node.args.length < 3) throw new Error("slopefield requires at least 3 arguments: equation, x, y");
                 const eq = args[0];
                 const xVar = args[1];
                 const yVar = args[2];

                 let expr = eq;
                 // If eq is Eq (y' = ...), extract RHS
                 if (eq instanceof Eq) {
                     expr = eq.right;
                 }

                 let min = -10, max = 10, yMin = -10, yMax = 10;
                 if (node.args.length >= 7) {
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
                 if (node.args.length < 3) throw new Error("vectorfield requires at least 3 arguments: vector, x, y");
                 const vec = args[0];
                 const xVar = args[1];
                 const yVar = args[2];

                 if (!(vec instanceof Vec) || vec.elements.length !== 2) throw new Error("Vector field must be 2D vector [u, v]");

                 let min = -10, max = 10, yMin = -10, yMax = 10;
                 if (node.args.length >= 7) {
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
                 if (node.args.length < 3) throw new Error("plotimplicit requires at least 3 arguments: equation, x, y");
                 const eq = args[0];
                 const xVar = args[1];
                 const yVar = args[2];

                 let xMin = -10, xMax = 10, yMin = -10, yMax = 10;
                 if (node.args.length >= 7) {
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


            if (node.funcName === 'molarMass') {
                if (node.args.length !== 1) throw new Error("molarMass requires 1 argument (string formula)");
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
                if (node.args.length !== 1) throw new Error("atomicWeight requires 1 argument (symbol)");
                let sym = "";
                if (args[0] instanceof Sym) sym = args[0].name;
                if (!sym) throw new Error("Invalid element symbol");
                return this._atomicWeight(sym);
            }

            if (node.funcName === 'balance') {
                if (node.args.length !== 1) throw new Error("balance requires 1 argument (equation string)");
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
                if (node.args.length !== 1) throw new Error("diagonalize requires 1 argument (matrix)");
                return this._diagonalize(args[0]);
            }

            if (node.funcName === 'tTest') {
                if (node.args.length !== 2) throw new Error("tTest requires 2 arguments: data, mu0");
                return this._tTest(args[0], args[1]);
            }

            if (node.funcName === 'zInterval') {
                 // zInterval(data, sigma, level)
                 if (node.args.length !== 3) throw new Error("zInterval requires 3 arguments: data, sigma, level");
                 return this._zInterval(args[0], args[1], args[2]);
            }

            if (node.funcName === 'tInterval') {
                 // tInterval(data, level)
                 if (node.args.length !== 2) throw new Error("tInterval requires 2 arguments: data, level");
                 return this._tInterval(args[0], args[1]);
            }

            if (node.funcName === 'propTest') {
                 // propTest(successes, n, p0)
                 if (node.args.length !== 3) throw new Error("propTest requires 3 arguments: x, n, p0");
                 return this._propTest(args[0], args[1], args[2]);
            }

            if (node.funcName === 'propTest2') {
                 // propTest2(x1, n1, x2, n2)
                 if (node.args.length !== 4) throw new Error("propTest2 requires 4 arguments: x1, n1, x2, n2");
                 return this._propTest2(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'tTest2') {
                if (node.args.length !== 2) throw new Error("tTest2 requires 2 arguments: data1, data2");
                return this._tTest2(args[0], args[1]);
            }

            if (node.funcName === 'chiSquareTest') {
                if (node.args.length !== 2) throw new Error("chiSquareTest requires 2 arguments: observed, expected");
                return this._chiSquareTest(args[0], args[1]);
            }

            if (node.funcName === 'kron') {
                if (node.args.length !== 2) throw new Error("kron requires 2 arguments");
                return this._kron(args[0], args[1]);
            }

            if (node.funcName === 'svd') {
                if (node.args.length !== 1) throw new Error("svd requires 1 argument");
                return this._svd(args[0]);
            }

            if (node.funcName === 'geoMean') {
                if (node.args.length !== 1) throw new Error("geoMean requires 1 argument");
                return this._geoMean(args[0]);
            }

            if (node.funcName === 'harmMean') {
                if (node.args.length !== 1) throw new Error("harmMean requires 1 argument");
                return this._harmMean(args[0]);
            }

            if (node.funcName === 'rms') {
                if (node.args.length !== 1) throw new Error("rms requires 1 argument");
                return this._rms(args[0]);
            }

            if (node.funcName === 'mad') {
                if (node.args.length !== 1) throw new Error("mad requires 1 argument");
                return this._mad(args[0]);
            }

            if (node.funcName === 'moment') {
                 if (node.args.length !== 2) throw new Error("moment requires 2 arguments: list, k");
                 return this._moment(args[0], args[1]);
            }

            if (node.funcName === 'isSubset') {
                if (node.args.length !== 2) throw new Error("isSubset requires 2 arguments (lists)");
                return this._isSubset(args[0], args[1]);
            }

            if (node.funcName === 'cartesianProduct') {
                if (node.args.length !== 2) throw new Error("cartesianProduct requires 2 arguments (lists)");
                return this._cartesianProduct(args[0], args[1]);
            }

            if (node.funcName === 'skewness') {
                 if (node.args.length !== 1) throw new Error("skewness requires 1 argument");
                 return this._skewness(args[0]);
            }

            if (node.funcName === 'kurtosis') {
                 if (node.args.length !== 1) throw new Error("kurtosis requires 1 argument");
                 return this._kurtosis(args[0]);
            }

            if (node.funcName === 'curvature') {
                if (node.args.length < 2) throw new Error("curvature requires at least 2 arguments: expr, var, [point]");
                const point = node.args.length > 2 ? args[2] : null;
                return this._curvature(args[0], args[1], point);
            }

            if (node.funcName === 'par' || node.funcName === 'parallel') {
                return this._parallel(args);
            }

            if (node.funcName === 'cis') {
                if (node.args.length !== 1) throw new Error("cis requires 1 argument (angle in degrees)");
                return this._cis(args[0]);
            }

            if (node.funcName === 'phasor') {
                if (node.args.length !== 2) throw new Error("phasor requires 2 arguments: magnitude, angle(deg)");
                return this._phasor(args[0], args[1]);
            }

            if (node.funcName === 'toPolar') {
                if (node.args.length !== 1) throw new Error("toPolar requires 1 argument");
                return this._toPolar(args[0]);
            }

            if (node.funcName === 'bernoulli') {
                if (node.args.length !== 1) throw new Error("bernoulli requires 1 argument (n)");
                return this._bernoulli(args[0]);
            }

            if (node.funcName === 'harmonic') {
                if (node.args.length !== 1) throw new Error("harmonic requires 1 argument (n)");
                return this._harmonic(args[0]);
            }

            if (node.funcName === 'lineEquation') {
                if (node.args.length !== 2) throw new Error("lineEquation requires 2 arguments: point1, point2");
                return this._lineEquation(args[0], args[1]);
            }

            if (node.funcName === 'circleEquation') {
                if (node.args.length !== 2) throw new Error("circleEquation requires 2 arguments: center, radius");
                return this._circleEquation(args[0], args[1]);
            }

            if (node.funcName === 'planeEquation') {
                if (node.args.length !== 3) throw new Error("planeEquation requires 3 arguments: p1, p2, p3");
                return this._planeEquation(args[0], args[1], args[2]);
            }

            if (node.funcName === 'root') {
                 if (node.args.length !== 2) throw new Error("root requires 2 arguments: x, n");
                 // root(x, n) -> x^(1/n)
                 return new Pow(args[0], new Div(new Num(1), args[1])).simplify();
            }

            if (node.funcName === 'cfrac' || node.funcName === 'cf' || node.funcName === 'continued_fraction') {
                 // cfrac(val, [depth])
                 const depth = node.args.length > 1 ? args[1].evaluateNumeric() : 15;
                 return this._cfrac(args[0], depth);
            }

            if (node.funcName === 'sturm') {
                 if (node.args.length !== 2) throw new Error("sturm requires 2 arguments: poly, var");
                 return this._sturm(args[0], args[1]);
            }

            if (node.funcName === 'num_real_roots' || node.funcName === 'numRealRoots') {
                 // num_real_roots(poly, var, [a, b])
                 if (node.args.length < 2) throw new Error("num_real_roots requires at least 2 arguments");
                 const a = node.args.length > 2 ? args[2] : new Mul(new Num(-1), new Sym('Infinity'));
                 const b = node.args.length > 3 ? args[3] : new Sym('Infinity');
                 return this._numRealRoots(args[0], args[1], a, b);
            }

            if (node.funcName === 'isSquare') {
                 if (node.args.length !== 1) throw new Error("isSquare requires 1 argument");
                 return this._isSquare(args[0]);
            }

            if (node.funcName === 'truthTable') {
                 if (node.args.length !== 2) throw new Error("truthTable requires 2 arguments: expr, vars");
                 return this._truthTable(args[0], args[1]);
            }

            if (node.funcName === 'legendre') {
                 if (node.args.length !== 2) throw new Error("legendre requires 2 arguments: n, x");
                 return this._legendre(args[0], args[1]);
            }
            if (node.funcName === 'hermite') {
                 if (node.args.length !== 2) throw new Error("hermite requires 2 arguments: n, x");
                 return this._hermite(args[0], args[1]);
            }
            if (node.funcName === 'chebyshev') {
                 if (node.args.length !== 2) throw new Error("chebyshev requires 2 arguments: n, x");
                 return this._chebyshev(args[0], args[1]);
            }
            if (node.funcName === 'laguerre') {
                 if (node.args.length !== 2) throw new Error("laguerre requires 2 arguments: n, x");
                 return this._laguerre(args[0], args[1]);
            }

            if (node.funcName === 'pinv') {
                if (node.args.length !== 1) throw new Error("pinv requires 1 argument (matrix)");
                return this._pinv(args[0]);
            }

            if (node.funcName === 'cond') {
                if (node.args.length !== 1) throw new Error("cond requires 1 argument (matrix)");
                return this._cond(args[0]);
            }

            if (node.funcName === 'isDiagonal') {
                 if (node.args.length !== 1) throw new Error("isDiagonal requires 1 argument");
                 return this._isDiagonal(args[0]);
            }

            if (node.funcName === 'isSymmetric') {
                 if (node.args.length !== 1) throw new Error("isSymmetric requires 1 argument");
                 return this._isSymmetric(args[0]);
            }

            if (node.funcName === 'isOrthogonal') {
                 if (node.args.length !== 1) throw new Error("isOrthogonal requires 1 argument");
                 return this._isOrthogonal(args[0]);
            }

            if (node.funcName === 'collect') {
                if (node.args.length !== 2) throw new Error("collect requires 2 arguments: expr, var");
                return this._collect(args[0], args[1]);
            }

            if (node.funcName === 'matrix') {
                // matrix(rows, cols, expr, var1, var2)
                // or matrix(rows, cols, func) if we supported lambda
                // or matrix(rows, cols, element) -> fill
                if (node.args.length === 5) {
                    return this._matrixGen(args[0], args[1], args[2], args[3], args[4]);
                }
                // Fallback to zeros/ones-like behavior or just generic constructor if 3 args?
                // matrix(r, c, val) -> fill
                if (node.args.length === 3) {
                    const r = args[0];
                    const c = args[1];
                    const val = args[2];
                    if (r instanceof Num && c instanceof Num) {
                        // Create matrix filled with val
                        const rows = r.value;
                        const cols = c.value;
                        const m = [];
                        for(let i=0; i<rows; i++) {
                            const row = [];
                            for(let j=0; j<cols; j++) row.push(val);
                            m.push(new Vec(row));
                        }
                        return new Vec(m);
                    }
                }
                // If just matrix([[1,2]]) -> return Vec
                if (node.args.length === 1 && args[0] instanceof Vec) return args[0];

                throw new Error("matrix command usage: matrix(rows, cols, expr, i, j) or matrix(rows, cols, fillValue)");
            }

            if (node.funcName === 'min' || node.funcName === 'max') {
                // Handle list argument
                if (node.args.length === 1 && args[0] instanceof Vec) {
                    return this._minMaxList(args[0], node.funcName);
                }
                // Check if any arg is a Vec, if so, map?
                // min([1,2], [3,4]) -> [1, 2] elementwise?
                // Standard CAS usually treats min(list) as reduction.
                // Call handles varargs min(a,b,c).
                return new Call(node.funcName, args).simplify();
            }

            if (node.funcName === 'ceil' || node.funcName === 'floor') {
                if (node.args.length !== 1) throw new Error(`${node.funcName} requires 1 argument`);
                return new Call(node.funcName, args).simplify();
            }

            if (node.funcName === 'round') {
                // round(x, [n])
                return new Call(node.funcName, args).simplify();
            }

            if (node.funcName === 'moebius') {
                if (node.args.length !== 1) throw new Error("moebius requires 1 argument (integer)");
                return this._moebius(args[0]);
            }

            if (node.funcName === 'sigma' || node.funcName === 'divisorSum') {
                if (node.args.length < 1) throw new Error("sigma requires at least 1 argument");
                const k = node.args.length > 1 ? args[1] : new Num(1);
                return this._sigma(args[0], k);
            }

            if (node.funcName === 'legendreSymbol') {
                if (node.args.length !== 2) throw new Error("legendreSymbol requires 2 arguments: a, p");
                return this._legendreSymbol(args[0], args[1]);
            }

            if (node.funcName === 'stirling1') {
                if (node.args.length !== 2) throw new Error("stirling1 requires 2 arguments: n, k");
                return this._stirling1(args[0], args[1]);
            }

            if (node.funcName === 'stirling2') {
                if (node.args.length !== 2) throw new Error("stirling2 requires 2 arguments: n, k");
                return this._stirling2(args[0], args[1]);
            }

            if (node.funcName === 'bell') {
                if (node.args.length !== 1) throw new Error("bell requires 1 argument: n");
                return this._bell(args[0]);
            }

            if (node.funcName === 'isPerfect') {
                if (node.args.length !== 1) throw new Error("isPerfect requires 1 argument");
                return this._isPerfect(args[0]);
            }

            if (node.funcName === 'matrixPow' || node.funcName === 'matrix_pow') {
                if (node.args.length !== 2) throw new Error("matrixPow requires 2 arguments: matrix, n");
                return this._matrixPow(args[0], args[1]);
            }

            if (node.funcName === 'kroneckerDelta') {
                if (node.args.length !== 2) throw new Error("kroneckerDelta requires 2 arguments");
                return this._kroneckerDelta(args[0], args[1]);
            }

            if (node.funcName === 'fft') {
                if (node.args.length !== 1) throw new Error("fft requires 1 argument (list)");
                return this._fft(args[0]);
            }

            if (node.funcName === 'ifft') {
                if (node.args.length !== 1) throw new Error("ifft requires 1 argument (list)");
                return this._ifft(args[0]);
            }

            if (node.funcName === 'cnf') {
                if (node.args.length !== 1) throw new Error("cnf requires 1 argument");
                return this._cnf(args[0]);
            }

            if (node.funcName === 'dnf') {
                if (node.args.length !== 1) throw new Error("dnf requires 1 argument");
                return this._dnf(args[0]);
            }

            if (node.funcName === 'annuity') {
                // annuity(rate, n, payment) -> FV
                if (node.args.length !== 3) throw new Error("annuity requires 3 arguments: rate, n, payment");
                return this._annuity(args[0], args[1], args[2]);
            }

            if (node.funcName === 'amortization') {
                // amortization(rate, n, principal) -> payment
                if (node.args.length !== 3) throw new Error("amortization requires 3 arguments: rate, n, principal");
                return this._amortization(args[0], args[1], args[2]);
            }

            if (node.funcName === 'codegen' || node.funcName === 'toCode') {
                if (node.args.length !== 2) throw new Error("codegen requires 2 arguments: expr, language");
                // Expect language as string identifier
                let lang = "";
                if (args[1] instanceof Sym) lang = args[1].name;
                else if (args[1] instanceof Num) lang = args[1].toString(); // unlikely
                else throw new Error("Language must be an identifier (e.g. python, js, c)");
                return this._codegen(args[0], lang);
            }

            if (node.funcName === 'convert') {
                if (node.args.length !== 3) throw new Error("convert requires 3 arguments: value, fromUnit, toUnit");
                // fromUnit, toUnit as symbols
                let fromUnit = "", toUnit = "";
                if (args[1] instanceof Sym) fromUnit = args[1].name;
                if (args[2] instanceof Sym) toUnit = args[2].name;
                // handle power units like m^2 represented as Pow in AST if parser allows,
                // but parser likely sees m^2 as Pow(Sym(m), Num(2)).
                // We need to convert AST to unit string if complex.
                // Simple helper to get unit string
                const getUnitStr = (n) => {
                    if (n instanceof Sym) return n.name;
                    if (n instanceof Pow) return getUnitStr(n.left) + "^" + getUnitStr(n.right);
                    if (n instanceof Num) return n.value.toString();
                    if (n instanceof Div) return getUnitStr(n.left) + "/" + getUnitStr(n.right);
                    if (n instanceof Mul) return getUnitStr(n.left) + "*" + getUnitStr(n.right); // rarely used in simple converters
                    return "";
                };

                if (!fromUnit) fromUnit = getUnitStr(args[1]);
                if (!toUnit) toUnit = getUnitStr(args[2]);

                return this._convert(args[0], fromUnit, toUnit);
            }

            if (node.funcName === 'mse') {
                if (node.args.length !== 1) throw new Error("mse requires 1 argument (list)");
                return this._mse(args[0]);
            }

            if (node.funcName === 'mae') {
                if (node.args.length !== 1) throw new Error("mae requires 1 argument (list)");
                return this._mae(args[0]);
            }

            if (node.funcName === 'convolution') {
                if (node.args.length !== 3) throw new Error("convolution requires 3 arguments: f, g, t");
                return this._convolution(args[0], args[1], args[2]);
            }

            if (node.funcName === 'conv') {
                if (node.args.length !== 2) throw new Error("conv requires 2 arguments (lists)");
                return this._conv(args[0], args[1]);
            }

            if (node.funcName === 'xcorr') {
                if (node.args.length !== 2) throw new Error("xcorr requires 2 arguments (lists)");
                return this._xcorr(args[0], args[1]);
            }

            if (node.funcName === 'pade') {
                if (node.args.length < 4) throw new Error("pade requires 4 arguments: expr, var, n, m");
                return this._pade(args[0], args[1], args[2], args[3]);
            }

            if (node.funcName === 'vandermonde') {
                if (node.args.length !== 1) throw new Error("vandermonde requires 1 argument (vector)");
                return this._vandermonde(args[0]);
            }

            if (node.funcName === 'hilbert') {
                if (node.args.length !== 1) throw new Error("hilbert requires 1 argument (size)");
                return this._hilbert(args[0]);
            }

            if (node.funcName === 'toeplitz') {
                if (node.args.length < 1) throw new Error("toeplitz requires at least 1 argument");
                const c = args[0];
                const r = args.length > 1 ? args[1] : null;
                return this._toeplitz(c, r);
            }

            if (node.funcName === 'wronskian') {
                if (node.args.length !== 2) throw new Error("wronskian requires 2 arguments: list_of_funcs, var");
                return this._wronskian(args[0], args[1]);
            }

            return new Call(node.funcName, args);
        }

        if (node instanceof BinaryOp) {
            const left = this._recursiveEval(node.left);
            const right = this._recursiveEval(node.right);

            // Special handling for Matrix Power with negative exponent (Inverse)
            // or large powers using diagonalization if needed
            if (node instanceof Pow && left instanceof Vec && right instanceof Num) {
                if (right.value === -1) {
                    return this._inv(left);
                }
                // Forward to matrix pow logic if integer
                if (Number.isInteger(right.value)) {
                    // Let Pow.simplify handle small positive powers, or use _matrixPow for large?
                    // Currently Pow.simplify handles small positive integers via loop.
                    // For consistency, we can leave it or enhance.
                }
            }

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
            return this._productSymbolic(expr, varNode, start, end);
        }

        let prod = new Num(1);
        for (let i = start.value; i <= end.value; i++) {
            const term = expr.substitute(varNode, new Num(i)).simplify();
            prod = new Mul(prod, term).simplify();
        }
        return prod;
    }

    _productSymbolic(expr, varNode, start, end) {
        expr = expr.expand().simplify();

        // product(c, k, 1, n) = c^n
        if (!this._dependsOn(expr, varNode)) {
            const count = new Add(new Sub(end, start), new Num(1)).simplify();
            return new Pow(expr, count).simplify();
        }

        // Normalize start to 1
        // product(f(k), k, a, b) -> product(f(j+a-1), j, 1, b-a+1)
        if (!(start instanceof Num && start.value === 1)) {
             const j = new Sym(varNode.name + "_idx");
             const shift = new Sub(start, new Num(1));
             const subExpr = expr.substitute(varNode, new Add(j, shift)).simplify();
             const newEnd = new Sub(end, shift).simplify();
             return this._productSymbolic(subExpr, j, new Num(1), newEnd);
        }

        // Now start is 1. Number of terms is 'end' (denoted n)
        const n = end;

        // product(A*B) = product(A) * product(B)
        if (expr instanceof Mul) {
            return new Mul(
                this._productSymbolic(expr.left, varNode, start, end),
                this._productSymbolic(expr.right, varNode, start, end)
            ).simplify();
        }

        // product(A/B) = product(A) / product(B)
        if (expr instanceof Div) {
            return new Div(
                this._productSymbolic(expr.left, varNode, start, end),
                this._productSymbolic(expr.right, varNode, start, end)
            ).simplify();
        }

        // c^k -> c^(sum k)
        if (expr instanceof Pow && !this._dependsOn(expr.left, varNode)) {
             const exponentSum = this._sumSymbolic(expr.right, varNode, start, end);
             return new Pow(expr.left, exponentSum).simplify();
        }

        // Check for Linear Term: c1*k + c0
        // product(c1*k + c0) = c1^n * gamma(n + 1 + c0/c1) / gamma(1 + c0/c1)
        try {
            const poly = this._getPolyCoeffs(expr, varNode);
            if (poly && poly.maxDeg === 1) {
                const c1 = poly.coeffs[1];
                const c0 = poly.coeffs[0] || new Num(0);

                // If c1 = 1: gamma(n + 1 + c0) / gamma(1 + c0)
                // If c1 != 1: c1^n * ...
                const B = new Div(c0, c1).simplify(); // c0/c1

                const gammaNum = new Call('gamma', [new Add(new Add(n, new Num(1)), B)]);
                const gammaDen = new Call('gamma', [new Add(new Num(1), B)]);
                const termGamma = new Div(gammaNum, gammaDen);

                if (c1 instanceof Num && c1.value === 1) {
                    return termGamma.simplify();
                } else {
                    const termC1 = new Pow(c1, n);
                    return new Mul(termC1, termGamma).simplify();
                }
            }
        } catch(e) {}

        return new Call("product", [expr, varNode, start, end]);
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

    _solveInequality(ineq, varNode) {
        // Normalize to P(x) > 0 (or >=, <, <=)
        // Move everything to Left Side: LHS - RHS
        const expr = new Sub(ineq.left, ineq.right).simplify();
        let type = ''; // gt, ge, lt, le

        if (ineq instanceof Gt) type = 'gt';
        else if (ineq instanceof Ge) type = 'ge';
        else if (ineq instanceof Lt) type = 'lt';
        else if (ineq instanceof Le) type = 'le';

        // 1. Try Linear Logic: ax + b > 0
        const poly = this._getPolyCoeffs(expr, varNode);
        if (poly && poly.maxDeg === 1) {
            const a = poly.coeffs[1];
            const b = poly.coeffs[0] || new Num(0);

            // Inequality: ax + b > 0  =>  ax > -b
            const negB = new Mul(new Num(-1), b).simplify();

            // We need to know the sign of 'a' to flip operator
            const aVal = a.evaluateNumeric();
            if (!isNaN(aVal)) {
                const rhs = new Div(negB, a).simplify();
                if (aVal > 0) {
                    // Sign preserved
                    if (type === 'gt') return new Gt(varNode, rhs);
                    if (type === 'ge') return new Ge(varNode, rhs);
                    if (type === 'lt') return new Lt(varNode, rhs);
                    if (type === 'le') return new Le(varNode, rhs);
                } else if (aVal < 0) {
                    // Sign flipped
                    if (type === 'gt') return new Lt(varNode, rhs);
                    if (type === 'ge') return new Le(varNode, rhs);
                    if (type === 'lt') return new Gt(varNode, rhs);
                    if (type === 'le') return new Ge(varNode, rhs);
                } else {
                    // a = 0. b > 0?
                    // Check logic of 0*x + b > 0 => b > 0
                    const bVal = b.evaluateNumeric();
                    if (!isNaN(bVal)) {
                        let holds = false;
                        if (type === 'gt') holds = bVal > 0;
                        if (type === 'ge') holds = bVal >= 0;
                        if (type === 'lt') holds = bVal < 0;
                        if (type === 'le') holds = bVal <= 0;
                        return holds ? new Sym('true') : new Sym('false');
                    }
                }
            } else {
                // Symbolic 'a'. We can't determine direction.
                // Return undefined or just the simplified form?
                // Standard CAS might return piecewise or warning.
                // Let's return strict form ax > -b (assuming positive?) No, dangerous.
                return new Call('solve', [ineq, varNode]);
            }
        }

        // 2. Polynomial / Rational Logic via Roots and Poles
        // Zeros: Roots of P(x) = 0
        // Poles: Roots of Q(x) = 0 where expr = P/Q
        // Critical Points = Zeros U Poles

        // Find Zeros
        const zerosRes = this._solve(expr, varNode);
        let zeros = [];
        if (zerosRes instanceof Call && zerosRes.funcName === 'set') {
            zeros = zerosRes.args;
        } else if (zerosRes instanceof Expr && !(zerosRes instanceof Call && zerosRes.funcName === 'solve')) {
            zeros = [zerosRes];
        }

        // Find Poles (Singularities)
        let poles = [];
        if (expr instanceof Div) {
            const polesRes = this._solve(expr.right, varNode);
            if (polesRes instanceof Call && polesRes.funcName === 'set') {
                poles = polesRes.args;
            } else if (polesRes instanceof Expr && !(polesRes instanceof Call && polesRes.funcName === 'solve')) {
                poles = [polesRes];
            }
        }

        // Collect all critical points with type
        const criticalPoints = [];

        for (const r of zeros) {
            const val = r.evaluateNumeric();
            if (!isNaN(val)) criticalPoints.push({ val: val, node: r, type: 'zero' });
        }
        for (const r of poles) {
            const val = r.evaluateNumeric();
            if (!isNaN(val)) criticalPoints.push({ val: val, node: r, type: 'pole' });
        }

        // Sort
        criticalPoints.sort((a, b) => a.val - b.val);

        // Deduplicate
        const uniquePoints = [];
        if (criticalPoints.length > 0) {
            uniquePoints.push(criticalPoints[0]);
            for(let i=1; i<criticalPoints.length; i++) {
                const diff = Math.abs(criticalPoints[i].val - uniquePoints[uniquePoints.length-1].val);
                if (diff < 1e-9) {
                    // Duplicate/Collision
                    // Pole overrides Zero
                    if (criticalPoints[i].type === 'pole') {
                        uniquePoints[uniquePoints.length-1] = criticalPoints[i];
                    }
                } else {
                    uniquePoints.push(criticalPoints[i]);
                }
            }
        }

        if (uniquePoints.length === 0) {
            // No roots (or only complex/unknown).
            // Test one point (0 or any)
            const testVal = expr.substitute(varNode, new Num(0)).evaluateNumeric();
            if (isNaN(testVal)) return new Call('solve', [ineq, varNode]);

            let holds = false;
            if (type === 'gt') holds = testVal > 0;
            if (type === 'ge') holds = testVal >= 0;
            if (type === 'lt') holds = testVal < 0;
            if (type === 'le') holds = testVal <= 0;

            return holds ? new Sym('true') : new Sym('false'); // Or Empty set / All reals
        }

        // Handle case where we have points but they are all poles/zeros.
        // If the only points are poles, we still test intervals.
        // But what if `poles` are not found because `expr` structure isn't `Div`?
        // Fallback check for Mul(..., Pow(..., -1))
        if (poles.length === 0 && expr instanceof Mul) {
             const findDenom = (node) => {
                 if (node instanceof Pow && node.right instanceof Num && node.right.value < 0) return node.left;
                 if (node instanceof Mul) {
                     const d1 = findDenom(node.left);
                     if (d1) return d1;
                     return findDenom(node.right);
                 }
                 return null;
             };
             const den = findDenom(expr);
             if (den) {
                 const polesRes = this._solve(den, varNode);
                 let morePoles = [];
                 if (polesRes instanceof Call && polesRes.funcName === 'set') {
                     morePoles = polesRes.args;
                 } else if (polesRes instanceof Expr && !(polesRes instanceof Call && polesRes.funcName === 'solve')) {
                     morePoles = [polesRes];
                 }
                 for (const r of morePoles) {
                     const val = r.evaluateNumeric();
                     if (!isNaN(val)) {
                         // Add to uniquePoints and re-sort
                         uniquePoints.push({ val: val, node: r, type: 'pole' });
                     }
                 }
                 uniquePoints.sort((a, b) => a.val - b.val);
                 // Dedupe again? (Simplified logic: assume adding poles didn't duplicate zeros exactly, or sort handles order)
             }
        }

        // Test Intervals
        // Intervals: (-inf, p0), (p0, p1), ..., (pn, inf)
        const intervals = [];
        // (-inf, p0)
        intervals.push({
            check: uniquePoints[0].val - 1,
            cond: (x) => new Lt(x, uniquePoints[0].node)
        });

        for(let i=0; i<uniquePoints.length - 1; i++) {
            const mid = (uniquePoints[i].val + uniquePoints[i+1].val) / 2;
            intervals.push({
                check: mid,
                cond: (x) => new And(new Gt(x, uniquePoints[i].node), new Lt(x, uniquePoints[i+1].node))
            });
        }

        // (pn, inf)
        intervals.push({
            check: uniquePoints[uniquePoints.length - 1].val + 1,
            cond: (x) => new Gt(x, uniquePoints[uniquePoints.length - 1].node)
        });

        // Collect valid intervals
        const validConds = [];
        for(const interval of intervals) {
            const testVal = expr.substitute(varNode, new Num(interval.check)).evaluateNumeric();
            let holds = false;
            if (type === 'gt') holds = testVal > 0;
            if (type === 'ge') holds = testVal >= 0;
            if (type === 'lt') holds = testVal < 0;
            if (type === 'le') holds = testVal <= 0;

            if (holds) {
                validConds.push(interval.cond(varNode));
            }
        }

        // Handle Equality at roots for >= and <=
        if (type === 'ge' || type === 'le') {
            // Only include points that are 'zero', not 'pole'
            for(const pt of uniquePoints) {
                if (pt.type === 'zero') {
                    validConds.push(new Eq(varNode, pt.node));
                }
            }
        }

        if (validConds.length === 0) return new Sym('false'); // No solution
        if (validConds.length === 1) return validConds[0];

        // Combine with OR
        let res = validConds[0];
        for(let i=1; i<validConds.length; i++) {
            res = new Or(res, validConds[i]);
        }
        return res;
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

    _implicitDiff(eq, depVar, indepVar) {
        if (!(depVar instanceof Sym) || !(indepVar instanceof Sym)) throw new Error("Variables must be symbols");

        // Normalize F(x, y) = 0
        let F = eq;
        if (eq instanceof Eq) {
            F = new Sub(eq.left, eq.right).simplify();
        }

        // dy/dx = - (dF/dx) / (dF/dy)
        // Partial derivatives
        const Fx = F.diff(indepVar).simplify();
        const Fy = F.diff(depVar).simplify();

        return new Mul(new Num(-1), new Div(Fx, Fy)).simplify();
    }

    _minimize(expr, varNode) {
        return this._optimize(expr, varNode, 'min');
    }

    _maximize(expr, varNode) {
        return this._optimize(expr, varNode, 'max');
    }

    _analyze(expr, varNode) {
        const results = [];

        // Helper to format
        const label = (txt, val) => {
            return new Eq(new Sym(txt), val);
        };

        // 1. Function
        results.push(label("Function", expr));

        // 2. Roots (Zeros)
        try {
            const roots = this._roots(expr, varNode);
            results.push(label("Roots", roots));
        } catch(e) {
            results.push(label("Roots", new Sym("Error")));
        }

        // 3. Y-Intercept
        try {
            const yInt = expr.substitute(varNode, new Num(0)).simplify();
            results.push(label("Y_Intercept", yInt));
        } catch(e) { }

        // 4. Extrema
        try {
            const extr = this._extrema(expr, varNode);
            results.push(label("Extrema", extr));
        } catch(e) {}

        // 5. Asymptotes
        try {
            const asym = this._asymptotes(expr, varNode);
            results.push(label("Asymptotes", asym));
        } catch(e) {}

        // 6. Inflection Points candidates
        try {
            const deriv = expr.diff(varNode).simplify();
            const deriv2 = deriv.diff(varNode).simplify();
            const inf = this._roots(deriv2, varNode);
            results.push(label("Inflection_Pts", inf));
        } catch(e) {}

        return new Vec(results);
    }

    _extrema(expr, varNode) {
        // Find critical points: f'(x) = 0
        const deriv = expr.diff(varNode).simplify();
        const roots = this._roots(deriv, varNode);
        const points = (roots instanceof Vec) ? roots.elements : [roots];

        const deriv2 = deriv.diff(varNode).simplify();
        const results = [];

        for (const pt of points) {
            try {
                // Classify
                const val2 = deriv2.substitute(varNode, pt).evaluateNumeric();
                let type = "unknown";
                if (!isNaN(val2)) {
                    if (val2 > 0) type = "min";
                    else if (val2 < 0) type = "max";
                    else type = "saddle"; // or inflection
                }

                // Get value
                const val = expr.substitute(varNode, pt).simplify();
                results.push(new Vec([pt, val, new Sym(type)]));
            } catch(e) {}
        }
        return new Vec(results);
    }

    _stationary_points(expr, varNode) {
        const deriv = expr.diff(varNode).simplify();
        return this._roots(deriv, varNode);
    }

    _asymptotes(expr, varNode) {
        const results = [];

        // Vertical Asymptotes: Roots of denominator
        if (expr instanceof Div) {
            const den = expr.right;
            const poles = this._roots(den, varNode);
            const pList = (poles instanceof Vec) ? poles.elements : [poles];

            for(const p of pList) {
                // Check limit?
                results.push(new Eq(varNode, p));
            }
        }

        // Horizontal Asymptotes: Limit at Infinity
        const limInf = this._limit(expr, varNode, new Sym("Infinity"));
        if (!(limInf instanceof Sym && (limInf.name === "Infinity" || limInf.name === "NaN" || limInf.name === "-Infinity"))) {
             const y = new Sym('y');
             results.push(new Eq(y, limInf));
        }

        const limNegInf = this._limit(expr, varNode, new Mul(new Num(-1), new Sym("Infinity")));
        if (!(limNegInf instanceof Sym && (limNegInf.name === "Infinity" || limNegInf.name === "NaN" || limNegInf.name === "-Infinity"))) {
             // Avoid duplicate if same as limInf
             if (limNegInf.toString() !== limInf.toString()) {
                 const y = new Sym('y');
                 results.push(new Eq(y, limNegInf));
             }
        }

        return new Vec(results);
    }

    _completeSquare(expr, varNode) {
        expr = expr.expand().simplify();
        const poly = this._getPolyCoeffs(expr, varNode);

        if (!poly || poly.maxDeg !== 2) {
            throw new Error("completeSquare requires a quadratic polynomial");
        }

        // ax^2 + bx + c
        const a = poly.coeffs[2];
        const b = poly.coeffs[1] || new Num(0);
        const c = poly.coeffs[0] || new Num(0);

        // a(x + b/(2a))^2 + (c - b^2/4a)

        const term1_inner = new Add(varNode, new Div(b, new Mul(new Num(2), a)));
        const term1 = new Mul(a, new Pow(term1_inner, new Num(2)));

        const b2 = new Pow(b, new Num(2));
        const four_a = new Mul(new Num(4), a);
        const term2 = new Sub(c, new Div(b2, four_a));

        return new Add(term1, term2).simplify();
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
        // Multivariate Newton-Raphson
        if (eq instanceof Vec && varNode instanceof Vec && guess instanceof Vec) {
            return this._fsolveSystem(eq, varNode, guess);
        }

        // Single Variable Newton-Raphson Method
        // x_{n+1} = x_n - f(x_n) / f'(x_n)

        if (!(varNode instanceof Sym)) throw new Error("Numeric solve requires a single variable for single equation");

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

    _fsolveSystem(eqsVec, varsVec, guessVec) {
        const n = eqsVec.elements.length;
        if (varsVec.elements.length !== n || guessVec.elements.length !== n) {
            throw new Error("Dimension mismatch in fsolve system");
        }

        // Normalize equations f(x) = 0
        const functions = eqsVec.elements.map(e => {
            if (e instanceof Eq) return new Sub(e.left, e.right).simplify();
            return e.simplify();
        });

        const vars = varsVec.elements;
        let currentX = guessVec.elements.map(g => g.evaluateNumeric());
        if (currentX.some(v => isNaN(v))) throw new Error("Initial guess must be numeric");

        // Jacobian Matrix J[i][j] = df_i / dx_j
        // Pre-compute symbolic derivatives
        const jacobianSymbolic = [];
        for(let i=0; i<n; i++) {
            const row = [];
            for(let j=0; j<n; j++) {
                row.push(functions[i].diff(vars[j]).simplify());
            }
            jacobianSymbolic.push(row);
        }

        const maxIter = 100;
        const tol = 1e-9;

        for (let iter = 0; iter < maxIter; iter++) {
            // Evaluate F(x)
            const F = [];
            for(let i=0; i<n; i++) {
                let val = functions[i];
                for(let k=0; k<n; k++) val = val.substitute(vars[k], new Num(currentX[k]));
                F.push(val.evaluateNumeric());
            }

            // Check convergence (norm of F)
            const normF = Math.sqrt(F.reduce((sum, v) => sum + v*v, 0));
            if (normF < tol) {
                return new Vec(currentX.map(v => new Num(v)));
            }

            // Evaluate Jacobian J(x)
            const J = [];
            for(let i=0; i<n; i++) {
                const row = [];
                for(let j=0; j<n; j++) {
                    let val = jacobianSymbolic[i][j];
                    for(let k=0; k<n; k++) val = val.substitute(vars[k], new Num(currentX[k]));
                    row.push(val.evaluateNumeric());
                }
                J.push(row);
            }

            // Solve J * deltaX = -F
            // Use Gaussian elimination for numerical system
            // Augmented matrix [J | -F]
            const M = [];
            for(let i=0; i<n; i++) {
                const row = [...J[i], -F[i]];
                M.push(row);
            }

            // Gaussian elimination
            for(let i=0; i<n; i++) {
                // Pivot
                let maxRow = i;
                for(let k=i+1; k<n; k++) {
                    if (Math.abs(M[k][i]) > Math.abs(M[maxRow][i])) maxRow = k;
                }
                [M[i], M[maxRow]] = [M[maxRow], M[i]];

                if (Math.abs(M[i][i]) < 1e-12) throw new Error("Jacobian is singular");

                for(let k=i+1; k<n; k++) {
                    const factor = M[k][i] / M[i][i];
                    for(let j=i; j<=n; j++) {
                        M[k][j] -= factor * M[i][j];
                    }
                }
            }

            // Back substitution
            const deltaX = new Array(n).fill(0);
            for(let i=n-1; i>=0; i--) {
                let sum = M[i][n];
                for(let j=i+1; j<n; j++) {
                    sum -= M[i][j] * deltaX[j];
                }
                deltaX[i] = sum / M[i][i];
            }

            // Update currentX
            for(let i=0; i<n; i++) {
                currentX[i] += deltaX[i];
            }

            // Check step size
            const stepSize = Math.sqrt(deltaX.reduce((sum, v) => sum + v*v, 0));
            if (stepSize < tol) {
                return new Vec(currentX.map(v => new Num(v)));
            }
        }

        throw new Error("System solver did not converge");
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

    _roots(poly, varNode) {
        const sol = this._solve(poly, varNode);
        const values = [];

        const process = (s) => {
            if (s instanceof Eq) values.push(s.right);
            else if (s instanceof Expr && !(s instanceof Call && s.funcName === 'solve')) {
                values.push(s);
            }
        };

        if (sol instanceof Call && sol.funcName === 'set') {
            sol.args.forEach(process);
        } else {
            process(sol);
        }

        return new Vec(values);
    }

    _solve(eq, varNode) {
        if (eq instanceof Vec && varNode instanceof Vec) {
            return this._solveSystem(eq, varNode);
        }

        // Handle Inequalities
        if (eq instanceof Lt || eq instanceof Gt || eq instanceof Le || eq instanceof Ge) {
            return this._solveInequality(eq, varNode);
        }

        // Handle Abs(A) = B  => A = B or A = -B
        // Check if eq is Eq(abs(A), B) or abs(A) = B (implicitly A-B=0 if input was expression)
        // If input was Eq, we have explicit Left/Right.
        if (eq instanceof Eq) {
            if (eq.left instanceof Call && eq.left.funcName === 'abs') {
                const A = eq.left.args[0];
                const B = eq.right;
                if (!this._dependsOn(B, varNode)) { // B should be constant w.r.t var (or at least handle splitting)
                    const eq1 = new Eq(A, B);
                    const eq2 = new Eq(A, new Mul(new Num(-1), B));
                    const sol1 = this._solve(eq1, varNode);
                    const sol2 = this._solve(eq2, varNode);
                    // Merge solutions
                    const sols = [];
                    const collect = (s) => {
                        if (s instanceof Call && s.funcName === 'set') s.args.forEach(collect);
                        else sols.push(s);
                    };
                    collect(sol1);
                    collect(sol2);
                    sols.sort((a, b) => {
                        const va = a.evaluateNumeric();
                        const vb = b.evaluateNumeric();
                        if (!isNaN(va) && !isNaN(vb)) return va - vb;
                        return a.toString().localeCompare(b.toString());
                    });
                    return new Call('set', sols);
                }
            }
             // B = abs(A)?
            if (eq.right instanceof Call && eq.right.funcName === 'abs') {
                 const A = eq.right.args[0];
                 const B = eq.left;
                 if (!this._dependsOn(B, varNode)) {
                    const eq1 = new Eq(A, B);
                    const eq2 = new Eq(A, new Mul(new Num(-1), B));
                    const sol1 = this._solve(eq1, varNode);
                    const sol2 = this._solve(eq2, varNode);
                    // Merge solutions
                    const sols = [];
                    const collect = (s) => {
                        if (s instanceof Call && s.funcName === 'set') s.args.forEach(collect);
                        else sols.push(s);
                    };
                    collect(sol1);
                    collect(sol2);
                    sols.sort((a, b) => {
                        const va = a.evaluateNumeric();
                        const vb = b.evaluateNumeric();
                        if (!isNaN(va) && !isNaN(vb)) return va - vb;
                        return a.toString().localeCompare(b.toString());
                    });
                    return new Call('set', sols);
                 }
            }
        } else {
            // Expression = 0
            // Check if expression contains abs
            // e.g. abs(x) - 5
            if (eq instanceof Sub) {
                if (eq.left instanceof Call && eq.left.funcName === 'abs') {
                    // abs(A) - B = 0 -> abs(A) = B
                    return this._solve(new Eq(eq.left, eq.right), varNode);
                }
            }
        }

        let expr;
        if (eq instanceof Eq) {
            expr = new Sub(eq.left, eq.right).simplify();
        } else {
            expr = eq.simplify();
        }

        // Handle Rational Equation P/Q = 0 => P = 0
        if (expr instanceof Div) {
            return this._solve(expr.left, varNode);
        }

        // Expand to ensure polynomial form for identification
        expr = expr.expand().simplify();

        // Re-check Div after expansion/simplification
        if (expr instanceof Div) {
            return this._solve(expr.left, varNode);
        }

        try {
            // Try factoring FIRST if degree > 2, before expanding?
            // No, expanding is needed to know degree.

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

                        const results = [sol1.simplify(), sol2.simplify()];
                        results.sort((a, b) => {
                            const va = a.evaluateNumeric();
                            const vb = b.evaluateNumeric();
                            if (!isNaN(va) && !isNaN(vb)) return va - vb;
                            return a.toString().localeCompare(b.toString());
                        });

                        return new Call("set", results);
                    }
                } else if (poly.maxDeg > 2) {
                    // Try factoring higher degree polynomials
                    // Note: _factor internally uses _factorRational which might use the same poly coeffs
                    const factored = this._factor(expr);

                    // If factorization successfully broke it down (not just returned the input wrapped)
                    // Important: compare strings or structure.
                    const isFactorCall = (factored instanceof Call && factored.funcName === 'factor');

                    if (!isFactorCall) {
                        const terms = [];
                        const collect = (n) => {
                            if (n instanceof Mul) {
                                collect(n.left);
                                collect(n.right);
                            } else {
                                terms.push(n);
                            }
                        };

                        // Handle 'factored' call wrapper if it exists (for integers, but maybe mixed)
                        if (factored instanceof Call && factored.funcName === 'factored') {
                            factored.args.forEach(a => terms.push(a));
                        } else {
                            collect(factored);
                        }

                        if (terms.length >= 1) {
                            let solutions = [];

                            for(const term of terms) {
                                // term could be Pow(base, exp)
                                let target = term;
                                if (term instanceof Pow) target = term.left;

                                // Solve factor = 0
                                if (this._dependsOn(target, varNode)) {
                                    const sol = this._solve(target, varNode);

                                    if (sol instanceof Call && sol.funcName === 'set') {
                                        sol.args.forEach(s => solutions.push(s));
                                    } else if (sol instanceof Expr && !(sol instanceof Call && sol.funcName === 'solve')) {
                                        solutions.push(sol);
                                    }
                                }
                            }

                            if (solutions.length > 0) {
                                // Deduplicate
                                const unique = [];
                                const seen = new Set();
                                solutions.forEach(s => {
                                    const str = s.toString();
                                    if (!seen.has(str)) {
                                        seen.add(str);
                                        unique.push(s);
                                    }
                                });

                                if (unique.length === 1) return unique[0];
                                // Sort solutions if possible
                                unique.sort((a, b) => {
                                    const va = a.evaluateNumeric();
                                    const vb = b.evaluateNumeric();
                                    if (!isNaN(va) && !isNaN(vb)) return va - vb;
                                    return a.toString().localeCompare(b.toString());
                                });
                                return new Call("set", unique);
                            }
                        }
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

            // Attempt Trigonometric Solving
            // sin(u) = c, cos(u) = c, tan(u) = c
            // expr is LHS-RHS, so sin(u) - c = 0 => sin(u) = c
            // We look for form Func(Arg) - Const or Const - Func(Arg)
            let funcTerm = null;
            let constTerm = new Num(0);
            let sign = 1;

            if (expr instanceof Sub) {
                if (expr.left instanceof Call && !this._dependsOn(expr.right, varNode)) {
                    funcTerm = expr.left;
                    constTerm = expr.right;
                } else if (expr.right instanceof Call && !this._dependsOn(expr.left, varNode)) {
                    funcTerm = expr.right;
                    constTerm = expr.left; // Const - Func = 0 => Func = Const
                }
            } else if (expr instanceof Add) {
                if (expr.left instanceof Call && !this._dependsOn(expr.right, varNode)) {
                    funcTerm = expr.left;
                    constTerm = new Mul(new Num(-1), expr.right).simplify(); // Func = -Const
                } else if (expr.right instanceof Call && !this._dependsOn(expr.left, varNode)) {
                    funcTerm = expr.right;
                    constTerm = new Mul(new Num(-1), expr.left).simplify();
                }
            } else if (expr instanceof Call) {
                funcTerm = expr; // Func = 0
            }

            if (funcTerm && (funcTerm.funcName === 'sin' || funcTerm.funcName === 'cos' || funcTerm.funcName === 'tan')) {
                const arg = funcTerm.args[0];
                if (this._dependsOn(arg, varNode)) {
                    // sin(arg) = C
                    const C = constTerm;
                    let solutions = [];
                    const pi = new Sym('pi');

                    if (funcTerm.funcName === 'sin') {
                        // arg = asin(C)
                        // arg = pi - asin(C)
                        const asinC = new Call('asin', [C]).simplify();
                        const sol1 = this._solve(new Eq(arg, asinC), varNode);
                        const sol2 = this._solve(new Eq(arg, new Sub(pi, asinC)), varNode);
                        solutions.push(sol1);
                        solutions.push(sol2);
                    } else if (funcTerm.funcName === 'cos') {
                        // arg = acos(C)
                        // arg = -acos(C)
                        const acosC = new Call('acos', [C]).simplify();
                        const sol1 = this._solve(new Eq(arg, acosC), varNode);
                        const sol2 = this._solve(new Eq(arg, new Mul(new Num(-1), acosC)), varNode);
                        solutions.push(sol1);
                        solutions.push(sol2);
                    } else if (funcTerm.funcName === 'tan') {
                        // arg = atan(C)
                        const atanC = new Call('atan', [C]).simplify();
                        const sol1 = this._solve(new Eq(arg, atanC), varNode);
                        solutions.push(sol1);
                    }

                    // Flatten solutions if they returned sets
                    const flatSols = [];
                    const process = (s) => {
                        if (s instanceof Call && s.funcName === 'set') s.args.forEach(process);
                        else if (s instanceof Expr && !(s instanceof Call && s.funcName === 'solve')) flatSols.push(s);
                    };
                    solutions.forEach(process);

                    if (flatSols.length > 0) {
                        if (flatSols.length === 1) return flatSols[0];
                        return new Call('set', flatSols);
                    }
                }
            }

        } catch (e) {
            console.error(e);
        }

        return new Call("solve", [eq, varNode]);
    }

    _solveInequality(ineq, varNode) {
        // Normalize to expr > 0 or < 0
        const expr = new Sub(ineq.left, ineq.right).simplify();

        // Find roots of expr = 0
        const rootsRes = this._solve(expr, varNode);

        // Check if solver failed (returned symbolic call)
        if (rootsRes instanceof Call && rootsRes.funcName === 'solve') {
            return new Call('solve', [ineq, varNode]);
        }

        let roots = [];
        if (rootsRes instanceof Call && rootsRes.funcName === 'set') {
            roots = rootsRes.args;
        } else if (rootsRes instanceof Expr) {
            roots = [rootsRes];
        }

        // Filter roots to real numeric values
        const realRoots = [];
        for(const r of roots) {
            const val = r.evaluateNumeric();
            if (!isNaN(val)) realRoots.push({ val: val, node: r });
        }

        if (realRoots.length === 0 && roots.length > 0) {
             // Symbolic roots only and we can't evaluate? Cannot reliably solve inequality.
             // Return symbolic call
             return new Call('solve', [ineq, varNode]);
        }

        // Sort
        realRoots.sort((a, b) => a.val - b.val);

        // Test Intervals
        const intervals = [];
        const points = [];

        if (realRoots.length === 0) {
            points.push({ val: 0, desc: 'all' }); // Test 0
        } else {
            // Pick safe test points avoiding roots
            points.push({ val: realRoots[0].val - 1, desc: '(-inf, r0)' });
            for(let i=0; i<realRoots.length - 1; i++) {
                points.push({ val: (realRoots[i].val + realRoots[i+1].val)/2, desc: 'mid' });
            }
            points.push({ val: realRoots[realRoots.length-1].val + 1, desc: '(inf)' });
        }

        const validIntervals = [];

        const check = (val) => {
            try {
                const res = expr.substitute(varNode, new Num(val)).evaluateNumeric();
                if (isNaN(res)) return false; // Undefined
                if (ineq instanceof Lt) return res < 0;
                if (ineq instanceof Gt) return res > 0;
                if (ineq instanceof Le) return res <= 0;
                if (ineq instanceof Ge) return res >= 0;
            } catch(e) {}
            return false;
        };

        // Helper to check boundaries for non-strict
        // If Le or Ge, we must include roots where expr=0.
        // Roots are where expr=0.
        // My logic below constructs intervals ( ) and joins them.
        // If it's non-strict, we should output [ ] intervals.

        for(let i=0; i<points.length; i++) {
            if (check(points[i].val)) {
                // Determine interval bounds
                // i corresponds to interval before roots[i] (if i < roots.length)
                // Actually my points array logic:
                // points[0] is before roots[0]
                // points[1] is between roots[0] and roots[1]
                // ...
                // points[k] is after roots[k-1] (last root)

                const strict = (ineq instanceof Lt || ineq instanceof Gt);

                if (realRoots.length === 0) {
                    // All reals satisfied
                    // Return "true"? Or (-inf, inf)
                    // Let's return logic: varNode > -inf and varNode < inf?
                    // Just 1 (true)
                    return new Num(1);
                }

                let lower, upper;
                let cond;

                if (i === 0) {
                    // x < root[0]
                    // Strict or non-strict depends on operator
                    cond = strict ? new Lt(varNode, realRoots[0].node) : new Le(varNode, realRoots[0].node);
                } else if (i === points.length - 1) {
                    // x > root[n-1]
                    cond = strict ? new Gt(varNode, realRoots[i-1].node) : new Ge(varNode, realRoots[i-1].node);
                } else {
                    // root[i-1] < x < root[i]
                    const c1 = strict ? new Gt(varNode, realRoots[i-1].node) : new Ge(varNode, realRoots[i-1].node);
                    const c2 = strict ? new Lt(varNode, realRoots[i].node) : new Le(varNode, realRoots[i].node);
                    cond = new And(c1, c2);
                }
                validIntervals.push(cond);
            }
        }

        if (validIntervals.length === 0) return new Num(0); // False
        if (validIntervals.length === 1) return validIntervals[0];

        // Join with Or
        let res = validIntervals[0];
        for(let k=1; k<validIntervals.length; k++) {
            res = new Or(res, validIntervals[k]);
        }
        return res;
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
        if (n < 0) return 0; // Or handle negative logic F(-n) = (-1)^(n+1)F(n)
        if (n === 0) return 0;
        if (n === 1) return 1;
        if (n === 2) return 1;

        // O(log n) Fast Doubling
        // F(2k) = F(k) * (2*F(k+1) - F(k))
        // F(2k+1) = F(k+1)^2 + F(k)^2

        // We compute pair (F(k), F(k+1))
        const fib = (k) => {
            if (k === 0) return [0, 1]; // F(0), F(1)
            const p = Math.floor(k / 2);
            const [fk, fk1] = fib(p);
            // c = F(2p) = F(p)*(2F(p+1) - F(p))
            const c = fk * (2 * fk1 - fk);
            // d = F(2p+1) = F(p+1)^2 + F(p)^2
            const d = fk1 * fk1 + fk * fk;

            if (k % 2 === 0) {
                return [c, d];
            } else {
                return [d, c + d];
            }
        };

        const [fn, fn1] = fib(n);
        return fn;
    }

    _gamma(z) {
        if (z < 0.5) {
            // Reflection formula for negative z or z < 0.5
            // Gamma(z) * Gamma(1-z) = pi / sin(pi*z)
            // Use this directly instead of logGamma to preserve sign for negative inputs
            return Math.PI / (Math.sin(Math.PI * z) * this._gamma(1 - z));
        }
        return Math.exp(this._logGamma(z));
    }

    _logGamma(z) {
        // Lanczos approximation for z >= 0.5
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

        if (z < 0.5) {
            // We shouldn't really use logGamma for negative z if we want the sign?
            // But if user calls logGamma explicitly?
            // log(Gamma(z)). If Gamma(z) < 0, this is complex.
            // For now, assume this is used for probability where z > 0 usually.
            // But if z is just small positive (< 0.5), we use reflection.
            return Math.log(Math.PI / (Math.sin(Math.PI * z) * this._gamma(1 - z)));
        }

        z -= 1;
        let x = C[0];
        for (let i = 1; i < g + 2; i++) x += C[i] / (z + i);

        const t = z + g + 0.5;
        return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
    }

    _gammainc(s, x) {
        // Regularized Incomplete Gamma P(s, x)
        // Uses Series for x < s+1, Continued Fraction for x >= s+1
        if (x < 0) return 0;
        if (x === 0) return 0;
        if (s <= 0) return 0; // Undefined usually

        if (x < s + 1.0) {
            // Series
            let ap = s;
            let sum = 1.0 / s;
            let del = sum;
            for (let n = 1; n < 100; n++) {
                ap += 1.0;
                del *= x / ap;
                sum += del;
                if (Math.abs(del) < Math.abs(sum) * 1e-14) {
                    return sum * Math.exp(-x + s * Math.log(x) - this._logGamma(s));
                }
            }
        } else {
            // Continued Fraction
            const gln = this._logGamma(s);
            let b = x + 1.0 - s;
            let c = 1.0 / 1e-30;
            let d = 1.0 / b;
            let h = d;
            for (let i = 1; i < 100; i++) {
                const an = -i * (i - s);
                b += 2.0;
                d = an * d + b;
                if (Math.abs(d) < 1e-30) d = 1e-30;
                c = b + an / c;
                if (Math.abs(c) < 1e-30) c = 1e-30;
                d = 1.0 / d;
                const del = d * c;
                h *= del;
                if (Math.abs(del - 1.0) < 1e-14) {
                    return 1.0 - Math.exp(-x + s * Math.log(x) - gln) * h;
                }
            }
        }
        return 0; // Failed
    }

    _betaNumeric(a, b) {
        // B(a, b) = Gamma(a)Gamma(b) / Gamma(a+b)
        return Math.exp(this._logGamma(a) + this._logGamma(b) - this._logGamma(a + b));
    }

    _betainc(x, a, b) {
        // Regularized Incomplete Beta I_x(a, b)
        if (x < 0 || x > 1) return 0; // Or NaN
        if (x === 0) return 0;
        if (x === 1) return 1;

        const bt = Math.exp(this._logGamma(a + b) - this._logGamma(a) - this._logGamma(b) + a * Math.log(x) + b * Math.log(1.0 - x));

        if (x < (a + 1.0) / (a + b + 2.0)) {
            return bt * this._betacf(x, a, b) / a;
        } else {
            return 1.0 - bt * this._betacf(1.0 - x, b, a) / b;
        }
    }

    _betacf(x, a, b) {
        // Continued Fraction for Beta
        const MAXIT = 100;
        const EPS = 1e-14;
        let qab = a + b;
        let qap = a + 1.0;
        let qam = a - 1.0;
        let c = 1.0;
        let d = 1.0 - qab * x / qap;
        if (Math.abs(d) < 1e-30) d = 1e-30;
        d = 1.0 / d;
        let h = d;

        for (let m = 1; m <= MAXIT; m++) {
            let m2 = 2 * m;
            let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1.0 + aa * d;
            if (Math.abs(d) < 1e-30) d = 1e-30;
            c = 1.0 + aa / c;
            if (Math.abs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d;
            h *= d * c;
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1.0 + aa * d;
            if (Math.abs(d) < 1e-30) d = 1e-30;
            c = 1.0 + aa / c;
            if (Math.abs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d;
            const del = d * c;
            h *= del;
            if (Math.abs(del - 1.0) < EPS) break;
        }
        return h;
    }

    _integrateTrig(expr, varNode) {
        if (expr instanceof Call) {
            // tan(x) -> -ln(cos(x))
            if (expr.funcName === 'tan' && expr.args.length === 1 && expr.args[0].toString() === varNode.toString()) {
                return new Mul(new Num(-1), new Call('ln', [new Call('cos', [varNode])])).simplify();
            }
            // cot(x) -> ln(sin(x))
            if (expr.funcName === 'cot' && expr.args.length === 1 && expr.args[0].toString() === varNode.toString()) {
                return new Call('ln', [new Call('sin', [varNode])]).simplify();
            }
            // sec(x) -> ln(sec(x) + tan(x))
            if (expr.funcName === 'sec' && expr.args.length === 1 && expr.args[0].toString() === varNode.toString()) {
                const sum = new Add(new Call('sec', [varNode]), new Call('tan', [varNode]));
                return new Call('ln', [sum]).simplify();
            }
            // csc(x) -> -ln(csc(x) + cot(x))
            if (expr.funcName === 'csc' && expr.args.length === 1 && expr.args[0].toString() === varNode.toString()) {
                const sum = new Add(new Call('csc', [varNode]), new Call('cot', [varNode]));
                return new Mul(new Num(-1), new Call('ln', [sum])).simplify();
            }
        }

        // Helper to match func(var)^2
        const isTrigSq = (node, func) => {
            if (node instanceof Pow && node.right instanceof Num && node.right.value === 2) {
                if (node.left instanceof Call && node.left.funcName === func) {
                    if (node.left.args.length === 1 && node.left.args[0].toString() === varNode.toString()) {
                        return true;
                    }
                }
            }
            return false;
        };

        // sec^2(x) -> tan(x)
        if (isTrigSq(expr, 'sec')) {
            return new Call('tan', [varNode]);
        }

        // csc^2(x) -> -cot(x)
        if (isTrigSq(expr, 'csc')) {
            return new Mul(new Num(-1), new Call('cot', [varNode])).simplify();
        }

        // sec(x)tan(x) -> sec(x)
        // csc(x)cot(x) -> -csc(x)
        if (expr instanceof Mul) {
            const l = expr.left;
            const r = expr.right;
            const isCall = (n, name) => (n instanceof Call && n.funcName === name && n.args.length === 1 && n.args[0].toString() === varNode.toString());

            if ((isCall(l, 'sec') && isCall(r, 'tan')) || (isCall(l, 'tan') && isCall(r, 'sec'))) {
                return new Call('sec', [varNode]);
            }

            if ((isCall(l, 'csc') && isCall(r, 'cot')) || (isCall(l, 'cot') && isCall(r, 'csc'))) {
                return new Mul(new Num(-1), new Call('csc', [varNode])).simplify();
            }
        }

        // 1/cos(x)^2 -> tan(x)
        if (expr instanceof Div && expr.left instanceof Num && expr.left.value === 1) {
            if (isTrigSq(expr.right, 'cos')) {
                return new Call('tan', [varNode]);
            }
            if (isTrigSq(expr.right, 'sin')) {
                return new Mul(new Num(-1), new Call('cot', [varNode])).simplify();
            }
        }
        if (expr instanceof Pow && expr.left instanceof Call && expr.left.funcName === 'sec' && expr.right instanceof Num && expr.right.value === 2) {
             if (expr.left.args[0].toString() === varNode.toString()) return new Call('tan', [varNode]);
        }
        // csc^2(x) -> -cot(x)
        if (expr instanceof Pow && expr.left instanceof Call && expr.left.funcName === 'csc' && expr.right instanceof Num && expr.right.value === 2) {
             if (expr.left.args[0].toString() === varNode.toString()) return new Mul(new Num(-1), new Call('cot', [varNode]));
        }
        // sec(x)*tan(x) -> sec(x)
        if (expr instanceof Mul) {
             const check = (func1, func2) => {
                 if (expr.left instanceof Call && expr.left.funcName === func1 && expr.right instanceof Call && expr.right.funcName === func2) {
                     if (expr.left.args[0].toString() === varNode.toString() && expr.right.args[0].toString() === varNode.toString()) return true;
                 }
                 if (expr.right instanceof Call && expr.right.funcName === func1 && expr.left instanceof Call && expr.left.funcName === func2) {
                     if (expr.right.args[0].toString() === varNode.toString() && expr.left.args[0].toString() === varNode.toString()) return true;
                 }
                 return false;
             };

             if (check('sec', 'tan')) return new Call('sec', [varNode]);
             if (check('csc', 'cot')) return new Mul(new Num(-1), new Call('csc', [varNode]));
        }

        return null;
    }

    _integrateAtan(expr, varNode) {
        // Look for 1/(Ax^2 + Bx + C)
        // expr should be Div(..., ...) or Mul(Div..., ...)

        let num, den;
        if (expr instanceof Div) {
            num = expr.left;
            den = expr.right;
        } else if (expr instanceof Mul) {
            // c * 1/den or num/den
            // This case is tricky without fully analyzing Mul structure
            // Let's assume Div case first or simple Mul(c, Div(1, den))
            // But simplification usually puts coeff in numerator if possible.
            // If expr is Mul(c, Div(1, Q)), let's treat it as Div(c, Q)
            if (expr.right instanceof Div && expr.right.left instanceof Num && expr.right.left.value === 1) {
                num = expr.left;
                den = expr.right.right;
            } else if (expr.left instanceof Div && expr.left.left instanceof Num && expr.left.left.value === 1) {
                num = expr.right;
                den = expr.left.right;
            } else {
                return null;
            }
        } else {
            return null;
        }

        if (this._dependsOn(num, varNode)) return null;

        // Check denominator for Ax^2 + Bx + C form
        const poly = this._getPolyCoeffs(den, varNode);
        if (poly && poly.maxDeg === 2) {
             const A = poly.coeffs[2];
             const B = poly.coeffs[1] || new Num(0);
             const C = poly.coeffs[0] || new Num(0);

             // Completing the square
             // Ax^2 + Bx + C = A(x + B/2A)^2 + (C - B^2/4A)
             // Let u = x + B/2A
             // K = C - B^2/4A
             // Integral 1 / (A(u^2 + K/A))
             // If K/A > 0, it is atan.
             // K/A = (4AC - B^2) / 4A^2.
             // Sign depends on 4AC - B^2 (Discriminant of Q).
             // If D = B^2 - 4AC < 0 => 4AC - B^2 > 0 => Atan.
             // If D > 0 => Partial fractions (ln).

             const discriminant = new Sub(new Pow(B, new Num(2)), new Mul(new Num(4), new Mul(A, C))).simplify();
             const discriminantVal = discriminant.evaluateNumeric();

             // Check if 4AC - B^2 > 0 (i.e. discriminant < 0)
             let isAtan = false;
             let isAtanh = false;

             if (!isNaN(discriminantVal)) {
                 if (discriminantVal < 0) isAtan = true; // D < 0 -> 4AC > B^2
                 else if (discriminantVal > 0) isAtanh = true; // D > 0 -> 4AC < B^2 (real roots)
             } else {
                 // Symbolic check
                 if (B instanceof Num && B.value === 0) {
                     const AC = new Mul(A, C).simplify();
                     const isNeg = (node) => {
                         if (node instanceof Num && node.value < 0) return true;
                         if (node instanceof Mul && node.left instanceof Num && node.left.value < 0) return true;
                         return false;
                     };
                     const ACVal = AC.evaluateNumeric();
                     if ((!isNaN(ACVal) && ACVal < 0) || (!isNaN(ACVal) ? false : isNeg(AC))) isAtanh = true;
                     else isAtan = true; // Assume positive (atan) if unknown
                 }
             }

             if (isAtan) {
                 // Formula: integral dx / (Ax^2 + Bx + C)
                 // let Delta = 4AC - B^2
                 // = 2/sqrt(Delta) * atan((2Ax + B)/sqrt(Delta))

                 const Delta = new Sub(new Mul(new Num(4), new Mul(A, C)), new Pow(B, new Num(2))).simplify();
                 const sqrtDelta = new Call('sqrt', [Delta]).simplify();

                 const argNum = new Add(new Mul(new Num(2), new Mul(A, varNode)), B).simplify();
                 const arg = new Div(argNum, sqrtDelta).simplify();

                 const coeff = new Div(new Num(2), sqrtDelta).simplify();

                 return new Mul(num, new Mul(coeff, new Call('atan', [arg]))).simplify();
             }

             if (isAtanh) {
                 // 1/(Ax^2+Bx+C) with D = B^2 - 4AC > 0.
                 // Formula: -2/sqrt(D) * atanh( (2Ax+B) / sqrt(D) )

                 const D = new Sub(new Pow(B, new Num(2)), new Mul(new Num(4), new Mul(A, C))).simplify();
                 const sqrtD = new Call('sqrt', [D]).simplify();

                 const argNum = new Add(new Mul(new Num(2), new Mul(A, varNode)), B).simplify();
                 const arg = new Div(argNum, sqrtD).simplify();

                 const coeff = new Div(new Num(-2), sqrtD).simplify();

                 return new Mul(num, new Mul(coeff, new Call('atanh', [arg]))).simplify();
             }
        }
        return null;
    }

    _integrateInverseHyperbolic(expr, varNode) {
        // Look for 1/sqrt(Ax^2 + Bx + C)
        // expr could be Div(1, sqrt(...)) or Div(c, sqrt(...)) or Mul(c, Div...)

        let num, den;
        if (expr instanceof Div) {
            num = expr.left;
            den = expr.right;
        } else if (expr instanceof Mul) {
             if (expr.right instanceof Div && expr.right.left instanceof Num && expr.right.left.value === 1) {
                num = expr.left;
                den = expr.right.right;
            } else if (expr.left instanceof Div && expr.left.left instanceof Num && expr.left.left.value === 1) {
                num = expr.right;
                den = expr.left.right;
            } else {
                return null;
            }
        } else {
            return null;
        }

        if (this._dependsOn(num, varNode)) return null;

        if (den instanceof Call && den.funcName === 'sqrt') {
            const inner = den.args[0];
            const poly = this._getPolyCoeffs(inner, varNode);

            if (poly && poly.maxDeg === 2) {
                const A = poly.coeffs[2];
                const B = poly.coeffs[1] || new Num(0);
                const C = poly.coeffs[0] || new Num(0);

                // Ax^2 + Bx + C
                // Cases based on sign of A and Discriminant D = B^2 - 4AC

                const AVal = A.evaluateNumeric();
                if (isNaN(AVal)) return null; // Need to know sign of A

                if (AVal > 0) {
                    // Form: 1/sqrt(A) * asinh( (2Ax+B) / sqrt(4AC-B^2) )  if 4AC - B^2 > 0
                    // Form: 1/sqrt(A) * acosh( (2Ax+B) / sqrt(B^2-4AC) )  if B^2 - 4AC > 0
                    // Or simpler: 1/sqrt(A) * ln( 2sqrt(A)*sqrt(Q) + 2Ax + B )

                    // Let's use the log form which covers both asinh/acosh usually
                    // Result = (1/sqrt(A)) * ln( 2*sqrt(A)*sqrt(Ax^2+Bx+C) + 2Ax + B )

                    const sqrtA = new Call('sqrt', [A]).simplify();

                    // Check if it simplifies to standard asinh/acosh first
                    if (B instanceof Num && B.value === 0) {
                        // Check if C is positive/negative (numerically or symbolically)
                        const isNeg = (node) => {
                            if (node instanceof Num && node.value < 0) return true;
                            if (node instanceof Mul && node.left instanceof Num && node.left.value < 0) return true;
                            return false;
                        };

                        const cVal = C.evaluateNumeric();
                        if ((!isNaN(cVal) && cVal > 0) || (!isNaN(cVal) ? false : !isNeg(C))) {
                            // 1/sqrt(Ax^2 + C) = 1/sqrt(A) * asinh( x * sqrt(A/C) )
                            const arg = new Mul(varNode, new Call('sqrt', [new Div(A, C)])).simplify();
                            return new Mul(new Div(num, sqrtA), new Call('asinh', [arg])).simplify();
                        } else {
                            // 1/sqrt(Ax^2 - |C|) = 1/sqrt(A) * acosh( x * sqrt(A/|C|) )
                            const negC = new Mul(new Num(-1), C).simplify();
                            const arg = new Mul(varNode, new Call('sqrt', [new Div(A, negC)])).simplify();
                            return new Mul(new Div(num, sqrtA), new Call('acosh', [arg])).simplify();
                        }
                    }

                    const term1 = new Mul(new Num(2), new Mul(sqrtA, den)).simplify(); // 2*sqrt(A)*sqrt(Q)
                    const term2 = new Add(new Mul(new Num(2), new Mul(A, varNode)), B).simplify(); // 2Ax + B

                    const arg = new Add(term1, term2).simplify();
                    const res = new Mul(new Div(num, sqrtA), new Call('ln', [arg])).simplify();

                    return res;

                } else if (AVal < 0) {
                    // A < 0. Must be asin.
                    // 1/sqrt(-A) * asin( (-2Ax - B) / sqrt(B^2 - 4AC) )
                    // Requires B^2 - 4AC > 0 for real solution range.

                    const negA = new Mul(new Num(-1), A).simplify();
                    const sqrtNegA = new Call('sqrt', [negA]).simplify();

                    const disc = new Sub(new Pow(B, new Num(2)), new Mul(new Num(4), new Mul(A, C))).simplify();
                    const sqrtDisc = new Call('sqrt', [disc]).simplify();

                    const argNum = new Sub(new Mul(new Mul(new Num(-2), A), varNode), B).simplify();
                    const arg = new Div(argNum, sqrtDisc).simplify();

                    return new Mul(new Div(num, sqrtNegA), new Call('asin', [arg])).simplify();
                }
            }
        }
        return null;
    }

    _integrateSubstitution(expr, varNode) {
        // Handle trig rewrites for substitution
        if (expr instanceof Call) {
            // tan/cot/tanh handled by _integrateTrig for simple cases now,
            // but keep rewrite for complex args e.g. tan(2x)
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

    _limit(expr, varNode, point, depth = 0, dir = 0) {
        // Special check for direction sensitivity on basic functions like sign/abs at discontinuity
        if (dir !== 0) {
             const val = expr.substitute(varNode, point).simplify();
             // If evaluation yields 0 for sign(x) but we have direction, check argument
             if (expr instanceof Call && expr.funcName === 'sign' && val.evaluateNumeric() === 0) {
                 const arg = expr.args[0];
                 // Check derivative of argument at point
                 const argDiff = arg.diff(varNode).substitute(varNode, point).evaluateNumeric();
                 if (!isNaN(argDiff) && argDiff !== 0) {
                     // If arg is increasing, sign(arg) approaches sign(dir)
                     // If arg is decreasing, sign(arg) approaches sign(-dir)
                     return new Num(Math.sign(argDiff * dir));
                 }
             }
        }
        if (depth > 5) return new Call("limit", [expr, varNode, point]);

        // Check for Indeterminate form NaN (from 0*Inf)
        // If expr is Mul, convert to Div
        if (expr instanceof Mul) {
            const val = expr.substitute(varNode, point).simplify();
            if (val instanceof Sym && val.name === 'NaN') {
                // 0 * Inf -> Transform to Div
                // f * g -> f / (1/g) or g / (1/f)
                // Prefer keeping log or simpler term in numerator
                // Heuristic: If one is exp or trig, keep in num?
                // Let's try 1st / (1/2nd)
                // If right is 1/x (Div(1,x)), 1/right is x.
                // 1/ (1/x) -> x.
                // if expr.right is Div(1, x).
                // den = new Div(1, Div(1, x)) -> simplifies to x.
                // so newExpr = num / den = x / x.

                const num = expr.left;
                const den = new Div(new Num(1), expr.right).simplify();
                const newExpr = new Div(num, den).simplify(); // Simplify here!
                return this._limit(newExpr, varNode, point, depth + 1, dir);
            }
        }

        if (expr instanceof Div) {
            let num = expr.left.substitute(varNode, point).simplify();
            let den = expr.right.substitute(varNode, point).simplify();

            // Handle Num or zero-value Num from simplification
            const isZero = (n) => (n instanceof Num && n.value === 0);
            const isInf = (n) => (n instanceof Sym && (n.name === 'Infinity' || n.name === 'infinity'));
            const isNegInf = (n) => (n instanceof Mul && n.left instanceof Num && n.left.value === -1 && isInf(n.right));
            const isInfinite = (n) => isInf(n) || isNegInf(n);

            if (isZero(num) && isZero(den)) {
                 const diffNum = expr.left.diff(varNode).simplify();
                 const diffDen = expr.right.diff(varNode).simplify();
                 // Use a temporary Division that doesn't eager simplify to float if possible,
                 // but Div constructor does not simplify unless operands are numbers.
                 // The issue is likely 'new Div(diffNum, diffDen)' if diffNum/diffDen are simple integers.
                 // Div.simplify() handles integer division. We want to avoid it if it returns float?
                 // But Div(1, 2).simplify() -> Div(1, 2) unless one is float.
                 // Check if simplify() was called on new Div inside the recursion?
                 // Yes, recursive _limit might simplify result.
                 return this._limit(new Div(diffNum, diffDen).simplify(), varNode, point, depth + 1, dir);
            }

            if (isInfinite(num) && isInfinite(den)) {
                 const diffNum = expr.left.diff(varNode).simplify();
                 const diffDen = expr.right.diff(varNode).simplify();
                 return this._limit(new Div(diffNum, diffDen).simplify(), varNode, point, depth + 1, dir);
            }

            if (num instanceof Num && den instanceof Num) {
                if (den.value !== 0) {
                    // Prefer exact division if integer
                    if (Number.isInteger(num.value) && Number.isInteger(den.value)) {
                        if (num.value % den.value === 0) return new Num(num.value / den.value);
                        // GCD reduction
                        const gcd = (a, b) => !b ? a : gcd(b, a % b);
                        const common = gcd(Math.abs(num.value), Math.abs(den.value));
                        if (common > 1) return new Div(new Num(num.value / common), new Num(den.value / common));
                        return new Div(num, den);
                    }
                    return new Num(num.value / den.value);
                }
                if (num.value !== 0 && den.value === 0) {
                     // Directional Limit Logic
                     if (dir !== 0) {
                         try {
                             const eps = (dir > 0 ? 1 : -1) * 1e-9;
                             const val = point.evaluateNumeric() + eps;
                             const denVal = expr.right.substitute(varNode, new Num(val)).evaluateNumeric();

                             if (!isNaN(denVal) && denVal !== 0) {
                                 const resultSign = (num.value > 0 ? 1 : -1) * (denVal > 0 ? 1 : -1);
                                 return (resultSign > 0) ? new Sym("Infinity") : new Mul(new Num(-1), new Sym("Infinity"));
                             }
                         } catch(e) {}
                     }
                     if (num.value < 0) return new Mul(new Num(-1), new Sym("Infinity"));
                     return new Sym("Infinity");
                }
            }
        }

        try {
            // Check for Infinity - Infinity form
            if (expr instanceof Sub || expr instanceof Add) {
                const L_left = this._limit(expr.left, varNode, point, depth + 1, dir);
                const L_right = this._limit(expr.right, varNode, point, depth + 1, dir);

                const isInf = (node) => node instanceof Sym && (node.name === 'Infinity' || node.name === 'infinity');
                const isNegInf = (node) => node instanceof Mul && node.left instanceof Num && node.left.value === -1 && isInf(node.right);

                const leftInf = isInf(L_left) || isNegInf(L_left);
                const rightInf = isInf(L_right) || isNegInf(L_right);

                if (leftInf && rightInf) {
                    // Possible Inf - Inf form
                    // Attempt to combine fractions or simplify
                    // E.g. 1/x - 1/sin(x) -> (sin(x)-x) / (x*sin(x))
                    // simplify() handles basic fraction combination if they are numeric
                    // But we need to force combination for symbolic fractions.
                    // Let's try combining manually if they are Divs.

                    // Check if expr can be combined
                    if (expr.left instanceof Div && expr.right instanceof Div) {
                        const n1 = expr.left.left;
                        const d1 = expr.left.right;
                        const n2 = expr.right.left;
                        const d2 = expr.right.right;

                        // (n1*d2 +/- n2*d1) / (d1*d2)
                        let num;
                        if (expr instanceof Add) num = new Add(new Mul(n1, d2), new Mul(n2, d1));
                        else num = new Sub(new Mul(n1, d2), new Mul(n2, d1));

                        const den = new Mul(d1, d2);
                        const combined = new Div(num, den).simplify(); // This might cancel terms or prepare for L'Hopital

                        if (combined instanceof Div) {
                             // Now evaluate limit of combined
                             return this._limit(combined, varNode, point, depth + 1);
                        }
                    }
                }
            }

            // Check for 1^Infinity form: lim (base^exp)
            if (expr instanceof Pow) {
                const L_base = this._limit(expr.left, varNode, point, depth + 1, dir);
                const L_exp = this._limit(expr.right, varNode, point, depth + 1, dir);

                // Check if base -> 1 and exp -> Infinity
                const isOne = (node) => node instanceof Num && Math.abs(node.value - 1) < 1e-9;
                const isInf = (node) => node instanceof Sym && (node.name === 'Infinity' || node.name === 'infinity');

                if (isOne(L_base) && isInf(L_exp)) {
                    // L = exp( lim (base - 1) * exp )
                    // For (1+1/x)^x -> base-1 = 1/x. exp = x. (1/x)*x = 1. lim=1. res=e.
                    const baseMinusOne = new Sub(expr.left, new Num(1)).simplify();
                    // Do NOT simplify Mul immediately if it collapses to constant 1 before limit?
                    // Actually, if it simplifies to 1, limit(1) is 1. That is correct.
                    // If it simplifies to 1, result is exp(1).
                    const newExpr = new Mul(baseMinusOne, expr.right).simplify();
                    const limNew = this._limit(newExpr, varNode, point, depth + 1, dir);
                    // Force return of exp(limNew)
                    // If limNew is 1, Call('exp', [1]).simplify() might return Num(2.718) or Num(1) if logic is wrong?
                    // Call.simplify for exp: if arg is 0 -> 1. If arg is 1 -> exp(1).
                    // Evaluate numeric might turn exp(1) -> 2.718...
                    return new Call('exp', [limNew]).simplify();
                }
            }

            const val = expr.substitute(varNode, point).simplify();

            // Check for directional limit at pole (Infinity or NaN)
            if (dir !== 0) {
                const isInf = (n) => n instanceof Sym && (n.name === 'Infinity' || n.name === 'infinity');
                const isUnsignedInf = (n) => isInf(n) || (n instanceof Sym && n.name === 'NaN');

                // If direct substitution gave Infinity (unsigned) or NaN (pole), we check direction
                if (isUnsignedInf(val)) {
                    // Test point approach
                    // Check sign of f(point + dir*delta)
                    // If point is numeric, use small delta.
                    const pVal = point.evaluateNumeric();
                    if (!isNaN(pVal)) {
                        const delta = dir * 1e-9;
                        const testVal = expr.substitute(varNode, new Num(pVal + delta)).evaluateNumeric();
                        if (!isNaN(testVal)) {
                            if (testVal > 1e6) return new Sym("Infinity");
                            if (testVal < -1e6) return new Mul(new Num(-1), new Sym("Infinity")).simplify();
                        }
                    }
                }
            }

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

        // Use Gaussian Elimination for Determinant O(n^3)
        // Check if matrix is numeric to safely use Gaussian without huge symbolic expression explosion
        // Even for symbolic, Gaussian is usually better than cofactor O(n!) for n > 4

        // Clone matrix
        const M = matrix.elements.map(row => [...row.elements]);
        let sign = 1;
        let det = new Num(1);

        for (let i = 0; i < rows; i++) {
            // Find pivot
            let pivotIdx = i;
            let found = false;

            // Search for non-zero pivot
            for(let k=i; k<rows; k++) {
                const val = M[k][i].simplify();
                // Check if zero
                let isZ = false;
                if (val instanceof Num && val.value === 0) isZ = true;

                if (!isZ) {
                    pivotIdx = k;
                    found = true;
                    break;
                }
            }

            if (!found) return new Num(0); // Zero column -> det 0

            if (pivotIdx !== i) {
                // Swap rows
                const temp = M[i];
                M[i] = M[pivotIdx];
                M[pivotIdx] = temp;
                sign *= -1;
            }

            const pivot = M[i][i];

            // Eliminate rows below
            for(let k=i+1; k<rows; k++) {
                const val = M[k][i];
                if (val instanceof Num && val.value === 0) continue;

                // factor = M[k][i] / pivot
                const factor = new Div(val, pivot).simplify();

                for(let j=i; j<cols; j++) {
                    // M[k][j] -= factor * M[i][j]
                    M[k][j] = new Sub(M[k][j], new Mul(factor, M[i][j])).simplify();
                }
            }
        }

        // Product of diagonal
        let prod = new Num(sign);
        for(let i=0; i<rows; i++) {
            prod = new Mul(prod, M[i][i]);
        }
        return prod.simplify();
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

    _union(l1, l2) {
        if (l1 instanceof Vec && l2 instanceof Vec) {
            const seen = new Set();
            const res = [];
            const add = (list) => {
                for(const e of list.elements) {
                    const s = e.toString();
                    if (!seen.has(s)) {
                        seen.add(s);
                        res.push(e);
                    }
                }
            };
            add(l1);
            add(l2);
            return this._sort(new Vec(res));
        }
        return new Call('union', [l1, l2]);
    }

    _intersect(l1, l2) {
        if (l1 instanceof Vec && l2 instanceof Vec) {
            const s1 = new Set();
            for(const e of l1.elements) s1.add(e.toString());

            const res = [];
            const seen = new Set(); // Prevent duplicates in result
            for(const e of l2.elements) {
                const s = e.toString();
                if (s1.has(s) && !seen.has(s)) {
                    seen.add(s);
                    res.push(e);
                }
            }
            return this._sort(new Vec(res));
        }
        return new Call('intersect', [l1, l2]);
    }

    _setdiff(l1, l2) {
        if (l1 instanceof Vec && l2 instanceof Vec) {
            const s2 = new Set();
            for(const e of l2.elements) s2.add(e.toString());

            const res = [];
            const seen = new Set();
            for(const e of l1.elements) {
                const s = e.toString();
                if (!s2.has(s) && !seen.has(s)) {
                    seen.add(s);
                    res.push(e);
                }
            }
            return this._sort(new Vec(res));
        }
        return new Call('setdiff', [l1, l2]);
    }

    _cross(v1, v2) {
        if (!(v1 instanceof Vec) || !(v2 instanceof Vec)) throw new Error("cross requires two vectors");

        let u = v1.elements;
        let v = v2.elements;

        if (u.length === 2 && v.length === 2) {
            u = [u[0], u[1], new Num(0)];
            v = [v[0], v[1], new Num(0)];
        }

        if (u.length !== 3 || v.length !== 3) throw new Error("cross requires 3D vectors (or 2D)");

        const a1 = u[0], a2 = u[1], a3 = u[2];
        const b1 = v[0], b2 = v[1], b3 = v[2];

        const c1 = new Sub(new Mul(a2, b3), new Mul(a3, b2));
        const c2 = new Sub(new Mul(a3, b1), new Mul(a1, b3));
        const c3 = new Sub(new Mul(a1, b2), new Mul(a2, b1));

        return new Vec([c1.simplify(), c2.simplify(), c3.simplify()]);
    }

    _primeFactors(n) {
        n = n.simplify();
        if (!(n instanceof Num)) throw new Error("primeFactors requires a number");
        let val = n.value;
        if (!Number.isInteger(val)) throw new Error("primeFactors requires an integer");
        if (val < 2) return new Vec([]);

        const factors = [];
        let d = 2;
        while (d * d <= val) {
            while (val % d === 0) {
                factors.push(new Num(d));
                val = Math.floor(val / d);
            }
            d++;
        }
        if (val > 1) factors.push(new Num(val));
        return new Vec(factors);
    }

    _angle(u, v) {
        // angle = acos( (u . v) / (|u|*|v|) )
        const dot = this._recursiveEval(new Call('dot', [u, v])).simplify();
        const normU = this._recursiveEval(new Call('norm', [u])).simplify();
        const normV = this._recursiveEval(new Call('norm', [v])).simplify();

        const cosTheta = new Div(dot, new Mul(normU, normV)).simplify();
        return new Call('acos', [cosTheta]).simplify();
    }

    _projection(u, v) {
        // proj_v(u) = (u . v) / (v . v) * v
        const dotUV = this._recursiveEval(new Call('dot', [u, v])).simplify();
        const dotVV = this._recursiveEval(new Call('dot', [v, v])).simplify();

        const factor = new Div(dotUV, dotVV).simplify();
        return new Mul(factor, v).simplify();
    }

    _toSpherical(v) {
        if (!(v instanceof Vec) || v.elements.length !== 3) throw new Error("toSpherical requires a 3D vector");
        const [x, y, z] = v.elements;

        // r = sqrt(x^2 + y^2 + z^2)
        const rSq = new Add(new Add(new Pow(x, new Num(2)), new Pow(y, new Num(2))), new Pow(z, new Num(2))).simplify();
        const r = new Call('sqrt', [rSq]).simplify();

        // theta = acos(z / r)  (polar angle from z-axis)
        const theta = new Call('acos', [new Div(z, r)]).simplify();

        // phi = arg(x + iy)
        const phi = this._arg(new Add(x, new Mul(new Sym('i'), y))).simplify();

        return new Vec([r, theta, phi]);
    }

    _toCylindrical(v) {
        if (!(v instanceof Vec) || v.elements.length !== 3) throw new Error("toCylindrical requires a 3D vector");
        const [x, y, z] = v.elements;

        // r = sqrt(x^2 + y^2)
        const rSq = new Add(new Pow(x, new Num(2)), new Pow(y, new Num(2))).simplify();
        const r = new Call('sqrt', [rSq]).simplify();

        // theta = arg(x + iy)
        const theta = this._arg(new Add(x, new Mul(new Sym('i'), y))).simplify();

        return new Vec([r, theta, z]);
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
             const ltR = pR.coeffs[degA] || new Num(0);
             const ltB = pB.coeffs[degB] || new Num(0);

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

         // Polynomial LCM: (A * B) / GCD(A, B)
         const gcdVal = this._gcd(a, b);
         if (!(gcdVal instanceof Call && gcdVal.funcName === 'gcd')) {
             const prod = new Mul(a, b).simplify();
             return new Div(prod, gcdVal).simplify();
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

    _studentTPDF(x, df) {
        x = x.simplify();
        df = df.simplify();

        if (x instanceof Num && df instanceof Num) {
            const xv = x.value;
            const v = df.value;
            if (v <= 0) return new Num(0);

            // Use Log Gamma to avoid overflow for large df
            // log(coeff) = logGamma((v+1)/2) - 0.5*log(v*pi) - logGamma(v/2)
            const logNum = this._logGamma((v + 1) / 2);
            const logDen = 0.5 * Math.log(v * Math.PI) + this._logGamma(v / 2);
            const logCoeff = logNum - logDen;

            const base = 1 + (xv * xv) / v;
            const power = -(v + 1) / 2;

            // result = exp(logCoeff) * base^power
            // result = exp(logCoeff + power * log(base))
            const logResult = logCoeff + power * Math.log(base);
            return new Num(Math.exp(logResult));
        }
        return new Call('studentTPDF', [x, df]);
    }

    _studentTCDF(x, df) {
        x = x.simplify();
        df = df.simplify();

        if (x instanceof Num && df instanceof Num) {
            const t = x.value;
            const v = df.value;
            if (v <= 0) return new Num(0);

            // F(t) = 1 - 0.5 * I_x(v/2, 1/2) where x = v / (v + t^2)
            const xt = v / (v + t * t);
            const p = 0.5 * this._betainc(xt, v / 2, 0.5);

            if (t > 0) return new Num(1 - p);
            return new Num(p);
        }
        return new Call('studentTCDF', [x, df]);
    }

    _fPDF(x, d1, d2) {
        x = x.simplify();
        d1 = d1.simplify();
        d2 = d2.simplify();

        if (x instanceof Num && d1 instanceof Num && d2 instanceof Num) {
            const xv = x.value;
            const v1 = d1.value;
            const v2 = d2.value;

            if (xv < 0) return new Num(0);
            if (v1 <= 0 || v2 <= 0) return new Num(0);

            // f(x) = sqrt( (v1*x)^v1 * v2^v2 / (v1*x + v2)^(v1+v2) ) / (x * B(v1/2, v2/2))
            // Logarithmic calculation for stability
            const logNum = (v1 / 2) * Math.log(v1 * xv) + (v2 / 2) * Math.log(v2);
            const logDen = ( (v1 + v2) / 2 ) * Math.log(v1 * xv + v2);
            const logBeta = this._logGamma(v1/2) + this._logGamma(v2/2) - this._logGamma((v1+v2)/2);

            // log(f) = logNum - logDen - log(x) - logBeta
            const logF = logNum - logDen - Math.log(xv) - logBeta;
            return new Num(Math.exp(logF));
        }
        return new Call('fPDF', [x, d1, d2]);
    }

    _fCDF(x, d1, d2) {
        x = x.simplify();
        d1 = d1.simplify();
        d2 = d2.simplify();

        if (x instanceof Num && d1 instanceof Num && d2 instanceof Num) {
            const xv = x.value;
            const v1 = d1.value;
            const v2 = d2.value;

            if (xv <= 0) return new Num(0);
            if (v1 <= 0 || v2 <= 0) return new Num(0);

            // F(x) = I_z(d1/2, d2/2) where z = d1*x / (d1*x + d2)
            const z = (v1 * xv) / (v1 * xv + v2);
            return new Num(this._betainc(z, v1/2, v2/2));
        }
        return new Call('fCDF', [x, d1, d2]);
    }

    _invT(area, df) {
        area = area.simplify();
        df = df.simplify();

        if (area instanceof Num && df instanceof Num) {
            const p = area.value;
            const v = df.value;
            if (p <= 0 || p >= 1) return new Sym("NaN");
            if (v <= 0) return new Sym("NaN");

            // Use binary search on CDF
            // Range roughly -10 to 10? T can be large.
            // Start with rough guess using Normal approx or just search range.
            let min = -100, max = 100;
            // Expand range if needed
            while (this._studentTCDF(new Num(min), df).value > p) min *= 2;
            while (this._studentTCDF(new Num(max), df).value < p) max *= 2;

            for (let i = 0; i < 100; i++) {
                const mid = (min + max) / 2;
                const cdf = this._studentTCDF(new Num(mid), df).value;
                if (Math.abs(cdf - p) < 1e-9) return new Num(mid);
                if (cdf < p) min = mid;
                else max = mid;
            }
            return new Num((min + max) / 2);
        }
        return new Call('invT', [area, df]);
    }

    _chisquareCDF(x, k) {
        x = x.simplify();
        k = k.simplify();

        if (x instanceof Num && k instanceof Num) {
            const xv = x.value;
            const kv = k.value;
            if (kv <= 0) return new Num(0);
            if (xv <= 0) return new Num(0);

            // P(k/2, x/2)
            // _gammainc returns P(s, x) (Regularized Incomplete Gamma)
            // Arguments for _gammainc are s, x.
            // We need P(k/2, x/2).
            return new Num(this._gammainc(kv / 2, xv / 2));
        }
        return new Call('chisquareCDF', [x, k]);
    }

    _invChiSquare(area, k) {
        area = area.simplify();
        k = k.simplify();

        if (area instanceof Num && k instanceof Num) {
            const p = area.value;
            const kv = k.value;
            if (p <= 0 || p >= 1) return new Sym("NaN");
            if (kv <= 0) return new Sym("NaN");

            // Binary search on CDF (x > 0)
            let min = 0, max = 100;
            while (this._chisquareCDF(new Num(max), k).value < p) max *= 2;

            for (let i = 0; i < 100; i++) {
                const mid = (min + max) / 2;
                const cdf = this._chisquareCDF(new Num(mid), k).value;
                if (Math.abs(cdf - p) < 1e-9) return new Num(mid);
                if (cdf < p) min = mid;
                else max = mid;
            }
            return new Num((min + max) / 2);
        }
        return new Call('invChiSquare', [area, k]);
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

    _conservative(field, vars) {
        if (!(field instanceof Vec) || !(vars instanceof Vec)) throw new Error("Arguments must be vectors");
        // Check dimension
        const n = field.elements.length;
        if (vars.elements.length !== n) throw new Error("Dimension mismatch");

        // Check partials: dFi/dxj = dFj/dxi
        for(let i=0; i<n; i++) {
            for(let j=i+1; j<n; j++) {
                const dFi_dxj = field.elements[i].diff(vars.elements[j]).simplify();
                const dFj_dxi = field.elements[j].diff(vars.elements[i]).simplify();
                // Check equality
                const diff = new Sub(dFi_dxj, dFj_dxi).simplify();
                if (diff instanceof Num && diff.value === 0) continue;
                // If symbolic difference is 0?
                // simplify() handles basic algebra.
                // If not 0, return 0 (false)
                return new Num(0);
            }
        }
        return new Num(1);
    }

    _jacobian(vector, vars) {
        if (!(vector instanceof Vec)) throw new Error("First argument to jacobian must be a vector of expressions");
        if (!(vars instanceof Vec)) throw new Error("Second argument to jacobian must be a vector of variables");

        const rows = [];
        for (const expr of vector.elements) {
            const row = [];
            for (const v of vars.elements) {
                row.push(expr.diff(v).simplify());
            }
            rows.push(new Vec(row));
        }
        return new Vec(rows);
    }

    _hessian(expr, vars) {
        if (!(vars instanceof Vec)) throw new Error("Second argument to hessian must be a vector of variables");

        const rows = [];
        for (const v1 of vars.elements) {
            const row = [];
            for (const v2 of vars.elements) {
                // d^2 f / (dx_i dx_j)
                const d1 = expr.diff(v1).simplify();
                const d2 = d1.diff(v2).simplify();
                row.push(d2);
            }
            rows.push(new Vec(row));
        }
        return new Vec(rows);
    }

    _potential(field, vars) {
        const isCons = this._conservative(field, vars);
        if (isCons.value === 0) throw new Error("Field is not conservative");

        const n = field.elements.length;
        let phi = new Num(0);

        // phi = int(F1, x1) + C(x2, x3...)
        // Then dphi/dx2 = F2 ...

        // We can do this iteratively.
        // Current phi accounts for first k variables.
        // We calculate partial derivative of phi w.r.t next variable, subtract from F_next.
        // Remainder should depend only on x_{k+1}...
        // Integrate remainder.

        for(let i=0; i<n; i++) {
             // Differentiate current phi w.r.t current var
             const dPhi = phi.diff(vars.elements[i]).simplify();
             // Remainder to integrate
             const rem = new Sub(field.elements[i], dPhi).simplify();
             // Integrate rem w.r.t current var
             const term = rem.integrate(vars.elements[i]).simplify();
             phi = new Add(phi, term).simplify();
        }
        return phi;
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

    _modInverse(a, m) {
        if (a instanceof Num && m instanceof Num) {
             let n = a.value;
             let mod = m.value;

             // Extended Euclidean Algorithm
             let t = 0, newt = 1;
             let r = mod, newr = n;

             while (newr !== 0) {
                 const quotient = Math.floor(r / newr);
                 [t, newt] = [newt, t - quotient * newt];
                 [r, newr] = [newr, r - quotient * newr];
             }

             if (r > 1) throw new Error("a is not invertible");
             if (t < 0) t = t + mod;

             return new Num(t);
        }
        return new Call('modInverse', [a, m]);
    }

    _modPow(base, exp, mod) {
        if (base instanceof Num && exp instanceof Num && mod instanceof Num) {
             let res = 1n;
             let b = BigInt(Math.floor(base.value));
             let e = BigInt(Math.floor(exp.value));
             let m = BigInt(Math.floor(mod.value));

             b = b % m;
             while (e > 0n) {
                 if (e % 2n === 1n) res = (res * b) % m;
                 e = e / 2n;
                 b = (b * b) % m;
             }
             return new Num(Number(res));
        }
        return new Call('modPow', [base, exp, mod]);
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

        // Check for 0 (undefined usually, or 0)
        if (z instanceof Num && z.value === 0) return new Num(0);

        if (z instanceof Num) {
             if (z.value >= 0) return new Num(0);
             return new Sym('pi');
        }

        // Handle complex parts
        const parts = this._getComplexParts(z);
        const re = parts.re;
        const im = parts.im;

        const reVal = re.evaluateNumeric();
        const imVal = im.evaluateNumeric();

        if (!isNaN(reVal) && !isNaN(imVal)) {
            const val = Math.atan2(imVal, reVal);
            return new Num(val);
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

    _rsolve(eq, term, conds) {
        // Linear Constant Coefficient Recurrence Relation
        // eq: a(n) - a(n-1) - a(n-2) = 0
        // term: a(n) (Call)
        // conds: [a(0)=0, a(1)=1] (Vec of Eqs)

        if (!(term instanceof Call)) throw new Error("Recurrence term must be a function call, e.g., a(n)");
        const funcName = term.funcName;
        const nVar = term.args[0]; // n

        // Normalize equation to Expression = 0
        let expr = eq;
        if (eq instanceof Eq) expr = new Sub(eq.left, eq.right).simplify();
        expr = expr.expand().simplify();

        // 1. Identify Shifted Terms
        // Find all occurrences of a(n+k) or a(n-k)
        // Map shift k -> coefficient
        const coeffs = {}; // integer shift -> Expr
        let minShift = 0;
        let maxShift = 0;
        let foundAny = false;

        const terms = [];
        const collect = (e, sign=1) => {
            if (e instanceof Add) { collect(e.left, sign); collect(e.right, sign); }
            else if (e instanceof Sub) { collect(e.left, sign); collect(e.right, -sign); }
            else terms.push({e, sign});
        };
        collect(expr);

        for(const item of terms) {
            // Check for coeff * a(...)
            // Separate coeff and call
            let coeff = new Num(item.sign);
            let callNode = null;

            const analyze = (node) => {
                if (node instanceof Call && node.funcName === funcName) return node;
                if (node instanceof Mul) {
                    const l = analyze(node.left);
                    if (l) { coeff = new Mul(coeff, node.right); return l; }
                    const r = analyze(node.right);
                    if (r) { coeff = new Mul(coeff, node.left); return r; }
                }
                return null;
            };

            callNode = analyze(item.e);

            if (callNode) {
                // Check argument n+k or n-k
                const arg = callNode.args[0].simplify();
                let shift = 0;
                // n -> 0
                // n+k -> k
                // n-k -> -k
                if (arg.toString() === nVar.toString()) shift = 0;
                else if (arg instanceof Add && arg.left.toString() === nVar.toString() && arg.right instanceof Num) {
                    shift = arg.right.value;
                }
                else if (arg instanceof Add && arg.right.toString() === nVar.toString() && arg.left instanceof Num) {
                    shift = arg.left.value;
                }
                else if (arg instanceof Sub && arg.left.toString() === nVar.toString() && arg.right instanceof Num) {
                    shift = -arg.right.value;
                }
                else {
                    // Complex shift or not linear shift
                    // throw new Error("Only constant shifts supported: " + arg);
                    return new Call('rsolve', [eq, term, conds]);
                }

                if (!foundAny) { minShift = shift; maxShift = shift; foundAny = true; }
                else {
                    if (shift < minShift) minShift = shift;
                    if (shift > maxShift) maxShift = shift;
                }

                coeffs[shift] = new Add(coeffs[shift] || new Num(0), coeff).simplify();
            } else {
                // Non-homogeneous term?
                // For now, assume homogeneous (=0) or ignore/error
                // Or subtract from 0?
                // We strictly solve homogeneous for now.
                const val = item.e.evaluateNumeric();
                if (isNaN(val) || Math.abs(val) > 1e-9) {
                    // return new Call('rsolve', [eq, term, conds]); // Fallback
                }
            }
        }

        if (!foundAny) return new Call('rsolve', [eq, term, conds]);

        // 2. Characteristic Polynomial
        // Shift indices so lowest is 0.
        // P(r) = sum c_k * r^(k - minShift)
        let charPoly = new Num(0);
        const rVar = new Sym('r');

        for(const s in coeffs) {
            const shift = parseInt(s);
            const power = shift - minShift;
            const c = coeffs[s];
            if (power === 0) charPoly = new Add(charPoly, c);
            else charPoly = new Add(charPoly, new Mul(c, new Pow(rVar, new Num(power))));
        }

        charPoly = charPoly.simplify();

        // 3. Solve for Roots
        const rootsRes = this._solve(charPoly, rVar);
        let roots = [];
        if (rootsRes instanceof Call && rootsRes.funcName === 'set') roots = rootsRes.args;
        else if (rootsRes instanceof Expr && !(rootsRes instanceof Call && rootsRes.funcName === 'solve')) roots = [rootsRes];
        else {
            // Unsolvable
            return new Call('rsolve', [eq, term, conds]);
        }

        // Count multiplicities
        const rootCounts = {};
        for(const r of roots) {
            const s = r.toString();
            rootCounts[s] = (rootCounts[s] || 0) + 1;
        }

        // 4. General Solution
        // sum (C_i * r^n)
        // repeated root r (m times): C1*r^n + C2*n*r^n + ... + Cm*n^(m-1)*r^n
        let solution = new Num(0);
        let constants = [];
        let cIdx = 1;

        for(const rStr in rootCounts) {
            // Find the actual root object
            const root = roots.find(r => r.toString() === rStr);
            const m = rootCounts[rStr];

            for(let i=0; i<m; i++) {
                const C = new Sym('C' + cIdx++);
                constants.push(C);
                let termSol;
                // n^i * r^n
                const rPowN = new Pow(root, nVar);
                if (i === 0) termSol = rPowN;
                else if (i === 1) termSol = new Mul(nVar, rPowN);
                else termSol = new Mul(new Pow(nVar, new Num(i)), rPowN);

                solution = new Add(solution, new Mul(C, termSol));
            }
        }
        solution = solution.simplify();

        // 5. Apply Initial Conditions
        if (conds && conds instanceof Vec) {
            const system = [];
            const vars = new Vec(constants);

            for(const cond of conds.elements) {
                // a(k) = val
                // Extract k and val
                let kVal = null;
                let val = null;
                if (cond instanceof Eq) {
                    if (cond.left instanceof Call && cond.left.funcName === funcName) {
                        kVal = cond.left.args[0].evaluateNumeric();
                        val = cond.right;
                    }
                }
                if (kVal !== null && val !== null) {
                    // Substitute n=k in solution -> eq
                    const eqSol = solution.substitute(nVar, new Num(kVal)).simplify();
                    system.push(new Eq(eqSol, val));
                }
            }

            if (system.length > 0) {
                const solConsts = this._solveSystem(new Vec(system), vars);
                if (solConsts instanceof Vec && solConsts.elements.length === constants.length) {
                    // Substitute constants back
                    let finalSol = solution;
                    for(let i=0; i<constants.length; i++) {
                        finalSol = finalSol.substitute(constants[i], solConsts.elements[i]);
                    }
                    return finalSol.simplify();
                }
            }
        }

        return solution;
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
                // Ensure order is 1 (node.args.length == 2 or args[2] == 1)
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
                    if (inner instanceof Call && inner.funcName === 'diff' && inner.node.args.length === 2) {
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
        // Improved logic: If solving leads to complex roots for a real quadratic, we might want to avoid splitting
        // IF we are aiming for real integration (atan).
        // But partfrac is strictly partial fraction decomposition.
        // In Complex field, 1/(x^2+1) -> 1/2i(x-i) - 1/2i(x+i).
        // If we want to support real integration, we should handle irreducible quadratics.
        // However, this function _partfrac is supposed to return the decomposition.
        // We will stick to full decomposition if possible, but handle integration separately.

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

        // Ensure maxDeg exists in coeffs
        if (!coeffs[trueMaxDeg]) coeffs[trueMaxDeg] = new Num(0);

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

    _ztrans(expr, n, z) {
        expr = expr.expand().simplify();

        if (expr instanceof Add) {
            return new Add(this._ztrans(expr.left, n, z), this._ztrans(expr.right, n, z)).simplify();
        }
        if (expr instanceof Sub) {
            return new Sub(this._ztrans(expr.left, n, z), this._ztrans(expr.right, n, z)).simplify();
        }
        if (expr instanceof Mul) {
            if (!this._dependsOn(expr.left, n)) return new Mul(expr.left, this._ztrans(expr.right, n, z)).simplify();
            if (!this._dependsOn(expr.right, n)) return new Mul(expr.right, this._ztrans(expr.left, n, z)).simplify();
        }

        // Table
        // 1 -> z / (z-1)
        if (expr instanceof Num && expr.value === 1) {
            return new Div(z, new Sub(z, new Num(1))).simplify();
        }
        if (!this._dependsOn(expr, n)) {
            // c -> c * z / (z-1)
            const Z = new Div(z, new Sub(z, new Num(1)));
            return new Mul(expr, Z).simplify();
        }

        // n -> z / (z-1)^2
        if (expr instanceof Sym && expr.name === n.name) {
            return new Div(z, new Pow(new Sub(z, new Num(1)), new Num(2))).simplify();
        }

        // a^n -> z / (z-a)
        if (expr instanceof Pow && expr.left instanceof Sym && expr.right instanceof Sym && expr.right.name === n.name) {
            const a = expr.left; // Assume a does not depend on n
            return new Div(z, new Sub(z, a)).simplify();
        }
        if (expr instanceof Pow && expr.right instanceof Sym && expr.right.name === n.name) {
             const a = expr.left;
             if (!this._dependsOn(a, n)) {
                 return new Div(z, new Sub(z, a)).simplify();
             }
        }

        // sin(an) -> z*sin(a) / (z^2 - 2z*cos(a) + 1)
        // cos(an) -> z(z-cos(a)) / (z^2 - 2z*cos(a) + 1)
        if (expr instanceof Call && (expr.funcName === 'sin' || expr.funcName === 'cos')) {
            const arg = expr.args[0];
            const poly = this._getPolyCoeffs(arg, n);
            let a = new Num(1);
            if (poly && poly.maxDeg === 1) {
                a = poly.coeffs[1]; // a*n
            } else {
                return new Call('ztrans', [expr, n, z]);
            }

            const den = new Add(new Sub(new Pow(z, new Num(2)), new Mul(new Num(2), new Mul(z, new Call('cos', [a])))), new Num(1));

            if (expr.funcName === 'sin') {
                const num = new Mul(z, new Call('sin', [a]));
                return new Div(num, den).simplify();
            }
            if (expr.funcName === 'cos') {
                const num = new Mul(z, new Sub(z, new Call('cos', [a])));
                return new Div(num, den).simplify();
            }
        }

        return new Call('ztrans', [expr, n, z]);
    }

    _iztrans(expr, z, n) {
        // Inverse Z-Transform using Partial Fractions
        // F(z) -> PFD F(z)/z -> sum A/(z-p)
        // inv(A * z / (z-p)) -> A * p^n

        // Try to factor out z? Or decompose F(z)/z.
        // My _partfrac handles expr w.r.t var.
        const F_div_z = new Div(expr, z).simplify();
        const pf = this._partfrac(F_div_z, z);

        if (pf instanceof Call && pf.funcName === 'partfrac' && pf.toString() === F_div_z.toString()) {
             // Failed decomposition
             return new Call('iztrans', [expr, z, n]);
        }

        // pf is sum of terms like A/(z-p)
        // We want to transform back to F(z) form: A*z / (z-p)
        // So multiply pf by z
        const F_expanded = new Mul(z, pf).expand().simplify();

        // Linearity
        const terms = [];
        const collect = (node) => {
            if (node instanceof Add) { collect(node.left); collect(node.right); }
            else if (node instanceof Sub) { collect(node.left); collect(new Mul(new Num(-1), node.right)); }
            else terms.push(node);
        };
        collect(F_expanded);

        let result = new Num(0);
        for(const term of terms) {
            // Check form A * z / (z - p)
            // term is simplified.
            // e.g. z / (z-1) -> 1
            // e.g. z / (z-a) -> a^n
            // e.g. k*z / (z-a) -> k*a^n

            // Extract coefficient independent of z?
            // This is hard to pattern match generally on simplified expression.
            // But usually we get k * z * (z-p)^-1
            // or k * z / (z-p)^m

            // Let's rely on basic patterns:
            // 1. k (constant) -> k * delta(n) (impulse) -- but standard z-trans usually 1 -> delta.
            // Wait, ztrans(1) = z/(z-1).
            // ztrans(delta(n)) = 1.
            // So if term is constant k, inv is k * dirac(n)?
            // Usually we deal with causal sequences.

            // Pattern: k * z / (z-a)^m
            // m=1: k * a^n
            // m=2: n * a^(n-1) ? z/(z-1)^2 -> n.  az/(z-a)^2 -> n a^n.
            // z/(z-a)^2 -> n * a^(n-1)

            // Let's handle simple case z/(z-a)
            // Check if term has z in numerator and (z-a) in denominator
            if (!this._dependsOn(term, z)) {
                // Constant C -> C * delta(n)
                // result += C * dirac(n)
                result = new Add(result, new Mul(term, new Call('dirac', [n])));
                continue;
            }

            // Simple pattern matching for Div
            if (term instanceof Div || (term instanceof Mul && (term.left instanceof Div || term.right instanceof Div))) {
                 // Try to match k * z / (z-p)
                 // divide term by z, see if it is k/(z-p)
                 const ratio = new Div(term, z).simplify(); // Should be k/(z-p)
                 // Check if ratio denominator is linear z-p
                 if (ratio instanceof Div) {
                     const num = ratio.left;
                     const den = ratio.right;
                     const polyDen = this._getPolyCoeffs(den, z);
                     if (polyDen && polyDen.maxDeg === 1 && !this._dependsOn(num, z)) {
                         // den = c1*z + c0 = c1(z + c0/c1)
                         // form: num / (c1(z - p)) -> (num/c1) / (z - p)
                         // p = -c0/c1
                         const c1 = polyDen.coeffs[1];
                         const c0 = polyDen.coeffs[0] || new Num(0);
                         const p = new Div(new Mul(new Num(-1), c0), c1).simplify();
                         const k = new Div(num, c1).simplify();

                         // inv( k * z / (z-p) ) = k * p^n
                         // if p=1, k
                         let inv;
                         if (p instanceof Num && p.value === 1) inv = k;
                         else inv = new Mul(k, new Pow(p, n));

                         result = new Add(result, inv);
                         continue;
                     }
                     // Check (z-p)^2
                     // z/(z-1)^2 -> n
                     if (polyDen && polyDen.maxDeg === 2) {
                         // Check square (z-p)^2 = z^2 - 2pz + p^2
                         // We can also use partial fraction on ratio (which failed earlier if we are here?)
                         // ratio is k/(z-p)^2 ?
                         // This requires identifying square.
                     }
                 }
            }

            // Fallback
            result = new Add(result, new Call('iztrans', [term, z, n]));
        }

        return result.simplify();
    }

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

    _tTest2(data1, data2) {
        // Two-sample t-test (independent, equal variance assumed for simplicity or Welch's)
        // Using Welch's t-test by default as it's safer
        if (!(data1 instanceof Vec) || !(data2 instanceof Vec)) throw new Error("Data must be lists");

        const n1 = new Num(data1.elements.length);
        const n2 = new Num(data2.elements.length);
        const m1 = this._mean(data1);
        const m2 = this._mean(data2);
        const v1 = this._variance(data1);
        const v2 = this._variance(data2);

        const num = new Sub(m1, m2).simplify();
        const denSq = new Add(new Div(v1, n1), new Div(v2, n2)).simplify();
        const t = new Div(num, new Call('sqrt', [denSq])).simplify();

        // Welch-Satterthwaite equation for df
        // df = (v1/n1 + v2/n2)^2 / [ (v1/n1)^2/(n1-1) + (v2/n2)^2/(n2-1) ]
        const term1 = new Div(v1, n1);
        const term2 = new Div(v2, n2);
        const dfNum = new Pow(new Add(term1, term2), new Num(2));
        const dfDen1 = new Div(new Pow(term1, new Num(2)), new Sub(n1, new Num(1)));
        const dfDen2 = new Div(new Pow(term2, new Num(2)), new Sub(n2, new Num(1)));
        const df = new Div(dfNum, new Add(dfDen1, dfDen2)).simplify();

        return new Vec([t, df]);
    }

    _zInterval(data, sigma, level) {
         if (!(data instanceof Vec)) throw new Error("Data must be a list");
         const n = new Num(data.elements.length);
         const mean = this._mean(data);
         const alpha = new Div(new Sub(new Num(1), level), new Num(2)).simplify();
         // critical z: area = 1 - alpha
         const area = new Sub(new Num(1), alpha).simplify();
         const zCrit = this._invNorm(area, new Num(0), new Num(1));

         const ME = new Mul(zCrit, new Div(sigma, new Call('sqrt', [n]))).simplify();
         return new Vec([new Sub(mean, ME).simplify(), new Add(mean, ME).simplify()]);
    }

    _tInterval(data, level) {
         if (!(data instanceof Vec)) throw new Error("Data must be a list");
         const n = new Num(data.elements.length);
         const mean = this._mean(data);
         const s = this._std(data);
         const df = new Sub(n, new Num(1));

         const alpha = new Div(new Sub(new Num(1), level), new Num(2)).simplify();
         const area = new Sub(new Num(1), alpha).simplify();
         const tCrit = this._invT(area, df);

         const ME = new Mul(tCrit, new Div(s, new Call('sqrt', [n]))).simplify();
         return new Vec([new Sub(mean, ME).simplify(), new Add(mean, ME).simplify()]);
    }

    _propTest(x, n, p0) {
        // One proportion z-test
        // z = (p_hat - p0) / sqrt(p0(1-p0)/n)
        const pHat = new Div(x, n).simplify();
        const num = new Sub(pHat, p0).simplify();
        const denSq = new Div(new Mul(p0, new Sub(new Num(1), p0)), n).simplify();
        const z = new Div(num, new Call('sqrt', [denSq])).simplify();

        // Calculate p-value (two-tailed)
        // p = 2 * (1 - normalCDF(abs(z)))
        const absZ = new Call('abs', [z]).simplify();
        const prob = this._normalCDF(absZ, new Num(0), new Num(1));
        const pVal = new Mul(new Num(2), new Sub(new Num(1), prob)).simplify();

        return new Vec([z, pVal]);
    }

    _propTest2(x1, n1, x2, n2) {
         // Two proportion z-test
         // p_pool = (x1+x2)/(n1+n2)
         // z = (p1 - p2) / sqrt( p_pool(1-p_pool)(1/n1 + 1/n2) )
         const p1 = new Div(x1, n1).simplify();
         const p2 = new Div(x2, n2).simplify();
         const pPool = new Div(new Add(x1, x2), new Add(n1, n2)).simplify();

         const num = new Sub(p1, p2).simplify();
         const qPool = new Sub(new Num(1), pPool).simplify();
         const factor = new Add(new Div(new Num(1), n1), new Div(new Num(1), n2)).simplify();
         const denSq = new Mul(new Mul(pPool, qPool), factor).simplify();
         const z = new Div(num, new Call('sqrt', [denSq])).simplify();

         const absZ = new Call('abs', [z]).simplify();
         const prob = this._normalCDF(absZ, new Num(0), new Num(1));
         const pVal = new Mul(new Num(2), new Sub(new Num(1), prob)).simplify();

         return new Vec([z, pVal]);
    }

    _chiSquareTest(observed, expected) {
        // Goodness of Fit
        // chi2 = sum( (O-E)^2 / E )
        if (!(observed instanceof Vec) || !(expected instanceof Vec)) throw new Error("Arguments must be lists");
        if (observed.elements.length !== expected.elements.length) throw new Error("Lists must be equal length");

        let sum = new Num(0);
        for(let i=0; i<observed.elements.length; i++) {
            const O = observed.elements[i];
            const E = expected.elements[i];
            const diff = new Sub(O, E);
            const term = new Div(new Pow(diff, new Num(2)), E);
            sum = new Add(sum, term).simplify();
        }
        return sum;
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

    _powerRegression(data) {
        if (!(data instanceof Vec)) throw new Error("Data must be list");
        // Model: y = A * x^B
        // Linearize: ln(y) = ln(A) + B * ln(x)
        // Y = a + b * X
        // Y = ln(y), X = ln(x), a = ln(A), b = B

        const logData = [];
        for(const pt of data.elements) {
             let x, y;
             if (pt instanceof Vec && pt.elements.length >= 2) {
                 x = pt.elements[0].evaluateNumeric();
                 y = pt.elements[1].evaluateNumeric();
             }

             if (isNaN(x) || isNaN(y)) throw new Error("Regression data must be numeric");
             if (x <= 0 || y <= 0) throw new Error("Power regression requires positive x and y values");

             const logX = Math.log(x);
             const logY = Math.log(y);
             logData.push({x: logX, y: logY});
        }

        let sumX=0, sumY=0, sumXY=0, sumXX=0, n=0;
        for(const pt of logData) {
             sumX += pt.x; sumY += pt.y; sumXY += pt.x*pt.y; sumXX += pt.x*pt.x;
             n++;
        }

        const den = n*sumXX - sumX*sumX;
        if(den===0) throw new Error("Vertical line");

        const B = (n*sumXY - sumX*sumY) / den; // Slope b
        const a = (sumY - B*sumX) / n; // Intercept a = ln(A)
        const A = Math.exp(a);

        // A * x^B
        return new Mul(new Num(A), new Pow(new Sym('x'), new Num(B))).simplify();
    }

    _logRegression(data) {
        if (!(data instanceof Vec)) throw new Error("Data must be list");
        // Model: y = A + B * ln(x)
        // Linearize: y = A + B * X
        // X = ln(x)

        const logData = [];
        for(const pt of data.elements) {
             let x, y;
             if (pt instanceof Vec && pt.elements.length >= 2) {
                 x = pt.elements[0].evaluateNumeric();
                 y = pt.elements[1].evaluateNumeric();
             }

             if (isNaN(x) || isNaN(y)) throw new Error("Regression data must be numeric");
             if (x <= 0) throw new Error("Logarithmic regression requires positive x values");

             const logX = Math.log(x);
             logData.push({x: logX, y: y});
        }

        let sumX=0, sumY=0, sumXY=0, sumXX=0, n=0;
        for(const pt of logData) {
             sumX += pt.x; sumY += pt.y; sumXY += pt.x*pt.y; sumXX += pt.x*pt.x;
             n++;
        }

        const den = n*sumXX - sumX*sumX;
        if(den===0) throw new Error("Vertical line");

        const B = (n*sumXY - sumX*sumY) / den; // Slope B
        const A = (sumY - B*sumX) / n; // Intercept A

        // A + B * ln(x)
        return new Add(new Num(A), new Mul(new Num(B), new Call('ln', [new Sym('x')]))).simplify();
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

    _moment(list, k) {
        if (!(list instanceof Vec)) throw new Error("First argument to moment must be a list");
        if (!(k instanceof Num)) throw new Error("Second argument to moment must be a number");

        const mean = this._mean(list);
        let sum = new Num(0);
        const n = list.elements.length;
        for(const x of list.elements) {
            const diff = new Sub(x, mean);
            sum = new Add(sum, new Pow(diff, k));
        }
        return new Div(sum.simplify(), new Num(n)).simplify();
    }

    _skewness(list) {
        // mu3 / sigma^3
        const m3 = this._moment(list, new Num(3));
        const m2 = this._moment(list, new Num(2)); // Variance (population)
        // skewness is usually sample skewness or population?
        // Moment gives population moment.
        // g1 = m3 / m2^(3/2)
        const sigma = new Call('sqrt', [m2]);
        const sigma3 = new Pow(sigma, new Num(3));
        return new Div(m3, sigma3).simplify();
    }

    _kurtosis(list) {
        // mu4 / sigma^4
        const m4 = this._moment(list, new Num(4));
        const m2 = this._moment(list, new Num(2));
        const sigma4 = new Pow(m2, new Num(2)); // (sigma^2)^2
        return new Div(m4, sigma4).simplify();
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

    _bernoulli(n) {
        n = n.simplify();
        if (n instanceof Num && Number.isInteger(n.value) && n.value >= 0) {
            const k = n.value;
            if (k === 0) return new Num(1);
            if (k === 1) return new Div(new Num(-1), new Num(2)); // B1 = -1/2
            if (k % 2 !== 0) return new Num(0); // Odd B_k are 0 for k > 1

            // Compute recursively: B_m = -1/(m+1) * sum(binom(m+1, k) * B_k) for k=0..m-1
            // Cache?
            const cache = { 0: new Num(1) };
            for(let m=1; m<=k; m++) {
                if (m === 1) { cache[1] = new Div(new Num(-1), new Num(2)); continue; }
                if (m % 2 !== 0) { cache[m] = new Num(0); continue; }

                let sum = new Num(0);
                for(let j=0; j<m; j++) {
                    const bin = this._nCr(new Num(m+1), new Num(j));
                    sum = new Add(sum, new Mul(bin, cache[j])).simplify();
                }
                const term = new Div(new Num(-1), new Num(m+1));
                cache[m] = new Mul(term, sum).simplify();
            }
            return cache[k];
        }
        return new Call('bernoulli', [n]);
    }

    _harmonic(n) {
        // H_n = sum(1/k, k, 1, n)
        n = n.simplify();
        if (n instanceof Num && Number.isInteger(n.value) && n.value >= 0) {
            if (n.value === 0) return new Num(0);
            let sum = new Num(0);
            for(let i=1; i<=n.value; i++) {
                sum = new Add(sum, new Div(new Num(1), new Num(i)));
            }
            return sum.simplify();
        }
        // Symbolic sum
        const k = new Sym('__k_harmonic__');
        return this._sum(new Div(new Num(1), k), k, new Num(1), n);
    }

    _lineEquation(p1, p2) {
        if (!(p1 instanceof Vec) || !(p2 instanceof Vec)) throw new Error("Points must be vectors [x, y]");
        if (p1.elements.length < 2 || p2.elements.length < 2) throw new Error("Points must be at least 2D");

        const x1 = p1.elements[0];
        const y1 = p1.elements[1];
        const x2 = p2.elements[0];
        const y2 = p2.elements[1];

        const dx = new Sub(x2, x1).simplify();
        const dy = new Sub(y2, y1).simplify();

        const x = new Sym('x');
        const y = new Sym('y');

        // Check if vertical line (dx = 0)
        let isVertical = false;
        if (dx instanceof Num && dx.value === 0) isVertical = true;

        if (isVertical) {
            // x = x1
            return new Eq(x, x1);
        }

        // slope m = dy/dx
        const m = new Div(dy, dx).simplify();

        // y - y1 = m(x - x1)
        // y = m(x - x1) + y1
        const rhs = new Add(new Mul(m, new Sub(x, x1)), y1).simplify();
        return new Eq(y, rhs);
    }

    _circleEquation(center, radius) {
        if (!(center instanceof Vec)) throw new Error("Center must be vector [x, y]");
        const cx = center.elements[0];
        const cy = center.elements[1];
        const r = radius.simplify();

        const x = new Sym('x');
        const y = new Sym('y');

        // (x - cx)^2 + (y - cy)^2 = r^2
        const lhs = new Add(new Pow(new Sub(x, cx), new Num(2)), new Pow(new Sub(y, cy), new Num(2))).simplify();
        const rhs = new Pow(r, new Num(2)).simplify();

        return new Eq(lhs, rhs);
    }

    _planeEquation(p1, p2, p3) {
        if (!(p1 instanceof Vec) || !(p2 instanceof Vec) || !(p3 instanceof Vec)) throw new Error("Points must be vectors");
        if (p1.elements.length !== 3) throw new Error("Points must be 3D");

        // Vectors u = p2 - p1, v = p3 - p1
        const u = new Vec([
            new Sub(p2.elements[0], p1.elements[0]),
            new Sub(p2.elements[1], p1.elements[1]),
            new Sub(p2.elements[2], p1.elements[2])
        ]).simplify();

        const v = new Vec([
            new Sub(p3.elements[0], p1.elements[0]),
            new Sub(p3.elements[1], p1.elements[1]),
            new Sub(p3.elements[2], p1.elements[2])
        ]).simplify();

        // Normal n = u x v
        const n = this._cross(u, v);

        // Equation: n . (r - p1) = 0 => nx*x + ny*y + nz*z = n . p1
        const nx = n.elements[0];
        const ny = n.elements[1];
        const nz = n.elements[2];

        const x = new Sym('x');
        const y = new Sym('y');
        const z = new Sym('z');

        const lhs = new Add(new Add(new Mul(nx, x), new Mul(ny, y)), new Mul(nz, z)).simplify();

        // rhs = n . p1
        const rhs = new Add(new Add(
            new Mul(nx, p1.elements[0]),
            new Mul(ny, p1.elements[1])),
            new Mul(nz, p1.elements[2])
        ).simplify();

        return new Eq(lhs, rhs);
    }
    _cfrac(x, depth) {
        x = x.simplify();
        if (x instanceof Num) {
            let val = x.value;
            const res = [];
            // Handle depth default
            const d = (depth instanceof Num) ? depth.value : 15;

            for(let i=0; i<d; i++) {
                const a = Math.floor(val);
                res.push(new Num(a));
                if (Math.abs(val - a) < 1e-9) break;
                val = 1 / (val - a);
            }
            return new Vec(res);
        }
        return new Call('cfrac', [x, depth]);
    }

    _isSquare(n) {
        n = n.simplify();
        if (n instanceof Num && Number.isInteger(n.value) && n.value >= 0) {
            const sqrt = Math.sqrt(n.value);
            if (Number.isInteger(sqrt)) return new Num(1);
            return new Num(0);
        }
        return new Call('isSquare', [n]);
    }

    _truthTable(expr, vars) {
        if (!(vars instanceof Vec)) throw new Error("Second argument to truthTable must be a list of variables");
        const variables = vars.elements; // Array of Sym
        const n = variables.length;
        const numRows = 1 << n; // 2^n

        const rows = [];
        for(let i=0; i<numRows; i++) {
            const row = [];
            const subMap = [];

            for(let j=0; j<n; j++) {
                const val = (i >> (n - 1 - j)) & 1;
                const numVal = new Num(val);
                row.push(numVal);
                subMap.push({ v: variables[j], val: numVal });
            }

            let current = expr;
            for(const pair of subMap) {
                current = current.substitute(pair.v, pair.val);
            }
            const res = current.simplify();

            let resVal = res;
            if (res instanceof Sym) {
                if (res.name === 'true') resVal = new Num(1);
                if (res.name === 'false') resVal = new Num(0);
            }

            row.push(resVal);
            rows.push(new Vec(row));
        }

        return new Vec(rows);
    }

    _legendre(n, x) {
        // P_n(x)
        // P_0 = 1, P_1 = x
        // (n+1)P_{n+1} = (2n+1)x P_n - n P_{n-1}
        n = n.simplify();
        if (!(n instanceof Num && Number.isInteger(n.value) && n.value >= 0)) {
            return new Call('legendre', [n, x]);
        }
        const val = n.value;
        if (val === 0) return new Num(1);
        if (val === 1) return x;

        let p_prev = new Num(1);
        let p_curr = x;

        for (let i = 1; i < val; i++) {
            // P_{i+1} = ((2i+1)x P_i - i P_{i-1}) / (i+1)
            const num = new Sub(
                new Mul(new Mul(new Num(2 * i + 1), x), p_curr),
                new Mul(new Num(i), p_prev)
            );
            const p_next = new Div(num, new Num(i + 1)).simplify();
            p_prev = p_curr;
            p_curr = p_next;
        }
        return p_curr;
    }

    _hermite(n, x) {
        // Physicists' Hermite polynomials H_n(x)
        // H_0 = 1, H_1 = 2x
        // H_{n+1} = 2x H_n - 2n H_{n-1}
        n = n.simplify();
        if (!(n instanceof Num && Number.isInteger(n.value) && n.value >= 0)) {
            return new Call('hermite', [n, x]);
        }
        const val = n.value;
        if (val === 0) return new Num(1);
        if (val === 1) return new Mul(new Num(2), x).simplify();

        let h_prev = new Num(1);
        let h_curr = new Mul(new Num(2), x).simplify();

        for (let i = 1; i < val; i++) {
            // H_{i+1} = 2x H_i - 2i H_{i-1}
            const term1 = new Mul(new Mul(new Num(2), x), h_curr);
            const term2 = new Mul(new Mul(new Num(2 * i), new Num(1)), h_prev); // i is number
            const h_next = new Sub(term1, term2).simplify();
            h_prev = h_curr;
            h_curr = h_next;
        }
        return h_curr;
    }

    _chebyshev(n, x) {
        // First kind T_n(x)
        // T_0 = 1, T_1 = x
        // T_{n+1} = 2x T_n - T_{n-1}
        n = n.simplify();
        if (!(n instanceof Num && Number.isInteger(n.value) && n.value >= 0)) {
            return new Call('chebyshev', [n, x]);
        }
        const val = n.value;
        if (val === 0) return new Num(1);
        if (val === 1) return x;

        let t_prev = new Num(1);
        let t_curr = x;

        for (let i = 1; i < val; i++) {
            const term1 = new Mul(new Mul(new Num(2), x), t_curr);
            const t_next = new Sub(term1, t_prev).simplify();
            t_prev = t_curr;
            t_curr = t_next;
        }
        return t_curr;
    }

    _laguerre(n, x) {
        // L_n(x)
        // L_0 = 1, L_1 = 1 - x
        // (n+1)L_{n+1} = (2n+1-x)L_n - nL_{n-1}
        n = n.simplify();
        if (!(n instanceof Num && Number.isInteger(n.value) && n.value >= 0)) {
            return new Call('laguerre', [n, x]);
        }
        const val = n.value;
        if (val === 0) return new Num(1);
        if (val === 1) return new Sub(new Num(1), x).simplify();

        let l_prev = new Num(1);
        let l_curr = new Sub(new Num(1), x).simplify();

        for (let i = 1; i < val; i++) {
            // L_{i+1} = ( (2i+1-x)L_i - i L_{i-1} ) / (i+1)
            const factor = new Sub(new Num(2 * i + 1), x);
            const term1 = new Mul(factor, l_curr);
            const term2 = new Mul(new Num(i), l_prev);
            const num = new Sub(term1, term2);
            const l_next = new Div(num, new Num(i + 1)).simplify();
            l_prev = l_curr;
            l_curr = l_next;
        }
        return l_curr;
    }

    _pinv(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("pinv requires a matrix");
        // Pseudo-inverse A^+ = V * S^+ * U^T
        const [U, S, V] = this._svd(matrix).elements; // svd returns [U, S, V] (V is V, not V^T in my implementation return)

        // S is diagonal matrix
        // Construct S^+ (reciprocal of non-zero elements, transposed - but S is square/diag so transpose is same)
        // Wait, if A is mxn, U is mxm, S is mxn, V is nxn.
        // My SVD implementation returns U, S, V.
        // S is a matrix.
        // We need to invert diagonal elements.

        const rows = S.elements.length;
        const cols = S.elements[0].elements.length;

        const S_plus_rows = [];
        // S^+ should be nxm (transpose dimensions)
        for(let r=0; r<cols; r++) { // New rows = old cols
            const row = [];
            for(let c=0; c<rows; c++) { // New cols = old rows
                if (r === c) { // Diagonal element
                     const val = S.elements[r].elements[c]; // S[r][r] since r==c
                     // Invert if non-zero
                     const numVal = val.evaluateNumeric();
                     if (!isNaN(numVal) && Math.abs(numVal) > 1e-9) {
                         row.push(new Div(new Num(1), val).simplify());
                     } else {
                         row.push(new Num(0));
                     }
                } else {
                     row.push(new Num(0));
                }
            }
            S_plus_rows.push(new Vec(row));
        }
        const S_plus = new Vec(S_plus_rows);

        // A^+ = V * S^+ * U^T
        // Check my SVD return: A = U * S * V^T
        // So A^+ = (V^T)^-1 * S^+ * U^-1 = V * S^+ * U^T (since U, V orthogonal)

        // My _svd returns V. Does it return V or V^T?
        // Code: "Result: [U, S, V] (V, not V^T, following standard svd returns)"
        // And "A = U * S * V^T" implies V contains eigenvectors of A^T A.
        // And result says "Columns are eigenvectors".
        // So V is indeed V.

        const UT = this._trans(U);
        const VSplus = new Mul(V, S_plus).simplify();
        const res = new Mul(VSplus, UT).simplify();

        return res;
    }

    _cond(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("cond requires a matrix");
        // Condition number = max(sigma) / min(sigma)
        // Use SVD
        const svd = this._svd(matrix);
        const S = svd.elements[1]; // Diagonal matrix

        // Extract singular values
        const sigmas = [];
        const rows = S.elements.length;
        const cols = S.elements[0].elements.length;

        const minDim = Math.min(rows, cols);

        for(let i=0; i<minDim; i++) {
             sigmas.push(S.elements[i].elements[i]);
        }

        if (sigmas.length === 0) return new Num(0);

        // SVD logic sorts eigenvalues descending => sigma0 is max.
        // Last one is min.

        const maxSigma = sigmas[0];
        let minSigma = sigmas[sigmas.length - 1];

        // Warning: minSigma might be 0.
        const minVal = minSigma.evaluateNumeric();
        if (Math.abs(minVal) < 1e-12) return new Sym('Infinity');

        return new Div(maxSigma, minSigma).simplify();
    }

    _isDiagonal(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("Argument must be a matrix");
        const rows = matrix.elements.length;
        if (rows === 0) return new Num(1);
        const cols = matrix.elements[0].elements.length;
        if (rows !== cols) return new Num(0); // Must be square

        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                if (i !== j) {
                    const val = matrix.elements[i].elements[j].evaluateNumeric();
                    if (isNaN(val) || Math.abs(val) > 1e-10) {
                        // Check symbolic zero?
                        const s = matrix.elements[i].elements[j].simplify();
                        if (!(s instanceof Num && s.value === 0)) return new Num(0);
                    }
                }
            }
        }
        return new Num(1);
    }

    _isSymmetric(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("Argument must be a matrix");
        const rows = matrix.elements.length;
        if (rows === 0) return new Num(1);
        const cols = matrix.elements[0].elements.length;
        if (rows !== cols) return new Num(0); // Must be square

        const T = this._trans(matrix);
        // Compare matrix and T
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                const diff = new Sub(matrix.elements[i].elements[j], T.elements[i].elements[j]).simplify();
                if (diff instanceof Num && diff.value === 0) continue;
                // If not numeric zero, maybe symbolic zero?
                return new Num(0);
            }
        }
        return new Num(1);
    }

    _isOrthogonal(matrix) {
        if (!(matrix instanceof Vec)) throw new Error("Argument must be a matrix");
        const rows = matrix.elements.length;
        if (rows === 0) return new Num(1);
        const cols = matrix.elements[0].elements.length;
        if (rows !== cols) return new Num(0); // Must be square

        // A^T * A = I
        const T = this._trans(matrix);
        const P = new Mul(T, matrix).simplify();
        const I = this._identity(new Num(rows));

        // Compare P and I
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                const diff = new Sub(P.elements[i].elements[j], I.elements[i].elements[j]).simplify();
                if (diff instanceof Num && Math.abs(diff.value) < 1e-9) continue;
                return new Num(0);
            }
        }
        return new Num(1);
    }

    _collect(expr, varNode) {
        if (!(varNode instanceof Sym)) throw new Error("collect: second argument must be a variable");
        const poly = this._getPolyCoeffs(expr, varNode);
        if (poly) {
            // Reconstruct polynomial from coefficients
            let res = new Num(0);
            // poly.coeffs keys are degrees (strings)
            for (const d in poly.coeffs) {
                const deg = parseInt(d);
                const coeff = poly.coeffs[d];
                let term;
                if (deg === 0) term = coeff;
                else if (deg === 1) term = new Mul(coeff, varNode);
                else term = new Mul(coeff, new Pow(varNode, new Num(deg)));
                res = new Add(res, term);
            }
            return res.simplify();
        }
        return new Call('collect', [expr, varNode]);
    }

    _matrixGen(rowsNode, colsNode, expr, var1, var2) {
        if (!(rowsNode instanceof Num) || !(colsNode instanceof Num)) throw new Error("matrix dimensions must be numbers");
        if (!(var1 instanceof Sym) || !(var2 instanceof Sym)) throw new Error("matrix variables must be symbols");

        const rows = rowsNode.value;
        const cols = colsNode.value;
        const mat = [];

        for(let i=0; i<rows; i++) {
            const row = [];
            for(let j=0; j<cols; j++) {
                // Substitute index variables (0-based)
                // Use temp substitution to avoid collision
                // But simplified substitution should work if i, j are standard
                let val = expr.substitute(var1, new Num(i)).substitute(var2, new Num(j));
                // If var1 == var2, first sub removes both? No, Expr.substitute is specific.
                val = val.simplify();
                row.push(val);
            }
            mat.push(new Vec(row));
        }
        return new Vec(mat);
    }

    _minMaxList(list, type) {
        // list is Vec
        if (list.elements.length === 0) return new Num(NaN);

        // If all numeric, return numeric min/max
        const allNumeric = list.elements.every(e => e instanceof Num);
        if (allNumeric) {
            const vals = list.elements.map(e => e.value);
            return new Num(type === 'min' ? Math.min(...vals) : Math.max(...vals));
        }

        // Return symbolic min(a,b,c...) unpacked
        return new Call(type, list.elements).simplify();
    }

    _volumeSolid(expr, varNode, a, b, axis) {
        if (!(varNode instanceof Sym)) throw new Error("Variable must be a symbol");

        let integrand;
        let axisName = 'x';
        if (axis instanceof Sym) axisName = axis.name;
        if (axis instanceof Num && axis.value === 0) axisName = 'x';

        if (axisName === 'x' || axisName === 'X') {
            integrand = new Mul(new Sym('pi'), new Pow(expr, new Num(2)));
        } else if (axisName === 'y' || axisName === 'Y') {
            integrand = new Mul(new Mul(new Num(2), new Sym('pi')), new Mul(varNode, expr));
        } else {
            integrand = new Mul(new Sym('pi'), new Pow(new Sub(expr, axis), new Num(2)));
        }

        return this.evaluate(new Call('integrate', [integrand, varNode, a, b])).simplify();
    }

    _surfaceArea(expr, varNode, a, b, axis) {
        if (!(varNode instanceof Sym)) throw new Error("Variable must be a symbol");
        // S = int 2*pi * r(x) * sqrt(1 + (f')^2) dx
        // if axis=x, r(x) = f(x) (or |f(x)|)
        // if axis=y, r(x) = x (or |x|)

        let axisName = 'x';
        if (axis instanceof Sym) axisName = axis.name;

        const deriv = expr.diff(varNode).simplify();
        const ds = new Call('sqrt', [new Add(new Num(1), new Pow(deriv, new Num(2)))]).simplify();

        let radius;
        if (axisName === 'x' || axisName === 'X') {
            radius = expr; // Should be abs(expr), but usually expr >= 0 assumed
        } else if (axisName === 'y' || axisName === 'Y') {
            radius = varNode;
        } else {
             // Rotate around y = k? radius = expr - k
             radius = new Sub(expr, axis);
        }

        const integrand = new Mul(new Mul(new Num(2), new Sym('pi')), new Mul(radius, ds));
        return this.evaluate(new Call('integrate', [integrand, varNode, a, b])).simplify();
    }

    _beta(x, y) {
        // B(x, y) = gamma(x)gamma(y) / gamma(x+y)
        const gx = this.evaluate(new Call('gamma', [x]));
        const gy = this.evaluate(new Call('gamma', [y]));
        const gxy = this.evaluate(new Call('gamma', [new Add(x, y)]));
        return new Div(new Mul(gx, gy), gxy).simplify();
    }

    _isSubset(list1, list2) {
        if (!(list1 instanceof Vec) || !(list2 instanceof Vec)) return new Call('isSubset', [list1, list2]);
        // Check if every element of list1 is in list2
        // Simplistic check using toString comparison for symbolic equality
        for (const el1 of list1.elements) {
            const s1 = el1.toString();
            let found = false;
            for (const el2 of list2.elements) {
                if (el2.toString() === s1) {
                    found = true;
                    break;
                }
            }
            if (!found) return new Num(0);
        }
        return new Num(1);
    }

    _cartesianProduct(list1, list2) {
        if (!(list1 instanceof Vec) || !(list2 instanceof Vec)) return new Call('cartesianProduct', [list1, list2]);
        const res = [];
        for (const el1 of list1.elements) {
            for (const el2 of list2.elements) {
                res.push(new Vec([el1, el2]));
            }
        }
        return new Vec(res);
    }

    _integrateSpecial(expr, varNode) {
        // exp(-a*x^2) -> sqrt(pi)/(2*sqrt(a)) * erf(sqrt(a)*x)
        // Check if expr is exp(...)
        if (expr instanceof Call && expr.funcName === 'exp') {
            const arg = expr.args[0]; // -a*x^2
            // Check if arg is quadratic in varNode with no linear term
            const poly = this._getPolyCoeffs(arg, varNode);
            // coeffs[2] = -a, coeffs[1] = 0, coeffs[0] = c (ignore c for now or exp(c) factor)
            if (poly && poly.maxDeg === 2) {
                const c2 = poly.coeffs[2];
                const c1 = poly.coeffs[1] || new Num(0);
                const c0 = poly.coeffs[0] || new Num(0);

                if (c1.evaluateNumeric() === 0 && c2.evaluateNumeric() < 0) {
                     // -a = c2 => a = -c2
                     const a = new Mul(new Num(-1), c2).simplify();
                     const sqrtA = new Call('sqrt', [a]).simplify();

                     // Integral of exp(-ax^2) = sqrt(pi)/(2*sqrt(a)) * erf(sqrt(a)*x)
                     const factor = new Div(new Call('sqrt', [new Sym('pi')]), new Mul(new Num(2), sqrtA));
                     const erfTerm = new Call('erf', [new Mul(sqrtA, varNode)]);

                     let result = new Mul(factor, erfTerm).simplify();

                     // Handle c0: exp(c0) factor
                     if (c0.evaluateNumeric() !== 0) {
                         const C = new Call('exp', [c0]).simplify();
                         result = new Mul(C, result).simplify();
                     }
                     return result;
                }
            }
        }
        return null;
    }

    _sturm(poly, varNode) {
        poly = poly.expand().simplify();
        const seq = [poly];

        let deriv = poly.diff(varNode).simplify();
        if (deriv instanceof Num && deriv.value === 0) return new Vec(seq);

        seq.push(deriv);

        let p_prev2 = poly;
        let p_prev1 = deriv;

        for(let i=0; i<20; i++) {
            if (!this._dependsOn(p_prev1, varNode)) break;

            const rem = this._polyRem(p_prev2, p_prev1, varNode);

            const p_next = new Mul(new Num(-1), rem).simplify();

            if (p_next instanceof Num && p_next.value === 0) break;

            seq.push(p_next);
            p_prev2 = p_prev1;
            p_prev1 = p_next;
        }

        return new Vec(seq);
    }

    _numRealRoots(poly, varNode, a, b) {
        const seq = this._sturm(poly, varNode);

        const signChanges = (val) => {
            let changes = 0;
            let lastSign = 0;

            for(const p of seq.elements) {
                let res;
                const isInf = (v) => v instanceof Sym && (v.name === 'Infinity' || v.name === 'infinity');
                const isNegInf = (v) => v instanceof Mul && v.left instanceof Num && v.left.value === -1 && isInf(v.right);

                if (isInf(val)) {
                    const info = this._getPolyCoeffs(p, varNode);
                    if (info) {
                        const lc = info.coeffs[info.maxDeg].evaluateNumeric();
                        if (!isNaN(lc)) {
                            res = lc > 0 ? 1 : -1;
                        }
                    }
                } else if (isNegInf(val)) {
                    const info = this._getPolyCoeffs(p, varNode);
                    if (info) {
                        const lc = info.coeffs[info.maxDeg].evaluateNumeric();
                        if (!isNaN(lc)) {
                            const signLC = lc > 0 ? 1 : -1;
                            const signPow = (info.maxDeg % 2 === 0) ? 1 : -1;
                            res = signLC * signPow;
                        }
                    }
                } else {
                    res = p.substitute(varNode, val).evaluateNumeric();
                }

                if (typeof res === 'number' && !isNaN(res)) {
                    const sign = res > 0 ? 1 : (res < 0 ? -1 : 0);
                    if (sign !== 0) {
                        if (lastSign !== 0 && sign !== lastSign) {
                            changes++;
                        }
                        lastSign = sign;
                    }
                }
            }
            return changes;
        };

        const Va = signChanges(a);
        const Vb = signChanges(b);

        return new Num(Math.abs(Va - Vb));
    }

    _getPrimeFactorsMap(n) {
        const factors = this._primeFactors(n); // Vec of Num
        const map = {};
        if (factors instanceof Vec) {
            factors.elements.forEach(f => {
                if (f instanceof Num) {
                    const val = f.value;
                    map[val] = (map[val] || 0) + 1;
                }
            });
        }
        return map;
    }

    _moebius(n) {
        n = n.simplify();
        if (!(n instanceof Num && Number.isInteger(n.value) && n.value > 0)) {
            return new Call('moebius', [n]);
        }
        if (n.value === 1) return new Num(1);

        const map = this._getPrimeFactorsMap(n);
        let count = 0;
        for (const p in map) {
            if (map[p] > 1) return new Num(0); // Square factor
            count++;
        }
        return new Num(count % 2 === 0 ? 1 : -1);
    }

    _sigma(n, k) {
        // sum of k-th powers of divisors
        // sigma_k(n) = prod ( (p^(k(e+1)) - 1) / (p^k - 1) )
        n = n.simplify();
        k = k.simplify();
        if (!(n instanceof Num && Number.isInteger(n.value) && n.value > 0) || !(k instanceof Num)) {
            return new Call('sigma', [n, k]);
        }

        const map = this._getPrimeFactorsMap(n);
        let res = new Num(1);
        const kVal = k.value;

        for (const pStr in map) {
            const p = parseInt(pStr);
            const e = map[pStr];
            if (kVal === 0) {
                // Count divisors: prod(e+1)
                res = new Mul(res, new Num(e + 1));
            } else {
                // Geometric series sum: 1 + p^k + ... + p^(ke)
                // = (p^(k(e+1)) - 1) / (p^k - 1)
                const num = Math.pow(p, kVal * (e + 1)) - 1;
                const den = Math.pow(p, kVal) - 1;
                res = new Mul(res, new Num(num / den));
            }
        }
        return res.simplify();
    }

    _legendreSymbol(a, p) {
        a = a.simplify();
        p = p.simplify();
        if (!(a instanceof Num && Number.isInteger(a.value)) || !(p instanceof Num && Number.isInteger(p.value))) {
            return new Call('legendreSymbol', [a, p]);
        }
        // (a/p) = a^((p-1)/2) mod p
        // returns 1, -1, 0
        const A = a.value;
        const P = p.value;
        if (P <= 2) throw new Error("legendreSymbol requires odd prime p > 2"); // Or p=2 trivial

        if (A % P === 0) return new Num(0);

        // Euler's criterion
        const exp = (P - 1) / 2;
        const res = this._modPow(new Num(A), new Num(exp), new Num(P));

        // Result is in [0, P-1]. If P-1, return -1.
        if (res.value === P - 1) return new Num(-1);
        return res;
    }

    _stirling1(n, k) {
        // Signed Stirling numbers of first kind s(n, k)
        // Coeff of x^k in x(x-1)...(x-n+1) (falling factorial x_n)
        // s(n, k) = s(n-1, k-1) - (n-1)s(n-1, k)
        n = n.simplify();
        k = k.simplify();
        if (!(n instanceof Num && Number.isInteger(n.value)) || !(k instanceof Num && Number.isInteger(k.value))) {
            return new Call('stirling1', [n, k]);
        }
        const N = n.value;
        const K = k.value;
        if (K < 0 || K > N) return new Num(0);
        if (N === 0 && K === 0) return new Num(1);
        if (N === 0 || K === 0) return new Num(0);

        // DP table
        // Map "n,k" -> val
        const memo = {};
        const s = (n, k) => {
            if (k === n) return 1;
            if (k === 0) return 0; // n>0
            const key = n + "," + k;
            if (memo[key] !== undefined) return memo[key];
            const res = s(n - 1, k - 1) - (n - 1) * s(n - 1, k);
            memo[key] = res;
            return res;
        };
        return new Num(s(N, K));
    }

    _stirling2(n, k) {
        // Stirling numbers of second kind S(n, k)
        // S(n, k) = S(n-1, k-1) + k*S(n-1, k)
        n = n.simplify();
        k = k.simplify();
        if (!(n instanceof Num && Number.isInteger(n.value)) || !(k instanceof Num && Number.isInteger(k.value))) {
            return new Call('stirling2', [n, k]);
        }
        const N = n.value;
        const K = k.value;
        if (K < 0 || K > N) return new Num(0);
        if (N === 0 && K === 0) return new Num(1);
        if (N === 0 || K === 0) return new Num(0);

        const memo = {};
        const S = (n, k) => {
            if (k === n) return 1;
            if (k === 0) return 0;
            const key = n + "," + k;
            if (memo[key] !== undefined) return memo[key];
            const res = S(n - 1, k - 1) + k * S(n - 1, k);
            memo[key] = res;
            return res;
        };
        return new Num(S(N, K));
    }

    _bell(n) {
        n = n.simplify();
        if (!(n instanceof Num && Number.isInteger(n.value) && n.value >= 0)) {
            return new Call('bell', [n]);
        }
        const N = n.value;
        let sum = new Num(0);
        for (let k = 0; k <= N; k++) {
            sum = new Add(sum, this._stirling2(new Num(N), new Num(k)));
        }
        return sum.simplify();
    }

    _isPerfect(n) {
        n = n.simplify();
        if (!(n instanceof Num)) return new Call('isPerfect', [n]);
        // sigma(n) == 2n
        const sig = this._sigma(n, new Num(1));
        const twoN = new Num(2 * n.value);
        return new Num(sig.value === twoN.value ? 1 : 0);
    }

    _matrixPow(matrix, n) {
        if (!(matrix instanceof Vec)) throw new Error("matrixPow requires a matrix");
        n = n.simplify();

        // Integer power handled by loop or binary exponentiation
        if (n instanceof Num && Number.isInteger(n.value)) {
            const exp = n.value;
            if (exp === 0) {
                // Identity
                return this._identity(new Num(matrix.elements.length));
            }
            if (exp === -1) return this._inv(matrix);
            if (exp < 0) return this._matrixPow(this._inv(matrix), new Num(-exp));

            // Binary Exponentiation
            let res = this._identity(new Num(matrix.elements.length));
            let base = matrix;
            let e = exp;
            while (e > 0) {
                if (e % 2 === 1) res = new Mul(res, base).simplify();
                base = new Mul(base, base).simplify();
                e = Math.floor(e / 2);
            }
            return res;
        }

        // Fractional/Symbolic power: Use diagonalization P D^n P^-1
        // Requires matrix to be diagonalizable
        try {
            const [P, D] = this._diagonalize(matrix).elements;
            const D_pow_rows = [];
            const rows = D.elements.length;
            for(let i=0; i<rows; i++) {
                const row = [];
                for(let j=0; j<rows; j++) {
                    if (i === j) {
                        const val = D.elements[i].elements[j];
                        row.push(new Pow(val, n).simplify());
                    } else {
                        row.push(new Num(0));
                    }
                }
                D_pow_rows.push(new Vec(row));
            }
            const D_pow = new Vec(D_pow_rows);

            const P_inv = this._inv(P);
            const PD = new Mul(P, D_pow).simplify();
            return new Mul(PD, P_inv).simplify();
        } catch (e) {
            return new Call('matrixPow', [matrix, n]);
        }
    }

    _kroneckerDelta(i, j) {
        i = i.simplify();
        j = j.simplify();
        if (i instanceof Num && j instanceof Num) {
            return new Num(i.value === j.value ? 1 : 0);
        }
        if (i.toString() === j.toString()) return new Num(1);
        return new Call('kroneckerDelta', [i, j]);
    }

    _fft(list) {
        if (!(list instanceof Vec)) throw new Error("fft requires a list");
        const n = list.elements.length;
        if ((n & (n - 1)) !== 0 || n === 0) throw new Error("fft requires length to be a power of 2");

        // Convert elements to complex {re, im}
        const input = list.elements.map(e => {
            const c = this._getComplexParts(e);
            return { re: c.re.evaluateNumeric(), im: c.im.evaluateNumeric() };
        });

        // Check if numeric
        if (input.some(x => isNaN(x.re) || isNaN(x.im))) {
            return new Call('fft', [list]);
        }

        const fftRec = (x) => {
            const N = x.length;
            if (N <= 1) return x;
            const even = fftRec(x.filter((_, i) => i % 2 === 0));
            const odd = fftRec(x.filter((_, i) => i % 2 !== 0));
            const T = new Array(N / 2);
            for (let k = 0; k < N / 2; k++) {
                const angle = -2 * Math.PI * k / N;
                const wk = { re: Math.cos(angle), im: Math.sin(angle) };
                // odd[k] * wk
                const t = {
                    re: wk.re * odd[k].re - wk.im * odd[k].im,
                    im: wk.re * odd[k].im + wk.im * odd[k].re
                };
                T[k] = t;
            }
            const res = new Array(N);
            for (let k = 0; k < N / 2; k++) {
                res[k] = { re: even[k].re + T[k].re, im: even[k].im + T[k].im };
                res[k + N / 2] = { re: even[k].re - T[k].re, im: even[k].im - T[k].im };
            }
            return res;
        };

        const output = fftRec(input);
        // Convert back to Expr
        return new Vec(output.map(c => {
            // Simplify near-zero
            const re = (Math.abs(c.re) < 1e-9) ? 0 : c.re;
            const im = (Math.abs(c.im) < 1e-9) ? 0 : c.im;
            if (im === 0) return new Num(re);
            if (re === 0) return new Mul(new Num(im), new Sym('i'));
            return new Add(new Num(re), new Mul(new Num(im), new Sym('i')));
        }));
    }

    _ifft(list) {
        if (!(list instanceof Vec)) throw new Error("ifft requires a list");
        const n = list.elements.length;
        if ((n & (n - 1)) !== 0 || n === 0) throw new Error("ifft requires length to be a power of 2");

        // IFFT(x) = 1/N * conj(FFT(conj(x)))
        const conjList = new Vec(list.elements.map(e => new Call('conj', [e]).simplify()));
        const fftRes = this._fft(conjList);

        if (fftRes instanceof Call) return new Call('ifft', [list]);

        const result = new Vec(fftRes.elements.map(e => {
            // conj(e) / N
            const c = new Call('conj', [e]).simplify();
            return new Div(c, new Num(n)).simplify();
        }));
        return result;
    }

    _cnf(expr) {
        // Convert to Conjunctive Normal Form (AND of ORs)
        // 1. Eliminate implications (already done by simplifier usually, but let's ensure)
        // 2. Move NOT inwards (De Morgan)
        // 3. Distribute OR over AND: A or (B and C) -> (A or B) and (A or C)

        let current = expr;

        // Helper: push NOT inwards
        const pushNot = (e) => {
            if (e instanceof Not) {
                const arg = e.arg;
                if (arg instanceof Not) return pushNot(arg.arg); // Double neg
                if (arg instanceof And) return new Or(pushNot(new Not(arg.left)), pushNot(new Not(arg.right)));
                if (arg instanceof Or) return new And(pushNot(new Not(arg.left)), pushNot(new Not(arg.right)));
                return e;
            }
            if (e instanceof BinaryOp) {
                // Recurse
                // Note: BinaryOp constructor is abstract, use e.constructor
                return new e.constructor(pushNot(e.left), pushNot(e.right));
            }
            return e;
        };

        // Helper: distribute OR over AND
        const distribute = (e) => {
            if (e instanceof And) {
                return new And(distribute(e.left), distribute(e.right));
            }
            if (e instanceof Or) {
                const l = distribute(e.left);
                const r = distribute(e.right);
                // (P and Q) or R -> (P or R) and (Q or R)
                if (l instanceof And) {
                    return new And(distribute(new Or(l.left, r)), distribute(new Or(l.right, r)));
                }
                // L or (P and Q) -> (L or P) and (L or Q)
                if (r instanceof And) {
                    return new And(distribute(new Or(l, r.left)), distribute(new Or(l, r.right)));
                }
                return new Or(l, r);
            }
            return e;
        };

        // Apply
        // Simplify first to handle Implies/Iff -> Or/And/Not
        let s = current.simplify();
        // Expand logical operators if they are not basic And/Or/Not?
        // simplify() handles implies/iff -> logic ops?
        // Implies.simplify() -> Implies (unless numeric).
        // We need transformation: A->B to !A or B.
        const eliminate = (e) => {
            if (e instanceof Implies) return new Or(new Not(eliminate(e.left)), eliminate(e.right));
            if (e instanceof Iff) {
                // (A->B) and (B->A)
                const A = eliminate(e.left);
                const B = eliminate(e.right);
                return new And(new Or(new Not(A), B), new Or(new Not(B), A));
            }
            if (e instanceof Xor) {
                // (A or B) and not (A and B)
                const A = eliminate(e.left);
                const B = eliminate(e.right);
                return new And(new Or(A, B), new Not(new And(A, B)));
            }
            if (e instanceof BinaryOp) return new e.constructor(eliminate(e.left), eliminate(e.right));
            if (e instanceof Not) return new Not(eliminate(e.arg));
            return e;
        };

        s = eliminate(s);
        s = pushNot(s); // simplify handles some, but be explicit
        // Repeat distribute until stable? Or just once recursive?
        // Distribute is recursive.
        return distribute(s);
    }

    _dnf(expr) {
        // Disjunctive Normal Form (OR of ANDs)
        // Distribute AND over OR: A and (B or C) -> (A and B) or (A and C)

        const eliminate = (e) => {
            if (e instanceof Implies) return new Or(new Not(eliminate(e.left)), eliminate(e.right));
            if (e instanceof Iff) {
                // (A and B) or (!A and !B)
                const A = eliminate(e.left);
                const B = eliminate(e.right);
                return new Or(new And(A, B), new And(new Not(A), new Not(B)));
            }
            if (e instanceof Xor) {
                // (A and !B) or (!A and B)
                const A = eliminate(e.left);
                const B = eliminate(e.right);
                return new Or(new And(A, new Not(B)), new And(new Not(A), B));
            }
            if (e instanceof BinaryOp) return new e.constructor(eliminate(e.left), eliminate(e.right));
            if (e instanceof Not) return new Not(eliminate(e.arg));
            return e;
        };

        const pushNot = (e) => {
            if (e instanceof Not) {
                const arg = e.arg;
                if (arg instanceof Not) return pushNot(arg.arg);
                if (arg instanceof And) return new Or(pushNot(new Not(arg.left)), pushNot(new Not(arg.right)));
                if (arg instanceof Or) return new And(pushNot(new Not(arg.left)), pushNot(new Not(arg.right)));
                return e;
            }
            if (e instanceof BinaryOp) return new e.constructor(pushNot(e.left), pushNot(e.right));
            return e;
        };

        const distribute = (e) => {
            if (e instanceof Or) {
                return new Or(distribute(e.left), distribute(e.right));
            }
            if (e instanceof And) {
                const l = distribute(e.left);
                const r = distribute(e.right);
                // A and (B or C) -> (A and B) or (A and C)
                if (l instanceof Or) {
                    return new Or(distribute(new And(l.left, r)), distribute(new And(l.right, r)));
                }
                if (r instanceof Or) {
                    return new Or(distribute(new And(l, r.left)), distribute(new And(l, r.right)));
                }
                return new And(l, r);
            }
            return e;
        };

        let s = eliminate(expr);
        s = pushNot(s);
        return distribute(s);
    }

    _annuity(r, n, p) {
        // FV = P * ((1+r)^n - 1) / r
        const onePlusR = new Add(new Num(1), r);
        const term = new Pow(onePlusR, n);
        const num = new Sub(term, new Num(1));
        const frac = new Div(num, r);
        return new Mul(p, frac).simplify();
    }

    _amortization(r, n, pv) {
        // P = r * PV / (1 - (1+r)^-n)
        const onePlusR = new Add(new Num(1), r);
        const term = new Pow(onePlusR, new Mul(new Num(-1), n));
        const den = new Sub(new Num(1), term);
        const num = new Mul(r, pv);
        return new Div(num, den).simplify();
    }

    _codegen(expr, lang) {
        lang = lang.toLowerCase();

        const mapOp = (op) => {
            if (lang === 'python' || lang === 'py') {
                if (op === '^') return '**';
                if (op === '&&') return 'and';
                if (op === '||') return 'or';
                if (op === '!') return 'not ';
            } else if (lang === 'js' || lang === 'javascript') {
                if (op === '^') return '**'; // ES6
                if (op === '&&') return '&&';
                if (op === '||') return '||';
                if (op === '!') return '!';
            } else if (lang === 'c' || lang === 'cpp') {
                if (op === '^') return ', '; // pow(a, b) needs comma
                if (op === '&&') return '&&';
                if (op === '||') return '||';
                if (op === '!') return '!';
            } else if (lang === 'matlab' || lang === 'octave') {
                if (op === '^') return '^';
                if (op === '&&') return '&&';
                if (op === '||') return '||';
                if (op === '!') return '~';
            }
            return op;
        };

        const mapFunc = (func) => {
            if (lang === 'python' || lang === 'py') {
                const np = ['sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh', 'exp', 'log', 'sqrt', 'abs', 'floor', 'ceil', 'erf', 'gamma'];
                if (np.includes(func)) return 'np.' + func;
                if (func === 'ln') return 'np.log';
                if (func === 'arctan') return 'np.arctan';
                if (func === 'arcsin') return 'np.arcsin';
                if (func === 'arccos') return 'np.arccos';
            } else if (lang === 'js' || lang === 'javascript') {
                const math = ['sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'exp', 'log', 'sqrt', 'abs', 'floor', 'ceil', 'round', 'min', 'max', 'pow', 'sinh', 'cosh', 'tanh', 'asinh', 'acosh', 'atanh', 'cbrt', 'sign', 'trunc'];
                if (math.includes(func)) return 'Math.' + func;
                if (func === 'ln') return 'Math.log';
            } else if (lang === 'c' || lang === 'cpp') {
                if (func === 'ln') return 'log';
                if (func === 'abs') return 'fabs';
                const std = ['sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh', 'exp', 'sqrt', 'floor', 'ceil', 'erf', 'tgamma'];
                if (func === 'gamma') return 'tgamma';
                if (std.includes(func)) return func;
            } else if (lang === 'matlab' || lang === 'octave') {
                if (func === 'ln') return 'log';
                if (func === 'log') return 'log10';
            }
            return func;
        };

        const mapConst = (c) => {
            if (lang === 'python' || lang === 'py') {
                if (c === 'pi') return 'np.pi';
                if (c === 'e') return 'np.e';
                if (c === 'i') return '1j';
                if (c === 'true') return 'True';
                if (c === 'false') return 'False';
                if (c === 'infinity' || c === 'Infinity') return 'np.inf';
            } else if (lang === 'js' || lang === 'javascript') {
                if (c === 'pi') return 'Math.PI';
                if (c === 'e') return 'Math.E';
                if (c === 'true') return 'true';
                if (c === 'false') return 'false';
                if (c === 'infinity' || c === 'Infinity') return 'Infinity';
            } else if (lang === 'c' || lang === 'cpp') {
                if (c === 'pi') return 'M_PI';
                if (c === 'e') return 'M_E';
                if (c === 'true') return '1';
                if (c === 'false') return '0';
                if (c === 'infinity' || c === 'Infinity') return 'INFINITY';
            }
            return c;
        };

        const rec = (node) => {
            if (node instanceof Num) return node.toString();
            if (node instanceof Sym) return mapConst(node.name);

            if (node instanceof Add) return `(${rec(node.left)} + ${rec(node.right)})`;
            if (node instanceof Sub) return `(${rec(node.left)} - ${rec(node.right)})`;
            if (node instanceof Mul) return `(${rec(node.left)} * ${rec(node.right)})`;
            if (node instanceof Div) return `(${rec(node.left)} / ${rec(node.right)})`;

            if (node instanceof Pow) {
                if ((lang === 'c' || lang === 'cpp')) {
                    return `pow(${rec(node.left)}, ${rec(node.right)})`;
                }
                return `(${rec(node.left)} ${mapOp('^')} ${rec(node.right)})`;
            }

            if (node instanceof Call) {
                const fn = mapFunc(node.funcName);
                const args = node.args.map(rec).join(', ');
                return `${fn}(${args})`;
            }

            if (node instanceof Eq) {
                if (lang === 'python' || lang === 'c' || lang === 'js') return `(${rec(node.left)} == ${rec(node.right)})`;
            }

            if (node instanceof Not) return `(${mapOp('!')}${rec(node.arg)})`;
            if (node instanceof And) return `(${rec(node.left)} ${mapOp('&&')} ${rec(node.right)})`;
            if (node instanceof Or) return `(${rec(node.left)} ${mapOp('||')} ${rec(node.right)})`;

            return node.toString();
        };

        return { type: 'info', text: rec(expr), toLatex: () => `\\texttt{${rec(expr)}}` };
    }

    _convert(val, from, to) {
        val = val.evaluateNumeric();
        if (isNaN(val)) return new Sym('NaN');

        // Check Temp
        if (from === 'C' || from === 'F' || from === 'K' || to === 'C' || to === 'F' || to === 'K') {
            let k = val;
            if (from === 'C') k = val + 273.15;
            else if (from === 'F') k = (val - 32) * 5/9 + 273.15;

            let res = k;
            if (to === 'C') res = k - 273.15;
            else if (to === 'F') res = (k - 273.15) * 9/5 + 32;

            return new Num(res);
        }

        const rateFrom = this.conversionRates[from];
        const rateTo = this.conversionRates[to];

        if (rateFrom === undefined || rateTo === undefined) {
            throw new Error(`Unknown unit conversion: ${from} -> ${to}`);
        }

        // Base value = val * rateFrom
        // Target value = Base / rateTo
        // Note: assumes same dimension check is implicit or user knows.
        // We could check if they share a base category but I flattened the map.

        return new Num(val * rateFrom / rateTo);
    }

    _fourierTransform(expr, t, w) {
        expr = expr.expand().simplify();

        if (expr instanceof Add) return new Add(this._fourierTransform(expr.left, t, w), this._fourierTransform(expr.right, t, w)).simplify();
        if (expr instanceof Sub) return new Sub(this._fourierTransform(expr.left, t, w), this._fourierTransform(expr.right, t, w)).simplify();
        if (expr instanceof Mul) {
            if (!this._dependsOn(expr.left, t)) return new Mul(expr.left, this._fourierTransform(expr.right, t, w)).simplify();
            if (!this._dependsOn(expr.right, t)) return new Mul(expr.right, this._fourierTransform(expr.left, t, w)).simplify();
        }

        // Table
        // delta(t - t0) -> exp(-i*w*t0)
        // exp(-a*t^2) -> sqrt(pi/a) * exp(-w^2/(4a))
        // exp(-a*abs(t)) -> 2a / (a^2 + w^2)
        // 1 -> 2*pi*delta(w) (Not handled nicely, but possible)

        // Gaussian: exp(-a*t^2)
        if (expr instanceof Call && expr.funcName === 'exp') {
             const arg = expr.args[0];
             // check -a*t^2
             const poly = this._getPolyCoeffs(arg, t);
             if (poly && poly.maxDeg === 2) {
                 const c2 = poly.coeffs[2];
                 // const c1 = poly.coeffs[1]; // Shift?
                 if (poly.coeffs[1] && poly.coeffs[1].evaluateNumeric() !== 0) {
                     // Shift theorem? F{f(t-t0)}? No, this is exp(-(t-t0)^2) maybe
                     return new Call('fourier_transform', [expr, t, w]);
                 }

                 // Check c2 is negative
                 const c2Val = c2.evaluateNumeric();
                 if (!isNaN(c2Val) && c2Val < 0) {
                     // a = -c2
                     const a = new Mul(new Num(-1), c2).simplify();

                     // sqrt(pi/a) * exp(-w^2 / 4a)
                     const factor = new Call('sqrt', [new Div(new Sym('pi'), a)]);
                     const expArg = new Div(new Mul(new Num(-1), new Pow(w, new Num(2))), new Mul(new Num(4), a));
                     return new Mul(factor, new Call('exp', [expArg])).simplify();
                 }
             }
        }

        // dirac(t - a)
        if (expr instanceof Call && expr.funcName === 'dirac') {
            const arg = expr.args[0];
            const poly = this._getPolyCoeffs(arg, t);
            if (poly && poly.maxDeg === 1 && poly.coeffs[1].evaluateNumeric() === 1) {
                // t + c0 = t - (-c0)
                const a = new Mul(new Num(-1), poly.coeffs[0] || new Num(0)).simplify();
                // e^(-i w a)
                return new Call('exp', [new Mul(new Mul(new Num(-1), new Sym('i')), new Mul(w, a))]).simplify();
            }
        }

        return new Call('fourier_transform', [expr, t, w]);
    }

    _inverseFourierTransform(expr, w, t) {
        // Just inverse of the above
        return new Call('inverse_fourier_transform', [expr, w, t]);
    }

    _mse(list) {
        if (!(list instanceof Vec)) throw new Error("mse requires a list");
        // mean squared error (from 0? or mean?)
        // If single list, usually MSE of errors implies elements ARE errors?
        // Or MSE of data from mean? That's variance (biased).
        // Let's assume list is errors. sum(x^2)/n
        // Or if list of points? No, usually mse(residuals).

        const n = list.elements.length;
        if (n === 0) return new Num(0);
        let sum = new Num(0);
        for(const e of list.elements) {
            sum = new Add(sum, new Pow(e, new Num(2)));
        }
        return new Div(sum.simplify(), new Num(n)).simplify();
    }

    _mae(list) {
        if (!(list instanceof Vec)) throw new Error("mae requires a list");
        // mean absolute error
        const n = list.elements.length;
        if (n === 0) return new Num(0);
        let sum = new Num(0);
        for(const e of list.elements) {
            sum = new Add(sum, new Call('abs', [e]));
        }
        return new Div(sum.simplify(), new Num(n)).simplify();
    }

    _convolution(f, g, t) {
        if (!(t instanceof Sym)) throw new Error("Third argument to convolution must be a variable (t)");
        // Integral from 0 to t of f(tau) * g(t - tau) d_tau
        const tau = new Sym('tau_' + Math.floor(Math.random() * 1000)); // Unique var
        const f_tau = f.substitute(t, tau);
        const t_minus_tau = new Sub(t, tau);
        const g_shifted = g.substitute(t, t_minus_tau);
        const integrand = new Mul(f_tau, g_shifted);

        return this.evaluate(new Call('integrate', [integrand, tau, new Num(0), t])).simplify();
    }

    _conv(u, v) {
        if (!(u instanceof Vec) || !(v instanceof Vec)) throw new Error("Arguments must be lists/vectors");
        const n = u.elements.length;
        const m = v.elements.length;
        const len = n + m - 1;
        const res = [];

        for(let k = 0; k < len; k++) {
            let sum = new Num(0);
            // sum_{j} u[j] * v[k-j]
            // range of j: max(0, k-(m-1)) to min(n-1, k)
            const lower = Math.max(0, k - (m - 1));
            const upper = Math.min(n - 1, k);

            for(let j = lower; j <= upper; j++) {
                sum = new Add(sum, new Mul(u.elements[j], v.elements[k-j])).simplify();
            }
            res.push(sum);
        }
        return new Vec(res);
    }

    _xcorr(u, v) {
        if (!(u instanceof Vec) || !(v instanceof Vec)) throw new Error("Arguments must be lists/vectors");
        // Cross-correlation: (f star g)[n] = sum f[m] * conj(g[m+n])
        // Or standard signal processing def?
        // xcorr(u, v) usually implies full correlation (lags -(N-1) to (M-1))?
        // Or valid/same?
        // Let's implement full cross-correlation via convolution.
        // xcorr(u, v) = conv(u, reverse(conj(v)))
        // but convolution index shift?
        // Matlab xcorr returns vector of length 2*max(N,M)-1 centered at lag 0.
        // Let's implement simple discrete correlation.
        // conv(u, reverse(v))
        const vRev = this._reverse(v);
        // Conj if complex?
        const vConjRev = new Vec(vRev.elements.map(e => new Call('conj', [e]).simplify()));
        return this._conv(u, vConjRev);
    }

    _pade(expr, varNode, n, m) {
        // Pade Approximation R(x) = P(x) / Q(x)
        // deg(P) <= n, deg(Q) <= m, Q(0)=1
        // expr - P/Q = O(x^(n+m+1)) => expr * Q - P = O(x^(n+m+1))
        if (!(n instanceof Num) || !(m instanceof Num)) throw new Error("Degrees must be numbers");
        const N = n.value;
        const M = m.value;
        const order = N + M;

        // 1. Taylor Series
        const taylor = this._taylor(expr, varNode, new Num(0), order);
        const poly = this._getPolyCoeffs(taylor, varNode); // coeffs[k] is coeff of x^k
        if (!poly) throw new Error("Could not compute Taylor series");

        const c = [];
        for(let i=0; i<=order; i++) c[i] = poly.coeffs[i] || new Num(0);

        // 2. Solve for Q coefficients q_1 ... q_m (q_0 = 1)
        // System: sum_{j=0}^M c_{k-j} * q_j = 0 for k = N+1 ... N+M
        // q_0 = 1
        // c_k + c_{k-1}q_1 + ... + c_{k-m}q_m = 0
        // Matrix A * [q_1 ... q_m]^T = B

        if (M > 0) {
            const rows = [];
            for(let i=0; i<M; i++) {
                const k = N + 1 + i;
                const row = [];
                for(let j=1; j<=M; j++) {
                    // coeff of q_j is c_{k-j}
                    const idx = k - j;
                    if (idx < 0) row.push(new Num(0));
                    else row.push(c[idx] || new Num(0));
                }
                rows.push(new Vec(row));
            }
            const A = new Vec(rows);

            const b = [];
            for(let i=0; i<M; i++) {
                const k = N + 1 + i;
                // RHS is -c_k * q_0 = -c_k
                b.push(new Mul(new Num(-1), c[k] || new Num(0)).simplify());
            }
            const B = new Vec(b.map(x => new Vec([x]))); // Col vector

            // Solve A*q = B
            // If A is singular, Pade might not exist or non-unique.
            let q_vec;
            try {
                // Use linear solve
                // _solveSystem or matrix solve?
                // A*x=B -> x = A^-1 * B or solve(Ax=B)
                // Use _inv?
                const invA = this._inv(A);
                q_vec = new Mul(invA, B).simplify();
            } catch (e) {
                // Fallback or singular
                throw new Error("Pade approximation failed (singular matrix)");
            }

            // Construct Q
            let Q = new Num(1);
            for(let j=1; j<=M; j++) {
                const qj = q_vec.elements[j-1].elements[0];
                Q = new Add(Q, new Mul(qj, new Pow(varNode, new Num(j)))).simplify();
            }

            // 3. Compute P
            // P(x) = (Q(x) * T(x)) mod x^(N+1)
            // p_k = sum_{j=0}^k c_{k-j} * q_j for k=0..N
            let P = new Num(0);
            for(let k=0; k<=N; k++) {
                let pk = new Num(0);
                for(let j=0; j<=k; j++) {
                    // q_j
                    let qj;
                    if (j===0) qj = new Num(1);
                    else if (j <= M) qj = q_vec.elements[j-1].elements[0];
                    else qj = new Num(0);

                    const ckj = c[k-j] || new Num(0);
                    pk = new Add(pk, new Mul(ckj, qj)).simplify();
                }
                P = new Add(P, new Mul(pk, new Pow(varNode, new Num(k)))).simplify();
            }

            return new Div(P, Q).simplify();
        } else {
            // M=0, just Taylor polynomial
            return taylor;
        }
    }

    _vandermonde(vec) {
        if (!(vec instanceof Vec)) throw new Error("vandermonde argument must be a vector");
        const n = vec.elements.length;
        const rows = [];
        for(let i=0; i<n; i++) {
            const row = [];
            const x = vec.elements[i];
            for(let j=0; j<n; j++) {
                // x^j
                row.push(new Pow(x, new Num(j)).simplify());
            }
            rows.push(new Vec(row));
        }
        return new Vec(rows);
    }

    _hilbert(n) {
        if (!(n instanceof Num)) throw new Error("hilbert size must be a number");
        const size = n.value;
        const rows = [];
        for(let i=1; i<=size; i++) {
            const row = [];
            for(let j=1; j<=size; j++) {
                // H_ij = 1 / (i + j - 1)
                row.push(new Div(new Num(1), new Num(i + j - 1)).simplify());
            }
            rows.push(new Vec(row));
        }
        return new Vec(rows);
    }

    _toeplitz(c, r) {
        if (!(c instanceof Vec)) throw new Error("Column argument must be a vector");
        if (r && !(r instanceof Vec)) throw new Error("Row argument must be a vector");
        // Symmetric Toeplitz if r is null (uses conj of c)
        // c defines first col. r defines first row.
        // c[0] must equal r[0] usually.

        const m = c.elements.length; // rows
        const n = r ? r.elements.length : m; // cols

        const rows = [];
        for(let i=0; i<m; i++) {
            const row = [];
            for(let j=0; j<n; j++) {
                if (j >= i) {
                    // Upper triangle (use r)
                    const idx = j - i;
                    if (r) {
                        row.push(r.elements[idx]);
                    } else {
                        // Symmetric/Hermitian: conj(c[idx])
                        row.push(new Call('conj', [c.elements[idx]]).simplify());
                    }
                } else {
                    // Lower triangle (use c)
                    const idx = i - j;
                    row.push(c.elements[idx]);
                }
            }
            rows.push(new Vec(row));
        }
        return new Vec(rows);
    }

    _wronskian(funcs, varNode) {
        if (!(funcs instanceof Vec)) throw new Error("wronskian expects a list of functions");
        if (!(varNode instanceof Sym)) throw new Error("wronskian expects a variable");

        const n = funcs.elements.length;
        const rows = [];

        // Row 0: functions themselves
        let currentDerivs = funcs.elements;
        rows.push(new Vec(currentDerivs));

        // Subsequent rows: derivatives
        for (let i = 1; i < n; i++) {
            const nextDerivs = currentDerivs.map(f => f.diff(varNode).simplify());
            rows.push(new Vec(nextDerivs));
            currentDerivs = nextDerivs;
        }

        const mat = new Vec(rows);
        return this._det(mat);
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
