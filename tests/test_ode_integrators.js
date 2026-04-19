// test_ode_integrators.js — ODE integrators: rk4 + rk45 (Dormand-Prince)
// Run: node tests/test_ode_integrators.js

const { CAS } = require('../js/cas.js');
const { Expr, Num, Sym, Call, Vec, Eq, Add, Sub, Mul, Div, Pow } = require('../js/expression.js');

const cas = new CAS();

let passed = 0, failed = 0;

function assertApprox(actual, expected, label, tolerance) {
    tolerance = tolerance || 1e-6;
    if (Math.abs(actual - expected) < tolerance) {
        console.log(`  ✓ ${label}: ${actual.toFixed(8)} (expected ≈ ${expected.toFixed(8)})`);
        passed++;
    } else {
        console.log(`  ✗ ${label}: ${actual.toFixed(8)} (expected ≈ ${expected.toFixed(8)}, diff=${Math.abs(actual-expected).toExponential(2)})`);
        failed++;
    }
}

function assertTrue(condition, label) {
    if (condition) {
        console.log(`  ✓ ${label}`);
        passed++;
    } else {
        console.log(`  ✗ ${label}`);
        failed++;
    }
}

console.log('\n=== Test: ODE Integrators (RK4 + RK45/Dormand-Prince) ===\n');

// ── 1. Exponential decay: y' = -0.5*y, y(0)=1 ───────────────────────────────
console.log('--- Problem 1: Exponential decay (y\' = -0.5*y, y(0)=1) ---');

{
    // RK4 with fixed step
    const y = new Sym('y');
    const t = new Sym('t');
    const ode = new Mul(new Num(-0.5), y);
    const sol4 = cas.evaluate(new Call('rk4', [ode, y, t, new Num(1), new Num(0), new Num(2), new Num(0.01)]));
    const last4 = sol4.elements[sol4.elements.length - 1];
    const lastT = last4.elements[0].value;
    const lastY = last4.elements[1].value;
    const exact = Math.exp(-1); // y(2) = e^(-1)
    assertApprox(lastY, exact, 'rk4 at t=2 (y≈e^-1≈0.1353)', 1e-4);
    assertApprox(lastT, 2.0,    'rk4 final t', 1e-8);
    console.log(`  Points: ${sol4.elements.length} (expected ~200)`);
}

{
    // RK45 adaptive
    const y = new Sym('y');
    const t = new Sym('t');
    const ode = new Mul(new Num(-0.5), y);
    const sol45 = cas.evaluate(new Call('rk45', [ode, y, t, new Num(1), new Num(0), new Num(2), new Num(0.1), new Num(1e-8), new Vec([])]));
    const meta  = sol45.elements[sol45.elements.length - 1];
    const last  = sol45.elements[sol45.elements.length - 2];
    const lastT = last.elements[0].value;
    const lastY = last.elements[1].value;
    const lastE = last.elements[2].value;
    const exact = Math.exp(-1);
    assertApprox(lastY, exact, 'rk45 at t=2 (y≈e^-1≈0.1353)', 1e-4);
    assertApprox(lastT, 2.0,   'rk45 final t', 1e-8);
    console.log(`  Steps=${meta.elements[1].value}, Rejected=${meta.elements[2].value}, Evals=${meta.elements[3].value}`);
    console.log(`  Last error flag: ${lastE.toFixed(8)} (0 = accepted step)`);
}

// ── 2. Harmonic oscillator: y'' = -y  →  y1'=y2, y2'=-y1, y1(0)=0, y2(0)=1 ──
console.log('\n--- Problem 2: Harmonic oscillator (y\'\'=-y, y(0)=0, y\'(0)=1) ---');

{
    const y1 = new Sym('y1');
    const y2 = new Sym('y2');
    const t  = new Sym('t');
    const odeVec = new Vec([y2, new Mul(new Num(-1), y1)]);
    const depVec = new Vec([y1, y2]);
    const initV  = new Vec([new Num(0), new Num(1)]);
    const sol4 = cas.evaluate(new Call('rk4', [odeVec, depVec, t, initV, new Num(0), new Num(10), new Num(0.001)]));
    const last = sol4.elements[sol4.elements.length - 1];
    const exp_y1 = Math.sin(10), exp_y2 = Math.cos(10);
    assertApprox(last.elements[1].value, exp_y1, 'rk4 y1(10)≈sin(10)', 1e-3);
    assertApprox(last.elements[2].value, exp_y2, 'rk4 y2(10)≈cos(10)', 1e-3);
    console.log(`  Points: ${sol4.elements.length}`);
}

{
    const y1 = new Sym('y1');
    const y2 = new Sym('y2');
    const t  = new Sym('t');
    const odeVec = new Vec([y2, new Mul(new Num(-1), y1)]);
    const depVec = new Vec([y1, y2]);
    const initV  = new Vec([new Num(0), new Num(1)]);
    const sol45 = cas.evaluate(new Call('rk45', [odeVec, depVec, t, initV, new Num(0), new Num(10), new Num(0.1), new Num(1e-6), new Vec([])]));
    const meta  = sol45.elements[sol45.elements.length - 1];
    const last  = sol45.elements[sol45.elements.length - 2];
    const exp_y1 = Math.sin(10), exp_y2 = Math.cos(10);
    assertApprox(last.elements[1].value, exp_y1, 'rk45 y1(10)', 1e-3);
    assertApprox(last.elements[2].value, exp_y2, 'rk45 y2(10)', 1e-3);
    console.log(`  Steps=${meta.elements[1].value}, Rejected=${meta.elements[2].value}`);
}

// ── 3. Damped oscillation: y'' = -0.5*y' - y, y(0)=1, y'(0)=0 ───────────────
console.log('\n--- Problem 3: Damped oscillation (y\'\' = -0.5*y\' - y) ---');

{
    const y1 = new Sym('y1');
    const y2 = new Sym('y2');
    const t  = new Sym('t');
    const odeVec = new Vec([
        new Add(new Mul(new Num(-0.5), y2), new Mul(new Num(-1), y1)),
        y2
    ]);
    const depVec = new Vec([y1, y2]);
    const initV  = new Vec([new Num(1), new Num(0)]);
    const sol4 = cas.evaluate(new Call('rk4', [odeVec, depVec, t, initV, new Num(0), new Num(10), new Num(0.005)]));
    const last = sol4.elements[sol4.elements.length - 1];
    const y1_val = last.elements[1].value;
    const y2_val = last.elements[2].value;
    const ampEst = Math.sqrt(y1_val * y1_val + y2_val * y2_val);
    console.log(`  y1(10)≈${y1_val.toFixed(6)}, y2(10)≈${y2_val.toFixed(6)}, Amp≈${ampEst.toFixed(6)}`);
    assertTrue(ampEst < 0.5, 'Amplitude at t=10 is small (damped), as expected');
}

{
    const y1 = new Sym('y1');
    const y2 = new Sym('y2');
    const t  = new Sym('t');
    const odeVec = new Vec([
        new Add(new Mul(new Num(-0.5), y2), new Mul(new Num(-1), y1)),
        y2
    ]);
    const depVec = new Vec([y1, y2]);
    const initV  = new Vec([new Num(1), new Num(0)]);
    const sol45 = cas.evaluate(new Call('rk45', [odeVec, depVec, t, initV, new Num(0), new Num(10), new Num(0.1), new Num(1e-10), new Vec([])]));
    const meta  = sol45.elements[sol45.elements.length - 1];
    const last  = sol45.elements[sol45.elements.length - 2];
    const y1_val = last.elements[1].value;
    const ampEst = Math.sqrt(y1_val * y1_val + last.elements[2].value * last.elements[2].value);
    console.log(`  rk45: y1(10)≈${y1_val.toFixed(6)}, Steps=${meta.elements[1].value}, Rejected=${meta.elements[2].value}`);
    assertTrue(ampEst < 0.5, 'Damped oscillation amplitude small (rk45)');
}

// ── 4. SIR model (epidemiology) ─────────────────────────────────────────────
// dS/dt = -β*S*I,  dI/dt = β*S*I - γ*I,  dR/dt = γ*I
// Parameters: β=0.3, γ=0.1,  S(0)=0.99, I(0)=0.01, R(0)=0
console.log('\n--- Problem 4: SIR model (S\'=-β*S*I, I\'=β*S*I-γ*I, R\'=γ*I) ---');

{
    const S = new Sym('S');
    const I = new Sym('I');
    const R = new Sym('R');
    const t = new Sym('t');
    const beta = new Num(0.3), gamma = new Num(0.1);
    const sir = new Vec([
        new Mul(new Num(-1), new Mul(beta, new Mul(S, I))),
        new Sub(new Mul(beta, new Mul(S, I)), new Mul(gamma, I)),
        new Mul(gamma, I)
    ]);
    const dep = new Vec([S, I, R]);
    const ini = new Vec([new Num(0.99), new Num(0.01), new Num(0)]);
    const sol = cas.evaluate(new Call('rk4', [sir, dep, t, ini, new Num(0), new Num(50), new Num(0.05)]));
    const last = sol.elements[sol.elements.length - 1];
    const Sval = last.elements[1].value;
    const Ival = last.elements[2].value;
    const Rval = last.elements[3].value;
    console.log(`  S(50)≈${Sval.toFixed(6)}, I(50)≈${Ival.toFixed(6)}, R(50)≈${Rval.toFixed(6)}`);
    const total = Sval + Ival + Rval;
    assertApprox(total, 1.0, 'S+I+R≈1 (conserved)', 1e-3);
    assertTrue(Ival < 0.1, 'I(50) < 0.1 (epidemic subsiding)');
}

{
    const S = new Sym('S');
    const I = new Sym('I');
    const R = new Sym('R');
    const t = new Sym('t');
    const beta = new Num(0.3), gamma = new Num(0.1);
    const sir = new Vec([
        new Mul(new Num(-1), new Mul(beta, new Mul(S, I))),
        new Sub(new Mul(beta, new Mul(S, I)), new Mul(gamma, I)),
        new Mul(gamma, I)
    ]);
    const dep = new Vec([S, I, R]);
    const ini = new Vec([new Num(0.99), new Num(0.01), new Num(0)]);
    const sol45 = cas.evaluate(new Call('rk45', [sir, dep, t, ini, new Num(0), new Num(50), new Num(0.1), new Num(1e-6), new Vec([])]));
    const meta  = sol45.elements[sol45.elements.length - 1];
    const last  = sol45.elements[sol45.elements.length - 2];
    const Ival = last.elements[2].value;
    console.log(`  rk45: I(50)≈${Ival.toFixed(6)}, Steps=${meta.elements[1].value}, Rejected=${meta.elements[2].value}`);
    assertTrue(Ival < 0.1, 'SIR epidemic subsiding (rk45)');
}

// ── 5. Error estimation check: rk45 should report error=0 for accepted steps ─
console.log('\n--- Problem 5: rk45 error estimation check ---');
{
    const y = new Sym('y');
    const t = new Sym('t');
    const ode = new Mul(new Num(-0.5), y);
    const sol45 = cas.evaluate(new Call('rk45', [ode, y, t, new Num(1), new Num(0), new Num(1), new Num(0.1), new Num(1e-6), new Vec([])]));
    // All non-meta rows should have err=0 (accepted steps)
    let allZeroErr = true;
    for (let i = 0; i < sol45.elements.length - 1; i++) {
        const err = sol45.elements[i].elements[2].value;
        if (err !== 0) { allZeroErr = false; break; }
    }
    assertTrue(allZeroErr, 'All accepted steps have error flag = 0');
    const meta = sol45.elements[sol45.elements.length - 1];
    console.log(`  Steps=${meta.elements[1].value}, Rejected=${meta.elements[2].value}, Evals=${meta.elements[3].value}`);
}

// ── Summary ──────────────────────────────────────────────────────────────────
// ── 6. odeplot check: ensure odeplot returns correct structure for list plotting ─
console.log('\n--- Problem 6: odeplot test ---');
{
    const y = new Sym('y');
    const t = new Sym('t');
    const ode = new Mul(new Num(-0.5), y);
    // Get sol from rk45
    const sol45 = cas.evaluate(new Call('rk45', [ode, y, t, new Num(1), new Num(0), new Num(1), new Num(0.1), new Num(1e-6), new Vec([])]));
    const plotObj = cas.evaluate(new Call('odeplot', [sol45]));

    assertTrue(plotObj.type === 'plot', 'odeplot returns a plot type');
    // It might return multiple plots, or single
    const pts = plotObj.scatter || (plotObj.plots && plotObj.plots[0].scatter);
    assertTrue(pts && pts.length > 0, 'odeplot has scatter points');

    let validPoints = true;
    for(let pt of pts) {
        if(typeof pt.x !== 'number' || typeof pt.y !== 'number') validPoints = false;
    }
    assertTrue(validPoints, 'odeplot scatter points are valid numeric coordinates');
}

console.log(`\n=== Results: ${passed} passed, ${failed} failed ===`);
if (failed > 0) process.exit(1);
