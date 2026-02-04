
const { CAS } = require('../js/cas');
const { Expr, Num, Sym, Call, Vec } = require('../js/expression');

const cas = new CAS();

console.log("=== Testing Graph Features ===");
// Adjacency Matrix for a triangle graph with weights: 0-1 (1), 1-2 (2), 0-2 (4)
const adj = new Vec([
    new Vec([new Num(0), new Num(1), new Num(4)]),
    new Vec([new Num(1), new Num(0), new Num(2)]),
    new Vec([new Num(4), new Num(2), new Num(0)])
]);

try {
    const path = cas.evaluate(new Call('shortestPath', [adj, new Num(0), new Num(2)]));
    console.log("Shortest Path 0->2:", path.toString());

    if (path.toString() === "[0, 1, 2]") {
        console.log("PASS: Shortest Path");
    } else {
        console.log("FAIL: Shortest Path");
        throw new Error("Shortest Path Failed");
    }

    const mst = cas.evaluate(new Call('mst', [adj]));
    console.log("MST Matrix:", mst.toString());
    const row0 = mst.elements[0].toString();
    if (row0.includes("0, 1, 0") || row0.includes("0, 1, 0")) {
         console.log("PASS: MST");
    } else {
         console.log("FAIL: MST");
         throw new Error("MST Failed");
    }

} catch (e) {
    console.log("Graph Error:", e.stack);
    process.exit(1);
}

console.log("\n=== Testing Statistics Features (ANOVA) ===");
const groups = new Vec([
    new Vec([new Num(1), new Num(2), new Num(3)]),
    new Vec([new Num(4), new Num(5), new Num(6)]),
    new Vec([new Num(7), new Num(8), new Num(9)])
]);

try {
    const res = cas.evaluate(new Call('anova', [groups]));
    console.log("ANOVA Result [F, p]:", res.toString());
    const F = res.elements[0].value;
    const p = res.elements[1].value;

    if (F > 10 && p < 0.05) {
        console.log("PASS: ANOVA significant");
    } else {
        console.log("FAIL: ANOVA insignificant?");
        throw new Error("ANOVA Failed");
    }
} catch (e) {
    console.log("Stats Error:", e.stack);
    process.exit(1);
}

console.log("\n=== Testing Optimization (Knapsack) ===");
const w = new Vec([new Num(10), new Num(20), new Num(30)]);
const v = new Vec([new Num(60), new Num(100), new Num(120)]);
const cap = new Num(50);

try {
    const res = cas.evaluate(new Call('knapsack', [w, v, cap]));
    console.log("Knapsack Result [Max, Items]:", res.toString());
    const maxVal = res.elements[0].value;
    const items = res.elements[1].toString();

    if (maxVal === 220 && (items.includes("1, 2") || items.includes("2, 1"))) {
        console.log("PASS: Knapsack");
    } else {
        console.log("FAIL: Knapsack");
        throw new Error("Knapsack Failed");
    }
} catch (e) {
    console.log("Optimization Error:", e.stack);
    process.exit(1);
}

console.log("\nAll new features passed.");
