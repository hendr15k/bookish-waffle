
const cas = new CAS();

// 1. Limit check
const lim1 = cas.evaluate(new Call('limit', [
    new Pow(new Add(new Num(1), new Div(new Num(1), new Sym('x'))), new Sym('x')),
    new Sym('x'),
    new Sym('Infinity')
]));
console.log("limit((1+1/x)^x, x, Infinity) =", lim1.toString());

// 2. Integration check
const int1 = cas.evaluate(new Call('integrate', [
    new Mul(new Pow(new Sym('x'), new Num(2)), new Call('exp', [new Sym('x')])),
    new Sym('x')
]));
console.log("integrate(x^2 * exp(x), x) =", int1.toString());

// 3. Partitions check
try {
    console.log("partitions(5) =", cas.evaluate(new Call('partitions', [new Num(5)])).toString());
} catch(e) {
    console.log("partitions(5) error:", e.message);
}
