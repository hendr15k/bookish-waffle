const { JSDOM } = require('jsdom');
const fs = require('fs');
const path = require('path');

const html = fs.readFileSync(path.join(__dirname, '..', 'index.html'), 'utf8');
const dom = new JSDOM(html, { runScripts: 'dangerously' });
global.window = dom.window;
global.document = dom.window.document;
global.navigator = dom.window.navigator;

const codeHelp = fs.readFileSync(path.join(__dirname, '..', 'js', 'help.js'), 'utf8');
const codeCas = fs.readFileSync(path.join(__dirname, '..', 'js', 'cas.js'), 'utf8');
dom.window.eval(codeHelp);
dom.window.eval(codeCas);

const cas = new dom.window.CAS();

function assertEqual(actual, expected, label) {
  if (actual !== expected) {
    throw new Error(`${label} failed: expected ${expected}, got ${actual}`);
  }
  console.log(`[PASS] ${label}`);
}

assertEqual(cas.evaluate('egcd(240,46)').toString(), cas.evaluate('xgcd(240,46)').toString(), 'egcd alias');
assertEqual(cas.evaluate('invmod(3,11)').toString(), cas.evaluate('modInverse(3,11)').toString(), 'invmod alias');
assertEqual(cas.evaluate('inverse(3,11)').toString(), cas.evaluate('modInverse(3,11)').toString(), 'inverse alias');
assertEqual(cas.evaluate('powmod(2,10,1000)').toString(), cas.evaluate('modPow(2,10,1000)').toString(), 'powmod alias');

const help = dom.window.HELP_DATA;
for (const key of ['egcd', 'invmod', 'inverse', 'powmod']) {
  if (!help[key]) throw new Error(`Missing help entry for ${key}`);
  console.log(`[PASS] help entry for ${key}`);
}

console.log('All alias tests passed');
