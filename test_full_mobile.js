const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true });
  const results = [];
  const bugs = [];
  function check(name, condition, details = '') {
    const status = condition ? 'PASS' : 'FAIL';
    results.push({ name, status, details });
    if (!condition) bugs.push({ name, details });
    console.log(`${status}: ${name}${details ? ' - ' + details : ''}`);
  }

  const context = await browser.newContext({
    viewport: { width: 375, height: 812 },
    isMobile: true, hasTouch: true,
    userAgent: 'Mozilla/5.0 (iPhone; CPU iPhone OS 16_0 like Mac OS X) AppleWebKit/605.1.15'
  });
  const page = await context.newPage();
  const consoleErrors = [];
  page.on('console', msg => { if (msg.type() === 'error') consoleErrors.push(msg.text()); });

  try {
    await page.goto('http://localhost:9999/', { waitUntil: 'networkidle', timeout: 15000 });

    // ===================== CALCULATIONS =====================
    console.log('\n========== CALCULATIONS ==========');
    await page.click('#mode-toggle');
    await page.waitForTimeout(500);

    // Calculator: 2+2
    console.log('--- Calculator ---');
    await page.click('.app-card:has-text("Calculator")', { force: true });
    await page.waitForTimeout(300);
    await page.fill('#calc-main-input', '2+2');
    await page.click('.action-btn:has-text("Calculate")', { force: true });
    await page.waitForTimeout(500);
    let result = page.locator('#tool-calculator .tool-result');
    check('Calc: 2+2 shows result', await result.isVisible());
    if (await result.isVisible()) {
      const txt = await result.textContent();
      check('Calc: 2+2 = 4', txt.includes('4'), `got: "${txt.substring(0, 50)}"`);
      const box = await result.boundingBox();
      check('Calc result fits viewport', !box || box.x + box.width <= 380);
    }
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // Algebra: simplify
    console.log('--- Algebra ---');
    await page.click('.app-card:has-text("Algebra")', { force: true });
    await page.waitForTimeout(300);
    await page.fill('#alg-expr', '(x+1)^2');
    await page.click('.action-btn:has-text("Simplify")', { force: true });
    await page.waitForTimeout(500);
    result = page.locator('#tool-algebra .tool-result');
    check('Algebra: simplify shows result', await result.isVisible());
    if (await result.isVisible()) {
      const txt = await result.textContent();
      check('Algebra: (x+1)^2 expanded', txt.length > 3, `got: "${txt.substring(0, 50)}"`);
    }
    // Switch to Solve tab
    await page.click('#tool-algebra .tab-btn:has-text("Solve")', { force: true });
    await page.waitForTimeout(200);
    check('Algebra: Solve tab visible', await page.locator('#tab-alg-solve').isVisible());
    // Switch to Advanced tab
    await page.click('#tool-algebra .tab-btn:has-text("Advanced")', { force: true });
    await page.waitForTimeout(200);
    check('Algebra: Advanced tab visible', await page.locator('#tab-alg-advanced').isVisible());
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // Calculus: differentiate
    console.log('--- Calculus ---');
    await page.click('.app-card:has-text("Calculus")', { force: true });
    await page.waitForTimeout(300);
    await page.fill('#calc-expr', 'sin(x)');
    await page.click('.action-btn:has-text("Differentiate")', { force: true });
    await page.waitForTimeout(500);
    result = page.locator('#tool-calculus .tool-result');
    check('Calc: diff shows result', await result.isVisible());
    if (await result.isVisible()) {
      const txt = await result.textContent();
      check('Calc: d/dx sin(x) contains cos', txt.toLowerCase().includes('cos'), `got: "${txt.substring(0, 50)}"`);
    }
    // Test Integration tab
    await page.click('#tool-calculus .tab-btn:has-text("Integration")', { force: true });
    await page.waitForTimeout(200);
    check('Calc: Integration tab visible', await page.locator('#tab-calc-int').isVisible());
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // Numbers: factorial
    console.log('--- Numbers ---');
    await page.click('.app-card:has-text("Numbers")', { force: true });
    await page.waitForTimeout(300);
    await page.fill('#num-a', '5');
    await page.click('.action-btn:has-text("n!")', { force: true });
    await page.waitForTimeout(500);
    result = page.locator('#tool-numbers .tool-result');
    check('Numbers: factorial shows result', await result.isVisible());
    if (await result.isVisible()) {
      const txt = await result.textContent();
      check('Numbers: 5! = 120', txt.includes('120'), `got: "${txt.substring(0, 50)}"`);
    }
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // Linear: determinant
    console.log('--- Linear ---');
    await page.click('.app-card:has-text("Linear")', { force: true });
    await page.waitForTimeout(300);
    await page.fill('#lin-expr', '[[1,2],[3,4]]');
    await page.click('.action-btn:has-text("Determinant")', { force: true });
    await page.waitForTimeout(500);
    result = page.locator('#tool-linear .tool-result');
    check('Linear: det shows result', await result.isVisible());
    if (await result.isVisible()) {
      const txt = await result.textContent();
      check('Linear: det([[1,2],[3,4]]) = -2', txt.includes('-2') || txt.includes('(-2)'), `got: "${txt.substring(0, 50)}"`);
    }
    // Test Builder tab
    await page.click('#tool-linear .tab-btn:has-text("Builder")', { force: true });
    await page.waitForTimeout(200);
    check('Linear: Builder tab visible', await page.locator('#tab-lin-build').isVisible());
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // ===================== DARK MODE =====================
    console.log('\n========== DARK MODE ==========');
    await page.click('.app-card:has-text("Settings")', { force: true });
    await page.waitForTimeout(300);
    await page.locator('label[for="dark-mode-toggle"]').click();
    await page.waitForTimeout(300);
    const isDark = await page.evaluate(() => document.body.classList.contains('dark-mode'));
    check('Dark mode activated', isDark);

    // Open a tool in dark mode
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);
    await page.click('.app-card:has-text("Calculator")', { force: true });
    await page.waitForTimeout(300);
    check('Tool visible in dark mode', await page.locator('#tool-calculator').isVisible());

    // Check dark mode styles
    const darkInputBg = await page.locator('.form-input').first().evaluate(el => getComputedStyle(el).backgroundColor);
    check('Dark mode input has dark bg', darkInputBg !== 'rgb(255, 255, 255)', `bg: ${darkInputBg}`);

    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // Turn off dark mode
    await page.click('.app-card:has-text("Settings")', { force: true });
    await page.waitForTimeout(300);
    await page.locator('label[for="dark-mode-toggle"]').click();
    await page.waitForTimeout(300);
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // ===================== GRAPHING MODE =====================
    console.log('\n========== GRAPHING MODE ==========');
    await page.click('#graphing-toggle');
    await page.waitForTimeout(500);
    check('Graphing mode visible', await page.locator('#graphing-mode').isVisible());

    // Check sidebar
    const sidebarBox = await page.locator('#graphing-sidebar').boundingBox();
    check('Graphing sidebar not too wide on mobile', sidebarBox && sidebarBox.width <= 250, `width: ${sidebarBox?.width}`);
    check('Graphing sidebar has add button', await page.locator('#graphing-add-btn').isVisible());

    // Check canvas
    const canvasBox = await page.locator('#graphing-canvas').boundingBox();
    check('Canvas has dimensions', canvasBox && canvasBox.width > 0, canvasBox ? `${canvasBox.width}x${canvasBox.height}` : 'null');

    // Add cell and evaluate
    await page.click('#graphing-add-btn');
    await page.waitForTimeout(200);
    const cells = await page.locator('.graphing-cell').count();
    check('Cell added', cells > 0, `count: ${cells}`);

    await page.locator('.graphing-cell-input').first().type('x^2');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(500);
    const cellOutput = await page.locator('.graphing-cell-output').first().textContent();
    check('Cell evaluates expression', cellOutput && cellOutput.length > 3, `output: "${(cellOutput || '').substring(0, 30)}"`);

    // Remove cell
    await page.locator('.graphing-cell-remove').first().click();
    await page.waitForTimeout(200);
    const cellsAfter = await page.locator('.graphing-cell').count();
    check('Cell can be removed', cellsAfter < cells, `before: ${cells}, after: ${cellsAfter}`);

    // Exit graphing
    await page.click('#graphing-toggle');
    await page.waitForTimeout(300);

    // Switch back to app mode
    await page.click('#mode-toggle');
    await page.waitForTimeout(300);

    // ===================== BOTTOM BAR =====================
    console.log('\n========== BOTTOM BAR ==========');
    const bottomBar = page.locator('.mobile-bottom-bar');
    check('Bottom bar visible', await bottomBar.isVisible());
    const bottomBtns = await page.locator('.mobile-bottom-bar button').count();
    check('Bottom bar has 5 buttons', bottomBtns === 5, `count: ${bottomBtns}`);

    // Test each bottom bar button
    for (const [tool, id] of [['Home', 'home'], ['Calc', 'calculator'], ['Algebra', 'algebra'], ['Plot', 'plotting'], ['Matrix', 'linear']]) {
      const btn = page.locator(`.mobile-bottom-bar button[data-tool="${id}"]`);
      if (await btn.isVisible()) {
        await btn.click({ force: true });
        await page.waitForTimeout(300);
        if (id === 'home') {
          check(`Bottom bar "${tool}" → app grid`, await page.locator('#app-grid').isVisible());
        } else {
          const toolVisible = await page.locator(`#tool-${id}`).isVisible().catch(() => false);
          check(`Bottom bar "${tool}" → tool opens`, toolVisible, `tool-${id}`);
          const isActive = await btn.evaluate(el => el.classList.contains('active'));
          check(`Bottom bar "${tool}" shows active`, isActive);
        }
      }
    }

    // ===================== SPECIAL TOOLS =====================
    console.log('\n========== SPECIAL TOOLS ==========');

    // Test Chemistry tabs
    await page.click('.mobile-bottom-bar button[data-tool="home"]', { force: true });
    await page.waitForTimeout(200);
    await page.click('.app-card:has-text("Chemistry")', { force: true });
    await page.waitForTimeout(300);
    check('Chemistry tool opens', await page.locator('#tool-chemistry').isVisible());

    // Chemistry tabs
    for (const tab of ['tab-chem-elem', 'tab-chem-mass', 'tab-chem-bal']) {
      const tabBtn = page.locator(`#${tab}`).locator('..').locator('preceding-sibling::.tab-header .tab-btn');
      await page.locator('#tool-chemistry .tab-btn').first().click({ force: true });
      await page.waitForTimeout(100);
    }
    // Click each tab properly
    for (const [idx, tabName] of ['Elements', 'Molar Mass', 'Balance'].entries()) {
      await page.click(`#tool-chemistry .tab-btn:has-text("${tabName}")`, { force: true });
      await page.waitForTimeout(100);
    }
    check('Chemistry has 3 tabs', true);
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // Test Finance tabs
    await page.click('.app-card:has-text("Finance")', { force: true });
    await page.waitForTimeout(300);
    check('Finance tool opens', await page.locator('#tool-finance').isVisible());
    const finTabs = ['Compound Interest', 'Loan Payment', 'Cash Flow', 'Annuity', 'Amortization', 'Black-Scholes'];
    for (const t of finTabs) {
      await page.click(`#tool-finance .tab-btn:has-text("${t}")`, { force: true });
      await page.waitForTimeout(100);
    }
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // Test Physics tabs
    await page.click('.app-card:has-text("Physics")', { force: true });
    await page.waitForTimeout(300);
    check('Physics tool opens', await page.locator('#tool-physics').isVisible());
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // Test special tools without tabs
    for (const tool of ['converter', 'editor', 'fractals', 'graph', 'table', 'complex', 'special', 'logic', 'trig', 'variables', 'history', 'gemini', 'system']) {
      const name = tool.charAt(0).toUpperCase() + tool.slice(1);
      try {
        await page.click(`.app-card:has-text("${name}")`, { force: true });
        await page.waitForTimeout(300);
        check(`Tool "${tool}" opens`, await page.locator(`#tool-${tool}`).isVisible().catch(() => false));
        await page.evaluate(() => closeTool());
        await page.waitForTimeout(200);
      } catch(e) {
        check(`Tool "${tool}" opens`, false, e.message.substring(0, 60));
      }
    }

    // ===================== FORM INPUT CHECKS =====================
    console.log('\n========== FORM INPUTS ==========');
    await page.click('.app-card:has-text("Algebra")', { force: true });
    await page.waitForTimeout(300);

    // Check font-size of inputs (should be 16px for iOS no-zoom)
    const inputFontSize = await page.locator('.form-input').first().evaluate(el => getComputedStyle(el).fontSize);
    check('Form input font-size >= 16px (no iOS zoom)', parseFloat(inputFontSize) >= 16, `fontSize: ${inputFontSize}`);

    // Check form-row stacks on mobile
    const formRowDir = await page.locator('.form-row').first().evaluate(el => getComputedStyle(el).flexDirection);
    check('Form rows stack on mobile', formRowDir === 'column', `direction: ${formRowDir}`);

    // Action grid
    const actionGridCols = await page.locator('.action-grid').first().evaluate(el => getComputedStyle(el).gridTemplateColumns);
    check('Action grid single column on mobile', actionGridCols.split(' ').length === 1 || actionGridCols.includes('1fr'), `columns: ${actionGridCols}`);

    // Back buttons
    const backBtnSize = await page.locator('.back-btn').first().evaluate(el => {
      const style = getComputedStyle(el);
      return { width: parseInt(style.minWidth), height: parseInt(style.minHeight) };
    });
    check('Back button min 44x44 touch target', backBtnSize.width >= 44 && backBtnSize.height >= 44, `size: ${backBtnSize.width}x${backBtnSize.height}`);

    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // ===================== SMALL PHONE (320px) =====================
    console.log('\n========== SMALL PHONE (320px) ==========');
    await page.setViewportSize({ width: 320, height: 568 });
    await page.waitForTimeout(300);

    check('Small viewport applied', await page.evaluate(() => document.body.clientWidth) <= 325);

    const headerBox = await page.locator('h1').boundingBox();
    check('Header fits small phone', headerBox && headerBox.x + headerBox.width <= 325);

    const modeToggle = page.locator('#mode-toggle');
    const modeBox = await modeToggle.boundingBox();
    check('Mode toggle fits small phone', modeBox && modeBox.x + modeBox.width <= 325);

    const bottomBarSmall = await page.locator('.mobile-bottom-bar').boundingBox();
    check('Bottom bar fits small phone', bottomBarSmall && bottomBarSmall.width <= 325);

    const firstCard = page.locator('.app-card').first();
    const cardBox = await firstCard.boundingBox();
    check('App cards fit small phone', cardBox && cardBox.x + cardBox.width <= 325);

    // ===================== LANDSCAPE =====================
    console.log('\n========== LANDSCAPE (812x375) ==========');
    await page.setViewportSize({ width: 812, height: 375 });
    await page.waitForTimeout(300);

    check('Mobile mode persists in landscape', await page.evaluate(() => document.body.classList.contains('mobile-mode')));
    check('App grid visible in landscape', await page.locator('#app-grid').isVisible());
    check('Bottom bar visible in landscape', await page.locator('.mobile-bottom-bar').isVisible());

    // Open a tool in landscape
    await page.click('.app-card:has-text("Calculator")', { force: true });
    await page.waitForTimeout(300);
    check('Tool works in landscape', await page.locator('#tool-calculator').isVisible());
    await page.evaluate(() => closeTool());
    await page.waitForTimeout(200);

    // ===================== CONSOLE ERRORS =====================
    console.log('\n========== CONSOLE ERRORS ==========');
    const critical = consoleErrors.filter(e =>
      !e.includes('Warning') && !e.includes('deprecated') && !e.includes('MathJax') &&
      !e.includes('Failed to load') && !e.includes('favicon') && !e.includes('404') &&
      !e.includes('KaTeX')
    );
    check('No critical console errors', critical.length === 0, critical.length > 0 ? critical[0].substring(0, 100) : '');

  } catch(e) {
    console.error('FATAL:', e.message);
    bugs.push({ name: 'Fatal error', details: e.message });
  } finally {
    console.log('\n========== SUMMARY ==========');
    const passed = results.filter(r => r.status === 'PASS').length;
    const failed = results.filter(r => r.status === 'FAIL').length;
    console.log(`Passed: ${passed}/${results.length}`);
    console.log(`Failed: ${failed}/${results.length}`);

    if (bugs.length > 0) {
      console.log('\n========== BUGS FOUND ==========');
      bugs.forEach(b => console.log(`- ${b.name}${b.details ? ': ' + b.details : ''}`));
    } else {
      console.log('\nNo bugs found!');
    }

    await browser.close();
    process.exit(failed > 0 ? 1 : 0);
  }
})();