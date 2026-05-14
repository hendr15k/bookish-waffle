const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true });
  const results = [];
  const bugs = [];
  function check(name, condition, details = '') {
    const status = condition ? 'PASS' : 'FAIL';
    results.push({ name, status, details });
    if (!condition) {
      bugs.push({ name, details });
    }
    console.log(`${status}: ${name}${details ? ' - ' + details : ''}`);
  }

  console.log('\n=== Comprehensive Mobile UI Bug Testing ===\n');

  const context = await browser.newContext({
    viewport: { width: 375, height: 812 },
    isMobile: true,
    hasTouch: true,
    userAgent: 'Mozilla/5.0 (iPhone; CPU iPhone OS 16_0 like Mac OS X) AppleWebKit/605.1.15'
  });
  const page = await context.newPage();
  const consoleErrors = [];
  page.on('console', msg => {
    if (msg.type() === 'error') consoleErrors.push(msg.text());
  });

  try {
    await page.goto('http://localhost:9999/', { waitUntil: 'networkidle', timeout: 15000 });

    // 1. App Mode Basic Tests
    console.log('--- App Mode Basics ---');
    await page.click('#mode-toggle');
    await page.waitForTimeout(500);

    const appGrid = page.locator('#app-grid');
    const gridVisible = await appGrid.isVisible();
    check('App grid visible', gridVisible);

    const cardCount = await page.locator('.app-card').count();
    check('App has multiple cards', cardCount >= 20, `count: ${cardCount}`);

    // 2. App Categories
    console.log('\n--- App Categories ---');
    const catBtns = await page.locator('#app-categories button').count();
    check('Has category filters', catBtns >= 3, `count: ${catBtns}`);

    await page.click('#app-categories button:has-text("Math")');
    await page.waitForTimeout(300);
    const mathVisible = await page.evaluate(() => {
      return Array.from(document.querySelectorAll('.app-card'))
        .filter(c => c.style.display !== 'none').length;
    });
    check('Math filter works', mathVisible < cardCount, `visible: ${mathVisible}`);
    await page.click('#app-categories button:has-text("All")');
    await page.waitForTimeout(200);

    // 3. App Search
    console.log('\n--- App Search ---');
    await page.fill('#app-search', 'calc');
    await page.waitForTimeout(300);
    const searchVisible = await page.evaluate(() => {
      return Array.from(document.querySelectorAll('.app-card'))
        .filter(c => c.style.display !== 'none').length;
    });
    check('Search filters apps', searchVisible < cardCount, `visible: ${searchVisible}`);
    await page.fill('#app-search', '');
    await page.waitForTimeout(200);

    // 4. Tool Navigation
    console.log('\n--- Tool Navigation ---');
    const tools = ['Calculator', 'Algebra', 'Calculus', 'Plotting', 'Linear', 'Statistics'];
    for (const tool of tools) {
      await page.click(`.app-card:has-text("${tool}")`);
      await page.waitForTimeout(300);
      const toolView = page.locator(`#tool-${tool.toLowerCase()}`);
      const visible = await toolView.isVisible();
      check(`${tool} tool opens`, visible, visible ? '' : 'FAILED');

      const backBtn = toolView.locator('.back-btn');
      if (await backBtn.isVisible()) {
        await backBtn.click();
        await page.waitForTimeout(200);
      }
    }

    // 5. Algebra Tabs
    console.log('\n--- Algebra Tabs ---');
    await page.click('.app-card:has-text("Algebra")');
    await page.waitForTimeout(300);

    const algTabs = ['Simplify & Expand', 'Solve', 'Numeric Solve', 'Advanced'];
    for (const tab of algTabs) {
      await page.click(`#tool-algebra .tab-btn:has-text("${tab}")`);
      await page.waitForTimeout(200);
      let tabId;
      if (tab === 'Simplify & Expand') tabId = 'tab-alg-basic';
      else if (tab === 'Numeric Solve') tabId = 'tab-alg-numsolve';
      else tabId = 'tab-alg-' + tab.toLowerCase().replace(/ & /g, '-').replace(/ /g, '-');
      check(`${tab} tab works`, await page.locator(`#${tabId}`).isVisible());
    }

    // 6. Calculus Tabs
    console.log('\n--- Calculus Tabs ---');
    await page.locator('#tool-algebra .back-btn').click();
    await page.waitForTimeout(200);
    await page.click('.app-card:has-text("Calculus")');
    await page.waitForTimeout(300);

    const tabHeader = page.locator('#tool-calculus .tab-header');
    const scrolls = await tabHeader.evaluate(el => el.scrollWidth > el.clientWidth);
    check('Calculus tabs scroll horizontally', scrolls);

    const calcTabTests = [
      { btn: 'Diff & Limit', id: 'tab-calc-diff' },
      { btn: 'Integration', id: 'tab-calc-int' },
      { btn: 'Diff Eq', id: 'tab-calc-desolve' },
      { btn: 'Taylor Series', id: 'tab-calc-series' },
      { btn: 'Riemann', id: 'tab-calc-riemann' },
      { btn: 'Fourier', id: 'tab-calc-fourier' },
      { btn: 'Tangent', id: 'tab-calc-tangent' },
      { btn: 'Vector Calc', id: 'tab-calc-vector' },
      { btn: 'Laplace', id: 'tab-calc-laplace' },
      { btn: 'Curve Tools', id: 'tab-calc-geom' },
      { btn: 'ODE Num.', id: 'tab-calc-ode-num' },
      { btn: 'Transforms', id: 'tab-calc-transforms' },
      { btn: 'Optimization', id: 'tab-calc-optimize' },
      { btn: 'Analysis', id: 'tab-calc-analyze' }
    ];

    for (const tab of calcTabTests) {
      await page.click(`#tool-calculus .tab-btn:has-text("${tab.btn}")`);
      await page.waitForTimeout(100);
      check(`${tab.btn} tab works`, await page.locator(`#${tab.id}`).isVisible());
    }

    // 7. Graphing Mode
    console.log('\n--- Graphing Mode ---');
    await page.locator('#tool-calculus .back-btn').click();
    await page.waitForTimeout(200);
    await page.click('#graphing-toggle');
    await page.waitForTimeout(500);

    check('Graphing mode visible', await page.locator('#graphing-mode').isVisible());

    const canvas = page.locator('#graphing-canvas');
    const canvasBox = await canvas.boundingBox();
    check('Graphing canvas has dimensions', canvasBox && canvasBox.width > 0, canvasBox ? `${canvasBox.width}x${canvasBox.height}` : 'null');

    await page.click('#graphing-add-btn');
    await page.waitForTimeout(200);

    const cellInput = page.locator('.graphing-cell-input').first();
    await cellInput.type('x^2');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(300);

    const output = await page.locator('.graphing-cell-output').first().textContent();
    check('Graphing cell evaluates expression', output && output.length > 5, `output: "${output.substring(0, 30)}..."`);

    // 8. Dark Mode - switch back to app mode first
    console.log('\n--- Dark Mode ---');
    await page.click('#graphing-toggle');
    await page.waitForTimeout(200);

    // We're now in desktop mode, switch to app mode
    await page.click('#mode-toggle');
    await page.waitForTimeout(300);

    await page.click('.app-card:has-text("Settings")');
    await page.waitForTimeout(300);

    check('Settings tool works', await page.locator('#tool-settings').isVisible());

    await page.locator('label[for="dark-mode-toggle"]').click();
    await page.waitForTimeout(300);
    const darkOn = await page.evaluate(() => document.body.classList.contains('dark-mode'));
    check('Dark mode toggle works', darkOn);

    // 9. Bottom Bar
    console.log('\n--- Mobile Bottom Bar ---');
    await page.locator('#tool-settings .back-btn').click();
    await page.waitForTimeout(200);

    const bottomBar = page.locator('.mobile-bottom-bar');
    check('Bottom bar visible', await bottomBar.isVisible());

    await page.click('.mobile-bottom-bar button[data-tool="calculator"]');
    await page.waitForTimeout(300);
    check('Calculator opens from bottom bar', await page.locator('#tool-calculator').isVisible());

    const calcActive = await page.evaluate(() =>
      document.querySelector('.mobile-bottom-bar button[data-tool="calculator"]').classList.contains('active')
    );
    check('Bottom bar button shows active state', calcActive);

    // 10. Small Phone (320px)
    console.log('\n--- Small Phone (320px) ---');
    await page.locator('#tool-calculator .back-btn').click();
    await page.waitForTimeout(200);
    await page.setViewportSize({ width: 320, height: 568 });
    await page.waitForTimeout(300);

    const smallBody = await page.evaluate(() => document.body.clientWidth);
    check('Small viewport applied', smallBody <= 325, `width: ${smallBody}`);

    const headerBox = await page.locator('h1').boundingBox();
    check('Header fits on small phone', headerBox && headerBox.x + headerBox.width <= 325, `right: ${headerBox?.x + headerBox?.width}`);

    const cardBox = await page.locator('.app-card').first().boundingBox();
    check('App card fits on small phone', cardBox && cardBox.x + cardBox.width <= 325, `right: ${cardBox?.x + cardBox?.width}`);

    // 11. Landscape
    console.log('\n--- Landscape Mode ---');
    await page.setViewportSize({ width: 812, height: 375 });
    await page.waitForTimeout(300);

    check('Mobile mode class persists in landscape', await page.evaluate(() => document.body.classList.contains('mobile-mode')));
    check('App grid visible in landscape', await page.locator('#app-grid').isVisible());

    // 12. Action Buttons on mobile
    console.log('\n--- Action Buttons ---');
    await page.setViewportSize({ width: 375, height: 812 });
    await page.waitForTimeout(200);

    await page.click('.app-card:has-text("Algebra")');
    await page.waitForTimeout(300);
    await page.click('#tool-algebra .tab-btn:has-text("Advanced")');
    await page.waitForTimeout(200);

    const actionBtns = await page.locator('#tab-alg-advanced .action-btn').all();
    for (let i = 0; i < Math.min(actionBtns.length, 3); i++) {
      const box = await actionBtns[i].boundingBox();
      check(`Action button ${i + 1} fits`, box && box.x + box.width <= 380, box ? `right: ${box.x + box.width}` : 'null');
    }

    // 13. Console Errors
    console.log('\n--- Console Errors ---');
    const criticalErrors = consoleErrors.filter(e =>
      !e.includes('Warning') && !e.includes('deprecated') && !e.includes('MathJax') && !e.includes('Failed to load')
    );
    check('No critical console errors', criticalErrors.length === 0, criticalErrors.length > 0 ? criticalErrors[0].substring(0, 100) : '');

  } catch (error) {
    console.error('Test failed:', error.message);
    bugs.push({ name: 'Test error', details: error.message });
  } finally {
    console.log('\n=== Summary ===');
    const passed = results.filter(r => r.status === 'PASS').length;
    const failed = results.filter(r => r.status === 'FAIL').length;
    console.log(`Passed: ${passed}/${results.length}`);
    console.log(`Failed: ${failed}/${results.length}`);

    if (bugs.length > 0) {
      console.log('\n=== Bugs Found ===');
      bugs.forEach(b => {
        console.log(`- ${b.name}${b.details ? ': ' + b.details : ''}`);
      });
    } else {
      console.log('\nNo bugs found!');
    }

    await browser.close();
    process.exit(failed > 0 ? 1 : 0);
  }
})();
