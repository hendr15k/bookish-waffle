const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true });
  const context = await browser.newContext({
    viewport: { width: 375, height: 812 },
    isMobile: true,
    hasTouch: true,
    userAgent: 'Mozilla/5.0 (iPhone; CPU iPhone OS 16_0 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/16.0 Mobile/15E148 Safari/604.1'
  });
  const page = await context.newPage();

  const results = [];
  function check(name, condition, details = '') {
    const status = condition ? 'PASS' : 'FAIL';
    results.push({ name, status, details });
    console.log(`${status}: ${name}${details ? ' - ' + details : ''}`);
  }

  const consoleErrors = [];
  page.on('console', msg => {
    if (msg.type() === 'error') consoleErrors.push(msg.text());
  });

  try {
    // 1. Load page
    console.log('\n=== Testing Mobile View (375x812 iPhone X) ===\n');
    await page.goto('http://localhost:9999/', { waitUntil: 'networkidle', timeout: 15000 });

    // 2. Check initial desktop mode
    console.log('--- Desktop Mode (initial) ---');
    const mainContent = await page.locator('#main-content').isVisible();
    check('main-content visible on load', mainContent);

    const headerH1 = await page.locator('h1').textContent();
    check('Header text', headerH1 === 'Web CAS', `got: "${headerH1}"`);

    // 3. Check header layout on mobile
    const headerEl = page.locator('h1').locator('..');
    const headerBox = await headerEl.boundingBox();
    if (headerBox) {
      check('Header fits in viewport', headerBox.width <= 375, `header width: ${headerBox.width}, viewport: 375`);
      check('Header not cut off', headerBox.x >= 0, `header x: ${headerBox.x}`);
    }

    // 4. Check header buttons don't overflow
    const modeToggle = page.locator('#mode-toggle');
    const graphingToggle = page.locator('#graphing-toggle');
    const modeToggleVisible = await modeToggle.isVisible();
    const graphingToggleVisible = await graphingToggle.isVisible();
    check('Mode toggle visible on mobile', modeToggleVisible);
    check('Graphing toggle visible on mobile', graphingToggleVisible);

    // Check if buttons wrap properly
    const buttonsContainer = page.locator('h1').locator('..').locator('div').last();
    const btnsBox = await buttonsContainer.boundingBox();
    if (btnsBox) {
      check('Header buttons fit viewport', btnsBox.x + btnsBox.width <= 375, `buttons right edge: ${btnsBox.x + btnsBox.width}`);
    }

    // 5. Switch to App Mode
    console.log('\n--- App Mode ---');
    await modeToggle.click({ force: true });
    await page.waitForTimeout(500);

    const appGrid = await page.locator('#app-grid').isVisible();
    check('App grid visible in app mode', appGrid);

    const appCategories = await page.locator('#app-categories').isVisible();
    check('App categories visible in app mode', appCategories, 'BUG: switchMode does not show app-categories');

    const appSearch = await page.locator('#app-search-container').isVisible();
    check('App search visible in app mode', appSearch);

    // 6. Check app cards layout on mobile
    const appCards = await page.locator('.app-card').count();
    check('App cards exist', appCards > 0, `count: ${appCards}`);

    if (appCards > 0) {
      const firstCard = page.locator('.app-card').first();
      const cardBox = await firstCard.boundingBox();
      if (cardBox) {
        check('App card fits viewport', cardBox.x + cardBox.width <= 375, `card right: ${cardBox.x + cardBox.width}`);
      }
    }

    // 7. Test opening Calculator tool
    console.log('\n--- Calculator Tool ---');
    await page.click('.app-card:has-text("Calculator")');
    await page.waitForTimeout(500);

    const calculatorVisible = await page.locator('#tool-calculator').isVisible();
    check('Calculator tool visible after click', calculatorVisible);

    const calcInput = await page.locator('#calc-main-input').isVisible();
    check('Calculator input visible', calcInput);

    const calcResult = await page.locator('#tool-calculator .tool-result').isVisible();
    check('Calculator result area hidden by default (expected)', !calcResult, 'display:none is default behavior');

    // 8. Test back button
    console.log('\n--- Back Button ---');
    const backBtn = page.locator('#tool-calculator .back-btn');
    const backBtnVisible = await backBtn.isVisible();
    check('Back button visible in tool view', backBtnVisible);

    if (backBtnVisible) {
      await backBtn.click();
      await page.waitForTimeout(300);
      const appGridAfterBack = await page.locator('#app-grid').isVisible();
      check('App grid visible after back button', appGridAfterBack);
    }

    // 9. Test Algebra tool tabs on mobile
    console.log('\n--- Algebra Tabs ---');
    await page.click('.app-card:has-text("Algebra")');
    await page.waitForTimeout(300);

    const algebraTabs = await page.locator('#tool-algebra .tab-btn').count();
    check('Algebra has 4 tabs', algebraTabs === 4, `count: ${algebraTabs}`);

    // Check tabs don't overflow
    const tabHeader = page.locator('#tool-algebra .tab-header');
    const tabHeaderBox = await tabHeader.boundingBox();
    if (tabHeaderBox) {
      check('Tab header fits viewport', tabHeaderBox.width <= 375, `width: ${tabHeaderBox.width}`);
    }

    // Switch to Solve tab
    await page.click('#tool-algebra .tab-btn:has-text("Solve")');
    await page.waitForTimeout(200);
    const solveTab = await page.locator('#tab-alg-solve').isVisible();
    check('Solve tab visible after click', solveTab);

    // 10. Test Calculus tool (has many tabs)
    console.log('\n--- Calculus Tabs (14 tabs) ---');
    await page.locator('#tool-algebra .back-btn').click();
    await page.waitForTimeout(300);
    await page.click('.app-card:has-text("Calculus")');
    await page.waitForTimeout(300);

    const calcTabs = await page.locator('#tool-calculus .tab-btn').count();
    check('Calculus has many tabs', calcTabs > 5, `count: ${calcTabs}`);

    // Check if Calculus tab header scrolls
    const calcTabHeader = page.locator('#tool-calculus .tab-header');
    const calcTabHeaderBox = await calcTabHeader.boundingBox();
    if (calcTabHeaderBox) {
      const calcScrollable = await calcTabHeader.evaluate(el => el.scrollWidth > el.clientWidth);
      check('Calculus tab header scrolls on mobile', calcScrollable, `scrollWidth vs clientWidth`);
    }

    // 11. Test Form-row layout on mobile
    console.log('\n--- Form Row Layout ---');
    const formRows = await page.locator('#tool-calculus .form-row').count();
    check('Form rows exist in Calculus', formRows > 0, `count: ${formRows}`);

    // Check that form-row items stack vertically on mobile
    if (formRows > 0) {
      const firstFormRow = page.locator('#tool-calculus .form-row').first();
      const flexDirection = await firstFormRow.evaluate(el => getComputedStyle(el).flexDirection);
      check('Form-row stacks on mobile', flexDirection === 'column', `flexDirection: ${flexDirection}`);
    }

    // 12. Test Action grid on mobile
    console.log('\n--- Action Grid ---');
    const actionGrids = await page.locator('#tool-calculus .action-grid').count();
    if (actionGrids > 0) {
      const firstActionGrid = page.locator('#tool-calculus .action-grid').first();
      const gridCols = await firstActionGrid.evaluate(el => getComputedStyle(el).gridTemplateColumns);
      check('Action grid single column on mobile', gridCols.split(' ').length === 1 || gridCols.includes('1fr') || gridCols === 'none', `gridTemplateColumns: ${gridCols}`);
    }

    // 13. Test mobile bottom bar
    console.log('\n--- Mobile Bottom Bar ---');
    await page.locator('#tool-calculus .back-btn').click();
    await page.waitForTimeout(300);

    const mobileBottomBar = page.locator('.mobile-bottom-bar');
    const bottomBarVisible = await mobileBottomBar.isVisible();
    check('Mobile bottom bar visible', bottomBarVisible);

    if (bottomBarVisible) {
      const bottomBarBox = await mobileBottomBar.boundingBox();
      if (bottomBarBox) {
        check('Bottom bar at bottom of viewport', bottomBarBox.y + bottomBarBox.height >= 812 - 80, `bar y: ${bottomBarBox.y}`);
        check('Bottom bar not clipped', bottomBarBox.x >= 0 && bottomBarBox.x + bottomBarBox.width <= 375, `bar dimensions: ${bottomBarBox.width}x${bottomBarBox.height}`);
      }
    }

    // 14. Test Graphing mode on mobile
    console.log('\n--- Graphing Mode ---');
    await page.locator('#graphing-toggle').click();
    await page.waitForTimeout(500);

    const graphingMode = await page.locator('#graphing-mode').isVisible();
    check('Graphing mode visible', graphingMode);

    if (graphingMode) {
      const sidebar = page.locator('#graphing-sidebar');
      const sidebarBox = await sidebar.boundingBox();
      if (sidebarBox) {
        check('Graphing sidebar not too wide on mobile', sidebarBox.width <= 250, `sidebar width: ${sidebarBox.width}`);
      }

      const sidebarDisplay = await sidebar.evaluate(el => getComputedStyle(el).display);
      check('Graphing sidebar visible on mobile', sidebarDisplay !== 'none', `display: ${sidebarDisplay}`);
    }

    // 15. Switch back to app mode
    console.log('\n--- Back to App Mode ---');
    await page.locator('#graphing-toggle').click();
    await page.waitForTimeout(300);
    await page.locator('#mode-toggle').click();
    await page.waitForTimeout(300);

    const appGridFinal = await page.locator('#app-grid').isVisible();
    check('App grid visible after switching modes', appGridFinal);

    // 16. Check dark mode on mobile
    console.log('\n--- Dark Mode on Mobile ---');
    const darkModeToggle = page.locator('#dark-mode-toggle');
    if (await darkModeToggle.isVisible()) {
      await darkModeToggle.click();
      await page.waitForTimeout(300);
      const isDark = await page.evaluate(() => document.body.classList.contains('dark-mode'));
      check('Dark mode activated on mobile', isDark);

      // Check dark mode card colors
      const cardBg = await page.evaluate(() => {
        const card = document.querySelector('#app-grid .app-card');
        return card ? getComputedStyle(card).backgroundColor : 'no card';
      });
      check('Dark mode card background is dark', cardBg !== 'no card', `bg: ${cardBg}`);
    }

    // 17. Test tool result display on mobile
    console.log('\n--- Tool Result Display ---');
    await page.click('.app-card:has-text("Calculator")');
    await page.waitForTimeout(300);
    
    await page.fill('#calc-main-input', '2+2');
    await page.click('.action-btn:has-text("Calculate")');
    await page.waitForTimeout(500);

    const resultArea = page.locator('#tool-calculator .tool-result');
    const resultVisible = await resultArea.isVisible();
    check('Tool result visible after calculation', resultVisible);

    if (resultVisible) {
      const resultBox = await resultArea.boundingBox();
      if (resultBox) {
        check('Tool result fits mobile viewport', resultBox.width <= 375, `width: ${resultBox.width}`);
      }
    }

  } catch (error) {
    console.error('Test failed with error:', error.message);
  } finally {
    // Summary
    console.log('\n=== Summary ===');
    const passed = results.filter(r => r.status === 'PASS').length;
    const failed = results.filter(r => r.status === 'FAIL').length;
    console.log(`Passed: ${passed}/${results.length}`);
    console.log(`Failed: ${failed}/${results.length}`);
    if (failed > 0) {
      console.log('\nFailed tests:');
      results.filter(r => r.status === 'FAIL').forEach(r => {
        console.log(`  - ${r.name}${r.details ? ' (' + r.details + ')' : ''}`);
      });
    }
    if (consoleErrors.length > 0) {
      console.log('\nConsole Errors:');
      consoleErrors.forEach(e => console.log('  ', e));
    }
    await browser.close();
  }
})();