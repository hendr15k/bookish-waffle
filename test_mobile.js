const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true });
  const context = await browser.newContext({
    viewport: { width: 375, height: 812 }, // iPhone X viewport
    isMobile: true,
    hasTouch: true
  });
  const page = await context.newPage();

  console.log('=== Testing Mobile View ===\n');

  // Capture console errors
  const errors = [];
  page.on('console', msg => {
    if (msg.type() === 'error') errors.push(msg.text());
  });

  try {
    // 1. Load the page
    console.log('1. Loading page...');
    await page.goto('http://localhost:9000/', { waitUntil: 'networkidle' });
    console.log('   Page loaded successfully\n');

    // 2. Check desktop mode elements
    console.log('2. Checking desktop mode elements...');
    const header = await page.locator('h1').textContent();
    console.log('   Header:', header);

    const mainContent = await page.locator('#main-content').isVisible();
    console.log('   Main content visible:', mainContent);

    const inputArea = await page.locator('#input-area').isVisible();
    console.log('   Input area visible:', inputArea);
    console.log('');

    // 3. Switch to App Mode
    console.log('3. Switching to App Mode...');
    await page.click('#mode-toggle');
    await page.waitForTimeout(500);

    const appGrid = await page.locator('#app-grid').isVisible();
    console.log('   App grid visible:', appGrid);

    const appCategories = await page.locator('#app-categories').isVisible();
    console.log('   App categories visible:', appCategories);
    console.log('');

    // 4. Check app grid layout
    console.log('4. Checking app grid layout...');
    const appCards = await page.locator('.app-card').count();
    console.log('   Number of app cards:', appCards);
    console.log('');

    // 5. Open Calculator tool
    console.log('5. Opening Calculator tool...');
    await page.click('.app-card:has-text("Calculator")');
    await page.waitForTimeout(300);

    const calculatorVisible = await page.locator('#tool-calculator').isVisible();
    console.log('   Calculator tool visible:', calculatorVisible);

    const calcInput = await page.locator('#calc-main-input').isVisible();
    console.log('   Calculator input visible:', calcInput);

    const functionGrid = await page.locator('.function-grid').first().isVisible();
    console.log('   Function grid visible:', functionGrid);
    console.log('');

    // 6. Close tool and go back to grid
    console.log('6. Closing tool...');
    await page.click('.back-btn');
    await page.waitForTimeout(300);

    const appGridAfterClose = await page.locator('#app-grid').isVisible();
    console.log('   App grid visible after close:', appGridAfterClose);
    console.log('');

    // 7. Test Algebra tool
    console.log('7. Testing Algebra tool...');
    await page.click('.app-card:has-text("Algebra")');
    await page.waitForTimeout(300);

    const algebraVisible = await page.locator('#tool-algebra').isVisible();
    console.log('   Algebra tool visible:', algebraVisible);

    // Check tabs
    const tabs = await page.locator('#tool-algebra .tab-btn').count();
    console.log('   Number of tabs:', tabs);
    console.log('');

    // 8. Test Calculus tool
    console.log('8. Testing Calculus tool...');
    await page.click('.back-btn');
    await page.waitForTimeout(300);
    await page.click('.app-card:has-text("Calculus")');
    await page.waitForTimeout(300);

    const calcTabs = await page.locator('#tool-calculus .tab-btn').count();
    console.log('   Calculus tabs count:', calcTabs);
    console.log('');

    // 9. Check mobile bottom bar
    console.log('9. Checking mobile bottom bar...');
    const bottomBar = await page.locator('.mobile-bottom-bar').isVisible();
    console.log('   Bottom bar visible:', bottomBar);
    console.log('');

    // 10. Test dark mode toggle
    console.log('10. Testing dark mode...');
    await page.click('.back-btn');
    await page.waitForTimeout(300);
    const darkModeToggle = await page.locator('#dark-mode-toggle');
    if (await darkModeToggle.isVisible()) {
      await darkModeToggle.click();
      await page.waitForTimeout(300);
      const isDark = await page.evaluate(() => document.body.classList.contains('dark-mode'));
      console.log('   Dark mode activated:', isDark);
    } else {
      console.log('   Dark mode toggle not visible on mobile');
    }
    console.log('');

    // 11. Test Graphing mode
    console.log('11. Testing Graphing mode...');
    const graphingToggle = await page.locator('#graphing-toggle');
    await graphingToggle.click();
    await page.waitForTimeout(500);

    const graphingMode = await page.locator('#graphing-mode').isVisible();
    console.log('   Graphing mode visible:', graphingMode);
    console.log('');

    // Report errors
    console.log('=== Console Errors ===');
    if (errors.length === 0) {
      console.log('No console errors detected!');
    } else {
      errors.forEach(e => console.log('ERROR:', e));
    }

  } catch (error) {
    console.error('Test failed:', error.message);
  } finally {
    await browser.close();
  }
})();