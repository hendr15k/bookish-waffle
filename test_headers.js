const { chromium } = require('playwright');

(async () => {
  const browser = await chromium.launch({ headless: true });
  
  // Test on iPhone SE (320x568)
  const context1 = await browser.newContext({ viewport: { width: 320, height: 568 }, isMobile: true, hasTouch: true });
  const page1 = await context1.newPage();
  await page1.goto('http://localhost:9999/', { waitUntil: 'networkidle', timeout: 15000 });

  const header1 = await page1.evaluate(() => {
    const header = document.body.firstElementChild;
    const rect = header.getBoundingClientRect();
    const style = getComputedStyle(header);
    return { height: rect.height, width: rect.width, y: rect.y, position: style.position, zIndex: style.zIndex };
  });
  console.log('iPhone SE (320x568) header:', JSON.stringify(header1));

  await context1.close();

  // Test on iPhone X (375x812)
  const context2 = await browser.newContext({ viewport: { width: 375, height: 812 }, isMobile: true, hasTouch: true });
  const page2 = await context2.newPage();
  await page2.goto('http://localhost:9999/', { waitUntil: 'networkidle', timeout: 15000 });

  const header2 = await page2.evaluate(() => {
    const header = document.body.firstElementChild;
    const rect = header.getBoundingClientRect();
    return { height: rect.height, width: rect.width, y: rect.y };
  });
  console.log('iPhone X (375x812) header:', JSON.stringify(header2));

  // Test smaller phone (280x560) - ultra small
  const context3 = await browser.newContext({ viewport: { width: 280, height: 560 }, isMobile: true, hasTouch: true });
  const page3 = await context3.newPage();
  await page3.goto('http://localhost:9999/', { waitUntil: 'networkidle', timeout: 15000 });

  const header3 = await page3.evaluate(() => {
    const header = document.body.firstElementChild;
    const rect = header.getBoundingClientRect();
    return { height: rect.height, width: rect.width, y: rect.y };
  });
  console.log('Ultra small (280x560) header:', JSON.stringify(header3));

  await context3.close();
  await browser.close();
})();