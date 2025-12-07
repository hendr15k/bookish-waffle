import os
from playwright.sync_api import sync_playwright

def run(playwright):
    browser = playwright.chromium.launch()
    page = browser.new_page()

    cwd = os.getcwd()
    url = f"file://{cwd}/index.html"
    page.goto(url)

    page.click("#mode-toggle")

    print("Testing Stats...")
    page.click(".app-card:has-text('Statistics')")
    page.click("button:has-text('Hypothesis Test')")
    page.fill("#test-data", "[10, 10, 10, 10]")
    page.fill("#test-mu", "9")
    page.fill("#test-sigma", "2")
    page.click("#tool-statistics button:has-text('Run Z-Test')")

    page.wait_for_function("document.querySelector('#tool-statistics .tool-result').innerText.length > 5")
    page.screenshot(path="verification/stats_ztest.png")

    page.click(".tool-view:visible .back-btn")

    print("Testing Geometry...")
    page.click(".app-card:has-text('Geometry')")
    page.click("button:has-text('Analytic')")
    page.fill("#geo-p1", "[0, 0]")
    page.fill("#geo-p2", "[3, 4]")
    page.click("#tool-geometry button:has-text('Distance')")

    page.wait_for_selector("#tool-geometry .tool-result", state="visible")
    page.wait_for_function("document.querySelector('#tool-geometry .tool-result').innerText.includes('5')")
    page.screenshot(path="verification/geo_analytic.png")

    page.click(".tool-view:visible .back-btn")

    print("Testing Linear...")
    page.click(".app-card:has-text('Linear')")
    page.screenshot(path="verification/linear_buttons.png")

    print("Testing Plot...")
    page.click(".tool-view:visible .back-btn")
    page.click(".app-card:has-text('Plotting')")
    page.fill("#plot-expr", "sin(x)")
    page.click("#tool-plotting button:has-text('Plot')")

    # Wait for result in the plotting tool
    page.wait_for_selector("#tool-plotting .tool-result canvas", state="visible")

    # The button is inside the tool-result (copied from history)
    # The text is "💾".
    page.wait_for_selector("#tool-plotting .tool-result button:has-text('💾')", state="visible")
    page.screenshot(path="verification/plot_save.png")

    browser.close()

with sync_playwright() as playwright:
    run(playwright)
