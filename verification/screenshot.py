
import os
from playwright.sync_api import sync_playwright

def run(playwright):
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page()

    # Get absolute path to index.html
    cwd = os.getcwd()
    url = f"file://{cwd}/index.html"
    print(f"Loading {url}")

    page.goto(url)

    # Enable Mobile Mode to see icons (optional, but good for visual context)
    # This makes the app grid visible
    page.click("#mode-toggle")

    # 1. Test Regression Plot (Statistics Tool)
    print("Testing Regression Plot...")

    # Navigate using JS for reliability
    page.evaluate("openTool('statistics')")

    # Wait for tool to appear
    page.wait_for_selector("#tool-statistics", state="visible")

    # Enter data
    page.fill("#reg-data", "[[0,0], [1,2], [2,4], [3,6]]")

    # Click Calculate & Plot
    page.click("button:has-text('Calculate & Plot')")

    # Wait for result
    page.wait_for_selector("#tool-statistics .tool-result canvas", state="visible")

    # Take screenshot of the tool view
    page.screenshot(path="verification/screenshot_regression.png")
    print("Screenshot saved to verification/screenshot_regression.png")

    # 2. Test Matrix Builder (Linear Algebra Tool)

    print("Testing Matrix Builder...")
    page.evaluate("openTool('linear')")

    page.wait_for_selector("#tool-linear", state="visible")

    # Click Generate Grid
    page.click("button:has-text('Generate Grid')")

    # Wait for grid
    page.wait_for_selector("#matrix-builder-container table")

    # Fill top-left cell (0,0)
    page.fill("input[data-row='0'][data-col='0']", "5")
    page.fill("input[data-row='0'][data-col='1']", "7")

    # Click Insert Matrix
    page.click("button:has-text('Insert Matrix')")

    # Check if value inserted into #lin-expr
    val = page.input_value("#lin-expr")
    print(f"Inserted Matrix: {val}")

    # Take screenshot
    page.screenshot(path="verification/screenshot_matrix.png")
    print("Screenshot saved to verification/screenshot_matrix.png")

    browser.close()

with sync_playwright() as playwright:
    run(playwright)
