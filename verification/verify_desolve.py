
from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        # Load local HTML file
        page.goto(f"file://{os.getcwd()}/index.html")

        # 1. Switch to App Mode
        page.click("#mode-toggle")

        # 2. Click Calculus App
        page.click(".app-card:has-text('Calculus')")

        # 3. Click Diff Eq Tab
        page.click("button:has-text('Diff Eq')")

        # 4. Enter equation
        # Note: diff(y,x) - y = 0 might be tricky to type if tokens not handled?
        # We can just type valid string.
        page.fill("#desolve-eq", "diff(y,x) - y")
        page.fill("#desolve-var", "y")

        # 5. Click Solve
        page.click("button:has-text('Solve ODE')")

        # 6. Wait for result in Calculus tool
        # Must be visible
        page.wait_for_selector("#tool-calculus .tool-result", state="visible")

        # 7. Screenshot
        page.screenshot(path="verification/desolve_test.png")

        browser.close()

if __name__ == "__main__":
    run()
