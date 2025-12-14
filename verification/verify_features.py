
from playwright.sync_api import sync_playwright
import time

def verify_features():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto("http://localhost:8080/index.html")

        # 1. Switch to App Mode
        page.click("#mode-toggle")

        # 2. Open Linear Algebra Tool and verify extended tabs
        print("Checking Linear Algebra...")
        page.click(".app-card:has-text('Linear')")
        # Check for 'Binary' tab instead of 'Vectors'
        page.wait_for_selector("button:has-text('Binary')")
        page.click("button:has-text('Binary')")
        page.wait_for_selector("button:has-text('Kronecker Prod')")
        # Check Operations for SVD
        page.click("button:has-text('Operations')")
        page.wait_for_selector("button:has-text('SVD')")
        page.screenshot(path="verification/linear_features.png")

        # Use specific back button for Linear Algebra tool
        # The selector .back-btn matches multiple elements (one for each tool).
        # We need the visible one.
        page.click("#tool-linear .back-btn")

        # 3. Open Calculus Tool and verify Transforms
        print("Checking Calculus...")
        page.click(".app-card:has-text('Calculus')")
        page.click("button:has-text('Transforms')")
        page.wait_for_selector("button:has-text('Laplace L{f}')")
        page.screenshot(path="verification/calculus_features.png")
        page.click("#tool-calculus .back-btn")

        # 4. Open Statistics Tool and verify Mode/Multi-Var
        print("Checking Statistics...")
        page.click(".app-card:has-text('Statistics')")
        page.wait_for_selector("button:has-text('Mode')")
        page.click("button:has-text('Multi-Var')")
        page.wait_for_selector("button:has-text('Covariance')")
        page.screenshot(path="verification/statistics_features.png")

        browser.close()

if __name__ == "__main__":
    verify_features()
