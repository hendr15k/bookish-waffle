from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        # Load local file
        url = "file://" + os.path.abspath("index.html")
        page.goto(url)

        # 1. Switch to App Mode
        page.click("#mode-toggle")

        # 2. Verify Table Tool
        # Open Table Tool
        page.click("text=Table") # Assuming card has text 'Table'

        # Fill inputs
        page.fill("#table-expr", "x^2")
        page.fill("#table-start", "0")
        page.fill("#table-end", "5")
        page.fill("#table-step", "1")

        # Generate - use specific selector or onclick since there are two "Generate Table" buttons (Logic and Table)
        # Table tool button calls runTable()
        page.click('button[onclick="runTable()"]')

        # Wait for result
        page.wait_for_selector("#tool-table .tool-result table")

        # Take screenshot of Table Tool
        page.screenshot(path="verification/table_tool.png")
        print("Table Tool Verified")

        # 3. Verify F-Distribution
        page.click("#tool-table .back-btn")

        # Open Statistics
        page.click("text=Statistics")

        # Switch to Distributions Tab
        page.click("text=Distributions")

        # Select F-Distribution
        page.select_option("#stat-dist-type", "fdist")

        # Check inputs visible
        if page.is_visible("#dist-inputs-fdist"):
            print("F-Distribution Inputs Visible")
        else:
            print("F-Distribution Inputs NOT Visible")

        # Take screenshot
        page.screenshot(path="verification/f_dist.png")

        # 4. Verify Matrix CSV
        page.click("#tool-statistics .back-btn")
        page.click("text=Linear")
        page.click("text=Builder")

        # Check buttons
        if page.is_visible("text=Import CSV") and page.is_visible("text=Export CSV"):
            print("CSV Buttons Visible")

        page.screenshot(path="verification/matrix_csv.png")

        # 5. Verify Session Save
        page.click("#tool-linear .back-btn")
        page.click("text=Variables")

        if page.is_visible("text=Save Variables") and page.is_visible("text=Load Variables"):
            print("Session Buttons Visible")

        page.screenshot(path="verification/session.png")

        browser.close()

if __name__ == "__main__":
    run()
