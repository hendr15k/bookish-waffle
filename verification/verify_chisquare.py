from playwright.sync_api import sync_playwright
import os

def run(playwright):
    browser = playwright.chromium.launch()
    page = browser.new_page()

    # Load the local index.html file
    # We use file:// protocol. Assuming current dir is root.
    cwd = os.getcwd()
    page.goto(f"file://{cwd}/index.html")

    # 1. Enable Mobile/App Mode to see the tools grid
    page.click("#mode-toggle")

    # 2. Click on Statistics Tool
    # Find the statistics app card.
    # It has text "Statistics".
    page.click("text=Statistics")

    # 3. Click on Distributions Tab
    # Button text "Distributions"
    page.click("text=Distributions")

    # 4. Select Chi-Squared from Dropdown
    page.select_option("#stat-dist-type", "chisquare")

    # 5. Wait for the plot to update/render
    # We can wait for the slider to be visible as a proxy that logic ran
    page.wait_for_selector("#dist-inputs-chisquare", state="visible")

    # Take a screenshot of the tool view
    page.screenshot(path="verification/chisquare_plot.png")

    browser.close()

with sync_playwright() as playwright:
    run(playwright)
