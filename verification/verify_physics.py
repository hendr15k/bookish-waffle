from playwright.sync_api import sync_playwright

def verify_physics(page):
    import os
    file_path = "file://" + os.path.abspath("index.html")

    page.on("console", lambda msg: print(f"CONSOLE: {msg.text}"))
    page.on("pageerror", lambda err: print(f"PAGE ERROR: {err}"))

    page.goto(file_path)

    page.click("#mode-toggle")
    page.wait_for_selector(".app-card:has-text(\"Physics\")")
    page.click(".app-card:has-text(\"Physics\")")

    page.wait_for_selector("#phys-u", state="visible")

    page.fill("#phys-u", "10")
    page.fill("#phys-a", "2")
    page.fill("#phys-t", "5")

    u_val = page.input_value("#phys-u")
    print(f"u: {u_val}")

    page.click("#tool-physics .action-btn:has-text(\"Calculate\")")

    page.wait_for_timeout(500)

    is_visible = page.is_visible("#tool-physics .tool-result")
    print(f"Result visible: {is_visible}")

    if not is_visible:
        # Check if element exists
        el = page.query_selector("#tool-physics .tool-result")
        if el:
            display = page.eval_on_selector("#tool-physics .tool-result", "el => el.style.display")
            print(f"Result display style: {display}")
        else:
            print("Element not found")

    result_text = page.text_content("#tool-physics .tool-result")
    print(f"Result Text: {result_text}")

    page.screenshot(path="verification/physics_tool.png")

with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    verify_physics(page)
    browser.close()
