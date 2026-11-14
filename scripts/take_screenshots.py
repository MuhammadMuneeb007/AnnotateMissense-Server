#!/usr/bin/env python3
"""Drive a running AnnotateMissense server with a headless browser and save screenshots.

    python scripts/take_screenshots.py --base http://127.0.0.1:8000 --out docs/screenshots

Requires: pip install playwright  (uses the bundled Chromium, or --channel msedge/chrome).
"""

import argparse
from pathlib import Path

from playwright.sync_api import sync_playwright

LIST = "22 36305985 G A\n22 19963748 G A\n22 21549460 C A,G\n22 36284163 C T\n22 36284474 C T\n22 36305985 GG G\n22 10000000 A G"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default="http://127.0.0.1:8000")
    parser.add_argument("--out", default="docs/screenshots")
    parser.add_argument("--channel", default=None, help="msedge or chrome to use an installed browser")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    base = args.base.rstrip("/")

    with sync_playwright() as p:
        browser = p.chromium.launch(channel=args.channel) if args.channel else p.chromium.launch()
        page = browser.new_page(viewport={"width": 1400, "height": 900}, device_scale_factor=1)
        page.set_default_timeout(90_000)

        def shot(name: str, full: bool = True, locator: str = None) -> None:
            page.wait_for_timeout(900)  # let Chart.js animations finish
            target = out / f"{name}.png"
            if locator:
                page.locator(locator).screenshot(path=str(target))
            else:
                page.screenshot(path=str(target), full_page=full)
            print("saved", target)

        # 1. Home page with an inline single-variant result
        page.goto(f"{base}/", wait_until="networkidle")
        page.click("[data-example='22|36305985|G|A']")
        page.wait_for_selector("#single-result-preview a")
        shot("01-home-search")

        # 2. Full variant report (MYH9 p.R702C)
        page.goto(f"{base}/results?chr=22&pos=36305985&ref=G&alt=A", wait_until="networkidle")
        shot("02-variant-report")

        # 3. Transcript context for a multi-transcript allele
        page.goto(f"{base}/results?chr=22&pos=21549460&ref=C&alt=A", wait_until="networkidle")
        page.wait_for_timeout(500)
        shot("03-transcript-context", locator="#section-transcripts")

        # 4. VCF annotation
        page.goto(f"{base}/annotate", wait_until="networkidle")
        page.click("#loadExampleVcf")
        page.click("#annotate-form button[type=submit]")
        page.wait_for_selector("#annotateBody tr td.font-mono")
        shot("04-annotate-vcf")

        # 5. Variant-list batch lookup
        page.goto(f"{base}/", wait_until="networkidle")
        page.click("[data-home-tab='batch']")
        page.fill("#variants_text", LIST)
        page.click("#batch-upload-form button[type=submit]")
        page.wait_for_url("**/batch_results")
        page.wait_for_load_state("networkidle")
        shot("05-batch-results")

        # 6. Gene view
        page.goto(f"{base}/gene/MYH9?limit=25", wait_until="networkidle")
        shot("06-gene-myh9")

        # 7. API reference with a live call
        page.goto(f"{base}/api_docs", wait_until="networkidle")
        block = page.locator("#variant")
        block.locator(".run-get-tester").click()
        page.wait_for_function("document.querySelector('#variant .tester-output').textContent.includes('annotatemissense')")
        shot("07-api-docs", locator="#variant")

        # 8. Compare
        page.goto(f"{base}/compare?v1=22:36305985:G:A&v2=22:19963748:G:A&v3=22:36284163:C:T", wait_until="networkidle")
        page.click("#compareSubmit")
        page.wait_for_selector("#compareTableBody tr")
        shot("08-compare")

        browser.close()


if __name__ == "__main__":
    main()
