#!/usr/bin/env python3

import argparse
import http.server
import socketserver
import webbrowser
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Serve pipeline HTML reports locally."
    )
    parser.add_argument(
        "--results-dir",
        default="results",
        help="Results directory to serve (default: results)",
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host to bind to (default: 127.0.0.1)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=4174,
        help="Port to use (default: 4174)",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not automatically open the browser",
    )

    args = parser.parse_args()

    results_dir = Path(args.results_dir).resolve()

    if not results_dir.is_dir():
        raise SystemExit(
            f"Results directory not found: {results_dir}"
        )

    dashboard = results_dir / "run_summary" / "run_summary.html"

    if not dashboard.is_file():
        raise SystemExit(
            f"Run summary not found: {dashboard}\n"
            "Generate the run summary before starting the report server."
        )

    handler = lambda *handler_args, **handler_kwargs: http.server.SimpleHTTPRequestHandler(
        *handler_args,
        directory=str(results_dir),
        **handler_kwargs,
    )

    url = (
        f"http://{args.host}:{args.port}/"
        "run_summary/run_summary.html"
    )

    print()
    print("Pipeline report server")
    print("======================")
    print(f"Serving: {results_dir}")
    print(f"Dashboard: {url}")
    print()
    print("Press Ctrl+C to stop.")
    print()

    class ReusableTCPServer(socketserver.TCPServer):
        allow_reuse_address = True

    with ReusableTCPServer((args.host, args.port), handler) as server:
        if not args.no_browser:
            webbrowser.open(url)

        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print("\nReport server stopped.")


if __name__ == "__main__":
    main()