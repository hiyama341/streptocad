How to run StreptoCAD locally or hack your own workflows.

Published paper: https://pubs.acs.org/doi/10.1021/acssynbio.5c00261

## Run locally

From the project root:

```bash
pip install -e .
python web_app/application.py
```

Then open the local URL shown in the terminal, usually `http://127.0.0.1:8050`.

## Hack workflows

- Main app entry point: `web_app/application.py`
- Workflow callbacks: `web_app/callbacks/`
- Workflow tabs/layouts: `web_app/tabs/`
- Shared styling and components: `web_app/styling.py`, `web_app/components.py`, `web_app/tooltip.py`

If you add or change a workflow, restart the app and refresh the browser.
