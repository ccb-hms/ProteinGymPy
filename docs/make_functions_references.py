import inspect
import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / 'src'
DOCS_REFERENCE_DIR = REPO_ROOT / 'docs' / 'reference'

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import proteingympy  # noqa: E402


def get_exported_functions():
    exported_names = getattr(proteingympy, '__all__', [])
    if exported_names:
        return [
            name for name in exported_names
            if inspect.isfunction(getattr(proteingympy, name, None))
        ]
    # Fallback to all module-level callables if __all__ is missing/empty
    return [name for name, obj in inspect.getmembers(proteingympy, inspect.isfunction)]


functions = get_exported_functions()

os.makedirs(DOCS_REFERENCE_DIR, exist_ok=True)

for func in functions:
    doc_path = DOCS_REFERENCE_DIR / f'{func}.md'
    with open(doc_path, 'w') as f:
        f.write(f'# {func}\n\n::: proteingympy.{func}\n')