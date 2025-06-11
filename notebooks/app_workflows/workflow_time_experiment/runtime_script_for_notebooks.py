import os
import time
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import pandas as pd
from tqdm import tqdm


def run_notebook(path: str, timeout: int = 600) -> float:
    """
    Execute a Jupyter notebook and measure its runtime.

    Args:
        path: Path to the .ipynb file.
        timeout: Maximum execution time in seconds for the notebook.

    Returns:
        Runtime in seconds.

    Raises:
        Exception: If execution fails.
    """
    with open(path, "r", encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)
    ep = ExecutePreprocessor(timeout=timeout, kernel_name="python3")
    start_time = time.time()
    ep.preprocess(nb, {"metadata": {"path": os.path.dirname(path) or "."}})
    end_time = time.time()
    return end_time - start_time


def measure_all_notebooks(
    directory: str, timeout: int = 600, show_progress: bool = True
) -> pd.DataFrame:
    """
    Walk through the given directory, execute all .ipynb files, and record runtimes.

    Args:
        directory: Root directory to search for notebooks.
        timeout: Execution timeout per notebook.
        show_progress: If True, display progress bar and current notebook being run.

    Returns:
        DataFrame indexed by notebook filename with columns ['runtime', 'error'].
    """
    records = []
    nb_paths = []
    for root, _, files in os.walk(directory):
        for fname in files:
            if fname.endswith(".ipynb"):
                nb_paths.append(os.path.join(root, fname))

    iterator = tqdm(nb_paths, desc="Running notebooks") if show_progress else nb_paths
    for fullpath in iterator:
        if show_progress:
            iterator.set_postfix(notebook=os.path.basename(fullpath))

        try:
            print(f"\n▶ Running {fullpath}")
            runtime = run_notebook(fullpath, timeout)
            print(f"✔ Completed in {runtime:.2f}s")
            error_msg = ""
        except Exception as e:
            runtime = None
            error_msg = str(e)
            print(f"✖ Error: {error_msg}")

        records.append(
            {
                "notebook": os.path.relpath(fullpath, directory),
                "runtime": runtime,
                "error": error_msg,
            }
        )

    df = pd.DataFrame(records).set_index("notebook")
    return df


if __name__ == "__main__":
    # This script is intended to be run from the command line.
    # example usage:
    # python runtime_script_for_notebooks.py /path/to/notebooks --timeout 6000
    import argparse

    print("Starting notebook execution script with progress reporting...")
    parser = argparse.ArgumentParser(
        description="Measure execution time of all Jupyter notebooks in a directory and save results to CSV."
    )
    parser.add_argument(
        "directory", type=str, help="Path to the directory containing notebooks"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="Max execution time per notebook in seconds (default: 600)",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bar and live output",
    )
    args = parser.parse_args()

    df = measure_all_notebooks(
        args.directory, timeout=args.timeout, show_progress=not args.no_progress
    )
    output_file = "runtime_experiment.csv"
    df.to_csv(output_file)
    print(f"Results saved to {output_file}")
