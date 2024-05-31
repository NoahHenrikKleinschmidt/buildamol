"""
Make a new version of biobuild
"""

import argparse
import os
import re
import subprocess

BASE_DIR = os.path.abspath(__file__)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(BASE_DIR)))
CODE_DIR = BASE_DIR + "/buildamol"
DOCS_DIR = BASE_DIR + "/docs"

FILES_TO_UPDATE = {
    CODE_DIR + "/utils/info.py": "inpackage",
    BASE_DIR + "/setup.py": "setup",
    DOCS_DIR + "/conf.py": "docs",
}


VERSION_PATTERN = r"\"(\d+\.\d+\.\d+-?(a|b|dev)?)\""

REPLACEMENT_PATTERNS = {
    "inpackage": "__version__ ?= ?{}\n",
    "setup": "version ?= ?{},",
    "docs": "release ?= ?{}",
}

DEFAULT_COMMIT_MESSAGE = "[automatic] Version update to {}"


def version_to_tuple(version: str):
    """
    Convert a version string to a tuple

    Parameters
    ----------
    version : str
        The version string

    Returns
    -------
    tuple
        The version tuple
    """
    return tuple(map(int, version.split("-")[0].split(".")))


def get_current_version(filename: str, pattern: str):
    """
    Get the current version of biobuild

    Parameters
    ----------
    filename : str
        The path to the file to open

    Returns
    -------
    str
        The current version
    """
    pattern = REPLACEMENT_PATTERNS[pattern].format(VERSION_PATTERN)
    pattern = re.compile(pattern)
    with open(filename) as f:
        for line in f:
            match = re.search(pattern, line)
            if match:
                return match.group(1)
    raise ValueError(f"Could not find version in {filename}")


def make_new_version(
    filename: str,
    pattern: str,
    args: argparse.Namespace,
) -> str:
    """
    Make a new version
    """
    current_version = get_current_version(filename, pattern)
    increments = parse_increment(args)
    if increments[0] is not None:
        new_version = increment_new_version(
            current_version,
            *increments,
        )
    else:
        new_version = args.version
    if args.alpha:
        new_version += "-a"
    elif args.beta:
        new_version += "-b"
    elif args.dev:
        new_version += "-dev"
    return new_version


def increment_new_version(
    current_version: str,
    top_increment: int,
    mid_increment: int,
    base_increment: int,
):
    """
    Get the new version of biobuild

    Parameters
    ----------
    current_version : str
        The current version
    top_increment : int
        The amount to increment the top by
    mid_increment : int
        The amount to increment the mid by
    base_increment : int
        The amount to increment the base by

    Returns
    -------
    str
        The new version
    """
    top, mid, base = version_to_tuple(current_version)
    if top_increment == -1:
        top = 0
    else:
        top += top_increment
    if mid_increment == -1:
        mid = 0
    else:
        mid += mid_increment
    if base_increment == -1:
        base = 0
    else:
        base += base_increment
    new = f"{top}.{mid}.{base}"
    return new


def update_file(filename: str, pattern: str, new_version: str):
    """
    Update a file with a new version

    Parameters
    ----------
    filename : str
        The path to the file to open
    pattern : str
        The pattern to search for
    new_version : str
        The new version to replace the old one with
    """
    pattern = REPLACEMENT_PATTERNS[pattern].format(VERSION_PATTERN)
    pattern = re.compile(pattern)
    with open(filename) as f:
        lines = f.readlines()
    with open(filename, "w") as f:
        for line in lines:
            match = re.search(pattern, line)
            if match:
                line = line.replace(match.group(1), new_version)
            f.write(line)


def setup():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "version",
        help="The new version designation. This can be either '+'=increment base by one (e.g. '1.1.1' -> '1.1.2'), '++' to update the mid by one and reset the base (e.g. '1.1.4' -> '1.2.0'), or '+++' to increment the top by one and reset mid and base (e.g. '1.4.2' -> '2.0.0'). Use '-' to leave the version unchanged (e.g. when you only want to remove an alpha or beta tag). Alternatively, a specific string can be provided as the version directly (e.g. '5.2.1').",
    )
    parser.add_argument(
        "-a",
        "--alpha",
        action="store_true",
        help="mark version as an alpha-release (add -a to the end of the version)",
    )
    parser.add_argument(
        "-b",
        "--beta",
        action="store_true",
        help="mark version as a beta-release (add -b to the end of the version)",
    )
    parser.add_argument(
        "-d",
        "--dev",
        action="store_true",
        help="mark version as a development release (add -dev to the end of the version)",
    )
    parser.add_argument(
        "-u",
        "--build-docs",
        action="store_true",
        help="build a new version of the documentation. This requires a directory named 'docs' wherein the sphinx conf.py and makefile are stored is located in the base directory (same as where the main package directory and setup.py are stored).",
    )
    parser.add_argument(
        "-i", "--install", action="store_true", help="pip-install the new version"
    )
    parser.add_argument(
        "-g", "--git", action="store_true", help="commit the updated files"
    )
    parser.add_argument(
        "-c",
        "--commit-message",
        default=None,
        help="the commit message to use when committing the updated files. A default message will be used if this is not provided.",
    )

    args = parser.parse_args()
    return args


def parse_increment(args):
    if args.version == "+++":
        return 1, -1, -1
    elif args.version == "++":
        return 0, 1, -1
    elif args.version == "+":
        return 0, 0, 1
    elif args.version == "-":
        return 0, 0, 0
    else:
        return None, None, None


def build_docs():
    """
    Build the documentation
    """
    import add_example_gallery

    os.chdir(DOCS_DIR)
    g = add_example_gallery.gallery("_static/gallery/")
    with open(DOCS_DIR + "/index.rst", "r") as f:
        contents = f.read()
    with open(DOCS_DIR + "/index.rstb", "w") as f:
        f.write(contents)
    with open(DOCS_DIR + "/index.rst", "w") as f:
        f.write(contents.replace(".. gallery", g))

    subprocess.run("make clean", shell=True)
    subprocess.run("make html", shell=True)
    os.remove(DOCS_DIR + "/index.rst")
    os.rename(DOCS_DIR + "/index.rstb", DOCS_DIR + "/index.rst")


def install():
    """
    Install the package
    """
    subprocess.run("pip install " + BASE_DIR, shell=True)


def commit(args, new_version):
    """
    Commit the changes
    """
    for filename in FILES_TO_UPDATE:
        subprocess.run("git add " + filename, shell=True)
    if args.commit_message is None:
        args.commit_message = DEFAULT_COMMIT_MESSAGE.format(new_version)
    subprocess.run(f"git commit -m '{args.commit_message}'", shell=True)


def main(args):
    for filename, pattern in FILES_TO_UPDATE.items():
        new_version = make_new_version(filename, pattern, args)
        update_file(filename, pattern, new_version)
        print(f"Updated {filename} to {new_version}")
    if args.install:
        install()
    if args.build_docs:
        build_docs()
    if args.git:
        commit(args, new_version)
    print("Update to {} complete".format(new_version))


if __name__ == "__main__":
    args = setup()
    main(args)

    # file = list(FILES_TO_UPDATE.keys())[0]
    # pattern = "inpackage"

    # version = get_current_version(file, pattern)

    # class Args:
    #     version = "-"
    #     alpha = False
    #     beta = False

    # args = Args()
    # new = make_new_version(file, pattern, args)
    # print(version, "->", new)
    # update_file(file, pattern, new)
