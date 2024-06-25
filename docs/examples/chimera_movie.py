"""
One script to rule them all!
... and make a chimeraX movie of an MD trajectory
... and also make a gif of the movie (outside of chimeraX)
"""

try:
    from chimerax.core.commands import run
except:
    run = None

from pathlib import Path

try:
    import imageio
except:
    imageio = None

DIRECTORY = "/Users/noahhk/Downloads/frames/"
DIRECTORY = Path(DIRECTORY)
DIRECTORY.mkdir(exist_ok=True)

if run is not None:

    run(session, "hide")

    for i in range(1, 101):
        run(session, f"select #1.{i}; show sel atoms; show sel surface; select clear")
        run(session, f"turn y {360/100}")
        run(session, f"save {DIRECTORY}/frame_{i}.png")
        run(session, f"select #1.{i}; hide sel atoms; hide sel surface")

if imageio is not None:

    images = DIRECTORY.glob("*.png")
    images = sorted(images, key=lambda x: int(x.stem.split("_")[-1]))
    images = [imageio.imread(str(image)) for image in images]
    imageio.mimsave(DIRECTORY / "movie.gif", images, duration=0.1)

print("done")
