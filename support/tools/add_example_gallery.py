import os

_style_str = """
.. raw:: html

    <style>
    .slideshow-container {
      max-width: 1000px;
      position: relative;
      margin: auto;
    }

   img {
      background-color: transparent !important;
   }

    .mySlides {
      display: none;
    }

    .fade {
      animation: fade <animdur>ms infinite;
    }

    @keyframes fade {
      from {opacity: .7} 
      to {opacity: 1}
    }
    </style>

.. raw:: html

    <script>
    var slideIndex = 0;
    showSlides();

    function showSlides() {
      var i;
      var slides = document.getElementsByClassName("mySlides");
      for (i = 0; i < slides.length; i++) {
        slides[i].style.display = "none";  
      }
      slideIndex++;
      if (slideIndex > slides.length) {slideIndex = 1}    
      slides[slideIndex-1].style.display = "block";  
      setTimeout(showSlides, <animspeed>); // Change image every N seconds
    }
    </script>
""".strip()


def gallery(directory):
    """
    Add a gallery to the documentation

    Parameters
    ----------
    directory : str
        The directory containing the images
    """
    images = [directory + f for f in os.listdir(directory) if f.endswith(".png")]
    return f"""
.. raw:: html

    <div class="slideshow-container">
        {"".join([slide(f) for f in images])}
    </div>

{style()}
""".strip()


def slide(fname):
    """
    Add a slide to the documentation

    Parameters
    ----------
    fname : str
        The filename of the image to use
    """
    return f"""
    <div class="mySlides fade">
        <img src="{fname}" style="width:100%">
    </div>
"""


def style(duration: int = 5000, speed: int = 5000):
    """
    Add a style to the documentation

    Parameters
    ----------
    duration : int, optional
        The duration of the fade, by default 2000
    speed : int, optional
        The speed of the animation, by default 5000
    """
    return _style_str.replace("<animdur>", str(duration)).replace(
        "<animspeed>", str(speed)
    )
