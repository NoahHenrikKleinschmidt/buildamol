.. <gallery>

.. raw:: html

    <div class="slideshow-container">
        
    <div class="mySlides fade">
        <img src="_static/gallery/sampling.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/complex.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/linear-polymer.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/conf_samples.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/polyphenylene.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/nanotube.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/xwing.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/prot-glyco.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/ring.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/MOF-top.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/MOF-angle.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/rotaxan-linear.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/rotaxan-small.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/peptide.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/drd2-ligand.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/MOF-detail.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/his20-2.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/his20.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/glyco-membrane.png" style="width:100%">
    </div>

    </div>

    <!-- Placeholder to maintain layout spacing -->
    <div class="slideshow-placeholder"></div>

.. raw:: html

    <style>
    /* Hide the "On this Page" sidebar on index page */
    .bd-toc {
      display: none !important;
    }

    /* Expand main content to use the full available width */
    .bd-main {
      grid-template-columns: 1fr !important;
    }

    /* Create a placeholder for the slideshow space */
    .slideshow-placeholder {
      height: 400px;
      width: 100%;
    }

    .slideshow-container {
      position: fixed !important;
      top: 60px !important; /* Account for navbar height */
      left: 0 !important;
      right: 0 !important;
      width: 100vw !important;
      height: 400px !important; /* Fixed height matching placeholder */
      max-width: none !important;
      margin: 0 !important;
      padding: 0 !important;
      z-index: 999999 !important;
      box-sizing: border-box;
      overflow: hidden;
      background-color: #000000bb !important;
      /* Additional properties to ensure it stays on top */
      isolation: isolate !important;
      transform: translateZ(0) !important; /* Force hardware acceleration */
    }

    /* Prevent any content from appearing above slideshow */
    .bd-main .bd-content {
      position: relative !important;
      z-index: 1 !important;
    }

    /* Ensure navbar doesn't interfere */
    .bd-header {
      z-index: 1000000 !important;
    }

    /* Style only slideshow images, not all images on the page */
    .slideshow-container img {
      width: 100% !important;
      height: 100% !important;
      object-fit: cover !important; /* Fill container, crop edges if needed */
      object-position: center !important; /* Center the image when cropping */
      display: block !important;
      background-color: transparent !important;
    }

    .mySlides {
      display: none;
      width: 100%;
      height: 100%;
      /* Center the content within each slide */
      display: flex;
      align-items: center;
      justify-content: center;
    }

    .fade {
      animation: fade 5000ms infinite;
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
      setTimeout(showSlides, 5000); // Change image every N seconds
    }
    </script>

.. <gallery>

.. title:: BuildAMol
   

.. ====================================
.. Welcome to BuildAMol's documentation
.. ====================================

.. image:: _resources/logo_large_light.png
   :class: only-light
   :width: 80%
   :align: center
   :alt: logo

.. image:: _resources/logo_large_dark.png
   :class: only-dark
   :width: 80%
   :align: center
   :alt: logo


`BuildAMol` (formerly Biobuild) is a fragment-based molecular assembly toolkit for the generation of atomic models for complex molecular structures.
It is designed to leverage the simplicity of python-coding and the power of fragment-based assembly to provide a slim and streamlined workflow.
Based on `biopython <http://biopython.org/wiki/Main_Page>`_ and accessible as a `python package`, `BuildAMol` not only offers a straightforward API to generate, manipulate, visualize, and export 3D structures of molecules, but also provides easy interfaces with other molecular modeling tools such as `RDKit <https://www.rdkit.org/docs/index.html>`_.



.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :hidden:

   whatfor
   installation
   tutorials
   documentation


.. grid:: 3


    .. grid-item-card::  Tutorials
        :link: tutorials
        :link-type: ref
        :link-alt: Tutorials

    .. grid-item-card::  API Documentation
        :link: apidocumentation
        :link-type: ref
        :link-alt: API Documentation    


    .. grid-item-card::  Extensions
        :link: buildamol.extensions
        :link-type: doc
        :link-alt: Extensions
