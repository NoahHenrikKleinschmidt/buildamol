.. only:: html

    .. title:: BuildAMol | Fragment-based Molecular Assembly Toolkit


.. raw:: html

    <div class="hero-section">
        <div class="hero-content">
            
            <h2>Fragment-based Molecular Assembly Toolkit</h2>
            <p>Generate complex molecular structures with Python. Simple, powerful, and intuitive.</p>
            <div class="hero-buttons">
                <a href="installation.html" class="btn btn-primary">Get Started</a>
                <a href="tutorials.html" class="btn btn-secondary">View Tutorials</a>
            </div>
        </div>
        <div class="hero-slide">
            <img src="_static/gallery/nanotube.png" alt="Molecular structures">
        </div>
    </div>

A Python toolkit for generating atomic models of complex molecular structures
================================================================================

`BuildAMol` (formerly Biobuild) is designed to leverage the simplicity of Python coding and the power of **fragment-based assembly** to provide a streamlined molecular modeling workflow.

**Key Features:**

- ✨ **Fragment-based assembly**: Build complex molecules from molecular fragments
- 🧬 **Biology support**: Dedicated extensions for proteins, glycosylation, and docking
- ⚛️ **Chemistry-aware**: Automatic handling of reactions, functional groups, and geometry
- 🔧 **Materials science**: Tools for polymers, nanotubes, and metal-organic frameworks
- 📊 **Integration**: Seamless interfaces with RDKit, OpenMM, and molecular dynamics tools
- 📓 **Documentation**: Interactive Jupyter notebooks and comprehensive API reference

.. raw:: html

    <div class="feature-cards">

.. grid:: 2

    .. grid-item-card:: 🎓 Learn by Example
        :link: tutorials
        :link-type: ref
        :class-title: card-title-emphasized

        Explore interactive tutorials organized by use case: basic usage, chemistry, biology, materials science, and molecular dynamics.

    .. grid-item-card:: 📖 API Reference
        :link: apidocumentation
        :link-type: ref
        :class-title: card-title-emphasized

        Comprehensive API documentation for the core library and all extensions, with detailed docstrings and examples.

    .. grid-item-card:: 🏗️ What is BuildAMol?
        :link: whatfor
        :link-type: doc
        :class-title: card-title-emphasized

        Understand the philosophy, core concepts, and design decisions behind BuildAMol.

    .. grid-item-card:: 📦 Installation
        :link: installation
        :link-type: doc
        :class-title: card-title-emphasized

        Quick setup guide with all requirements and optional dependencies.

.. raw:: html

    </div>

Showcase
========

BuildAMol can model a diverse range of molecular systems:

.. raw:: html

        <div class="showcase-rotator-wrap">
                <img id="showcase-rotator" class="showcase-rotator-image" src="_static/gallery/complex.png" alt="BuildAMol showcase image" />
        </div>

.. raw:: html

        <script>
        (function () {
            const images = [
                "_static/gallery/complex.png",
                "_static/gallery/peptide.png",
                "_static/gallery/nanotube.png",
                "_static/gallery/MOF-top.png"
            ];
            const el = document.getElementById("showcase-rotator");
            if (!el) return;
            let idx = 0;
            setInterval(function () {
                idx = (idx + 1) % images.length;
                el.classList.add("is-fading");
                setTimeout(function () {
                    el.src = images[idx];
                    el.classList.remove("is-fading");
                }, 180);
            }, 4200);
        })();
        </script>


.. toctree::
   :maxdepth: 2
   :caption: Documentation
   :hidden:

   whatfor
   installation
   tutorials
   Main library <documentation>
   Extensions <buildamol.extensions>
