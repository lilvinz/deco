deco
====

This is a python implementation of the original gfdeco.f fortran ZHL16B based decompression model by Eric Baker.
The main idea is to have an implementation which matches the original fortran implementation without any difference. Most decompression programs don't match the original fortran implementation and thus you never know if a program is erroneous or not.

I used ipython as it is a very intuitive and easy to use environment for doing scientific calculations.  
See: [ipython.org](http://ipython.org/notebook.html)

IPython notebooks can be rendered online:
* Validation against original fortran code: [zhl16b_validate.ipynb](http://nbviewer.ipython.org/github/lilvinz/deco/blob/master/zhl16b_validate.ipynb)
* Example dive profile analysis: [zhl16b_analyze_uddf.ipynb](http://nbviewer.ipython.org/github/lilvinz/deco/blob/master/zhl16b_analyze_uddf.ipynb)
* Example dive profile uing [bokeh](http://bokeh.pydata.org): [zhl16b_bokeh.ipynb](http://nbviewer.ipython.org/github/lilvinz/deco/blob/master/zhl16b_bokeh.ipynb)

The input file format is [uddf](http://uddf.org)

Note that the current state of this is work in progress and the main goal of matching gfdeco.f results has NOT been reached. This is due to the impossibility to choose floating point precision in python. The only possible solution i currently see is to call into an external c library to do the actual math which is really not satisfactory.
