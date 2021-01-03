# shwfs
This project explores the data-reduction pipeline for processing spot-field
images acquired from a Shack-Hartmann wavefront sensor (WFS). All images
are from Roordalab's custom-built [WFS](http://roorda.vision.berkeley.edu/iq_res.htm#SHWS)

## Eye Model
The Arizona Eye Model serves as the foundation

| Effect | Assumption |
| ---------- | ------- |
| Corneal diameter | Dependent on unaccommodated ACD |
| Foveal diameter | 1 mm |
| Scleral thickness | Equal to difference in corneal sag at corneal edge |