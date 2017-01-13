# grandpopo
Python scripts for Grand Popo related processing.

This set of scripts handles the **processing of Grand Popo** shore-based video data:
- frame extraction from video file;
- rectification;
- filtering;
- output frames as images;
- motion estimation through Typhoon;
- output estimate as time-series and/or fields;
- create the figures.

## Requirements
- Scientific Python distribution including numpy, scipy, matplotlib, pandas. E.g. [Anaconda](https://www.continuum.io/downloads), [Canopy](https://www.enthought.com/products/canopy/);
- FFMPEG for frames extraction;
- [Typhoon estimator](http://www.pierrederian.net/typhoon.html).

## Experiment configuration
![Grand Popo experiment configuration][resources/grandpopo_config.jpg]





