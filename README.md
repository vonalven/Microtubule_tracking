# Bioimage informatics - EB3Track Project

The aim of our porject is to develop a tool, which detects, tracks, and measures the behavior of fluorescently labeled microtubules imaged by fluorescent microscopy techniques.

The ImageJ plugin that we have developed is a tri-modular tool. The first optional module permits to preprocess the input stack image. The second module detects microtubule seeds, tracks the growing microtubule ends and creates trajectories. The third module clusters these trajectories and fit models of dynamic behavior and is able to compute population statistics such as length distributions.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

For development purpose, the following extra libraries are required.

```
commons-math3-3.6.1.jar
```

### Installing

 - Download the **EB3T_Plugin.jar** release from [here](https://github.com/vonalven/Microtubule_tracking/releases)
 - Place it into the jars folder located into the ImageJ's *plugins* folder 
 - Restart ImageJ (or refresh menus)
 - Import a stack image and run *EB3Track Plugin*

## Running the plugin

Once launched, the plugin will ask for a de-noising pre-process step. 

Next steps require parameters the following parameters to be specified.

### Tracking - Parameters

 * *Threshold* - defines the minimum intensity with which to detect the points.
 * *Radius spot* - defines the minimum size with which to detect the points.
 * *Holes threshold* - factor which, weighted to the average intensity of a trajectory, allows the separation of points of intersection between two microtubules, characterized by an increase in brightness (in non-saturated images).
 * *Velocity max* - maximum distance to the next point to be considered in the same trajectory. It can be seen as the maximum velocity that a microtubule can have.  


### Clustering - Parameters

 * *Max gap* - defines the maximum gap in frames that can exist between the end and the start of 2 trajectories to be connected
 * *Min sub-track length* - defines the minimal length in terms of frames that a trajectory need to have. Otherwise it is discarded
 * *Forward cone half-opening* - defines the half-angle of the forward searching cone. The cone start from the last point of a reference trajectory; if some starting point of other trajectories is inside this cone, the corresponding trajectory is considered as a potential candidate for linking to the reference trajectory.
 * *Sedentary factor* - defines the sedentarity threshold of trajectories. If it is high, sedentary trajectories are discared, if it is zero, any trajectory is discarded by sedentarity. By sedentarity it is intended the spatial mouvement range covered by a trajectory during its life [in all the frames where it is defined].



## Built With

* [Java 1.8](http://www.dropwizard.io/1.0.2/docs/) - Java Compiler v1.8
* [Maven](https://maven.apache.org/) - Dependency Management

## Students

* **Pietro Airaghi** - [pietroairaghi](https://github.com/pietroairaghi)
* **Giacomo von Alvensleben** - [vonalven](https://github.com/vonalven)

## Professors

* **Arne Seitz** - [People@EPFL](https://people.epfl.ch/arne.seitz)
* **Daniel Sage** - [People@EPFL](https://people.epfl.ch/daniel.sage)

